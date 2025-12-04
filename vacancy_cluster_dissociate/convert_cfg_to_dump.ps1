<#
convert_cfg_to_dump.ps1

Batch-convert all `*.cfg` in the current directory to LAMMPS dump files.
Conversion rule: cart = H0_diag * (frac - 0.5)
Outputs are written to `./corrected/<basename>.dump`.

Usage:
  powershell -ExecutionPolicy Bypass -File .\convert_cfg_to_dump.ps1

#>

[CmdletBinding()]
param()

function ToDouble([string]$s){
    $inv = [System.Globalization.CultureInfo]::InvariantCulture
    $d = 0.0
    if ([double]::TryParse($s, [System.Globalization.NumberStyles]::Float, $inv, [ref]$d)) { return $d }
    return $null
}

$cwd = Get-Location
$cfgFiles = Get-ChildItem -Path $cwd -Filter '*.cfg' -File -ErrorAction SilentlyContinue
if (-not $cfgFiles -or $cfgFiles.Count -eq 0) {
    Write-Host "No .cfg files found in $cwd" -ForegroundColor Yellow
    return
}

$outDir = Join-Path $cwd 'corrected'
if (-not (Test-Path $outDir)) { New-Item -Path $outDir -ItemType Directory | Out-Null }

foreach ($file in $cfgFiles) {
    Write-Host "Processing $($file.Name) ..."
    $text = Get-Content $file.FullName -Raw

    # Number of particles
    $mN = [regex]::Match($text, 'Number of particles\s*=\s*(\d+)', 'IgnoreCase')
    if (-not $mN.Success) {
        Write-Warning "Could not find 'Number of particles' in $($file.Name) - skipping"
        continue
    }
    $N = [int]$mN.Groups[1].Value

    # H0 diagonal elements
    $m11 = [regex]::Match($text, 'H0\(1,1\)\s*=\s*([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)', 'IgnoreCase')
    $m22 = [regex]::Match($text, 'H0\(2,2\)\s*=\s*([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)', 'IgnoreCase')
    $m33 = [regex]::Match($text, 'H0\(3,3\)\s*=\s*([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)', 'IgnoreCase')
    if (-not ($m11.Success -and $m22.Success -and $m33.Success)) {
        Write-Warning "Could not parse H0 diagonal for $($file.Name) - skipping"
        continue
    }
    $h0x = ToDouble($m11.Groups[1].Value)
    $h0y = ToDouble($m22.Groups[1].Value)
    $h0z = ToDouble($m33.Groups[1].Value)

    if ($h0x -eq $null -or $h0y -eq $null -or $h0z -eq $null) {
        Write-Warning "H0 parse produced nulls for $($file.Name) - skipping"
        continue
    }

    # Parse coordinates robustly by scanning lines: for each line, if it contains
    # at least three floats, take the last three floats on that line as a fractional
    # coordinate triple. This avoids picking up unrelated numeric triples elsewhere
    # in the file and preserves per-line ordering.
    $lines = Get-Content $file.FullName
    $floatPattern = '([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)'
    $coords = @()
    $bad = @()
    foreach ($ln in $lines) {
        $matches = [regex]::Matches($ln, $floatPattern)
        if ($matches.Count -ge 3) {
            # prefer the FIRST three floats on the line (many .cfg files store
            # fractional coords at the start of the atom line, followed by
            # indices/charges). Using the first three avoids picking up those
            # trailing integers.
            $a = ToDouble($matches[0].Value)
            $b = ToDouble($matches[1].Value)
            $c = ToDouble($matches[2].Value)
            if ($a -ne $null -and $b -ne $null -and $c -ne $null) {
                # Accept triples that look like fractional coordinates (tolerant range)
                if (($a -ge -1.0 -and $a -le 2.0) -and ($b -ge -1.0 -and $b -le 2.0) -and ($c -ge -1.0 -and $c -le 2.0)) {
                    $coords += ,@($a,$b,$c)
                } else {
                    $bad += ,@{line=$ln; triple=@($a,$b,$c)}
                }
            }
        }
    }

    if ($coords.Count -lt $N) {
        Write-Warning "Found only $($coords.Count) coordinate triples in $($file.Name), expected $N - skipping"
        # Write a small debug file containing the first few rejected numeric triples to help diagnose
        if ($bad.Count -gt 0) {
            $dbgPath = Join-Path $outDir ($file.BaseName + '.debug.txt')
            $dbg = @()
            $dbg += "Rejected numeric triples (first 50):"
            $i = 0
            foreach ($b in $bad) {
                $dbg += ("LINE: " + $b.line)
                $dbg += ("TRIPLE: " + ($b.triple -join ", "))
                $i++
                if ($i -ge 50) { break }
            }
            $dbg | Set-Content -Path $dbgPath -Encoding UTF8
            Write-Host "Wrote debug file: $dbgPath" -ForegroundColor Yellow
        }
        continue
    }

    # Take the last N coordinate lines (most cfg files have coordinates at the end)
    $selected = $coords[-$N..-1]

    # Convert fractional -> cartesian using cart = H0 * (frac - 0.5)
    $inv = [System.Globalization.CultureInfo]::InvariantCulture
    $outLines = @()
    # attempt to infer timestep from filename (e.g. voids-relax.100000 -> 100000)
    $timestep = 0
    $mstep = [regex]::Match($file.BaseName, '(\d+)$')
    if ($mstep.Success) { $timestep = [int]$mstep.Groups[1].Value }
    $outLines += 'ITEM: TIMESTEP'
    $outLines += $timestep.ToString()
    # record source filename for traceability
    $outLines += ("# source: " + $file.Name)
    $outLines += 'ITEM: NUMBER OF ATOMS'
    $outLines += $N.ToString()
    $outLines += 'ITEM: BOX BOUNDS pp pp pp'
    $outLines += ((-1.0 * $h0x/2.0).ToString('F10',$inv) + ' ' + ($h0x/2.0).ToString('F10',$inv))
    $outLines += ((-1.0 * $h0y/2.0).ToString('F10',$inv) + ' ' + ($h0y/2.0).ToString('F10',$inv))
    $outLines += ((-1.0 * $h0z/2.0).ToString('F10',$inv) + ' ' + ($h0z/2.0).ToString('F10',$inv))
    $outLines += 'ITEM: ATOMS x y z'

    foreach ($t in $selected) {
        $fx = [double]$t[0]; $fy = [double]$t[1]; $fz = [double]$t[2]
        $cx = ($fx - 0.5) * $h0x
        $cy = ($fy - 0.5) * $h0y
        $cz = ($fz - 0.5) * $h0z
        $outLines += ($cx.ToString('F10',$inv) + ' ' + $cy.ToString('F10',$inv) + ' ' + $cz.ToString('F10',$inv))
    }

    $outPath = Join-Path $outDir ($file.BaseName + '.dump')
    $outLines | Set-Content -Path $outPath -Encoding UTF8
    Write-Host "Wrote $outPath (atoms=$N)"
}

Write-Host 'Done.' -ForegroundColor Green