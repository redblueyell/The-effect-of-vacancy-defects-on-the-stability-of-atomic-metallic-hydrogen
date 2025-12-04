<#
convert_voidsrelax_every_k_steps.ps1

Select files named like `voids-relax.<step>.cfg` in a target directory,
convert those where step % K == 0 to LAMMPS dump (cart = H0_diag*(frac-0.5)),
and write outputs to a `corrected_voids_relax` subdirectory under the target directory.

Usage:
  powershell -ExecutionPolicy Bypass -File .\convert_voidsrelax_every_k_steps.ps1 -dir "C:\path\to\df" -step 1000
#>

[CmdletBinding()]
param(
    [Parameter(Mandatory=$true)][string]$dir,
    [Parameter(Mandatory=$false)][int]$step = 1000
)

function ToDouble([string]$s){
    $inv = [System.Globalization.CultureInfo]::InvariantCulture
    $d = 0.0
    if ([double]::TryParse($s, [System.Globalization.NumberStyles]::Float, $inv, [ref]$d)) { return $d }
    return $null
}

if (-not (Test-Path $dir)) { throw "Directory not found: $dir" }

$outDir = Join-Path $dir 'corrected_voids_relax'
if (-not (Test-Path $outDir)) { New-Item -Path $outDir -ItemType Directory | Out-Null }

$files = Get-ChildItem -Path $dir -File | Where-Object { $_.Name -match '^voids-relax\.(\d+)(?:\.cfg)?$' }
if (-not $files -or $files.Count -eq 0) { Write-Host "No voids-relax files found in $dir" -ForegroundColor Yellow; return }

$processed = @()
foreach ($f in $files) {
    $m = [regex]::Match($f.Name, '^voids-relax\.(\d+)(?:\.cfg)?$')
    if (-not $m.Success) { continue }
    $idx = [int]$m.Groups[1].Value
    if (($idx % $step) -ne 0) { continue }

    Write-Host "Converting $($f.Name) (step=$idx) ..."
    $text = Get-Content $f.FullName -Raw

    # Number of particles
    $mN = [regex]::Match($text, 'Number of particles\s*=\s*(\d+)', 'IgnoreCase')
    if (-not $mN.Success) { Write-Warning "No particle count in $($f.Name) - skipping"; continue }
    $N = [int]$mN.Groups[1].Value

    # H0 diagonal
    $m11 = [regex]::Match($text, 'H0\(1,1\)\s*=\s*([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)', 'IgnoreCase')
    $m22 = [regex]::Match($text, 'H0\(2,2\)\s*=\s*([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)', 'IgnoreCase')
    $m33 = [regex]::Match($text, 'H0\(3,3\)\s*=\s*([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)', 'IgnoreCase')
    if (-not ($m11.Success -and $m22.Success -and $m33.Success)) { Write-Warning "H0 diag missing in $($f.Name) - skipping"; continue }
    $h0x = ToDouble($m11.Groups[1].Value); $h0y = ToDouble($m22.Groups[1].Value); $h0z = ToDouble($m33.Groups[1].Value)
    if ($h0x -eq $null -or $h0y -eq $null -or $h0z -eq $null) { Write-Warning "H0 parse null in $($f.Name) - skipping"; continue }

    # parse numeric triples
    $coords = @()
    foreach ($m2 in [regex]::Matches($text, '([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)\s+([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)')) {
        $x = ToDouble($m2.Groups[1].Value); $y = ToDouble($m2.Groups[2].Value); $z = ToDouble($m2.Groups[3].Value)
        if ($x -ne $null -and $y -ne $null -and $z -ne $null) { $coords += ,@($x,$y,$z) }
    }
    if ($coords.Count -lt $N) { Write-Warning "Found $($coords.Count) triples, expected $N in $($f.Name) - skipping"; continue }
    $selected = $coords[-$N..-1]

    # convert fractional -> cartesian using cart = H0_diag * (frac - 0.5)
    $inv = [System.Globalization.CultureInfo]::InvariantCulture
    $outLines = @()
    $outLines += 'ITEM: TIMESTEP'; $outLines += '0'
    $outLines += 'ITEM: NUMBER OF ATOMS'; $outLines += $N.ToString()
    $outLines += 'ITEM: BOX BOUNDS pp pp pp'
    $outLines += ((-1.0 * $h0x/2.0).ToString('F10',$inv) + ' ' + ($h0x/2.0).ToString('F10',$inv))
    $outLines += ((-1.0 * $h0y/2.0).ToString('F10',$inv) + ' ' + ($h0y/2.0).ToString('F10',$inv))
    $outLines += ((-1.0 * $h0z/2.0).ToString('F10',$inv) + ' ' + ($h0z/2.0).ToString('F10',$inv))
    $outLines += 'ITEM: ATOMS x y z'
    foreach ($t in $selected) {
        $fx = [double]$t[0]; $fy = [double]$t[1]; $fz = [double]$t[2]
        $cx = ($fx - 0.5) * $h0x; $cy = ($fy - 0.5) * $h0y; $cz = ($fz - 0.5) * $h0z
        $outLines += ($cx.ToString('F10',$inv) + ' ' + $cy.ToString('F10',$inv) + ' ' + $cz.ToString('F10',$inv))
    }

    $outName = "voids-relax.$idx.dump"
    $outPath = Join-Path $outDir $outName
    $outLines | Set-Content -Path $outPath -Encoding UTF8
    Write-Host "Wrote $outPath"
    $processed += $outPath
}

Write-Host "Done. Generated $($processed.Count) files in $outDir" -ForegroundColor Green
