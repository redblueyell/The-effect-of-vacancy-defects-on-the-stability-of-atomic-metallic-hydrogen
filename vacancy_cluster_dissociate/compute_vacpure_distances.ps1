<#
compute_vacpure_distances.ps1

Compute pairwise distances for the 6 atoms in each vac-pure dump file
located in ./corrected_s1. Produces per-file `<basename>_distances.txt` files
with lines `distance i-j` (distance rounded to 2 decimals) and a summary
`dist_counts_per_1000.txt` that lists for each step (multiple of 1000)
how many of the 15 distances are greater than 3.0.

Usage:
  powershell -ExecutionPolicy Bypass -File .\compute_vacpure_distances.ps1

#>

[CmdletBinding()]
param(
    [string]$dir = (Join-Path (Get-Location) 'corrected_s1'),
    [double]$threshold = 3.0
)

function ToDouble([string]$s){
    $inv = [System.Globalization.CultureInfo]::InvariantCulture
    $d = 0.0
    if ([double]::TryParse($s, [System.Globalization.NumberStyles]::Float, $inv, [ref]$d)) { return $d }
    return $null
}

if (-not (Test-Path $dir)) { Write-Error "Directory not found: $dir"; exit 1 }

$inv = [System.Globalization.CultureInfo]::InvariantCulture
$files = Get-ChildItem -Path $dir -File | Where-Object { $_.Name -match '^vac-pure\.(\d+)\.dump$' } | Sort-Object Name
if (-not $files) { Write-Host "No vac-pure dump files found in $dir"; exit 0 }

$summary = @()

foreach ($f in $files) {
    $m = [regex]::Match($f.Name, '^vac-pure\.(\d+)\.dump$')
    if (-not $m.Success) { continue }
    $step = [int]$m.Groups[1].Value

    $lines = Get-Content -Path $f.FullName
    # robustly find the ITEM: ATOMS line index (case-insensitive)
    $startIndex = -1
    for ($i = 0; $i -lt $lines.Count; $i++) {
        if ($lines[$i] -match '^\s*ITEM:\s*ATOMS') { $startIndex = $i + 1; break }
    }
    if ($startIndex -lt 0) {
        Write-Warning "No 'ITEM: ATOMS' found in $($f.Name) - skipping"
        continue
    }
    # lines after ITEM: ATOMS
    $atomLines = @()
    if ($startIndex -lt $lines.Count) { $atomLines = $lines[$startIndex..($lines.Count - 1)] | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne '' } }

    # parse first 6 coordinate triples (assume they are the 6 target atoms in appearance order)
    $coords = @()
    foreach ($ln in $atomLines) {
        $parts = -split $ln
        if ($parts.Count -ge 3) {
            $x = ToDouble($parts[0]); $y = ToDouble($parts[1]); $z = ToDouble($parts[2])
            if ($x -ne $null -and $y -ne $null -and $z -ne $null) { $coords += ,@($x,$y,$z) }
        }
        if ($coords.Count -ge 6) { break }
    }
    if ($coords.Count -ne 6) { Write-Warning "File $($f.Name) parsed $($coords.Count) atoms (expected 6) - skipping"; continue }

    # compute pairwise distances
    $pairs = @()
    for ($i=0; $i -lt 6; $i++) {
        for ($j=$i+1; $j -lt 6; $j++) {
            $dx = $coords[$i][0] - $coords[$j][0]
            $dy = $coords[$i][1] - $coords[$j][1]
            $dz = $coords[$i][2] - $coords[$j][2]
            $d = [math]::Sqrt($dx*$dx + $dy*$dy + $dz*$dz)
            $pairs += [pscustomobject]@{ Distance = [double]$d; Pair = "{0}-{1}" -f ($i+1),($j+1) }
        }
    }

    # sort ascending by distance
    $pairs = $pairs | Sort-Object Distance

    # write per-file distances: include coordinates + 15 distances lines
    $outLines = @()
    $outLines += "# File: $($f.Name)"
    $outLines += "# Step: $step"
    $outLines += "# Coordinates (index x y z)"
    for ($i=0; $i -lt 6; $i++) { $outLines += "{0} {1} {2} {3}" -f ($i+1), $coords[$i][0].ToString('F10',$inv), $coords[$i][1].ToString('F10',$inv), $coords[$i][2].ToString('F10',$inv) }
    $outLines += "# Distances (sorted ascending): format -> distance i-j"
    foreach ($p in $pairs) { $outLines += ($p.Distance.ToString('F2',$inv) + ' ' + $p.Pair) }

    $outPath = Join-Path $dir (($f.BaseName) + '_distances.txt')
    $outLines | Set-Content -Path $outPath -Encoding UTF8

    # for summary: only count files where step % 1000 == 0
    if (($step % 1000) -eq 0) {
        $countGt = ($pairs | Where-Object { $_.Distance -gt $threshold }).Count
        $summary += [pscustomobject]@{ Step = $step; CountGreaterThanThreshold = $countGt }
    }
}

# write summary sorted by step
$summary = $summary | Sort-Object Step
$summaryPath = Join-Path $dir 'dist_counts_per_1000.txt'
$hdr = 'step count_gt_3'
$out = @($hdr)
foreach ($s in $summary) { $out += ("{0} {1}" -f $s.Step, $s.CountGreaterThanThreshold) }
$out | Set-Content -Path $summaryPath -Encoding UTF8

Write-Host "Done. Wrote per-file distances and summary to $summaryPath" -ForegroundColor Green
