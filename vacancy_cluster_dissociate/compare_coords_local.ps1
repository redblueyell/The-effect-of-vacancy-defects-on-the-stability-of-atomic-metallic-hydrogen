param(
    [Parameter(Mandatory=$true)][string]$canonicalPath,
    [Parameter(Mandatory=$true)][string]$correctedPath
)
function ParseDump($path){
    if (-not (Test-Path $path)) { Write-Error "File not found: $path"; return @() }
    $lines = Get-Content -LiteralPath $path
    $start = -1
    for ($li=0;$li -lt $lines.Count;$li++){
        if ($lines[$li] -match '^ITEM:\s*ATOMS') { $start = $li; break }
    }
    if ($start -lt 0) { return @() }
    $coords = @()
    for ($i = $start + 1; $i -lt $lines.Count; $i++){
        $ln = $lines[$i].Trim()
        if ([string]::IsNullOrWhiteSpace($ln)) { continue }
        $parts = $ln -split '\s+'
        $x = 0.0; $y = 0.0; $z = 0.0
        $ok1 = $false; $ok2 = $false; $ok3 = $false
        if ($parts.Count -ge 1) { $ok1 = [double]::TryParse($parts[0],[ref]$x) }
        if ($parts.Count -ge 2) { $ok2 = [double]::TryParse($parts[1],[ref]$y) }
        if ($parts.Count -ge 3) { $ok3 = [double]::TryParse($parts[2],[ref]$z) }
        if ($ok1 -and $ok2 -and $ok3) { $coords += ,@($x,$y,$z) }
    }
    return $coords
}
$c = ParseDump $canonicalPath
$r = ParseDump $correctedPath
Write-Host "Canonical parsed count: $($c.Count)"
Write-Host "Corrected parsed count: $($r.Count)"
Write-Host "Printing first 12 canonical coords:"
for ($i=0;$i -lt [Math]::Min(12,$c.Count);$i++){
    $t=$c[$i]
    Write-Host ($i.ToString() + ': ' + ('{0:F10}' -f $t[0]) + ' ' + ('{0:F10}' -f $t[1]) + ' ' + ('{0:F10}' -f $t[2]))
}
Write-Host ""
Write-Host "Printing first 24 corrected coords:"
for ($i=0;$i -lt [Math]::Min(24,$r.Count);$i++){
    $t=$r[$i]
    Write-Host ($i.ToString() + ': ' + ('{0:F10}' -f $t[0]) + ' ' + ('{0:F10}' -f $t[1]) + ' ' + ('{0:F10}' -f $t[2]))
}
Write-Host ""
$min = [Math]::Min($c.Count,$r.Count)
if ($min -gt 0) {
    $maxDiff = 0.0
    $sumDiff = 0.0
    for ($i=0;$i -lt [Math]::Min(100,$min);$i++){
        $tc = $c[$i]; $tr = $r[$i]
        $dx = [math]::Abs($tc[0]-$tr[0]); $dy=[math]::Abs($tc[1]-$tr[1]); $dz=[math]::Abs($tc[2]-$tr[2]);
        $d = [math]::Sqrt($dx*$dx+$dy*$dy+$dz*$dz)
        if ($d -gt $maxDiff) { $maxDiff = $d }
        $sumDiff += $d
    }
    $avg = $sumDiff / [Math]::Min(100,$min)
    Write-Host "Compared first $([Math]::Min(100,$min)) atoms: max dist = $([string]::Format('{0:F6}',$maxDiff)), avg = $([string]::Format('{0:F6}',$avg))"
}
