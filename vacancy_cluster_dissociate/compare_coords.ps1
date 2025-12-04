param()
$canonical = 'C:\Users\Administrator\Documents\xwechat_files\wxid_ccnsk401f7tu22_5a73\msg\file\2025-11\df\voids.0'
$corrected = 'C:\Users\Administrator\Documents\xwechat_files\wxid_ccnsk401f7tu22_5a73\msg\file\2025-11\df\corrected_voids_relax\voids-relax.0.dump'
function ParseDump($path){
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
        # try parse first three tokens as doubles (declare targets first)
        $x = 0.0; $y = 0.0; $z = 0.0
        $ok1 = $false; $ok2 = $false; $ok3 = $false
        if ($parts.Count -ge 1) { $ok1 = [double]::TryParse($parts[0],[ref]$x) }
        if ($parts.Count -ge 2) { $ok2 = [double]::TryParse($parts[1],[ref]$y) }
        if ($parts.Count -ge 3) { $ok3 = [double]::TryParse($parts[2],[ref]$z) }
        if ($ok1 -and $ok2 -and $ok3) { $coords += ,@($x,$y,$z) }
    }
    return $coords
}
$c = ParseDump $canonical
$r = ParseDump $corrected
Write-Host "Canonical parsed count: $($c.Count)"
Write-Host "Corrected parsed count: $($r.Count)"
$min = [Math]::Min($c.Count,$r.Count)
Write-Host "Printing first 12 canonical coords:"
for ($i=0;$i -lt [Math]::Min(12,$c.Count);$i++){
    $t=$c[$i]
    $line = $i.ToString() + ': ' + [string]::Format('{0:F10}',$t[0]) + ' ' + [string]::Format('{0:F10}',$t[1]) + ' ' + [string]::Format('{0:F10}',$t[2])
    Write-Host $line
}
Write-Host ""
Write-Host "Printing first 24 corrected coords:"
for ($i=0;$i -lt [Math]::Min(24,$r.Count);$i++){
    $t=$r[$i]
    $line = $i.ToString() + ': ' + [string]::Format('{0:F10}',$t[0]) + ' ' + [string]::Format('{0:F10}',$t[1]) + ' ' + [string]::Format('{0:F10}',$t[2])
    Write-Host $line
}
Write-Host ""
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
