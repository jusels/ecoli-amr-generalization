if (Test-Path data\PanPred) {
    Write-Host "data\PanPred already exists, skipping download"
} else {
    New-Item -ItemType Directory -Force data | Out-Null
    git clone --depth 1 https://github.com/DaneshMoradigaravand/PanPred.git data/PanPred
}