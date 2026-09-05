$ErrorActionPreference = "Stop"

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$julia = "C:\Users\Stadnichuk\AppData\Local\Programs\Julia-1.12.6\bin\julia.exe"
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$runRoot = Join-Path $PSScriptRoot (Join-Path "runs" $stamp)
$hndpOutput = Join-Path $runRoot "network_design"
$bobilibOutput = Join-Path $runRoot "bobilib\results.csv"
$launcherLog = Join-Path $runRoot "sequential_launcher.log"

New-Item -ItemType Directory -Force -Path $hndpOutput, (Split-Path $bobilibOutput) | Out-Null

Push-Location $repoRoot
try {
    "Starting HNDP at $(Get-Date -Format o)" | Tee-Object -FilePath $launcherLog
    $env:JUBIC_HNDP_OUTPUT = $hndpOutput
    & $julia --project=. benchmark/network_design/run_experiment.jl 2>&1 |
        Tee-Object -FilePath $launcherLog -Append
    if ($LASTEXITCODE -ne 0) {
        throw "HNDP experiment exited with code $LASTEXITCODE. BOBILib was not started."
    }

    "HNDP completed at $(Get-Date -Format o); starting full BOBILib at $(Get-Date -Format o)" |
        Tee-Object -FilePath $launcherLog -Append
    $env:JUBIC_BOBILIB_OUTPUT = $bobilibOutput
    Remove-Item Env:JUBIC_BOBILIB_MAX_INSTANCES -ErrorAction SilentlyContinue
    & $julia --project=. benchmark/bobilib/run_experiment.jl 2>&1 |
        Tee-Object -FilePath $launcherLog -Append
    if ($LASTEXITCODE -ne 0) {
        throw "BOBILib experiment exited with code $LASTEXITCODE."
    }
    "BOBILib completed at $(Get-Date -Format o)" | Tee-Object -FilePath $launcherLog -Append
}
finally {
    Remove-Item Env:JUBIC_HNDP_OUTPUT -ErrorAction SilentlyContinue
    Remove-Item Env:JUBIC_BOBILIB_OUTPUT -ErrorAction SilentlyContinue
    Pop-Location
}
