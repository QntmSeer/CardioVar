# Setup Windows Task Scheduler for Daily Backup
# This script creates a scheduled task to run the backup script every day

# Configuration
$TaskName = "CardioVar Daily Backup"
$ScriptPath = "C:\Users\Gebruiker\Documents\Cardiovar\daily_backup.ps1"
$BackupTime = "02:00"  # 2:00 AM daily

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "Setting Up Daily Backup Automation" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

# Check if script exists
if (-not (Test-Path $ScriptPath)) {
    Write-Host "[ERROR] Backup script not found at: $ScriptPath" -ForegroundColor Red
    exit 1
}

# Remove existing task if it exists
$ExistingTask = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
if ($ExistingTask) {
    Write-Host "Removing existing scheduled task..." -ForegroundColor Yellow
    Unregister-ScheduledTask -TaskName $TaskName -Confirm:$false
}

# Create the scheduled task action
$Action = New-ScheduledTaskAction -Execute "PowerShell.exe" `
    -Argument "-NoProfile -ExecutionPolicy Bypass -File `"$ScriptPath`""

# Create the trigger (daily at specified time)
$Trigger = New-ScheduledTaskTrigger -Daily -At $BackupTime

# Create task settings
$Settings = New-ScheduledTaskSettingsSet `
    -AllowStartIfOnBatteries `
    -DontStopIfGoingOnBatteries `
    -StartWhenAvailable `
    -RunOnlyIfNetworkAvailable:$false

# Create the principal (run with highest privileges)
$Principal = New-ScheduledTaskPrincipal -UserId "$env:USERDOMAIN\$env:USERNAME" `
    -LogonType Interactive -RunLevel Highest

# Register the scheduled task
try {
    Register-ScheduledTask -TaskName $TaskName `
        -Action $Action `
        -Trigger $Trigger `
        -Settings $Settings `
        -Principal $Principal `
        -Description "Automatically backs up CardioVar project daily at $BackupTime" | Out-Null
    
    Write-Host "`n[SUCCESS] Daily backup scheduled!" -ForegroundColor Green
    Write-Host "`nTask Details:" -ForegroundColor Cyan
    Write-Host "  Task Name: $TaskName" -ForegroundColor White
    Write-Host "  Schedule: Daily at $BackupTime" -ForegroundColor White
    Write-Host "  Script: $ScriptPath" -ForegroundColor White
    Write-Host "  Backup Location: C:\Users\Gebruiker\Documents\CardioVar_Backups" -ForegroundColor White
    Write-Host "  Retention: Last 30 backups" -ForegroundColor White
    
    Write-Host "`nYou can manage this task in:" -ForegroundColor Yellow
    Write-Host "  Task Scheduler > Task Scheduler Library > $TaskName" -ForegroundColor White
    
}
catch {
    Write-Host "`n[ERROR] Failed to create scheduled task: $_" -ForegroundColor Red
    exit 1
}

# Test the backup script now (optional)
Write-Host "`n========================================" -ForegroundColor Cyan
$Response = Read-Host "Would you like to run a test backup now? (Y/N)"
if ($Response -eq 'Y' -or $Response -eq 'y') {
    Write-Host "`nRunning test backup..." -ForegroundColor Green
    & $ScriptPath
}
else {
    Write-Host "`nSkipping test backup. First backup will run at $BackupTime" -ForegroundColor Yellow
}

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "Setup Complete!" -ForegroundColor Green
Write-Host "========================================`n" -ForegroundColor Cyan
