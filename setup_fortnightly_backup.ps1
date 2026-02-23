# Setup Windows Task Scheduler for Fortnightly Backup
# This script creates a scheduled task to run the backup script every 2 weeks

# Configuration
$OldTaskName = "CardioVar Daily Backup"
$NewTaskName = "CardioVar Fortnightly Backup"
$ScriptPath = "C:\Users\Gebruiker\Documents\Cardiovar\daily_backup.ps1"
$BackupTime = "02:00"  # 2:00 AM

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "Setting Up Fortnightly Backup Automation" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

# Check if script exists
if (-not (Test-Path $ScriptPath)) {
    Write-Host "[ERROR] Backup script not found at: $ScriptPath" -ForegroundColor Red
    exit 1
}

# 1. Clean up old daily task
$OldTask = Get-ScheduledTask -TaskName $OldTaskName -ErrorAction SilentlyContinue
if ($OldTask) {
    Write-Host "Removing old daily backup task..." -ForegroundColor Yellow
    Unregister-ScheduledTask -TaskName $OldTaskName -Confirm:$false
}

# 2. Clean up existing fortnightly task (to update it)
$ExistingTask = Get-ScheduledTask -TaskName $NewTaskName -ErrorAction SilentlyContinue
if ($ExistingTask) {
    Write-Host "Removing existing fortnightly task..." -ForegroundColor Yellow
    Unregister-ScheduledTask -TaskName $NewTaskName -Confirm:$false
}

# Create the scheduled task action
$Action = New-ScheduledTaskAction -Execute "PowerShell.exe" `
    -Argument "-NoProfile -ExecutionPolicy Bypass -File `"$ScriptPath`""

# Create the trigger (Weekly, every 2 weeks, on Sunday)
$Trigger = New-ScheduledTaskTrigger -Weekly -WeeksInterval 2 -DaysOfWeek Sunday -At $BackupTime

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
    Register-ScheduledTask -TaskName $NewTaskName `
        -Action $Action `
        -Trigger $Trigger `
        -Settings $Settings `
        -Principal $Principal `
        -Description "Automatically backs up CardioVar project every 2 weeks at $BackupTime" | Out-Null
    
    Write-Host "`n[SUCCESS] Fortnightly backup scheduled!" -ForegroundColor Green
    Write-Host "`nTask Details:" -ForegroundColor Cyan
    Write-Host "  Task Name: $NewTaskName" -ForegroundColor White
    Write-Host "  Schedule: Every 2 weeks on Sunday at $BackupTime" -ForegroundColor White
    Write-Host "  Script: $ScriptPath" -ForegroundColor White
    Write-Host "  Retention: Last 30 backups" -ForegroundColor White
    
    Write-Host "`nYou can manage this task in:" -ForegroundColor Yellow
    Write-Host "  Task Scheduler > Task Scheduler Library > $NewTaskName" -ForegroundColor White
    
}
catch {
    Write-Host "`n[ERROR] Failed to create scheduled task: $_" -ForegroundColor Red
    exit 1
}

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "Setup Complete!" -ForegroundColor Green
Write-Host "========================================`n" -ForegroundColor Cyan
