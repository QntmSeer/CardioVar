# CardioVar Daily Backup Script
# Automatically creates timestamped backups of the entire project

# Configuration
$ProjectPath = "C:\Users\Gebruiker\Documents\Cardiovar"
$BackupRootPath = "C:\Users\Gebruiker\Documents\CardioVar_Backups"
$MaxBackups = 30  # Keep last 30 days of backups

# Create backup root directory if it doesn't exist
if (-not (Test-Path $BackupRootPath)) {
    New-Item -ItemType Directory -Path $BackupRootPath -Force | Out-Null
    Write-Host "Created backup directory: $BackupRootPath" -ForegroundColor Green
}

# Generate timestamp for backup folder
$Timestamp = Get-Date -Format "yyyy-MM-dd_HHmm"
$BackupPath = Join-Path $BackupRootPath "Cardiovar_Backup_$Timestamp"

# Create the backup
Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "CardioVar Daily Backup" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Source: $ProjectPath" -ForegroundColor Yellow
Write-Host "Destination: $BackupPath" -ForegroundColor Yellow
Write-Host "`nStarting backup..." -ForegroundColor Green

try {
    # Copy entire project folder
    Copy-Item -Path $ProjectPath -Destination $BackupPath -Recurse -Force
    
    # Get backup size
    $BackupSize = (Get-ChildItem -Path $BackupPath -Recurse | Measure-Object -Property Length -Sum).Sum
    $BackupSizeMB = [math]::Round($BackupSize / 1MB, 2)
    
    Write-Host "`n[SUCCESS] Backup completed!" -ForegroundColor Green
    Write-Host "Backup size: $BackupSizeMB MB" -ForegroundColor Cyan
    Write-Host "Location: $BackupPath" -ForegroundColor Cyan
    
    # Create backup log
    $LogPath = Join-Path $BackupPath "backup_info.txt"
    $LogContent = @"
CardioVar Backup Information
============================
Backup Date: $(Get-Date -Format "yyyy-MM-dd HH:mm:ss")
Source: $ProjectPath
Backup Size: $BackupSizeMB MB
Status: SUCCESS

Files Backed Up:
$(Get-ChildItem -Path $BackupPath -Recurse -File | Select-Object -First 20 | ForEach-Object { "  - $($_.FullName.Replace($BackupPath, '.'))" })
... and more

Total Files: $((Get-ChildItem -Path $BackupPath -Recurse -File).Count)
Total Folders: $((Get-ChildItem -Path $BackupPath -Recurse -Directory).Count)
"@
    
    $LogContent | Out-File -FilePath $LogPath -Encoding UTF8
    
} catch {
    Write-Host "`n[ERROR] Backup failed: $_" -ForegroundColor Red
    exit 1
}

# Clean up old backups (keep only last $MaxBackups)
Write-Host "`nCleaning up old backups..." -ForegroundColor Yellow
$AllBackups = Get-ChildItem -Path $BackupRootPath -Directory | 
              Where-Object { $_.Name -like "Cardiovar_Backup_*" } | 
              Sort-Object Name -Descending

if ($AllBackups.Count -gt $MaxBackups) {
    $BackupsToDelete = $AllBackups | Select-Object -Skip $MaxBackups
    
    foreach ($OldBackup in $BackupsToDelete) {
        Write-Host "  Removing old backup: $($OldBackup.Name)" -ForegroundColor Gray
        Remove-Item -Path $OldBackup.FullName -Recurse -Force
    }
    
    Write-Host "Removed $($BackupsToDelete.Count) old backup(s)" -ForegroundColor Green
} else {
    Write-Host "No old backups to remove (keeping last $MaxBackups)" -ForegroundColor Green
}

# Summary
Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "Backup Summary" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Current backups: $($AllBackups.Count)" -ForegroundColor Cyan
Write-Host "Latest backup: $BackupPath" -ForegroundColor Green
Write-Host "`nBackup completed successfully!" -ForegroundColor Green
Write-Host "========================================`n" -ForegroundColor Cyan
