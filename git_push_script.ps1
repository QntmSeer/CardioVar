# Git Push Script
Write-Host "Adding all changes..."
git add .

Write-Host "Committing changes..."
git commit -m "feat: Integrate Enformer Deep Learning model and validate against clinical variants"

Write-Host "Pushing to remote..."
git push origin main

Write-Host "Done!"
