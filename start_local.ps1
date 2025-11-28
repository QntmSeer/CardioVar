# Start FastAPI in the background
Write-Host "Starting FastAPI backend..."
$apiProcess = Start-Process -FilePath "uvicorn" -ArgumentList "api:app", "--host", "127.0.0.1", "--port", "8000" -PassThru -NoNewWindow

# Wait for API to be ready
Write-Host "Waiting for API to be ready..."
Start-Sleep -Seconds 5

# Start Streamlit
Write-Host "Starting Streamlit frontend..."
$streamlitProcess = Start-Process -FilePath "streamlit" -ArgumentList "run", "dashboard.py", "--server.port", "8501" -PassThru -NoNewWindow

Write-Host "Application started!"
Write-Host "API: http://127.0.0.1:8000"
Write-Host "Dashboard: http://localhost:8501"
Write-Host "Press any key to stop..."
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")

# Stop processes
Stop-Process -Id $apiProcess.Id -Force
Stop-Process -Id $streamlitProcess.Id -Force
Write-Host "Stopped."
