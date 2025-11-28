#!/bin/bash

# Start FastAPI in the background
echo "Starting FastAPI backend..."
uvicorn api:app --host 0.0.0.0 --port 8000 &

# Wait for FastAPI to start
echo "Waiting for API to be ready..."
timeout 30s bash -c 'until curl -s http://localhost:8000/status > /dev/null; do sleep 1; done'

# Start Streamlit
echo "Starting Streamlit frontend..."
streamlit run dashboard.py --server.port 7860 --server.address 0.0.0.0
