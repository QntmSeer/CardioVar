# Use Python 3.9 as base image
FROM python:3.9-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    software-properties-common \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements file
COPY requirements.txt .

# Install Python dependencies
# Note: enformer-pytorch and torch might be heavy, so we install them carefully
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Create start script
COPY start.sh .
RUN chmod +x start.sh

# Expose port 7860 for Hugging Face Spaces
EXPOSE 7860

# Run the start script
CMD ["./start.sh"]
