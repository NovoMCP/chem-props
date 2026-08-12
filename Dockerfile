FROM python:3.11-slim

WORKDIR /app

# Install system dependencies for RDKit
RUN apt-get update && apt-get install -y \
    curl \
    gcc \
    g++ \
    libxrender1 \
    libx11-6 \
    libcairo2 \
    libfreetype6 \
    libboost-system-dev \
    libboost-thread-dev \
    libboost-serialization-dev \
    libboost-python-dev \
    libboost-regex-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN pip install --upgrade pip

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY main.py .
COPY qed_v2.py .

# Expose service port
EXPOSE 8003

# Run the service
CMD ["python", "main.py"]