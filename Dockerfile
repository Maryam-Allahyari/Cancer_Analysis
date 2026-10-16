FROM python:3.13-slim

# libgomp1 is needed by XGBoost, procps lets Nextflow collect task stats
RUN apt-get update \
    && apt-get install -y --no-install-recommends libgomp1 procps \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# matplotlib needs a writable folder for its cache when running as a non-root user
ENV MPLCONFIGDIR=/tmp/matplotlib

LABEL description="Python environment for the cancer origin pipeline"
