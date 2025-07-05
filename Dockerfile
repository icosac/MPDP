FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
    python3 \
    build-essential \
    cmake 



