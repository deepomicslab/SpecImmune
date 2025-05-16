# Use Ubuntu as the base image
FROM ubuntu:latest

# Install necessary packages and dependencies
RUN apt-get update && apt-get install -y \
    wget \
    software-properties-common \
    curl \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /SpecImmune

# Define SHELL
SHELL ["/bin/bash", "-c"] 

# Copy code and scripts into the container
COPY . /SpecImmune

# Download and install Miniconda
RUN wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh

# Add Conda to the PATH
ENV PATH="/opt/conda/bin:$PATH"

# Update Conda
RUN conda update -n base -c defaults conda -y  \
    && conda clean --all \
    && conda env create -p /opt/conda/envs/SpecImmune -f /SpecImmune/environment.yml \ 
    && source /opt/conda/envs/SpecImmune/bin/activate \    
    && echo "source /opt/conda/bin/activate SpecImmune" >> ~/.bashrc

# Set the ENTRYPOINT
ENTRYPOINT ["/bin/bash", "-c", "echo \"Passed arguments: $@\"; source /opt/conda/bin/activate SpecImmune && exec python /SpecImmune/scripts/main.py \"$@\"", "--"]
