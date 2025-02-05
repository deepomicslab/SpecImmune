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

# Copy code and scripts into the container
COPY . /SpecImmune

# Download and install Miniconda
RUN wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh

# Add Conda to the PATH
ENV PATH="/opt/conda/bin:$PATH"

# Update Conda
RUN conda update -n base -c defaults conda -y

# Install mamba
#RUN conda install -n base -c conda-forge mamba -y

# Create a Conda environment
RUN conda env create -n SpecImmune -f /SpecImmune/environment.yml

# Initialize Conda and activate the environment
RUN echo "source /opt/conda/bin/activate SpecImmune" >> ~/.bashrc

# Delete unnecessary directories. You can also directly ignore them using the .dockerignore file.
#RUN rm -rf /SpecImmune/test \
#    /SpecImmune/pharmvar-6.2 \
#    /SpecImmune/simulation \
#    /SpecImmune/evaluation \
#    /SpecImmune/SpecImmune

# Set the ENTRYPOINT
ENTRYPOINT ["/bin/bash", "-c", "echo \"Passed arguments: $@\"; source /opt/conda/bin/activate SpecImmune && exec python /SpecImmune/scripts/main.py \"$@\"", "--"]