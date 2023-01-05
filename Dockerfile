FROM condaforge/mambaforge

ENV DEBIAN_FRONTEND=noninteractive

# Install ssh (missing dependency to run conda envs)
RUN apt-get update && \
    apt-get install -y ssh build-essential libgl1-mesa-glx mesa-utils

# Upgrade mamba
RUN mamba upgrade -y mamba

# Copy environment and requirements files into docker env
COPY environment.yml .

# Update environment file with new environment name
RUN mamba env update --file environment.yml --name dockerenv

RUN echo "source activate dockerenv" > ~/.bashrc

# Fix issues with VMTK
RUN sed -i "s/file[[:space:]]=[[:space:]]open(self.OutputFileName,'r')/file = open(self.OutputFileName,'rb')/g" /opt/conda/envs/dockerenv/lib/python3.10/site-packages/vmtk/vmtkmeshwriter.py
RUN sed -i "s/gzfile[[:space:]]=[[:space:]]gzip.open(self.OutputFileName,'w')/gzfile = gzip.open(self.OutputFileName,'wb')/g" /opt/conda/envs/dockerenv/lib/python3.10/site-packages/vmtk/vmtkmeshwriter.py


SHELL ["mamba", "run", "-n", "dockerenv", "/bin/bash", "-c"]
ENV MESA_LOADER_DRIVER_OVERRIDE=""
