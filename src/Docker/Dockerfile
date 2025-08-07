# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

FROM python:3.9-slim

WORKDIR /app

COPY src/Docker/requirements.txt requirements.txt
RUN pip install --trusted-host pypi.org --trusted-host pypi.python.org --trusted-host=files.pythonhosted.org --upgrade pip && \
    pip install --trusted-host pypi.org --trusted-host pypi.python.org --trusted-host=files.pythonhosted.org --prefer-binary -r requirements.txt

RUN	mkdir -p /usr/share/man/man1 && \
	apt-get update && \
	apt-get install -y \
		openjdk-17-jre \
		unzip \
        nano \
        vim

COPY DockerImport/libOPERA2.9_Py.tar.gz libOPERA2.9_Py.tar.gz
ENV OPERA_EXE_PATH=/usr/local/bin/OPERA/application/run_OPERA.sh
ENV MATLAB_RUNTIME_PATH=/usr/local/MATLAB/MATLAB_Runtime/v912
ENV LD_LIBRARY_PATH=.:/usr/local/MATLAB/MATLAB_Runtime/v912/runtime/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v912/bin/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v912/sys/os/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v912/sys/opengl/lib/glnxa64
ENV XAPPLRESDIR=/usr/local/MATLAB/MATLAB_Runtime/v912/X11/app-defaults

RUN	tar xzf libOPERA2.9_Py.tar.gz
RUN	rm libOPERA2.9_Py.tar.gz
RUN cd libOPERA2_Py/ && \
	./OPERA2.9_Py_mcr.install -mode silent -agreeToLicense yes
    
# Temporarily downgrade setuptools to allow OPERA's invalid version string
RUN pip install --trusted-host pypi.org --trusted-host pypi.python.org --trusted-host=files.pythonhosted.org "setuptools<60.0.0"
RUN cd /usr/local/bin/OPERA/application/ && \
	python setup.py install

RUN rm /usr/local/MATLAB/MATLAB_Runtime/v912/java/jarext/log4j.*

ENV FLASK_APP=App/App.py
EXPOSE 5005

COPY DockerImport/volumes/etcd /etcd_data
COPY DockerImport/volumes/milvus /milvus_data
COPY DockerImport/volumes/minio /minio_data

# Define the volume mount points
VOLUME ["/etcd_data", "/milvus_data", "/minio_data"]

# Copy application files
COPY src/Comparison src/Comparison
COPY src/App src/App
COPY src/main.py src/main.py
COPY src/utils src/utils
COPY src/OPERA /app/OPERA
COPY src/Docker/startup.sh /app/startup.sh
COPY tests tests

# Clean up any existing files and create all required directories
RUN rm -rf src/Comparison/LocalIO/* || true
RUN rm -rf src/App/static/LocalIO/* || true
RUN mkdir -p output logs \
    && mkdir -p src/Comparison/LocalIO \
    && mkdir -p src/App/static/LocalIO \
    && mkdir -p src/App/static/single_jobs \
    && mkdir -p src/Comparison/LocalIO/single_jobs \
    && touch src/Comparison/LocalIO/data.csv \
    && touch src/App/static/LocalIO/graph.png \
    && touch src/App/static/LocalIO/Property-Predictions.csv \
    && touch src/App/static/LocalIO/results.csv \
    && chmod -R 777 src/Comparison/LocalIO \
    && chmod -R 777 src/App/static/LocalIO \
    && chmod -R 777 output \
    && chmod -R 777 logs \
    && chmod +x /app/startup.sh \
    && chmod +x tests/run_tests.sh

# Create a placeholder no-graph image
RUN echo '<svg width="800" height="400" xmlns="http://www.w3.org/2000/svg"><rect width="100%" height="100%" fill="#f0f0f0"/><text x="50%" y="50%" font-family="Arial" font-size="24" text-anchor="middle" fill="#666">No graph data available</text></svg>' > src/App/static/assets/no-graph.svg \
    && echo '<svg width="800" height="400" xmlns="http://www.w3.org/2000/svg"><rect width="100%" height="100%" fill="#f0f0f0"/><text x="50%" y="50%" font-family="Arial" font-size="24" text-anchor="middle" fill="#666">No graph data available</text></svg>' > src/App/static/LocalIO/graph.png

LABEL org.opencontainers.image.source="https://github.com/sandialabs/chemical-recommender-system"

# Use the startup script instead of direct command
CMD ["/app/startup.sh"]