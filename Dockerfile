FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install essential dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libboost-all-dev \
    libcgal-dev \
    libeigen3-dev \
    libomp-dev \
    libceres-dev \
    libmatio-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install yaml-cpp 
WORKDIR /tmp
RUN git clone https://github.com/jbeder/yaml-cpp.git \
    && cd yaml-cpp \
    && mkdir build \
    && cd build \
    && cmake -DYAML_BUILD_SHARED_LIBS=ON .. \
    && make -j$(nproc) \
    && make install

# Install MatioCpp
WORKDIR /tmp
RUN git clone https://github.com/ami-iit/matio-cpp.git \
    && cd matio-cpp \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make -j$(nproc) \
    && make install

# Install Catch2
WORKDIR /tmp
RUN git clone https://github.com/catchorg/Catch2.git \
    && cd Catch2 \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make -j$(nproc) \
    && make install

# Update library path
RUN ldconfig

# Create app directory and set it as the working directory
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app/

# Create build directory and compile the project
RUN rm -rf /app/build/Release && mkdir -p /app/build/Release && \
    cd /app/build/Release && \
    cmake ../.. && \
    make -j12

# Set the default command
CMD ["/bin/bash"]