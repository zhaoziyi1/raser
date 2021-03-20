#
# A Dockerfile for TRACS v1.1 image
#

# Use Debian 9 ( "stretch" ) image as a base 
FROM debian:stretch

# Execute by "root" user
USER root

# Set directory to /tmp for build processes
WORKDIR /tmp

# Set password for the root user
RUN echo "root:tracs_dev" | chpasswd

# Install necessary commands/libraries ( # some of them are already included in build-essential ... )
RUN apt-get update && \
    apt-get -y install  \
            build-essential  \
            cmake \
            doxygen \
            emacs \
            evince \
            gimp \
            git \
            less \
            locate \
            locales \
            sudo \
            vim \
            wget \
            xterm && \
    apt-get clean  && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# Set locale environment
ENV  LC_ALL        C.UTF-8  
ENV  LANG           C.UTF-8  
ENV  LANGUAGE  C.UTF-8
    

# Install python3 environment
RUN apt-get update && \
    apt-get -y install \
            python3 \
            python3-dev \
            python3-pip   && \
#    ln -s /usr/bin/python3 /usr/bin/python  && \   
    apt-get clean  && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# Install pip & pip modules
RUN wget https://bootstrap.pypa.io/get-pip.py  && \
    python3 get-pip.py  && \
    python3 -m pip install numpy  && \
    python3 -m pip install sympy  && \
    python3 -m pip install ply  && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*



# Install libraries ahead of the CERN ROOT installation
RUN apt-get update && \
    apt-get -y install \
            libx11-dev libxpm-dev libxft-dev libxext-dev  \
            libtiff5-dev libgif-dev libgsl-dev  \
            libkrb5-dev libxml2-dev libssl-dev \
            default-libmysqlclient-dev libpq-dev libqt4-opengl-dev libgl2ps-dev libpcre-ocaml-dev \
            libgraphviz-dev libdpm-dev unixodbc-dev libsqlite3-dev libfftw3-dev libcfitsio-dev  \
            dcap-dev libldap2-dev libavahi-compat-libdnssd-dev  \
            libboost-all-dev  && \
    apt-get clean  && \            
    rm -rf  /tmp/* 
            

# Install CERN ROOT libraries
RUN wget https://root.cern/download/root_v6.12.04.source.tar.gz  && \
    mkdir root && \
    tar -xzf root_v6.12.04.source.tar.gz -C root --strip-components 1  && \
    cd root  &&  mkdir obj  &&  cd obj  && \
    cmake  -DCMAKE_INSTALL_PREFIX="/usr/local/root"  \
           -DPYTHON_EXECUTABLE="/usr/bin/python3" \
           -DPYTHON_LIBRARIES="/usr/lib/python3.5/config-3.5m-x86_64-linux-gnu/libpython3.5.so" \
           -Dbuiltin_cfitsio="ON" \
           -Dbuiltin_fftw3="ON"  \
           -Dbuiltin_gsl="ON"  \
           -Dbuiltin_xrootd="ON"  \
           -Dgnuinstall="ON"  \
           -Dmathmore="ON"  \
           -Dminuit2="ON"  \
           -Droofit="ON"  \
           -Dtable="ON"  \
           -Dx11="ON"  \
           -DCMAKE_INSTALL_INCLUDEDIR:PATH="include"  \
           -DCMAKE_INSTALL_LIBDIR:PATH="lib"  \
           ../               && \
    make   && \
    make install  && \
    rm -rf  /tmp/*


# Install QT4
RUN apt-get update && \
    apt-get -y install qt4-default  && \
    apt-get clean  && \    
    rm -rf  /tmp/*

# Install Eigen3
RUN wget https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz  && \
    mkdir eigen  && \
    tar -xzf eigen-3.3.7.tar.gz -C eigen --strip-components 1  && \
    cd eigen  &&  mkdir build  && cd build  && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr  ../   && \
    make install   && \
    rm -rf  /tmp/*


### For installation of the dolfin ###

# Install Swig
RUN wget https://sourceforge.net/projects/swig/files/swig/swig-3.0.12/swig-3.0.12.tar.gz  && \
    mkdir swig  &&  \
    tar -xzf swig-3.0.12.tar.gz  -C swig --strip-components 1  && \
    cd swig  && \
    ./configure  && \
    make  &&  make install  && \
    rm -rf  /tmp/*

# Install libraries
RUN apt-get update && \
    apt-get -y install mpich  libxt-dev liblapack-dev libopenblas-dev && \
    apt-get clean  && \    
    rm -rf  /tmp/*

# Install VTK
RUN wget https://www.vtk.org/files/release/7.1/VTK-7.1.1.tar.gz  && \
    mkdir vtk  && \
    tar -xzf VTK-7.1.1.tar.gz  -C vtk --strip-components 1  && \
    cd vtk  && \
    mkdir build  &&  cd build  && \
    cmake ../ -DPYTHON_EXECUTABLE=/usr/bin/python3 -DVTK_Group_Qt:BOOL=ON -DCMAKE_INSTALL_PREFIX=/usr/local/vtk/  && \
    make  &&  make install   && \
    rm -rf  /tmp/*

# Install UMFPACK
RUN wget https://people.engr.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.4.tar.gz   && \
    mkdir suitesparse  && \
    tar -xzf SuiteSparse-4.4.4.tar.gz  -C suitesparse --strip-components 1   && \
    cd suitesparse   &&  \
    make   &&  make install  && \
    rm -rf  /tmp/*


# Install FEniCS libraries
RUN wget https://bitbucket.org/fenics-project/fiat/downloads/fiat-1.5.0.tar.gz   && \
    mkdir fiat  && \
    tar -xzf fiat-1.5.0.tar.gz  -C fiat --strip-components 1  && \
    cd fiat  &&  python3 -m pip install .   && \
    rm -rf  /tmp/*

RUN wget https://bitbucket.org/fenics-project/dijitso/downloads/dijitso-2016.1.0.tar.gz  && \
    mkdir dijitso  && \
    tar -xzf dijitso-2016.1.0.tar.gz  -C dijitso --strip-components 1  && \
    cd dijitso  &&  python3 -m pip install .   && \
    rm -rf  /tmp/*

RUN wget https://bitbucket.org/fenics-project/ufl/downloads/ufl-1.5.0.tar.gz   && \
    mkdir ufl  && \
    tar -xzf ufl-1.5.0.tar.gz  -C ufl --strip-components 1   && \
    cd ufl  &&  python3 -m pip install .   && \
    rm -rf  /tmp/*

RUN wget https://bitbucket.org/fenics-project/ffc/downloads/ffc-1.5.0.tar.gz   && \
    mkdir ffc  && \
    tar -xzf ffc-1.5.0.tar.gz  -C ffc --strip-components 1   && \
    cd ffc  &&  python3 -m pip install .   && \
    rm -rf  /tmp/*

RUN wget https://bitbucket.org/fenics-project/dolfin/downloads/dolfin-1.5.0.tar.gz   && \
    mkdir dolfin   && \
    tar -xzf dolfin-1.5.0.tar.gz  -C dolfin --strip-components 1   && \
    cd dolfin  &&  \
    mkdir build  && cd build   && \
    cmake -DPYTHON_EXECUTABLE=/usr/bin/python3  ../     && \
    make install   && \
    rm -rf  /tmp/*

### End the install process for dolfin library ###


# User add.
# username : tracs
# pass : tracs
# group : -> sudo ( to use sudo command )
RUN useradd -m -s /bin/bash -G sudo tracs &&  \
    echo "tracs:tracs" | chpasswd && \
    echo "tracs ALL=(ALL) NOPASSWD: ALL \n" >> /etc/sudoers


# Visit to home directory
WORKDIR /home/tracs

# Switch to "tracs" user
USER tracs

RUN mkdir /home/tracs/work