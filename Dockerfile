FROM rootproject/root:6.24.00-ubuntu20.04
#FROM rootproject/root

# Run the following commands as super user (root):
USER root

# Set directory to /tmp for build processes
WORKDIR /tmp
# WORKDIR /var/lib/apt/lists/
# WORKDIR /var/tmp/
# Set directory to /tmp for build processes

# 2021.5.14 jia ning
# Set locale environment
# ENV  LC_ALL        C.UTF-8  
# ENV  LANG           C.UTF-8  
# ENV  LANGUAGE  C.UTF-8


# # Install pip & pip modules
# RUN wget https://bootstrap.pypa.io/get-pip.py  && \
#     python3 get-pip.py  && \
#     python3 -m pip install numpy  && \
#     python3 -m pip install sympy  && \
#     python3 -m pip install ply  && \
#     rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# # Install Eigen3
# RUN wget https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz  && \
#     mkdir eigen  && \
#     tar -xzf eigen-3.3.7.tar.gz -C eigen --strip-components 1  && \
#     cd eigen  &&  mkdir build  && cd build  && \
#     cmake -DCMAKE_INSTALL_PREFIX=/usr  ../   && \
#     make install   && \
#     rm -rf  /tmp/*

# # Install libraries ahead of the CERN ROOT installation
# RUN apt-get update && \
#     apt-get -y install \
#             libboost-all-dev   &&\
#     apt-get clean  && \
#     rm -rf  /tmp/* 

# ### For installation of the dolfin ###

# # Install Swig
# RUN wget https://sourceforge.net/projects/swig/files/swig/swig-3.0.12/swig-3.0.12.tar.gz  && \
#     mkdir swig  &&  \
#     tar -xzf swig-3.0.12.tar.gz  -C swig --strip-components 1  && \
#     cd swig  && \
#     ./configure  && \
#     make  &&  make install  && \
#     rm -rf  /tmp/*

# # Install libraries
# RUN apt-get update && \
#     apt-get -y install mpich  libxt-dev liblapack-dev libopenblas-dev && \
#     apt-get clean  && \    
#     rm -rf  /tmp/*


# # Install UMFPACK
# RUN wget https://people.engr.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.4.tar.gz   && \
#     mkdir suitesparse  && \
#     tar -xzf SuiteSparse-4.4.4.tar.gz  -C suitesparse --strip-components 1   && \
#     cd suitesparse   &&  \
#     make   &&  make install  && \
#     rm -rf  /tmp/*


# # Install FEniCS libraries
# RUN wget https://bitbucket.org/fenics-project/fiat/downloads/fiat-1.5.0.tar.gz   && \
#     mkdir fiat  && \
#     tar -xzf fiat-1.5.0.tar.gz  -C fiat --strip-components 1  && \
#     cd fiat  &&  python3 -m pip install .   && \
#     rm -rf  /tmp/*

# RUN wget https://bitbucket.org/fenics-project/dijitso/downloads/dijitso-2016.1.0.tar.gz  && \
#     mkdir dijitso  && \
#     tar -xzf dijitso-2016.1.0.tar.gz  -C dijitso --strip-components 1  && \
#     cd dijitso  &&  python3 -m pip install .   && \
#     rm -rf  /tmp/*

# RUN wget https://bitbucket.org/fenics-project/ufl/downloads/ufl-1.5.0.tar.gz   && \
#     mkdir ufl  && \
#     tar -xzf ufl-1.5.0.tar.gz  -C ufl --strip-components 1   && \
#     cd ufl  &&  python3 -m pip install .   && \
#     rm -rf  /tmp/*

# RUN wget https://bitbucket.org/fenics-project/ffc/downloads/ffc-1.5.0.tar.gz   && \
#     mkdir ffc  && \
#     tar -xzf ffc-1.5.0.tar.gz  -C ffc --strip-components 1   && \
#     cd ffc  &&  python3 -m pip install .   && \
#     rm -rf  /tmp/*

# RUN wget https://bitbucket.org/fenics-project/dolfin/downloads/dolfin-1.5.0.tar.gz   && \
#     mkdir dolfin   && \
#     tar -xzf dolfin-1.5.0.tar.gz  -C dolfin --strip-components 1   && \
#     cd dolfin  &&  \
#     mkdir build  && cd build   && \
#     cmake -DPYTHON_EXECUTABLE=/usr/bin/python3  ../     && \
#     make install   && \
#     rm -rf  /tmp/*

# Set locale environment
RUN apt update 
RUN apt -y install --no-install-recommends fenics
RUN apt -y install tk-dev
RUN apt -y install python3-tk
RUN apt -y install vim

RUN wget https://bootstrap.pypa.io/get-pip.py
RUN python3 get-pip.py
RUN python3 -m pip install matplotlib

# Create a user that does not have root privileges 
ARG username=physicist
RUN useradd --create-home --home-dir /home/${username} ${username}
ENV HOME /home/${username}

# Switch to our newly created user
USER ${username}

# Our working directory will be in our home directory where we have permissions
WORKDIR /home/${username}





