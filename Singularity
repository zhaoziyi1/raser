Bootstrap: docker
From: rootproject/root:6.24.00-ubuntu20.04 

%post
    apt update
    apt -y install --no-install-recommends fenics
    apt -y install tk-dev
    apt -y install python3-tk
    apt -y install vim 

    wget https://bootstrap.pypa.io/get-pip.py
    python3 get-pip.py
    python3 -m pip install matplotlib

    
    echo "-----build geant4-----" 
    apt -y install build-essential libboost-all-dev  
    apt -y install libgl1-mesa-dev libglu1-mesa-dev libxt-dev libxmu-dev libxi-dev zlib1g-dev libgl2ps-dev libexpat1-dev libxerces-c-dev
    
    cd / 
    mkdir -p /opt/geant4/src
    mkdir -p /opt/geant4/build 
    mkdir -p /opt/geant4/install
    mkdir -p /opt/geant4/data 
    wget https://geant4-data.web.cern.ch/releases/geant4.10.07.p02.tar.gz
    tar xzf geant4.10.07.p02.tar.gz -C /opt/geant4/src 
    rm /geant4.10.07.p02.tar.gz 
    cd /opt/geant4/build 
    cmake -DCMAKE_INSTALL_PREFIX=/opt/geant4/install \
          -DGEANT4_INSTALL_DATA=ON \
	  -DGEANT4_INSTALL_DATADIR=/opt/geant4/data ../src/geant4.10.07.p02  
    make -j24 
    make install 

    

