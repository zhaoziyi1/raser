Bootstrap: docker
From: rootproject/root:6.24.00-ubuntu20.04 

%environment 
    #GEANT4_INSTALL=/opt/geant4/install
    #export GEANT4_INSTALL 


%post
    apt update
    apt -y install --no-install-recommends fenics
    apt -y install tk-dev
    apt -y install python3-tk
    apt -y install vim emacs 

    #wget https://bootstrap.pypa.io/get-pip.py
    #python3 get-pip.py
    #python3 -m pip install matplotlib

    
    #echo "-----build geant4 environment -----" 
    apt -y install build-essential libboost-all-dev  
    apt -y install libgl1-mesa-dev libglu1-mesa-dev libxt-dev libxmu-dev libxi-dev zlib1g-dev libgl2ps-dev libexpat1-dev libxerces-c-dev
    
    #cd / 
    #mkdir -p /opt/geant4/src
    #mkdir -p /opt/geant4/build 
    #mkdir -p /opt/geant4/install
    #mkdir -p /opt/geant4/data 
    #wget https://geant4-data.web.cern.ch/releases/geant4.10.07.p02.tar.gz
    #tar xzf geant4.10.07.p02.tar.gz -C /opt/geant4/src 
    #rm /geant4.10.07.p02.tar.gz 
    #cd /opt/geant4/build 

    #echo 'export GEANT4_INSTALL=/opt/geant4/install'  >> $SINGULARITY_ENVIRONMENT 
    #echo ">>>> Print Environment >>>>"
    #printenv
    #echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>"      
    #cmake -DCMAKE_INSTALL_PREFIX=/opt/geant4/install ../src/geant4.10.07.p02
    #cmake -DCMAKE_INSTALL_PREFIX=/opt/geant4/install -DGEANT4_BUILD_TLS_MODEL=global-dynamic -DGEANT4_USE_PYTHON=ON  ../src/geant4.10.07.p02   
          #-DGEANT4_INSTALL_DATA=ON \
	  #-DGEANT4_INSTALL_DATADIR=/opt/geant4/data
	  #-DGEANT4_BUILD_TLS_MODEL=global-dynamic \
	  # -DGEANT4_USE_PYTHON=ON
    #make -j24 
    #make install 




