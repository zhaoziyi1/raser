FROM rootproject/root:6.24.00-ubuntu20.04
#FROM rootproject/root

# Run the following commands as super user (root):
USER root

# Set directory to /tmp for build processes
WORKDIR /tmp

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





