FROM rootproject/root

# Run the following commands as super user (root):
USER root

# Create a user that does not have root privileges 
ARG username=physicist
RUN useradd --create-home --home-dir /home/${username} ${username}
ENV HOME /home/${username}

# Switch to our newly created user
USER ${username}

# Our working directory will be in our home directory where we have permissions
WORKDIR /home/${username}


