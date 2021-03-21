# Setup the env for Mac

ip=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}')
xhost + $ip




