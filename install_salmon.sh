#! bin/bash

# This only works up to the v1.10.0. 

wget -q https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
tar zxf salmon-1.10.0_linux_x86_64.tar.gz
export PATH=$PATH:$PWD/salmon-latest_linux_x86_64/bin
salmon
