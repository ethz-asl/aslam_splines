#!/bin/bash -e
if [[ $(uname) == "Linux" ]]; then
  sudo apt-get install -y python-scipy
else
  echo "Platform $(uname) is not supported!"
fi
