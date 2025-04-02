#!/bin/bash

PACKAGE_FILE=multiDoE_0.9.3.tar.gz

if [[ ! -f $PACKAGE_FILE ]]; then
  echo "File not found: $PACKAGE_FILE. Run 'R CMD build .' first. Then check the file name."
  exit 1
fi

# CRAN suggests to check the tar.gz and not the folder
R CMD check --as-cran $PACKAGE_FILE
