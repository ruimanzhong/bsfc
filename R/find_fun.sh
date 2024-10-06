#!/bin/bash

# Get a list of all cccd functions from the R environment
functions=$(Rscript -e 'library(spdep); cat(ls("package:spdep"), sep=" ")')
# Search for these functions in your current directory
for func in $functions; do
  echo "Searching for $func..."
  grep -r "$func" .
done

