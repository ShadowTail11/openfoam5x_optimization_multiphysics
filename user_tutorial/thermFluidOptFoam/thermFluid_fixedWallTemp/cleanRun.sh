#!/bin/bash

find . -type d -regex './[1-9][0-9]*' -exec rm -rf {} +
blockMesh
thermFluidOptFoam
