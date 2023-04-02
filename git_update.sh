#!/bin/bash

// This includes all functions required even if not used
// git status
git add --all
git commit -m "Change heat generation input to W/m^3 in thermFluidOptFoam"
git push -u origin main
