#!/bin/bash

// This includes all functions required even if not used
// git status
git add --all
git commit -m "Expand zone_test to include temperature and pressure drop in thermFluidOptFoam"
git push -u origin main
