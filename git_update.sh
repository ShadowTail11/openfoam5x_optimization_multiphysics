#!/bin/bash

// This includes all functions required even if not used
// git status
git add --all
git commit -m "Add fluidOptFoam solver for optimizing only pressure drop (via power loss)"
git push -u origin main
