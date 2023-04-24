/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

// Provide help information
static char help[] = "Structure topology optimization solver\n";

// Include necessary libraries
#include "fvCFD.H"                              // Establishes basic CFD capabilities
#include "simpleControl.H"                      // Creates SIMPLE control
#include <math.h>                               // Enables standard math functions
//#include <fstream>
#include <iostream>
#include <iosfwd>
#include <stdio.h>
//#include "petsc.h"
//#include "petscvec.h"
#include "MMA/MMA.h"                            // Enables use of Method of Moving Asymptotes
//#include <MMA.h>
//#include <mpi.h>
#include <delta_gamma_filter.c>                 // Function that differentiates Heaviside filter by pseudo density

// Run main program
int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"                    // Sets root case
    #include "createTime.H"                     // Creates time
    #include "createMesh.H"                     // Creates mesh
    #include "createControl.H"                  // Enables controls
    #include "createFields.H"                   // Create parameters & fields
    #include "readMechanicalProperties.H" 
    #include "initializeSIMP.H"

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        #include "updateOptProps.h"
        #include "primalSolver.H"
        #include "costFunction.H"
        #include "sensitivity.H"
//        #include "updateOptProps.h"
        #include "output.h"
    }
    return 0;
}
