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
static char help[] = "Thermal-fluid-structure topology optimization solver\n";

// Include necessary libraries
#include "fvCFD.H"                              // Establishes basic CFD capabilities
#include "singlePhaseTransportModel.H"          // Establishes single phase transport model
#include "turbulentTransportModel.H"            // Establishes turbulent model
#include "simpleControl.H"                      // Creates SIMPLE control
#include "fvOptions.H"                          // Establishes fvOptions
#include <math.h>                               // Enables standard math functions
#include "MMA/MMA.h"                            // Enables use of Method of Moving Asymptotes
#include <delta_gamma_filter.c>                 // Function that differentiates Heaviside filter by pseudo density

// Begin main program
int main(int argc, char *argv[])
{
    #include "setRootCase.H"                    // Sets root case
    #include "createTime.H"                     // Creates time
    #include "createMesh.H"                     // Creates mesh
    #include "createControl.H"                  // Enables controls
    #include "createFvOptions.H"                // Creates fvOptions for relaxation factors
    #include "createFields.H"                   // Create parameters & fields
    #include "initContinuityErrs.H"             // Enables continuity error tracking
    #include "initializeParameters.H"           // Initialize other parameters

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        #include "updateOptProperties.H"        // Update optimization-related properties

        #include "primalFlowSolver.H"           // Run primal solver for U and p
        #include "primalThermalSolver.H"        // Run primal solver for T
        #include "LinearElasticity.H"           // Run primal solver for D
        #include "updatePrimalProperties.h"     // Update primal properties that are based on temperature

        #include "adjointThermalSolverT.H"      // Run adjoint solver for T_adj
        #include "adjointThermalSolverU.H"      // Run adjoint solver for U_adj_T
        #include "adjointFlowSolverU.H"         // Run adjoint solver for U_adj_U

        #include "costFunction.H"               // Calculate the cost function & convergence properties
        #include "sensitivity.H"                // Calculate sensitivity and update pseudo density accordingly
        #include "output.h"                     // Check convergence and output monitoring/state variables
    }
    Info << "End\n" << endl;
    return 0;
}
