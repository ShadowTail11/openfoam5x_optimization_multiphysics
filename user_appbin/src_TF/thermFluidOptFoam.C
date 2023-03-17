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

// Include necessary libraries

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include <math.h>
#include <fstream>
#include <iostream>
#include <iosfwd>
#include <stdio.h>
#include "petsc.h"
#include "petscvec.h"
#include <MMA.h>
#include <mpi.h>

// Initialize classes and cells

// Function sets the value of a zone of cells
template<class Type>
void setCells(
    
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells,
    double value
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = value;
    }
}

double fun(double gamma[], double del, double eta, int allcells);
static char help[] = "topology optimization \n";

// Begin main program

int main(int argc, char *argv[])
{
    PetscInitialize(&argc,&argv,PETSC_NULL,help);
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"       // Read & create parameters & fields
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "SIMP_initialize.H"    // Initialize other parameters

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        #include "primal_flow_solver.H"
        #include "primal_thermal_solver.H"
        #include "update_primal_properties.H"

        #include "adjoint_thermal_Tb_solver.H"
        #include "adjoint_flow_T.H"
        #include "adjoint_flow_U.H"

        #include "costfunction.H"
        #include "sensitivity.H"
        
        if(runTime.writeTime())
        {
            gamma.write();
            T.write();
            U.write();
            p.write();
            nu_eff.write();
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

        // Exit if convergence has been achieved
        if (dmeanT_RMS < converge_tol && opt > 20) {
            gamma.write();
            T.write();
            U.write();
            p.write();
            nu_eff.write();
            Info << "Convergence criterion (<" << converge_tol * 100 << "%) has been met after "
                 << opt - 1 << " iterations!" << endl << "Program ending!" << endl;

            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
            
            break;
        }
    }

    delete mma;
    PetscFinalize();

    Info << "End\n" << endl;
    
    return 0;
}


//*********************************************************//function

double fun(double gamma[], double del, double eta, int allcells)
{
    int i;
    double z = 0;
    double *fg = new double[allcells];

    for(i = 0; i < allcells; i++)
    {
        if(gamma[i] <= eta)
        {
            fg[i] = eta * (Foam::exp(-del * (1 - gamma[i] / eta)) - (1 - gamma[i] / eta) * Foam::exp(-del));
        }
        else
        {
            fg[i] = eta + (1 - eta) * (1 - Foam::exp(-del * (gamma[i] - eta) / (1 - eta)) + (gamma[i] - eta) * Foam::exp(-del) / (1 - eta));
        }
    }

    for(i = 0; i < allcells; i++)
    {
        z = z + gamma[i] - fg[i];
    }

    delete fg;

    return {z};
}
