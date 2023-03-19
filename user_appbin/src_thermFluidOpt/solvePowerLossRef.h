//
// This function is a pre-run that calculated a fluid power dissipation reference value that is used to
// constrain the power dissipation of the final design by means of a power dissipation ratio.
//

#ifndef THERMFLUIDOPT_SOLVEPOWERLOSSREF_H
#define THERMFLUIDOPT_SOLVEPOWERLOSSREF_H

if (test_area == 1)
{
    Info << "Initial power_loss_ref = " << power_loss_ref << nl << endl;

    #include "primalFlowSolver.H"           // Run primal solver for U and p

    // Calculate the power loss reference value based on the test area
    power_loss_ref = 0;
    for(i = 0; i < nObjPatch; i++)
    {
    power_loss_ref = power_loss_ref - sum(
            phi.boundaryField()[conPatchList[i]]
            * (p.boundaryField()[conPatchList[i]] + 0.5 * magSqr(U.boundaryField()[conPatchList[i]])));
    }

    if(geo_dim == 2)
    {
        power_loss_ref = power_loss_ref / len;
    }

    Info << "Updated power_loss_ref = " << power_loss_ref << nl << endl;

    // Reset domain for the beginning of the optimization run
    forAll(gamma, i){gamma[i] = set_vol_frac;}
    setCells(gamma, cells1, 0);
    if(solid_area == 1){setCells(gamma, cells_solid, 0);}
    if(fluid_area == 1){setCells(gamma, cells_fluid, 1);}
}

#endif //THERMFLUIDOPT_SOLVEPOWERLOSSREF_H
