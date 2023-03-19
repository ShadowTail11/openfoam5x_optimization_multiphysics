//
// This function is a pre-run that calculated a fluid power dissipation reference value that is used to
// constrain the power dissipation of the final design by means of a power dissipation ratio.
//

#ifndef THERMFLUIDOPT_SOLVEPOWERLOSSREF_H
#define THERMFLUIDOPT_SOLVEPOWERLOSSREF_H

if (test_area == 1)
{
    // Iterate primal solver to calculate velocity U and kinematic pressure p under reference conditions
    for(i = 0 ; i < 400; i++)
    {
        Info << "\nSolve for U_test\n" << endl;

        tmp<fvVectorMatrix> tU_testEqn(
                fvm::div(phi_test, U_test)                // Divergent flow term
                + MRF.DDt(U_test)                    // Material derivative (not necessary for steady-state solver)
    //            + turbulence->divDevReff(U_test)   // Used for turbulent flow
                - fvm::laplacian(nu_eff, U_test)     // Used specifically for laminar flow
                + fvm::Sp(alpha_U_test, U_test)           // Creates a pseudo-porosity-based body force with a ramp function
                ==
                fvOptions(U_test)                    // Enables OpenFOAM options for velocity such as relaxation factors
        );

        fvVectorMatrix& U_testEqn = tU_testEqn.ref();

        U_testEqn.relax();

        fvOptions.constrain(U_testEqn);

        solve(U_testEqn == -fvc::grad(p_test));

        fvOptions.correct(U_test);

        //****************************************

        volScalarField rAU_test(1.0 / U_testEqn.A());
        volVectorField HbyA_test("HbyA_test", U_test);
        HbyA_test = rAU_test * U_testEqn.H();
        tU_testEqn.clear();
        surfaceScalarField phiHbyA_test("phiHbyA_test", fvc::flux(HbyA_test));
        adjustPhi(phiHbyA_test, U_test, p_test);

        tmp<volScalarField> rAtU_test(rAU_test);

        if(simple.consistent())
        {
            rAtU_test = 1.0 / (1.0 / rAU_test - U_testEqn.H1());
            phiHbyA_test +=
            fvc::interpolate(rAtU_test() - rAU_test) * fvc::snGrad(p_test) * mesh.magSf();
            HbyA_test -= (rAU_test - rAtU_test()) * fvc::grad(p_test);
        }

        tU_testEqn.clear();

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p_test, U_test, phiHbyA_test, rAtU_test(), MRF);

        // Non-orthogonal pressure corrector loop
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix p_testEqn(
                    fvm::laplacian(rAU_test, p_test) == fvc::div(phiHbyA_test)
                    );

            p_testEqn.setReference(p_testRefCell, p_testRefValue);
            p_testEqn.solve();

            if (simple.finalNonOrthogonalIter())
            {
                phi_test = phiHbyA_test - p_testEqn.flux();
            }
        }

        //#include "adjointContinuityErrs.H"

        // Explicitly relax pressure for adjoint momentum corrector
        p_test.relax();

        // Adjoint momentum corrector
        U_test = HbyA_test - rAU_test * fvc::grad(p_test);
        U_test.correctBoundaryConditions();
        fvOptions.correct(U_test);

        U_test.storePrevIter();
        p_test.storePrevIter();
        phi_test.storePrevIter();
    }

    Info << "Initial power_loss_ref = " << power_loss_ref << nl << endl;

    // Calculate the power loss reference value based on the test area
    power_loss_ref = 0;
    for(i = 0; i < nObjPatch; i++)
    {
    power_loss_ref = power_loss_ref - sum(
            phi_test.boundaryField()[conPatchList[i]]
            * (p_test.boundaryField()[conPatchList[i]] + 0.5 * magSqr(U_test.boundaryField()[conPatchList[i]])));
    }

    if(geo_dim == 2)
    {
        power_loss_ref = power_loss_ref / len;
    }

    Info << "Updated power_loss_ref = " << power_loss_ref << nl << endl;

    // Turn off test
    test_area = 0;
}

#endif //THERMFLUIDOPT_SOLVEPOWERLOSSREF_H
