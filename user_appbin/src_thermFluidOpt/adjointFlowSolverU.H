for(i = 0 ; i < 40; i++)
{
    Info << "\nSolve for adjoint flow-U\n" << endl;

    volVectorField adjointTransposeConvectionc((fvc::grad(U_adj_U) & U));

    tmp<fvVectorMatrix> tUcEqn(

        fvm::div(-phi, U_adj_U)
        - adjointTransposeConvectionc
//        + turbulence->divDevReff(U_adj_U)
        - fvm::laplacian(nu_eff, U_adj_U)
        + fvm::Sp(alpha_U, U_adj_U)
        ==
        fvOptions(U_adj_U)
    );

    fvVectorMatrix& UcEqn = tUcEqn.ref();

    UcEqn.relax();

    fvOptions.constrain(UcEqn);

    solve(UcEqn == -fvc::grad(p_adj_U));

    fvOptions.correct(U_adj_U);

    //****************************************

    volScalarField rAUc(1.0 / UcEqn.A());
    volVectorField HbyAc("HbyAc", U_adj_U);
    HbyAc = rAUc * UcEqn.H();
    tUcEqn.clear();
    surfaceScalarField phiHbyAc("phiHbyAc", fvc::flux(HbyAc));
    adjustPhi(phiHbyAc, U_adj_U, p_adj_U);

    tmp<volScalarField> rAtUc(rAUc);

    if(simple.consistent())
    {
        rAtUc = 1.0 / (1.0 / rAUc - UcEqn.H1());
        phiHbyAc +=
            fvc::interpolate(rAtUc() - rAUc) * fvc::snGrad(p_adj_U) * mesh.magSf();
        HbyAc -= (rAUc - rAtUc()) * fvc::grad(p_adj_U);
    }

    tUcEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_adj_U, U_adj_U, phiHbyAc, rAtUc(), MRF);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pcEqn(

                fvm::laplacian(rAUc, p_adj_U) == fvc::div(phiHbyAc)
        );

        pcEqn.setReference(p_adj_URefCell, p_adj_URefValue);
        pcEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi_adj_U = phiHbyAc - pcEqn.flux();
        }
    }

    //#include "adjointContinuityErrs.H"

    // Explicitly relax pressure for adjoint momentum corrector
    p_adj_U.relax();

    // Adjoint momentum corrector
    U_adj_U = HbyAc - rAUc * fvc::grad(p_adj_U);
    U_adj_U.correctBoundaryConditions();
    fvOptions.correct(U_adj_U);

    U_adj_U.storePrevIter();
    p_adj_U.storePrevIter();
    phi_adj_U.storePrevIter();
}
