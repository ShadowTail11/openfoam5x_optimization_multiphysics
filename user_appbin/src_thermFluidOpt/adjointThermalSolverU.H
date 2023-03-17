
if (opt == 1)
{
    n_loops = 100;
}
else
{
    n_loops = 40;
}

for (i = 0; i < n_loops; i++)
{
    Info << "\nSolve for adjoint flow-T\n" << endl;

    volVectorField adjointTransposeConvectionb((fvc::grad(U_adj_T) & U));

    tmp<fvVectorMatrix> tUbEqn
    (
            fvm::div(-phi, U_adj_T)
            - adjointTransposeConvectionb
//            + turbulence->divDevReff(U_adj_T)
            - fvm::laplacian(nu_eff, U_adj_T)
            + fvm::Sp(alpha_U, U_adj_T)
            + cp_eff * T_adj * fvc::grad(T) * unit_correct_T
            ==
            fvOptions(U_adj_T)
    );

    fvVectorMatrix& UbEqn = tUbEqn.ref();

    UbEqn.relax();

    fvOptions.constrain(UbEqn);

    solve(UbEqn == -fvc::grad(p_adj_T));

    fvOptions.correct(U_adj_T);

    //****************************************

    volScalarField rAUb(1.0 / UbEqn.A());
    volVectorField HbyAb("HbyAb", U_adj_T);
    HbyAb = rAUb * UbEqn.H();
    tUbEqn.clear();
    surfaceScalarField phiHbyAb("phiHbyAb", fvc::flux(HbyAb));
    adjustPhi(phiHbyAb, U_adj_T, p_adj_T);

    tmp<volScalarField> rAtUb(rAUb);

    if(simple.consistent())
    {
        rAtUb = 1.0 / (1.0 / rAUb - UbEqn.H1());
        phiHbyAb +=
        fvc::interpolate(rAtUb() - rAUb) * fvc::snGrad(p_adj_T) * mesh.magSf();
        HbyAb -= (rAUb - rAtUb()) * fvc::grad(p_adj_T);
    }

    tUbEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_adj_T, U_adj_T, phiHbyAb, rAtUb(), MRF);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pbEqn
        (
                fvm::laplacian(rAUb, p_adj_T) == fvc::div(phiHbyAb)
        );

        pbEqn.setReference(p_adj_TRefCell, p_adj_TRefValue);
        pbEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi_adj_T = phiHbyAb - pbEqn.flux();
        }
    }

//    #include "adjointContinuityErrs.H"

    // Explicitly relax pressure for adjoint momentum corrector
    p_adj_T.relax();

    // Adjoint momentum corrector
    U_adj_T = HbyAb - rAUb * fvc::grad(p_adj_T);
    U_adj_T.correctBoundaryConditions();
    fvOptions.correct(U_adj_T);

    U_adj_T.storePrevIter();
    p_adj_T.storePrevIter();
    phi_adj_T.storePrevIter();
}
