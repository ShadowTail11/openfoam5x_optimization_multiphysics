/*
 * This function implements the primal equation for fluid flow by means of the Navier Stokes equation.
 * The special considerations for this function is that it is only viable for laminar flow, and it
 * uses a ramp function (alpha_U) to slow the velocity in the pseudo porous media, which allows the design
 * to gradually transition between fluid and solid.
 */

{
    // Iterate primal solver to calculate velocity U and kinematic pressure p
    for (i = 0; i < n_flow_solve; i++)
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        // Correct velocity boundary conditions
        MRF.correctBoundaryVelocity(U);

        // Define velocity portion of primal flow equation
        tmp<fvVectorMatrix> tUEqn(

            fvm::div(phi, U)                // Divergent flow term
            + MRF.DDt(U)                    // Material derivative (not necessary for steady-state solver)
            + turbulence->divDevReff(U)     // Used for turbulent flow
            + fvm::Sp(alpha_U, U)           // Creates a pseudo-porosity-based body force with a ramp function
            ==
            fvOptions(U)                    // Enables OpenFOAM options for velocity such as relaxation factors
        );

        fvVectorMatrix& UEqn = tUEqn.ref(); // Creates a virtual pointer to the previous equation

        UEqn.relax();                       // Relaxes the velocity equation to improve stability

        fvOptions.constrain(UEqn);          // Applies constraints to the equation

        // Incorporates pressure into primal equation
        if(simple.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));

            fvOptions.correct(U);
        }

        //***************** Employs SIMPLE algorithm to solve primal flow equation***********************//

        volScalarField rAU(1.0 / UEqn.A());
        volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));
        surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
        MRF.makeRelative(phiHbyA);
        adjustPhi(phiHbyA, U, p);

        tmp<volScalarField> rAtU(rAU);

        if (simple.consistent())
        {
            rAtU = 1.0 / (1.0 / rAU - UEqn.H1());
            phiHbyA +=
                fvc::interpolate(rAtU() - rAU) * fvc::snGrad(p) * mesh.magSf();
            HbyA -= (rAU - rAtU()) * fvc::grad(p);
        }

        tUEqn.clear();

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, U, phiHbyA, rAtU(), MRF);

        // Non-orthogonal pressure corrector loop
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(

                fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
            );

            pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve();

            if (simple.finalNonOrthogonalIter())
            {
                phi = phiHbyA - pEqn.flux();
            }
        }

        // Calculate and print continuity errors for convergence tracking
        #include "continuityErrs.H"

        // Explicitly relax pressure for momentum corrector
        p.relax();

        // Momentum corrector
        U = HbyA - rAtU() * fvc::grad(p);
        U.correctBoundaryConditions();
        fvOptions.correct(U);

        // Store the resulting values
        U.storePrevIter();
        p.storePrevIter();
        phi.storePrevIter();
    }
}
