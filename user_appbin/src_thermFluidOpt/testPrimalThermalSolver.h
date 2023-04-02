/*
 * This function implements the primal equation for heat transfer by means of the energy equation.
 * The function is optionally applied to the test_zone in order to get a reference condition for optimization.
 * This equation assumes that the thermal diffusion is only happening in the fluid region (gamma=1)
 * and that there is heat generation Q in the solid region only (gamma=0).
 */

{
    // Iterate thermal solver to calculate temperature T
    for (i = 0; i < 2; i++)
    {
        Info << "\nSolve for T_test\n" << endl;

        // Define primal energy equation
        fvScalarMatrix T_testEqn
        (
                fvm::div(phi_test, T_test)                  // Thermal diffusivity
                - fvm::laplacian(alpha_T_eff, T_test)       // Thermal conduction term
                ==
                fvOptions(T_test) + q_gen * (1 - gamma)     // Heat generation in solid region
//                fvOptions(T_test) + q_gen                 // Heat generation in whole domain
        );

        // Relax and solve for temperature T_test
        T_testEqn.relax();
        fvOptions.constrain(T_testEqn);
        T_testEqn.solve();
        fvOptions.correct(T_test);

        // Correct temperature so that minimum does not go below T_min
        forAll(T_test, i)
        {
            T_test[i] = Foam::max(T_test[i], T_min.value());
        }
    }
}
