/*
 * This function implements the primal equation for heat transfer by means of the energy equation.
 * The function is optionally applied to the test_zone in order to get a reference condition for optimization.
 * This equation assumes that the thermal diffusion is only happening in the fluid region (gamma=1)
 * and that there is heat generation Q in the solid region only (gamma=0).
 */

{
    // Implement thermal solver to calculate temperature T
    Info << "\nSolve for T_test\n" << endl;

    // Compute effective parameters that transition between fluid and solid using the RAMP function for the test zone
    ramp_test = (q_factor * (1 - gamma_test) / (q_factor + gamma_test));    // RAMP Function
    rho_eff_test = (rho_fluid + (rho_solid - rho_fluid) * ramp_test);           // Effective density
    cp_eff_test = (cp_fluid + (cp_solid - cp_fluid) * ramp_test);               // Effective heat capacity
    k_eff_test = (k_fluid + (k_solid - k_fluid) * ramp_test);                   // Effective thermal conductivity
    alpha_T_eff_test = k_eff_test / (rho_eff_test * cp_eff_test);               // Effective thermal diffusivity

    // Define primal energy equation
    fvScalarMatrix T_testEqn
    (
            fvm::div(phi_test, T_test)                      // Thermal diffusivity
            - fvm::laplacian(alpha_T_eff_test, T_test)      // Thermal conduction term
            ==
            fvOptions(T_test) + q_gen * (1 - gamma_test)  // Heat generation in solid region
//                fvOptions(T_test) + q_gen                     // Heat generation in whole domain
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
