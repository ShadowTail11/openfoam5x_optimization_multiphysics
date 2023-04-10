/*
 * This function implements the adjoint equation for heat transfer by means of the energy equation.
 * This equation includes heat generation as a factor in the function.
 */

{
    // Iterate adjoint thermal solver to calculate adjoint temperature T_adj
    Info << "\nSolving for adjoint temperature T_adj\n" << endl;

    fvScalarMatrix TEqnb
            (
                    fvm::div(-phi, T_adj)                   // Thermal diffusion
                    - fvm::laplacian(alpha_T_eff, T_adj)    // Thermal conduction
                    ==
                    fvOptions(T_adj) - q_gen * (1 - gamma)  // Heat generation
            );

    // Relax and solve for temperature T_adj
    TEqnb.relax();
    fvOptions.constrain(TEqnb);
    TEqnb.solve();
    fvOptions.correct(T_adj);
}
