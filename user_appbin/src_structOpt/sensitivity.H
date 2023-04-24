/*
 * This header calculates the sensitivity and constraint functions, performs the optimization, and applies Heaviside filtering
 */

{
    Info << "Performing sensitivity analysis" << endl;

    cost_sens_comp0 = (
            -gradD && (rho * (3 * gamma * gamma * E * (1 - eps) / (2.0 * (1.0 + Po)) * twoSymm(gradD)
            + (3 * gamma * gamma * Po * E * (1 - eps) / ((1.0 + Po) * (1.0 - 2.0 * Po)) * I) * tr(gradD)))
    );

    #include "heavisideSens.H"

    // Store current iteration's constraint values
    gx[0] = cost_vol_frac;

    // Store optimization parameters for use in MMA
    for (i = 0; i < n_cells; i++)
    {
        xmma[i] = x[i];
        dfdx[i] = cost_sens_comp[i] / n_cells;    //sensitivity of objective function
        dgdx[0][i] = cost_sens_vol_frac[i] / n_cells;    //sensitivity of constraint function
    }

    Info << "\nRunning Method of Moving Asymptotes Solver\n" << endl;
    mma.MMAsolver(xmma, dfdx, gx, dgdx);

    for (i = 0; i < n_cells; i++)
    {
        x[i] = xmma[i];
    }

    #include "heavisideRho.H"
}