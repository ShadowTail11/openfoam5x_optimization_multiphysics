/*
 * This header calculates the sensitivity and constraint functions, performs the optimization, and applies Heaviside filtering
 */

{
    Info << "Performing sensitivity analysis" << endl;

    // Update sensitivity function
    mu = -qd * (qd + 1.0) / ((qd + gamma) * (qd + gamma)) * (Esp - Esp_min) / (2.0 * (1.0 + Po));
    lambda = -qd * (qd + 1.0) / ((qd + gamma) * (qd + gamma)) * Po * (Esp - Esp_min) / ((1.0 + Po) * (1.0 - 2.0 * Po));

    cost_sens_comp0 = -gradD && (rho_solid * (mu * twoSymm(gradD) + (lambda * I) * tr(gradD)));

    #include "HeavisideSensExternal.h"

    // Store current iteration's constraint values
    gx[0] = cost_vol_frac;
    gx[1] = 0;

    // Store optimization parameters for use in MMA
    for(i = 0; i < n_cells; i++)
    {
        xmma[i] = x[i];
        dfdx[i] = weight_sens_tot * weight_sens_comp * cost_sens_comp[i] / n_cells / vol_frac_frozen;           // Sensitivity of objective function
        dgdx[0][i] = weight_sens_tot * weight_sens_vf * cost_sens_vol_frac[i] / n_cells / vol_frac_frozen;      // Sensitivity of constraint function (volume fraction)
        dgdx[1][i] = 0;                                     // Null
    }

    Info << "\nRunning Method of Moving Asymptotes Solver\n" << endl;
    mma.MMAsolver(xmma, dfdx, gx, dgdx);

    // Update pseudo density gamma and track change between previous and current iteration
    gamma_switch_prev = gamma_switch;
    gamma_switch = 0;
    for(i = 0; i < n_cells; i++)
    {
        gamma_switch += abs(std::round(xmma[i]) - std::round(gamma[i]));
        x[i] = xmma[i];
    }
    gamma_switch = gamma_switch / n_cells;
    Info << "\n--> The percentage of cells that switch material is: " << gamma_switch * 100 << "%\n" << endl;

    dgamma_switch = abs((gamma_switch - gamma_switch_prev) / (gamma_switch + SMALL));
    dgamma_switch_ave = (1.0 - 1.0 / n_ave) * dgamma_switch_ave + dgamma_switch / n_ave;

    #include "HeavisideRhoExternal.h"
}