/*
 * This header calculates the sensitivity and constraint functions, performs the optimization, and applies Heaviside filtering
 */

{
    Info << "Performing sensitivity analysis" << endl;

    // Update sensitivity functions
    adjoint_sensitivity_init = (
            -alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_U)
    );

    // Calculate cost sensitivity for compliance
    mu = -(1 + qd) * qd / ((qd + gamma) * (qd + gamma)) * (E - Emin) / (2.0 * (1.0 + Po));
    lambda = -(1 + qd) * qd / ((qd + gamma) * (qd + gamma)) * Po * (E - Emin) / ((1.0 + Po) * (1.0 - 2.0 * Po));
    gsenshCompliance = -gradD && (rho * (mu * twoSymm(gradD) + (lambda * I) * tr(gradD)));

    #include "HeavisideSens.H"

    // Store current iteration's constraint values
    gx[0] = cost_vol_frac;

    // Store optimization parameters for use in MMA
    for(i = 0; i < n_cells; i++)
    {
        xmma[i] = x[i];
        dfdx[i] = weight_sens * adjoint_sensitivity[i] / n_cells;
        dgdx[2][i] = gsensCompliance[i] / n_cells;
        dgdx[0][i] = cost_sensitivity[i] / n_cells;
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
    dgamma_switch_ave = (1 - 1 / n_ave) * dgamma_switch_ave + dgamma_switch / n_ave;

    Info << "\n--> The percentage change of cells that switch material is: " << dgamma_switch_ave * 100 << "%\n" << endl;

    #include "HeavisideRho.H"
}
