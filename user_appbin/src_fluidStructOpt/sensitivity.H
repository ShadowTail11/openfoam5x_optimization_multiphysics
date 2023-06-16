/*
 * This header calculates the sensitivity and constraint functions, performs the optimization, and applies Heaviside filtering
 */

{
    Info << "Performing sensitivity analysis" << endl;

    // Update sensitivity functions
    cost_sens_power_loss0 = (
            -alpha_scale * alpha_U_max * (1 + q_ramp) * q_ramp / ((q_ramp + gamma) * (q_ramp + gamma)) * (U & U_adj_U)
    );

    mu = -(1 + qd) * qd / ((qd + gamma) * (qd + gamma)) * (Esp - Esp_min) / (2.0 * (1.0 + Po));
    lambda = -(1 + qd) * qd / ((qd + gamma) * (qd + gamma)) * Po * (Esp - Esp_min) / ((1.0 + Po) * (1.0 - 2.0 * Po));
    cost_sens_comp0 = -gradD && (rho_solid * (mu * twoSymm(gradD) + (lambda * I) * tr(gradD)));

    #include "HeavisideSens.H"

    // Store current iteration's constraint values
    gx[0] = cost_compliance;
    gx[1] = cost_vol_frac;

    // Store optimization parameters for use in MMA
    for(i = 0; i < n_cells; i++)
    {
        xmma[i] = x[i];
        dfdx[i] = weight_sens_tot * weight_sens_U * cost_sens_power_loss[i] / n_cells;    // Sensitivity of objective function (power loss)
        dgdx[0][i] = weight_sens_tot * weight_sens_comp * cost_sens_comp[i] / n_cells;    // Sensitivity of constraint function (compliance)
        dgdx[1][i] = weight_sens_tot * weight_sens_vf * cost_sens_vol_frac[i] / n_cells;  // Sensitivity of constraint function (volume fraction)
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

    #include "HeavisideRho.H"
}