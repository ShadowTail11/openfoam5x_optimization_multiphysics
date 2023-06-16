/*
 * This header calculates the sensitivity and constraint functions, performs the optimization, and applies Heaviside filtering
 */

{
    Info << "Performing sensitivity analysis" << endl;

    // Update sensitivity functions
    obj_sens_T0 = (
            - weight_sens_TU * alpha_U_max * (1 + q_ramp) * q_ramp / ((q_ramp + gamma) * (q_ramp + gamma)) * (U & U_adj_T) * unit_correct_T / unit_correct
            + weight_sens_TK * fvc::laplacian((k_fluid-k_solid)*(1+q_ramp)*q_ramp/((q_ramp+gamma)*(q_ramp+gamma)), T_adj)/(rho_eff*cp_eff)
            + weight_sens_TQ * q_gen
            - weight_sens_TT * fvc::div(T_adj * U)
    );

    cost_sens_power_loss0 = (
            -alpha_U_max * (1 + q_ramp) * q_ramp / ((q_ramp + gamma) * (q_ramp + gamma)) * (U & U_adj_U)
    );

    mu = -(1 + qd) * qd / ((qd + gamma) * (qd + gamma)) * (Esp - Esp_min) / (2.0 * (1.0 + Po));
    lambda = -(1 + qd) * qd / ((qd + gamma) * (qd + gamma)) * Po * (Esp - Esp_min) / ((1.0 + Po) * (1.0 - 2.0 * Po));
    cost_sens_comp0 = -gradD && (rho_solid * (mu * twoSymm(gradD) + (lambda * I) * tr(gradD)));

    #include "HeavisideSens.H"

    // Store current iteration's constraint value (volume fraction)
    gx[0] = cost_vol_frac;

    // Store optimization parameters for use in MMA
    for(i = 0; i < n_cells; i++)
    {
        xmma[i] = x[i];
        dfdx[i] = ((weight_sens_T * obj_sens_T[i]
                + weight_sens_U * cost_sens_power_loss[i]
                + weight_sens_comp * cost_sens_comp[i]) * weight_sens_tot / n_cells);              // Sensitivities of objective function (temperature)
        dgdx[0][i] = weight_sens_tot * weight_sens_vf * cost_sens_vol_frac[i] / n_cells;  // Sensitivities of constraint function (volume fraction)
    }

    i_cost = 0;
    if (set_power_ratio != -1){
        i_cost++;
        gx[i_cost] = cost_power_loss_ratio;
        for(i = 0; i < n_cells; i++)
        {
            dgdx[i_cost][i] = weight_sens_tot * weight_sens_U * cost_sens_power_loss[i] / n_cells; // Sensitivities of constraint function (power loss)
            dfdx[i] = dfdx[i] - dgdx[i_cost][i];
        }
    }
    if (set_compliance != -1){
        i_cost++;
        gx[i_cost] = cost_compliance;
        for(i = 0; i < n_cells; i++)
        {
            dgdx[i_cost][i] = weight_sens_tot * weight_sens_comp * cost_sens_comp[i] / n_cells;    // Sensitivities of constraint function (compliance)
            dfdx[i] = dfdx[i] - dgdx[i_cost][i];
        }
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