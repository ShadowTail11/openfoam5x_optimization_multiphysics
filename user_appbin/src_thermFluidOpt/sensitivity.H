/*
 * This header calculates the sensitivity and constraint functions, performs the optimization, and applies Heaviside filtering
 */

{
    Info << "Performing sensitivity analysis" << endl;

    // Update sensitivity functions
    f_sens_TU = (
            - weight_sens_TU * alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_T) * unit_correct_T / unit_correct
    );

    f_sens_TK = (
            - weight_sens_TK * fvc::laplacian((k_fluid-k_solid)*(1+q_factor)*q_factor/((q_factor+gamma)*(q_factor+gamma)), T_adj)/(rho_eff*cp_eff)
    );

    f_sens_TT = (
            - weight_sens_TT * fvc::div(T_adj * U)
    );

    f_sens_T_init = (
            - weight_sens_TU * alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_T) * unit_correct_T / unit_correct
            + weight_sens_TK * fvc::laplacian((k_fluid-k_solid)*(1+q_factor)*q_factor/((q_factor+gamma)*(q_factor+gamma)), T_adj)/(rho_eff*cp_eff)
            + weight_sens_TQ * q_gen
            - weight_sens_TT * fvc::div(T_adj * U)
    );

    g_sens_power_loss_init = (
            -alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_U)
    );

    #include "HeavisideSens.H"

    // Store current iteration's constraint values
    gx[0] = cost_power_loss_ratio;
    gx[1] = cost_vol_frac;

    // Store optimization parameters for use in MMA
    for(i = 0; i < n_cells; i++)
    {
        xmma[i] = x[i];
        dfdx[i] = weight_sens_T * f_sens_T[i] / n_cells;
        dgdx[0][i] = weight_cost_U * g_sens_power_loss[i] / n_cells;
        dgdx[1][i] = g_sens_vol_frac[i] / n_cells;
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
