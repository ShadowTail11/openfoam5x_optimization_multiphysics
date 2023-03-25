{
    Info << "Perform sensitivity analysis" << endl;

    f_sens_T_init = (
            - alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_T) * unit_correct_T / unit_correct
            + fvc::laplacian((k_fluid-k_solid)*(1+q_factor)*q_factor/((q_factor+gamma)*(q_factor+gamma)), T_adj)/(rho_eff*cp_eff)
            + q_gen
            - weight_TU * fvc::div(T_adj * U)
    );

    g_sens_power_loss_init = (
            -alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_U)
    );

    #include "HeavisideSens.H"

    gx[0] = cost_power_loss_ratio;
    gx[1] = cost_vol_frac;

    VecGetArray(dfdx, &dfdx_array);
    VecGetArrays(dgdx, m, &dgdx_array);

    for(i = 0; i < n_cells; i++)
    {
        dfdx_array[i] = f_sens_T[i] / n_cells;
        dgdx_array[0][i] = g_sens_power_loss[i] / n_cells;
        dgdx_array[1][i] = g_sens_vol_frac[i] / n_cells;
    }

    VecRestoreArray(dfdx, &dfdx_array);
    VecRestoreArrays(dgdx, m, &dgdx_array);

    Info << "\nRun mmasub\n" << endl;

    mma->SetOuterMovelimit(Xmin, Xmax, mma_limit, xmma, xmin, xmax);
    mma->Update(xmma, dfdx, gx, dgdx, xmin, xmax, raa0, mma_dec, mma_inc);

    VecGetArray(xmma, &xmma_array);

    // Update values

    gamma_switch_prev = gamma_switch;
    gamma_switch = 0;
    for(i = 0; i < n_cells; i++)
    {
        gamma_switch += abs(std::round(xmma_array[i]) - std::round(gamma[i]));
        x[i] = xmma_array[i];
    }
    gamma_switch = gamma_switch / n_cells;
    Info << "\n--> The percentage of cells that switch material is: " << gamma_switch << endl;

    dgamma_switch = abs((gamma_switch - gamma_switch_prev) / (gamma_switch + SMALL));
    dgamma_switch_ave = (1 - 1 / n_ave) * dgamma_switch_ave + dgamma_switch / n_ave;

    VecRestoreArray(xmma, &xmma_array);

    #include "HeavisideRho.H"

    opt++;
}