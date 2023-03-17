{
    Info << "Perform sensitivity analysis" << endl;

    f_sens_T = (
        -alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_T)
        + unit_correct * (k_fluid - k_solid) * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma))
        * (fvc::grad(T) & fvc::grad(T_adj)) / (rho_eff * cp_eff)
        + unit_correct * T_adj * q_gen
    );

    g_sens_power_loss = (
            -alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_U)
    );

    #include "HeavisideSens.H"

    gx[0] = DissPower / Foam::max(set_power_ratio * 3 - opt * 0.2, set_power_ratio) - 1.0;
    gx[1] = V;

    VecGetArray(dfdx, &dfdx_array);
    VecGetArrays(dgdx, m, &dgdx_array);

    for(i = 0; i < allcells; i++)
    {
      dfdx_array[i] = f_sens_T[i] / scalef / nallcells;
      dgdx_array[0][i] = g_sens_power_loss[i] / scaleg / nallcells;
      dgdx_array[1][i] = gsens[i];
    }

    VecRestoreArray(dfdx, &dfdx_array);
    VecRestoreArrays(dgdx, m, &dgdx_array);

    Info << "\nRun mmasub\n" << endl;

    mma->SetOuterMovelimit(Xmin, Xmax, mma_limit, xmma, xmin, xmax);
    mma->Update(xmma, dfdx, gx, dgdx, xmin, xmax, raa0, mma_dec, mma_inc);

    VecGetArray(xmma, &xmma_array);

    // Update values

    for(i = 0; i < allcells; i++)
    {
      gamma[i] = xmma_array[i];
    }

    VecRestoreArray(xmma, &xmma_array);

    #include "HeavisideRho.H"

    opt++;

    if(opt <= 63)
    {
      alpha_U_max = alpha_U_max * (opt / 7.0 + 1.0);
    }

    if(opt > 100)
    {
      alpha_U_max = alpha_U_max * 1.05;
    }

    alpha_U_max.value() = Foam::min(alpha_U_max.value(), 10000000.0);

    if(opt >= 80)
    {
      q_factor = 0.005 + (opt - 80) * 1e-4;
    }

    q_factor = Foam::min(q_factor, 0.01);
    alpha_U = alpha_U_max * q_factor * (1 - gamma) / (q_factor + gamma);
    alpha_T_eff = (k_solid + (k_fluid - k_solid) * gamma * (1 + q_factor) / (q_factor + gamma)) / (rho_eff * cp_eff);
}