/*
 * This header calculates the sensitivity and constraint functions, performs the optimization, and applies Heaviside filtering
 */

{
    Info << "Performing sensitivity analysis" << endl;

    // Update sensitivity functions
    fsenshMeanT = (
            - weight_sens_TU * alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_T) * unit_correct_T / unit_correct
            + weight_sens_TK * fvc::laplacian((k_fluid-k_solid)*(1+q_factor)*q_factor/((q_factor+gamma)*(q_factor+gamma)), T_adj)/(rho_eff*cp_eff)
            + weight_sens_TQ * q_gen
            - weight_sens_TT * fvc::div(T_adj * U)
    );

    gsenshPowerDiss = (
            -alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_U)
    );

    mu=-(1+qd)*qd/((qd+gamma)*(qd+gamma))*(E-Emin)/(2.0*(1.0+Po));
    lambda=-(1+qd)*qd/((qd+gamma)*(qd+gamma))*Po*(E-Emin)/((1.0+Po)*(1.0-2.0*Po));
    gsenshCompliance=-gradD && (rho*(mu*twoSymm(gradD)+ (lambda*I)*tr(gradD)));

    //链式求导计算对物理密度的灵敏度
    #include "HeavisideSens.H"

    //往MMA中传入数据
    gx[0] = cost_power_loss_ratio;
    gx[1] = (C/Foam::max(CRelax-opt*0.2,CMax)-1.0);
    gx[2]= cost_vol_frac;

    for(i=0;i<n_cells;i++)
    {
     xmma[i] = x[i];
     dfdx[i] = weight_sens_T * fsensMeanT[i] / n_cells;//sensitivities of objective function
     dgdx[0][i] = weight_cost_U * gsensPowerDiss[i] / n_cells;//sensitivities of constraint function
     dgdx[1][i] = gsensCompliance[i] / n_cells;//sensitivities of constraint function
     dgdx[2][i] = gsensVol[i] / n_cells;//sensitivities of constraint function
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