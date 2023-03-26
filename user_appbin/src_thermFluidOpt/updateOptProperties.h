//
// This function updates any properties or variables that are related
// to the optimization progressing.
//

Info << "\nUpdating parameters related to optimization process " << endl;

// Reset permanent solid and fluid regions, respectively
if (solid_area == 1)
{
    setCells(gamma, cells_solid, 0);
}

if (fluid_area == 1)
{
    setCells(gamma, cells_fluid, 1);
}

// Ensure rapid increase in flow resistance in early iterations
if (opt <= 63)
{
    alpha_U_max = alpha_U_max * (opt / 7.0 + 1.0);
}

// Mark threshold for when within 2x of convergence threshold
if (opt_threshold == -1 && opt > 50 && dgamma_switch_ave < 2 * gamma_tol && power_loss_conv < 2 * merit_tol && dT_drop_ave < 2 * merit_tol)
{
    opt_threshold = opt;
    Info << "\nNearing convergence. Optimization threshold updated to: " << opt_threshold << endl;
    ofstream outfile8("CheckThreshold.txt", std::ios::app);
    outfile8 << opt_threshold << "\n";
    outfile8.close();
}

// Increase flow resistance and RAMP function as temperature drop approaches convergence
if (opt_threshold > 0 && opt > opt_threshold)
{
    alpha_U_max = alpha_U_max * 1.05;
    q_factor = 0.005 + (opt - opt_threshold) * 1e-4;
}

// Limit flow resistance to a maximum threshold
alpha_U_max.value() = Foam::min(alpha_U_max.value(), 10000000.0);

// Update boundary conditions and parameters based on changes in flow resistance and RAMP function
gamma.correctBoundaryConditions();
q_factor = Foam::min(q_factor, 0.01);
alpha_U = alpha_U_max * q_factor * (1 - gamma) / (q_factor + gamma);
alpha_T_eff = (k_solid + (k_fluid - k_solid) * gamma * (1 + q_factor) / (q_factor + gamma)) / (rho_eff * cp_eff);
