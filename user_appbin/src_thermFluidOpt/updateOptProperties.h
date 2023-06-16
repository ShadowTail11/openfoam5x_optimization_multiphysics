//
// This function updates any properties or variables that are related
// to the optimization progressing.
//

Info << "\nUpdating parameters related to optimization process " << endl;

// Reset permanent solid and fluid regions, respectively
if (solid_area == 1)
{
    forAll(cells_solid, i)
    {
        gamma[cells_solid[i]] = 0;
    }
}

if (fluid_area == 1)
{
    forAll(cells_fluid, i)
    {
        gamma[cells_fluid[i]] = 1;
    }
}

if (test_area == 1)
{
    forAll(cells_fluid, i)
    {
        gamma[cells_fluid[i]] = 1;
    }
}

// Ensure rapid increase in flow resistance
alpha_scale = Foam::min(1.0, alpha_scale * pow(dalpha / 100 + 1.0, opt));

// Mark threshold for when within 2x of convergence threshold
if (opt_threshold == -1 && opt > 50 && dgamma_switch_ave < n_rapid * gamma_tol && power_loss_conv < n_rapid * merit_tol && dT_drop_ave < n_rapid * merit_tol)
{
    opt_threshold = opt;
    Info << "\nConvergence threshold reached. Optimization threshold updated to: " << opt_threshold << endl;
    Info << "RAMP function shape factor will begin increasing at prescribed rate." << endl;
}

// Increase flow resistance and RAMP function as temperature drop approaches convergence
if (opt_threshold > 0 && opt > opt_threshold)
{
    q_factor *= dq_factor / 100 + 1.0;
}

// Update boundary conditions and parameters based on changes in flow resistance and RAMP function
gamma.correctBoundaryConditions();
q_factor = Foam::min(q_factor, q_factor_limit);

// RAMP factor
ramp =  q_factor * (1 - gamma) / (q_factor + gamma);
rho_eff = rho_fluid + (rho_solid - rho_fluid) * ramp;
cp_eff = cp_fluid + (cp_solid - cp_fluid) * ramp;
k_eff = k_fluid + (k_solid - k_fluid) * ramp;
alpha_T_eff = k_eff / (rho_eff * cp_eff);

// Calculate flow resistance
alpha_U = alpha_scale * alpha_U_max * ramp;
