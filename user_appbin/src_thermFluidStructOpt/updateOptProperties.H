//
// This function updates any properties or variables that are related
// to the optimization progressing.
//

Info << "\nUpdating parameters related to optimization process " << endl;

Info << "\nCurrent optimization loop count is: " << opt << endl;

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

// Ensure rapid increase in flow resistance
alpha_U_max = alpha_U_max_init * pow(dalpha / 100 + 1.0, opt);

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
    q_ramp = q_ramp_init * pow(dq_ramp / 100 + 1.0, (opt - opt_threshold));
}

// Limit flow resistance metrics by respective maximum thresholds
alpha_U_max.value() = Foam::min(alpha_U_max.value(), alpha_U_limit.value());
q_ramp = Foam::min(q_ramp, q_ramp_limit);

// Update boundary conditions and parameters based on changes in flow resistance and RAMP function
gamma.correctBoundaryConditions();
ramp = q_ramp * (1 - gamma) / (q_ramp + gamma);
alpha_U = alpha_U_max * ramp;
rho_eff = rho_fluid + (rho_solid - rho_fluid) * ramp;
cp_eff = cp_fluid + (cp_solid - cp_fluid) * ramp;
k_eff = k_fluid + (k_solid - k_fluid) * ramp;
alpha_T_eff = k_eff / (rho_eff * cp_eff);

// Update linear elastic properties that use the RAMP function
ramp_d = qd * (1 - gamma) / (qd + gamma);
mu = (ramp_d * (Esp - Esp_min) + Esp_min) / (2.0 * (1.0 + Po));
lambda = (ramp_d * Po * (Esp - Esp_min) + Po * Esp_min) / ((1.0 + Po) * (1.0 - 2.0 * Po));
if (planeStress)
{
    lambda = (ramp_d * Po * (Esp - Esp_min) + Po * Esp_min) / ((1.0 + Po) * (1.0 - Po));
}
