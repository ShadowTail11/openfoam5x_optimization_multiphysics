//
// This function updates any properties or variables that are related
// to the optimization progressing.
//

Info << "\nUpdating parameters related to optimization process " << endl;

Info << "\nCurrent optimization loop count is: " << opt << endl;

// Reset permanent solid and fluid regions, respectively
if (solid_area)
{
    forAll(cells_solid, i)
    {
        gamma[cells_solid[i]] = 0.0;
    }
}

if (fluid_area)
{
    forAll(cells_fluid, i)
    {
        gamma[cells_fluid[i]] = 1.0;
    }
}

if (test_area)
{
    forAll(cells_test_fluid, i)
    {
        gamma[cells_test_fluid[i]] = 1.0;
    }

    forAll(cells_test_solid, i)
    {
        gamma[cells_test_solid[i]] = 0.0;
    }
}

// Ensure rapid increase in flow resistance
alpha_scale = Foam::min(1.0, alpha_scale * pow(dalpha / 100 + 1.0, opt));

// Mark threshold for when within 2x of convergence threshold
if (opt_threshold == -1 && opt > 50 && dgamma_switch_ave < n_rapid * gamma_tol && power_loss_conv < n_rapid * merit_tol)
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

// Update boundary conditions and parameters based on changes in flow resistance and RAMP function
gamma.correctBoundaryConditions();
q_ramp = Foam::min(q_ramp, q_ramp_limit);

// Update boundary conditions and parameters based on changes in flow resistance and RAMP function
ramp = q_ramp * (1 - gamma) / (q_ramp + gamma);
rho_eff = rho_fluid + (rho_solid - rho_fluid) * ramp;

// Calculate flow resistance
alpha_U = alpha_scale * alpha_U_max * ramp;

// Update linear elastic properties that use the RAMP function
mu = (qd * (1 - gamma) / (qd + gamma) * (Esp - Esp_min) + Esp_min) / (2.0 * (1.0 + Po));
lambda = (qd * (1 - gamma) / (qd + gamma) * Po * (Esp - Esp_min) + Po * Esp_min) / ((1.0 + Po) * (1.0 - 2.0 * Po));
if (planeStress)
{
    lambda = (qd * (1 - gamma) / (qd + gamma) * Po * (Esp - Esp_min) + Po * Esp_min) / ((1.0 + Po) * (1.0 - Po));
}
