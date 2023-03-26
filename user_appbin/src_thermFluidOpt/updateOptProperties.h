//
// This function updates any properties or variables that are related
// to the optimization progressing.
//

if (solid_area == 1)
{
    setCells(gamma, cells_solid, 0);
}

if (fluid_area == 1)
{
    setCells(gamma, cells_fluid, 1);
}

if (opt <= 63)
{
    alpha_U_max = alpha_U_max * (opt / 7.0 + 1.0);
}

if (opt > 300) // Previously 100
{
    alpha_U_max = alpha_U_max * 1.05;
}

alpha_U_max.value() = Foam::min(alpha_U_max.value(), 10000000.0);

if (opt >= 250) // previously 80 here and in (opt-x) below
{
    q_factor = 0.005 + (opt - 250) * 1e-4;
}

gamma.correctBoundaryConditions();
q_factor = Foam::min(q_factor, 0.01);
alpha_U = alpha_U_max * q_factor * (1 - gamma) / (q_factor + gamma);
alpha_T_eff = (k_solid + (k_fluid - k_solid) * gamma * (1 + q_factor) / (q_factor + gamma)) / (rho_eff * cp_eff);
