/*
 * This header implements the heaviside and PDE filter into the adjoint and cost function sensitivities
 */

Info << "\nCalculating Heaviside_rho sensitivity and implementing PDE filter\n" << endl;

// Normalize sensitivity functions by maximum value of respective functions
cost_sens_comp0 = cost_sens_comp0 / gMax(mag(cost_sens_comp0.primitiveField()));
cost_sens_vol_frac0.primitiveFieldRef() = 1.0;

// Update fields based on material zones
if (fluid_area == 1)
{
    forAll(cells_fluid, i)
    {
        cost_sens_comp0[cells_fluid[i]] = 0;
        cost_sens_vol_frac0[cells_fluid[i]] = 0;
    }
}

// Use Heaviside density filter to update sensitivities
if (opt > 1)
{
    for (i = 0; i < n_cells; i++)
    {
        if (gamma_filter[i] <= eta5)
        {
            drho[i] = del * Foam::exp(-del * (1.0 - gamma_filter[i] / eta5)) + Foam::exp(-del);
        }
        else
        {
            drho[i] = del * Foam::exp(-del * (gamma_filter[i] - eta5) / (1.0 - eta5)) + Foam::exp(-del);
        }
    }
    cost_sens_comp0 = cost_sens_comp0 * drho;
    cost_sens_vol_frac0 = cost_sens_vol_frac0 * drho;
}

// Use PDE filter to update sensitivities
solve(fvm::laplacian(cost_sens_vol_frac) - fvm::Sp(b, cost_sens_vol_frac) + cost_sens_vol_frac0 * b);
solve(fvm::laplacian(cost_sens_comp) - fvm::Sp(b, cost_sens_comp) + cost_sens_comp0 * b);

// Rescale sensitivities
cost_sens_comp.primitiveFieldRef() = cost_sens_comp.primitiveFieldRef() * mesh.V() / gMax(mesh.V());
cost_sens_vol_frac.primitiveFieldRef() = cost_sens_vol_frac.primitiveFieldRef() * mesh.V() / gMax(mesh.V());
