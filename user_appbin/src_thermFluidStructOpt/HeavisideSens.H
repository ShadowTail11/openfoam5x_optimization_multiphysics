/*
 * This header implements the heaviside and PDE filter into the adjoint and cost function sensitivities
 */

Info << "\nCalculating Heaviside_rho sensitivity and implementing PDE filter\n" << endl;

// Normalize sensitivity functions by maximum value of respective functions
obj_sens_T0 = obj_sens_T0 / gMax(mag(obj_sens_T0.primitiveField()));
cost_sens_power_loss0 = cost_sens_power_loss0 / gMax(mag(cost_sens_power_loss0.primitiveField()));
cost_sens_comp0 = cost_sens_comp0 / gMax(mag(cost_sens_comp0.primitiveField()));
cost_sens_vol_frac0.primitiveFieldRef() = 1.0;

// Update fields based on material zones
if (solid_area)
{
    forAll(cells_solid, i)
    {
        obj_sens_T0[cells_solid[i]] = 0;
        cost_sens_power_loss0[cells_solid[i]] = 0;
        cost_sens_comp0[cells_solid[i]] = 0;
        cost_sens_vol_frac0[cells_solid[i]] = 0;
    }
}
if (fluid_area)
{
    forAll(cells_fluid, i)
    {
        obj_sens_T0[cells_fluid[i]] = 0;
        cost_sens_power_loss0[cells_fluid[i]] = 0;
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
    obj_sens_T0 = obj_sens_T0 * drho;
    cost_sens_power_loss0 = cost_sens_power_loss0 * drho;
    cost_sens_comp0 = cost_sens_comp0 * drho;
    cost_sens_vol_frac0 = cost_sens_vol_frac0 * drho;
}

// Use PDE filter to update sensitivities
solve(fvm::laplacian(obj_sens_T) - fvm::Sp(b, obj_sens_T) + obj_sens_T0 * b);
solve(fvm::laplacian(cost_sens_power_loss) - fvm::Sp(b, cost_sens_power_loss) + cost_sens_power_loss0 * b);
solve(fvm::laplacian(cost_sens_vol_frac) - fvm::Sp(b, cost_sens_vol_frac) + cost_sens_vol_frac0 * b);
solve(fvm::laplacian(cost_sens_comp) - fvm::Sp(b, cost_sens_comp) + cost_sens_comp0 * b);

// Rescale sensitivities
obj_sens_T.primitiveFieldRef() = obj_sens_T.primitiveFieldRef() * mesh.V() / gMax(mesh.V());
cost_sens_power_loss.primitiveFieldRef() = cost_sens_power_loss.primitiveFieldRef() * mesh.V() / gMax(mesh.V());
cost_sens_comp.primitiveFieldRef() = cost_sens_comp.primitiveFieldRef() * mesh.V() / gMax(mesh.V());
cost_sens_vol_frac.primitiveFieldRef() = cost_sens_vol_frac.primitiveFieldRef() * mesh.V() / gMax(mesh.V());
