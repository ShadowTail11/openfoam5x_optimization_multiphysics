/*
 * This header implements the heaviside and PDE filter into the adjoint and cost function sensitivities
 */

Info << "\nCalculating Heaviside_rho sensitivity and implementing PDE filter\n" << endl;

// Normalize sensitivity functions by maximum value of respective functions
f_sens_T_init = f_sens_T_init / gMax(mag(f_sens_T_init.primitiveField()));
g_sens_power_loss_init = g_sens_power_loss_init / gMax(mag(g_sens_power_loss_init.primitiveField()));
g_sens_vol_frac_norm.primitiveFieldRef() = 1.0;

// Use Heaviside density filter to update sensitivities
if (opt > 1)
{
    for (i = 0; i < n_cells; i++)
    {
        if (gamma_densityfilter[i] <= eta5)
        {
            drho[i] = del * Foam::exp(-del * (1.0 - gamma_densityfilter[i] / eta5)) + Foam::exp(-del);
        }
        else
        {
            drho[i] = del * Foam::exp(-del * (gamma_densityfilter[i] - eta5) / (1.0 - eta5)) + Foam::exp(-del);
        }
    }
    f_sens_T_init = f_sens_T_init * drho;
    g_sens_power_loss_init = g_sens_power_loss_init * drho;
    g_sens_vol_frac_norm = g_sens_vol_frac_norm * drho;
}

// Use PDE filter to update sensitivities
solve(fvm::laplacian(f_sens_T) - fvm::Sp(b, f_sens_T) + f_sens_T_init * b);
solve(fvm::laplacian(g_sens_power_loss) - fvm::Sp(b, g_sens_power_loss) + g_sens_power_loss_init * b);
solve(fvm::laplacian(g_sens_vol_frac) - fvm::Sp(b, g_sens_vol_frac) + g_sens_vol_frac_norm * b);

f_sens_T.primitiveFieldRef() = f_sens_T.primitiveFieldRef() * mesh.V() / gMax(mesh.V());
g_sens_power_loss.primitiveFieldRef() = g_sens_power_loss.primitiveFieldRef() * mesh.V() / gMax(mesh.V());
g_sens_vol_frac.primitiveFieldRef() = g_sens_vol_frac.primitiveFieldRef() * mesh.V() / gMax(mesh.V());
