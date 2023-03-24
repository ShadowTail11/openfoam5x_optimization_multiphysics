{
Info<< "\nImplement Heaviside_rho function\n" << endl;

//***********************************//PDE filter
solve(fvm::laplacian(gamma_densityfilter) - fvm::Sp(b, gamma_densityfilter) + x * b);


//***********************************//Heaviside_filter
del = Foam::min(0.2 * opt, 100.0);
eta0 = 0.0001;
eta1 = 1.0;
y0 = fun(gamma_densityfilter, mesh.V(), del, eta0, n_cells);

do
{
eta5 = (eta0 + eta1) / 2.0;
y5 = fun(gamma_densityfilter, mesh.V(), del, eta5, n_cells);

if(y0 * y5 < 0)
{
eta1 = eta5;
}
else
{
eta0 = eta5;
y0 = y5;
}
} while ((eta1 - eta0) > 0.0001);

for(i = 0; i < n_cells; i++)
{
if(gamma_densityfilter[i] <= eta5)
{
gamma[i] = eta5 * (Foam::exp(-del * (1.0 - gamma_densityfilter[i] / eta5)) - (1.0 - gamma_densityfilter[i] / eta5) * Foam::exp(-del));
}
else
{
gamma[i] = eta5 + (1.0 - eta5) * (1.0 - Foam::exp(-del * (gamma_densityfilter[i] - eta5) / (1.0 - eta5)) + (gamma_densityfilter[i] - eta5) * Foam::exp(-del) / (1.0 - eta5));
}
}

}
