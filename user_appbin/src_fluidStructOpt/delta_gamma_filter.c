//
// This function returns the difference between the pseudo density gamma
// and the implementation of the heaviside filter on gamma.
//

scalar delta_gamma_filter(volScalarField &gamma,
                          const scalarField &cell_volume,
                          double del,
                          double eta,
                          int cell_count)
{
    int i;
    scalar z = 0;
    double *x = new double[cell_count];

    for (i = 0; i < cell_count; i++)
    {
        if (gamma[i] <= eta)
        {
            x[i] = (eta * (Foam::exp(-del * (1 - gamma[i] / eta))
                    - (1 - gamma[i] / eta) * Foam::exp(-del)));
        }
        else
        {
            x[i] = (eta + (1 - eta) * (1 - Foam::exp(-del * (gamma[i] - eta) / (1 - eta))
                    + (gamma[i] - eta) * Foam::exp(-del) / (1 - eta)));
        }
    }
    for (i = 0; i < cell_count; i++)
    {
        z = z + (gamma[i] - x[i]) * cell_volume[i];
    }
    delete x;
    return {z};
}
