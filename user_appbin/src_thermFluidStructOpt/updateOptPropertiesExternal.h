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
        gamma[cells_solid[i]] = 0.0;
    }
}

if (fluid_area == 1)
{
    forAll(cells_fluid, i)
    {
        gamma[cells_fluid[i]] = 1.0;
    }
}

// Update boundary conditions and parameters based on changes in flow resistance and RAMP function
gamma.correctBoundaryConditions();

// Update linear elastic properties that use the RAMP function
mu = (qd * (1 - gamma) / (qd + gamma) * (Esp - Esp_min) + Esp_min) / (2.0 * (1.0 + Po));
lambda = (qd * (1 - gamma) / (qd + gamma) * Po * (Esp - Esp_min) + Po * Esp_min) / ((1.0 + Po) * (1.0 - 2.0 * Po));
if (planeStress)
{
    lambda = (qd * (1 - gamma) / (qd + gamma) * Po * (Esp - Esp_min) + Po * Esp_min) / ((1.0 + Po) * (1.0 - Po));
}
