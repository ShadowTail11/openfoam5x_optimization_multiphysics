Info << "Reading pseudo density field gamma\n" << endl;

// Pseudo density gamma
volScalarField gamma
(
        IOobject
        (
                "gamma",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        mesh,
        scalar(set_vol_frac),
        zeroGradientFvPatchScalarField::typeName
);

// Field used to facilitate density filtering for the pseudo density after optimization
volScalarField gamma_filter
(
        IOobject
        (
                "gamma_filter",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        mesh,
        scalar(set_vol_frac),
        zeroGradientFvPatchScalarField::typeName
);

volScalarField x(gamma_filter);                                     // Field used for Heaviside sensitivity calculation
volScalarField drho(x);                                             // Density factor used in Heaviside function
volScalarField ramp(q_ramp * (1 - gamma) / (q_ramp + gamma));       // RAMP Function

// Set pseudo density values in fixed material regions
labelList cells_solid,cells_fluid,cells_test;
if (solid_area)
{
    word zoneName="solid_area";
    label zoneID=mesh.cellZones().findZoneID(zoneName);
    cells_solid=mesh.cellZones()[zoneID];
    forAll(cells_solid, i)
    {
        x[cells_solid[i]] = 0;
    }
}
if (fluid_area)
{
    word zoneName="fluid_area";
    label zoneID=mesh.cellZones().findZoneID(zoneName);
    cells_fluid=mesh.cellZones()[zoneID];
    forAll(cells_fluid, i)
    {
        x[cells_fluid[i]] = 1.0;
    }
}
if (test_area)
{
    word zoneName="test_area";
    label zoneID=mesh.cellZones().findZoneID(zoneName);
    cells_test=mesh.cellZones()[zoneID];
    forAll(cells_test, i)
    {
        x[cells_test[i]] = 1.0;
        gamma[cells_test[i]] = 1.0;
    }
}
volScalarField cost_sens_vol_frac
(
        IOobject
        (
                "cost_sens_vol_frac",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        gamma,
        zeroGradientFvPatchScalarField::typeName
);
volScalarField cost_sens_vol_frac0(gamma);
