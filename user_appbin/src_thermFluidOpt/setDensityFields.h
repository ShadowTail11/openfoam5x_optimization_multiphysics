Info << "Reading pseudo density field gamma\n" << endl;

// Pseudo density gamma
volScalarField gamma
(
        IOobject
        (
                "gamma",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        mesh,
        scalar(set_vol_frac),
        zeroGradientFvPatchScalarField::typeName
);

// Field used to facilitate density filtering for the pseudo density after optimization
volScalarField gamma_densityfilter
(
        IOobject
        (
                "gamma_densityfilter",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        gamma,
        zeroGradientFvPatchScalarField::typeName
);

// Density field filter used for guiding external solver (freezes fluid region)
volScalarField gamma_conv
(
        IOobject
        (
                "gamma_conv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        gamma,
        zeroGradientFvPatchScalarField::typeName
);

volScalarField x(gamma_densityfilter);                              // Field used for Heaviside sensitivity calculation
volScalarField drho(x);                                             // Density factor used in Heaviside function
volScalarField ramp(q_factor * (1 - gamma) / (q_factor + gamma));   // RAMP Function

// Set pseudo density values in fixed material regions
labelList cells_solid,cells_fluid,cells_test;
if (solid_area)
{
    word zoneName="zone_solid";
    label zoneID=mesh.cellZones().findZoneID(zoneName);
    cells_solid=mesh.cellZones()[zoneID];
    forAll(cells_solid, i)
    {
        //        gamma_test[cells_solid[i]] = 0;
        gamma[cells_solid[i]] = 0;
    }
}
if (fluid_area)
{
    word zoneName="zone_fluid";
    label zoneID=mesh.cellZones().findZoneID(zoneName);
    cells_fluid=mesh.cellZones()[zoneID];
    forAll(cells_fluid, i)
    {
        //        gamma_test[cells_fluid[i]] = 1;
        gamma[cells_fluid[i]] = 1;
    }
}
if (test_area)
{
    word zoneName="zone_test";
    label zoneID=mesh.cellZones().findZoneID(zoneName);
    cells_test=mesh.cellZones()[zoneID];
    forAll(cells_test, i)
    {
        gamma[cells_test[i]] = 1.0;
        x[cells_test[i]] = 1.0;
    }
    q_factor = q_factor_limit;
}
