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

// Density field filter used for guiding external solver (initializes wall)
volScalarField gamma_wall
(
        IOobject
        (
                "gamma_wall",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        gamma,
        zeroGradientFvPatchScalarField::typeName
);

volScalarField x(gamma_filter);                                     // Field used for Heaviside sensitivity calculation
volScalarField drho(x);                                             // Density factor used in Heaviside function
volScalarField ramp(q_ramp * (1 - gamma) / (q_ramp + gamma));       // RAMP Function

// A field that sets the initial solid region
volScalarField vol_set_solid
(
        IOobject
        (
                "vol_set_solid",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        gamma * 0.0,
        zeroGradientFvPatchScalarField::typeName
);

// A field that sets the initial fluid region
volScalarField vol_set_fluid
(
        IOobject
        (
                "vol_set_fluid",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        gamma * 0.0,
        zeroGradientFvPatchScalarField::typeName
);

// A field that sets the initial fluid region
volScalarField vol_set_test
(
        IOobject
        (
                "vol_set_test",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        gamma * 0.0,
        zeroGradientFvPatchScalarField::typeName
);

// A field that allows groovyBC to access the solver_ext Boolean
volScalarField p_operating_field
(
        IOobject
        (
                "p_operating_field",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        gamma * 0.0,
        zeroGradientFvPatchScalarField::typeName
);

// Set pseudo density values in fixed material regions
labelList cells_solid,cells_fluid,cells_test;
if (solid_area)
{
    word zoneName="solid_area";
    label zoneID=mesh.cellZones().findZoneID(zoneName);
    cells_solid=mesh.cellZones()[zoneID];
    forAll(cells_solid, i)
    {
        x[cells_solid[i]] = 0.0;
        vol_set_solid[cells_solid[i]] = 1.0;
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
        vol_set_fluid[cells_solid[i]] = 1.0;
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
        vol_set_test[cells_solid[i]] = 1.0;
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
