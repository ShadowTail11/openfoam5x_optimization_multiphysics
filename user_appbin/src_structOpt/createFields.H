
Info << "Reading optimization parameters\n" << endl;
IOdictionary optimizationProperties
(
    IOobject
    (
        "optimizationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//// Domain Controls
// Domain definition
scalar solid_area(readScalar(optimizationProperties.lookup("solid_area")));             // Is solid region set? (1=yes)
scalar fluid_area(readScalar(optimizationProperties.lookup("fluid_area")));             // Is fluid region used? (1=yes)
scalar test_area(readScalar(optimizationProperties.lookup("test_area")));               // Is test region used? (1=yes)
scalar geo_dim(readScalar(optimizationProperties.lookup("geo_dim")));                   // Dimensions of solver (2D vs. 3D)

//// Method of Moving Asymptotes (MMA) Controls
// Optimization properties
scalar mma_init(readScalar(optimizationProperties.lookup("mma_init")));                 // Initial optimization increment
scalar mma_dec(readScalar(optimizationProperties.lookup("mma_dec")));                   // Optimization decrement
scalar mma_inc(readScalar(optimizationProperties.lookup("mma_inc")));                   // Optimization increment
scalar mma_limit(readScalar(optimizationProperties.lookup("mma_limit")));               // Optimization limit
scalar raa0(readScalar(optimizationProperties.lookup("raa0")));                         // MMA smoothing parameter
int r_filter(round(readScalar(optimizationProperties.lookup("r_filter"))));             // Density filter radius
dimensionedScalar b("b", dimensionSet(0,-2,0,0,0,0,0),1.0);                             // Filter radius of the PDE filter

//// Adjoint Optimization Controls
// Cost function constraints
scalar compliance_ref(readScalar(optimizationProperties.lookup("compliance_ref")));     // Compliance reference
scalar set_vol_frac(readScalar(optimizationProperties.lookup("set_vol_frac")));         // Target volume fraction

// Convergence controls
scalar merit_tol(readScalar(optimizationProperties.lookup("merit_tol")));               // Convergence tolerance for merit parameters

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

volScalarField x(gamma_filter);                                     // Field used for Heaviside sensitivity calculation
volScalarField drho(x);                                             // Density factor used in Heaviside function

// Set pseudo density values in fixed material regions
labelList cells_solid,cells_fluid,cells_test;
if (solid_area == 1)
{
    word zoneName="solid_area";
    label zoneID=mesh.cellZones().findZoneID(zoneName);
    cells_solid=mesh.cellZones()[zoneID];
//    forAll(cells_solid, i)
//    {
//        x[cells_solid[i]] = 1.0;
//    }
}
if (fluid_area == 1)
{
    word zoneName="fluid_area";
    label zoneID=mesh.cellZones().findZoneID(zoneName);
    cells_fluid=mesh.cellZones()[zoneID];
    forAll(cells_fluid, i)
    {
        x[cells_fluid[i]] = 0;
    }
}
if (test_area == 1)
{
    word zoneName="test_area";
    label zoneID=mesh.cellZones().findZoneID(zoneName);
    cells_test=mesh.cellZones()[zoneID];
    forAll(cells_test, i)
    {
        x[cells_test[i]] = 0;
        gamma[cells_test[i]] = 0;
    }
}
