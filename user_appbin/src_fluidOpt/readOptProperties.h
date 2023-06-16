//
// This header file reads user-defined properties from the optimizationProperties file in the application's
// 'constant' directory
//

Info << "Reading optimization properties\n" << endl;

// Establishing link to optimizationProperties file
IOdictionary optimizationProperties(

        IOobject(

                "optimizationProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
        )
);

// Domain definition
scalar solid_area(readScalar(optimizationProperties.lookup("solid_area")));             // Is solid region set? (1=yes)
scalar fluid_area(readScalar(optimizationProperties.lookup("fluid_area")));             // Is fluid region used? (1=yes)
scalar test_area(readScalar(optimizationProperties.lookup("test_area")));               // Is test region used? (1=yes)
scalar geo_dim(readScalar(optimizationProperties.lookup("geo_dim")));                   // Dimensions of solver (2D vs. 3D)

// Optimization properties
scalar mma_init(readScalar(optimizationProperties.lookup("mma_init")));                 // Initial optimization increment
scalar mma_dec(readScalar(optimizationProperties.lookup("mma_dec")));                   // Optimization decrement
scalar mma_inc(readScalar(optimizationProperties.lookup("mma_inc")));                   // Optimization increment
scalar mma_limit(readScalar(optimizationProperties.lookup("mma_limit")));               // Optimization limit
scalar raa0(readScalar(optimizationProperties.lookup("raa0")));                         // MMA smoothing parameter
int r_filter(round(readScalar(optimizationProperties.lookup("r_filter"))));             // Density filter radius
dimensionedScalar b("b", dimensionSet(0,-2,0,0,0,0,0),1.0);                             // Filter radius of the PDE filter

// Adjoint solver controls
scalar vol_frac_ref(readScalar(optimizationProperties.lookup("vol_frac_ref")));         // Reference volume fraction
scalar vol_frac_ratio(readScalar(optimizationProperties.lookup("vol_frac_ratio")));     // Ratio used to calculate the set volume fraction
scalar power_loss_ref(readScalar(optimizationProperties.lookup("power_loss_ref")));     // Power loss reference
scalar weight_sens(readScalar(optimizationProperties.lookup("weight_sens")));           // Weight factor for sensitivity function
scalar set_vol_frac = vol_frac_ref * vol_frac_ratio;                                    // Target volume fraction

// Convergence controls
int n_flow_solve(round(readScalar(optimizationProperties.lookup("n_flow_solve"))));     // Number of iterations used for solver per optimization loop
scalar merit_tol(readScalar(optimizationProperties.lookup("merit_tol")));               // Convergence tolerance for merit parameters
scalar gamma_tol(readScalar(optimizationProperties.lookup("gamma_tol")));               // Convergence tolerance for pseudo density gamma
