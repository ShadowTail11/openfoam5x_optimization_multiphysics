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

//// Domain Controls
// Domain definition
int solve_ext(round(readScalar(optimizationProperties.lookup("solve_ext"))));                   // Boolean used to switch between internal and external optimizer
scalar solid_area(readScalar(optimizationProperties.lookup("solid_area")));                     // Is solid region set? (1=yes)
scalar fluid_area(readScalar(optimizationProperties.lookup("fluid_area")));                     // Is fluid region used? (1=yes)
scalar test_area(readScalar(optimizationProperties.lookup("test_area")));                       // Is test region used? (1=yes)
scalar geo_dim(readScalar(optimizationProperties.lookup("geo_dim")));                           // Dimensions of solver (2D vs. 3D)


//// Method of Moving Asymptotes (MMA) Controls
// Optimization properties
scalar mma_init(readScalar(optimizationProperties.lookup("mma_init")));                         // Initial optimization increment
scalar mma_dec(readScalar(optimizationProperties.lookup("mma_dec")));                           // Optimization decrement
scalar mma_inc(readScalar(optimizationProperties.lookup("mma_inc")));                           // Optimization increment
scalar mma_limit(readScalar(optimizationProperties.lookup("mma_limit")));                       // Optimization limit
scalar raa0(readScalar(optimizationProperties.lookup("raa0")));                                 // MMA smoothing parameter
int r_filter(round(readScalar(optimizationProperties.lookup("r_filter"))));                     // Density filter radius
dimensionedScalar b("b", dimensionSet(0,-2,0,0,0,0,0),1.0);                                     // Filter radius of the PDE filter

// RAMP function constant (increases transition rate between fluid and solid)
scalar q_ramp_init(readScalar(optimizationProperties.lookup("q_ramp_init")));                   // Initial shape factor used to scale ramp function
scalar q_ramp = q_ramp_init;                                                                    // Shape factor used to scale ramp function
scalar q_ramp_limit(readScalar(optimizationProperties.lookup("q_ramp_limit")));                 // Shape factor limit
scalar dq_ramp(readScalar(optimizationProperties.lookup("dq_ramp")));                           // Percent increase rate for shape factor (near convergence)
scalar n_rapid(round(readScalar(optimizationProperties.lookup("n_rapid"))));                    // Factor that controls when to increase q_ramp based on proximity to convergence


//// Adjoint Optimization Controls
// Cost function constraints
scalar power_loss_ref(readScalar(optimizationProperties.lookup("power_loss_ref")));             // Power loss reference
scalar compliance_ref(readScalar(optimizationProperties.lookup("compliance_ref")));             // Compliance reference
scalar set_vol_frac(readScalar(optimizationProperties.lookup("set_vol_frac")));                 // Target volume fraction
scalar set_vol_frac_solid(readScalar(optimizationProperties.lookup("set_vol_frac_solid")));     // Target volume fraction for region remaining after initial optimization
scalar set_power_ratio(readScalar(optimizationProperties.lookup("set_power_ratio")));           // Target power loss ratio
scalar set_compliance(readScalar(optimizationProperties.lookup("set_compliance")));             // Target compliance ratio

// Optimization weights
scalar weight_cost_U(readScalar(optimizationProperties.lookup("weight_cost_U")));               // Sensitivity weight for velocity term

// Adjoint unit correctors (used because adjoint variables have inverse units of parent, but OpenFOAM uses same units)
dimensionedScalar unit_correct_U("unit_correct_U", dimensionSet(0,-1,1,0,0,0,0),1.0);           // Adjoint U unit correction

// Convergence controls
scalar merit_tol(readScalar(optimizationProperties.lookup("merit_tol")));                       // Convergence tolerance for merit parameters
scalar gamma_tol(readScalar(optimizationProperties.lookup("gamma_tol")));                       // Convergence tolerance for pseudo density gamma
