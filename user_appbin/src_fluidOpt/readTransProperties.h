//
// This header file reads user-defined properties from the transportProperties file in the application's
// 'constant' directory
//

Info << "Reading transport properties\n" << endl;

// Establishing link to transportProperties file
IOdictionary transportProperties(

        IOobject(

                "transportProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
        )
);

// Material densities
dimensionedScalar rho_fluid(transportProperties.lookup("rho_fluid"));               // Fluid density
dimensionedScalar rho_solid(transportProperties.lookup("rho_solid"));               // Solid density

//// Porosity controls
// Flow resistance (rapidly slows fluid velocity as pseudo density transitions to 0 -> solid)
dimensionedScalar alpha_U_max(transportProperties.lookup("alpha_U_max"));           // Flow resistance max reference
dimensionedScalar alpha_U_limit(transportProperties.lookup("alpha_U_limit"));       // Flow resistance limit
int n_rapid(round(readScalar(transportProperties.lookup("n_rapid"))));              // Number of initial iterations under rapid change
scalar dalpha_init(readScalar(transportProperties.lookup("dalpha_init")));          // Initial percent increase rate for flow resistance (before n_rapid iterations)
                                                                                    //  (this factor is multiplied by loop count for rapid increase)
scalar dalpha_end(readScalar(transportProperties.lookup("dalpha_end")));            // Percent increase rate for flow resistance (near convergence)

// RAMP function constant (increases transition rate between fluid and solid)
scalar q_factor(readScalar(transportProperties.lookup("q_factor")));                // Shape factor used to scale ramp function
scalar q_factor_limit(readScalar(transportProperties.lookup("q_factor_limit")));    // Shape factor limit
scalar dq_factor(readScalar(transportProperties.lookup("dq_factor")));              // Percent increase rate for shape factor (near convergence)
