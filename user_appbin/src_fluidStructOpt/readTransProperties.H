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
dimensionedScalar rho_solid(transportProperties.lookup("rho_solid"));

//// Porosity Controls
// Flow resistance (rapidly slows fluid velocity as pseudo density transitions to 0 -> solid)
dimensionedScalar alpha_U_max_init(transportProperties.lookup("alpha_U_max_init")); // Initial flow resistance max reference
dimensionedScalar alpha_U_max = alpha_U_max_init;                                   // Flow resistance max reference
dimensionedScalar alpha_U_limit(transportProperties.lookup("alpha_U_limit"));       // Flow resistance limit
scalar dalpha(readScalar(transportProperties.lookup("dalpha")));                    // Percent increase rate for flow resistance

// Pressure controls
scalar P_operating(readScalar(transportProperties.lookup("P_operating")));          // Operating pressure of the fluid
scalar P_ref(readScalar(transportProperties.lookup("P_ref")));                      // Reference pressure
