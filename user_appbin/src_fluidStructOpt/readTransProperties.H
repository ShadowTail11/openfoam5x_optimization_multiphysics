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

//// Porosity Controls
// Flow resistance (rapidly slows fluid velocity as pseudo density transitions to 0 -> solid)
dimensionedScalar l_char(transportProperties.lookup("l_char"));                     // Characteristic length
scalar darcy(readScalar(transportProperties.lookup("darcy")));                      // Darcy number
scalar alpha_scale(readScalar(transportProperties.lookup("alpha_scale")));          // Scale factor for early implementation
scalar dalpha(readScalar(transportProperties.lookup("dalpha")));                    // Percent increase rate for flow resistance

// Pressure controls
scalar P_operating(readScalar(transportProperties.lookup("P_operating")));          // Operating pressure of the fluid
scalar P_ref(readScalar(transportProperties.lookup("P_ref")));                      // Reference pressure
