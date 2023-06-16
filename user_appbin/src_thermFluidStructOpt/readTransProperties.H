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

// Viscosity properties
scalar nu_temp_dep(readScalar(transportProperties.lookup("nu_temp_dep")));  // Use temperature dependent viscosity?(y=1)
dimensionedScalar nu_k(transportProperties.lookup("nu_k"));                 // Viscosity at nu_T_base
dimensionedScalar nu_n(transportProperties.lookup("nu_n"));                 // Fluid behavior index (1 for Newtonian)
dimensionedScalar nu_slope(transportProperties.lookup("nu_slope"));         // Slope of viscosity vs. temperature
dimensionedScalar nu_T_ref(transportProperties.lookup("nu_T_ref"));         // Reference temperature
dimensionedScalar nu_min(transportProperties.lookup("nu_min"));             // Viscosity minimum
dimensionedScalar nu_max(transportProperties.lookup("nu_max"));             // Viscosity maximum

//// Porosity Controls
// Flow resistance (rapidly slows fluid velocity as pseudo density transitions to 0 -> solid)
dimensionedScalar l_char(transportProperties.lookup("l_char"));             // Characteristic length
scalar darcy(readScalar(transportProperties.lookup("darcy")));              // Darcy number
scalar alpha_scale(readScalar(transportProperties.lookup("alpha_scale")));  // Scale factor for early implementation
scalar dalpha(readScalar(transportProperties.lookup("dalpha")));            // Percent increase rate for flow resistance

// Pressure controls
scalar P_operating(readScalar(transportProperties.lookup("P_operating")));  // Operating pressure of the fluid
scalar P_ref(readScalar(transportProperties.lookup("P_ref")));              // Reference pressure
