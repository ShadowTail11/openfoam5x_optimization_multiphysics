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

// Porosity properties
dimensionedScalar alpha_U_max(transportProperties.lookup("alpha_U_max"));   // Flow resistance max reference
scalar q_factor(readScalar(transportProperties.lookup("q_factor")));        // Shape factor used to scale ramp function