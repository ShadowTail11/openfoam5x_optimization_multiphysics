//
// This header file reads user-defined properties from the thermalProperties file in the application's
// 'constant' directory
//

Info << "Reading thermal properties\n" << endl;

// Establishing link to thermalProperties file
IOdictionary thermalProperties(

        IOobject(

                "thermalProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
        )
);

// Fluid thermal properties
dimensionedScalar k_fluid(thermalProperties.lookup("k_fluid"));     // Fluid thermal conductivity
dimensionedScalar rho_fluid(thermalProperties.lookup("rho_fluid")); // Fluid density
dimensionedScalar cp_fluid(thermalProperties.lookup("cp_fluid"));   // Fluid specific heat capacity

// Solid thermal properties
dimensionedScalar k_solid(thermalProperties.lookup("k_solid"));     // Solid thermal conductivity
dimensionedScalar rho_solid(thermalProperties.lookup("rho_solid")); // Solid density
dimensionedScalar cp_solid(thermalProperties.lookup("cp_solid"));   // Solid specific heat capacity

// Heat transfer properties
dimensionedScalar q_gen(thermalProperties.lookup("q_gen"));         // Specific heat generation Q/(rho*cp)
dimensionedScalar T_min(thermalProperties.lookup("T_min"));         // Minimum temperature of domain
