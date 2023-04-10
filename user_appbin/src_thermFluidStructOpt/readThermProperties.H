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
dimensionedScalar Q(thermalProperties.lookup("Q"));                 // Heat generation [W/m^3]
dimensionedScalar T_min(thermalProperties.lookup("T_min"));         // Minimum temperature of domain

// Compute effective parameters that transition between fluid and solid using the RAMP function
volScalarField rho_eff(rho_fluid + (rho_solid - rho_fluid) * ramp);                 // Effective density
volScalarField cp_eff(cp_fluid + (cp_solid - cp_fluid) * ramp);                     // Effective heat capacity
volScalarField k_eff(k_fluid + (k_solid - k_fluid) * ramp);                         // Effective thermal conductivity

// Specific heat generation within the solid region
dimensionedScalar q_gen = Q / (rho_solid * cp_solid);

// Effective thermal diffusivity
volScalarField alpha_T_eff
(
        IOobject
        (
                "alpha_T_eff",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        k_eff / (rho_eff * cp_eff)
);

Info << "Reading temperature field T\n" << endl;

volScalarField T(

        IOobject(

                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
        ),
        mesh
);

Info << "Reading adjoint temperature field T_adj\n" << endl;

volScalarField T_adj(

        IOobject(

                "T_adj",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        mesh
);

// Effective viscosity with optional temperature dependent capabilities
volScalarField nu_eff(

        IOobject(
                "nu_eff",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        ((1 - nu_temp_dep) * nu_const
         + max(nu_min,
               min(nu_max,
                   (nu_k - nu_slope * (T - nu_T_ref)) * pow(
                           max(
                                   dimensionedScalar("one", dimTime, 1.0) * std::sqrt(2.0)*mag(symm(fvc::grad(U))),
                                   dimensionedScalar("VSMALL", dimless, VSMALL)
                           ),
                           nu_n - scalar(1)
                   )
               )
        ) * nu_temp_dep)
);

volScalarField obj_sens_T
(
    IOobject
    (
        "obj_sens_T",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    - weight_sens_TU * alpha_U_max * (1 + q_ramp) * q_ramp / ((q_ramp + gamma) * (q_ramp + gamma)) * (U & U_adj_T) * unit_correct_T / unit_correct
    + weight_sens_TK * fvc::laplacian((k_fluid - k_solid) * (1 + q_ramp) * q_ramp / ((q_ramp + gamma) * (q_ramp + gamma)), T_adj)/(rho_eff * cp_eff)
    + weight_sens_TQ * q_gen
    - weight_sens_TT * fvc::div(T_adj * U),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField obj_sens_T0(obj_sens_T);
