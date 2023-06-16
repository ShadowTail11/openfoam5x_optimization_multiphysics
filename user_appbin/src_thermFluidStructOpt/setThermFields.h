//
// This header file sets fields related to thermal properties
//

Info << "Setting thermal fields\n" << endl;

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


// Maximum flow resistance
volScalarField alpha_U_max
(
        IOobject
        (
                "alpha_U_max",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        nu_eff / (l_char * l_char * darcy)
);

// Calculate flow resistance
volScalarField alpha_U(alpha_scale * alpha_U_max * ramp);
