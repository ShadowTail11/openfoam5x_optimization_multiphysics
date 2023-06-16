//
// This header file sets fields related to adjoint sensitivity
//

Info << "Setting sensitivity fields\n" << endl;

volScalarField g_sens_vol_frac
(
        IOobject
        (
                "g_sens_vol_frac",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        gamma,
        zeroGradientFvPatchScalarField::typeName
);
volScalarField g_sens_vol_frac_norm(gamma);

// Secondary adjoint sensitivity equation for flow optimization
volScalarField g_sens_power_loss(

        IOobject(
                "g_sens_power_loss",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        -alpha_scale * alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_U),
        zeroGradientFvPatchScalarField::typeName
);

volScalarField g_sens_power_loss_init(g_sens_power_loss); // Initialize the normalization of g_sense_power_loss

// Primary adjoint sensitivity equation for thermal optimization (with broken up into subcomponents and whole term)
volScalarField f_sens_TU(

        IOobject(

                "f_sens_TU",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        - weight_sens_TU * alpha_scale * alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_T) * unit_correct_T / unit_correct
);

volScalarField f_sens_TK(

        IOobject(

                "f_sens_TK",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        - weight_sens_TK * fvc::laplacian((k_fluid-k_solid)*(1+q_factor)*q_factor/((q_factor+gamma)*(q_factor+gamma)), T_adj)/(rho_eff*cp_eff)
);

volScalarField f_sens_TT(

        IOobject(

                "f_sens_TT",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        - weight_sens_TT * fvc::div(T_adj * U)
);

volScalarField f_sens_T(

        IOobject(

                "f_sens_T",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        - weight_sens_TU * alpha_scale * alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_T) * unit_correct_T / unit_correct
        + weight_sens_TK * fvc::laplacian((k_fluid-k_solid)*(1+q_factor)*q_factor/((q_factor+gamma)*(q_factor+gamma)), T_adj)/(rho_eff*cp_eff)
        + weight_sens_TQ * q_gen
        - weight_sens_TT * fvc::div(T_adj * U),
        zeroGradientFvPatchScalarField::typeName
);

volScalarField f_sens_T_init(f_sens_T);  // Initialize the normalization of f_sens_T
