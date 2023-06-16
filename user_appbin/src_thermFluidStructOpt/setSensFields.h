//
// This header file sets fields related to adjoint sensitivity
//

Info << "Setting sensitivity fields\n" << endl;

volScalarField cost_sens_vol_frac
(
        IOobject
        (
                "cost_sens_vol_frac",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        gamma,
        zeroGradientFvPatchScalarField::typeName
);
volScalarField cost_sens_vol_frac0(gamma);

volScalarField cost_sens_power_loss
        (
                IOobject
                        (
                                "cost_sens_power_loss",
                                runTime.timeName(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                        ),
                -alpha_scale * alpha_U_max * (1 + q_ramp) * q_ramp / ((q_ramp + gamma) * (q_ramp + gamma)) * (U & U_adj_U),
                zeroGradientFvPatchScalarField::typeName
        );
volScalarField cost_sens_power_loss0(cost_sens_power_loss);

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
        - weight_sens_TU * alpha_scale * alpha_U_max * (1 + q_ramp) * q_ramp / ((q_ramp + gamma) * (q_ramp + gamma)) * (U & U_adj_T) * unit_correct_T / unit_correct
        + weight_sens_TK * fvc::laplacian((k_fluid - k_solid) * (1 + q_ramp) * q_ramp / ((q_ramp + gamma) * (q_ramp + gamma)), T_adj)/(rho_eff * cp_eff)
        + weight_sens_TQ * q_gen
        - weight_sens_TT * fvc::div(T_adj * U),
        zeroGradientFvPatchScalarField::typeName
);
volScalarField obj_sens_T0(obj_sens_T);

volScalarField cost_sens_comp
(
    IOobject
    (
        "cost_sens_comp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -gradD && (rho_solid * (3.0 * gamma * gamma * (Esp - Esp_min) / (2.0 * (1.0 + Po)) * twoSymm(gradD))),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField cost_sens_comp0(cost_sens_comp);
