//
// This header file manages the calculation and creation of fields.
//

#include "readOptProperties.h"      // Read user-defined optimization properties
#include "readThermProperties.h"    // Read user-defined thermal properties
#include "readTransProperties.h"    // Read user-defined transport properties

#include "createMRF.H"

// Pseudo density gamma
Info << "Reading pseudo density field gamma\n" << endl;

volScalarField gamma(

        IOobject(

                "gamma",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        mesh,
        scalar(set_vol_frac),
        zeroGradientFvPatchScalarField::typeName
);

// Field used to facilitate density filtering for the pseudo density after optimization
volScalarField gamma_densityfilter
        (
                IOobject
                        (
                                "gamma_densityfilter",
                                runTime.timeName(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                        ),
                mesh,
                scalar(set_vol_frac),
                zeroGradientFvPatchScalarField::typeName
        );

volScalarField x(gamma_densityfilter);      // Field used to facilitate Heaviside sensitivity calculation
volScalarField drho(x);                     // Density factor used in Heaviside function

// Pseudo density field initialized to get a reference power loss value
volScalarField gamma_test(

        IOobject(

                "gamma_test",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        mesh,
        scalar(0),
        zeroGradientFvPatchScalarField::typeName
);

// Effective density
volScalarField rho_eff(
        rho_solid + (rho_fluid - rho_solid) * gamma * (1 + q_factor) / (q_factor + gamma)
);

// Effective heat capacity
volScalarField cp_eff(
        cp_solid + (cp_fluid - cp_solid) * gamma * (1 + q_factor) / (q_factor + gamma)
);

// Effective thermal conductivity
volScalarField k_eff(
        k_solid + (k_fluid - k_solid) * gamma * (1 + q_factor) / (q_factor + gamma)
);

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
        (k_solid + (k_fluid - k_solid) * gamma * (1 + q_factor) / (q_factor + gamma)) / (rho_eff * cp_eff)
);

// Fields related to primal and adjoint functions

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

Info << "Reading kinematic pressure field p\n" << endl;

volScalarField p(

    IOobject(

        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading velocity field U\n" << endl;

volVectorField U(

    IOobject(

        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

// Sets reference point for pressure
label pRefCell = 0;
scalar pRefValue = 0.0;

setRefCell(

        p,
        simple.dict(),
        pRefCell,
        pRefValue
);

// Set turbulence model and create pointer to model
singlePhaseTransportModel laminarTransport(U, phi);
autoPtr<incompressible::turbulenceModel> turbulence(

    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

Info << "Reading thermal adjoint pressure field p_adj_T\n" << endl;

volScalarField p_adj_T(

    IOobject(

        "p_adj_T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field U_adj_T\n" << endl;

volVectorField U_adj_T(

    IOobject(

        "U_adj_T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading/calculating adjoint potential velocity field phi_adj_T\n" << endl;

surfaceScalarField phi_adj_T(

    IOobject(

        "phi_adj_T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U_adj_T) & mesh.Sf()
);

// Set adjoint pressure reference
label p_adj_TRefCell = 0;
scalar p_adj_TRefValue = 0.0;

setRefCell(

    p_adj_T,
    simple.dict(),
    p_adj_TRefCell,
    p_adj_TRefValue
);

mesh.setFluxRequired(p_adj_T.name());

Info << "Reading adjoint kinematic pressure (flow) field p_adj_U\n" << endl;

volScalarField p_adj_U(

        IOobject(

                "p_adj_U",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        mesh
);

Info << "Reading adjoint velocity field U_adj_U\n" << endl;

volVectorField U_adj_U(

        IOobject(

                "U_adj_U",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        mesh
);

Info << "Reading/calculating face flux field phi_adj_U\n" << endl;

surfaceScalarField phi_adj_U(

        IOobject(

                "phi_adj_U",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        linearInterpolate(U_adj_U) & mesh.Sf()
);

// Set adjoint kinematic pressure reference (flow)
label p_adj_URefCell = 0;
scalar p_adj_URefValue = 0.0;

setRefCell(

        p_adj_U,
        simple.dict(),
        p_adj_URefCell,
        p_adj_URefValue
);

mesh.setFluxRequired(p_adj_U.name());

// Primal fields used to get power loss reference

Info << "Reading test pressure field p_test\n" << endl;

volScalarField p_test(

        IOobject(

                "p_test",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        mesh
);

Info << "Reading test velocity field U_test\n" << endl;

volVectorField U_test(

        IOobject(

                "U_test",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        mesh
);

Info << "Reading/calculating test face flux field phi_test\n" << endl;

surfaceScalarField phi_test(

        IOobject(

                "phi_test",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        linearInterpolate(U_test) & mesh.Sf()
);

// Set adjoint kinematic pressure reference (flow)
label p_testRefCell = 0;
scalar p_testRefValue = 0.0;

setRefCell(

        p_test,
        simple.dict(),
        p_testRefCell,
        p_testRefValue
);

mesh.setFluxRequired(p_test.name());

// Set pseudo density values in fixed material regions

labelList cells1, cells_solid, cells_fluid;

if (test_area == 1)
{
    word zoneName1 = "zone_test";
    label zoneID1 = mesh.cellZones().findZoneID(zoneName1);
    cells1 = mesh.cellZones()[zoneID1];
    setCells(gamma_test, cells1, 1);
}

if (solid_area == 1)
{
    word zoneName2 = "zone_solid";
    label zoneID2 = mesh.cellZones().findZoneID(zoneName2);
    cells_solid = mesh.cellZones()[zoneID2];
    setCells(gamma, cells_solid, 0);
    setCells(gamma_test, cells_solid, 0);
}

if (fluid_area == 1)
{
    word zoneName3 = "zone_fluid";
    label zoneID3 = mesh.cellZones().findZoneID(zoneName3);
    cells_fluid = mesh.cellZones()[zoneID3];
    setCells(gamma, cells_fluid, 1);
    setCells(gamma_test, cells_fluid, 1);
}

// Calculate flow resistance
volScalarField alpha_U(alpha_U_max * q_factor * (1 - gamma) / (q_factor + gamma));
volScalarField alpha_U_test(alpha_U_max * q_factor * (1 - gamma_test) / (q_factor + gamma_test));

// Primary adjoint sensitivity equation for thermal optimization
volScalarField f_sens_TU(

        IOobject(

                "f_sens_TU",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        - weight_TU * alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_T) * unit_correct_T / unit_correct
);

volScalarField f_sens_TK(

        IOobject(

                "f_sens_TK",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        - weight_TK * fvc::laplacian((k_fluid-k_solid)*(1+q_factor)*q_factor/((q_factor+gamma)*(q_factor+gamma)), T_adj)/(rho_eff*cp_eff)
);

volScalarField f_sens_TT(

        IOobject(

                "f_sens_TT",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
        - weight_TT * fvc::div(T_adj * U)
);

volScalarField f_sens_T(

        IOobject(

                "f_sens_T",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        - weight_TU * alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_T) * unit_correct_T / unit_correct
        + weight_TK * fvc::laplacian((k_fluid-k_solid)*(1+q_factor)*q_factor/((q_factor+gamma)*(q_factor+gamma)), T_adj)/(rho_eff*cp_eff)
        + weight_TQ * q_gen
        - weight_TT * fvc::div(T_adj * U),
        zeroGradientFvPatchScalarField::typeName
);

volScalarField f_sens_T_init(f_sens_T);  // Initialize the normalization of f_sens_T

// Secondary adjoint sensitivity equation for flow optimization
volScalarField g_sens_power_loss(

        IOobject(
                "g_sens_power_loss",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        -alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_U),
        zeroGradientFvPatchScalarField::typeName
);

volScalarField g_sens_power_loss_init(g_sens_power_loss); // Initialize the normalization of g_sense_power_loss

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

volScalarField nu_const = turbulence->nu();

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

// Establish path reference for inlet and outlet
dictionary Fun1 = mesh.solutionDict().subDict("power_dissipation");
const int nObjPatch = Fun1.lookupOrDefault<scalar>("numberConstraintPatches",2);    // Number of patches
wordList conPatchNames = Fun1.lookup("constraintPatchesNames");                     // Function name
label conPatchList[nObjPatch];                                                      // Patch list
int i_loop;
for (i_loop = 0; i_loop < nObjPatch; i_loop++)
{
    conPatchList[i_loop] = mesh.boundaryMesh().findPatchID(conPatchNames[i_loop]);
}
