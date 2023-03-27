//
// This header file manages the calculation and creation of fields.
//

#include "readOptProperties.h"      // Read user-defined optimization properties
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

// Fields related to primal and adjoint functions

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

labelList cells_test, cells_solid, cells_fluid;

if (test_area == 1)
{
    word zoneName1 = "zone_test";
    label zoneID1 = mesh.cellZones().findZoneID(zoneName1);
    cells_test = mesh.cellZones()[zoneID1];
    forAll(cells_test, i)
    {
        gamma_test[cells_test[i]] = 1;
    }
}

if (solid_area == 1)
{
    word zoneName2 = "zone_solid";
    label zoneID2 = mesh.cellZones().findZoneID(zoneName2);
    cells_solid = mesh.cellZones()[zoneID2];
    forAll(cells_solid, i)
    {
        gamma_test[cells_solid[i]] = 0;
        gamma[cells_solid[i]] = 0;
    }
}

if (fluid_area == 1)
{
    word zoneName3 = "zone_fluid";
    label zoneID3 = mesh.cellZones().findZoneID(zoneName3);
    cells_fluid = mesh.cellZones()[zoneID3];
    forAll(cells_fluid, i)
    {
        gamma_test[cells_fluid[i]] = 1;
        gamma[cells_fluid[i]] = 1;
    }
}

// Calculate flow resistance
volScalarField alpha_U(alpha_U_max * q_factor * (1 - gamma) / (q_factor + gamma));
volScalarField alpha_U_test(alpha_U_max * q_factor * (1 - gamma_test) / (q_factor + gamma_test));

// Primary adjoint sensitivity equation for flow optimization
volScalarField adjoint_sensitivity(

        IOobject(

                "adjoint_sensitivity",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        -alpha_U_max * (1 + q_factor) * q_factor / ((q_factor + gamma) * (q_factor + gamma)) * (U & U_adj_U),
        zeroGradientFvPatchScalarField::typeName
);

volScalarField adjoint_sensitivity_init(adjoint_sensitivity);  // Initialize the normalization of f_sens_T

// Secondary adjoint sensitivity equation for flow optimization
volScalarField cost_sensitivity
(
        IOobject
                (
                        "cost_sensitivity",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                ),
        gamma,
        zeroGradientFvPatchScalarField::typeName
);
volScalarField cost_sensitivity_init(gamma);

// Establish path reference for inlet and outlet
dictionary inlet_outlet = mesh.solutionDict().subDict("inlet_outlet");                  // Reference subDict in fvSolution
const int nObjPatch = inlet_outlet.lookupOrDefault<int>("numberConstraintPatches",2);   // Number of patches
wordList conPatchNames = inlet_outlet.lookup("constraintPatchesNames");                 // Function name
label conPatchList[nObjPatch];                                                          // Patch list
for (int i_loop = 0; i_loop < nObjPatch; i_loop++)
{
    conPatchList[i_loop] = mesh.boundaryMesh().findPatchID(conPatchNames[i_loop]);
}
