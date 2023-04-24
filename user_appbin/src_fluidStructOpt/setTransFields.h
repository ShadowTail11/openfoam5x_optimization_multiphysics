//
// This header file sets fields related to the transport properties
//

Info << "Setting transport fields\n" << endl;

// Compute effective parameters that transition between fluid and solid using the RAMP function
volScalarField rho_eff(rho_fluid + (rho_solid - rho_fluid) * ramp);                 // Effective density

// Calculate flow resistance
volScalarField alpha_U(alpha_U_max * ramp);

#include "createMRF.H"

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

volScalarField nu_const = turbulence->nu();

// Establish path reference for inlet and outlet
dictionary inlet_outlet = mesh.solutionDict().subDict("inlet_outlet");                  // Reference subDict in fvSolution
const int nObjPatch = inlet_outlet.lookupOrDefault<int>("numberConstraintPatches",2);   // Number of patches
wordList conPatchNames = inlet_outlet.lookup("constraintPatchesNames");                 // Function name
label conPatchList[nObjPatch];                                                          // Patch list
for (int i_loop = 0; i_loop < nObjPatch; i_loop++)
{
conPatchList[i_loop] = mesh.boundaryMesh().findPatchID(conPatchNames[i_loop]);
}

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
                -alpha_U_max * (1 + q_ramp) * q_ramp / ((q_ramp + gamma) * (q_ramp + gamma)) * (U & U_adj_U),
                zeroGradientFvPatchScalarField::typeName
        );
volScalarField cost_sens_power_loss0(cost_sens_power_loss);
