//******************************************************************//pri-fluid

#include "createMRF.H"

Info << "Reading optimization properties\n" << endl;

IOdictionary optimizationProperties(

        IOobject(

                "optimizationProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
        )
);

// Domain definition
scalar solid_area(readScalar(optimizationProperties.lookup("solid_area")));             // Is solid region set? (1=yes)
scalar fluid_area(readScalar(optimizationProperties.lookup("fluid_area")));             // Is fluid region used? (1=yes)
scalar test_area(readScalar(optimizationProperties.lookup("test_area")));               // Is test region used? (1=yes)
scalar geo_dim(readScalar(optimizationProperties.lookup("geo_dim")));                   // Dimensions of solver (2D vs. 3D)

// Optimization properties
scalar mma_init(readScalar(optimizationProperties.lookup("mma_init")));                 // Initial optimization increment
scalar mma_dec(readScalar(optimizationProperties.lookup("mma_dec")));                   // Optimization decrement
scalar mma_inc(readScalar(optimizationProperties.lookup("mma_inc")));                   // Optimization increment
scalar mma_limit(readScalar(optimizationProperties.lookup("mma_limit")));               // Optimization limit
scalar raa0(readScalar(optimizationProperties.lookup("raa0")));                         // MMA smoothing parameter
scalar r_filter_scalar(readScalar(optimizationProperties.lookup("r_filter")));          // Density filter radius
int r_filter = r_filter_scalar;
dimensionedScalar b("b", dimensionSet(0,-2,0,0,0,0,0),1.0);                             // Filter radius of the PDE filter

// Cost function constraints
scalar set_vol_frac(readScalar(optimizationProperties.lookup("set_vol_frac")));         // Target volume fraction
scalar set_power_ratio(readScalar(optimizationProperties.lookup("set_power_ratio")));   // Target power loss ratio
scalar power_loss_ref(readScalar(optimizationProperties.lookup("power_loss_ref")));     // Power loss reference

// Optimization weights
scalar weight_TU(readScalar(optimizationProperties.lookup("weight_TU")));               // Sensitivity weight for thermal velocity term
scalar weight_TK(readScalar(optimizationProperties.lookup("weight_TK")));               // Sensitivity weight for thermal conductivity term
scalar weight_TQ(readScalar(optimizationProperties.lookup("weight_TQ")));               // Sensitivity weight for heat generation term
scalar weight_TT(readScalar(optimizationProperties.lookup("weight_TT")));               // Sensitivity weight for thermal diffusivity term
scalar weight_T(readScalar(optimizationProperties.lookup("weight_T")));                 // Sensitivity weight for thermal term
scalar weight_U(readScalar(optimizationProperties.lookup("weight_U")));                 // Sensitivity weight for velocity term

// Adjoint unit correctors (used because adjoint variables have inverse units of parent, but OpenFOAM uses same units)
dimensionedScalar unit_correct("unit_correct", dimensionSet(0,2,-2,-2,0,0,0),1.0);      // Adjoint (U/T)^2 unit correction
dimensionedScalar unit_correct_U("unit_correct_U", dimensionSet(0,-1,1,0,0,0,0),1.0);   // Adjoint U unit correction
dimensionedScalar unit_correct_T("unit_correct_T", dimensionSet(0,0,0,-1,0,0,0),1.0);   // Adjoint T unit correction

// Convergence controls
scalar merit_tol(readScalar(optimizationProperties.lookup("merit_tol")));               // Convergence tolerance for merit parameters
scalar gamma_tol(readScalar(optimizationProperties.lookup("gamma_tol")));               // Convergence tolerance for pseudo density gamma

Info << "Reading thermal properties\n" << endl;

// Access thermalProperties file to extract relevant properties
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

Info << "Reading transport properties\n" << endl;

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

volScalarField x(gamma_densityfilter);
volScalarField drho(x);

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

//***************************************************************************//pri-heat

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

//************************************************************************//

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

//************************************************************************//

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

//************************************************************************//

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

//************************************************************************//


labelList cells1, cells_solid, cells_fluid;

// Set fixed regions
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

// Create the parallel run pseudo density variable Allgamma

//List<scalarField> AllrList(Pstream::nProcs());
//scalarField localr(gamma.internalField());
//AllrList[Pstream::myProcNo()] = localr;
//Pstream::gatherList(AllrList);
//Pstream::scatterList(AllrList);
//
//scalarField Allgamma(
//
//    ListListOps::combine<scalarField>(
//
//        AllrList,
//        accessOp<scalarField>()
//    )
//);

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
