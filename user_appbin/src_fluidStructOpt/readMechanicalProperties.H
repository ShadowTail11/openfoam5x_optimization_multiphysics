//
// This header file reads user-defined properties from the mechanicalProperties file in the application's
// 'constant' directory
//

Info << "Reading mechanical properties\n" << endl;

IOdictionary mechanicalProperties
(
    IOobject
    (
        "mechanicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

int iCorr, nCorr = 150;                                                 // Initialize iterater and max iterations for primal structural solver
scalar convergenceTolerance(1e-5),initialResidual(0);                   // Initialize convergence for primal structural solver

Switch planeStress(mechanicalProperties.lookup("planeStress"));         // Read plane stress variable
dimensionedScalar E("E",
                    dimMass/dimLength/sqr(dimTime),
                    readScalar(mechanicalProperties.lookup("E")));      // Read Young's modulus
scalar Po(readScalar(mechanicalProperties.lookup("Po")));               // Read Poisson's ratio
dimensionedScalar Esp(E / rho_solid);                                   // Calculate specific modulus
dimensionedScalar Esp_min(Esp * 1e-9);                                  // Calculate minimum specific modulus threshold
double qd = 0.1;                                                        // Intialize RAMP factor for solid region
