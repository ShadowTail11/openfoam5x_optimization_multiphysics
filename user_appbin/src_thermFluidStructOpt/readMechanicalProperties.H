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

int iCorr, nCorr = 150;
scalar convergenceTolerance(1e-5),initialResidual(0);

Switch planeStress(mechanicalProperties.lookup("planeStress"));
dimensionedScalar E("E", dimMass/dimLength/sqr(dimTime),readScalar(mechanicalProperties.lookup("E")));
scalar Po(readScalar(mechanicalProperties.lookup("Po")));
dimensionedScalar Esp(E / rho_solid);
dimensionedScalar Esp_min(Esp * 1e-9);
double qd = 0.1;
//ramp_d = qd * (1 - gamma) / (qd + gamma);
