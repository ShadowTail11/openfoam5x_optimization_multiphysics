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
volScalarField ramp_d(qd * (1 - gamma) / (qd + gamma));
volScalarField mu((ramp_d * (Esp - Esp_min) + Esp_min) / (2.0 * (1.0 + Po)));
volScalarField lambda((ramp_d * Po * (Esp - Esp_min) + Po * Esp_min) / ((1.0 + Po) * (1.0 - 2.0 * Po)));

Info << "Reading field D\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField gradD
(
    IOobject
    (
        "gradD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
      fvc::grad(D)
);

Info << "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mu * twoSymm(gradD) + (lambda * I) * tr(gradD)
);

Info << "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);

Switch compactNormalStress("yes");

if (compactNormalStress)
{
    divSigmaExp -= fvc::laplacian(2.0 * mu + lambda, D, "laplacian(DD,D)");
}
else
{
    divSigmaExp -= fvc::div((2.0 * mu + lambda) * fvc::grad(D), "div(sigmaD)");
}
mesh.setFluxRequired(D.name());

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
