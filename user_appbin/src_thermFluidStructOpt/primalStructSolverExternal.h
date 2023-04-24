Info << "\nSolving for the linear elasticity equation\n" << endl;

initialResidual = 0;
iCorr = 0;
do
{
    {
        fvVectorMatrix D_extEqn(

                fvm::d2dt2(D_ext)
                ==
                fvm::laplacian(2 * mu + lambda, D_ext, "laplacian(DD,D_ext)")
                + divSigmaExp_ext
        //              - fvc::grad(threeKalpha*T)
        );

        initialResidual = D_extEqn.solve().max().initialResidual();

        if (!compactNormalStress)
        {
            divSigmaExp_ext = fvc::div(D_extEqn.flux());
        }
    }

    {
        volTensorField gradD_ext(fvc::grad(D_ext));
        sigmaD_ext = mu * twoSymm(gradD_ext) + (lambda * I) * tr(gradD_ext);
//        sigmaD = sigmaD_ext;

        if (compactNormalStress)
        {
            divSigmaExp_ext = fvc::div
            (
                    sigmaD_ext - (2 * mu + lambda) * gradD_ext,
                    "div(sigmaD_ext)"
            );
        }
        else
        {
            divSigmaExp_ext += fvc::div(sigmaD_ext);
        }
    }

} while (initialResidual > convergenceTolerance && ++iCorr < nCorr);
