{
    Info<< "Calculating the Cost Function" << endl;
    meanT_prev = meanT;
//    V=0;
//    meanT=0;
//
//    for(i = 0; i < allcells; i++)
//    {
//      V += gamma[i];
//      meanT += T[i] * gamma[i];
//    }

    reduce(meanT, sumOp<scalar>());
    reduce(V, sumOp<scalar>());
    meanT = fvc::domainIntegrate(T * gamma).value() / (allcells * set_vol_frac);
    V = fvc::domainIntegrate(gamma).value() / area - set_vol_frac;

    dmeanT = (meanT - meanT_prev) / (nallcells * set_vol_frac);
    if(opt == 1)
    {
        dmeanT = 0;
    }
    dmeanT_sum = dmeanT_sum * converge_smooth + dmeanT * dmeanT;
    dmeanT_RMS = std::sqrt(dmeanT_sum / opt) / std::abs(meanT / (nallcells * set_vol_frac) + 0.0000000005);

//DissPower
    DissPower = 0;
    for(i = 0; i < nObjPatch; i++)
    {
        DissPower = DissPower - sum(
        phi.boundaryField()[conPatchList[i]] * (p.boundaryField()[conPatchList[i]]
        + 0.5 * magSqr(U.boundaryField()[conPatchList[i]])));
    }

    DissPower = DissPower / power_loss_ref;

    if(geo_dim == 2)
    {
      DissPower = DissPower / len;
    }

    reduce(DissPower, sumOp<scalar>());
    Info << "Constraint function (DissPower): " << DissPower << endl;

    pressure_drop = (sum(phi.boundaryField()[conPatchList[1]] * p.boundaryField()[conPatchList[1]])
            / sum(phi.boundaryField()[conPatchList[1]])
            -sum(phi.boundaryField()[conPatchList[0]] * p.boundaryField()[conPatchList[0]])
            / sum(phi.boundaryField()[conPatchList[0]])) * rho_fluid.value();

    temperature_drop = (sum(phi.boundaryField()[conPatchList[1]] * T.boundaryField()[conPatchList[1]])
            / sum(phi.boundaryField()[conPatchList[1]])
            -sum(phi.boundaryField()[conPatchList[0]] * T.boundaryField()[conPatchList[0]])
            / sum(phi.boundaryField()[conPatchList[0]]));

}
