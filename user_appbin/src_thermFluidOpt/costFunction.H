{
    Info<< "\nCalculating the monitoring variables" << endl;

    // Temperature
    T_ave = fvc::domainIntegrate(T * gamma).value() / (area * set_vol_frac);

    Info << "\n--> Average fluid temperature is: " << T_ave << endl;

    T_drop_prev = T_drop;
    T_drop = (sum(phi.boundaryField()[conPatchList[1]] * T.boundaryField()[conPatchList[1]])
            / sum(phi.boundaryField()[conPatchList[1]])
            -sum(phi.boundaryField()[conPatchList[0]] * T.boundaryField()[conPatchList[0]])
            / sum(phi.boundaryField()[conPatchList[0]]));

    Info << "\n--> Temperature drop between inlet and outlet is: " << T_drop << endl;

    dT_drop = (T_drop - T_drop_prev) / T_drop;
    dT_drop_ave = (1 - 1 / n_ave) * dT_drop_ave + dT_drop / n_ave;

    Info << "\n--> Percentage change in temperature drop is: " << dT_drop * 100 << "%" << endl;

    // Volume Fraction
    reduce(V, sumOp<scalar>());
    V = fvc::domainIntegrate(gamma).value() / area - set_vol_frac;

    Info<< "\nCalculating the Cost Function" << endl;
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

}
