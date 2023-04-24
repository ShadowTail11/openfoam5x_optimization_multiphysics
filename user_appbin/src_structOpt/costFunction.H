/*
 * This header calculates cost functions and convergence tracking variables.
 */

{
    //// Volume Fraction (used for monitoring, cost, & convergence)
    vol_frac = fvc::domainIntegrate(gamma).value() / area;

    Info << "\n--> Volume fraction is: " << vol_frac << endl;

    cost_vol_frac = vol_frac - set_vol_frac;
    vol_frac_conv = abs(cost_vol_frac / set_vol_frac);

    Info << "\n--> Percentage difference between volume fraction and target is: " << vol_frac_conv * 100 << "%" << endl;

//    //// Compliance (used for monitoring, cost, & convergence)
//    gradD = fvc::grad(D);
//    volSymmTensorField sigma(rho * (mu * twoSymm(gradD) + (lambda * I) * tr(gradD)));
//    volScalarField Energy((gradD + gradD.T()) && sigma);
//
//    compliance_prev = compliance_ratio;
//    compliance_ratio = fvc::domainIntegrate(Energy).value() / compliance_ref;
//
//    Info << "\n--> Compliance ratio is: " << compliance_ratio << endl;
//
//    compliance_conv = abs((compliance_ratio - compliance_prev) / (compliance_ratio + SMALL));
//
//    Info << "\n--> Percentage change in compliance ratio is: " << compliance_conv * 100 << "%\n" << endl;
}