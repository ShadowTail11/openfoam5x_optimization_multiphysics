/*
 * This function updates properties that change based on the results of the primal solve
 * such as temperature dependent variables so that they can be used in the adjoint solvers.
 */

{
    // Temperature dependent viscosity (with an override option to make it constant)
    nu_eff = ((1 - nu_temp_dep) * nu_const
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
                                      ) * nu_temp_dep);
}
