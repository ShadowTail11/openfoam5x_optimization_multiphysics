//
// This function handles output from the optimization algorithm including
// monitoring values and mesh data (such as velocity at each node).
//

{
    // Output monitoring variables
    if (Pstream::master())
    {
        ofstream outfile1("Monitor_Time.txt", std::ios::app);
        outfile1 << runTime.elapsedClockTime() << "\n";
        outfile1.close();

        ofstream outfile4("Monitor_PowerLoss.txt", std::ios::app);
        outfile4 << power_loss_ratio << "\n";
        outfile4.close();

        ofstream outfile5("Monitor_PressureDrop.txt", std::ios::app);
        outfile5 << P_drop << "\n";
        outfile5.close();

        ofstream outfile6("Monitor_VolumeFraction.txt", std::ios::app);
        outfile6 << vol_frac << "\n";
        outfile6.close();

        ofstream outfile7("Monitor_Compliance.txt",std::ios::app);
        outfile7 << compliance_ratio << "\n";
        outfile7.close();

        ofstream outfile8("Monitor_MaxStress.txt",std::ios::app);
        outfile8 << gMax(mag(sigmaD.primitiveField())) << "\n";
        outfile8.close();

        ofstream outfile9("Monitor_PseudoDensityDelta.txt", std::ios::app);
        outfile9 << gamma_switch << "\n";
        outfile9.close();

        ofstream outfile10("Monitor_Iterations.txt", std::ios::app);
        outfile10 << opt << "\n";
        outfile10.close();
    }

    // Exit if convergence has been achieved
    if (opt > 0.75 * runTime.endTime().value() || (opt > 50 && dgamma_switch_ave < gamma_tol && power_loss_conv < merit_tol) )
//    if (opt == 10)
    {
        // Create pseudo density field masks to freeze the fluid region and create an initial fluid wall
        forAll(gamma, i)
        {
            gamma_conv[i] = (
                    std::round(Foam::min(1.0, gamma[i] / 0.2))
                    * std::round(Foam::min(1.0, mag(U.primitiveField()[i]) / gMax(mag(U.boundaryField()[conPatchList[0]])) / 0.2))
                    );
            gamma_wall[i] = (
                    std::round(Foam::min(1.0, gamma[i] / 0.2))
                    * std::round(Foam::min(1.0, mag(U.primitiveField()[i]) / gMax(mag(U.boundaryField()[conPatchList[0]])) / 0.01))
                    - gamma_conv[i]);
        }

        // Calculate the frozen region
        vol_frac_frozen = fvc::domainIntegrate(gamma_conv).value() / area;
        n_frozen = vol_frac_frozen * n_cells;

        // Calculate the frozen region with gamma factored in
        scalar vol_frac_frozen_gamma = fvc::domainIntegrate(gamma_conv * gamma).value() / area;

        // Calculate the solid region
        scalar vol_frac_solid = fvc::domainIntegrate(vol_set_solid).value() / area;

        Info << "The number of cells that are frozen are: " << n_frozen << endl
        << "Which represents: " << vol_frac_frozen * 100 << "% of the total design space." << endl;

        set_vol_frac = set_vol_frac_solid * (vol_frac_frozen + vol_frac_solid) + vol_frac_frozen_gamma;

        gamma = (gamma * gamma_conv + (1.0 - set_vol_frac_solid) * (1.0 - gamma_conv)) * (1.0 - gamma_wall);
//        gamma = gamma * gamma_conv + (1.0 - set_vol_frac_solid) * (1.0 - gamma_conv);
        x = gamma;
        if (solid_area)
        {
            forAll(cells_solid, i)
            {
                x[cells_solid[i]] = 0;
            }
        }

        gamma.write();
        gamma_conv.write();
        gamma_wall.write();
        U.write();
        p.write();
        D.write();
        sigmaD.write();

        cost_sens_power_loss.write();
        cost_sens_comp.write();
        cost_sens_vol_frac.write();

        // Remaining outputs used to initialize continued run
        p_adj_U.write();
        U_adj_U.write();

        Info << "Convergence criterion (<" << merit_tol * 100 << "%) has been met after "
        << opt << " iterations!" << endl << "Thermal-Fluid solver ending!" << endl;

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

        solve_ext = 1;
        p_operating_field = p_operating_field + (P_operating - P_ref) / rho_fluid.value();
        opt++;
    }

    // Write mesh data
    if (runTime.writeTime())
    {
        // Eventual output (commented out until changes can be finalized)
        gamma.write();
        U.write();
        p.write();
        D.write();
        sigmaD.write();

        cost_sens_power_loss.write();
        cost_sens_comp.write();
        cost_sens_vol_frac.write();

        // Remaining outputs used to initialize continued run

        p_adj_U.write();
        U_adj_U.write();
    }

    // Print elapsed time
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;

    opt++;
}