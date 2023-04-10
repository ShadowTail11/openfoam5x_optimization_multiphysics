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

        ofstream outfile2("Monitor_AverageTemperature.txt", std::ios::app);
        outfile2 << T_ave << "\n";
        outfile2.close();

        ofstream outfile3("Monitor_TemperatureDrop.txt", std::ios::app);
        outfile3 << T_drop << "\n";
        outfile3.close();

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
        outfile7 << C << "\n";
        outfile7.close();

        ofstream outfile8("Monitor_PseudoDensityDelta.txt", std::ios::app);
        outfile8 << gamma_switch << "\n";
        outfile8.close();
    }

    // Write mesh data
    if (runTime.writeTime())
    {
        // Eventual output (commented out until changes can be finalized)
        gamma.write();
        T.write();
        U.write();
        p.write();
        D.write();
        sigmaD.write();
        nu_eff.write();

        obj_sens_T.write();
        cost_sens_power_loss.write();
        cost_sens_comp.write();
        cost_sens_vol_frac.write();
    }

    // Print elapsed time
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;

    // Exit if convergence has been achieved
    if (opt > 50 && dgamma_switch_ave < gamma_tol && power_loss_conv < merit_tol && dT_drop_ave < merit_tol)
    {
        // Round pseudo density so each node is either solid (gamma=0) or fluid (gamma=1) rather than porous (0<gamma<1)
        forAll(gamma, i)
        {
            gamma[i] = std::round(gamma[i]);
        }

        gamma.write();
        T.write();
        U.write();
        p.write();
        D.write();
        sigmaD.write();
        nu_eff.write();

        obj_sens_T.write();
        cost_sens_power_loss.write();
        cost_sens_comp.write();
        cost_sens_vol_frac.write();

        Info << "Convergence criterion (<" << merit_tol * 100 << "%) has been met after "
        << opt - 1 << " iterations!" << endl << "Program ending!" << endl;

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

        break;
    }

    opt++;
}