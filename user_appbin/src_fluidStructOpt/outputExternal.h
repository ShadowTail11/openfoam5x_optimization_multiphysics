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
    if (opt > 50 && dgamma_switch < gamma_tol / 100 && compliance_conv < merit_tol / 100 && vol_frac_conv < merit_tol / 100)
//    if (opt == 300)
    {
        // Round pseudo density so each node is either solid (gamma=0) or fluid (gamma=1) rather than porous (0<gamma<1)
        forAll(gamma, i)
        {
            gamma_conv[i] = (
                    std::round(gamma[i] / 0.2)
                    * std::round(Foam::min(1.0, mag(U.primitiveField()[i]) / gMax(mag(U.primitiveField())) / 0.02))
                    );
        }

        gamma.write();
        gamma_conv.write();
        U.write();
        p.write();
        D.write();
        sigmaD.write();

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

    // Write mesh data
    if (runTime.writeTime())
    {
        gamma.write();
        U.write();
        p.write();
        D.write();
        sigmaD.write();

        cost_sens_power_loss.write();
        cost_sens_comp.write();
        cost_sens_vol_frac.write();
    }

    // Print elapsed time
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;

    opt++;
}