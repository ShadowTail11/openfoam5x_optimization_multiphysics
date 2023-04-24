//
// This function handles output from the optimization algorithm including
// monitoring values and mesh data (such as velocity at each node).
//

{
    // Output monitoring variables
    if (Pstream::master())
    {
        ofstream outfile1("Monitor_Volume_Fraction.txt", std::ios::app);
        outfile1 << vol_frac << "\n";
        outfile1.close();

        ofstream outfile2("Monitor_Compliance_Ratio.txt", std::ios::app);
        outfile2 << compliance_ratio << "\n";
        outfile2.close();

        ofstream outfile3("Monitor_Time.txt", std::ios::app);
        outfile3 << runTime.elapsedClockTime() << "\n";
        outfile3.close();

        ofstream outfile4("Monitor_Iterations.txt", std::ios::app);
        outfile4 << opt << "\n";
        outfile4.close();
    }

    // Exit if convergence has been achieved
    if (opt > 50 && compliance_conv < merit_tol && vol_frac_conv < merit_tol)
    {
        Info << "Convergence criterion (<" << merit_tol * 100 << "%) has been met after "
        << opt - 1 << " iterations!" << endl << "Program ending!" << endl;

        // Write mesh data
        gamma.write();
        D.write();
        sigmaD.write();

        // Print elapsed time
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

        break;
    }

    // Write mesh data
    if (runTime.writeTime())
    {
        gamma.write();
        D.write();
        sigmaD.write();
    }

    // Print elapsed time
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;
}
