//
// This function handles output from the optimization algorithm including
// monitoring values and mesh data (such as velocity at each node).
//

#ifndef THERMFLUIDOPT_OUTPUT_H
#define THERMFLUIDOPT_OUTPUT_H

// Output monitoring variables
if (Pstream::master())
{
    ofstream outfile1("Monitor_Time.txt", std::ios::app);
    outfile1 << runTime.elapsedClockTime() << "\n";
    outfile1.close();

    ofstream outfile2("Monitor_PseudoDensityDelta.txt", std::ios::app);
    outfile2 << std::sqrt(gamma_rms) << "\n";
    outfile2.close();

    ofstream outfile3("Monitor_AverageTemperature.txt", std::ios::app);
    outfile3 << T_ave << "\n";
    outfile3.close();

    ofstream outfile4("Monitor_TemperatureDrop.txt", std::ios::app);
    outfile4 << T_drop << "\n";
    outfile4.close();

    ofstream outfile5("Monitor_PowerDissipation.txt", std::ios::app);
    outfile5 << DissPower << "\n";
    outfile5.close();

    ofstream outfile6("Monitor_PressureDrop.txt", std::ios::app);
    outfile6 << pressure_drop << "\n";
    outfile6.close();

    ofstream outfile7("Monitor_VolumeFraction.txt", std::ios::app);
    outfile7 << fvc::domainIntegrate(gamma).value() / area << "\n";
    outfile7.close();
}

// Write mesh data
if (runTime.writeTime())
{
    gamma.write();
    T.write();
    U.write();
    p.write();
    nu_eff.write();
}

// Print elapsed time
Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
<< nl << endl;

// Exit if convergence has been achieved
if (dT_drop_ave < converge_tol && opt > 50) {
    gamma.write();
    T.write();
    U.write();
    p.write();
    nu_eff.write();
    Info << "Convergence criterion (<" << converge_tol * 100 << "%) has been met after "
    << opt - 1 << " iterations!" << endl << "Program ending!" << endl;

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;

    break;
}

#endif //THERMFLUIDOPT_OUTPUT_H
