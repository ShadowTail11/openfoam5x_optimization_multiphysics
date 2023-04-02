//
// This function is a pre-run that calculated a fluid power dissipation reference value that is used to
// constrain the power dissipation of the final design by means of a power dissipation ratio.
//

{
    if (test_area == 1)
    {
        // Solve for primal fields in test region
        #include "testPrimalFlowSolver.h"
        #include "testPrimalThermalSolver.h"

        // Calculate the temperature drop reference value based on the test area
        scalar T_drop_ref;
        T_drop_ref = (sum(phi_test.boundaryField()[conPatchList[1]] * T_test.boundaryField()[conPatchList[1]])
                / sum(phi_test.boundaryField()[conPatchList[1]])
                -sum(phi_test.boundaryField()[conPatchList[0]] * T_test.boundaryField()[conPatchList[0]])
                / sum(phi_test.boundaryField()[conPatchList[0]]));

        Info << "\n--> Temperature drop for between inlet and outlet for the test zone is: " << T_drop_ref << endl;

        // Calculate the pressure drop reference value based on the test area
        scalar P_drop_ref;
        P_drop_ref = (sum(phi_test.boundaryField()[conPatchList[1]] * p_test.boundaryField()[conPatchList[1]])
                / sum(phi_test.boundaryField()[conPatchList[1]])
                -sum(phi_test.boundaryField()[conPatchList[0]] * p_test.boundaryField()[conPatchList[0]])
                / sum(phi_test.boundaryField()[conPatchList[0]])) * rho_fluid.value();

        Info << "\n--> Pressure drop between inlet and outlet for the test zone is: " << P_drop_ref << endl;

        Info << "\nInitial power_loss_ref = " << power_loss_ref << nl << endl;

        // Calculate the power loss reference value based on the test area
        scalar power_loss_ref(0);
        for(i = 0; i < nObjPatch; i++)
        {
            power_loss_ref = power_loss_ref - sum(
                    phi_test.boundaryField()[conPatchList[i]]
                    * (p_test.boundaryField()[conPatchList[i]] + 0.5 * magSqr(U_test.boundaryField()[conPatchList[i]])));
        }

        if(geo_dim == 2)
        {
            power_loss_ref = power_loss_ref / len;
        }

        Info << "Updated power_loss_ref = " << power_loss_ref << nl << endl;

        // Turn off test
        test_area = 0;

        // Write parameters to
        U_test.write();
        p_test.write();
        T_test.write();
        gamma_test.write();

        // Output reference values
        ofstream outfile8("Reference_Test_Zone.txt", std::ios::app);
        outfile8 << "Temperature Drop Reference = " << T_drop_ref << "\n";
        outfile8 << "Pressure Drop Reference = " << P_drop_ref << "\n";
        outfile8 << "Power Loss Reference = " << power_loss_ref << "\n";
        outfile8.close();
    }
}
