// Remove monitoring files from previous run
std::remove("Monitor_Time.txt");
std::remove("Monitor_PseudoDensityDelta.txt");
std::remove("Monitor_AverageTemperature.txt");
std::remove("Monitor_TemperatureDrop.txt");
std::remove("Monitor_PowerDissipation.txt");
std::remove("Monitor_PressureDrop.txt");
std::remove("Monitor_VolumeFraction.txt");

// Initialize variables
int opt(1);                                     // Optimization loop counter
int i, n_loops;                                 // Initialize loop counters (generic and solver)
scalar n_cells(gamma.size());                   // The number of cells in the domain
n_cells = static_cast<int>(n_cells);            // Force n_cells to be an integer

scalar area = gSum(mesh.V());                   // Domain area
scalar len = std::pow(mesh.V()[0], 1.0 / 3);    // Characteristic length of domain
scalar n_ave = 30;                              // Number of samples to average over in moving average method

scalar T_ave;                                   // Average temperature of domain
scalar T_drop(0);                               // Temperature drop between inlet and outlet
scalar T_drop_prev;                             // Previous temperature drop
scalar dT_drop(0);                              // Percentage change in temperature drop
scalar dT_drop_ave(dT_drop);                    // Moving average percentage change in temperature drop

scalar vol_frac;                                // Volume fraction of domain
scalar cost_vol_frac;                           // Cost variable for volume fraction
scalar vol_frac_conv;                           // Percentage difference between volume fraction and target

scalar P_drop(0);                               // Pressure drop between inlet and outlet
scalar P_drop_prev;                             // Previous pressure drop
scalar dP_drop(0);                              // Percentage change in pressure drop

scalar power_loss;                              // Power loss of the fluid
scalar power_loss_ratio;                        // Power loss ratio of the fluid
scalar cost_power_loss_ratio;                   // Cost variable for power loss ratio of the fluid
scalar power_loss_conv;                         // Percentage difference between power loss ratio and target

scalar gamma_switch(0);                         // Number of cells that switch between fluid and solid
scalar gamma_switch_prev;                       // Previous number of cells that switch between fluid and solid
scalar dgamma_switch(0);                        // Percentage change in number of cells that switch material
scalar dgamma_switch_ave(0);                    // Moving average of percent change in number of material switched cells

// Heaviside filter variables
double eta0, eta1, eta5 = 0, del = 0.1, y0, y1, y5;
b.value()= 1.0 / (r_filter * len / 3.464) / (r_filter * len / 3.464);

// PETSc variables

PetscInt m = 2;                 // Number of cost function constraints
PetscInt n_dof = n_cells;
PetscInt N_dof = n_cells;

PetscScalar Xmin = 0;           // Min. value of design variables
PetscScalar Xmax = 1;           // Max. value of design variables

// Initialize Method of Moving Asymptotes (MMA) bounds and communications
Vec xmma, xold, xmin, xmax;
PetscScalar *xmma_array;
VecCreate(PETSC_COMM_WORLD, &xmma);
VecSetSizes(xmma, n_dof, N_dof);
VecSetFromOptions(xmma);
VecGetArray(xmma, &xmma_array);

for(i = 0 ; i < n_cells; i++)
{
  xmma_array[i] = gamma[i];
}
VecRestoreArray(xmma, &xmma_array);

VecDuplicate(xmma, &xold);
VecDuplicate(xmma, &xmin);
VecDuplicate(xmma, &xmax);

PetscScalar aMMA[m];
PetscScalar cMMA[m];
PetscScalar dMMA[m];

for (PetscInt i = 0; i < m; i++)
{
  aMMA[i] = 0.0;
  dMMA[i] = 1.0;
  cMMA[i] = 1000.0;
}

MMA *mma;
mma = new MMA(N_dof, m, xmma, aMMA, cMMA, dMMA, mma_init, mma_dec, mma_inc);

// Sensitivities vectors of objective
Vec dfdx;
VecDuplicate(xmma, &dfdx);
PetscScalar *dfdx_array;
VecGetArray(dfdx, &dfdx_array);
VecRestoreArray(dfdx, &dfdx_array);

// Array with cost function constraint values
PetscScalar *gx = new PetscScalar[m];

// Sensitivities of cost function constraints
Vec *dgdx;
VecDuplicateVecs(xmma, m, &dgdx);
PetscScalar **dgdx_array;
VecGetArrays(dgdx, m, &dgdx_array);
VecRestoreArrays(dgdx, m, &dgdx_array);
