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
int i, n_loops;                                 // Initialize loop counters (generic and solver, respectively)
scalar n_cells(gamma.size());                   // The number of cells in the domain
n_cells = static_cast<int>(n_cells);            // Force n_cells to be an integer
int opt_threshold = -1;                         // Updated to iteration count when merit value is within 2x of tolerance

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
scalar vol_frac_conv(0);                        // Percentage difference between volume fraction and target

scalar P_drop(0);                               // Pressure drop between inlet and outlet
scalar P_drop_prev;                             // Previous pressure drop
scalar dP_drop(0);                              // Percentage change in pressure drop

scalar power_loss;                              // Power loss of the fluid
scalar power_loss_ratio;                        // Power loss ratio of the fluid
scalar cost_power_loss_ratio;                   // Cost variable for power loss ratio of the fluid
scalar power_loss_conv(0);                      // Percentage difference between power loss ratio and target

scalar gamma_switch(0);                         // Number of cells that switch between fluid and solid
scalar gamma_switch_prev;                       // Previous number of cells that switch between fluid and solid
scalar dgamma_switch(0);                        // Percentage change in number of cells that switch material
scalar dgamma_switch_ave(0);                    // Moving average of percent change in number of material switched cells

// Heaviside filter variables
double eta0, eta1, eta5 = 0, del = 0.1, y0, y5;
b.value()= 1.0 / (r_filter * len / 3.464) / (r_filter * len / 3.464);

// Initialize Method of Moving Asymptotes (MMA)
int m_constraints = 2;                                      // Number of constraints
MMA mma(n_cells, m_constraints);                            // Set MMA optimization size
mma.asyminit = mma_init;                                    // Set initial movement
mma.asymdec = mma_dec;                                      // Set decrement size
mma.asyminc = mma_inc;                                      // Set increment size
mma.raa0 = raa0;                                            // Set smoothing factor
mma.movelimit = mma_limit;                                  // Set movement limit
std::vector<double> xmma(n_cells);                          // Vector that stores optimized gamma
std::vector<double> dfdx(n_cells);                          // Vector that stores sensitivity function
std::vector<double> gx(m_constraints);                      // Vector that stores constraint values
std::vector<std::vector<double>> dgdx(m_constraints);       // Vector that stores constrain functions
for (int i = 0; i < m_constraints; i++)                     // Resize constraint vector to size of domain
{
    dgdx[i].resize(n_cells);
}
