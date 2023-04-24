//
// This header file initialized parameters not initialized previously.
//

// Remove monitoring files from previous run if starting a fresh run
Info << "Start Time value is: " << runTime.startTime().value() << endl;

if (runTime.startTime().value() == 0)
{
    std::remove("Monitor_Iterations.txt");
    std::remove("Monitor_Volume_Fraction.txt");
    std::remove("Monitor_Time.txt");
    std::remove("Monitor_Compliance_Ratio.txt");
}

// Initialize variables
int opt(runTime.startTime().value() + 1);       // Optimization loop counter
int i;                                          // Initialize generic loop counter
int n_cells = gamma.size();                     // The number of cells in the domain
scalar len = std::pow(mesh.V()[0], 1.0 / 3);    // Characteristic length of domain
scalar area = gSum(mesh.V());                   // Domain area

scalar compliance_ratio = 1;                    // Compliance ratio of the solid
scalar compliance_prev;                         // Previous compliance ratio
scalar compliance_conv;                         // Percentage change in compliance ratio

scalar vol_frac;                                // Volume fraction of domain
scalar cost_vol_frac;                           // Cost variable for volume fraction
scalar vol_frac_conv;                           // Percentage difference between volume fraction and target
scalar eps(1e-8);                               // Small reference value for Young's Modulus

// Heaviside filter variables
double eta0, eta1, eta5 = 0, del = 0.0, y0, y5;
b.value()= 1.0 / (r_filter * len / 3.464) / (r_filter * len / 3.464);

// Initialize Method of Moving Asymptotes (MMA)
int m_constraints = 1;                                      // Number of constraints
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
for (i = 0; i < m_constraints; i++)                         // Resize constraint vector to size of domain
{
    dgdx[i].resize(n_cells);
}
