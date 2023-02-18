
gradD = fvc::grad(D);
volSymmTensorField sigma(rho * (mu * twoSymm(gradD) + (lambda * I) * tr(gradD)));  
volScalarField Energy((gradD + gradD.T())&&sigma); 

// Volume use
V = 0;     
for(i = 0; i < allcells; i++)
{
  V = gamma[i] + V;
}

reduce(V, sumOp<scalar>()); 
V = V / nallcells - voluse;

// Compliance
C_prev = C;
for(i = 0; i < allcells; i++)
{
  C = C + Energy[i];
}

reduce(C, sumOp<scalar>());
C = C / 2.0 * len * len * len;
dC = C - C_prev;
if(opt==1)
{
    dC = 0;
}
dC_sum = dC_sum * opt_converge_smooth + dC * dC;
dC_RMS = std::sqrt(dC_sum / opt) / std::abs(C + 0.0000000005);

if(Pstream::master())
{
    ofstream outfile1("Voluse.txt", std::ios::app);
    outfile1 << V << "\n";
    outfile1.close();

    ofstream outfile2("Compliance.txt", std::ios::app);
    outfile2 << C << "\n";
    outfile2.close();

    ofstream outfile3("Time.txt", std::ios::app);
    outfile3 << runTime.elapsedClockTime() << "\n";
    outfile3.close();

    ofstream outfile4("dC_RMS.txt", std::ios::app);
    outfile4 << dC_RMS << "\n";
    outfile4.close();
}
