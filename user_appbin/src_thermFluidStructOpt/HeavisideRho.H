
//***********************************//PDE filter
fvScalarMatrix Eqn4(fvm::laplacian(gamma_filter)-fvm::Sp(b,gamma_filter)+x*b);
Eqn4.solve();

//***********************************//Heaviside_filter
del=Foam::min(0.2*opt,100);
eta0=0.0001;
eta1=1.0;
y0=diff(gamma_filter, mesh.V(), del, eta0, n_cells);
reduce(y0, sumOp<scalar>());
do
{
  eta5=(eta0+eta1)/2.0;
  y5=diff(gamma_filter, mesh.V(), del, eta5, n_cells);
  reduce(y5, sumOp<scalar>());
  if(y0*y5<0)
  { 
     eta1=eta5;
  }
  else
  { 
     eta0=eta5;
     y0=y5;
  }
} while ((eta1-eta0)>0.0001);

for(i=0;i<n_cells;i++)
{
   if(gamma_filter[i]<=eta5)
   {
        gamma[i]=eta5*(Foam::exp(-del*(1-gamma_filter[i]/eta5))-(1-gamma_filter[i]/eta5)*Foam::exp(-del));
   }
   else
   {
        gamma[i]=eta5+(1-eta5)*(1-Foam::exp(-del*(gamma_filter[i]-eta5)/(1-eta5))+(gamma_filter[i]-eta5)*Foam::exp(-del)/(1-eta5));
   }   
}

gamma.correctBoundaryConditions();
