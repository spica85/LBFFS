// Initialization
for(int k = 0; k < nz; k++)
{        
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            int ic = index1d(i,j,k,nx,ny);
        
            rho[ic] = rho0;
            u[ic] = 0.0;
            v[ic] = 0.0;
            w[ic] = 0.0;                 

            for(int q = 0; q < 19; q++)
            {
                float uSqr =u[ic]*u[ic]+v[ic]*v[ic]+w[ic]*w[ic];
                float uDotC = u[ic]*cx[q]+v[ic]*cy[q]+w[ic]*cz[q];
                float feq = (1.0+3.0*uDotC +4.5*uDotC*uDotC -1.5*uSqr)*wt[q]*rho[ic];
            
                int icf = idf(q, ic, nx, ny, nz);                       
                f[icf] = feq;
                ftmp[icf] = f[icf];                        
            }
        }
    }
}