// Initialization
for(int k = 0; k < nz; k++)
{        
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            int ic = index1d(i,j,k,nx,ny);
        
            p[ic] = rho0/3.0;
            u[ic] = 0.0;
            v[ic] = 0.0;
            w[ic] = 0.0;
        }
    }
}