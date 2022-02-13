// rho, u, v
for(int k = 0; k < nz; k++)
{    
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            int ic = index1d(i,j,k,nx,ny);
                            
            // Update of rho
            rho[ic] = 0.0;
            for(int q = 0; q < 19; q++)
            {
                int qic = idf(q,ic,nx,ny,nz);
                rho[ic] += f[qic];
            }
                                    
            // Update of u, v
            u[ic] = 0.0;
            v[ic] = 0.0;
            w[ic] = 0.0;                
            for(int q = 0; q < 19; q++)
            {
                int qic = idf(q,ic,nx,ny,nz);
                u[ic] += f[qic]*cx[q];
                v[ic] += f[qic]*cy[q];
                w[ic] += f[qic]*cz[q];
            }
            u[ic] /= rho[ic];
            v[ic] /= rho[ic];
            w[ic] /= rho[ic];                
        }
    }
}

