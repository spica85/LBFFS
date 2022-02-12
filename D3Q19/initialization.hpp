// Initialization
for(int k = 0; k < nz; k++)
{        
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            int ic = index1d(i,j,k,nx,ny);
        
            float rho = rho0;
            float u = 0.0;
            float v = 0.0;
            float w = 0.0;
            for(int q = 0; q < 19; q++)
            {
                float uSqr =u*u+v*v+w*w;
                float uDotC = u*cx[q]+v*cy[q]+w*cz[q];
                float feq = (1.0+3.0*uDotC +4.5*uDotC*uDotC -1.5*uSqr)*wt[q]*rho;

                int icf = idf(q,ic,nx,ny,nz);
                f[icf] = feq;
                ftmp[icf] = f[icf];
            }
        }
    }
}