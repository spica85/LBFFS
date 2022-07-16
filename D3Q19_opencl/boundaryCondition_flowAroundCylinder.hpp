// Setting for boundaries
for(int k = 0; k < nz; k++)
{                                    
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            int ic = index1d(i,j,k,nx,ny);

            if(k == 0)
            {
                boundary3[ic] = 0;
                v0[ic] = 0.0;
            }
            else if(k == nz-1)
            {
                boundary3[ic] = 0;
                v0[ic] = 0.0;
            }

            if(j == 0)
            {
                boundary2[ic] = 1;
                u0[ic] = uMax/c;
            }
            else if(j == ny-1)
            {
                boundary2[ic] = 1;
                u0[ic] = uMax/c;
            }

            if(i == 0)
            {
                boundary1[ic] = 1;
                u0[ic] = uMax/c;
            }
            else if(i == nx-1)
            {
                boundary1[ic] = 3;
                u0[ic] = uMax/c;
            }

            if(k == 0 && j == 0)
            {
                boundary2[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(k == 0 && j == ny-1)
            {
                boundary2[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(k == nz-1 && j == 0)
            {
                boundary2[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(k == nz-1 && j == ny-1)
            {
                boundary2[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }            


            if(i == 0 && j == 0)
            {
                boundary1[ic] = 1;
                boundary2[ic] = 1;
                u0[ic] = uMax/c;
            }
            if(i == 0 && j == ny-1)
            {
                boundary1[ic] = 1;
                boundary2[ic] = 1;
                u0[ic] = uMax/c;
            }
            if(i == nx-1 && j == 0)
            {
                boundary1[ic] = 3;
                boundary2[ic] = 1;
                u0[ic] = uMax/c;
            }
            if(i == nx-1 && j == ny-1)
            {
                boundary1[ic] = 3;
                boundary2[ic] = 1;
                u0[ic] = uMax/c;
            }


            if(i == 0 && k == 0)
            {
                boundary1[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(i == 0 && k == nz-1)
            {
                boundary1[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(i == nx-1 && k == 0)
            {
                boundary1[ic] = 3;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(i == nx-1 && k == nz-1)
            {
                boundary1[ic] = 3;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }

            if(i == 0 && j == 0 && k == 0)
            {
                boundary1[ic] = 1;
                boundary2[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(i == nx-1 && j == 0 && k == 0)
            {
                boundary1[ic] = 3;
                boundary2[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(i == 0 && j == ny-1 && k == 0)
            {
                boundary1[ic] = 1;
                boundary2[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(i == 0 && j == 0 && k == nz-1)
            {
                boundary1[ic] = 1;
                boundary2[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(i == nx-1 && j == ny-1 && k == 0)
            {
                boundary1[ic] = 3;
                boundary2[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(i == 0 && j == ny-1 && k == nz-1)
            {
                boundary1[ic] = 1;
                boundary2[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(i == nx-1 && j == 0 && k == nz-1)
            {
                boundary1[ic] = 3;
                boundary2[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
            if(i == nx-1 && j == ny-1 && k == nz-1)
            {
                boundary1[ic] = 3;
                boundary2[ic] = 1;
                boundary3[ic] = 0;
                u0[ic] = uMax/c;
            }
        }
    }
}