// Setting for boundaries
for(int k = 0; k < nz; k++)
{                                    
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            int ic = index1d(i,j,k,nx,ny);

            if(i == 0)
            {
                boundary1[ic] = nxMinBoundary1;
                u0[ic] = nxMinU0;
                v0[ic] = nxMinV0;
                w0[ic] = nxMinW0;
            }
            else if(i == nx-1)
            {
                boundary1[ic] = nxMaxBoundary1;
                u0[ic] = nxMaxU0;
                v0[ic] = nxMaxV0;
                w0[ic] = nxMaxW0;
            }

            if(j == 0)
            {
                boundary2[ic] = nyMinBoundary2;
                u0[ic] = nyMinU0;
                v0[ic] = nyMinV0;
                w0[ic] = nyMinW0;
            }
            else if(j == ny-1)
            {
                boundary2[ic] = nyMaxBoundary2;
                u0[ic] = nyMaxU0;
                v0[ic] = nyMaxV0;
                w0[ic] = nyMaxW0;
            }

            if(k == 0)
            {
                boundary3[ic] = nzMinBoundary3;
                u0[ic] = nzMinU0;
                v0[ic] = nzMinV0;
                w0[ic] = nzMinW0;
            }
            else if(k == nz-1)
            {
                boundary3[ic] = nzMaxBoundary3;
                u0[ic] = nzMaxU0;
                v0[ic] = nzMaxV0;
                w0[ic] = nzMaxW0;
            }


            if(i == 0 && j == 0)
            {
                u0[ic] = (nxMinU0*nyMinU0 == 0.f) ? 0.f : nxMinU0;
                v0[ic] = (nxMinV0*nyMinV0 == 0.f) ? 0.f : nxMinV0;
                w0[ic] = (nxMinW0*nyMinW0 == 0.f) ? 0.f : nxMinW0;
            }
            if(i == 0 && j == ny-1)
            {
                u0[ic] = (nxMinU0*nyMaxU0 == 0.f) ? 0.f : nxMinU0;
                v0[ic] = (nxMinV0*nyMaxV0 == 0.f) ? 0.f : nxMinV0;
                w0[ic] = (nxMinW0*nyMaxW0 == 0.f) ? 0.f : nxMinW0;
            }
            if(i == nx-1 && j == 0)
            {
                u0[ic] = (nxMaxU0*nyMinU0 == 0.f) ? 0.f : nxMaxU0;
                v0[ic] = (nxMaxV0*nyMinV0 == 0.f) ? 0.f : nxMaxV0;
                w0[ic] = (nxMaxW0*nyMinW0 == 0.f) ? 0.f : nxMaxW0;
            }
            if(i == nx-1 && j == ny-1)
            {
                u0[ic] = (nxMaxU0*nyMaxU0 == 0.f) ? 0.f : nxMaxU0;
                v0[ic] = (nxMaxV0*nyMaxV0 == 0.f) ? 0.f : nxMaxV0;
                w0[ic] = (nxMaxW0*nyMaxW0 == 0.f) ? 0.f : nxMaxW0;
            }


            if(i == 0 && k == 0)
            {
                u0[ic] = (nz-1 != 0 && nxMinU0*nzMinU0 == 0.f) ? 0.f : nxMinU0;
                v0[ic] = (nz-1 != 0 && nxMinV0*nzMinV0 == 0.f) ? 0.f : nxMinV0;
                w0[ic] = (nz-1 != 0 && nxMinW0*nzMinW0 == 0.f) ? 0.f : nxMinW0;
            }
            if(i == 0 && k == nz-1)
            {
                u0[ic] = (nz-1 != 0 && nxMinU0*nzMaxU0 == 0.f) ? 0.f : nxMinU0;
                v0[ic] = (nz-1 != 0 && nxMinV0*nzMaxV0 == 0.f) ? 0.f : nxMinV0;
                w0[ic] = (nz-1 != 0 && nxMinW0*nzMaxW0 == 0.f) ? 0.f : nxMinW0;
            }
            if(i == nx-1 && k == 0)
            {
                u0[ic] = (nz-1 != 0 && nxMaxU0*nzMinU0 == 0.f) ? 0.f : nxMaxU0;
                v0[ic] = (nz-1 != 0 && nxMaxV0*nzMinV0 == 0.f) ? 0.f : nxMaxV0;
                w0[ic] = (nz-1 != 0 && nxMaxW0*nzMinW0 == 0.f) ? 0.f : nxMaxW0;
            }
            if(i == nx-1 && k == nz-1)
            {
                u0[ic] = (nz-1 != 0 && nxMaxU0*nzMaxU0 == 0.f) ? 0.f : nxMaxU0;
                v0[ic] = (nz-1 != 0 && nxMaxV0*nzMaxV0 == 0.f) ? 0.f : nxMaxV0;
                w0[ic] = (nz-1 != 0 && nxMaxW0*nzMaxW0 == 0.f) ? 0.f : nxMaxW0;
            }
            

            if(j == 0 && k == 0)
            {
                u0[ic] = (nz-1 != 0 && nyMinU0*nzMinU0 == 0.f) ? 0.f : nyMinU0;
                v0[ic] = (nz-1 != 0 && nyMinV0*nzMinV0 == 0.f) ? 0.f : nyMinV0;
                w0[ic] = (nz-1 != 0 && nyMinW0*nzMinW0 == 0.f) ? 0.f : nyMinW0;
            }
            if(j == 0 && k == nz-1)
            {
                u0[ic] = (nz-1 != 0 && nyMinU0*nzMaxU0 == 0.f) ? 0.f : nyMinU0;
                v0[ic] = (nz-1 != 0 && nyMinV0*nzMaxV0 == 0.f) ? 0.f : nyMinV0;
                w0[ic] = (nz-1 != 0 && nyMinW0*nzMaxW0 == 0.f) ? 0.f : nyMinW0;
            }
            if(j == ny-1 && k == 0)
            {
                u0[ic] = (nz-1 != 0 && nyMaxU0*nzMinU0 == 0.f) ? 0.f : nyMaxU0;
                v0[ic] = (nz-1 != 0 && nyMaxV0*nzMinV0 == 0.f) ? 0.f : nyMaxV0;
                w0[ic] = (nz-1 != 0 && nyMaxW0*nzMinW0 == 0.f) ? 0.f : nyMaxW0;
            }
            if(j == ny-1 && k == nz-1)
            {
                u0[ic] = (nz-1 != 0 && nyMaxU0*nzMaxU0 == 0.f) ? 0.f : nyMaxU0;
                v0[ic] = (nz-1 != 0 && nyMaxV0*nzMaxV0 == 0.f) ? 0.f : nyMaxV0;
                w0[ic] = (nz-1 != 0 && nyMaxW0*nzMaxW0 == 0.f) ? 0.f : nyMaxW0;
            }            


            

            if(i == 0 && j == 0 && k == 0)
            {
                u0[ic] = (nz-1 != 0 && nxMinU0*nyMinU0*nzMinU0 == 0.f) ? 0.f : nxMinU0;
                v0[ic] = (nz-1 != 0 && nxMinV0*nyMinV0*nzMinV0 == 0.f) ? 0.f : nxMinV0;
                w0[ic] = (nz-1 != 0 && nxMinW0*nyMinW0*nzMinW0 == 0.f) ? 0.f : nxMinW0;
            }
            if(i == nx-1 && j == 0 && k == 0)
            {
                u0[ic] = (nz-1 != 0 && nxMaxU0*nyMinU0*nzMinU0 == 0.f) ? 0.f : nxMaxU0;
                v0[ic] = (nz-1 != 0 && nxMaxV0*nyMinV0*nzMinV0 == 0.f) ? 0.f : nxMaxV0;
                w0[ic] = (nz-1 != 0 && nxMaxW0*nyMinW0*nzMinW0 == 0.f) ? 0.f : nxMaxW0;
            }
            if(i == 0 && j == ny-1 && k == 0)
            {
                u0[ic] = (nz-1 != 0 && nxMinU0*nyMaxU0*nzMinU0 == 0.f) ? 0.f : nxMinU0;
                v0[ic] = (nz-1 != 0 && nxMinV0*nyMaxV0*nzMinV0 == 0.f) ? 0.f : nxMinV0;
                w0[ic] = (nz-1 != 0 && nxMinW0*nyMaxW0*nzMinW0 == 0.f) ? 0.f : nxMinW0;
            }
            if(i == 0 && j == 0 && k == nz-1)
            {
                u0[ic] = (nz-1 != 0 && nxMinU0*nyMinU0*nzMaxU0 == 0.f) ? 0.f : nxMinU0;
                v0[ic] = (nz-1 != 0 && nxMinV0*nyMinV0*nzMaxV0 == 0.f) ? 0.f : nxMinV0;
                w0[ic] = (nz-1 != 0 && nxMinW0*nyMinW0*nzMaxW0 == 0.f) ? 0.f : nxMinW0;
            }
            if(i == nx-1 && j == ny-1 && k == 0)
            {
                u0[ic] = (nz-1 != 0 && nxMaxU0*nyMaxU0*nzMinU0 == 0.f) ? 0.f : nxMaxU0;
                v0[ic] = (nz-1 != 0 && nxMaxV0*nyMaxV0*nzMinV0 == 0.f) ? 0.f : nxMaxV0;
                w0[ic] = (nz-1 != 0 && nxMaxW0*nyMaxW0*nzMinW0 == 0.f) ? 0.f : nxMaxW0;
            }
            if(i == 0 && j == ny-1 && k == nz-1)
            {
                u0[ic] = (nz-1 != 0 && nxMinU0*nyMaxU0*nzMaxU0 == 0.f) ? 0.f : nxMinU0;
                v0[ic] = (nz-1 != 0 && nxMinV0*nyMaxV0*nzMaxV0 == 0.f) ? 0.f : nxMinV0;
                w0[ic] = (nz-1 != 0 && nxMinW0*nyMaxW0*nzMaxW0 == 0.f) ? 0.f : nxMinW0;
            }
            if(i == nx-1 && j == 0 && k == nz-1)
            {
                u0[ic] = (nz-1 != 0 && nxMaxU0*nyMinU0*nzMaxU0 == 0.f) ? 0.f : nxMaxU0;
                v0[ic] = (nz-1 != 0 && nxMaxV0*nyMinV0*nzMaxV0 == 0.f) ? 0.f : nxMaxV0;
                w0[ic] = (nz-1 != 0 && nxMaxW0*nyMinW0*nzMaxW0 == 0.f) ? 0.f : nxMaxW0;
            }
            if(i == nx-1 && j == ny-1 && k == nz-1)
            {
                u0[ic] = (nz-1 != 0 && nxMaxU0*nyMaxU0*nzMaxU0 == 0.f) ? 0.f : nxMaxU0;
                v0[ic] = (nz-1 != 0 && nxMaxV0*nyMaxV0*nzMaxV0 == 0.f) ? 0.f : nxMaxV0;
                w0[ic] = (nz-1 != 0 && nxMaxW0*nyMaxW0*nzMaxW0 == 0.f) ? 0.f : nxMaxW0;
            }
        }
    }
}