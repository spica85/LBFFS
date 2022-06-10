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
                obst[ic].boundary3 = 0;
                obst[ic].normal = {0,0,-1};
                obst[ic].inner = false;
                obst[ic].v0 = 0.0;
            }
            else if(k == nz-1)
            {
                obst[ic].boundary3 = 0;
                obst[ic].normal = {0,0,1};
                obst[ic].inner = false;
                obst[ic].v0 = 0.0;
            }

            if(j == 0)
            {
                obst[ic].boundary2 = 1;
                obst[ic].normal = {0,-1,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            else if(j == ny-1)
            {
                obst[ic].boundary2 = 1;
                obst[ic].normal = {0,1,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }

            if(i == 0)
            {
                obst[ic].boundary1 = 1;
                obst[ic].normal = {-1,0,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            else if(i == nx-1)
            {
                obst[ic].boundary1 = 3;
                obst[ic].normal = {1,0,0};
                obst[ic].inner = false;                                  
                obst[ic].u0 = uMax/c;
            }

            if(k == 0 && j == 0)
            {
                obst[ic].boundary2 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {0,-1,-1};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(k == 0 && j == ny-1)
            {
                obst[ic].boundary2 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {0,1,-1};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(k == nz-1 && j == 0)
            {
                obst[ic].boundary2 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {0,-1,1};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(k == nz-1 && j == ny-1)
            {
                obst[ic].boundary2 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {0,1,1};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }            


            if(i == 0 && j == 0)
            {
                obst[ic].boundary1 = 1;
                obst[ic].boundary2 = 1;
                obst[ic].normal = {-1,-1,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == 0 && j == ny-1)
            {
                obst[ic].boundary1 = 1;
                obst[ic].boundary2 = 1;
                obst[ic].normal = {-1,1,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == nx-1 && j == 0)
            {
                obst[ic].boundary1 = 3;
                obst[ic].boundary2 = 1;
                obst[ic].normal = {1,-1,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == nx-1 && j == ny-1)
            {
                obst[ic].boundary1 = 3;
                obst[ic].boundary2 = 1;
                obst[ic].normal = {1,1,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }


            if(i == 0 && k == 0)
            {
                obst[ic].boundary1 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {-1,-1,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == 0 && k == nz-1)
            {
                obst[ic].boundary1 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {-1,1,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == nx-1 && k == 0)
            {
                obst[ic].boundary1 = 3;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {1,-1,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == nx-1 && k == nz-1)
            {
                obst[ic].boundary1 = 3;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {1,1,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }

            if(i == 0 && j == 0 && k == 0)
            {
                obst[ic].boundary1 = 1;
                obst[ic].boundary2 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {-1,-1,-1};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == nx-1 && j == 0 && k == 0)
            {
                obst[ic].boundary1 = 3;
                obst[ic].boundary2 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {1,-1,-1};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == 0 && j == ny-1 && k == 0)
            {
                obst[ic].boundary1 = 1;
                obst[ic].boundary2 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {-1,1,-1};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == 0 && j == 0 && k == nz-1)
            {
                obst[ic].boundary1 = 1;
                obst[ic].boundary2 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {-1,-1,1};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == nx-1 && j == ny-1 && k == 0)
            {
                obst[ic].boundary1 = 3;
                obst[ic].boundary2 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {1,1,-1};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == 0 && j == ny-1 && k == nz-1)
            {
                obst[ic].boundary1 = 1;
                obst[ic].boundary2 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {-1,1,1};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == nx-1 && j == 0 && k == nz-1)
            {
                obst[ic].boundary1 = 3;
                obst[ic].boundary2 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {1,-1,1};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            if(i == nx-1 && j == ny-1 && k == nz-1)
            {
                obst[ic].boundary1 = 3;
                obst[ic].boundary2 = 1;
                obst[ic].boundary3 = 0;
                obst[ic].normal = {1,1,1};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
        }
    }
}