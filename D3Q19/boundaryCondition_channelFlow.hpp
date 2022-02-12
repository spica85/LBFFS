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
                obst[ic].boundary = 0;
                obst[ic].normal = -1;
                obst[ic].inner = false;
            }
            else if(i == nx-1)
            {
                obst[ic].boundary = 0;
                obst[ic].normal = 1;
                obst[ic].inner = false;                                  
            }
            // else if(j == 0)
            // {
            //     obst[ic].boundary = 1;
            //     obst[ic].normal = -2;
            //     obst[ic].inner = false;
            // }
            // else if(j == ny-1)
            // {
            //     obst[ic].boundary = 1;
            //     obst[ic].normal = 2;
            //     obst[ic].inner = false;
            // }                    
            else if(k == 0)
            {
                obst[ic].boundary = 0;
                obst[ic].normal = -3;
                obst[ic].inner = false;
            }
            else if(k == nz-1)
            {
                obst[ic].boundary = 0;
                obst[ic].normal = 3;
                obst[ic].inner = false;
            }                
            else
            {
                obst[ic].boundary = 0;
                obst[ic].normal = 0;
                obst[ic].inner = false;                                                       
            }
            if(j == 0)
            {
                obst[ic].boundary = 1;
                obst[ic].normal = -2;
                obst[ic].inner = false;
            }
            else if(j == ny-1)
            {
                obst[ic].boundary = 1;
                obst[ic].normal = 2;
                obst[ic].inner = false;
            }                    
        }
    }
}