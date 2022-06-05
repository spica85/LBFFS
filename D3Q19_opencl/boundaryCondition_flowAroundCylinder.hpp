// // Setting for boundaries
// for(int k = 0; k < nz; k++)
// {                                    
//     for(int j = 0; j < ny; j++)
//     {
//         for(int i = 0; i < nx; i++)
//         {
//             int ic = index1d(i,j,k,nx,ny);

//             obst[ic].boundary = 0;
//             obst[ic].normal = {0,0,0};
//             obst[ic].inner = false;
//             obst[ic].u0 = 0.0f;
//             obst[ic].v0 = 0.0f;
//             obst[ic].w0 = 0.0f;

//             if(k == 0)
//             {
//                 obst[ic].boundary = 0;
//                 obst[ic].normal = {0,0,-1};
//                 obst[ic].inner = false;
//                 obst[ic].v0 = 0.0;
//             }
//             else if(k == nz-1)
//             {
//                 obst[ic].boundary = 0;
//                 obst[ic].normal = {0,0,1};
//                 obst[ic].inner = false;
//                 obst[ic].v0 = 0.0;
//             }

//             if(i == 0)
//             {
//                 obst[ic].boundary = 0;
//                 obst[ic].normal = {-1,0,0};
//                 obst[ic].inner = false;
//             }
//             else if(i == nx-1)
//             {
//                 obst[ic].boundary = 0;
//                 obst[ic].normal = {1,0,0};
//                 obst[ic].inner = false;                                  
//             }
//             // else if(j == 0)
//             // {
//             //     obst[ic].boundary = 1;
//             //     obst[ic].normal = -2;
//             //     obst[ic].inner = false;
//             // }
//             // else if(j == ny-1)
//             // {
//             //     obst[ic].boundary = 1;
//             //     obst[ic].normal = 2;
//             //     obst[ic].inner = false;
//             // }                    
//             // if(k == 0)
//             // {
//             //     obst[ic].boundary = 0;
//             //     obst[ic].normal = {0,0,-1};
//             //     obst[ic].inner = false;
//             // }
//             // else if(k == nz-1)
//             // {
//             //     obst[ic].boundary = 0;
//             //     obst[ic].normal = {0,0,1};
//             //     obst[ic].inner = false;
//             // }                
            
//             if(j == 0)
//             {
//                 obst[ic].boundary = 1;
//                 obst[ic].normal = {0,-1,0};
//                 obst[ic].inner = false;
//             }
//             else if(j == ny-1)
//             {
//                 obst[ic].boundary = 1;
//                 obst[ic].normal = {0,1,0};
//                 obst[ic].inner = false;
//             }                    
//         }
//     }
// }

// // Setting for boundaries
for(int k = 0; k < nz; k++)
{                                    
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            int ic = index1d(i,j,k,nx,ny);

            obst[ic].boundary = 0;
            obst[ic].normal = {0,0,0};
            obst[ic].inner = false;

            if(k == 0)
            {
                obst[ic].boundary = 0;
                obst[ic].normal = {0,0,-1};
                obst[ic].inner = false;
                obst[ic].v0 = 0.0;
            }
            else if(k == nz-1)
            {
                obst[ic].boundary = 0;
                obst[ic].normal = {0,0,1};
                obst[ic].inner = false;
                obst[ic].v0 = 0.0;
            }

            if(j == 0)
            {
                obst[ic].boundary = 1;
                obst[ic].normal = {0,-1,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            else if(j == ny-1)
            {
                obst[ic].boundary = 1;
                obst[ic].normal = {0,1,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }

            if(i == 0)
            {
                obst[ic].boundary = 1;
                obst[ic].normal = {-1,0,0};
                obst[ic].inner = false;
                obst[ic].u0 = uMax/c;
            }
            else if(i == nx-1)
            {
                obst[ic].boundary = 0;
                obst[ic].normal = {1,0,0};
                obst[ic].inner = false;                                  
                obst[ic].u0 = uMax/c;
            }

            // if(k == 0 && j == 0)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {0,-1,-1};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(k == 0 && j == ny-1)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {0,1,-1};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(k == nz-1 && j == 0)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {0,-1,1};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(k == nz-1 && j == ny-1)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {0,1,1};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }            


            // if(i == 0 && j == 0)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {-1,-1,0};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(i == 0 && j == ny-1)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {-1,1,0};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(i == nx-1 && j == 0)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {1,-1,0};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(i == nx-1 && j == ny-1)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {1,1,0};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }

            // if(i == 0 && j == 0 && k == 0)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {-1,-1,-1};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(i == nx-1 && j == 0 && k == 0)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {1,-1,-1};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(i == 0 && j == ny-1 && k == 0)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {-1,1,-1};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(i == 0 && j == 0 && k == nz-1)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {-1,-1,1};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(i == nx-1 && j == ny-1 && k == 0)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {1,1,-1};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(i == 0 && j == ny-1 && k == nz-1)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {-1,1,1};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(i == nx-1 && j == 0 && k == nz-1)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {1,-1,1};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }
            // if(i == nx-1 && j == ny-1 && k == nz-1)
            // {
            //     obst[ic].boundary = 0;
            //     obst[ic].normal = {1,1,1};
            //     obst[ic].inner = false;
            //     obst[ic].u0 = 0.0;
            // }

        }
    }
}