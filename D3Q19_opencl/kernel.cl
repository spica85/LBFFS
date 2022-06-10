// -- D3Q19 functions
inline int ic2i(int ic, int nx, int ny)
{
    return (ic%(nx*ny))%nx;
}

inline int ic2j(int ic, int nx, int ny)
{
    return ic%(nx*ny)/nx;
}

inline int ic2k(int ic, int nx, int ny)
{
    return ic/(nx*ny);
}

inline int index1d(int i, int j, int k, int nx, int ny)
{
    return nx*ny*k+nx*j+i;
}

inline int index1df(int q, int i, int j, int k, int nx, int ny, int nz)
{
    return q*nx*ny*nz+nx*ny*k+nx*j+i;
}

inline int idf(int q, int i, int nx, int ny, int nz)
{
    return q*nx*ny*nz+i;
}

inline int reflectQ(const int q)
{
    if(q == 0)
    {
        return 0;
    }
    else if(q == 1)
    {
        return 2;
    }
    else if(q == 2)
    {
        return 1;
    }
    else if(q == 3)
    {
        return 4;
    }
    else if(q == 4)
    {
        return 3;
    }
    else if(q == 5)
    {
        return 6;
    }
    else if(q == 6)
    {
        return 5;
    }
    else if(q == 7)
    {
        return 8;
    }
    else if(q == 8)
    {
        return 7;
    }
    else if(q == 9)
    {
        return 10;
    }
    else if(q == 10)
    {
        return 9;
    }
    else if(q == 11)
    {
        return 12;
    }
    else if(q == 12)
    {
        return 11;
    }
    else if(q == 13)
    {
        return 14;
    }
    else if(q == 14)
    {
        return 13;
    }
    else if(q == 15)
    {
        return 16;
    }
    else if(q == 16)
    {
        return 15;
    }
    else if(q == 17)
    {
        return 18;
    }
    else if(q == 18)
    {
        return 17;
    }
    // std::cerr << "Error q = [0:18], q = " << q << std::endl;
    // exit(EXIT_FAILURE);    
}

inline int upwindID_B(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz, const int boundary1, const int boundary2, const int boundary3)
{
    int ic = index1d(i, j, k, nx, ny);
    if(q == 0)
    {
        return index1d(i,j,k,nx,ny);
    }
    else if(q == 1)
    {
        return i != 0 ? index1d(i-1,j,k,nx,ny) : (boundary1 != 1 ? index1d(nx-1,j,k,nx,ny) : -1);
    }
    else if(q == 2)
    {
        return i != nx-1 ? index1d(i+1,j,k,nx,ny) : (boundary1 != 1 ? index1d(0,j,k,nx,ny) : -1);
    }
    else if(q == 3)
    {
        if(j != 0)
        {
            return index1d(i,j-1,k,nx,ny);
        }
        // else if(i != 0 && i != nx-1)
        // {
        //     return boundary1 != 1 ? index1d(i,ny-1,k,nx,ny) : -1;
        // }
        else
        {
            return boundary2 != 1 ? index1d(i,ny-1,k,nx,ny) : -1;
        }
    }
    else if(q == 4)
    {
        if(j != ny-1)
        {
            return index1d(i,j+1,k,nx,ny);
        }
        // else if(i != 0 && i != nx-1)
        // {
        //     return boundary1 != 1 ? index1d(i,0,k,nx,ny) : -1;
        // }
        else
        {
            return boundary2 != 1 ? index1d(i,0,k,nx,ny) : -1;
        }
    }
    else if(q == 5)
    {
        if(k != 0)
        {
            return index1d(i,j,k-1,nx,ny);
        }
        // else if(i != 0 && i != nx-1)
        // {
        //     if(j != 0 && j != ny-1)
        //     {
        //         return boundary1 != 1 ? index1d(i,j,nz-1,nx,ny) : -1;
        //     }
        //     else
        //     {
        //         return boundary2 != 1 ? index1d(i,j,nz-1,nx,ny) : -1;
        //     }
        // }
        else
        {
            // if(j != 0 && j != ny-1)
            // {
            //     return boundary2 != 1 ? index1d(i,j,nz-1,nx,ny) : -1;
            // }
            // else
            // {
                return boundary3 != 1 ? index1d(i,j,nz-1,nx,ny) : -1;
            // }
        }
    }
    else if(q == 6)
    {
        if(k != nz-1)
        {
            return index1d(i,j,k+1,nx,ny);
        }
        // else if(i != 0 && i != nx-1)
        // {
        //     if(j != 0 && j != ny-1)
        //     {
        //         return boundary1 != 1 ? index1d(i,j,0,nx,ny) : -1;
        //     }
        //     else
        //     {
        //         return boundary2 != 1 ? index1d(i,j,0,nx,ny) : -1;
        //     }
        // }
        else
        {
            // if(j != 0 && j != ny-1)
            // {
            //     return boundary2 != 1 ? index1d(i,j,0,nx,ny) : -1;
            // }
            // else
            // {
                return boundary3 != 1 ? index1d(i,j,0,nx,ny) : -1;
            // }
        }
    }
    else if(q == 7)
    {
        return (i != 0 && j != 0) ? index1d(i-1, j-1, k, nx, ny) :
               (i == 0 && j != 0) ? (boundary1 != 1 ? index1d(nx-1, j-1, k, nx, ny) : -1) :
               (i != 0 && j == 0) ? (boundary2 != 1 ? index1d(i-1, ny-1, k, nx, ny) : -1) :
               ((boundary1 != 1 && boundary2 != 1) ? index1d(nx-1, ny-1, k, nx, ny) : -1);
    }
    else if(q == 8)
    {
        return (i != nx-1 && j != ny-1) ? index1d(i+1, j+1, k, nx, ny) :
               (i == nx-1 && j != ny-1) ? (boundary1 != 1 ? index1d(0, j+1, k, nx, ny) : -1) :
               (i != nx-1 && j == ny-1) ? (boundary2 != 1 ? index1d(i+1, 0, k, nx, ny) : -1) :
               ((boundary1 != 1 && boundary2 != 1) ? index1d(0, 0, k, nx, ny) : -1);
    }
    else if(q == 9)
    {
        return (i != 0 && j != ny-1) ? index1d(i-1, j+1, k, nx, ny) :
               (i == 0 && j != ny-1) ? (boundary1 != 1 ? index1d(nx-1, j+1, k, nx, ny) : -1):
               (i != 0 && j == ny-1) ? (boundary2 != 1 ? index1d(i-1, 0, k, nx, ny) : -1) :
               ((boundary1 != 1 && boundary2 != 1) ? index1d(nx-1, 0, k, nx, ny) : -1);
    }
    else if(q == 10)
    {
        return (i != nx-1 && j != 0) ? index1d(i+1, j-1, k, nx, ny) :
               (i == nx-1 && j != 0) ? (boundary1 != 1 ? index1d(0, j-1, k, nx, ny) : -1) :
               (i != nx-1 && j == 0) ? (boundary2 != 1 ? index1d(i+1, ny-1, k, nx, ny) : -1) :
               ((boundary1 != 1 && boundary2 != 1) ? index1d(0, ny-1, k, nx, ny) : -1);
    }
    else if(q == 11)
    {
        if(i != 0 && k != 0)
        {
            return index1d(i-1, j, k-1, nx, ny);
        }
        else if(i == 0 && k != 0)
        {
            return boundary1 != 1 ? index1d(nx-1, j, k-1, nx, ny) : -1;
        }
        else if(i != 0 && k == 0)
        {
            // if(j != 0  && j != ny-1)
            // {
            //     return boundary1 != 1 ? index1d(i-1, j, nz-1, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(i-1, j, nz-1, nx, ny) : -1;
            // }
            return boundary3 != 1 ? index1d(i-1, j, nz-1, nx, ny) : -1;
        }
        else
        {
            // if(j != 0  && j != ny-1)
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(nx-1, j, nz-1, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1 && boundary3 != 1) ? index1d(nx-1, j, nz-1, nx, ny) : -1;
            // }
            return (boundary1 != 1 && boundary3 != 1) ? index1d(nx-1, j, nz-1, nx, ny) : -1;
        }
    }
    else if(q == 12)
    {
        if(i != nx-1 && k != nz-1)
        {
            return index1d(i+1, j, k+1, nx, ny);
        }
        else if(i == nx-1 && k != nz-1)
        {
            return boundary1 != 1 ? index1d(0, j, k+1, nx, ny) : -1;
        }
        else if(i != nx-1 && k == nz-1)
        {
            // if(j != 0  && j != ny-1)
            // {
            //     return boundary1 != 1 ? index1d(i+1, j, 0, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(i+1, j, 0, nx, ny) : -1;
            // }
            return boundary3 != 1 ? index1d(i+1, j, 0, nx, ny) : -1;
        }
        else
        {
            // if(j != 0  && j != ny-1)
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(0, j, 0, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1 && boundary3 != 1) ? index1d(0, j, 0, nx, ny) : -1;
            // }
            return (boundary1 != 1 && boundary3 != 1) ? index1d(0, j, 0, nx, ny) : -1;
        }
    }
    else if(q == 13)
    {
        if(i != 0 && k != nz-1)
        {
            return index1d(i-1, j, k+1, nx, ny);
        }
        else if(i == 0 && k != nz-1)
        {
            return boundary1 != 1 ? index1d(nx-1, j, k+1, nx, ny) : -1;
        }
        else if(i != 0 && k == nz-1)
        {
            // if(j != 0  && j != ny-1)
            // {
            //     return boundary1 != 1 ? index1d(i-1, j, 0, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(i-1, j, 0, nx, ny) : -1;
            // }
            return boundary3 != 1 ? index1d(i-1, j, 0, nx, ny) : -1;
        }
        else
        {
            // if(j != 0  && j != ny-1)
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(nx-1, j, 0, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1 && boundary3 != 1) ? index1d(nx-1, j, 0, nx, ny) : -1;
            // }
            return (boundary1 != 1 && boundary3 != 1) ? index1d(nx-1, j, 0, nx, ny) : -1;
        }
    }
    else if(q == 14)
    {
        if(i != nx-1 && k != 0)
        {
            return index1d(i+1, j, k-1, nx, ny);
        }
        else if(i == nx-1 && k != 0)
        {
            return boundary1 != 1 ? index1d(0, j, k-1, nx, ny) : -1;
        }
        else if(i != nx-1 && k == 0)
        {
            // if(j != 0  && j != ny-1)
            // {
            //     return boundary1 != 1 ? index1d(i+1, j, nz-1, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(i+1, j, nz-1, nx, ny) : -1;
            // }
            return boundary3 != 1 ? index1d(i+1, j, nz-1, nx, ny) : -1;
        }
        else
        {
            // if(j != 0  && j != ny-1)
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(0, j, nz-1, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1 && boundary3 != 1) ? index1d(0, j, nz-1, nx, ny) : -1;
            // }
            return (boundary1 != 1 && boundary3 != 1) ? index1d(0, j, nz-1, nx, ny) : -1;
        }
    }
    else if(q == 15)
    {
        if(j != 0 && k != 0)
        {
            return index1d(i, j-1, k-1, nx, ny);
        }
        else if(j == 0 && k != 0)
        {
            // return boundary1 != 1 ? index1d(i, ny-1, k-1, nx, ny) : -1;
            return boundary2 != 1 ? index1d(i, ny-1, k-1, nx, ny) : -1;
        }
        else if(j != 0 && k == 0)
        {
            // if(i != 0 && i != nx-1)
            // {
            //     return boundary1 != 1 ? index1d(i, j-1, nz-1, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(i, j-1, nz-1, nx, ny) : -1;
            // }
            return boundary3 != 1 ? index1d(i, j-1, nz-1, nx, ny) : -1;
        }
        else
        {
            // if(i != 0 && i != nx-1)
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(i, ny-1, nz-1, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1 && boundary3 != 1) ? index1d(i, ny-1, nz-1, nx, ny) : -1;
            // }
            return (boundary2 != 1 && boundary3 != 1) ? index1d(i, ny-1, nz-1, nx, ny) : -1;
        }
    }
    else if(q == 16)
    {
        if(j != ny-1 && k != nz-1)
        {
            return index1d(i, j+1, k+1, nx, ny);
        }
        else if(j == ny-1 && k != nz-1)
        {
            // return boundary1 != 1 ? index1d(i, 0, k+1, nx, ny) : -1;
            return boundary2 != 1 ? index1d(i, 0, k+1, nx, ny) : -1;
        }
        else if(j != ny-1 && k == nz-1)
        {
            // if(i != 0 && i != nx-1)
            // {
            //     return boundary1 != 1 ? index1d(i, j+1, 0, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(i, j+1, 0, nx, ny) : -1;
            // }
            return boundary3 != 1 ? index1d(i, j+1, 0, nx, ny) : -1;
        }
        else
        {
            // if(i != 0 && i != nx-1)
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(i, 0, 0, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1 && boundary3 != 1) ? index1d(i, 0, 0, nx, ny) : -1;
            // }
            return (boundary2 != 1 && boundary3 != 1) ? index1d(i, 0, 0, nx, ny) : -1;
        }
    }
    else if(q == 17)
    {
        if(j != 0 && k != nz-1)
        {
            return index1d(i, j-1, k+1, nx, ny);
        }
        else if(j == 0 && k != nz-1)
        {
            // return boundary1 != 1 ? index1d(i, ny-1, k+1, nx, ny) : -1;
            return boundary2 != 1 ? index1d(i, ny-1, k+1, nx, ny) : -1;
        }
        else if(j != 0 && k == nz-1)
        {
            // if(i != 0 && i != nx-1)
            // {
            //     return boundary1 != 1 ? index1d(i, j-1, 0, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(i, j-1, 0, nx, ny) : -1;
            // }
            return boundary3 != 1 ? index1d(i, j-1, 0, nx, ny) : -1;
        }
        else
        {
            // if(i != 0 && i != nx-1)
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(i, ny-1, 0, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1 && boundary3 != 1) ? index1d(i, ny-1, 0, nx, ny) : -1;
            // }
            return (boundary2 != 1 && boundary3 != 1) ? index1d(i, ny-1, 0, nx, ny) : -1;
        }
    }
    else if(q == 18)
    {
        if(j != ny-1 && k != 0)
        {
            return index1d(i, j+1, k-1, nx, ny);
        }
        else if(j == ny-1 && k != 0)
        {
            // return boundary1 != 1 ? index1d(i, 0, k-1, nx, ny) : -1;
            return boundary2 != 1 ? index1d(i, 0, k-1, nx, ny) : -1;
        }
        else if(j != ny-1 && k == 0)
        {
            // if(i != 0 && i != nx-1)
            // {
            //     return boundary1 != 1 ? index1d(i, j+1, nz-1, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(i, j+1, nz-1, nx, ny) : -1;
            // }
            return boundary3 != 1 ? index1d(i, j+1, nz-1, nx, ny) : -1;
        }
        else
        {
            // if(i != 0 && i != nx-1)
            // {
            //     return (boundary1 != 1 && boundary2 != 1) ? index1d(i, 0, nz-1, nx, ny) : -1;
            // }
            // else
            // {
            //     return (boundary1 != 1 && boundary2 != 1 && boundary3 != 1) ? index1d(i, 0, nz-1, nx, ny) : -1;
            // }
            return (boundary2 != 1 && boundary3 != 1) ? index1d(i, 0, nz-1, nx, ny) : -1;
        }
    }
    else
    {
        return 0;
    }
}

__kernel void k_streamingCollision // Pull
(
   __global float* f, __global float* fTmp,
   __global int* boundary1, __global int* boundary2, __global int* boundary3,
   __global float* sdf, __global unsigned char* solid, __global unsigned char* neiSolid, 
   __global float* u0, __global float* v0, __global float* w0,
   __global float* Fwx, __global float* Fwy, __global float* Fwz,
   const unsigned elements,
   const float omega,
   const float dpdx,
   const float rho_av,
   const int nx, const int ny, const int nz
)
{
    int ic = get_global_id(0);
    if(solid[ic] == 0)
    {
        int i = ic2i(ic,nx,ny);
        int j = ic2j(ic,nx,ny);
        int k = ic2k(ic,nx,ny);

        float wt[19] = {1.0f/3.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f};

        //                 0     1      2     3      4     5      6     7      8      9     10    11     12     13     14    15     16     17     18
        float cx[19] = {0.0f, 1.0f, -1.0f, 0.0f,  0.0f, 0.0f,  0.0f, 1.0f, -1.0f,  1.0f, -1.0f, 1.0f, -1.0f,  1.0f, -1.0f, 0.0f,  0.0f,  0.0f,  0.0f};
        float cy[19] = {0.0f, 0.0f,  0.0f, 1.0f, -1.0f, 0.0f,  0.0f, 1.0f, -1.0f, -1.0f,  1.0f, 0.0f,  0.0f,  0.0f,  0.0f, 1.0f, -1.0f,  1.0f, -1.0f};
        float cz[19] = {0.0f, 0.0f,  0.0f, 0.0f,  0.0f, 1.0f, -1.0f, 0.0f,  0.0f,  0.0f,  0.0f, 1.0f, -1.0f, -1.0f,  1.0f, 1.0f, -1.0f, -1.0f,  1.0f};

        float ft[19];
        int upID[19];
        ft[0] = f[ic];
        for(int q = 1; q < 19; q++)
        {
            int qic = q*elements +ic;

            upID[q] = upwindID_B(q,i,j,k,nx,ny,nz,boundary1[ic],boundary2[ic],boundary3[ic]);
            if(upID[q] != -1) // Streaming
            {   
                int upQID = idf(q, upID[q], nx, ny, nz);
                ft[q] = f[upQID];
            }
            else // Bounce-Back for boundary wall
            {
                int qbb = reflectQ(q);
                int bbQID = idf(qbb, ic, nx, ny, nz);
                const float rhow = rho_av;
                if(q == 1)
                {
                    ft[q] = f[bbQID] +6.0f*rhow*u0[ic]/18.0f;
                }
                else if(q == 2)
                {
                    ft[q] = f[bbQID] -6.0f*rhow*u0[ic]/18.0f;
                }
                else if(q == 3)
                {
                    ft[q] = f[bbQID] +6.0f*rhow*v0[ic]/18.0f;
                }
                else if(q == 4)
                {
                    ft[q] = f[bbQID] -6.0f*rhow*v0[ic]/18.0f;
                }
                else if(q == 5)
                {
                    ft[q] = f[bbQID] +6.0f*rhow*w0[ic]/18.0f;
                }
                else if(q == 6)
                {
                    ft[q] = f[bbQID] -6.0f*rhow*w0[ic]/18.0f;
                }
                else if(q == 7)
                {
                    ft[q] = f[bbQID] +6.0f*rhow*(u0[ic]+v0[ic])/36.0f;
                }
                else if(q == 8)
                {
                    ft[q] = f[bbQID] -6.0f*rhow*(u0[ic]+v0[ic])/36.0f;
                }
                else if(q == 9)
                {
                    ft[q] = f[bbQID] +6.0f*rhow*(u0[ic]-v0[ic])/36.0f;
                }
                else if(q == 10)
                {
                    ft[q] = f[bbQID] -6.0f*rhow*(u0[ic]-v0[ic])/36.0f;
                }
                else if(q == 11)
                {
                    ft[q] = f[bbQID] +6.0f*rhow*(u0[ic]+w0[ic])/36.0f;
                }
                else if(q == 12)
                {
                    ft[q] = f[bbQID] -6.0f*rhow*(u0[ic]+w0[ic])/36.0f;
                }
                else if(q == 13)
                {
                    ft[q] = f[bbQID] +6.0f*rhow*(u0[ic]-w0[ic])/36.0f;
                }
                else if(q == 14)
                {
                    ft[q] = f[bbQID] -6.0f*rhow*(u0[ic]-w0[ic])/36.0f;
                }
                else if(q == 15)
                {
                    ft[q] = f[bbQID] +6.0f*rhow*(v0[ic]+w0[ic])/36.0f;
                }
                else if(q == 16)
                {
                    ft[q] = f[bbQID] -6.0f*rhow*(v0[ic]+w0[ic])/36.0f;
                }
                else if(q == 17)
                {
                    ft[q] = f[bbQID] +6.0f*rhow*(v0[ic]-w0[ic])/36.0f;
                }
                else if(q == 18)
                {
                    ft[q] = f[bbQID] -6.0f*rhow*(v0[ic]-w0[ic])/36.0f;
                }
                // ft[q] = f[bbQID] -6.0f*rhow*wt[qbb]*(cx[qbb]*u0[ic]+cy[qbb]*v0[ic]+cz[qbb]*w0[ic]);
            }
        }

        // Bounce-Back for internal walls
        {
            float rho = 0.0f;
            float u = 0.0f;
            float v = 0.0f;
            float w = 0.0f;
            for(int q = 0; q < 19; q++)
            {
                int qic = q*elements +ic;
                rho += f[qic];

                u += f[qic]*cx[q];
                v += f[qic]*cy[q];
                w += f[qic]*cz[q];
            }
            u /= rho;
            v /= rho;
            w /= rho;
            float p = rho/3.f;

            Fwx[ic] = 0.f;
            Fwy[ic] = 0.f;
            Fwz[ic] = 0.f;

            for(int q = 1; q < 19; q++)
            {
                if(neiSolid[ic] == 1)
                {
                    if(solid[upID[q]] == 1)
                    {
                        const float sdf0 = sdf[ic];
                        const float sdf1 = sdf[upID[q]];
                        const float qf = fabs(sdf0)/(fabs(sdf0)+fabs(sdf1));

                        int qbb = reflectQ(q);
                        int bbQID = idf(qbb, ic, nx, ny, nz);
                        int upQID = idf(q, upID[qbb], nx, ny, nz);
                        int upQBBID = idf(qbb, upID[qbb], nx, ny, nz);

                        float uSqr =u*u+v*v+w*w;
                        float uDotC = -u*cx[q]-v*cy[q]-w*cz[q];
                        float feq = (rho+3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q];
                        
                        if(qf <= 0.5f)
                        {
                            // ft[q] = (1.f -2.f*qf)*ft[qbb] +(qf*f[bbQID])*2.f; // Bouzidi et al.'s Interpolated Bounce-Back
                            // ft[q] = (1.f -2.f*qf)*f[upQBBID] +(qf*f[bbQID])*2.f; // Bouzidi et al.'s Interpolated Bounce-Back (local)
                            // ft[q] = f[bbQID]; // Simple Bounce-Back

                            float chi = omega*(2.f*qf -1.f)/(1.f-omega);
                            ft[q] = (1.f -chi)*f[bbQID] +chi*feq; // Filippova & Hanel's Interpolated Bounce-Back (physically local)
                        }
                        else
                        {
                            // ft[q] = (1.f -0.5f/qf)*ft[upQID] +(0.5f/qf)*f[bbQID]; // Bouzidi et al.'s Interpolated Bounce-Back
                            // ft[q] = (1.f -0.5f/qf)*f[q] +(0.5f/qf)*f[bbQID]; // Bouzidi et al.'s Interpolated Bounce-Back (local)
                            // ft[q] = f[bbQID]; // Simple Bounce-Back

                            uSqr *= (1.f -1.f/qf)*(1.f -1.f/qf);
                            uDotC *= (1.f -1.f/qf);
                            float chi = omega*(2.f*qf -1.f);
                            ft[q] = (1.f -chi)*f[bbQID] +chi*feq; // Filippova & Hanel's Interpolated Bounce-Back (physically local)
                        }
                        Fwx[ic] += -(f[bbQID] + ft[q])*cx[q];
                        Fwy[ic] += -(f[bbQID] + ft[q])*cy[q];
                        Fwz[ic] += -(f[bbQID] + ft[q])*cz[q];
                    }
                }
            }
        }

        // Outflow Boundary (Geier et al., Comput. Math. Appl. (2015), Appendix F)
        for(int q = 0; q < 19; q++)
        {
            if(boundary1[ic] == 3 || boundary2[ic] == 3 || boundary3[ic] == 3)
            {
                int i = ic2i(ic,nx,ny);
                int j = ic2j(ic,nx,ny);
                int k = ic2k(ic,nx,ny);

                if(i == 0)
                {
                    if(boundary1[ic] == 3)
                    {
                        if(q == 1 || q == 7 || q == 9 || q == 11 || q == 13)
                        {                            
                            int innerID = index1d(1,j,k,nx,ny);
                            int qinic = idf(q, innerID, nx, ny, nz);
                            int qic = idf(q, ic, nx, ny, nz);
                            ft[q] = sqrt(3.f)*f[qinic] +(1.f -sqrt(3.f))*f[qic];
                        }
                    }
                }
                if(i == nx-1)
                {
                    if(boundary1[ic] == 3)
                    {
                        if(q == 2 || q == 8 || q == 10 || q == 12 || q == 14)
                        {
                            int innerID = index1d(nx-2,j,k,nx,ny);
                            int qinic = idf(q, innerID, nx, ny, nz);
                            int qic = idf(q, ic, nx, ny, nz);
                            ft[q] = f[qinic]/sqrt(3.f) +(1.f -1.f/sqrt(3.f))*f[qic];
                        }
                    }
                }
                if(j == 0)
                {
                    if(boundary2[ic] == 3)
                    {
                        if(q == 3 || q == 7 || q == 10 || q == 15 || q == 17)
                        {
                            int innerID = index1d(i,1,k,nx,ny);
                            int qinic = idf(q, innerID, nx, ny, nz);
                            int qic = idf(q, ic, nx, ny, nz);
                            ft[q] = f[qinic]/sqrt(3.f) +(1.f -1.f/sqrt(3.f))*f[qic];
                        }
                    }
                }
                if(j == ny-1)
                {
                    if(boundary2[ic] == 3)
                    {
                        if(q == 4 || q == 8 || q == 9 || q == 16 || q == 18)
                        {
                            int innerID = index1d(i,ny-2,k,nx,ny);
                            int qinic = idf(q, innerID, nx, ny, nz);
                            int qic = idf(q, ic, nx, ny, nz);
                            ft[q] = f[qinic]/sqrt(3.f) +(1.f -1.f/sqrt(3.f))*f[qic];
                        }
                    }
                }
                if(k == 0)
                {
                    if(boundary3[ic] == 3)
                    {
                        if(q == 5 || q == 11 || q == 14 || q == 15 || q == 18)
                        {
                            int innerID = index1d(i,j,1,nx,ny);
                            int qinic = idf(q, innerID, nx, ny, nz);
                            int qic = idf(q, ic, nx, ny, nz);
                            ft[q] = f[qinic]/sqrt(3.f) +(1.f -1.f/sqrt(3.f))*f[qic];
                        }
                    }
                }
                if(k == nz-1)
                {
                    if(boundary3[ic] == 3)
                    {
                        if(q == 6 || q == 12 || q == 13 || q == 16 || q == 17)
                        {
                            int innerID = index1d(i,j,nz-2,nx,ny);
                            int qinic = idf(q, innerID, nx, ny, nz);
                            int qic = idf(q, ic, nx, ny, nz);
                            ft[q] = f[qinic]/sqrt(3.f) +(1.f -1.f/sqrt(3.f))*f[qic];
                        }
                    }
                }
            }
        }
        
        // Collision
        {
            float rho = 0.0f;
            float u = 0.0f;
            float v = 0.0f;
            float w = 0.0f;
            for(int q = 0; q < 19; q++)
            {
                rho += ft[q];

                u += ft[q]*cx[q];
                v += ft[q]*cy[q];
                w += ft[q]*cz[q];
            }
            u /= rho;
            v /= rho;
            w /= rho;

            for(int q = 0; q < 19; q++)
            {
                float uSqr =u*u+v*v+w*w;
                float uDotC = u*cx[q]+v*cy[q]+w*cz[q];
                float feq = (1.0f+3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rho;

                int qic = q*elements +ic;

                fTmp[qic] = (1.0f -omega)*ft[q] + omega *feq +rho*wt[q]*3.0f*dpdx*cx[q]; // Pull

                // Equilibrium Boundary
                if(boundary1[ic] == 2 || boundary2[ic] == 2 || boundary3[ic] == 2)
                {
                    int i = ic2i(ic,nx,ny);
                    int j = ic2j(ic,nx,ny);
                    int k = ic2k(ic,nx,ny);

                    const float rhow = rho_av;
                    const float u = u0[ic];
                    const float v = v0[ic];
                    const float w = w0[ic];

                    if(i == 0)
                    {
                        if(boundary1[ic] == 2)
                        {
                            if(q == 1 || q == 7 || q == 9 || q == 11 || q == 13)
                            {                            
                                float uSqr = u*u +v*v +w*w;
                                float uDotC = u*cx[q]+v*cy[q]+w*cz[q];

                                ft[q] = (1.f +3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rhow;
                            }
                        }
                    }
                    if(i == nx-1)
                    {
                        if(boundary1[ic] == 2)
                        {
                            if(q == 2 || q == 8 || q == 10 || q == 12 || q == 14)
                            {
                                float uSqr = u*u +v*v +w*w;
                                float uDotC = u*cx[q]+v*cy[q]+w*cz[q];

                                ft[q] = (1.f +3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rhow;
                            }
                        }
                    }
                    if(j == 0)
                    {
                        if(boundary2[ic] == 2)
                        {
                            if(q == 3 || q == 7 || q == 10 || q == 15 || q == 17)
                            {
                                float uSqr = u*u +v*v +w*w;
                                float uDotC = u*cx[q]+v*cy[q]+w*cz[q];

                                ft[q] = (1.f +3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rhow;
                            }
                        }
                    }
                    if(j == ny-1)
                    {
                        if(boundary2[ic] == 2)
                        {
                            if(q == 4 || q == 8 || q == 9 || q == 16 || q == 18)
                            {
                                float uSqr = u*u +v*v +w*w;
                                float uDotC = u*cx[q]+v*cy[q]+w*cz[q];

                                ft[q] = (1.f +3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rhow;
                            }
                        }
                    }
                    if(k == 0)
                    {
                        if(boundary3[ic] == 2)
                        {
                            if(q == 5 || q == 11 || q == 14 || q == 15 || q == 18)
                            {
                                float uSqr = u*u +v*v +w*w;
                                float uDotC = u*cx[q]+v*cy[q]+w*cz[q];

                                ft[q] = (1.f +3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rhow;
                            }
                        }
                    }
                    if(k == nz-1)
                    {
                        if(boundary3[ic] == 2)
                        {
                            if(q == 6 || q == 12 || q == 13 || q == 16 || q == 17)
                            {
                                float uSqr = u*u +v*v +w*w;
                                float uDotC = u*cx[q]+v*cy[q]+w*cz[q];

                                ft[q] = (1.f +3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rhow;
                            }
                        }
                    }
                }
            }
        }
    }
}
