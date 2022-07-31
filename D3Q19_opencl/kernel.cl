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
        else
        {
            return boundary3 != 1 ? index1d(i,j,nz-1,nx,ny) : -1;
        }
    }
    else if(q == 6)
    {
        if(k != nz-1)
        {
            return index1d(i,j,k+1,nx,ny);
        }
        else
        {
            return boundary3 != 1 ? index1d(i,j,0,nx,ny) : -1;
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
            return boundary3 != 1 ? index1d(i-1, j, nz-1, nx, ny) : -1;
        }
        else
        {
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
            return boundary3 != 1 ? index1d(i+1, j, 0, nx, ny) : -1;
        }
        else
        {
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
            return boundary3 != 1 ? index1d(i-1, j, 0, nx, ny) : -1;
        }
        else
        {
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
            return boundary3 != 1 ? index1d(i+1, j, nz-1, nx, ny) : -1;
        }
        else
        {
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
            return boundary2 != 1 ? index1d(i, ny-1, k-1, nx, ny) : -1;
        }
        else if(j != 0 && k == 0)
        {
            return boundary3 != 1 ? index1d(i, j-1, nz-1, nx, ny) : -1;
        }
        else
        {
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
            return boundary2 != 1 ? index1d(i, 0, k+1, nx, ny) : -1;
        }
        else if(j != ny-1 && k == nz-1)
        {
            return boundary3 != 1 ? index1d(i, j+1, 0, nx, ny) : -1;
        }
        else
        {
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
            return boundary2 != 1 ? index1d(i, ny-1, k+1, nx, ny) : -1;
        }
        else if(j != 0 && k == nz-1)
        {
            return boundary3 != 1 ? index1d(i, j-1, 0, nx, ny) : -1;
        }
        else
        {
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
            return boundary2 != 1 ? index1d(i, 0, k-1, nx, ny) : -1;
        }
        else if(j != ny-1 && k == 0)
        {
            return boundary3 != 1 ? index1d(i, j+1, nz-1, nx, ny) : -1;
        }
        else
        {
            return (boundary2 != 1 && boundary3 != 1) ? index1d(i, 0, nz-1, nx, ny) : -1;
        }
    }
    else
    {
        return 0;
    }
}

int qicNear(int q, int ic, int iNear, int nx, int ny, int nz)
{
    int i = ic2i(ic,nx,ny);
    int j = ic2j(ic,nx,ny);
    int k = ic2k(ic,nx,ny);

    if(i == 0 || i == nx-1 || j == 0 || j == ny-1 || k == 0 || k == nz-1)
    {
        printf("Error: Missing near points!\n");
    }

    if(iNear == 0)
    {
        return idf(q, ic, nx, ny, nz);
    }
    else if(iNear == 1)
    {
        return index1df(q, i+1, j, k, nx, ny, nz);
    }
    else if(iNear == 2)
    {
        return index1df(q, i-1, j, k, nx, ny, nz);
    }
    else if(iNear == 3)
    {
        return index1df(q, i, j+1, k, nx, ny, nz);
    }
    else if(iNear == 4)
    {
        return index1df(q, i, j-1, k, nx, ny, nz);
    }
    else if(iNear == 5)
    {
        return index1df(q, i, j, k+1, nx, ny, nz);
    }
    else if(iNear == 6)
    {
        return index1df(q, i+1, j, k-1, nx, ny, nz);
    }
    else if(iNear == 7)
    {
        return index1df(q, i+1, j+1, k, nx, ny, nz);
    }
    else if(iNear == 8)
    {
        return index1df(q, i-1, j-1, k, nx, ny, nz);
    }
    else if(iNear == 9)
    {
        return index1df(q, i+1, j-1, k, nx, ny, nz);
    }
    else if(iNear == 10)
    {
        return index1df(q, i-1, j+1, k, nx, ny, nz);
    }
    else if(iNear == 11)
    {
        return index1df(q, i+1, j, k+1, nx, ny, nz);
    }
    else if(iNear == 12)
    {
        return index1df(q, i-1, j, k-1, nx, ny, nz);
    }
    else if(iNear == 13)
    {
        return index1df(q, i+1, j, k-1, nx, ny, nz);
    }
    else if(iNear == 14)
    {
        return index1df(q, i-1, j, k+1, nx, ny, nz);
    }
    else if(iNear == 15)
    {
        return index1df(q, i, j+1, k+1, nx, ny, nz);
    }
    else if(iNear == 16)
    {
        return index1df(q, i, j-1, k-1, nx, ny, nz);
    }
    else if(iNear == 17)
    {
        return index1df(q, i, j+1, k-1, nx, ny, nz);
    }
    else if(iNear == 18)
    {
        return index1df(q, i, j-1, k+1, nx, ny, nz);
    }
    else if(iNear == 19)
    {
        return index1df(q, i+1, j+1, k+1, nx, ny, nz);
    }
    else if(iNear == 20)
    {
        return index1df(q, i-1, j-1, k-1, nx, ny, nz);
    }
    else if(iNear == 21)
    {
        return index1df(q, i+1, j-1, k+1, nx, ny, nz);
    }
    else if(iNear == 22)
    {
        return index1df(q, i-1, j+1, k-1, nx, ny, nz);
    }
    else if(iNear == 23)
    {
        return index1df(q, i+1, j+1, k-1, nx, ny, nz);
    }
    else if(iNear == 24)
    {
        return index1df(q, i-1, j-1, k+1, nx, ny, nz);
    }
    else if(iNear == 25)
    {
        return index1df(q, i+1, j-1, k-1, nx, ny, nz);
    }
    else if(iNear == 26)
    {
        return index1df(q, i-1, j+1, k+1, nx, ny, nz);
    }
}

void cal_rhoUVW(const float* f, float* rho, float* u, float* v, float* w)
{
    float cx[19] = {0.0f, 1.0f, -1.0f, 0.0f,  0.0f, 0.0f,  0.0f, 1.0f, -1.0f,  1.0f, -1.0f, 1.0f, -1.0f,  1.0f, -1.0f, 0.0f,  0.0f,  0.0f,  0.0f};
    float cy[19] = {0.0f, 0.0f,  0.0f, 1.0f, -1.0f, 0.0f,  0.0f, 1.0f, -1.0f, -1.0f,  1.0f, 0.0f,  0.0f,  0.0f,  0.0f, 1.0f, -1.0f,  1.0f, -1.0f};
    float cz[19] = {0.0f, 0.0f,  0.0f, 0.0f,  0.0f, 1.0f, -1.0f, 0.0f,  0.0f,  0.0f,  0.0f, 1.0f, -1.0f, -1.0f,  1.0f, 1.0f, -1.0f, -1.0f,  1.0f};

    *rho = 0.0f;
    *u = 0.0f;
    *v = 0.0f;
    *w = 0.0f;

    for(int q = 0; q < 19; q++)
    {
        *rho += f[q];
        *u += f[q]*cx[q];
        *v += f[q]*cy[q];
        *w += f[q]*cz[q];
    }
    *u /= *rho;
    *v /= *rho;
    *w /= *rho;
}

__kernel void k_streamingCollision // Pull
(
   __global float* f, __global float* fTmp,
   __global int* boundary1, __global int* boundary2, __global int* boundary3,
   __global float* sdf, __global unsigned char* solid, __global unsigned char* neiSolid, 
   __global float* u0, __global float* v0, __global float* w0,
   __global float* Fwx, __global float* Fwy, __global float* Fwz,
   __global float* tauSGS,
   const unsigned elements,
   const float omega,
   const float dpdx,
   const float rho_av,
   const int nx, const int ny, const int nz,
   const float LES
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

        //-- Bounce-Back for internal walls
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
                        // float feq = (1.0f+3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*rho*wt[q];

                        float tau = 1.f/omega;
                        float omegaEff = 1.f/(tau +tauSGS[ic]);
                        
                        if(qf <= 0.5f)
                        {
                            // ft[q] = (1.f -2.f*qf)*ft[qbb] +(qf*f[bbQID])*2.f; // Bouzidi et al.'s Interpolated Bounce-Back
                            // ft[q] = (1.f -2.f*qf)*f[upQBBID] +(qf*f[bbQID])*2.f; // Bouzidi et al.'s Interpolated Bounce-Back (local)
                            // ft[q] = f[bbQID]; // Simple Bounce-Back

                            float chi = omegaEff*(2.f*qf -1.f)/(1.f-omegaEff);
                            ft[q] = (1.f -chi)*f[bbQID] +chi*feq; // Filippova & Hanel's Interpolated Bounce-Back (physically local)
                        }
                        else
                        {
                            // ft[q] = (1.f -0.5f/qf)*ft[upQID] +(0.5f/qf)*f[bbQID]; // Bouzidi et al.'s Interpolated Bounce-Back
                            // ft[q] = (1.f -0.5f/qf)*f[q] +(0.5f/qf)*f[bbQID]; // Bouzidi et al.'s Interpolated Bounce-Back (local)
                            // ft[q] = f[bbQID]; // Simple Bounce-Back

                            uSqr *= (1.f -1.f/qf)*(1.f -1.f/qf);
                            uDotC *= (1.f -1.f/qf);
                            float chi = omegaEff*(2.f*qf -1.f);
                            ft[q] = (1.f -chi)*f[bbQID] +chi*feq; // Filippova & Hanel's Interpolated Bounce-Back (physically local)
                        }
                        Fwx[ic] += -(f[bbQID] + ft[q])*cx[q];
                        Fwy[ic] += -(f[bbQID] + ft[q])*cy[q];
                        Fwz[ic] += -(f[bbQID] + ft[q])*cz[q];
                    }
                }
            }
        }
        //--

        //-- Outflow Boundary (Geier et al., Comput. Math. Appl. (2015), Appendix F)
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
        //--

        //-- IBM
        {
            float rhoNear[27];
            float uNear[27];
            float vNear[27];
            float wNear[27];
            float fNear[27][19];

            int i = ic2i(ic,nx,ny);
            int j = ic2j(ic,nx,ny);
            int k = ic2k(ic,nx,ny);


            for(int iNear; iNear < 27; iNear++)
            {
                for(int q = 0; q < 19; q++)
                {
                    int qicN = qicNear(q, ic, iNear, nx, ny, nz);
                    fNear[iNear][q] = f[qicN];
                }
                cal_rhoUVW(fNear[iNear], &rhoNear[iNear], &uNear[iNear], &vNear[iNear], &wNear[iNear]);
            }
        }

        
        //-- Collision
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

            //-- LES viscosity
            float tau = 1.f/omega;
            
            float PIxx = 0.f;
            float PIxy = 0.f;
            float PIxz = 0.f;
            float PIyy = 0.f;
            float PIyz = 0.f;
            float PIzz = 0.f;

            for(int q = 0; q < 19; q++)
            {
                float uSqr =u*u+v*v+w*w;
                float uDotC = u*cx[q]+v*cy[q]+w*cz[q];
                // float feq = (1.0f+3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rho;
                float feq = (rho+3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q];
                
                PIxx += cx[q]*cx[q]*(ft[q] -feq);
                PIxy += cx[q]*cy[q]*(ft[q] -feq);
                PIxz += cx[q]*cz[q]*(ft[q] -feq);
                PIyy += cy[q]*cy[q]*(ft[q] -feq);
                PIyz += cy[q]*cz[q]*(ft[q] -feq);
                PIzz += cz[q]*cz[q]*(ft[q] -feq);
            }
            float sqrtPIPI = sqrt(PIxx*PIxx+PIyy*PIyy+PIzz*PIzz+2.f*(PIxy*PIxy+PIxz*PIxz+PIyz*PIyz));
            
            float Cs = 0.1f;// 0.1--0.2
            // float Cs = 0.2f;// 0.1--0.2
            // float Cs = 0.33f;// 0.1--0.2

            tauSGS[ic] = LES*0.5f*(-tau +sqrt(tau*tau +18.f*sqrt(2.f)*Cs*Cs*sqrtPIPI/rho));
            // tauSGS[ic] = 3.f*(Cs*Cs)*sqrt(2.f)*sqrtPIPI*0.5f/rho*3.0f/tau;

            //-- Damping of nuSGS (tauSGS) near wall
            const float y = sdf[ic];
            if(y != 10000.f && y > 0.f)
            {
                const float kappa = 0.41f;
                tauSGS[ic] *= pow(min(1.f, kappa*y/Cs),2.f); // Treatment in Fluent

                // //-- van Driest damping function
                // const float Aplus = 26.f;
                // const float Cdelta = 0.158f;
                // const float nu = (tau -0.5f)/3.f;
                // const float E = 9.8f;

                // float yPlus = 1.f;
                // float yPlusN = yPlus;
                // const float epsilon = 0.001f;
                // int iNew = 0;
                // do
                // {
                //     float U = sqrt(u*u+v*v+w*w);
                //     yPlusN = yPlus;
                //     yPlus = ((kappa*U*y/nu) +yPlus)/(1.f +log(E*yPlus));
                //     iNew++;
                // }while(fabs(yPlus -yPlusN) > epsilon && iNew < 10);
                
                // tauSGS[ic] *= pow(min(1.f, (kappa/Cs)*((1.f +1e-10f)-exp(-yPlus/Aplus))*y),2.f);
                // //--
            }

            //-- spongeZone
            int icX0 = index1d(0,ny/2,nz/2,nx,ny);
            int icXE = index1d(nx-1,ny/2,nz/2,nx,ny);
            int icY0 = index1d(nx/2,0,nz/2,nx,ny);
            int icYE = index1d(nx/2,ny-1,nz/2,nx,ny);
            int icZ0 = index1d(nx/2,ny/2,0,nx,ny);
            int icZE = index1d(nx/2,ny/2,nz-1,nx,ny);
            const float spzWidth = 0.1f;
            
            if(boundary1[icX0] == 3)
            {
                if(i < nx*spzWidth)
                {
                    const float fx = i/(nx*spzWidth);
                    tau = (1.f -fx)*(1.f -tauSGS[ic]) +fx*tau;
                }
            }
            if(boundary1[icXE] == 3)
            {
                if(i > nx*(1.f -spzWidth))
                {
                    const float fx = (nx-i)/(nx*spzWidth);
                    tau = (1.f -fx)*(1.f -tauSGS[ic]) +fx*tau;
                }
            }
            if(boundary2[icY0] == 3)
            {
                if(j < ny*spzWidth)
                {
                    const float fx = j/(ny*spzWidth);
                    tau = (1.f -fx)*(1.f -tauSGS[ic]) +fx*tau;
                }
            }
            if(boundary2[icYE] == 3)
            {
                if(j > ny*(1.f -spzWidth))
                {
                    const float fx = (ny-j)/(ny*spzWidth);
                    tau = (1.f -fx)*(1.f -tauSGS[ic]) +fx*tau;
                }
            }
            if(boundary3[icZ0] == 3)
            {
                if(k < nz*spzWidth)
                {
                    const float fx = k/(nz*spzWidth);
                    tau = (1.f -fx)*(1.f -tauSGS[ic]) +fx*tau;
                }
            }
            if(boundary3[icZE] == 3)
            {
                if(k > nz*(1.f -spzWidth))
                {
                    const float fx = (nz-k)/(nz*spzWidth);
                    tau = (1.f -fx)*(1.f -tauSGS[ic]) +fx*tau;
                }
            }
            //--
            
            float omegaEff = 1.f/(tau +tauSGS[ic]);

            const float sqrCs = 1.f/3.f;
            const float quadCs = 1.f/9.f;
            
            //-- BGK model
            // for(int q = 0; q < 19; q++)
            // {
            //     float uSqr =u*u+v*v+w*w;
            //     float uDotC = u*cx[q]+v*cy[q]+w*cz[q];
            //     float feq = (1.0f+3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rho;

            //     int qic = q*elements +ic;

            //     fTmp[qic] = (1.0f -omegaEff)*ft[q] + omegaEff *feq +rho*wt[q]*3.0f*dpdx*cx[q]; // Pull
            // }
            //--

            
            // //-- Cumulant model
            // float K200 = 0.f;
            // float K020 = 0.f;
            // float K002 = 0.f;
            // float K110 = 0.f;
            // float K101 = 0.f;
            // float K011 = 0.f;
            
            // float K210 = 0.f;
            // float K201 = 0.f;
            // float K021 = 0.f;
            // float K120 = 0.f;
            // float K102 = 0.f;
            // float K012 = 0.f;

            // float K220 = 0.f;
            // float K202 = 0.f;
            // float K022 = 0.f;

            // // float Keq200 = sqrCs;
            // // float Keq020 = sqrCs;
            // // float Keq002 = sqrCs;
            // // float Keq110 = 0.f;
            // // float Keq101 = 0.f;
            // // float Keq011 = 0.f;
            
            // // float Keq210 = 0.f;
            // // float Keq201 = 0.f;
            // // float Keq021 = 0.f;
            // // float Keq120 = 0.f;
            // // float Keq102 = 0.f;
            // // float Keq012 = 0.f;

            // // float Keq220 = 0.f;
            // // float Keq202 = 0.f;
            // // float Keq022 = 0.f;

            // // float invRho = 1.f/rho;
            // // // Order 4
            // // K220 = invRho * (ft[8] + ft[10] + ft[7] + ft[9]);
            // // K202 = invRho * (ft[12] + ft[14] + ft[11] + ft[13]);
            // // K022 = invRho * (ft[16] + ft[18] + ft[15] + ft[17]);
            // // // Order 2
            // // K200 = invRho * (ft[2] + ft[1]) + K220 + K202;
            // // K020 = invRho * (ft[4] + ft[3]) + K220 + K022;
            // // K002 = invRho * (ft[6] + ft[5]) + K202 + K022;
        
            // // K110 = K220 - 2.*invRho * (ft[10] + ft[9]);
            // // K101 = K202 - 2.*invRho * (ft[14] + ft[13]);
            // // K011 = K022 - 2.*invRho * (ft[18] + ft[17]);
            // // // Order 3
            // // K210 = K220 - 2.*invRho * (ft[8] + ft[9]);
            // // K201 = K202 - 2.*invRho * (ft[12] + ft[13]);
            // // K021 = K022 - 2.*invRho * (ft[16] + ft[17]);
            // // K120 = K220 - 2.*invRho * (ft[8] + ft[10]);
            // // K102 = K202 - 2.*invRho * (ft[12] + ft[14]);
            // // K012 = K022 - 2.*invRho * (ft[16] + ft[18]);

            // // // Compute central moments from raw moments using binomial formulas
            // // double ux2 = u*u;
            // // double uy2 = v*v;
            // // double uz2 = w*w;
            // // double uxy = u*v;
            // // double uxz = u*w;
            // // double uyz = v*w;

            // // K200 -= ux2;
            // // K020 -= uy2;
            // // K002 -= uz2;
            
            // // K110 -= uxy;
            // // K101 -= uxz;
            // // K011 -= uyz;

            // // K210 -= (v*K200 + 2.*u*K110 + ux2*v);
            // // K201 -= (w*K200 + 2.*u*K101 + ux2*w);
            // // K021 -= (w*K020 + 2.*v*K011 + uy2*w);
            // // K120 -= (u*K020 + 2.*v*K110 + u*uy2);
            // // K102 -= (u*K002 + 2.*w*K101 + u*uz2);
            // // K012 -= (v*K002 + 2.*w*K011 + v*uz2);
            
            // // K220 -= (2.*v*K210 + 2.*u*K120 + uy2*K200 + ux2*K020 + 4.*uxy*K110 + ux2*uy2);
            // // K202 -= (2.*w*K201 + 2.*u*K102 + uz2*K200 + ux2*K002 + 4.*uxz*K101 + ux2*uz2);
            // // K022 -= (2.*w*K021 + 2.*v*K012 + uz2*K020 + uy2*K002 + 4.*uyz*K011 + uy2*uz2);


            // for(int q = 0; q < 19; q++)
            // {
            //     float cxq = cx[q] -u;
            //     float cyq = cy[q] -v;
            //     float czq = cz[q] -w;

            //     K200 += cxq*cxq*ft[q];
            //     K020 += cyq*cyq*ft[q];
            //     K002 += czq*czq*ft[q];
            //     K110 += cxq*cyq*ft[q];
            //     K101 += cxq*czq*ft[q];
            //     K011 += cyq*czq*ft[q];
                
            //     K210 += cxq*cxq*cyq*ft[q];
            //     K201 += cxq*cxq*czq*ft[q];
            //     K021 += cyq*cyq*czq*ft[q];
            //     K120 += cxq*cyq*cyq*ft[q];
            //     K102 += cxq*czq*czq*ft[q];
            //     K012 += cyq*czq*czq*ft[q];

            //     K220 += cxq*cxq*cyq*cyq*ft[q];
            //     K202 += cxq*cxq*czq*czq*ft[q];
            //     K022 += cyq*cyq*czq*czq*ft[q];
            // }

            // float invRho = 1.f/rho;
            // K200 *= invRho; 
            // K020 *= invRho; 
            // K002 *= invRho; 
            // K110 *= invRho; 
            // K101 *= invRho; 
            // K011 *= invRho; 

                
            // K210 *= invRho; 
            // K201 *= invRho; 
            // K021 *= invRho; 
            // K120 *= invRho; 
            // K102 *= invRho; 
            // K012 *= invRho; 

                
            // K220 *= invRho; 
            // K202 *= invRho; 
            // K022 *= invRho;

            // K220 -= (K200*K020 +2.f*K110*K110);
            // K202 -= (K200*K002 +2.f*K101*K101);
            // K022 -= (K020*K002 +2.f*K011*K011);

            // // float omegaB = 1.f;
            // // float omegaM = (omegaB - omegaEff)/3.f;
            // // float omegaP = omegaM +omegaEff;
            // // float Kcoll200 = K200 -omegaP*(K200 -sqrCs) -omegaM*(K020 -sqrCs) -omegaM*(K002 -sqrCs);
            // // float Kcoll020 = K020 -omegaM*(K200 -sqrCs) -omegaP*(K020 -sqrCs) -omegaM*(K002 -sqrCs);
            // // float Kcoll002 = K002 -omegaM*(K200 -sqrCs) -omegaM*(K020 -sqrCs) -omegaP*(K002 -sqrCs);

            // float Kcoll200 = (1.f -omegaEff)*K200 +omegaEff*sqrCs;
            // float Kcoll020 = (1.f -omegaEff)*K020 +omegaEff*sqrCs;
            // float Kcoll002 = (1.f -omegaEff)*K002 +omegaEff*sqrCs;

            // float omega2 = omegaEff;
            // float omega3 = omegaEff;
            // float omega4 = omegaEff;

            // float Kcoll110 = (1.f -omega2)*K110;
            // float Kcoll101 = (1.f -omega2)*K101;
            // float Kcoll011 = (1.f -omega2)*K011;
            
            // float Kcoll210 = (1.f -omega3)*K210;
            // float Kcoll201 = (1.f -omega3)*K201;
            // float Kcoll021 = (1.f -omega3)*K021;
            // float Kcoll120 = (1.f -omega3)*K120;
            // float Kcoll102 = (1.f -omega3)*K102;
            // float Kcoll012 = (1.f -omega3)*K012;

            // float Kcoll220 = (1.f -omega4)*K220;
            // float Kcoll202 = (1.f -omega4)*K202;
            // float Kcoll022 = (1.f -omega4)*K022;
            

            // float CMcoll200 = Kcoll200;
            // float CMcoll020 = Kcoll020;
            // float CMcoll002 = Kcoll002;
            // float CMcoll110 = Kcoll110;
            // float CMcoll101 = Kcoll101;
            // float CMcoll011 = Kcoll011;
        
            // float CMcoll210 = Kcoll210;
            // float CMcoll201 = Kcoll201;
            // float CMcoll021 = Kcoll021;
            // float CMcoll120 = Kcoll120;
            // float CMcoll102 = Kcoll102;
            // float CMcoll012 = Kcoll012;

            // float CMcoll220 = Kcoll220 +Kcoll200*Kcoll020 +2.f*Kcoll110*Kcoll110;
            // float CMcoll202 = Kcoll202 +Kcoll200*Kcoll002 +2.f*Kcoll101*Kcoll101;
            // float CMcoll022 = Kcoll022 +Kcoll020*Kcoll002 +2.f*Kcoll011*Kcoll011;


            // float u2 = u*u;
            // float v2 = v*v;
            // float w2 = w*w;
            // float uv = u*v;
            // float uw = u*w;
            // float vw = v*w;


            // float RMcoll200 = CMcoll200 +u2;
            // float RMcoll020 = CMcoll020 +v2;
            // float RMcoll002 = CMcoll002 +w2;
            // float RMcoll110 = CMcoll110 +uv;
            // float RMcoll101 = CMcoll101 +uw;
            // float RMcoll011 = CMcoll011 +vw;
            
            // float RMcoll210 = CMcoll210 +v*CMcoll200 +2.f*u*CMcoll110 +u2*v;
            // float RMcoll201 = CMcoll201 +w*CMcoll200 +2.f*u*CMcoll101 +u2*w;
            // float RMcoll021 = CMcoll021 +w*CMcoll020 +2.f*v*CMcoll011 +v2*w;
            // float RMcoll120 = CMcoll120 +u*CMcoll020 +2.f*v*CMcoll110 +u*v2;
            // float RMcoll102 = CMcoll102 +u*CMcoll002 +2.f*w*CMcoll101 +u*w2;
            // float RMcoll012 = CMcoll012 +v*CMcoll002 +2.f*w*CMcoll011 +v*w2;

            // float RMcoll220 = CMcoll220 +2.f*v*CMcoll210 +2.f*u*CMcoll120 +v2*CMcoll200 +u2*CMcoll020 +4.f*uv*CMcoll110 +u2*v2;
            // float RMcoll202 = CMcoll202 +2.f*w*CMcoll201 +2.f*u*CMcoll102 +w2*CMcoll200 +u2*CMcoll002 +4.f*uw*CMcoll101 +u2*w2;
            // float RMcoll022 = CMcoll022 +2.f*w*CMcoll021 +2.f*v*CMcoll012 +w2*CMcoll020 +v2*CMcoll002 +4.f*vw*CMcoll011 +v2*w2;


            // fTmp[ 0*elements +ic] = rho*(1.f -RMcoll200 -RMcoll020 -RMcoll002 +RMcoll220 +RMcoll202 +RMcoll022);
            // fTmp[ 1*elements +ic] = 0.5f*rho*(u +RMcoll200 -RMcoll120 -RMcoll102 -RMcoll220 -RMcoll202);
            // fTmp[ 2*elements +ic] = rho*(-u +RMcoll120 +RMcoll102) +0.5f*rho*(u +RMcoll200 -RMcoll120 -RMcoll102 -RMcoll220 -RMcoll202);
            // fTmp[ 3*elements +ic] = 0.5f*rho*(v +RMcoll020 -RMcoll210 -RMcoll012 -RMcoll220 -RMcoll022);
            // fTmp[ 4*elements +ic] = rho*(-v +RMcoll210 +RMcoll012) +0.5f*rho*(v +RMcoll020 -RMcoll210 -RMcoll012 -RMcoll220 -RMcoll022);
            // fTmp[ 5*elements +ic] = 0.5f*rho*(w +RMcoll002 -RMcoll201 -RMcoll021 -RMcoll202 -RMcoll022);
            // fTmp[ 6*elements +ic] = rho*(-w +RMcoll201 +RMcoll021) +0.5f*rho*(w +RMcoll002 -RMcoll201 -RMcoll021 -RMcoll202 -RMcoll022);
            // fTmp[ 7*elements +ic] = 0.25f*rho*(RMcoll110 +RMcoll210 +RMcoll120 +RMcoll220);
            // fTmp[ 8*elements +ic] = 0.5f*rho*(-RMcoll210 -RMcoll120) +0.25f*rho*(RMcoll110 +RMcoll210 +RMcoll120 +RMcoll220);
            // fTmp[ 9*elements +ic] = 0.5f*rho*(-RMcoll110 -RMcoll210) +0.25f*rho*(RMcoll110 +RMcoll210 +RMcoll120 +RMcoll220);
            // fTmp[10*elements +ic] = 0.5f*rho*(-RMcoll110 -RMcoll120) +0.25f*rho*(RMcoll110 +RMcoll210 +RMcoll120 +RMcoll220);
            // fTmp[11*elements +ic] = 0.25f*rho*(RMcoll101 +RMcoll201 +RMcoll102 +RMcoll202);
            // fTmp[12*elements +ic] = 0.5f*rho*(-RMcoll201 -RMcoll102) +0.25f*rho*(RMcoll101 +RMcoll201 +RMcoll102 +RMcoll202);
            // fTmp[13*elements +ic] = 0.5f*rho*(-RMcoll101 -RMcoll201) +0.25f*rho*(RMcoll101 +RMcoll201 +RMcoll102 +RMcoll202);
            // fTmp[14*elements +ic] = 0.5f*rho*(-RMcoll101 -RMcoll102) +0.25f*rho*(RMcoll101 +RMcoll201 +RMcoll102 +RMcoll202);
            // fTmp[15*elements +ic] = 0.25f*rho*(RMcoll011 +RMcoll021 +RMcoll012 +RMcoll022);
            // fTmp[16*elements +ic] = 0.5f*rho*(-RMcoll021 -RMcoll012) +0.25f*rho*(RMcoll011 +RMcoll021 +RMcoll012 +RMcoll022);
            // fTmp[17*elements +ic] = 0.5f*rho*(-RMcoll011 -RMcoll021) +0.25f*rho*(RMcoll011 +RMcoll021 +RMcoll012 +RMcoll022);
            // fTmp[18*elements +ic] = 0.5f*rho*(-RMcoll011 -RMcoll012) +0.25f*rho*(RMcoll011 +RMcoll021 +RMcoll012 +RMcoll022);
            //--

            
            //-- Recursive-regularized model
            float RR200 = 0.f;
            float RR020 = 0.f;
            float RR002 = 0.f;
            float RR110 = 0.f;
            float RR101 = 0.f;
            float RR011 = 0.f;
            
            float RReq200 = u*u;
            float RReq020 = v*v;
            float RReq002 = w*w;
            float RReq110 = u*v;
            float RReq101 = u*w;
            float RReq011 = v*w;

            float RReq210 = RReq200*v;
            float RReq201 = RReq200*w;
            float RReq021 = RReq020*w;
            float RReq120 = RReq020*u;
            float RReq102 = RReq002*u;
            float RReq012 = RReq002*v;

            float RReq220 = RReq200*RReq020;
            float RReq202 = RReq200*RReq002;
            float RReq022 = RReq020*RReq002;
            
            
            for(int q = 0; q < 19; q++)
            {
                float Hxx = cx[q]*cx[q] -sqrCs;
                float Hyy = cy[q]*cy[q] -sqrCs;
                float Hzz = cz[q]*cz[q] -sqrCs;

                RR200 += Hxx*ft[q];
                RR020 += Hyy*ft[q];
                RR002 += Hzz*ft[q];
                RR110 += cx[q]*cy[q]*ft[q];
                RR101 += cx[q]*cz[q]*ft[q];
                RR011 += cy[q]*cz[q]*ft[q];
            }
            float invRho = 1.f/rho;
            RR200 *= invRho;
            RR020 *= invRho;
            RR002 *= invRho;
            RR110 *= invRho;
            RR101 *= invRho;
            RR011 *= invRho;

            float RRneq200 = RR200 -RReq200;
            float RRneq020 = RR020 -RReq020;
            float RRneq002 = RR002 -RReq002;
            float RRneq110 = RR110 -RReq110;
            float RRneq101 = RR101 -RReq101;
            float RRneq011 = RR011 -RReq011;

            float RRneq210 = v*RRneq200 +2.f*u*RRneq110;
            float RRneq201 = w*RRneq200 +2.f*u*RRneq101;
            float RRneq021 = w*RRneq020 +2.f*v*RRneq011;
            float RRneq120 = u*RRneq020 +2.f*v*RRneq110;
            float RRneq102 = u*RRneq002 +2.f*w*RRneq101;
            float RRneq012 = v*RRneq002 +2.f*w*RRneq011;

            float RRneq220 = RReq020*RRneq200 +RReq200*RRneq020 +4.f*RReq110*RRneq110;
            float RRneq202 = RReq002*RRneq200 +RReq200*RRneq002 +4.f*RReq101*RRneq101;
            float RRneq022 = RReq002*RRneq020 +RReq020*RRneq002 +4.f*RReq011*RRneq011;

            // float omegaB = 0.985f;
            float omegaB = 1.f;
            float omegaM = (omegaB - omegaEff)/3.f;
            float omegaP = omegaM +omegaEff;

            float RRcoll200 = RR200 -omegaP*RRneq200 -omegaM*RRneq020 -omegaM*RRneq002;
            float RRcoll020 = RR020 -omegaM*RRneq200 -omegaP*RRneq020 -omegaM*RRneq002;
            float RRcoll002 = RR002 -omegaM*RRneq200 -omegaM*RRneq020 -omegaP*RRneq002;


            // float RRcoll200 = RR200 -omegaEff*RRneq200;
            // float RRcoll020 = RR020 -omegaEff*RRneq020;
            // float RRcoll002 = RR002 -omegaEff*RRneq002;

            float omega2 = omegaEff;
            float omega3 = omegaEff;
            float omega4 = omegaEff;

            float RRcoll110 = RR110 -omega2*RRneq110;
            float RRcoll101 = RR101 -omega2*RRneq101;
            float RRcoll011 = RR011 -omega2*RRneq011;

            float RRcoll210 = RReq210 +(1.f -omega3)*RRneq210;
            float RRcoll201 = RReq201 +(1.f -omega3)*RRneq201;
            float RRcoll021 = RReq021 +(1.f -omega3)*RRneq021;
            float RRcoll120 = RReq120 +(1.f -omega3)*RRneq120;
            float RRcoll102 = RReq102 +(1.f -omega3)*RRneq102;
            float RRcoll012 = RReq012 +(1.f -omega3)*RRneq012;

            float RRcoll220 = RReq220 +(1.f -omega4)*RRneq220;
            float RRcoll202 = RReq202 +(1.f -omega4)*RRneq202;
            float RRcoll022 = RReq022 +(1.f -omega4)*RRneq022;

            float RMcoll200 = RRcoll200 +sqrCs;
            float RMcoll020 = RRcoll020 +sqrCs;
            float RMcoll002 = RRcoll002 +sqrCs;
            float RMcoll110 = RRcoll110;
            float RMcoll101 = RRcoll101;
            float RMcoll011 = RRcoll011;

            float RMcoll210 = RRcoll210 +sqrCs*v;
            float RMcoll201 = RRcoll201 +sqrCs*w;
            float RMcoll021 = RRcoll021 +sqrCs*w;
            float RMcoll120 = RRcoll120 +sqrCs*u;
            float RMcoll102 = RRcoll102 +sqrCs*u;
            float RMcoll012 = RRcoll012 +sqrCs*v;

            float RMcoll220 = RRcoll220 +sqrCs*(RRcoll200 +RRcoll020) +quadCs;
            float RMcoll202 = RRcoll202 +sqrCs*(RRcoll200 +RRcoll002) +quadCs;
            float RMcoll022 = RRcoll022 +sqrCs*(RRcoll020 +RRcoll002) +quadCs;

            fTmp[ 0*elements +ic] = rho*(1.f -RMcoll200 -RMcoll020 -RMcoll002 +RMcoll220 +RMcoll202 +RMcoll022) +rho*wt[0]*3.0f*dpdx*cx[0];
            fTmp[ 1*elements +ic] = 0.5f*rho*(u +RMcoll200 -RMcoll120 -RMcoll102 -RMcoll220 -RMcoll202) +rho*wt[1]*3.0f*dpdx*cx[1];
            fTmp[ 2*elements +ic] = rho*(-u +RMcoll120 +RMcoll102) +0.5f*rho*(u +RMcoll200 -RMcoll120 -RMcoll102 -RMcoll220 -RMcoll202) +rho*wt[2]*3.0f*dpdx*cx[2];
            fTmp[ 3*elements +ic] = 0.5f*rho*(v +RMcoll020 -RMcoll210 -RMcoll012 -RMcoll220 -RMcoll022) +rho*wt[3]*3.0f*dpdx*cx[3];
            fTmp[ 4*elements +ic] = rho*(-v +RMcoll210 +RMcoll012) +0.5f*rho*(v +RMcoll020 -RMcoll210 -RMcoll012 -RMcoll220 -RMcoll022) +rho*wt[4]*3.0f*dpdx*cx[4];
            fTmp[ 5*elements +ic] = 0.5f*rho*(w +RMcoll002 -RMcoll201 -RMcoll021 -RMcoll202 -RMcoll022) +rho*wt[5]*3.0f*dpdx*cx[5];
            fTmp[ 6*elements +ic] = rho*(-w +RMcoll201 +RMcoll021) +0.5f*rho*(w +RMcoll002 -RMcoll201 -RMcoll021 -RMcoll202 -RMcoll022) +rho*wt[6]*3.0f*dpdx*cx[4];
            fTmp[ 7*elements +ic] = 0.25f*rho*(RMcoll110 +RMcoll210 +RMcoll120 +RMcoll220) +rho*wt[7]*3.0f*dpdx*cx[7];
            fTmp[ 8*elements +ic] = 0.5f*rho*(-RMcoll210 -RMcoll120) +0.25f*rho*(RMcoll110 +RMcoll210 +RMcoll120 +RMcoll220) +rho*wt[8]*3.0f*dpdx*cx[8];
            fTmp[ 9*elements +ic] = 0.5f*rho*(-RMcoll110 -RMcoll210) +0.25f*rho*(RMcoll110 +RMcoll210 +RMcoll120 +RMcoll220) +rho*wt[9]*3.0f*dpdx*cx[9];
            fTmp[10*elements +ic] = 0.5f*rho*(-RMcoll110 -RMcoll120) +0.25f*rho*(RMcoll110 +RMcoll210 +RMcoll120 +RMcoll220) +rho*wt[10]*3.0f*dpdx*cx[10];
            fTmp[11*elements +ic] = 0.25f*rho*(RMcoll101 +RMcoll201 +RMcoll102 +RMcoll202) +rho*wt[11]*3.0f*dpdx*cx[11];
            fTmp[12*elements +ic] = 0.5f*rho*(-RMcoll201 -RMcoll102) +0.25f*rho*(RMcoll101 +RMcoll201 +RMcoll102 +RMcoll202) +rho*wt[12]*3.0f*dpdx*cx[12];
            fTmp[13*elements +ic] = 0.5f*rho*(-RMcoll101 -RMcoll201) +0.25f*rho*(RMcoll101 +RMcoll201 +RMcoll102 +RMcoll202) +rho*wt[13]*3.0f*dpdx*cx[13];
            fTmp[14*elements +ic] = 0.5f*rho*(-RMcoll101 -RMcoll102) +0.25f*rho*(RMcoll101 +RMcoll201 +RMcoll102 +RMcoll202) +rho*wt[14]*3.0f*dpdx*cx[14];
            fTmp[15*elements +ic] = 0.25f*rho*(RMcoll011 +RMcoll021 +RMcoll012 +RMcoll022) +rho*wt[15]*3.0f*dpdx*cx[15];
            fTmp[16*elements +ic] = 0.5f*rho*(-RMcoll021 -RMcoll012) +0.25f*rho*(RMcoll011 +RMcoll021 +RMcoll012 +RMcoll022) +rho*wt[16]*3.0f*dpdx*cx[16];
            fTmp[17*elements +ic] = 0.5f*rho*(-RMcoll011 -RMcoll021) +0.25f*rho*(RMcoll011 +RMcoll021 +RMcoll012 +RMcoll022) +rho*wt[17]*3.0f*dpdx*cx[17];
            fTmp[18*elements +ic] = 0.5f*rho*(-RMcoll011 -RMcoll012) +0.25f*rho*(RMcoll011 +RMcoll021 +RMcoll012 +RMcoll022) +rho*wt[18]*3.0f*dpdx*cx[18];
            //--
        }

        //-- Equilibrium Boundary
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
                    for(int q = 0; q < 19; q++)
                    {
                        if(q == 1 || q == 7 || q == 9 || q == 11 || q == 13)
                        {                            
                            float uSqr = u*u +v*v +w*w;
                            float uDotC = u*cx[q]+v*cy[q]+w*cz[q];

                            ft[q] = (1.f +3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rhow;
                        }
                    }
                }
            }
            if(i == nx-1)
            {
                if(boundary1[ic] == 2)
                {   for(int q = 0; q < 19; q++)
                    {
                        if(q == 2 || q == 8 || q == 10 || q == 12 || q == 14)
                        {
                            float uSqr = u*u +v*v +w*w;
                            float uDotC = u*cx[q]+v*cy[q]+w*cz[q];

                            ft[q] = (1.f +3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rhow;
                        }
                    }
                }
            }
            if(j == 0)
            {
                if(boundary2[ic] == 2)
                {
                    for(int q = 0; q < 19; q++)
                    {
                        if(q == 3 || q == 7 || q == 10 || q == 15 || q == 17)
                        {
                            float uSqr = u*u +v*v +w*w;
                            float uDotC = u*cx[q]+v*cy[q]+w*cz[q];

                            ft[q] = (1.f +3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rhow;
                        }
                    }
                }
            }
            if(j == ny-1)
            {
                if(boundary2[ic] == 2)
                {
                    for(int q = 0; q < 19; q++)
                    {
                        if(q == 4 || q == 8 || q == 9 || q == 16 || q == 18)
                        {
                            float uSqr = u*u +v*v +w*w;
                            float uDotC = u*cx[q]+v*cy[q]+w*cz[q];

                            ft[q] = (1.f +3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rhow;
                        }
                    }
                }
            }
            if(k == 0)
            {
                if(boundary3[ic] == 2)
                {
                    for(int q = 0; q < 19; q++)
                    {
                        if(q == 5 || q == 11 || q == 14 || q == 15 || q == 18)
                        {
                            float uSqr = u*u +v*v +w*w;
                            float uDotC = u*cx[q]+v*cy[q]+w*cz[q];

                            ft[q] = (1.f +3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rhow;
                        }
                    }
                }
            }
            if(k == nz-1)
            {
                if(boundary3[ic] == 2)
                {
                    for(int q = 0; q < 19; q++)
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
        //--
    }
}
