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

inline int upwindID(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz)
{
    if(q == 0)
    {
        return index1d(i,j,k,nx,ny);
    }
    else if(q == 1)
    {
        return i != 0 ? index1d(i-1,j,k,nx,ny) : index1d(nx-1,j,k,nx,ny);
    }
    else if(q == 2)
    {
        return i != nx-1 ? index1d(i+1,j,k,nx,ny) : index1d(0,j,k,nx,ny);
    }
    else if(q == 3)
    {
        return j != 0 ? index1d(i,j-1,k,nx,ny) : index1d(i,ny-1,k,nx,ny);
    }
    else if(q == 4)
    {
        return j != ny-1 ? index1d(i,j+1,k,nx,ny) : index1d(i,0,k,nx,ny);
    }
    else if(q == 5)
    {
        return k != 0 ? index1d(i,j,k-1,nx,ny) : index1d(i,j,nz-1,nx,ny);
    }
    else if(q == 6)
    {
        return k != nz-1 ? index1d(i,j,k+1,nx,ny) : index1d(i,j,0,nx,ny);
    }
    else if(q == 7)
    {
        return (i != 0 && j != 0) ? index1d(i-1, j-1, k, nx, ny) :
               (i == 0 && j != 0) ? index1d(nx-1, j-1, k, nx, ny) :
               (i != 0 && j == 0) ? index1d(i-1, ny-1, k, nx, ny) :
               index1d(nx-1, ny-1, k, nx, ny);
    }
    else if(q == 8)
    {
        return (i != nx-1 && j != ny-1) ? index1d(i+1, j+1, k, nx, ny) :
               (i == nx-1 && j != ny-1) ? index1d(0, j+1, k, nx, ny) :
               (i != nx-1 && j == ny-1) ? index1d(i+1, 0, k, nx, ny) :
               index1d(0, 0, k, nx, ny);
    }
    else if(q == 9)
    {
        return (i != 0 && j != ny-1) ? index1d(i-1, j+1, k, nx, ny) :
               (i == 0 && j != ny-1) ? index1d(nx-1, j+1, k, nx, ny) :
               (i != 0 && j == ny-1) ? index1d(i-1, 0, k, nx, ny) :
               index1d(nx-1, 0, k, nx, ny);
    }
    else if(q == 10)
    {
        return (i != nx-1 && j != 0) ? index1d(i+1, j-1, k, nx, ny) :
               (i == nx-1 && j != 0) ? index1d(0, j-1, k, nx, ny) :
               (i != nx-1 && j == 0) ? index1d(i+1, ny-1, k, nx, ny) :
               index1d(0, ny-1, k, nx, ny);
    }
    else if(q == 11)
    {
        return (i != 0 && k != 0) ? index1d(i-1, j, k-1, nx, ny) :
               (i == 0 && k != 0) ? index1d(nx-1, j, k-1, nx, ny) :
               (i != 0 && k == 0) ? index1d(i-1, j, nz-1, nx, ny) :
               index1d(nx-1, j, nz-1, nx, ny);
    }
    else if(q == 12)
    {
        return (i != nx-1 && k != nz-1) ? index1d(i+1, j, k+1, nx, ny) :
               (i == nx-1 && k != nz-1) ? index1d(0, j, k+1, nx, ny) :
               (i != nx-1 && k == nz-1) ? index1d(i+1, j, 0, nx, ny) :
               index1d(0, j, 0, nx, ny);
    }
    else if(q == 13)
    {
        return (i != 0 && k != nz-1) ? index1d(i-1, j, k+1, nx, ny) :
               (i == 0 && k != nz-1) ? index1d(nx-1, j, k+1, nx, ny) :
               (i != 0 && k == nz-1) ? index1d(i-1, j, 0, nx, ny) :
               index1d(nx-1, j, 0, nx, ny);
    }
    else if(q == 14)
    {
        return (i != nx-1 && k != 0) ? index1d(i+1, j, k-1, nx, ny) :
               (i == nx-1 && k != 0) ? index1d(0, j, k-1, nx, ny) :
               (i != nx-1 && k == 0) ? index1d(i+1, j, nz-1, nx, ny) :
               index1d(0, j, nz-1, nx, ny);
    }
    else if(q == 15)
    {
        return (j != 0 && k != 0) ? index1d(i, j-1, k-1, nx, ny) :
               (j == 0 && k != 0) ? index1d(i, ny-1, k-1, nx, ny) :
               (j != 0 && k == 0) ? index1d(i, j-1, nz-1, nx, ny) :
               index1d(i, ny-1, nz-1, nx, ny);
    }
    else if(q == 16)
    {
        return (j != ny-1 && k != nz-1) ? index1d(i, j+1, k+1, nx, ny) :
               (j == ny-1 && k != nz-1) ? index1d(i, 0, k+1, nx, ny) :
               (j != ny-1 && k == nz-1) ? index1d(i, j+1, 0, nx, ny) :
               index1d(i, 0, 0, nx, ny);
    }
    else if(q == 17)
    {
        return (j != 0 && k != nz-1) ? index1d(i, j-1, k+1, nx, ny) :
               (j == 0 && k != nz-1) ? index1d(i, ny-1, k+1, nx, ny) :
               (j != 0 && k == nz-1) ? index1d(i, j-1, 0, nx, ny) :
               index1d(i, ny-1, 0, nx, ny);
    }
    else if(q == 18)
    {
        return (j != ny-1 && k != 0) ? index1d(i, j+1, k-1, nx, ny) :
               (j == ny-1 && k != 0) ? index1d(i, 0, k-1, nx, ny) :
               (j != ny-1 && k == 0) ? index1d(i, j+1, nz-1, nx, ny) :
               index1d(i, 0, nz-1, nx, ny);
    }
    else
    {
        return 0;
    }
}

inline int upwindID_B(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz, const int* normal)
{
    int ic = index1d(i, j, k, nx, ny);
    if(q == 0)
    {
        return index1d(i,j,k,nx,ny);
    }
    else if(q == 1)
    {
        return i != 0 ? index1d(i-1,j,k,nx,ny) : (normal[ic] == 0 ? index1d(nx-1,j,k,nx,ny) : -1);
    }
    else if(q == 2)
    {
        return i != nx-1 ? index1d(i+1,j,k,nx,ny) : (normal[ic] == 0 ? index1d(0,j,k,nx,ny) : -1);
    }
    else if(q == 3)
    {
        return j != 0 ? index1d(i,j-1,k,nx,ny) : (normal[ic] == 0 ? index1d(i,ny-1,k,nx,ny) : -1);
    }
    else if(q == 4)
    {
        return j != ny-1 ? index1d(i,j+1,k,nx,ny) : (normal[ic] == 0 ? index1d(i,0,k,nx,ny) : -1);
    }
    else if(q == 5)
    {
        return k != 0 ? index1d(i,j,k-1,nx,ny) : (normal[ic] == 0 ? index1d(i,j,nz-1,nx,ny) : -1);
    }
    else if(q == 6)
    {
        return k != nz-1 ? index1d(i,j,k+1,nx,ny) : (normal[ic] == 0 ? index1d(i,j,0,nx,ny) : -1);
    }
    else if(q == 7)
    {
        return (i != 0 && j != 0) ? index1d(i-1, j-1, k, nx, ny) :
               (i == 0 && j != 0) ? (normal[ic] == 0 ? index1d(nx-1, j-1, k, nx, ny) : -1) :
               (i != 0 && j == 0) ? (normal[ic] == 0 ? index1d(i-1, ny-1, k, nx, ny) : -1) :
               (normal[ic] == 0 ? index1d(nx-1, ny-1, k, nx, ny) : -1);
    }
    else if(q == 8)
    {
        return (i != nx-1 && j != ny-1) ? index1d(i+1, j+1, k, nx, ny) :
               (i == nx-1 && j != ny-1) ? (normal[ic] == 0 ? index1d(0, j+1, k, nx, ny) : -1) :
               (i != nx-1 && j == ny-1) ? (normal[ic] == 0 ? index1d(i+1, 0, k, nx, ny) : -1) :
               (normal[ic] == 0 ? index1d(0, 0, k, nx, ny) : -1);
    }
    else if(q == 9)
    {
        return (i != 0 && j != ny-1) ? index1d(i-1, j+1, k, nx, ny) :
               (i == 0 && j != ny-1) ? (normal[ic] == 0 ? index1d(nx-1, j+1, k, nx, ny) : -1):
               (i != 0 && j == ny-1) ? (normal[ic] == 0 ? index1d(i-1, 0, k, nx, ny) : -1) :
               (normal[ic] == 0 ? index1d(nx-1, 0, k, nx, ny) : -1);
    }
    else if(q == 10)
    {
        return (i != nx-1 && j != 0) ? index1d(i+1, j-1, k, nx, ny) :
               (i == nx-1 && j != 0) ? (normal[ic] == 0 ? index1d(0, j-1, k, nx, ny) : -1) :
               (i != nx-1 && j == 0) ? (normal[ic] == 0 ? index1d(i+1, ny-1, k, nx, ny) : -1) :
               (normal[ic] == 0 ? index1d(0, ny-1, k, nx, ny) : -1);
    }
    else if(q == 11)
    {
        return (i != 0 && k != 0) ? index1d(i-1, j, k-1, nx, ny) :
               (i == 0 && k != 0) ? (normal[ic] == 0 ? index1d(nx-1, j, k-1, nx, ny) : -1) :
               (i != 0 && k == 0) ? (normal[ic] == 0 ? index1d(i-1, j, nz-1, nx, ny) : -1) :
               (normal[ic] == 0 ? index1d(nx-1, j, nz-1, nx, ny) : -1);
    }
    else if(q == 12)
    {
        return (i != nx-1 && k != nz-1) ? index1d(i+1, j, k+1, nx, ny) :
               (i == nx-1 && k != nz-1) ? (normal[ic] == 0 ? index1d(0, j, k+1, nx, ny) : -1) :
               (i != nx-1 && k == nz-1) ? (normal[ic] == 0 ? index1d(i+1, j, 0, nx, ny) : -1) :
               (normal[ic] == 0 ? index1d(0, j, 0, nx, ny) : -1);
    }
    else if(q == 13)
    {
        return (i != 0 && k != nz-1) ? index1d(i-1, j, k+1, nx, ny) :
               (i == 0 && k != nz-1) ? (normal[ic] == 0 ? index1d(nx-1, j, k+1, nx, ny) : -1) :
               (i != 0 && k == nz-1) ? (normal[ic] == 0 ? index1d(i-1, j, 0, nx, ny) : -1) :
               (normal[ic] == 0 ? index1d(nx-1, j, 0, nx, ny) : -1);
    }
    else if(q == 14)
    {
        return (i != nx-1 && k != 0) ? index1d(i+1, j, k-1, nx, ny) :
               (i == nx-1 && k != 0) ? (normal[ic] == 0 ? index1d(0, j, k-1, nx, ny) : -1) :
               (i != nx-1 && k == 0) ? (normal[ic] == 0 ? index1d(i+1, j, nz-1, nx, ny) : -1) :
               (normal[ic] == 0 ? index1d(0, j, nz-1, nx, ny) : -1);
    }
    else if(q == 15)
    {
        return (j != 0 && k != 0) ? index1d(i, j-1, k-1, nx, ny) :
               (j == 0 && k != 0) ? (normal[ic] == 0 ? index1d(i, ny-1, k-1, nx, ny) : -1) :
               (j != 0 && k == 0) ? (normal[ic] == 0 ? index1d(i, j-1, nz-1, nx, ny) : -1) :
               (normal[ic] == 0 ? index1d(i, ny-1, nz-1, nx, ny) : -1);
    }
    else if(q == 16)
    {
        return (j != ny-1 && k != nz-1) ? index1d(i, j+1, k+1, nx, ny) :
               (j == ny-1 && k != nz-1) ? (normal[ic] == 0 ? index1d(i, 0, k+1, nx, ny) : -1) :
               (j != ny-1 && k == nz-1) ? (normal[ic] == 0 ? index1d(i, j+1, 0, nx, ny) : -1) :
               (normal[ic] == 0 ? index1d(i, 0, 0, nx, ny) : -1);
    }
    else if(q == 17)
    {
        return (j != 0 && k != nz-1) ? index1d(i, j-1, k+1, nx, ny) :
               (j == 0 && k != nz-1) ? (normal[ic] == 0 ? index1d(i, ny-1, k+1, nx, ny) : -1) :
               (j != 0 && k == nz-1) ? (normal[ic] == 0 ? index1d(i, j-1, 0, nx, ny) : -1) :
               (normal[ic] == 0 ? index1d(i, ny-1, 0, nx, ny) : -1);
    }
    else if(q == 18)
    {
        return (j != ny-1 && k != 0) ? index1d(i, j+1, k-1, nx, ny) :
               (j == ny-1 && k != 0) ? (normal[ic] == 0 ? index1d(i, 0, k-1, nx, ny) : -1) :
               (j != ny-1 && k == 0) ? (normal[ic] == 0 ? index1d(i, j+1, nz-1, nx, ny) : -1) :
               (normal[ic] == 0 ? index1d(i, 0, nz-1, nx, ny) : -1);
    }
    else
    {
        return 0;
    }
}

inline int downwindID(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz)
{
    if(q == 0)
    {
        return index1d(i,j,k,nx,ny);
    }
    else if(q == 2)
    {
        return i != 0 ? index1d(i-1,j,k,nx,ny) : index1d(nx-1,j,k,nx,ny);
    }
    else if(q == 1)
    {
        return i != nx-1 ? index1d(i+1,j,k,nx,ny) : index1d(0,j,k,nx,ny);
    }
    else if(q == 4)
    {
        return j != 0 ? index1d(i,j-1,k,nx,ny) : index1d(i,ny-1,k,nx,ny);
    }
    else if(q == 3)
    {
        return j != ny-1 ? index1d(i,j+1,k,nx,ny) : index1d(i,0,k,nx,ny);
    }
    else if(q == 6)
    {
        return k != 0 ? index1d(i,j,k-1,nx,ny) : index1d(i,j,nz-1,nx,ny);
    }
    else if(q == 5)
    {
        return k != nz-1 ? index1d(i,j,k+1,nx,ny) : index1d(i,j,0,nx,ny);
    }
    else if(q == 8)
    {
        return (i != 0 && j != 0) ? index1d(i-1, j-1, k, nx, ny) :
               (i == 0 && j != 0) ? index1d(nx-1, j-1, k, nx, ny) :
               (i != 0 && j == 0) ? index1d(i-1, ny-1, k, nx, ny) :
               index1d(nx-1, ny-1, k, nx, ny);
    }
    else if(q == 7)
    {
        return (i != nx-1 && j != ny-1) ? index1d(i+1, j+1, k, nx, ny) :
               (i == nx-1 && j != ny-1) ? index1d(0, j+1, k, nx, ny) :
               (i != nx-1 && j == ny-1) ? index1d(i+1, 0, k, nx, ny) :
               index1d(0, 0, k, nx, ny);
    }
    else if(q == 10)
    {
        return (i != 0 && j != ny-1) ? index1d(i-1, j+1, k, nx, ny) :
               (i == 0 && j != ny-1) ? index1d(nx-1, j+1, k, nx, ny) :
               (i != 0 && j == ny-1) ? index1d(i-1, 0, k, nx, ny) :
               index1d(nx-1, 0, k, nx, ny);
    }
    else if(q == 9)
    {
        return (i != nx-1 && j != 0) ? index1d(i+1, j-1, k, nx, ny) :
               (i == nx-1 && j != 0) ? index1d(0, j-1, k, nx, ny) :
               (i != nx-1 && j == 0) ? index1d(i+1, ny-1, k, nx, ny) :
               index1d(0, ny-1, k, nx, ny);
    }
    else if(q == 12)
    {
        return (i != 0 && k != 0) ? index1d(i-1, j, k-1, nx, ny) :
               (i == 0 && k != 0) ? index1d(nx-1, j, k-1, nx, ny) :
               (i != 0 && k == 0) ? index1d(i-1, j, nz-1, nx, ny) :
               index1d(nx-1, j, nz-1, nx, ny);
    }
    else if(q == 11)
    {
        return (i != nx-1 && k != nz-1) ? index1d(i+1, j, k+1, nx, ny) :
               (i == nx-1 && k != nz-1) ? index1d(0, j, k+1, nx, ny) :
               (i != nx-1 && k == nz-1) ? index1d(i+1, j, 0, nx, ny) :
               index1d(0, j, 0, nx, ny);
    }
    else if(q == 14)
    {
        return (i != 0 && k != nz-1) ? index1d(i-1, j, k+1, nx, ny) :
               (i == 0 && k != nz-1) ? index1d(nx-1, j, k+1, nx, ny) :
               (i != 0 && k == nz-1) ? index1d(i-1, j, 0, nx, ny) :
               index1d(nx-1, j, 0, nx, ny);
    }
    else if(q == 13)
    {
        return (i != nx-1 && k != 0) ? index1d(i+1, j, k-1, nx, ny) :
               (i == nx-1 && k != 0) ? index1d(0, j, k-1, nx, ny) :
               (i != nx-1 && k == 0) ? index1d(i+1, j, nz-1, nx, ny) :
               index1d(0, j, nz-1, nx, ny);
    }
    else if(q == 16)
    {
        return (j != 0 && k != 0) ? index1d(i, j-1, k-1, nx, ny) :
               (j == 0 && k != 0) ? index1d(i, ny-1, k-1, nx, ny) :
               (j != 0 && k == 0) ? index1d(i, j-1, nz-1, nx, ny) :
               index1d(i, ny-1, nz-1, nx, ny);
    }
    else if(q == 15)
    {
        return (j != ny-1 && k != nz-1) ? index1d(i, j+1, k+1, nx, ny) :
               (j == ny-1 && k != nz-1) ? index1d(i, 0, k+1, nx, ny) :
               (j != ny-1 && k == nz-1) ? index1d(i, j+1, 0, nx, ny) :
               index1d(i, 0, 0, nx, ny);
    }
    else if(q == 18)
    {
        return (j != 0 && k != nz-1) ? index1d(i, j-1, k+1, nx, ny) :
               (j == 0 && k != nz-1) ? index1d(i, ny-1, k+1, nx, ny) :
               (j != 0 && k == nz-1) ? index1d(i, j-1, 0, nx, ny) :
               index1d(i, ny-1, 0, nx, ny);
    }
    else if(q == 17)
    {
        return (j != ny-1 && k != 0) ? index1d(i, j+1, k-1, nx, ny) :
               (j == ny-1 && k != 0) ? index1d(i, 0, k-1, nx, ny) :
               (j != ny-1 && k == 0) ? index1d(i, j+1, nz-1, nx, ny) :
               index1d(i, 0, nz-1, nx, ny);
    }
    else
    {
        return 0;
    }
}
// ---

__kernel void k_collision
(
   __global float* f, __global float* fTmp,

   const unsigned elements,
   const float omega
)
{
    int i = get_global_id(0);

    float wt[19] = {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    float cx[19] = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};
    float cy[19] = {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0};
    float cz[19] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0};

    float rho = 0.0;
    float u = 0.0;
    float v = 0.0;
    float w = 0.0;
    for(int q = 0; q < 19; q++)
    {
        int qic = q*elements +i;

        rho += f[qic];

        u += f[qic]*cx[q];
        v += f[qic]*cy[q];
        w += f[qic]*cz[q];
    }
    u /= rho;
    v /= rho;
    w /= rho;

    for(int q = 0; q < 19; q++)
    {
        float uSqr =u*u+v*v+w*w;
        float uDotC = u*cx[q]+v*cy[q]+w*cz[q];
        float feq = (1.0+3.0*uDotC +4.5*uDotC*uDotC -1.5*uSqr)*wt[q]*rho;

        int qic = q*elements +i;

        fTmp[qic] = f[qic] - omega *(f[qic] -feq);
        // printf("feq: %.1f", feq);
    }
}

__kernel void k_externalForce
(
   __global float* f, __global float* fTmp,
   const unsigned elements,
   const float dpdx
)
{
    int i = get_global_id(0);

    float wt[19] = {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    float cx[19] = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};

    float rho = 0.0;
    for(int q = 0; q < 19; q++)
    {
        int qic = q*elements +i;

        rho += f[qic];
    }

    for(int q = 0; q < 19; q++)
    {
        int qic = q*elements +i;

        fTmp[qic] += rho*wt[q]*3.0f*dpdx*cx[q];
    }
}

__kernel void k_streaming
(
   __global float* f, __global float* fTmp,
   __global unsigned* upQID
)
{
    int qi = get_global_id(0);
    f[qi] = fTmp[upQID[qi]];
}

__kernel void k_collisionStreaming // Push
(
   __global float* f, __global float* fTmp,
   const unsigned elements,
   const float omega,
   const float dpdx,
   const int nx, const int ny, const int nz
)
{
    int ic = get_global_id(0);

    float wt[19] = {1.0f/3.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f};
    float cx[19] = {0.0f, 1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    float cy[19] = {0.0f, 0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f, 1.0f, -1.0f, -1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -1.0f, 1.0f, -1.0f};
    float cz[19] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -1.0f, -1.0f, 1.0f, 1.0f, -1.0f, -1.0f, 1.0f};

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

    for(int q = 0; q < 19; q++)
    {
        float uSqr =u*u+v*v+w*w;
        float uDotC = u*cx[q]+v*cy[q]+w*cz[q];
        float feq = (1.0f+3.0f*uDotC +4.5f*uDotC*uDotC -1.5f*uSqr)*wt[q]*rho;

        int qic = q*elements +ic;
        int i = ic2i(ic,nx,ny);
        int j = ic2j(ic,nx,ny);
        int k = ic2k(ic,nx,ny);

        int downID = downwindID(q,i,j,k,nx,ny,nz);
        int downQID = idf(q, downID, nx, ny, nz);

        fTmp[downQID] = (1.0f -omega)*f[qic] + omega *feq +rho*wt[q]*3.0f*dpdx*cx[q]; // Push
    }
}

__kernel void k_bounceBack
(
   __global float* f,
   __global int* BBID, __global int* BBQID,
   const unsigned elements, const unsigned nBB
)
{
    int i = get_global_id(0);
    for(int q = 0; q < 19; q++)
    {
        int qic = q*elements +BBID[i];
        int qid = q*nBB +i;
        if(BBQID[qid] >= 0)
        {
            f[qic] = f[BBQID[qid]];
        }
        else
        {
            f[qic] = 1.0f/36.0f;
        }
    }   
}

__kernel void k_bounceBackMovingWall
(
    __global float* f,
    __global int* BBmovingWallID,
    __global float* u0, __global float* v0, __global float* w0,
    __global int* normal,
    const unsigned elements
)
{
    int i = get_global_id(0);
    int ic = BBmovingWallID[i];
    float u = u0[ic];
    float v = v0[ic];
    float w = w0[ic];
    int nVec = normal[ic];  

    float* Vc    = &f[0*elements +ic];//(0)
    float* Vin   = &f[1*elements +ic];//(+x)
    float* Vip   = &f[2*elements +ic];//(-x)
    float* Vjn   = &f[3*elements +ic];//(+y)
    float* Vjp   = &f[4*elements +ic];//(-y)
    float* Vkn   = &f[5*elements +ic];//(+z)
    float* Vkp   = &f[6*elements +ic];//(-z)
    float* Vinjn = &f[7*elements +ic];//(+x,+y)
    float* Vipjp = &f[8*elements +ic];//(-x,-y)
    float* Vinjp = &f[9*elements +ic];//(+x,-y)
    float* Vipjn = &f[10*elements +ic];//(-x,+y)
    float* Vinkn = &f[11*elements +ic];//(+x,+z)
    float* Vipkp = &f[12*elements +ic];//(-x,-z)
    float* Vinkp = &f[13*elements +ic];//(+x,-z)
    float* Vipkn = &f[14*elements +ic];//(-x,+z)
    float* Vjnkn = &f[15*elements +ic];//(+y,+z)
    float* Vjpkp = &f[16*elements +ic];//(-y,-z)
    float* Vjnkp = &f[17*elements +ic];//(+y,-z)
    float* Vjpkn = &f[18*elements +ic];//(-y,+z)

    if(nVec == -1)
    {
        float rho =( (*Vc)+(*Vjn)+(*Vjp)+(*Vkn)+(*Vkp)+(*Vjnkn)+(*Vjpkp)+(*Vjnkp)+(*Vjpkn) + 2.0f*((*Vip)+(*Vipjn)+(*Vipjp)+(*Vipkn)+(*Vipkp)))/( 1.0f-u);
        (*Vin) = (*Vip) +rho*u/ 3.0f;
        float Nyx = -rho*v/ 3.0f + 0.5f*((*Vjn)+(*Vjnkn)+(*Vjnkp)-((*Vjp)+(*Vjpkn)+(*Vjpkp)));
        float Nzx = -rho*w/ 3.0f + 0.5f*((*Vkn)+(*Vjnkn)+(*Vjpkn)-((*Vkp)+(*Vjnkp)+(*Vjpkp)));
        (*Vinjn) = (*Vipjp) +rho*(u+v)/ 6.0f -Nyx;
        (*Vinjp) = (*Vipjn) +rho*(u-v)/ 6.0f +Nyx;
        (*Vinkn) = (*Vipkp) +rho*(u+w)/ 6.0f -Nzx;
        (*Vinkp) = (*Vipkn) +rho*(u-w)/ 6.0f +Nzx;
    }
    else if(nVec == 1)
    {
        float rho =( (*Vc)+(*Vjn)+(*Vjp)+(*Vkn)+(*Vkp)+(*Vjnkn)+(*Vjpkp)+(*Vjnkp)+(*Vjpkn) + 2.0f*((*Vin)+(*Vinjn)+(*Vinjp)+(*Vinkn)+(*Vinkp)))/( 1.0f+u);
        (*Vip) = (*Vin) -rho*u/ 3.0f;
        float Nyx = -rho*v/ 3.0f + 0.5f*((*Vjn)+(*Vjnkn)+(*Vjnkp)-((*Vjp)+(*Vjpkn)+(*Vjpkp)));
        float Nzx = -rho*w/ 3.0f + 0.5f*((*Vkn)+(*Vjnkn)+(*Vjpkn)-((*Vkp)+(*Vjnkp)+(*Vjpkp)));
        (*Vipjp) = (*Vinjn) -rho*(u+v)/ 6.0f +Nyx;
        (*Vipjn) = (*Vinjp) -rho*(u-v)/ 6.0f -Nyx;
        (*Vipkp) = (*Vinkn) -rho*(u+w)/ 6.0f +Nzx;
        (*Vipkn) = (*Vinkp) -rho*(u-w)/ 6.0f -Nzx;
    }
    else if(nVec == -2)
    {
        float rho = ( (*Vc)+(*Vin)+(*Vip)+(*Vkn)+(*Vkp)+(*Vinkn)+(*Vipkp)+(*Vipkn)+(*Vinkp) + 2.0f*((*Vjp)+(*Vinjp)+(*Vipjp)+(*Vjpkn)+(*Vjpkp)))/( 1.0f-v);
        (*Vjn) = (*Vjp) +rho*v/ 3.0f;
        float Nxy = -rho*u/ 3.0f + 0.5f*((*Vin)+(*Vinkn)+(*Vinkp)-((*Vip)+(*Vipkn)+(*Vipkp)));
        float Nzy = -rho*w/ 3.0f + 0.5f*((*Vkn)+(*Vinkn)+(*Vipkn)-((*Vkp)+(*Vinkp)+(*Vipkp)));
        (*Vinjn) = (*Vipjp) +rho*(v+u)/ 6.0f -Nxy;
        (*Vipjn) = (*Vinjp) +rho*(v-u)/ 6.0f +Nxy;
        (*Vjnkn) = (*Vjpkp) +rho*(v+w)/ 6.0f -Nzy;
        (*Vjnkp) = (*Vjpkn) +rho*(v-w)/ 6.0f +Nzy;
    }
    else if(nVec == 2)
    {
        float rho =( (*Vc)+(*Vin)+(*Vip)+(*Vkn)+(*Vkp)+(*Vinkn)+(*Vipkp)+(*Vipkn)+(*Vinkp) + 2.0f*((*Vjn)+(*Vinjn)+(*Vipjn)+(*Vjnkn)+(*Vjnkp)))/( 1.0f+v);
        (*Vjp) = (*Vjn) -rho*v/ 3.0f;
        float Nxy = -rho*u/ 3.0f + 0.5f*((*Vin)+(*Vinkn)+(*Vinkp)-((*Vip)+(*Vipkn)+(*Vipkp)));
        float Nzy = -rho*w/ 3.0f + 0.5f*((*Vkn)+(*Vinkn)+(*Vipkn)-((*Vkp)+(*Vinkp)+(*Vipkp)));
        (*Vipjp) = (*Vinjn)  -rho*(v+u)/ 6.0f +Nxy;
        (*Vinjp) = (*Vipjn)  -rho*(v-u)/ 6.0f -Nxy;
        (*Vjpkp) = (*Vjnkn)  -rho*(v+w)/ 6.0f +Nzy;
        (*Vjpkn) = (*Vjnkp)  -rho*(v-w)/ 6.0f -Nzy;
    }
    else if(nVec == -3)
    {
        float rho = ( (*Vc)+(*Vin)+(*Vip)+(*Vjn)+(*Vjp)+(*Vinjn)+(*Vipjp)+(*Vinjp)+(*Vipjn)+ 2.0f*((*Vkp)+(*Vipkp)+(*Vinkp)+(*Vjpkp)+(*Vjnkp)))/( 1.0f-w);
        (*Vkn) = (*Vkp) +rho*w/ 3.0f;
        float Nxz = -rho*u/ 3.0f + 0.5f*((*Vin)+(*Vinjn)+(*Vinjp)-((*Vip)+(*Vipjn)+(*Vipjp)));
        float Nyz = -rho*v/ 3.0f + 0.5f*((*Vjn)+(*Vinjn)+(*Vipjp)-((*Vjp)+(*Vinjp)+(*Vipjp)));
        (*Vinkn) = (*Vipkp) +rho*(w+u)/ 6.0f -Nxz;
        (*Vipkn) = (*Vinkp) +rho*(w-u)/ 6.0f +Nxz;
        (*Vjnkn) = (*Vjpkp) +rho*(w+v)/ 6.0f -Nyz;
        (*Vjpkn) = (*Vjnkn) +rho*(w-v)/ 6.0f +Nyz;
    }
    else if(nVec == 3)
    {
        float rho = ( (*Vc)+(*Vin)+(*Vip)+(*Vjn)+(*Vjp)+(*Vinjn)+(*Vipjp)+(*Vinjp)+(*Vipjn)+ 2.0f*((*Vkn)+(*Vipkn)+(*Vinkn)+(*Vjpkn)+(*Vjnkn)))/( 1.0f+w);
        (*Vkp) = (*Vkn) -rho*w/ 3.0f;
        float Nxz = -rho*u/ 3.0f + 0.5f*((*Vin)+(*Vinjn)+(*Vinjp)-((*Vip)+(*Vipjn)+(*Vipjp)));
        float Nyz = -rho*v/ 3.0f + 0.5f*((*Vjn)+(*Vinjn)+(*Vipjp)-((*Vjp)+(*Vinjp)+(*Vipjp)));
        (*Vipkp) = (*Vinkn) -rho*(w+u)/ 6.0f +Nxz;
        (*Vinkp) = (*Vipkn) -rho*(w-u)/ 6.0f -Nxz;
        (*Vjpkp) = (*Vjnkn) -rho*(w+v)/ 6.0f +Nyz;
        (*Vjnkp) = (*Vjpkn) -rho*(w-v)/ 6.0f -Nyz;
    }
}


__kernel void k_streamingCollision // Pull
(
   __global float* f, __global float* fTmp,
   __global int* normal,
   __global float* sdf, __global unsigned char* solid, __global unsigned char* neiSolid, 
   __global float* u0, __global float* v0, __global float* w0,
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
        float cx[19] = {0.0f, 1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        float cy[19] = {0.0f, 0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f, 1.0f, -1.0f, -1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -1.0f, 1.0f, -1.0f};
        float cz[19] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -1.0f, -1.0f, 1.0f, 1.0f, -1.0f, -1.0f, 1.0f};

        float ft[19];
        int upID[19];
        for(int q = 0; q < 19; q++)
        {
            int qic = q*elements +ic;

            upID[q] = upwindID_B(q,i,j,k,nx,ny,nz,normal);
            if(upID[q] != -1)
            {   
                int upQID = idf(q, upID[q], nx, ny, nz);

                ft[q] = f[upQID];
                // // -- For Poiseuille flow (j == 0 and j == ny-1 are wall)
                // if(j == 0  && k == 0)
                // {
                //     int qbb = reflectQ(q);
                //     int bbQID = idf(qbb, ic, nx, ny, nz);
                //     if(q == 3 || q == 7 || q == 10 || q == 15 || q == 17)
                //     {
                //         // printf("i: %d, j: %d, k: %d, q: %d, upID: %d", i, j, k, q, upID);
                //         ft[q] = f[bbQID];
                //     }
                // }
                // if(j == ny-1  && k == 0)
                // {
                //     int qbb = reflectQ(q);
                //     int bbQID = idf(qbb, ic, nx, ny, nz);
                //     if(q == 4 || q == 8 || q == 9 || q == 16 || q == 18)
                //     {
                //         // printf("i: %d, j: %d, k: %d, q: %d, upID: %d", i, j, k, q, upID);
                //         ft[q] = f[bbQID];
                //     }
                // }
                // if(j == 0  && k == nz-1)
                // {
                //     int qbb = reflectQ(q);
                //     int bbQID = idf(qbb, ic, nx, ny, nz);
                //     if(q == 3 || q == 7 || q == 10 || q == 15 || q == 17)
                //     {
                //         // printf("i: %d, j: %d, k: %d, q: %d, upID: %d", i, j, k, q, upID);
                //         ft[q] = f[bbQID];
                //     }
                // }
                // if(j == ny-1  && k == nz-1)
                // {
                //     int qbb = reflectQ(q);
                //     int bbQID = idf(qbb, ic, nx, ny, nz);
                //     if(q == 4 || q == 8 || q == 9 || q == 16 || q == 18)
                //     {
                //         // printf("i: %d, j: %d, k: %d, q: %d, upID: %d", i, j, k, q, upID);
                //         ft[q] = f[bbQID];
                //     }
                // }
                // // --
            }
            else
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

        for(int q = 1; q < 19; q++)
        {
            if(neiSolid[ic] == 1)
            {
                if(solid[upID[q]] == 1)
                {
                    // const float sdf0 = sdf[ic];
                    // const float sdf1 = sdf[upID[q]];
                    // const float qf = fabs(sdf0)/(fabs(sdf0)+fabs(sdf1));
                    int qbb = reflectQ(q);
                    int bbQID = idf(qbb, ic, nx, ny, nz);
                    
                    // if(qf < 0.5f)
                    // {
                    //     ft[q] = (1.f -2.f*qf)*ft[qbb] +(qf*f[bbQID])*2.f;
                    //     // ft[q] = f[bbQID];
                    // }
                    // else
                    // {
                        // ft[q] = (1.f -0.5f/qf)*f[q] +(0.5f/qf)*f[bbQID];
                        ft[q] = f[bbQID];        
                    // }
                }
            }
        }
        
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
        }
    }
}
