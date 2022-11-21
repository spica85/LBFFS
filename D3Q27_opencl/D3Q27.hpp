#ifndef D3Q27_H
#define D3Q27_H

#include <vector>
#include "base.hpp"

//D3Q27
inline const std::vector<float> setWt()
{
    return {8.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 2.0/27.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/54.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0};
}

inline const std::vector<float> setCx()
{
    return {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0};
}

inline const std::vector<float> setCy()
{
    return {0.0, 0.0,  0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0};
}

inline const std::vector<float> setCz()
{
    return {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0};
}

void cal_rhoUVW(int ic, int nx, int ny, int nz, const std::vector<float>& f, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, float& rho, float& u, float& v, float& w, float& dpdx)
{
    rho = 0.0;
    u = 0.0;
    v = 0.0;
    w = 0.0;

    for(int q = 0; q < 27; q++)
    {
        int qic = idf(q,ic,nx,ny,nz);
        rho += f[qic];
        u += f[qic]*cx[q];
        v += f[qic]*cy[q];
        w += f[qic]*cz[q];
    }
    u /= rho;
    v /= rho;
    w /= rho;                   
    u += 0.5f*dpdx/rho;
}

inline int upwindID(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz)
{
    const std::vector<float> cx = setCx();
    const std::vector<float> cy = setCy();
    const std::vector<float> cz = setCz();

    int upi = i -cx[q];
    int upj = j -cy[q];
    int upk = k -cz[q];

    upi = (upi == -1) ? nx-1 : ((upi == nx) ? 0 : upi);
    upj = (upj == -1) ? ny-1 : ((upj == ny) ? 0 : upj);
    upk = (upk == -1) ? nz-1 : ((upk == nz) ? 0 : upk);

    return index1d(upi,upj,upk,nx,ny);

    // if(q == 0)
    // {
    //     return index1d(i,j,k,nx,ny);
    // }
    // else if(q == 1)
    // {
    //     return i != 0 ? index1d(i-1,j,k,nx,ny) : index1d(nx-1,j,k,nx,ny);
    // }
    // else if(q == 2)
    // {
    //     return i != nx-1 ? index1d(i+1,j,k,nx,ny) : index1d(0,j,k,nx,ny);
    // }
    // else if(q == 3)
    // {
    //     return j != 0 ? index1d(i,j-1,k,nx,ny) : index1d(i,ny-1,k,nx,ny);
    // }
    // else if(q == 4)
    // {
    //     return j != ny-1 ? index1d(i,j+1,k,nx,ny) : index1d(i,0,k,nx,ny);
    // }
    // else if(q == 5)
    // {
    //     return k != 0 ? index1d(i,j,k-1,nx,ny) : index1d(i,j,nz-1,nx,ny);
    // }
    // else if(q == 6)
    // {
    //     return k != nz-1 ? index1d(i,j,k+1,nx,ny) : index1d(i,j,0,nx,ny);
    // }
    // else if(q == 7)
    // {
    //     return (i != 0 && j != 0) ? index1d(i-1, j-1, k, nx, ny) :
    //            (i == 0 && j != 0) ? index1d(nx-1, j-1, k, nx, ny) :
    //            (i != 0 && j == 0) ? index1d(i-1, ny-1, k, nx, ny) :
    //            index1d(nx-1, ny-1, k, nx, ny);
    // }
    // else if(q == 8)
    // {
    //     return (i != nx-1 && j != ny-1) ? index1d(i+1, j+1, k, nx, ny) :
    //            (i == nx-1 && j != ny-1) ? index1d(0, j+1, k, nx, ny) :
    //            (i != nx-1 && j == ny-1) ? index1d(i+1, 0, k, nx, ny) :
    //            index1d(0, 0, k, nx, ny);
    // }
    // else if(q == 9)
    // {
    //     return (i != 0 && j != ny-1) ? index1d(i-1, j+1, k, nx, ny) :
    //            (i == 0 && j != ny-1) ? index1d(nx-1, j+1, k, nx, ny) :
    //            (i != 0 && j == ny-1) ? index1d(i-1, 0, k, nx, ny) :
    //            index1d(nx-1, 0, k, nx, ny);
    // }
    // else if(q == 10)
    // {
    //     return (i != nx-1 && j != 0) ? index1d(i+1, j-1, k, nx, ny) :
    //            (i == nx-1 && j != 0) ? index1d(0, j-1, k, nx, ny) :
    //            (i != nx-1 && j == 0) ? index1d(i+1, ny-1, k, nx, ny) :
    //            index1d(0, ny-1, k, nx, ny);
    // }
    // else if(q == 11)
    // {
    //     return (i != 0 && k != 0) ? index1d(i-1, j, k-1, nx, ny) :
    //            (i == 0 && k != 0) ? index1d(nx-1, j, k-1, nx, ny) :
    //            (i != 0 && k == 0) ? index1d(i-1, j, nz-1, nx, ny) :
    //            index1d(nx-1, j, nz-1, nx, ny);
    // }
    // else if(q == 12)
    // {
    //     return (i != nx-1 && k != nz-1) ? index1d(i+1, j, k+1, nx, ny) :
    //            (i == nx-1 && k != nz-1) ? index1d(0, j, k+1, nx, ny) :
    //            (i != nx-1 && k == nz-1) ? index1d(i+1, j, 0, nx, ny) :
    //            index1d(0, j, 0, nx, ny);
    // }
    // else if(q == 13)
    // {
    //     return (i != 0 && k != nz-1) ? index1d(i-1, j, k+1, nx, ny) :
    //            (i == 0 && k != nz-1) ? index1d(nx-1, j, k+1, nx, ny) :
    //            (i != 0 && k == nz-1) ? index1d(i-1, j, 0, nx, ny) :
    //            index1d(nx-1, j, 0, nx, ny);
    // }
    // else if(q == 14)
    // {
    //     return (i != nx-1 && k != 0) ? index1d(i+1, j, k-1, nx, ny) :
    //            (i == nx-1 && k != 0) ? index1d(0, j, k-1, nx, ny) :
    //            (i != nx-1 && k == 0) ? index1d(i+1, j, nz-1, nx, ny) :
    //            index1d(0, j, nz-1, nx, ny);
    // }
    // else if(q == 15)
    // {
    //     return (j != 0 && k != 0) ? index1d(i, j-1, k-1, nx, ny) :
    //            (j == 0 && k != 0) ? index1d(i, ny-1, k-1, nx, ny) :
    //            (j != 0 && k == 0) ? index1d(i, j-1, nz-1, nx, ny) :
    //            index1d(i, ny-1, nz-1, nx, ny);
    // }
    // else if(q == 16)
    // {
    //     return (j != ny-1 && k != nz-1) ? index1d(i, j+1, k+1, nx, ny) :
    //            (j == ny-1 && k != nz-1) ? index1d(i, 0, k+1, nx, ny) :
    //            (j != ny-1 && k == nz-1) ? index1d(i, j+1, 0, nx, ny) :
    //            index1d(i, 0, 0, nx, ny);
    // }
    // else if(q == 17)
    // {
    //     return (j != 0 && k != nz-1) ? index1d(i, j-1, k+1, nx, ny) :
    //            (j == 0 && k != nz-1) ? index1d(i, ny-1, k+1, nx, ny) :
    //            (j != 0 && k == nz-1) ? index1d(i, j-1, 0, nx, ny) :
    //            index1d(i, ny-1, 0, nx, ny);
    // }
    // else if(q == 18)
    // {
    //     return (j != ny-1 && k != 0) ? index1d(i, j+1, k-1, nx, ny) :
    //            (j == ny-1 && k != 0) ? index1d(i, 0, k-1, nx, ny) :
    //            (j != ny-1 && k == 0) ? index1d(i, j+1, nz-1, nx, ny) :
    //            index1d(i, 0, nz-1, nx, ny);
    // }
    // else
    // {
    //     return 0;
    // }
}

#endif