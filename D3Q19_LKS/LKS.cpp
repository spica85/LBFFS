//D3Q19
#include "LKS.hpp"
#include <iostream>
#include <string>
#include <vector>


char* asciiToBinary(char* str, const float x)
{
    str[0] = ((char*)&x)[3 + 0];
    str[1] = ((char*)&x)[2 + 0];
    str[2] = ((char*)&x)[1 + 0];
    str[3] = ((char*)&x)[3 + 0];

    return str;    
}


void boundaryConditionsP(obstructure& obst, std::vector<double>& p, int i, int j, int k, int nx, int ny, int nz)
{
    int ic = index1d(i,j,k,nx,ny);

    // Pressure
    if(obst.boundary == 1 || obst.boundary == 2) //wall
    {
        int innerID = index1d(i-obst.normal[0],j-obst.normal[1],k-obst.normal[2],nx,ny);
        p[ic] = p[innerID];
    }
}

void boundaryConditionsU(obstructure& obst, std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, int i, int j, int k, int nx, int ny, int nz)
{
    int ic = index1d(i,j,k,nx,ny);

    // Velocity
    if(obst.boundary == 1) //fixed wall
    {
        u[ic] = 0.0;
        v[ic] = 0.0;
        w[ic] = 0.0;
    }
    else if(obst.boundary == 2) //wall velocity
    {
        u[ic] = obst.u0;
        v[ic] = obst.v0;
        w[ic] = obst.w0;
    }
}

void boundaryConditions(obstructure& obst, std::vector<double>& p, std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, int i, int j, int k, int nx, int ny, int nz)
{
    boundaryConditionsP(obst, p, i, j, k, nx, ny, nz);
    boundaryConditionsU(obst, u, v, w, i, j, k, nx, ny, nz);
}

double f_eq_in(const double cx, const double cy, const double cz, const double wt, const double p, const double u, const double v, const double w, const double Fx, const double Fy, const double Fz)
{
    double uSqr =u*u+v*v+w*w;
    double uDotC = u*cx+v*cy+w*cz;
    double FDotC = Fx*cx+Fy*cy+Fz*cz;

    return wt*(3.0*p +3.0*uDotC +4.5*uDotC*uDotC -1.5*uSqr +3.0*FDotC);
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

double updateP(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Ap, const double Fx, const double Fy, const double Fz)
{
    int ic = index1d(i,j,k,nx,ny);
    double pNew = 0.0;
    for(int q = 0; q < 19; q++)
    {
        int iUpwind = upwindID(q,i,j,k,nx,ny,nz);

        const double pUpwind = p[iUpwind];
        const double uUpwind = u[iUpwind];
        const double vUpwind = v[iUpwind];
        const double wUpwind = w[iUpwind];

        pNew += f_eq_in(cx[q], cy[q], cz[q], wt[q], pUpwind, uUpwind, vUpwind, wUpwind, Fx, Fy, Fz)
            + 3.0*Ap*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind));
    }
    pNew /= 3.0;

    return pNew;
}

double updateU(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const double Fx, const double Fy, const double Fz)
{
    int ic = index1d(i,j,k,nx,ny);
    double uNew = 0.0;
    std::vector<int> qx = {1, 2, 7, 8, 9, 10, 11, 12, 13, 14};
    for(int qxi = 0; qxi < 10; qxi++)
    {
        int q = qx[qxi];
        int iUpwind = upwindID(q,i,j,k,nx,ny,nz);

        const double pUpwind = p[iUpwind];
        const double uUpwind = u[iUpwind];
        const double vUpwind = v[iUpwind];
        const double wUpwind = w[iUpwind];

        uNew += cx[q]*(f_eq_in(cx[q], cy[q], cz[q], wt[q], pUpwind, uUpwind, vUpwind, wUpwind, Fx, Fy, Fz)
            + 3.0*Au*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind)));
    }

    return uNew;
}

double updateV(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const double Fx, const double Fy, const double Fz)
{
    int ic = index1d(i,j,k,nx,ny);
    double vNew = 0.0;
    std::vector<int> qy = {3, 4, 7, 8, 9, 10, 15, 16, 17, 18};
    for(int qyi = 0; qyi < 10; qyi++)
    {
        int q = qy[qyi];
        int iUpwind = upwindID(q,i,j,k,nx,ny,nz);

        const double pUpwind = p[iUpwind];
        const double uUpwind = u[iUpwind];
        const double vUpwind = v[iUpwind];
        const double wUpwind = w[iUpwind];

        vNew += cy[q]*(f_eq_in(cx[q], cy[q], cz[q], wt[q], pUpwind, uUpwind, vUpwind, wUpwind, Fx, Fy, Fz)
            + 3.0*Au*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind)));
    }

    return vNew;
}

double updateW(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const double Fx, const double Fy, const double Fz)
{
    int ic = index1d(i,j,k,nx,ny);
    double wNew = 0.0;
    std::vector<int> qz = {5, 6, 11, 12, 13, 14, 15, 16, 17, 18};
    for(int qzi = 0; qzi < 10; qzi++)
    {
        int q = qz[qzi];
        int iUpwind = upwindID(q,i,j,k,nx,ny,nz);

        const double pUpwind = p[iUpwind];
        const double uUpwind = u[iUpwind];
        const double vUpwind = v[iUpwind];
        const double wUpwind = w[iUpwind];

        wNew += cz[q]*(f_eq_in(cx[q], cy[q], cz[q], wt[q], pUpwind, uUpwind, vUpwind, wUpwind, Fx, Fy, Fz)
            + 3.0*Au*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind)));
    }

    return wNew;
}


double updateP(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Ap, const std::vector<int>& upID, const double Fx, const double Fy, const double Fz)
{
    double pNew = 0.0;
    for(int q = 0; q < 19; q++)
    {
        int qic = idf(q,ic,nx,ny,nz);
        int iUpwind = upID[qic];

        const double pUpwind = p[iUpwind];
        const double uUpwind = u[iUpwind];
        const double vUpwind = v[iUpwind];
        const double wUpwind = w[iUpwind];

        pNew += f_eq_in(cx[q], cy[q], cz[q], wt[q], pUpwind, uUpwind, vUpwind, wUpwind, Fx, Fy, Fz)
            + 3.0*Ap*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind));
    }
    pNew /= 3.0;

    return pNew;
}

double updateU(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const std::vector<int>& upID, const double Fx, const double Fy, const double Fz)
{    
    double uNew = 0.0;
    std::vector<int> qx = {1, 2, 7, 8, 9, 10, 11, 12, 13, 14};
    for(int qxi = 0; qxi < 10; qxi++)
    {
        int q = qx[qxi];
        int qic = idf(q,ic,nx,ny,nz);
        int iUpwind = upID[qic];

        const double pUpwind = p[iUpwind];
        const double uUpwind = u[iUpwind];
        const double vUpwind = v[iUpwind];
        const double wUpwind = w[iUpwind];

        uNew += cx[q]*(f_eq_in(cx[q], cy[q], cz[q], wt[q], pUpwind, uUpwind, vUpwind, wUpwind, Fx, Fy, Fz)
            + 3.0*Au*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind)));
    }

    return uNew;
}

double updateV(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const std::vector<int>& upID, const double Fx, const double Fy, const double Fz)
{
    double vNew = 0.0;
    std::vector<int> qy = {3, 4, 7, 8, 9, 10, 15, 16, 17, 18};
    for(int qyi = 0; qyi < 10; qyi++)
    {
        int q = qy[qyi];
        int qic = idf(q,ic,nx,ny,nz);
        int iUpwind = upID[qic];

        const double pUpwind = p[iUpwind];
        const double uUpwind = u[iUpwind];
        const double vUpwind = v[iUpwind];
        const double wUpwind = w[iUpwind];

        vNew += cy[q]*(f_eq_in(cx[q], cy[q], cz[q], wt[q], pUpwind, uUpwind, vUpwind, wUpwind, Fx, Fy, Fz)
            + 3.0*Au*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind)));
    }

    return vNew;
}

double updateW(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const std::vector<int>& upID, const double Fx, const double Fy, const double Fz)
{
    double wNew = 0.0;
    std::vector<int> qz = {5, 6, 11, 12, 13, 14, 15, 16, 17, 18};
    for(int qzi = 0; qzi < 10; qzi++)
    {
        int q = qz[qzi];
        int qic = idf(q,ic,nx,ny,nz);
        int iUpwind = upID[qic];

        const double pUpwind = p[iUpwind];
        const double uUpwind = u[iUpwind];
        const double vUpwind = v[iUpwind];
        const double wUpwind = w[iUpwind];

        wNew += cz[q]*(f_eq_in(cx[q], cy[q], cz[q], wt[q], pUpwind, uUpwind, vUpwind, wUpwind, Fx, Fy, Fz)
            + 3.0*Au*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind)));
    }

    return wNew;
}

void updateUVW(double& uNew, double& vNew, double& wNew, const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const std::vector<int>& upID, const double Fx, const double Fy, const double Fz)
{
    uNew = 0.0;
    vNew = 0.0;
    wNew = 0.0;

    for(int q = 0; q < 19; q++)
    {
        int qic = idf(q,ic,nx,ny,nz);
        int iUpwind = upID[qic];

        const double pUpwind = p[iUpwind];
        const double uUpwind = u[iUpwind];
        const double vUpwind = v[iUpwind];
        const double wUpwind = w[iUpwind];

        const double feq = f_eq_in(cx[q], cy[q], cz[q], wt[q], pUpwind, uUpwind, vUpwind, wUpwind, Fx, Fy, Fz);
        const double cDotDeltaU = cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind);
        const double f = feq + 3.0*Au*wt[q]*cDotDeltaU;

        uNew += cx[q]*f;
        vNew += cy[q]*f;
        wNew += cz[q]*f;
    }
}

void updateUVW(const std::vector<double>& feq, double& uNew, double& vNew, double& wNew, const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const std::vector<int>& upID, const double Fx, const double Fy, const double Fz)
{
    uNew = 0.0;
    vNew = 0.0;
    wNew = 0.0;

    for(int q = 0; q < 19; q++)
    {
        int qic = idf(q,ic,nx,ny,nz);
        int iUpwind = upID[qic];

        const double pUpwind = p[iUpwind];
        const double uUpwind = u[iUpwind];
        const double vUpwind = v[iUpwind];
        const double wUpwind = w[iUpwind];

        const double cDotDeltaU = cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind);
        const double f = feq[qic] + 3.0*Au*wt[q]*cDotDeltaU;

        uNew += cx[q]*f;
        vNew += cy[q]*f;
        wNew += cz[q]*f;
    }
}

double updatePfeq(std::vector<double>& feq, const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Ap, const std::vector<int>& upID, const double Fx, const double Fy, const double Fz)
{
    double pNew = 0.0;
    for(int q = 0; q < 19; q++)
    {
        int qic = idf(q,ic,nx,ny,nz);
        int iUpwind = upID[qic];

        const double pUpwind = p[iUpwind];
        const double uUpwind = u[iUpwind];
        const double vUpwind = v[iUpwind];
        const double wUpwind = w[iUpwind];

        feq[qic] = f_eq_in(cx[q], cy[q], cz[q], wt[q], pUpwind, uUpwind, vUpwind, wUpwind, Fx, Fy, Fz);

        pNew += feq[qic]
            + 3.0*Ap*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind));
    }
    pNew /= 3.0;

    return pNew;
}