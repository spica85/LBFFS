#ifndef LKS_H
#define LKS_H

#include <vector>

//D3Q19
inline int ic2i(int ic, int nx, int ny)
{
    return int((ic%(nx*ny))%nx);
}

inline int ic2j(int ic, int nx, int ny)
{
    return int(ic%(nx*ny)/nx);
}

inline int ic2k(int ic, int nx, int ny)
{
    return int(ic/(nx*ny));
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

class obstructure
{
    private:
    public:
    int boundary;
    std::vector<int> normal;
    bool inner;
    double u0;
    double v0;
    double w0;
    double rho0;
    obstructure()
    {
        boundary = 0;
        normal = {0, 0, 0};
        inner = false;

        u0 = 0.0;
        v0 = 0.0;
        w0 = 0.0;
        rho0 = 1.0;
    }
};

char* asciiToBinary(char* str, const float x);

inline double getAu(const double nu)
{
    return 1.0 -6.0*nu;
}

inline double getAp(const double Cs, const int nl)
{
    return 1.0 -3.0*Cs*Cs/double(nl);
}

inline const std::vector<double> setWt()
{
    return {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
}

inline const std::vector<double> setCx()
{
    return {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};
}

inline const std::vector<double> setCy()
{
    return {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0};
}

inline const std::vector<double> setCz()
{
    return {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0};
}

void boundaryConditionsP(obstructure& obst, std::vector<double>& p, int i, int j, int k, int nx, int ny, int nz);

void boundaryConditionsU(obstructure& obst, std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, int i, int j, int k, int nx, int ny, int nz);

void boundaryConditions(obstructure& obst, std::vector<double>& p, std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, int i, int j, int k, int nx, int ny, int nz);

double f_eq_in(const int q, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double p, const double u, const double v, const double w, const double Fx = 0.0, const double Fy = 0.0, const double Fz = 0.0);

int upwindID(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz);

double updateP(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Ap, const double Fx = 0.0, const double Fy = 0.0, const double Fz = 0.0);

double updateU(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const double Fx = 0.0, const double Fy = 0.0, const double Fz = 0.0);

double updateV(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const double Fx = 0.0, const double Fy = 0.0, const double Fz = 0.0);

double updateW(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const double Fx = 0.0, const double Fy = 0.0, const double Fz = 0.0);


double updateP(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Ap, const std::vector<int>& upID, const double Fx = 0.0, const double Fy = 0.0, const double Fz = 0.0);

double updateU(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const std::vector<int>& upID, const double Fx = 0.0, const double Fy = 0.0, const double Fz = 0.0);

double updateV(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const std::vector<int>& upID, const double Fx = 0.0, const double Fy = 0.0, const double Fz = 0.0);

double updateW(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const std::vector<int>& upID, const double Fx = 0.0, const double Fy = 0.0, const double Fz = 0.0);

void updateUVW(double& uNew, double& vNew, double& wNew, const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const std::vector<int>& upID, const double Fx = 0.0, const double Fy = 0.0, const double Fz = 0.0);

void updateUVW(const std::vector<double>& feq, double& uNew, double& vNew, double& wNew, const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au, const std::vector<int>& upID, const double Fx = 0.0, const double Fy = 0.0, const double Fz = 0.0);

double updatePfeq(std::vector<double>& feq, const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double>& p, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Ap, const std::vector<int>& upID, const double Fx = 0.0, const double Fy = 0.0, const double Fz = 0.0);


#endif