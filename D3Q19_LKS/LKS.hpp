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

void boundaryConditionsP(obstructure& obst, std::vector<double>& p, int i, int j, int k, int nx, int ny, int nz);

void boundaryConditionsU(obstructure& obst, std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, int i, int j, int k, int nx, int ny, int nz);

void boundaryConditions(obstructure& obst, std::vector<double>& p, std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, int i, int j, int k, int nx, int ny, int nz);

double f_eq_in(const int q, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double p, const double u, const double v, const double w);

int upwindID(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz);

double updateP(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double> p, const std::vector<double> u, const std::vector<double> v, const std::vector<double> w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Ap);

double updateU(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double> p, const std::vector<double> u, const std::vector<double> v, const std::vector<double> w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au);

double updateV(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double> p, const std::vector<double> u, const std::vector<double> v, const std::vector<double> w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au);

double updateW(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double> p, const std::vector<double> u, const std::vector<double> v, const std::vector<double> w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au);

void externalForce(const double dpdx, const int ic, const int nx, const int ny, const int nz, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, std::vector<double>& f);


#endif