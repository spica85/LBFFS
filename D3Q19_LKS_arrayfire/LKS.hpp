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
    float u0;
    float v0;
    float w0;
    float rho0;
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

inline float getAu(const float nu)
{
    return 1.0 -6.0*nu;
}

inline float getAp(const float Cs, const int nl)
{
    return 1.0 -3.0*Cs*Cs/float(nl);
}

inline const std::vector<float> setWt()
{
    return {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
}

inline const std::vector<float> setCx()
{
    return {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};
}

inline const std::vector<float> setCy()
{
    return {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0};
}

inline const std::vector<float> setCz()
{
    return {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0};
}

void boundaryConditionsP(obstructure& obst, std::vector<float>& p, int i, int j, int k, int nx, int ny, int nz);

void boundaryConditionsU(obstructure& obst, std::vector<float>& u, std::vector<float>& v, std::vector<float>& w, int i, int j, int k, int nx, int ny, int nz);

void boundaryConditions(obstructure& obst, std::vector<float>& p, std::vector<float>& u, std::vector<float>& v, std::vector<float>& w, int i, int j, int k, int nx, int ny, int nz);

float f_eq_in(const int q, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, const float p, const float u, const float v, const float w, const float Fx = 0.0, const float Fy = 0.0, const float Fz = 0.0);

int upwindID(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz);

float updateP(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<float>& p, const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& w, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, const float Ap, const float Fx = 0.0, const float Fy = 0.0, const float Fz = 0.0);

float updateU(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<float>& p, const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& w, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, const float Au, const float Fx = 0.0, const float Fy = 0.0, const float Fz = 0.0);

float updateV(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<float>& p, const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& w, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, const float Au, const float Fx = 0.0, const float Fy = 0.0, const float Fz = 0.0);

float updateW(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<float>& p, const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& w, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, const float Au, const float Fx = 0.0, const float Fy = 0.0, const float Fz = 0.0);


float updateP(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<float>& p, const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& w, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, const float Ap, const std::vector<int>& upID, const float Fx = 0.0, const float Fy = 0.0, const float Fz = 0.0);

float updateU(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<float>& p, const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& w, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, const float Au, const std::vector<int>& upID, const float Fx = 0.0, const float Fy = 0.0, const float Fz = 0.0);

float updateV(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<float>& p, const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& w, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, const float Au, const std::vector<int>& upID, const float Fx = 0.0, const float Fy = 0.0, const float Fz = 0.0);

float updateW(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<float>& p, const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& w, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, const float Au, const std::vector<int>& upID, const float Fx = 0.0, const float Fy = 0.0, const float Fz = 0.0);

void updateUVW(float& uNew, float& vNew, float& wNew, const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<float>& p, const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& w, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, const float Au, const std::vector<int>& upID, const float Fx = 0.0, const float Fy = 0.0, const float Fz = 0.0);

void updateUVW(const std::vector<float>& feq, float& uNew, float& vNew, float& wNew, const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<float>& p, const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& w, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, const float Au, const std::vector<int>& upID, const float Fx = 0.0, const float Fy = 0.0, const float Fz = 0.0);

float updatePfeq(std::vector<float>& feq, const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<float>& p, const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& w, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, const float Ap, const std::vector<int>& upID, const float Fx = 0.0, const float Fy = 0.0, const float Fz = 0.0);


#endif