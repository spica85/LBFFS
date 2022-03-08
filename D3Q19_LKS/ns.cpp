#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <omp.h>
#include <chrono>

//D3Q19
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
    int normal;
    bool inner;
    double u0;
    double v0;
    double w0;
    double rho0;
    obstructure()
    {
        boundary = 0;
        normal = 0;
        inner = false;

        u0 = 0.0;
        v0 = 0.0;
        w0 = 0.0;
        rho0 = 1.0;
    }
};

char* asciiToBinary(char* str, const float x)
{
    str[0] = ((char*)&x)[3 + 0];
    str[1] = ((char*)&x)[2 + 0];
    str[2] = ((char*)&x)[1 + 0];
    str[3] = ((char*)&x)[3 + 0];

    return str;    
}


void boundaryConditions(obstructure& obst, std::vector<double>& p, std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, int i, int j, int k, int nx, int ny, int nz)
{
    int ic = index1d(i,j,k,nx,ny);

    // Pressure
    if(obst.boundary == 1 || obst.boundary == 2) //wall
    {
        int innerID = 0;
        if(obst.normal == -1) //i=0
        {
            innerID = index1d(i+1,j,k,nx,ny);
        }
        else if(obst.normal == 1) //i=nx-1
        {
            innerID = index1d(i-1,j,k,nx,ny);
        }
        else if(obst.normal == -2) //j=0
        {
            innerID = index1d(i,j+1,k,nx,ny);
        }
        else if(obst.normal == 2) //j=ny-1
        {
            innerID = index1d(i,j-1,k,nx,ny);
        }
        else if(obst.normal == -3) //k=0
        {
            innerID = index1d(i,j,k+1,nx,ny);
        }
        else if(obst.normal == 3) //k=nz-1
        {
            innerID = index1d(i,j,k-1,nx,ny);
        }
        p[ic] = p[innerID];
    }

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

double f_eq_in(const int q, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double p, const double u, const double v, const double w)
{
    double uSqr =u*u+v*v+w*w;
    double uDotC = u*cx[q]+v*cy[q]+w*cz[q];

    // return wt[q]*3.0*p*(1.0 +3.0*uDotC +4.5*uDotC*uDotC -1.5*uSqr);
    return wt[q]*(3.0*p +3.0*uDotC +4.5*uDotC*uDotC -1.5*uSqr);
}

int upwindID(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz)
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

double updateP(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double> p, const std::vector<double> u, const std::vector<double> v, const std::vector<double> w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Ap)
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

        pNew += f_eq_in(q, cx, cy, cz, wt, pUpwind, uUpwind, vUpwind, wUpwind)
            + 3.0*Ap*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind));
    }
    pNew /= 3.0;

    return pNew;
}

double updateU(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double> p, const std::vector<double> u, const std::vector<double> v, const std::vector<double> w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au)
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

        uNew += cx[q]*(f_eq_in(q, cx, cy, cz, wt, pUpwind, uUpwind, vUpwind, wUpwind)
            + 3.0*Au*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind)));
    }

    return uNew;
}

double updateV(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double> p, const std::vector<double> u, const std::vector<double> v, const std::vector<double> w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au)
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

        vNew += cy[q]*(f_eq_in(q, cx, cy, cz, wt, pUpwind, uUpwind, vUpwind, wUpwind)
            + 3.0*Au*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind)));
    }

    return vNew;
}

double updateW(const int i, const int j, const int k,const int nx, const int ny, const int nz, const std::vector<double> p, const std::vector<double> u, const std::vector<double> v, const std::vector<double> w, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, const double Au)
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

        wNew += cz[q]*(f_eq_in(q, cx, cy, cz, wt, pUpwind, uUpwind, vUpwind, wUpwind)
            + 3.0*Au*wt[q]*(cx[q]*(u[ic]-uUpwind) +cy[q]*(v[ic]-vUpwind) +cz[q]*(w[ic]-wUpwind)));
    }

    return wNew;
}

void externalForce(const double dpdx, const int ic, const int nx, const int ny, const int nz, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, std::vector<double>& f)
{
    double rho = 0.0;
    for(int q = 0; q < 19; q++)
    {
        int qic = idf(q,ic,nx,ny,nz);
        rho += f[qic];
    }
    for(int q = 0; q < 19; q++)
    {
        int qic = idf(q,ic,nx,ny,nz);
        f[qic] += rho*wt[q]*3.0*dpdx*cx[q];
    }
}

template <typename Type = std::string>
Type returnWrapper(std::string arg) {
    return arg;
}

template <>
int returnWrapper<int>(std::string arg) {
  return std::stoi(arg);
}

template <>
bool returnWrapper<bool>(std::string arg) {
  return arg == "true" ? true : false;
}

template <>
float returnWrapper<float>(std::string arg) {
  return std::stof(arg);
}

template <>
double returnWrapper<double>(std::string arg) {
  return std::stod(arg);
}


template<typename Type> Type lookup(std::vector<std::string>& lines, std::string& str)
{
    for(auto itr = std::begin(lines); itr != std::end(lines); ++itr)
    {
        if(*itr == str)
        {
            ++itr;
            Type i = returnWrapper<Type>(*itr);
            --itr;
            ++itr;
            std::cout <<
                str << ": " << i << std::endl;
            --itr;
            return i;
        }
    }
    std::cerr << "Could not find " << str << std::endl;

    exit(EXIT_FAILURE);
}




int main()
{
    std::string inputFileName("input.txt");
    std::vector<std::string> lines;
    std::string line;
    std::ifstream inputFile(inputFileName);

    if(!inputFile.is_open())
    {
        std::cerr << "Could not open the file - '"
            << inputFileName << "'" << std::endl;
        exit(EXIT_FAILURE);
    }


    while(std::getline(inputFile, line))
    {
        std::vector<std::string> list_string;
        boost::split(list_string,line, boost::is_space());
        for(auto& str: list_string)
        {
            lines.push_back(str);
        }
    }

    std::string nThreadsStr("nThreads");
    omp_set_num_threads(lookup<int>(lines, nThreadsStr));

    std::string restartStr("restart");
    bool restart = lookup<bool>(lines, restartStr);

    std::string FwriteStr("Fwrite");
    bool Fwrite = lookup<bool>(lines, FwriteStr);

    std::string writeBinaryStr("writeBinary");
    bool writeBinary = lookup<bool>(lines, writeBinaryStr);
    
    std::string startTimeStepStr("startTimeStep");
    int startTimeStep = lookup<int>(lines, startTimeStepStr);
    
    std::string endTimeStepStr("endTimeStep");
    const int endTimeStep = lookup<int>(lines, endTimeStepStr);

    int nextOutTime = startTimeStep;

    std::string outIntervalStr("outInterval");
    const int outInterval = lookup<int>(lines, outIntervalStr);
    std::cout << std::endl;

    std::string nxStr("nx");
    const int nx = lookup<int>(lines, nxStr);

    std::string nyStr("ny");
    const int ny = lookup<int>(lines, nyStr);

    std::string nzStr("nz");
    const int nz = lookup<int>(lines, nzStr);
    std::cout << std::endl;

    std::string u0Str("u0");
    const double u0 = lookup<double>(lines, u0Str);

    std::string rho0Str("rho0");
    const double rho0 = lookup<double>(lines, rho0Str);

    std::string ReStr("Re");
    const double Re = lookup<double>(lines, ReStr);
    std::cout << std::endl;

    inputFile.close();

    // Improved Lattice Kinetic Scheme model

    //For cavity flow
    double nu = std::abs(u0)*double(nx)/Re;
    double dpdx = 0.0;

    double Au = 1.0 -6.0*nu;
    int nl = 5;
    double Ma = 0.1;
    double Cs = u0/Ma;
    double Ap = 1.0 -3.0*Cs*Cs/double(nl);
        

    //For channel flow
    // double Retau = 10;
    // double Retau = 1;
    // double utau = 0.005;
    // double nu = utau*0.5*ny/Retau;
    // double dpdx = utau*utau/(0.5*ny);

    //For Poiseuille flow
    // double umax = 0.1;
    // double h = 1.0;
    // double nu = umax*h/Re;
    // double dpdx = umax/(h*h)*8.0*nu/(ny-1);
    // double dpdx = 0.00001;

    // std::cout << "dpdx = " << dpdx << std::endl;

    // std::cout << "nu = " << nu << std::endl;
    
    // nu = nu*(ny-1);

    // double omega = 1.0/(3.0*nu +0.5);
    // double omega = 1.0/0.56;

    // std::cout << "tau = " << 1.0/omega << std::endl;

    const std::vector<double> wt = {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

    // D3Q19 model
    const std::vector<double> cx = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};
    const std::vector<double> cy = {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0};
    const std::vector<double> cz = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0};
    
    std::vector<double> p(nx*ny*nz);
    std::vector<double> u(nx*ny*nz);
    std::vector<double> v(nx*ny*nz);
    std::vector<double> w(nx*ny*nz);

    std::vector<obstructure> obst(nx*ny*nz);
    

    // Setting conditions
    #include "boundaryCondition.hpp"
    // #include "boundaryCondition_channelFlow.hpp"
    // if(restart)
    // {
    //     #include "restart.hpp"
    // }
    // else
    // {
        #include "initialization.hpp"        
    // }

    std::chrono::system_clock::time_point start;
    std::chrono::system_clock::time_point end;
    start = std::chrono::system_clock::now();

    // Time marching
    for(int nt = startTimeStep; nt <= endTimeStep; nt++)
    {   
        #include "write.hpp"

        std::vector<double> uTmp = u;
        std::vector<double> vTmp = v;
        std::vector<double> wTmp = w;
        std::vector<double> pTmp = p;

        #pragma omp parallel for
        for(int k = 0; k < nz; k++)
        {            
            for(int j = 0; j < ny; j++)                    
            {
                for(int i = 0; i < nx; i++)                        
                {
                    int ic = index1d(i,j,k,nx,ny);
                    for(int l = 0; l < nl; l++)
                    {
                        p[ic] = updateP(i,j,k,nx,ny,nz,pTmp,uTmp,vTmp,wTmp,cx,cy,cz,wt,Ap);
                    }
                }
            }
        }


        #pragma omp parallel for
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 0; i < nx; i++)
                {
                    int ic = index1d(i,j,k,nx,ny);

                    u[ic] = updateU(i,j,k,nx,ny,nz,p,uTmp,vTmp,wTmp,cx,cy,cz,wt,Au);
                    v[ic] = updateV(i,j,k,nx,ny,nz,p,uTmp,vTmp,wTmp,cx,cy,cz,wt,Au);
                    w[ic] = updateW(i,j,k,nx,ny,nz,p,uTmp,vTmp,wTmp,cx,cy,cz,wt,Au);
                }
            }
        }

        #pragma omp parallel for
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 0; i < nx; i++)
                {
                    int ic = index1d(i,j,k,nx,ny);
                    boundaryConditions(obst[ic], p, u, v, w, i, j, k, nx, ny, nz);
                }
            }
        }
    }
    end = std::chrono::system_clock::now();
    double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() *1e-6);
    std::cout << "Execution time: " << time << " (s)" << std::endl;
    std::cout << "Speed: " << double(endTimeStep)*double(nx*ny*nz)/time << " (LUPS)" << std::endl;

    return EXIT_SUCCESS;
}
