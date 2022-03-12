#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <omp.h>
#include <chrono>
#include "LKS.hpp"
#include "input.hpp"

int main()
{
    // Improved Lattice Kinetic Scheme model

    bool restart;
    bool Fwrite;
    bool writeBinary;
    int startTimeStep;
    int endTimeStep;
    int nextOutTime;
    int outInterval;
    int nx;
    int ny;
    int nz;

    input(restart, Fwrite, writeBinary, startTimeStep, endTimeStep, nextOutTime, outInterval, nx, ny, nz);

    //For cavity flow
    const double u0 = 0.1;
    const double rho0 = 1.0;
    const double Re = 100.0;
    double nu = std::abs(u0)*double(nx-1)/Re;
    double dpdx = 0.0;

    double Au = 1.0 -6.0*nu;
    int nl = 5;
    double Ma = 0.08;
    double Cs = u0/Ma;
    double Ap = 1.0 -3.0*Cs*Cs/double(nl);
    Ap = 0.0;

    std::cout << "Au: " << Au << ", Ap: " << Ap << std::endl;
        

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


        
        for(int l = 0; l < nl; l++)
        {
            #pragma omp parallel for                
            for(int ic = 0; ic < nx*ny*nz; ic++)
            {
                int i = ic2i(ic, nx, ny);
                int j = ic2j(ic, nx, ny);
                int k = ic2k(ic, nx, ny);

                p[ic] = updateP(i,j,k,nx,ny,nz,pTmp,uTmp,vTmp,wTmp,cx,cy,cz,wt,Ap);
            }
            pTmp = p;
            
            #pragma omp parallel for                
            for(int ic = 0; ic < nx*ny*nz; ic++)
            {
                int i = ic2i(ic, nx, ny);
                int j = ic2j(ic, nx, ny);
                int k = ic2k(ic, nx, ny);

                boundaryConditionsP(obst[ic], p, i, j, k, nx, ny, nz);
            }
        }

        #pragma omp parallel for
        for(int ic = 0; ic < nx*ny*nz; ic++)
        {

            int i = ic2i(ic, nx, ny);
            int j = ic2j(ic, nx, ny);
            int k = ic2k(ic, nx, ny);
            
            u[ic] = updateU(i,j,k,nx,ny,nz,p,uTmp,vTmp,wTmp,cx,cy,cz,wt,Au);
            v[ic] = updateV(i,j,k,nx,ny,nz,p,uTmp,vTmp,wTmp,cx,cy,cz,wt,Au);
            w[ic] = updateW(i,j,k,nx,ny,nz,p,uTmp,vTmp,wTmp,cx,cy,cz,wt,Au);
        }

        #pragma omp parallel for
        for(int ic = 0; ic < nx*ny*nz; ic++)                
        {
            int i = ic2i(ic, nx, ny);
            int j = ic2j(ic, nx, ny);
            int k = ic2k(ic, nx, ny);

            boundaryConditionsU(obst[ic], u, v, w, i, j, k, nx, ny, nz);
        }
    }
    end = std::chrono::system_clock::now();
    double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() *1e-6);
    std::cout << "Execution time: " << time << " (s)" << std::endl;
    std::cout << "Speed: " << double(endTimeStep)*double(nx*ny*nz)/time << " (LUPS)" << std::endl;

    return EXIT_SUCCESS;
}
