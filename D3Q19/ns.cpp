#include <boost/algorithm/string.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <omp.h>
#include <chrono>

#include "D3Q19.hpp"
#include "input.hpp"

int main()
{
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
    float uMax;
    float rho0;
    float Re;

    input(restart, Fwrite, writeBinary, startTimeStep, endTimeStep, nextOutTime, outInterval, nx, ny, nz, uMax, rho0, Re);

    // Single Relaxation Time model
    // // -- For cavity flow --
    const float a = 1.0; //(m)
    const float L = a/float(nx-1);

    const float u0 = 0.1;
    const float c = uMax/u0;

    const float deltaT = L/c;

    float nu = uMax*a/Re;
    nu = nu/(L*c);
    float dpdx = 0.0;


    //For channel flow
    // float Retau = 10;
    // float Retau = 1;
    // float utau = 0.005;
    // float nu = utau*0.5*ny/Retau;
    // float dpdx = utau*utau/(0.5*ny);

    //For Poiseuille flow
    // float umax = 0.1;
    // float h = 1.0;
    // float nu = umax*h/Re;
    // float dpdx = umax/(h*h)*8.0*nu/(ny-1);
    // float dpdx = 0.00001;

    // std::cout << "dpdx = " << dpdx << std::endl;

    // std::cout << "nu = " << nu << std::endl;
    
    // nu = nu*(ny-1);

    float omega = 1.0/(3.0*nu +0.5);
    // float omega = 1.0/0.56;

    std::cout << "tau = " << 1.0/omega << std::endl;

    const std::vector<float> wt = setWt();

    // D3Q19 model
    const std::vector<float> cx = setCx();
    const std::vector<float> cy = setCy();
    const std::vector<float> cz = setCz();

    std::vector<float> f(19*nx*ny*nz);
    std::vector<float> ftmp(19*nx*ny*nz);
    
    std::vector<float> rho(nx*ny*nz);
    std::vector<float> u(nx*ny*nz);
    std::vector<float> v(nx*ny*nz);
    std::vector<float> w(nx*ny*nz);

    std::vector<obstructure> obst(nx*ny*nz);
    

    // Setting conditions
    #include "boundaryCondition_cavityFlow2d.hpp"
    // #include "boundaryCondition_cavityFlow3dDiagonal.hpp"
    // #include "boundaryCondition_channelFlow.hpp"
    if(restart)
    {
        #include "restart.hpp"
    }
    else
    {
        #include "initialization.hpp"        
    }

    std::vector<int> upID(19*nx*ny*nz);
    for(int ic = 0; ic < nx*ny*nz; ic++)
    {
        int i = ic2i(ic,nx,ny);
        int j = ic2j(ic,nx,ny);
        int k = ic2k(ic,nx,ny);
        for(int q = 0; q < 19; q++)
        {
            int qic = idf(q,ic,nx,ny,nz);
            upID[qic] = upwindID(q,i,j,k,nx,ny,nz);
        }
    }

    std::vector<int> internalID;
    std::vector<int> boundaryID;
    for(int ic = 0; ic < nx*ny*nz; ic++)
    {
        if(obst[ic].boundary == 0)
        {
            internalID.push_back(ic);
        }
        else
        {
            boundaryID.push_back(ic);
        }
    }

    std::chrono::system_clock::time_point start;
    std::chrono::system_clock::time_point end;
    start = std::chrono::system_clock::now();

    // Time marching
    for(int nt = startTimeStep; nt <= endTimeStep; nt++)
    {
        #pragma omp parallel for
        for(int ic =0; ic <nx*ny*nz; ic++)
        {
            cal_rhoUVW(ic, nx, ny, nz, f, cx, cy, cz, rho[ic], u[ic], v[ic], w[ic]);
        }
        #include "write.hpp"

        #pragma omp parallel for
        for(int ic = 0; ic < nx*ny*nz; ic++)
        {
            collision(omega,ic,nx,ny,nz,cx,cy,cz,wt,f);
            externalForce(dpdx,ic,nx,ny,nz,cx,cy,cz,wt,f);
            for(int q = 0; q < 19; q++)
            {
                int qic = idf(q,ic,nx,ny,nz);
                ftmp[qic] = f[qic];
            }
        }

        #pragma omp parallel for
        for(int ic = 0; ic < nx*ny*nz; ic++)
        {
            int i = ic2i(ic,nx,ny);
            int j = ic2j(ic,nx,ny);
            int k = ic2k(ic,nx,ny);

            streaming(ic,i,j,k,nx,ny,nz,ftmp,f,upID);
        }

        #pragma omp parallel for
        for(int bID = 0; bID < boundaryID.size(); bID++)
        {
            int ic = boundaryID[bID];
            boundaryConditions(obst[ic],f,ic,nx,ny,nz);
        }
    }
    end = std::chrono::system_clock::now();
    double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() *1e-6);
    std::cout << "Execution time: " << time << " (s)" << std::endl;
    std::cout << "Speed: " << double(endTimeStep)*double(nx*ny*nz)/time*1e-6 << " (MLUPS)" << std::endl;

    return EXIT_SUCCESS;
}
