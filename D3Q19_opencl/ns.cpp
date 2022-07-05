#include <boost/algorithm/string.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <chrono>

#include <omp.h>

#include "D3Q19.hpp"
#include "input.hpp"

//-- For OpenCL
#ifdef __APPLE__
#define CL_SILENCE_DEPRECATION
#endif
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#include <CL/cl2.hpp>

#include "util.hpp"

#include "err_code.h"

//-- pick up device type from compiler command line or from the default type
#ifndef DEVICE
#define DEVICE CL_DEVICE_TYPE_DEFAULT
#endif
//--
//-- 

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
    float U0;

    input(restart, Fwrite, writeBinary, startTimeStep, endTimeStep, nextOutTime, outInterval, nx, ny, nz, uMax, rho0, Re, U0);

    //- For cavity flow
    // const float a = 1.0; //Dimensional length of system (m)
    // const float L = a/float(nx); //Representative length (-)  
    // float nu = uMax*a/Re; //Dimensional kinematic viscosity
    // float dpdx = 0.0;//Dimensionless external force
    //--

    //For channel flow
    // float Retau = 10;
    // float Retau = 1;
    // float utau = 0.005;
    // float nu = utau*0.5*ny/Retau;
    // float dpdx = utau*utau/(0.5*ny);

    //- For Poiseuille flow
    // float h = 1.f;
    // float nu = uMax*h/Re;
    // const float L = h/float(ny);
    // float dpdx = 8.0*nu*uMax/(h*h);
    //--
    

    //- For flow around cylinder
    // float h = 1.f;
    // // float d = 6.f;
    // float d = 0.1f;
    // float nu = uMax*d/Re;
    // const float L = h/float(ny);
    // float dpdx = 0.f;
    //--

    //- For backward facing step flow
    float h = 1.f;
    float Ly = 3.f*h; // 5.f*h// Flow around Digital Science
    float Lx = 40.f*h;
    float nu = uMax*h/Re;
    const float L = Ly/float(ny);
    float dpdx = 0.f;
    //--

    const float c = uMax/U0; //Representative velocity (m/s)
    nu = nu/(L*c);
    dpdx *= L/(c*c);

    std::cout << "Size: "
              << "(" << L/2 << ", " << L/2 << ", " << L/2 << "), "
              <<  "(" << L/2+L*(nx-1) << ", " << L/2+L*(ny-1) << ", " << L/2+L*(nz-1) << ")"
              << std::endl;

    std::cout << "dpdx = " << dpdx << std::endl;
    std::cout << "nu = " << nu << std::endl;


    const float deltaT = L/c; //Dimensional time step (s)
    const float omega = 1.0/(3.0*nu +0.5);

    std::cout << "tau = " << 1.0/omega << ", taulim = " << 0.5f +uMax/c/8.0f << std::endl;
    std::cout << "Maximum Ma = " << (uMax/c)*sqrt(3.f) << std::endl;

    const std::vector<float> wt = setWt();

    //-- D3Q19 model
    const std::vector<float> cx = setCx();
    const std::vector<float> cy = setCy();
    const std::vector<float> cz = setCz();

    std::vector<float> f(19*nx*ny*nz);
    std::vector<float> fTmp(19*nx*ny*nz);
    //--
    
    std::vector<float> rho(nx*ny*nz);
    std::vector<float> u(nx*ny*nz);
    std::vector<float> v(nx*ny*nz);
    std::vector<float> w(nx*ny*nz);
    float rho_av = 1.0;

    std::vector<float> tauSGS(nx*ny*nz);

    std::vector<obstructure> obst(nx*ny*nz);
    std::vector<float> u0(nx*ny*nz,0.0f);
    std::vector<float> v0(nx*ny*nz,0.0f);
    std::vector<float> w0(nx*ny*nz,0.0f);
    std::vector<int> boundary1(nx*ny*nz,0);
    std::vector<int> boundary2(nx*ny*nz,0);
    std::vector<int> boundary3(nx*ny*nz,0);

    const unsigned elements = nx*ny*nz;
    const unsigned qElements = 19*elements;

    // Setting conditions
    // #include "boundaryCondition_cavityFlow.hpp"
    // #include "boundaryCondition_PoiseuilleFlow.hpp"
    // #include "boundaryCondition_flowAroundCylinder.hpp"
    #include "boundaryCondition_backStepFlow.hpp"
    if(restart)
    {
        #include "restart.hpp"
    }
    else
    {
        #include "initialization.hpp"        
    }

    for(int ic = 0; ic < nx*ny*nz; ic++)
    {
        u0[ic] = obst[ic].u0;
        v0[ic] = obst[ic].v0;
        w0[ic] = obst[ic].w0;

        boundary1[ic] = obst[ic].boundary1;
        boundary2[ic] = obst[ic].boundary2;
        boundary3[ic] = obst[ic].boundary3;
    }

    //-- Reading and settings of STL 
    const std::string STLname("walls.stl");
    std::vector<std::vector<float> > STLnormal(3);
    std::vector<std::vector<float> > STLv0(3);
    std::vector<std::vector<float> > STLv1(3);
    std::vector<std::vector<float> > STLv2(3);

    std::vector<float> Fwx(nx*ny*nz,0.0f);
    std::vector<float> Fwy(nx*ny*nz,0.0f);
    std::vector<float> Fwz(nx*ny*nz,0.0f);

    readSTL(STLname, STLnormal, STLv0, STLv1, STLv2);
    const int nSTL = STLnormal[0].size();
    std::cout << "Number of elements of STL: " << nSTL << "\n" << std::endl;

    std::vector<std::vector<float> > STLc(3, std::vector<float>(nSTL));
    for(int i = 0; i < nSTL; i++)
    {
        STLv0[0][i] = STLv0[0][i]/L -0.5;
        STLv0[1][i] = STLv0[1][i]/L -0.5;
        STLv0[2][i] = STLv0[2][i]/L -0.5;
        STLv1[0][i] = STLv1[0][i]/L -0.5;
        STLv1[1][i] = STLv1[1][i]/L -0.5;
        STLv1[2][i] = STLv1[2][i]/L -0.5;
        STLv2[0][i] = STLv2[0][i]/L -0.5;
        STLv2[1][i] = STLv2[1][i]/L -0.5;
        STLv2[2][i] = STLv2[2][i]/L -0.5;

        STLc[0][i] = (STLv0[0][i]+STLv1[0][i]+STLv2[0][i])/3.f;
        STLc[1][i] = (STLv0[1][i]+STLv1[1][i]+STLv2[1][i])/3.f;
        STLc[2][i] = (STLv0[2][i]+STLv1[2][i]+STLv2[2][i])/3.f;
    }
    //--

    //-- Calculation of sdf for the boundaries defined by STL
    const float sdfIni = 10000.f;
    const float qfIni = 0.5f;
    std::vector<float> sdf(nx*ny*nz,sdfIni);
    std::vector<float> qf(19*nx*ny*nz,qfIni);
    std::vector<unsigned char> solid(nx*ny*nz,0);
    std::vector<unsigned char> neiSolid(nx*ny*nz,0);
    float dr = 10.f;
    const float p = 7.f;
    int drn = 10;
    
    std::vector<int> nearSTL;
    for(int iSTL = 0; iSTL < nSTL; iSTL++)
    {
        int i = int(STLc[0][iSTL]);
        int j = int(STLc[1][iSTL]);
        int k = int(STLc[2][iSTL]);
        // int k = 0;

        if(0 < i && i < nx-1 && 0 < j && j < ny-1 && 0 < k && k < nz-1)
        {
            for(int ii = -drn; ii <= drn; ii++)
            {
                int iNear = i + ii;
                if(0 <= iNear && iNear <= nx-1)
                {
                    for(int jj = -drn; jj <= drn; jj++)
                    {
                        int jNear = j + jj;
                        if(0 <= jNear && jNear <= ny-1)
                        {
                            for(int kk = -drn; kk <= drn; kk++)
                            {
                                int kNear = k + kk;
                                if(0 <= kNear && kNear <= nz-1)
                                {
                                    int ic = index1d(iNear,jNear,kNear,nx,ny);
                                    nearSTL.push_back(ic);
                                }
                            }
                        }
                    }
                }   
            }
        }
    }
    std::sort(nearSTL.begin(), nearSTL.end());
    nearSTL.erase(std::unique(nearSTL.begin(), nearSTL.end()), nearSTL.end());


    std::cout << "Number of nearSTL: " << nearSTL.size() << std::endl;
    for(int iNSTL = 0; iNSTL < nearSTL.size(); iNSTL++)
    {
        int ic = nearSTL[iNSTL];
        sdf[ic] = 0.f;

        int i = ic2i(ic,nx,ny);
        int j = ic2j(ic,nx,ny);
        int k = ic2k(ic,nx,ny);
        
        float sumD = 0.f;
        for(int iSTL = 0; iSTL < nSTL; iSTL++)
        {
            const float r = sqrt
                            (
                                pow(float(i) -STLc[0][iSTL],2.f)
                                +pow(float(j) -STLc[1][iSTL],2.f)
                                +pow(float(k) -STLc[2][iSTL],2.f)
                            );
            if(r < dr)
            {
                const float sd = (float(i) -STLc[0][iSTL])*STLnormal[0][iSTL]
                                +(float(j) -STLc[1][iSTL])*STLnormal[1][iSTL]
                                +(float(k) -STLc[2][iSTL])*STLnormal[2][iSTL];
                sumD += pow(r,-p);
                sdf[ic] += sd*pow(r,-p);
            }
        }
        sdf[ic] = sumD != 0.f ? sdf[ic]/sumD : sdfIni;
        if(sdf[ic] < 0.f)
        {
            solid[ic] = 1;
        }
    }

    for(int j = 0; j < ny; j++)
    {
        for(int k = 0; k < nz; k++)
        {
            std::vector<int> tmp;
            for(int i = 0; i < nx; i++)
            {
                int ic = index1d(i,j,k,nx,ny);
                if(sdf[ic] == sdfIni)
                {
                    tmp.push_back(ic);
                }
                else if(tmp.size() != 0 && (sdf[ic] < 0.0 || i == nx-1))
                {
                    for(int id = 0; id < tmp.size(); id++)
                    {
                        solid[tmp[id]] = 1;
                    }
                    tmp.clear();
                }
                else if(tmp.size() != 0 && sdf[ic] != sdfIni && sdf[ic] > 0.0)
                {
                    tmp.clear();
                }
            }
        }
    }

    for(int ic = 0; ic < nx*ny*nz; ic++)
    {
        int i = ic2i(ic,nx,ny);
        int j = ic2j(ic,nx,ny);
        int k = ic2k(ic,nx,ny);
        
        if(i != 0 && i != nx-1 && j != 0 && j != ny-1 && k != 0 && k != nz-1)
        {
            for(int q = 0; q < 19; q++)
            {
                int qic = q*elements +ic;
                const float sdf0 = sdf[ic];
                const int upID = upwindID(q,i,j,k,nx,ny,nz);
                
                const float sdf1 = sdf[upID];
                if(solid[ic] == 0 && solid[upID] == 1)
                {
                    neiSolid[ic] = 1;
                    qf[qic] = abs(sdf0)/(abs(sdf0)+abs(sdf1));
                }
            }
        }
    }
    //--

    //-- Calculation of sdf for the boundaries of the calculation region
    dr = 3.f;
    drn = 3;
    std::vector<std::vector<float> > bWallsC(3);
    std::vector<std::vector<float> > bWallsNormal(3);
    for(int j = 0; j < ny; j++)
    {
        for(int k = 0; k < nz; k++)
        {
            int icStart = index1d(0,j,k,nx,ny);
            int icEnd = index1d(nx-1,j,k,nx,ny);
            if(boundary1[icStart] == 1)
            {
                bWallsC[0].push_back(-0.5f);
                bWallsC[1].push_back(j);
                bWallsC[2].push_back(k);
                bWallsNormal[0].push_back(1.f);
                bWallsNormal[1].push_back(0.f);
                bWallsNormal[2].push_back(0.f);
            }
            if(boundary1[icEnd] == 1)
            {
                bWallsC[0].push_back(nx-0.5f);
                bWallsC[1].push_back(j);
                bWallsC[2].push_back(k);
                bWallsNormal[0].push_back(-1.f);
                bWallsNormal[1].push_back(0.f);
                bWallsNormal[2].push_back(0.f);
            }
        }
    }
    for(int i = 0; i < nx; i++)
    {
        for(int k = 0; k < nz; k++)
        {
            int icStart = index1d(i,0,k,nx,ny);
            int icEnd = index1d(i,ny-1,k,nx,ny);
            if(boundary2[icStart] == 1)
            {
                bWallsC[0].push_back(i);
                bWallsC[1].push_back(-0.5f);
                bWallsC[2].push_back(k);
                bWallsNormal[0].push_back(0.f);
                bWallsNormal[1].push_back(1.f);
                bWallsNormal[2].push_back(0.f);
            }
            if(boundary2[icEnd] == 1)
            {
                bWallsC[0].push_back(i);
                bWallsC[1].push_back(ny-0.5f);
                bWallsC[2].push_back(k);
                bWallsNormal[0].push_back(0.f);
                bWallsNormal[1].push_back(-1.f);
                bWallsNormal[2].push_back(0.f);
            }
        }
    }
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            int icStart = index1d(i,j,0,nx,ny);
            int icEnd = index1d(i,j,nz-1,nx,ny);
            if(boundary3[icStart] == 1)
            {
                bWallsC[0].push_back(i);
                bWallsC[1].push_back(j);
                bWallsC[2].push_back(-0.5f);
                bWallsNormal[0].push_back(0.f);
                bWallsNormal[1].push_back(0.f);
                bWallsNormal[2].push_back(1.f);
            }
            if(boundary3[icEnd] == 1)
            {
                bWallsC[0].push_back(i);
                bWallsC[1].push_back(j);
                bWallsC[2].push_back(nz-0.5f);
                bWallsNormal[0].push_back(0.f);
                bWallsNormal[1].push_back(0.f);
                bWallsNormal[2].push_back(-1.f);
            }
        }
    }

    int nBwalls = bWallsC[0].size();
    std::cout << "nBwalls: " << nBwalls << std::endl;
    std::vector<int> nearBwalls;
    for(int iBwalls = 0; iBwalls < nBwalls; iBwalls++)
    {
        int i = int(bWallsC[0][iBwalls]);
        int j = int(bWallsC[1][iBwalls]);
        int k = int(bWallsC[2][iBwalls]);
        
        for(int ii = -drn; ii <= drn; ii++)
        {
            int iNear = i + ii;
            if(0 <= iNear && iNear <= nx-1)
            {
                for(int jj = -drn; jj <= drn; jj++)
                {
                    int jNear = j + jj;
                    if(0 <= jNear && jNear <= ny-1)
                    {
                        for(int kk = -drn; kk <= drn; kk++)
                        {
                            int kNear = k + kk;
                            if(0 <= kNear && kNear <= nz-1)
                            {
                                int ic = index1d(iNear,jNear,kNear,nx,ny);
                                nearBwalls.push_back(ic);
                            }
                        }
                    }
                }
            }   
        }
    }
    std::sort(nearBwalls.begin(), nearBwalls.end());
    nearBwalls.erase(std::unique(nearBwalls.begin(), nearBwalls.end()), nearBwalls.end());


    std::cout << "Number of nearBwalls: " << nearBwalls.size() << std::endl;
    for(int iNBwalls = 0; iNBwalls < nearBwalls.size(); iNBwalls++)
    {
        int ic = nearBwalls[iNBwalls];
        sdf[ic] = 0.f;

        int i = ic2i(ic,nx,ny);
        int j = ic2j(ic,nx,ny);
        int k = ic2k(ic,nx,ny);
        
        float sumD = 0.f;
        for(int iBwalls = 0; iBwalls < nBwalls; iBwalls++)
        {
            const float r = sqrt
                            (
                                pow(float(i) -bWallsC[0][iBwalls],2.f)
                                +pow(float(j) -bWallsC[1][iBwalls],2.f)
                                +pow(float(k) -bWallsC[2][iBwalls],2.f)
                            );
            if(r < dr)
            {
                const float sd = (float(i) -bWallsC[0][iBwalls])*bWallsNormal[0][iBwalls]
                                +(float(j) -bWallsC[1][iBwalls])*bWallsNormal[1][iBwalls]
                                +(float(k) -bWallsC[2][iBwalls])*bWallsNormal[2][iBwalls];
                sumD += pow(r,-p);
                sdf[ic] += sd*pow(r,-p);
            }
        }
        sdf[ic] = sumD != 0.f ? sdf[ic]/sumD : sdfIni;
    }
    //--



    try
    {
    // Create a context
    cl::Context context(DEVICE);

    cl::Device device = context.getInfo<CL_CONTEXT_DEVICES>()[0];

    std::cout << std::endl << "Using OpenCL device: "
                << device.getInfo<CL_DEVICE_NAME>() << std::endl
                << device.getInfo<CL_DEVICE_VERSION>() << std::endl
                << device.getInfo<CL_DEVICE_OPENCL_C_VERSION>() << std::endl << std::endl;
    

    // Load in kernel source, creating and building a program object for the context
    cl::Program program(context, util::loadProgram("kernel.cl"), true);

    
    cl::CommandQueue queue(context);

    cl::Buffer f_d(context, f.begin(), f.end(), true);
    
    cl::Buffer fTmp_d(context, fTmp.begin(), fTmp.end(), true);
    
    cl::Buffer u0_d(context, u0.begin(), u0.end(), true);
    cl::Buffer v0_d(context, v0.begin(), v0.end(), true);
    cl::Buffer w0_d(context, w0.begin(), w0.end(), true);
    cl::Buffer Fwx_d(context, Fwx.begin(), Fwx.end(), true);
    cl::Buffer Fwy_d(context, Fwy.begin(), Fwy.end(), true);
    cl::Buffer Fwz_d(context, Fwz.begin(), Fwz.end(), true);
    cl::Buffer boundary1_d(context, boundary1.begin(), boundary1.end(), true);
    cl::Buffer boundary2_d(context, boundary2.begin(), boundary2.end(), true);
    cl::Buffer boundary3_d(context, boundary3.begin(), boundary3.end(), true);
    cl::Buffer rho_d(context, rho.begin(), rho.end(), true);
    cl::Buffer sdf_d(context, sdf.begin(), sdf.end(), true);
    cl::Buffer solid_d(context, solid.begin(), solid.end(), true);
    cl::Buffer neiSolid_d(context, neiSolid.begin(), neiSolid.end(), true);

    cl::Buffer tauSGS_d(context, tauSGS.begin(), tauSGS.end(), true);
    
    // Create the kernel functor of streamingCollision
    cl::KernelFunctor
    <
        cl::Buffer, cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer,
        const unsigned,
        const float,
        const float,
        const float,
        const int, const int, const int
    > k_streamingCollision(program, "k_streamingCollision");  

    std::chrono::system_clock::time_point start;
    std::chrono::system_clock::time_point end;
    start = std::chrono::system_clock::now();

    // Time marching
    for(int nt = startTimeStep; nt <= endTimeStep; nt++)
    {
        #include "write.hpp"
        {
        util::Timer timer;
        k_streamingCollision
        (
            cl::EnqueueArgs(queue,cl::NDRange(elements)),
            f_d, fTmp_d,
            boundary1_d, boundary2_d, boundary3_d,
            sdf_d, solid_d, neiSolid_d,
            u0_d, v0_d, w0_d,
            Fwx_d, Fwy_d, Fwz_d,
            tauSGS_d,
            elements,
            omega,
            dpdx,
            rho_av,
            nx, ny, nz
        );
        queue.finish();
        double rtime = static_cast<double>(timer.getTimeMicroseconds()) / 1000.0;
        // printf("\nThe kernel of streamingCollision ran in %lf m seconds\n", rtime);
        }
        
        std::swap(fTmp_d,f_d);
    }
    end = std::chrono::system_clock::now();
    double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() *1e-6);
    std::cout << "Execution time: " << time << " (s)" << std::endl;
    std::cout << "Speed: " << double(endTimeStep)*double(nx*ny*nz)/time*1e-6 << " (MLUPS)" << std::endl;
    }
    catch (cl::BuildError error)
    {
      std::string log = error.getBuildLog()[0].second;
      std::cerr << std::endl << "Build failed:" << std::endl << log << std::endl;
    }
    catch (cl::Error err) {
        std::cout << "Exception\n";
        std::cerr
            << "ERROR: "
            << err.what()
            << "("
            << err_code(err.err())
           << ")"
           << std::endl;
    }

    return EXIT_SUCCESS;
}
