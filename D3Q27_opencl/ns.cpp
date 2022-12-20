#include <boost/algorithm/string.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <chrono>

#define _USE_MATH_DEFINES
#include <cmath>

#include "D3Q27.hpp"
#include "input.hpp"
#include "walls.hpp"

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
    float Lx;
    float uMax;
    float rho0;
    float Re;
    float nu;
    float U0;
    float dpdx;
    float LES;
    bool forceCoeffs;
    float Dref;
    float omegaB;
    float spzWidth;
    float clim;
    bool invertFluidSolid;
    float uIni;
    float vIni;
    float wIni;

    input(restart, Fwrite, writeBinary, startTimeStep, endTimeStep, nextOutTime, outInterval, nx, ny, nz, Lx, uMax, rho0, U0, nu, dpdx, LES, forceCoeffs, Dref, omegaB, spzWidth, clim, invertFluidSolid, uIni, vIni, wIni);

    //-- D3Q27 model
    const std::vector<float> wt = setWt();

    const std::vector<float> cx = setCx();
    const std::vector<float> cy = setCy();
    const std::vector<float> cz = setCz();

    std::vector<float> f(27*nx*ny*nz);
    std::vector<float> fTmp(27*nx*ny*nz);
    //--
    
    std::vector<float> rho(nx*ny*nz);
    std::vector<float> u(nx*ny*nz);
    std::vector<float> v(nx*ny*nz);
    std::vector<float> w(nx*ny*nz);
    float rho_av = 1.0;
    std::vector<float> uMean(nx*ny*nz,0.f);
    std::vector<float> vMean(nx*ny*nz,0.f);
    std::vector<float> wMean(nx*ny*nz,0.f);
    float timeMean = 0.f;

    std::vector<float> tauSGS(nx*ny*nz);

    std::vector<float> u0(nx*ny*nz,0.0f);
    std::vector<float> v0(nx*ny*nz,0.0f);
    std::vector<float> w0(nx*ny*nz,0.0f);
    std::vector<int> boundary(nx*ny*nz,0);
    std::vector<int> boundary1(nx*ny*nz,0);
    std::vector<int> boundary2(nx*ny*nz,0);
    std::vector<int> boundary3(nx*ny*nz,0);

    const unsigned elements = nx*ny*nz;
    const unsigned qElements = 27*elements;
    //--
    const float L = Lx/float(nx-1);
    const float c = uMax/U0; //Representative velocity (m/s)
    std::cout << "L: " << L << ", c: " << c << std::endl;
    
    Re = uMax*Lx/nu;
    std::cout << "Re: " << Re << std::endl;

    nu = nu/(L*c);
    dpdx *= L/(c*c);
    uIni /= c;
    vIni /= c;
    wIni /= c;

    std::cout << "Size: "
              << "(" << L/2 << ", " << L/2 << ", " << L/2 << "), "
              <<  "(" << L/2+L*(nx-1) << ", " << L/2+L*(ny-1) << ", " << L/2+L*(nz-1) << ")"
              << std::endl;

    std::cout << "dpdx = " << dpdx << std::endl;
    std::cout << "nu = " << nu << std::endl;


    const float deltaT = L/c; //Dimensional time step (s)
    const float omega0 = 1.0/(3.0*nu +0.5);
    std::vector<float> omega(nx*ny*nz,omega0);

    std::cout << "tau = " << 1.0/omega0 << ", taulim = " << 0.5f +uMax/c/8.0f << std::endl;
    std::cout << "Maximum Ma = " << (uMax/c)*sqrt(3.f) << std::endl;

    //-- Settings for boundary conditions
    int nxMinBoundary1 = 0; int nxMinBoundary2 = 0; int nxMinBoundary3 = 0;
    float nxMinU0 = 0.f; float nxMinV0 = 0.f; float nxMinW0 = 0.f;
    int nxMaxBoundary1 = 0; int nxMaxBoundary2 = 0; int nxMaxBoundary3 = 0;
    float nxMaxU0 = 0.f; float nxMaxV0 = 0.f; float nxMaxW0 = 0.f;
    int nyMinBoundary1 = 0; int nyMinBoundary2 = 0; int nyMinBoundary3 = 0;
    float nyMinU0 = 0.f; float nyMinV0 = 0.f; float nyMinW0 = 0.f;
    int nyMaxBoundary1 = 0; int nyMaxBoundary2 = 0; int nyMaxBoundary3 = 0;
    float nyMaxU0 = 0.f; float nyMaxV0 = 0.f; float nyMaxW0 = 0.f;
    int nzMinBoundary1 = 0; int nzMinBoundary2 = 0; int nzMinBoundary3 = 0;
    float nzMinU0 = 0.f; float nzMinV0 = 0.f; float nzMinW0 = 0.f;
    int nzMaxBoundary1 = 0; int nzMaxBoundary2 = 0; int nzMaxBoundary3 = 0;
    float nzMaxU0 = 0.f; float nzMaxV0 = 0.f; float nzMaxW0 = 0.f;
    inputBoundaryConditions(
        nxMinBoundary1, nxMinBoundary2, nxMinBoundary3,
        nxMinU0, nxMinV0, nxMinW0,
        nxMaxBoundary1, nxMaxBoundary2, nxMaxBoundary3,
        nxMaxU0, nxMaxV0, nxMaxW0,
        nyMinBoundary1, nyMinBoundary2, nyMinBoundary3,
        nyMinU0, nyMinV0, nyMinW0,
        nyMaxBoundary1, nyMaxBoundary2, nyMaxBoundary3,
        nyMaxU0, nyMaxV0, nyMaxW0,
        nzMinBoundary1, nzMinBoundary2, nzMinBoundary3,
        nzMinU0, nzMinV0, nzMinW0,
        nzMaxBoundary1, nzMaxBoundary2, nzMaxBoundary3,
        nzMaxU0, nzMaxV0, nzMaxW0,
        c
        );
    
    #include "setBoundaryConditions.hpp"

    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                if(i == 0 || i == nx-1 || j == 0 || j == ny-1 || (nz-1 != 0 && (k == 0 || k == nz-1)))
                {
                    int ic = index1d(i,j,k,nx,ny);
                    boundary[ic] = 1;
                }
            }
        }   
    }

    if(restart)
    {
        #include "restart.hpp"
    }
    else
    {
        #include "initialization.hpp"        
    }

    //-- Reading and settings of STL 
    STL wallSTL("walls.stl",L);
    std::vector<float> Fwx(nx*ny*nz,0.0f);
    std::vector<float> Fwy(nx*ny*nz,0.0f);
    std::vector<float> Fwz(nx*ny*nz,0.0f);

    std::vector<unsigned char> solid(elements,0);
    setSolid(solid,wallSTL,nx,ny,nz,invertFluidSolid);

    std::vector<unsigned char> neiSolid(nx*ny*nz,0);
    setNeiSolid(neiSolid, solid, nx, ny, nz, 27);

    //-- Calculations of sdf, solid, and qf for the boundaries defined by STL
    const float sdfIni = 10000.f;
    std::vector<float> sdf(nx*ny*nz,sdfIni);
    float dr = 10.f;
    const float p = 7.f;

    for(int iPatch = 0; iPatch < wallSTL.patch.size(); iPatch++)
    {
        const std::unique_ptr<STLpatch>& patch = wallSTL.patch[iPatch];
        if(LES)
        {
            setSDF(sdf, sdfIni, dr, p, patch, nx, ny, nz, false, neiSolid);
        }
    }

    const float qfIni = 0.5f;
    std::vector<float> qf(27*nx*ny*nz,qfIni);
    setQf(qf, neiSolid, sdf, nx, ny, nz, 27);
    //--

    //-- Calculation of sdf for the boundaries of the calculation region
    dr = 3.f;
    if(LES)
    {
        setSDFnearBwalls(sdf,dr,boundary1,boundary2,boundary3,nx,ny,nz);
    }
    //--

    //-- Reading and settings of moving STL
    const std::string movingSTLname("movingWalls.stl");
    std::vector<std::vector<float> > movingSTLnormal(3);
    std::vector<std::vector<float> > movingSTLv0(3);
    std::vector<std::vector<float> > movingSTLv1(3);
    std::vector<std::vector<float> > movingSTLv2(3);
    int nMovingSTL;
    std::vector<std::vector<float> > movingSTLc(3);

    const int isReadMovingWalls = readSTL(movingSTLname, movingSTLnormal, movingSTLv0, movingSTLv1, movingSTLv2, nMovingSTL, movingSTLc, L);
    // -- Exception handling to avoid zero buffer size
    if(isReadMovingWalls == 0)
    {
        nMovingSTL = 1;
        movingSTLc[0].push_back(0.f);
        movingSTLc[1].push_back(0.f);
        movingSTLc[2].push_back(0.f);
    }
    // --

    std::vector<float> movingSTLcList(3*nMovingSTL);
    for(int iMSTL = 0; iMSTL < nMovingSTL; iMSTL++)
    {
        movingSTLcList[iMSTL*3] = movingSTLc[0][iMSTL];
        movingSTLcList[iMSTL*3+1] = movingSTLc[1][iMSTL];
        movingSTLcList[iMSTL*3+2] = movingSTLc[2][iMSTL];
    }
    std::vector<float> GxMovingWall(nMovingSTL);
    std::vector<float> GyMovingWall(nMovingSTL);
    std::vector<float> GzMovingWall(nMovingSTL);

    std::vector<float> GxIBM(elements,0.f);
    std::vector<float> GyIBM(elements,0.f);
    std::vector<float> GzIBM(elements,0.f);

    float uMovingTrans = 0.f;
    float vMovingTrans = 0.f;
    float wMovingTrans = 0.f;

    float rotOmega = 0.f;
    float rotX = 0.f;
    float rotY = 0.f;
    float rotZ = 0.f;
    float rotAxisX = 0.f;
    float rotAxisY = 0.f;
    float rotAxisZ = 0.f;
    if(nMovingSTL != 1)
    {
        readMotions(uMovingTrans,vMovingTrans,wMovingTrans,rotOmega,rotX,rotY,rotZ,rotAxisX,rotAxisY,rotAxisZ,c,L);
    }

    // --
        STL inletSTL("inlets.stl",L);
        for(int iPatch = 0; iPatch < inletSTL.patch.size(); iPatch++)
        {
            const std::unique_ptr<STLpatch>& patch = inletSTL.patch[iPatch];
            int patchBoundary1;
            int patchBoundary2;
            int patchBoundary3;
            float patchU0;
            float patchV0;
            float patchW0;
            std::vector<int> points;

            inputBoundaryCondition(patch->patchName,patchBoundary1,patchBoundary2,patchBoundary3,patchU0,patchV0,patchW0,c);
            int spzWidth = inputSpzWidth(patch->patchName, L);
            float tau0 = 0.55f;

            patch->findPoints(points,nx,ny,nz);

            for(int id = 0; id < points.size(); id++)
            {
                int ic = points[id];
                int i = ic2i(ic,nx,ny);
                int j = ic2j(ic,nx,ny);
                int k = ic2k(ic,nx,ny);

                boundary1[ic] = patchBoundary1;
                boundary2[ic] = patchBoundary2;
                boundary3[ic] = patchBoundary3;
                u0[ic] = patchU0;
                v0[ic] = patchV0;
                w0[ic] = patchW0;

                if(patchBoundary1 != BCnameToNum("Cyclic"))
                {
                    for(int iSpz = fmax(i-spzWidth,0); iSpz < fmin(i+spzWidth, nx); iSpz++)
                    {
                        int icSpz = index1d(iSpz,j,k,nx,ny);
                        float fx = (i-spzWidth < 0) ? float(iSpz)/spzWidth : float((nx-1)-iSpz)/spzWidth;
                        fx = 1.f -cos(0.5f*M_PI*fx);
                        float tau = fx*(1.f/omega0) + (1.f-fx)*tau0;
                        omega[icSpz] = 1.f/tau;
                    }
                }
                if(patchBoundary2 != BCnameToNum("Cyclic"))
                {   
                    for(int jSpz = fmax(j-spzWidth,0); jSpz < fmin(j+spzWidth, ny); jSpz++)
                    {
                        int icSpz = index1d(i,jSpz,k,nx,ny);
                        float fx = (j-spzWidth < 0) ? float(jSpz)/spzWidth : float((ny-1)-jSpz)/spzWidth;
                        fx = 1.f -cos(0.5f*M_PI*fx);
                        float tau = fx*(1.f/omega0) + (1.f-fx)*tau0;
                        omega[icSpz] = 1.f/tau;
                    }
                }
                if(patchBoundary3 != BCnameToNum("Cyclic"))
                {
                    for(int kSpz = fmax(k-spzWidth,0); kSpz < fmin(k+spzWidth, nz); kSpz++)
                    {
                        int icSpz = index1d(i,j,kSpz,nx,ny);
                        float fx = (k-spzWidth < 0) ? float(kSpz)/spzWidth : float((nz-1)-kSpz)/spzWidth;
                        fx = 1.f -cos(0.5f*M_PI*fx);
                        float tau = fx*(1.f/omega0) + (1.f-fx)*tau0;
                        omega[icSpz] = 1.f/tau;
                    }
                }
            }
        }
    // --

    // --
    STL outletSTL("outlets.stl",L);
    for(int iPatch = 0; iPatch < outletSTL.patch.size(); iPatch++)
    {
        const std::unique_ptr<STLpatch>& patch = outletSTL.patch[iPatch];
        int patchBoundary1;
        int patchBoundary2;
        int patchBoundary3;
        float patchU0;
        float patchV0;
        float patchW0;
        std::vector<int> points;

        inputBoundaryCondition(patch->patchName,patchBoundary1,patchBoundary2,patchBoundary3,patchU0,patchV0,patchW0,c);
        int spzWidth = inputSpzWidth(patch->patchName, L);
        float tau0 = 0.55f;

        patch->findPoints(points,nx,ny,nz);

        for(int id = 0; id < points.size(); id++)
        {
            int ic = points[id];
            int i = ic2i(ic,nx,ny);
            int j = ic2j(ic,nx,ny);
            int k = ic2k(ic,nx,ny);

            boundary1[ic] = patchBoundary1;
            boundary2[ic] = patchBoundary2;
            boundary3[ic] = patchBoundary3;
            u0[ic] = patchU0;
            v0[ic] = patchV0;
            w0[ic] = patchW0;

            if(patchBoundary1 != BCnameToNum("Cyclic"))
            {
                for(int iSpz = fmax(i-spzWidth,0); iSpz < fmin(i+spzWidth, nx); iSpz++)
                {
                    int icSpz = index1d(iSpz,j,k,nx,ny);
                    float fx = (i-spzWidth < 0) ? float(iSpz)/spzWidth : float((nx-1)-iSpz)/spzWidth;
                    fx = 1.f -cos(0.5f*M_PI*fx);
                    float tau = fx*(1.f/omega0) + (1.f-fx)*tau0;
                    omega[icSpz] = 1.f/tau;
                }
            }
            if(patchBoundary2 != BCnameToNum("Cyclic"))
            {   
                for(int jSpz = fmax(j-spzWidth,0); jSpz < fmin(j+spzWidth, ny); jSpz++)
                {
                    int icSpz = index1d(i,jSpz,k,nx,ny);
                    float fx = (j-spzWidth < 0) ? float(jSpz)/spzWidth : float((ny-1)-jSpz)/spzWidth;
                    fx = 1.f -cos(0.5f*M_PI*fx);
                    float tau = fx*(1.f/omega0) + (1.f-fx)*tau0;
                    omega[icSpz] = 1.f/tau;
                }
            }
            if(patchBoundary3 != BCnameToNum("Cyclic"))
            {
                for(int kSpz = fmax(k-spzWidth,0); kSpz < fmin(k+spzWidth, nz); kSpz++)
                {
                    int icSpz = index1d(i,j,kSpz,nx,ny);
                    float fx = (k-spzWidth < 0) ? float(kSpz)/spzWidth : float((nz-1)-kSpz)/spzWidth;
                    fx = 1.f -cos(0.5f*M_PI*fx);
                    float tau = fx*(1.f/omega0) + (1.f-fx)*tau0;
                    omega[icSpz] = 1.f/tau;
                }
            }
        }
    }
    // --



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
    // cl::Program program(context, util::loadProgram("kernel.cl"), true);
    cl::Program program(context, util::loadProgram("kernel.cl"), false);
    // program.build("-cl-mad-enable -cl-fast-relaxed-math -cl-no-signed-zeros");
    program.build();
    
    

    
    cl::CommandQueue queue(context);

    cl::Buffer f_d(context, f.begin(), f.end(), true);
    
    cl::Buffer fTmp_d(context, fTmp.begin(), fTmp.end(), true);
    
    cl::Buffer u0_d(context, u0.begin(), u0.end(), true);
    cl::Buffer v0_d(context, v0.begin(), v0.end(), true);
    cl::Buffer w0_d(context, w0.begin(), w0.end(), true);
    cl::Buffer Fwx_d(context, Fwx.begin(), Fwx.end(), true);
    cl::Buffer Fwy_d(context, Fwy.begin(), Fwy.end(), true);
    cl::Buffer Fwz_d(context, Fwz.begin(), Fwz.end(), true);
    cl::Buffer boundary_d(context, boundary.begin(), boundary.end(), true);
    cl::Buffer boundary1_d(context, boundary1.begin(), boundary1.end(), true);
    cl::Buffer boundary2_d(context, boundary2.begin(), boundary2.end(), true);
    cl::Buffer boundary3_d(context, boundary3.begin(), boundary3.end(), true);
    cl::Buffer sdf_d(context, sdf.begin(), sdf.end(), true);
    cl::Buffer solid_d(context, solid.begin(), solid.end(), true);
    cl::Buffer neiSolid_d(context, neiSolid.begin(), neiSolid.end(), true);

    cl::Buffer tauSGS_d(context, tauSGS.begin(), tauSGS.end(), true);

    cl::Buffer rho_d(context, rho.begin(), rho.end(), true);
    cl::Buffer u_d(context, u.begin(), u.end(), true);
    cl::Buffer v_d(context, v.begin(), v.end(), true);
    cl::Buffer w_d(context, w.begin(), w.end(), true);
    cl::Buffer omega_d(context, omega.begin(), omega.end(), true);
    cl::Buffer movingSTLcList_d(context, movingSTLcList.begin(), movingSTLcList.end(), true);
    cl::Buffer GxMovingWall_d(context, GxMovingWall.begin(), GxMovingWall.end(), true);
    cl::Buffer GyMovingWall_d(context, GyMovingWall.begin(), GyMovingWall.end(), true);
    cl::Buffer GzMovingWall_d(context, GzMovingWall.begin(), GzMovingWall.end(), true);
    cl::Buffer GxIBM_d(context, GxIBM.begin(), GxIBM.end(), true);
    cl::Buffer GyIBM_d(context, GyIBM.begin(), GyIBM.end(), true);
    cl::Buffer GzIBM_d(context, GzIBM.begin(), GzIBM.end(), true);
    


    // Create the kernel functor of streamingCollision
    cl::KernelFunctor
    <
        cl::Buffer, cl::Buffer,
        cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer,
        cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer,
        const unsigned,
        const float,
        const float,
        const int, const int, const int,
        const float,
        const int,
        const float,
        const float,
        const float
    > k_streamingCollision(program, "k_streamingCollision");  

    // Create the kernel functor of Gwall
    cl::KernelFunctor
    <
        cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        const int,
        const int, const int, const int,
        const float, const float, const float,
        const float, const float, const float,
        const float, const float, const float,
        const float
    > k_Gwall(program, "k_Gwall");

    // Create the kernel functor of Gibm
    cl::KernelFunctor
    <
        cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        const int,
        const int, const int, const int,
        const float, const float, const float,
        const float, const float, const float,
        const float, const float, const float,
        const float
    > k_Gibm(program, "k_Gibm");

    // Create the kernel functor of Force
    cl::KernelFunctor
    <
        cl::Buffer,
        cl::Buffer,
        cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        const unsigned
    > k_Force(program, "k_Force");

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
            boundary_d,
            boundary1_d, boundary2_d, boundary3_d,
            sdf_d, solid_d, neiSolid_d,
            u0_d, v0_d, w0_d,
            Fwx_d, Fwy_d, Fwz_d,
            tauSGS_d,
            rho_d,
            u_d, v_d, w_d,
            omega_d,
            GxIBM_d, GyIBM_d, GzIBM_d,
            elements,
            dpdx,
            rho_av,
            nx, ny, nz,
            LES,
            isReadMovingWalls,
            omegaB,
            spzWidth,
            clim
        );
        queue.finish();
        double rtime = static_cast<double>(timer.getTimeMicroseconds()) / 1000.0;
        // printf("\nThe kernel of streamingCollision ran in %lf m seconds\n", rtime);
        }
    
        if(nMovingSTL > 1)
        {
            {
            util::Timer timer;
            k_Gwall
            (
                cl::EnqueueArgs(queue,cl::NDRange(nMovingSTL)),
                rho_d,
                u_d, v_d, w_d,
                movingSTLcList_d,
                GxMovingWall_d, GyMovingWall_d, GzMovingWall_d,
                nMovingSTL,
                nx, ny, nz,
                uMovingTrans, vMovingTrans, wMovingTrans,
                rotX, rotY, rotZ,
                rotAxisX, rotAxisY, rotAxisZ,
                rotOmega
            );
            queue.finish();
            double rtime = static_cast<double>(timer.getTimeMicroseconds()) / 1000.0;
            // printf("\nThe kernel of Gwall ran in %lf m seconds\n", rtime);
            }

            {
            util::Timer timer;
            k_Gibm
            (
                cl::EnqueueArgs(queue,cl::NDRange(nMovingSTL)),
                rho_d,
                GxIBM_d, GyIBM_d, GzIBM_d,
                movingSTLcList_d,
                GxMovingWall_d, GyMovingWall_d, GzMovingWall_d,
                nMovingSTL,
                nx, ny, nz,
                uMovingTrans, vMovingTrans, wMovingTrans,
                rotX, rotY, rotZ,
                rotAxisX, rotAxisY, rotAxisZ,
                rotOmega
            );
            queue.finish();
            double rtime = static_cast<double>(timer.getTimeMicroseconds()) / 1000.0;
            // printf("\nThe kernel of Gibm ran in %lf m seconds\n", rtime);
            }

            {
            util::Timer timer;
            k_Force
            (
                cl::EnqueueArgs(queue,cl::NDRange(elements)),
                fTmp_d,
                solid_d,
                rho_d,
                GxIBM_d, GyIBM_d, GzIBM_d,
                elements
            );
            queue.finish();
            double rtime = static_cast<double>(timer.getTimeMicroseconds()) / 1000.0;
            // printf("\nThe kernel of Force ran in %lf m seconds\n", rtime);
            }
        }
        std::swap(fTmp_d,f_d);
    }
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
