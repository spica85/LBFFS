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
    int nSTL;
    std::vector<std::vector<float> > STLc(3);

    std::vector<float> Fwx(nx*ny*nz,0.0f);
    std::vector<float> Fwy(nx*ny*nz,0.0f);
    std::vector<float> Fwz(nx*ny*nz,0.0f);

    readSTL(STLname, STLnormal, STLv0, STLv1, STLv2, nSTL, STLc, L);
    std::cout << "Number of elements of STL: " << nSTL << std::endl;
    //--

    //-- Calculations of sdf, solid, and qf for the boundaries defined by STL
    const float sdfIni = 10000.f;
    std::vector<float> sdf(nx*ny*nz,sdfIni);
    float dr = 10.f;
    const float p = 7.f;
    setSDF(sdf, sdfIni, dr, p, STLc, STLnormal, nx, ny, nz, false);
    
    std::vector<unsigned char> solid(elements,0);
    setSolid(solid, sdf, sdfIni, nx, ny, nz);

    const float qfIni = 0.5f;
    std::vector<float> qf(19*nx*ny*nz,qfIni);
    std::vector<unsigned char> neiSolid(nx*ny*nz,0);
    setQf(qf, neiSolid, sdf, solid, nx, ny, nz);
    //--

    //-- Calculation of sdf for the boundaries of the calculation region
    dr = 3.f;
    std::vector<std::vector<float> > bWallsC(3);
    std::vector<std::vector<float> > bWallsNormal(3);
    setBwalls(bWallsC, bWallsNormal, boundary1, boundary2, boundary3, nx, ny, nz);
    setSDF(sdf, sdfIni, dr, p, bWallsC, bWallsNormal, nx, ny, nz, true);
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
