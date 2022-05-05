#include <boost/algorithm/string.hpp>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <chrono>

#include "LKS.hpp"
#include "input.hpp"

// For OpenCL
#ifdef __APPLE__
#define CL_SILENCE_DEPRECATION
#endif
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#include <CL/cl2.hpp>

#include "util.hpp"

#include "err_code.h"

// pick up device type from compiler command line or from the default type
#ifndef DEVICE
#define DEVICE CL_DEVICE_TYPE_DEFAULT
#endif
// 


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

    // D3Q19 model
    std::vector<float> wt = setWt();
        
    std::vector<float> cx = setCx();
    std::vector<float> cy = setCy();
    std::vector<float> cz = setCz();
    
    std::vector<float> p(nx*ny*nz);
    std::vector<float> u(nx*ny*nz);
    std::vector<float> v(nx*ny*nz);
    std::vector<float> w(nx*ny*nz);

    std::vector<obstructure> obst(nx*ny*nz);

    const unsigned elements = nx*ny*nz;

    // // -- For cavity flow --
    const float u0 = 0.1f;
    const float rho0 = 1.0f;
    const float Re = 1000.0f;
    float nu = std::abs(u0)*float(nx-1)/Re;
    float dpdx = 0.0f;
    #include "boundaryCondition_cavityFlow.hpp"
    // // --

    //For channel flow
    // float Retau = 10;
    // float Retau = 1;
    // float utau = 0.005;
    // float nu = utau*0.5*ny/Retau;
    // float dpdx = utau*utau/(0.5*ny);
    // #include "boundaryCondition_channelFlow.hpp"

    // -- For Poiseuille flow --
    // float rho0 = 1.0;
    // float u0 = 0.005;
    // float h = 1.0;
    // const float Re = 100.0;
    // float nu = u0*h/Re;
    // float dpdx = u0/(h*h)*8.0*nu/(ny-1);
    // nu = nu*(ny-1);
    // #include "boundaryCondition_poiseuilleFlow.hpp"
    // --
     
    if(restart)
    {
        #include "restart.hpp"
    }
    else
    {
        #include "initialization.hpp"        
    }


    float Au = getAu(nu);
    unsigned nl = 1;
    float Ma = 0.1f;
    // float Cs = u0/Ma;
    // float Ap = getAp(Cs, nl);
    float Ap = 0.0f;
    float Cs = sqrt((nl/3.0)*(1.0f-Ap));
    Ma = u0/Cs;

    std::cout << "Au: " << Au << ", Ap: " << Ap << std::endl;    
    std::cout << "Ma: " << Ma << std::endl;



    std::chrono::system_clock::time_point start;
    std::chrono::system_clock::time_point end;
    start = std::chrono::system_clock::now();

    std::vector<unsigned> upID(19*elements);
    for(int ic = 0; ic < nx*ny*nz; ic++)
    {
        int i = ic2i(ic, nx, ny);
        int j = ic2j(ic, nx, ny);
        int k = ic2k(ic, nx, ny);

        for(int q = 0; q < 19; q++)
        {
            int qic = idf(q,ic,nx,ny,nz);
            upID[qic] = upwindID(q,i,j,k,nx,ny,nz);
        }
    }

    std::vector<unsigned> internalID;
    std::vector<unsigned> boundaryID;
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

    std::vector<unsigned> innerID;
    std::vector<unsigned> wallID;
    for(int bID = 0; bID < boundaryID.size(); bID++)
    {
        int ic = boundaryID[bID];
        int i = ic2i(ic, nx, ny);
        int j = ic2j(ic, nx, ny);
        int k = ic2k(ic, nx, ny);

        if(obst[ic].boundary == 1 || obst[ic].boundary == 2)
        {
            wallID.push_back(ic);
            innerID.push_back(index1d(i-obst[ic].normal[0],j-obst[ic].normal[1],k-obst[ic].normal[2],nx,ny));
        }
    }

    std::vector<unsigned> fixedWallID;
    std::vector<unsigned> movingWallID;
    for(int bID = 0; bID < boundaryID.size(); bID++)
    {
        int ic = boundaryID[bID];
        int i = ic2i(ic, nx, ny);
        int j = ic2j(ic, nx, ny);
        int k = ic2k(ic, nx, ny);

        if(obst[ic].boundary == 1)
        {
            fixedWallID.push_back(ic);
        }
        else if(obst[ic].boundary == 2)
        {
            movingWallID.push_back(ic);
        }
    }

    // Create a context
    cl::Context context(DEVICE);

    cl::Device device = context.getInfo<CL_CONTEXT_DEVICES>()[0];

    std::cout << std::endl << "Using OpenCL device: "
                << device.getInfo<CL_DEVICE_NAME>() << std::endl;
    

    // Load in kernel source, creating and building a program object for the context
    cl::Program program(context, util::loadProgram("kernel.cl"), true);

    
    cl::CommandQueue queue(context);

    cl::Buffer p_d(context, p.begin(), p.end(), true);
    cl::Buffer u_d(context, u.begin(), u.end(), true);
    cl::Buffer v_d(context, v.begin(), v.end(), true);
    cl::Buffer w_d(context, w.begin(), w.end(), true);

    cl::Buffer pTmp_d(context, p.begin(), p.end(), true);
    cl::Buffer uTmp_d(context, u.begin(), u.end(), true);
    cl::Buffer vTmp_d(context, v.begin(), v.end(), true);
    cl::Buffer wTmp_d(context, w.begin(), w.end(), true);

    cl::Buffer upID_d(context, upID.begin(), upID.end(), true);
    cl::Buffer internalID_d(context, internalID.begin(), internalID.end(), true);
    cl::Buffer boundaryID_d(context, boundaryID.begin(), boundaryID.end(), true);
    cl::Buffer innerID_d(context, innerID.begin(), innerID.end(), true);
    cl::Buffer wallID_d(context, wallID.begin(), wallID.end(), true);
    cl::Buffer fixedWallID_d(context, fixedWallID.begin(), fixedWallID.end(), true);
    cl::Buffer movingWallID_d(context, movingWallID.begin(), movingWallID.end(), true);

    // Create the kernel functor of updateP
    cl::KernelFunctor
    <
        cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer,
        cl::Buffer, 
        const unsigned,
        const float
    > updateP(program, "updateP");

    // Create the kernel functor of BC_P
    cl::KernelFunctor
    <
        cl::Buffer, 
        cl::Buffer, cl::Buffer
    > BC_P(program, "BC_P");

    // {
    // cl::Kernel ko_BC_P(program,"BC_P");
    // size_t work_group_size = ko_BC_P.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(cl::Device::getDefault());
    // std::cout << "work_group_size: " << work_group_size << std::endl;
    // }


    // Create the kernel functor of updateU
    cl::KernelFunctor
    <
        cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer, 
        const unsigned,
        const float
    > updateU(program, "updateU");

    // {
    // cl::Kernel ko_updateU(program,"updateU");
    // size_t work_group_size = ko_updateU.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(cl::Device::getDefault());
    // std::cout << "work_group_size: " << work_group_size << std::endl;
    // }

    // Create the kernel functor of BC_U_fixedWall
    cl::KernelFunctor
    <
        cl::Buffer, cl::Buffer, cl::Buffer, 
        cl::Buffer
    > BC_U_fixedWall(program, "BC_U_fixedWall");

    // Create the kernel functor of BC_U_movingWall
    cl::KernelFunctor
    <
        cl::Buffer, cl::Buffer, cl::Buffer, 
        cl::Buffer,
        const float
    > BC_U_movingWall(program, "BC_U_movingWall");

    // Time marching
    unsigned nt = startTimeStep;
    for(int nt = startTimeStep; nt <= endTimeStep; nt++)
    {   
        #include "write.hpp"

        {
        // util::Timer timer;
        updateP
        (
            cl::EnqueueArgs(queue,cl::NDRange(elements)),
            p_d, u_d, v_d, w_d,
            pTmp_d,
            upID_d,
            elements,
            Ap
        );
        queue.finish();
        // double rtime = static_cast<double>(timer.getTimeMicroseconds()) / 1000.0;
        // printf("\nThe kernel of updateP ran in %lf m seconds\n", rtime);
        }

        {
        // util::Timer timer;
        BC_P
        (
            cl::EnqueueArgs(queue,cl::NDRange(wallID.size())),
            pTmp_d,
            innerID_d, wallID_d
        );
        queue.finish();
        // double rtime = static_cast<double>(timer.getTimeMicroseconds()) / 1000.0;
        // printf("\nThe kernels of BC_P ran in %lf m seconds\n", rtime);
        }

        {
        // util::Timer timer;
        updateU
        (
            cl::EnqueueArgs(queue,cl::NDRange(elements)),
            pTmp_d, u_d, v_d, w_d,
            uTmp_d, vTmp_d, wTmp_d,
            upID_d,
            elements,
            Au
        );
        // queue.finish();
        // double rtime = static_cast<double>(timer.getTimeMicroseconds()) / 1000.0;
        // printf("\nThe kernel of updateU ran in %lf m seconds\n", rtime);
        }

        {
        // util::Timer timer;
        BC_U_fixedWall
        (
            cl::EnqueueArgs(queue,cl::NDRange(fixedWallID.size())),
            uTmp_d, vTmp_d, wTmp_d,
            fixedWallID_d
        );
        // queue.finish();
        // double rtime = static_cast<double>(timer.getTimeMicroseconds()) / 1000.0;
        // printf("\nThe kernels of BC_U_fixedWall ran in %lf m seconds\n", rtime);
        }

        {
        // util::Timer timer;
        BC_U_movingWall
        (
            cl::EnqueueArgs(queue,cl::NDRange(movingWallID.size())),
            uTmp_d, vTmp_d, wTmp_d,
            movingWallID_d,
            u0
        );
        queue.finish();
        // double rtime = static_cast<double>(timer.getTimeMicroseconds()) / 1000.0;
        // printf("\nThe kernels of BC_U_movingWall ran in %lf m seconds\n", rtime);
        }
        std::swap(pTmp_d,p_d);
        std::swap(uTmp_d,u_d);
        std::swap(vTmp_d,v_d);
        std::swap(wTmp_d,w_d);
    }

    end = std::chrono::system_clock::now();
    float time = static_cast<float>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() *1e-6);
    std::cout << "Execution time: " << time << " (s)" << std::endl;
    std::cout << "Speed: " << float(endTimeStep)*float(nx*ny*nz)/time*1e-6 << " (MLUPS)" << std::endl;

    return EXIT_SUCCESS;
}
