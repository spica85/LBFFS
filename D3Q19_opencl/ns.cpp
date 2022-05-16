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
    const float U0 = 0.05; //Dimensionless maximum velocity


    // -- For cavity flow --
    // const float a = 1.0; //Dimensional length of system (m)
    // const float L = a/float(nx); //Representative length (-)
    // const float c = uMax/U0; //Representative velocity (m/s)
    

    // float nu = uMax*a/Re; //Dimensional kinematic viscosity
    // nu = nu/(L*c); //Dimensionless kinematic viscosity
    // float dpdx = 0.0;//Dimensionless external force


    //For channel flow
    // float Retau = 10;
    // float Retau = 1;
    // float utau = 0.005;
    // float nu = utau*0.5*ny/Retau;
    // float dpdx = utau*utau/(0.5*ny);

    //For Poiseuille flow
    float h = 1.f;
    float nu = uMax*h/Re;
    const float c = uMax/U0;
    const float L = h/float(ny);
    nu = nu/(L*c);
    // const float dpdx = uMax/(h*h)*8.0*nu/(ny-1);
    float dpdx = 8.0*(nu*(L*c))*uMax/(h*h);
    // float dpdx = 1.f-5;
    dpdx *= L/(c*c);

    std::cout << "dpdx = " << dpdx << std::endl;

    std::cout << "nu = " << nu << std::endl;


    const float deltaT = L/c; //Dimensional time step (s)
    const float omega = 1.0/(3.0*nu +0.5);
    // float omega = 1.0/0.56;

    std::cout << "tau = " << 1.0/omega << std::endl;

    const std::vector<float> wt = setWt();

    // D3Q19 model
    const std::vector<float> cx = setCx();
    const std::vector<float> cy = setCy();
    const std::vector<float> cz = setCz();

    std::vector<float> f(19*nx*ny*nz);
    std::vector<float> fTmp(19*nx*ny*nz);
    
    std::vector<float> rho(nx*ny*nz);
    std::vector<float> u(nx*ny*nz);
    std::vector<float> v(nx*ny*nz);
    std::vector<float> w(nx*ny*nz);
    float rho_av = 1.0;

    std::vector<obstructure> obst(nx*ny*nz);
    std::vector<float> u0(nx*ny*nz,0.0f);
    std::vector<float> v0(nx*ny*nz,0.0f);
    std::vector<float> w0(nx*ny*nz,0.0f);
    std::vector<int> normal(nx*ny*nz,0);

    const unsigned elements = nx*ny*nz;
    const unsigned qElements = 19*elements;

    // Setting conditions
    // #include "boundaryCondition_cavityFlow2d.hpp"
    // #include "boundaryCondition_cavityFlow3dDiagonal.hpp"
    #include "boundaryCondition_channelFlow.hpp"
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

        const std::vector<int> normal_in = {-1,0,0};
        const std::vector<int> normal_ip = {1,0,0};
        const std::vector<int> normal_jn = {0,-1,0};
        const std::vector<int> normal_jp = {0,1,0};
        const std::vector<int> normal_kn = {0,0,-1};
        const std::vector<int> normal_kp = {0,0,1};

        const std::vector<int> normal_injn = {-1,-1,0};
        const std::vector<int> normal_injp = {-1,1,0};
        const std::vector<int> normal_ipjn = {1,-1,0};
        const std::vector<int> normal_ipjp = {1,1,0};
        const std::vector<int> normal_inkn = {-1,0,-1};
        const std::vector<int> normal_inkp = {-1,0,1};
        const std::vector<int> normal_ipkn = {1,0,-1};
        const std::vector<int> normal_ipkp = {1,0,1};
        const std::vector<int> normal_jnkn = {0,-1,-1};
        const std::vector<int> normal_jnkp = {0,-1,1};
        const std::vector<int> normal_jpkn = {0,1,-1};
        const std::vector<int> normal_jpkp = {0,1,1};

        const std::vector<int> normal_injnkn = {-1,-1,-1};
        const std::vector<int> normal_injpkn = {-1,1,-1};
        const std::vector<int> normal_injnkp = {-1,-1,1};
        const std::vector<int> normal_injpkp = {-1,1,1};
        const std::vector<int> normal_ipjnkn = {1,-1,-1};
        const std::vector<int> normal_ipjpkn = {1,1,-1};
        const std::vector<int> normal_ipjnkp = {1,-1,1};
        const std::vector<int> normal_ipjpkp = {1,1,1};

        if(obst[ic].boundary > 0)
        {
            if(obst[ic].normal == normal_in)
            {
                normal[ic] = -1;
            }
            else if(obst[ic].normal == normal_ip)
            {
                normal[ic] = 1;
            }
            else if(obst[ic].normal == normal_jn)
            {
                normal[ic] = -2;
            }
            else if(obst[ic].normal == normal_jp)
            {
                normal[ic] = 2;
            }
            else if(obst[ic].normal == normal_kn)
            {
                normal[ic] = -3;
            }
            else if(obst[ic].normal == normal_kp)
            {
                normal[ic] = 3;
            }
            else if(obst[ic].normal == normal_injn)
            {
                normal[ic] = -4;
            }
            else if(obst[ic].normal == normal_ipjp)
            {
                normal[ic] = 4;
            }
            else if(obst[ic].normal == normal_ipjn)
            {
                normal[ic] = -5;
            }
            else if(obst[ic].normal == normal_injp)
            {
                normal[ic] = 5;
            }
            else if(obst[ic].normal == normal_inkn)
            {
                normal[ic] = -6;
            }
            else if(obst[ic].normal == normal_ipkp)
            {
                normal[ic] = 6;
            }
            else if(obst[ic].normal == normal_ipkn)
            {
                normal[ic] = -7;
            }
            else if(obst[ic].normal == normal_inkp)
            {
                normal[ic] = 7;
            }
            else if(obst[ic].normal == normal_jnkn)
            {
                normal[ic] = -8;
            }
            else if(obst[ic].normal == normal_jpkp)
            {
                normal[ic] = 8;
            }
            else if(obst[ic].normal == normal_jpkn)
            {
                normal[ic] = -9;
            }
            else if(obst[ic].normal == normal_jnkp)
            {
                normal[ic] = 9;
            }
            else if(obst[ic].normal == normal_injnkn)
            {
                normal[ic] = -10;
            }
            else if(obst[ic].normal == normal_ipjpkp)
            {
                normal[ic] = 10;
            }
            else if(obst[ic].normal == normal_ipjnkn)
            {
                normal[ic] = -11;
            }
            else if(obst[ic].normal == normal_injpkp)
            {
                normal[ic] = 11;
            }
            else if(obst[ic].normal == normal_injpkn)
            {
                normal[ic] = -12;
            }
            else if(obst[ic].normal == normal_ipjnkp)
            {
                normal[ic] = 12;
            }
            else if(obst[ic].normal == normal_injnkp)
            {
                normal[ic] = -13;
            }
            else if(obst[ic].normal == normal_ipjpkn)
            {
                normal[ic] = 13;
            }
            else
            {
                normal[ic] = 0;
            }
        }
        else
        {
            normal[ic] = 0;
        }
    }

    std::vector<int> upID(19*nx*ny*nz);
    std::vector<unsigned> upQID(19*nx*ny*nz);
    std::vector<unsigned> downQID(19*nx*ny*nz);
    for(int ic = 0; ic < nx*ny*nz; ic++)
    {
        int i = ic2i(ic,nx,ny);
        int j = ic2j(ic,nx,ny);
        int k = ic2k(ic,nx,ny);
        for(int q = 0; q < 19; q++)
        {
            int qic = idf(q,ic,nx,ny,nz);
            upID[qic] = upwindID(q,i,j,k,nx,ny,nz);
            upQID[qic] = idf(q, upID[qic], nx, ny, nz);
            int downID = downwindID(q,i,j,k,nx,ny,nz);
            downQID[qic] = idf(q, downID, nx, ny, nz);
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

    // std::vector<int> BBID;
    // std::vector<int> BBmovingWallID;
    // for(int bID = 0; bID < boundaryID.size(); bID++)
    // {
    //     int ic = boundaryID[bID];
    //     if(obst[ic].boundary == 1)
    //     {
    //         BBID.push_back(ic);
    //     }
    //     else if(obst[ic].boundary == 2)
    //     {
    //         BBmovingWallID.push_back(ic);
    //     }
    // }
    // const unsigned nBB = BBID.size();
    // const unsigned nBBMW = BBmovingWallID.size();

    // std::vector<int> BBQID(19*BBID.size());
    // for(int bID = 0; bID < BBID.size(); bID++)
    // {
    //     int ic = BBID[bID];

    //     int i = ic2i(ic,nx,ny);
    //     int j = ic2j(ic,nx,ny);
    //     int k = ic2k(ic,nx,ny);
    //     for(int q = 0; q < 19; q++)
    //     {
    //         int qic = q*BBID.size()+bID;
    //         BBQID[qic] = bounceBackQID(obst[ic],q,i,j,k,nx,ny,nz);
    //     }
    // }


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
    // cl::Buffer u_d(context, u.begin(), u.end(), true);
    // cl::Buffer v_d(context, v.begin(), v.end(), true);
    // cl::Buffer w_d(context, w.begin(), w.end(), true);

    cl::Buffer fTmp_d(context, fTmp.begin(), fTmp.end(), true);
    // cl::Buffer uTmp_d(context, u.begin(), u.end(), true);
    // cl::Buffer vTmp_d(context, v.begin(), v.end(), true);
    // cl::Buffer wTmp_d(context, w.begin(), w.end(), true);

    cl::Buffer upQID_d(context, upQID.begin(), upQID.end(), true);
    cl::Buffer downQID_d(context, downQID.begin(), downQID.end(), true);
    // cl::Buffer BBID_d(context, BBID.begin(), BBID.end(), true);
    // cl::Buffer BBQID_d(context, BBQID.begin(), BBQID.end(), true);
    // cl::Buffer BBmovingWallID_d(context, BBmovingWallID.begin(), BBmovingWallID.end(), true);
    cl::Buffer u0_d(context, u0.begin(), u0.end(), true);
    cl::Buffer v0_d(context, v0.begin(), v0.end(), true);
    cl::Buffer w0_d(context, w0.begin(), w0.end(), true);
    cl::Buffer normal_d(context, normal.begin(), normal.end(), true);
    cl::Buffer rho_d(context, rho.begin(), rho.end(), true);
    

    // Create the kernel functor of collision
    cl::KernelFunctor
    <
        cl::Buffer, cl::Buffer,
        const unsigned,
        const float
    > k_collision(program, "k_collision");

    // Create the kernel functor of externalForce
    cl::KernelFunctor
    <
        cl::Buffer, cl::Buffer,
        const unsigned,
        const float
    > k_externalForce(program, "k_externalForce");

    // Create the kernel functor of streaming
    cl::KernelFunctor
    <
        cl::Buffer, cl::Buffer,
        cl::Buffer
    > k_streaming(program, "k_streaming");

    // Create the kernel functor of collisionStreaming
    cl::KernelFunctor
    <
        cl::Buffer, cl::Buffer,
        const unsigned,
        const float,
        const float,
        const int, const int, const int
    > k_collisionStreaming(program, "k_collisionStreaming");

    // Create the kernel functor of bounceBack
    cl::KernelFunctor
    <
        cl::Buffer,
        cl::Buffer, cl::Buffer,
        const unsigned, const unsigned
    > k_bounceBack(program, "k_bounceBack");

    // Create the kernel functor of bounceBackMovingWall
    cl::KernelFunctor
    <
        cl::Buffer,
        cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer,
        const unsigned
    > k_bounceBackMovingWall(program, "k_bounceBackMovingWall");

    // Create the kernel functor of streamingCollision
    cl::KernelFunctor
    <
        cl::Buffer, cl::Buffer,
        cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer,
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
            normal_d,
            u0_d, v0_d, w0_d,
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
