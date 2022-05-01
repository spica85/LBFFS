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
#include <arrayfire.h>
#include "LKS.hpp"
#include "input.hpp"

af::array mult(const af::array &lhs, const af::array &rhs)
{
    return lhs * rhs;
}

int main()
{
    af::info();

    // int backends = af::getAvailableBackends();
    // std::cout << backends << std::endl;
    // std::vector<af_backend> possible_backends;
    // if(backends & AF_BACKEND_CUDA)
    //     possible_backends.push_back(AF_BACKEND_CUDA);
    // if(backends & AF_BACKEND_OPENCL)
    //     possible_backends.push_back(AF_BACKEND_OPENCL);
    // if(backends & AF_BACKEND_CPU)
    //     possible_backends.push_back(AF_BACKEND_CPU);
    // if(possible_backends.size() == 0){
    //     std::cerr << "No backends available\n";
    //     return -1;
    // }
    // for(size_t i = 0; i < possible_backends.size(); ++i){
    //     af::setBackend (possible_backends[i]);
    //     af::setDevice (0);
    //     if(af::isDoubleAvailable (0)){
    //         std::cout << "Backend " << af::getActiveBackend () << " supports double\n";
    //         break;
    //     }
    // }
    // std::cout << "Backend " << af::getActiveBackend () << " is in use\n";
    

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
    const std::vector<float> wt = setWt();
    
    const std::vector<float> cx = setCx();
    const std::vector<float> cy = setCy();
    const std::vector<float> cz = setCz();

    const af::array wt_af = transpose(af::array(wt.size(), wt.data()));
    const af::array cx_af = transpose(af::array(cx.size(), cx.data()));
    const af::array cy_af = transpose(af::array(cy.size(), cy.data()));
    const af::array cz_af = transpose(af::array(cz.size(), cz.data()));

    
    std::vector<float> p(nx*ny*nz);
    std::vector<float> u(nx*ny*nz);
    std::vector<float> v(nx*ny*nz);
    std::vector<float> w(nx*ny*nz);

    std::vector<obstructure> obst(nx*ny*nz);

    const unsigned elements = nx*ny*nz;

    // // -- For cavity flow --
    const float u0 = 0.1f;
    const float rho0 = 1.0f;
    const float Re = 200.0f;
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

    af::array p_af(p.size(), p.data());
    af::array u_af(u.size(), u.data());
    af::array v_af(v.size(), v.data());
    af::array w_af(w.size(), w.data());

    float Au = getAu(nu);
    int nl = 1;
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
    const af::array upID_af(upID.size(), upID.data());

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
    const af::array internalID_af(internalID.size(), internalID.data());
    const af::array boundaryID_af(boundaryID.size(), boundaryID.data());

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
    const af::array innerID_af(innerID.size(), innerID.data());
    const af::array wallID_af(wallID.size(), wallID.data());

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
    const af::array fixedWallID_af(fixedWallID.size(), fixedWallID.data());
    const af::array movingWallID_af(movingWallID.size(), movingWallID.data());

    std::vector<float> uTmp(nx*ny*nz);
    std::vector<float> vTmp(nx*ny*nz);
    std::vector<float> wTmp(nx*ny*nz);
    std::vector<float> pTmp(nx*ny*nz);

    af::array pTmp_af = p_af;
    af::array uTmp_af = u_af;
    af::array vTmp_af = v_af;
    af::array wTmp_af = w_af;

    af::array pUp_af(elements*19);
    af::array uUp_af(elements*19);
    af::array vUp_af(elements*19);
    af::array wUp_af(elements*19);
    af::array uUpDotC(elements*19);
    af::array uSqr(elements*19);
    af::array f_eq_in(elements*19);
    af::array f_eq_in_u(elements*19);
    af::array uDotC(elements*19);
    af::array f(elements*19);
    af::array f_Delta_u(elements*19);
    af::array fcx(elements*19);
    af::array fcy(elements*19);
    af::array fcz(elements*19);

    // Time marching
    for(int nt = startTimeStep; nt <= endTimeStep; nt++)
    {   
        #include "write.hpp"

        // af::timer::start();

        uTmp_af = u_af;
        vTmp_af = v_af;
        wTmp_af = w_af;
        
        // af::sync();
        // auto stop_clock = af::timer::stop();
        // std::cout << "Tmp copy: " << stop_clock << std::endl;
        
        for(int l = 0; l < nl; l++)        
        {            
            

            pTmp_af = p_af;

            // af::timer::start();
            pUp_af = (tile(pTmp_af,19))(upID_af);
            // pUp_af = moddims(pUp_af,nx*ny*nz,19);
            // af::sync();
            // stop_clock = af::timer::stop();
            // std::cout << "Create pUp: " << stop_clock << std::endl;

            // af::timer::start();
            uUp_af = (tile(uTmp_af,19))(upID_af);
            // uUp_af = moddims(uUp_af,nx*ny*nz,19);
            // af::sync();
            // stop_clock = af::timer::stop();
            // std::cout << "Create uUp: " << stop_clock << std::endl;

            // af::timer::start();
            vUp_af = (tile(vTmp_af,19))(upID_af);
            // vUp_af = moddims(vUp_af,nx*ny*nz,19);
            // af::sync();
            // stop_clock = af::timer::stop();
            // std::cout << "Create vUp: " << stop_clock << std::endl;

            // af::timer::start();
            wUp_af = (tile(wTmp_af,19))(upID_af);
            // wUp_af = moddims(wUp_af,nx*ny*nz,19);
            // af::sync();
            // stop_clock = af::timer::stop();
            // std::cout << "Create wUp: " << stop_clock << std::endl;
            
            // af::timer::start();
            uUpDotC = flat(batchFunc(cx_af,moddims(uUp_af,elements,19),mult)) +flat(batchFunc(cy_af,moddims(vUp_af,elements,19),mult)) +flat(batchFunc(cz_af,moddims(wUp_af,elements,19),mult));
            // af::sync();
            // stop_clock = af::timer::stop();
            // std::cout << "Calculate uUpDotC: " << stop_clock << std::endl;
            
            // af::timer::start();
            uSqr = uUp_af*uUp_af +vUp_af*vUp_af +wUp_af*wUp_af;
            // af::sync();
            // stop_clock = af::timer::stop();
            // std::cout << "Calculate uSqr: " << stop_clock << std::endl;
            
            // af::timer::start();
            f_eq_in_u = flat(batchFunc(wt_af,moddims(3.0f*uUpDotC+4.5f*uUpDotC*uUpDotC-1.5*uSqr,elements,19),mult));
            f_eq_in = 3.0f*flat(batchFunc(wt_af,moddims(pUp_af,elements,19),mult)) +f_eq_in_u;
            // af::sync();
            // stop_clock = af::timer::stop();
            // std::cout << "Calculate f_eq_in: " << stop_clock << std::endl;
            
            // af::timer::start();
            uDotC = flat(batchFunc(cx_af,u_af,mult)) +flat(batchFunc(cy_af,v_af,mult)) +flat(batchFunc(cz_af,w_af,mult));
            // af::sync();
            // stop_clock = af::timer::stop();
            // std::cout << "Calculate uDotC: " << stop_clock << std::endl;
            
            // af::timer::start();
            f_Delta_u = flat(batchFunc(wt_af,moddims(uDotC-uUpDotC,elements,19),mult));
            f = f_eq_in +3.0f*Ap*f_Delta_u;
            // af::sync();
            // stop_clock = af::timer::stop();
            // std::cout << "Calculate f: " << stop_clock << std::endl;

            // af::timer::start();
            p_af = sum(moddims(f,elements,19),1)/3.0f;
            
            // af::sync();
            // stop_clock = af::timer::stop();
            // std::cout << "Update p: " << stop_clock << std::endl;

            // af::timer::start();
            p_af(wallID_af) = pTmp_af(innerID_af);
            // af::sync();
            // stop_clock = af::timer::stop();
            // std::cout << "Update BC of p: " << stop_clock << std::endl;

            // af::sync();
            // stop_clock = af::timer::stop();
            // std::cout << "Update p: " << stop_clock << std::endl;
        }
    
        // af::timer::start();

        pUp_af = (tile(p_af,19))(upID_af);
        // pUp_af = moddims(pUp_af,nx*ny*nz,19);
        
        f_eq_in = 3.0f*flat(batchFunc(wt_af,moddims(pUp_af,elements,19),mult)) + f_eq_in_u;

        f = f_eq_in +3.0f*Au*f_Delta_u;
        fcx = flat(batchFunc(cx_af, moddims(f,elements,19), mult));
        fcy = flat(batchFunc(cy_af, moddims(f,elements,19), mult));
        fcz = flat(batchFunc(cz_af, moddims(f,elements,19), mult));

        u_af = sum(moddims(fcx,elements,19),1);
        v_af = sum(moddims(fcy,elements,19),1);
        w_af = sum(moddims(fcz,elements,19),1);

        u_af(fixedWallID_af) = 0.0f;
        v_af(fixedWallID_af) = 0.0f;
        w_af(fixedWallID_af) = 0.0f;

        u_af(movingWallID_af) = u0;
        v_af(movingWallID_af) = 0.0f;
        w_af(movingWallID_af) = 0.0f;

        // af::sync();
        // stop_clock = af::timer::stop();
        // std::cout << "Update U: " << stop_clock << std::endl;
    }

    end = std::chrono::system_clock::now();
    float time = static_cast<float>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() *1e-6);
    std::cout << "Execution time: " << time << " (s)" << std::endl;
    std::cout << "Speed: " << float(endTimeStep)*float(nx*ny*nz)/time << " (LUPS)" << std::endl;

    return EXIT_SUCCESS;
}
