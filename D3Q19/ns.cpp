#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <omp.h>

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
    float u0;
    float v0;
    float w0;
    float rho0;
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

void streaming(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, std::vector<float>* pf, std::vector<float>& V)
{
    float& Vc    = V[0]; //(0)
    float& Vin   = V[1]; //(+x)
    float& Vip   = V[2]; //(-x)
    float& Vjn   = V[3]; //(+y)
    float& Vjp   = V[4]; //(-y)
    float& Vkn   = V[5]; //(+z)
    float& Vkp   = V[6]; //(-z)
    float& Vinjn = V[7]; //(+x,+y)
    float& Vipjp = V[8]; //(-x,-y)
    float& Vipjn = V[9]; //(-x,+y)
    float& Vinjp = V[10];//(+x,-y)
    float& Vinkn = V[11];//(+x,+z)
    float& Vipkp = V[12];//(-x,-z)
    float& Vjnkn = V[13];//(+y,+z)
    float& Vjpkp = V[14];//(-y,-z)
    float& Vipkn = V[15];//(-x,+z)
    float& Vinkp = V[16];//(+x,-z)
    float& Vjpkn = V[17];//(-y,+z)
    float& Vjnkp = V[18];//(+y,-z)

    int in = (i != 0)    ? index1d(i-1, j, k, nx, ny) : index1d(nx-1, j, k, nx, ny);
    int ip = (i != nx-1) ? index1d(i+1, j, k, nx, ny) : index1d(0, j, k, nx, ny);
    int jn = (j != 0 )   ? index1d(i, j-1, k, nx, ny) : index1d(i, ny-1, k, nx, ny);
    int jp = (j != ny-1) ? index1d(i, j+1, k, nx, ny) : index1d(i, 0, k, nx, ny);
    int kn = (k != 0 )   ? index1d(i, j, k-1, nx, ny) : index1d(i, j, nz-1, nx, ny);
    int kp = (k != nz-1) ? index1d(i, j, k+1, nx, ny) : index1d(i, j, 0, nx, ny);

    int injn = (i != 0)    ? jn -1  : jn +(nx-1);
    int injp = (i != 0)    ? jp -1  : jp +(nx-1);
    int inkn = (i != 0)    ? kn -1  : kn +(nx-1);
    int inkp = (i != 0)    ? kp -1  : kp +(nx-1);

    int ipjn = (i != nx-1) ? jn +1  : jn -(nx-1);
    int ipjp = (i != nx-1) ? jp +1  : jp -(nx-1);
    int ipkn = (i != nx-1) ? kn +1  : kn -(nx-1);
    int ipkp = (i != nx-1) ? kp +1  : kp -(nx-1);
    int jnkn = (j != 0)    ? kn -nx : kn +nx*(ny-1);
    int jnkp = (j != 0)    ? kp -nx : kp +nx*(ny-1);
    int jpkn = (j != ny-1) ? kn +nx : kn -nx*(ny-1);
    int jpkp = (j != ny-1) ? kp +nx : kp -nx*(ny-1);

    Vc    = (*pf)[idf(0, ic  ,nx,ny,nz)];
    Vin   = (*pf)[idf(1, in  ,nx,ny,nz)];
    Vip   = (*pf)[idf(2, ip  ,nx,ny,nz)];
    Vjn   = (*pf)[idf(3, jn  ,nx,ny,nz)];
    Vjp   = (*pf)[idf(4, jp  ,nx,ny,nz)];
    Vkn   = (*pf)[idf(5, kn  ,nx,ny,nz)];
    Vkp   = (*pf)[idf(6, kp  ,nx,ny,nz)];
    Vinjn = (*pf)[idf(7, injn,nx,ny,nz)];
    Vipjp = (*pf)[idf(8, ipjp,nx,ny,nz)];
    Vipjn = (*pf)[idf(9, ipjn,nx,ny,nz)];
    Vinjp = (*pf)[idf(10,injp,nx,ny,nz)];
    Vinkn = (*pf)[idf(11,inkn,nx,ny,nz)];
    Vipkp = (*pf)[idf(12,ipkp,nx,ny,nz)];
    Vjnkn = (*pf)[idf(13,jnkn,nx,ny,nz)];
    Vjpkp = (*pf)[idf(14,jpkp,nx,ny,nz)];
    Vipkn = (*pf)[idf(15,ipkn,nx,ny,nz)];
    Vinkp = (*pf)[idf(16,inkp,nx,ny,nz)];
    Vjpkn = (*pf)[idf(17,jpkn,nx,ny,nz)];
    Vjnkp = (*pf)[idf(18,jnkp,nx,ny,nz)];
}

void boundaryConditions(obstructure& obst, std::vector<float>& V)
{
    float& Vc    = V[0]; //(0)
    float& Vin   = V[1]; //(+x)
    float& Vip   = V[2]; //(-x)
    float& Vjn   = V[3]; //(+y)
    float& Vjp   = V[4]; //(-y)
    float& Vkn   = V[5]; //(+z)
    float& Vkp   = V[6]; //(-z)
    float& Vinjn = V[7]; //(+x,+y)
    float& Vipjp = V[8]; //(-x,-y)
    float& Vipjn = V[9]; //(-x,+y)
    float& Vinjp = V[10];//(+x,-y)
    float& Vinkn = V[11];//(+x,+z)
    float& Vipkp = V[12];//(-x,-z)
    float& Vjnkn = V[13];//(+y,+z)
    float& Vjpkp = V[14];//(-y,-z)
    float& Vipkn = V[15];//(-x,+z)
    float& Vinkp = V[16];//(+x,-z)
    float& Vjpkn = V[17];//(-y,+z)
    float& Vjnkp = V[18];//(+y,-z)

    if(obst.boundary == 1) //bounce back
    {
        if(obst.normal == -1)
        {
            Vin   = Vip;
            Vinjn = Vipjp;
            Vinjp = Vipjn;
            Vinkn = Vipkp;
            Vinkp = Vipkn;
        }
        else if(obst.normal == 1)
        {
            Vip   = Vin;
            Vipjp = Vinjn;
            Vipjn = Vinjp;
            Vipkp = Vinkn;
            Vipkn = Vinkp;
        }
        else if(obst.normal == -2)
        {
            Vjn  = Vjp;
            Vinjn = Vipjp;
            Vipjn = Vinjp;
            Vjnkn = Vjpkp;
            Vjnkp = Vjpkn;
        }
        else if(obst.normal == 2)
        {
            Vjp  = Vjn;
            Vinjp = Vipjn;
            Vipjp = Vinjn;
            Vjpkn = Vjnkp;
            Vjpkp = Vjnkn;
        }
        else if(obst.normal == -3)
        {
            Vkn   = Vkp;
            Vinkn = Vipkp;
            Vipkn = Vinkp;
            Vjnkn = Vjpkp;
            Vjpkn = Vjnkp;
        }
        else if(obst.normal == 3)
        {
            Vkp   = Vkn;
            Vinkp = Vipkn;
            Vipkp = Vinkn;
            Vjnkp = Vjpkn;
            Vjpkp = Vjnkn;
        }
    }
    else if(obst.boundary == 2) //wall velocity by bounce back
    {
        float u = obst.u0;
        float v = obst.v0;
        float w = obst.w0;

        if(obst.normal == -1)
        {
            float rho =(V[0]+V[3]+V[4]+V[5]+V[6]+V[13]+V[14]+V[17]+V[18]+2.0*(V[2]+V[8]+V[9]+V[12]+V[15]))/(1.0+u);
            V[1]  = V[2]  +rho*u/3.0;
            V[10] = V[9]  +rho*u/6.0 +0.5*(V[3]+V[13]+V[18] -(V[4]+V[14]+V[17])) -rho*v/2.0;
            V[7]  = V[8]  +rho*u/6.0 -0.5*(V[3]+V[13]+V[18] -(V[4]+V[14]+V[17])) +rho*v/2.0;
            V[11] = V[12] +rho*u/6.0 -0.5*(V[5]+V[9]+V[13]  -(V[6]+V[14]+V[18])) +rho*w/2.0;
            V[16] = V[15] +rho*u/6.0 +0.5*(V[5]+V[9]+V[13]  -(V[6]+V[14]+V[18])) -rho*w/2.0;
        }
        else if(obst.normal == 1)
        {
            float rho =(V[0]+V[3]+V[4]+V[5]+V[6]+V[13]+V[14]+V[17]+V[18]+2.0*(V[1]+V[7]+V[10]+V[11]+V[16]))/(1.0-u);
            V[2]  = V[1]  -rho*u/3.0;
            V[9]  = V[10] -rho*u/6.0 -0.5*(V[3]+V[13]+V[18] -(V[4]+V[14]+V[17])) +rho*v/2.0;
            V[8]  = V[7]  -rho*u/6.0 +0.5*(V[3]+V[13]+V[18] -(V[4]+V[14]+V[17])) -rho*v/2.0;
            V[12] = V[11] -rho*u/6.0 +0.5*(V[5]+V[9]+V[13]  -(V[6]+V[14]+V[18])) -rho*w/2.0;
            V[15] = V[16] -rho*u/6.0 -0.5*(V[5]+V[9]+V[13]  -(V[6]+V[14]+V[18])) +rho*w/2.0;
        }
        else if(obst.normal == -2)
        {
            float rho =(V[0]+V[1]+V[2]+V[5]+V[6]+V[11]+V[12]+V[15]+V[16] +2.0*(V[4]+V[8]+V[10]+V[14]+V[17]))/(1.0-v);
            V[3]  = V[4]  +rho*v/3.0;
            V[7]  = V[8]  +rho*v/6.0 -0.5*(V[1]+V[11]+V[16] -(V[2]+V[12]+V[15])) +rho*u/2.0;
            V[9]  = V[10] +rho*v/6.0 +0.5*(V[1]+V[11]+V[16] -(V[2]+V[12]+V[15])) -rho*u/2.0;
            V[13] = V[14] +rho*v/6.0 -0.5*(V[5]+V[11]+V[15] -(V[6]+V[12]+V[16])) +rho*w/2.0;
            V[18] = V[17] +rho*v/6.0 +0.5*(V[5]+V[11]+V[15] -(V[6]+V[12]+V[16])) -rho*w/2.0;
        }
        else if(obst.normal == 2)
        {
            float rho =(V[0]+V[1]+V[2]+V[5]+V[6]+V[11]+V[12]+V[15]+V[16] +2.0*(V[3]+V[7]+V[9]+V[13]+V[18]))/(1.0+v);
            V[4]  = V[3]  -rho*v/3.0;
            V[8]  = V[7]  -rho*v/6.0 +0.5*(V[1]+V[11]+V[16] -(V[2]+V[12]+V[15])) -rho*u/2.0;
            V[10] = V[9]  -rho*v/6.0 -0.5*(V[1]+V[11]+V[16] -(V[2]+V[12]+V[15])) +rho*u/2.0;
            V[14] = V[13] -rho*v/6.0 +0.5*(V[5]+V[11]+V[15] -(V[6]+V[12]+V[16])) -rho*w/2.0;
            V[17] = V[18] -rho*v/6.0 -0.5*(V[5]+V[11]+V[15] -(V[6]+V[12]+V[16])) +rho*w/2.0;
        }
        else if(obst.normal == -3)
        {
            float rho = (V[0]+V[1]+V[2]+V[3]+V[4]+V[7]+V[8]+V[9]+V[10]+2.0*(V[6]+V[12]+V[14]+V[16]+V[18]))/(1.0-w);
            V[5]  = V[6]  +rho*w/3.0;
            V[11] = V[12] +rho*w/6.0 -0.5*(V[1]+V[7]+V[10] -(V[2]+V[9]+V[8]))  +rho*u/2.0;
            V[15] = V[16] +rho*w/6.0 +0.5*(V[1]+V[7]+V[10] -(V[2]+V[9]+V[8]))  -rho*u/2.0;
            V[13] = V[14] +rho*w/6.0 -0.5*(V[3]+V[7]+V[9]  -(V[4]+V[8]+V[10])) +rho*v/2.0;
            V[17] = V[18] +rho*w/6.0 +0.5*(V[3]+V[7]+V[9]  -(V[4]+V[8]+V[10])) -rho*v/2.0;
        }
        else if(obst.normal == 3)
        {
            float rho = (V[0]+V[1]+V[2]+V[3]+V[4]+V[7]+V[8]+V[9]+V[10]+2.0*(V[5]+V[11]+V[13]+V[15]+V[17]))/(1.0+w);
            V[6]  = V[5]  -rho*w/3.0;
            V[16] = V[15] -rho*w/6.0 -0.5*(V[1]+V[7]+V[10] -(V[2]+V[8]+V[9]))  +rho*u/2.0;
            V[12] = V[11] -rho*w/6.0 +0.5*(V[1]+V[7]+V[10] -(V[2]+V[8]+V[9]))  -rho*u/2.0;
            V[18] = V[17] -rho*w/6.0 -0.5*(V[3]+V[7]+V[9]  -(V[4]+V[8]+V[10])) +rho*v/2.0;
            V[14] = V[13] -rho*w/6.0 +0.5*(V[3]+V[7]+V[9]  -(V[4]+V[8]+V[10])) -rho*v/2.0;
        }
    }
}

void collision(const float omega, const int ic, const int nx, const int ny, const int nz, const std::vector<float>& V, const std::vector<float>& cx, const std::vector<float>& cy, const std::vector<float>& cz, const std::vector<float>& wt, std::vector<float>* pftmp)
{
    float rho = 0.0;
    float u = 0.0;
    float v = 0.0;
    float w = 0.0;
    for(int q = 0; q < 19; q++)
    {
        rho += V[q];

        u += V[q]*cx[q];
        v += V[q]*cy[q];
        w += V[q]*cz[q];
    }
    u /= rho;
    v /= rho;
    w /= rho;

    for(int q = 0; q < 19; q++)
    {
        float uSqr =u*u+v*v+w*w;
        float uDotC = u*cx[q]+v*cy[q]+w*cz[q];
        float feq = (1.0+3.0*uDotC +4.5*uDotC*uDotC -1.5*uSqr)*wt[q]*rho;

        (*pftmp)[idf(q,ic,nx,ny,nz)] = (V[q] -omega *(V[q] -feq));
    }
}

int main()
{
    // omp_set_num_threads(4);
    
    bool restart = false;
    /* bool restart = true; */
    
    bool Fwrite = true;
    // bool Fwrite = false;
    bool writeBinary = true;
    /* bool writeBinary = false; */

    int startTimeStep = 0;
    const int endTimeStep = 1000;

    int nextOutTime = startTimeStep;
    int outInterval = 200;
    
    const int nx = 40; // number of x points
    const int ny = 40;// number of y points
    const int nz = 40;
    
    float u0 = 0.1;
    float rho0 = 1.0;
    float Re = 100.0;

    // Single Relaxation Time model
    float nu = std::abs(u0)*float(nx)/Re;
    std::cout << "nu = " << nu << "\n";
    float omega = 1.0/(3.0*nu +0.5);

    const std::vector<float> wt = {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

    // D3Q19 model
    const std::vector<float> cx = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0};
    const std::vector<float> cy = {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0, 1.0};
    const std::vector<float> cz = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
    

    std::vector<float> f(19*nx*ny*nz);
    std::vector<float> ftmp(19*nx*ny*nz);

    std::vector<float>* pf = &f;
    std::vector<float>* pftmp = &ftmp;
    std::vector<float>* ptmp;
    
    std::vector<float> rho(nx*ny*nz);
    std::vector<float> u(nx*ny*nz);
    std::vector<float> v(nx*ny*nz);
    std::vector<float> w(nx*ny*nz);

    std::vector<obstructure> obst(nx*ny*nz);
    

    // Setting conditions
    #include "boundaryCondition.hpp"
    if(restart)
    {
        #include "restart.hpp"
    }
    else
    {
        #include "initialization.hpp"        
    }

    // Time marching
    for(int nt = startTimeStep; nt <= endTimeStep; nt++)
    {
        #include "cal_rho_u_v.hpp"            
        #include "write.hpp"

        #pragma omp parallel for                
        for(int k = 0; k < nz; k++)
        {            
            for(int j = 0; j < ny; j++)                    
            {
                for(int i = 0; i < nx; i++)                        
                {
                    int ic = index1d(i,j,k,nx,ny);                    
                    std::vector<float> V(19);

                    streaming(ic,i,j,k,nx,ny,nz,pf,V);
                    boundaryConditions(obst[ic],V);
                    collision(omega,ic,nx,ny,nz,V,cx,cy,cz,wt,pftmp);
                }
            }
        }

        ptmp = pftmp;
        pftmp = pf;
        pf = ptmp;
    }

    return 0;
}
