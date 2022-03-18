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

inline int ic2i(int ic, int nx, int ny)
{
    return int((ic%(nx*ny))%nx);
}

inline int ic2j(int ic, int nx, int ny)
{
    return int(ic%(nx*ny)/nx);
}

inline int ic2k(int ic, int nx, int ny)
{
    return int(ic/(nx*ny));
}

class obstructure
{
    private:
    public:
    int boundary;
    std::vector<int> normal;
    bool inner;
    double u0;
    double v0;
    double w0;
    double rho0;
    obstructure()
    {
        boundary = 0;
        normal = {0, 0, 0};
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

inline int upwindID(const int q, const int i, const int j, const int k, const int nx, const int ny, const int nz)
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

void streaming(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, std::vector<double>& ftmp, std::vector<double>& f)
{
    for(int q = 0; q < 19; q++)
    {
        int qic = idf(q,ic,nx,ny,nz);
        int upID = upwindID(q,i,j,k,nx,ny,nz);
        f[qic] = ftmp[idf(q,upID,nx,ny,nz)];
    }
}

void streaming(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, std::vector<double>& ftmp, std::vector<double>& f, const std::vector<int>& upID)
{
    for(int q = 0; q < 19; q++)
    {
        int qic = idf(q,ic,nx,ny,nz);
        f[qic] = ftmp[idf(q,upID[qic],nx,ny,nz)];
    }
}

void boundaryConditions(obstructure& obst, std::vector<double>& f, int ic, int nx, int ny, int nz)
{
    double& Vc    = f[idf(0, ic,nx,ny,nz)];//(0)
    double& Vin   = f[idf(1, ic,nx,ny,nz)];//(+x)
    double& Vip   = f[idf(2, ic,nx,ny,nz)];//(-x)
    double& Vjn   = f[idf(3, ic,nx,ny,nz)];//(+y)
    double& Vjp   = f[idf(4, ic,nx,ny,nz)];//(-y)
    double& Vkn   = f[idf(5, ic,nx,ny,nz)];//(+z)
    double& Vkp   = f[idf(6, ic,nx,ny,nz)];//(-z)
    double& Vinjn = f[idf(7, ic,nx,ny,nz)];//(+x,+y)
    double& Vipjp = f[idf(8, ic,nx,ny,nz)];//(-x,-y)
    double& Vinjp = f[idf(9, ic,nx,ny,nz)];//(+x,-y)
    double& Vipjn = f[idf(10,ic,nx,ny,nz)];//(-x,+y)
    double& Vinkn = f[idf(11,ic,nx,ny,nz)];//(+x,+z)
    double& Vipkp = f[idf(12,ic,nx,ny,nz)];//(-x,-z)
    double& Vinkp = f[idf(13,ic,nx,ny,nz)];//(+x,-z)
    double& Vipkn = f[idf(14,ic,nx,ny,nz)];//(-x,+z)
    double& Vjnkn = f[idf(15,ic,nx,ny,nz)];//(+y,+z)
    double& Vjpkp = f[idf(16,ic,nx,ny,nz)];//(-y,-z)
    double& Vjnkp = f[idf(17,ic,nx,ny,nz)];//(+y,-z)
    double& Vjpkn = f[idf(18,ic,nx,ny,nz)];//(-y,+z)

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

    if(obst.boundary == 1) //bounce back
    {
        if(obst.normal == normal_in) //i=0
        {
            Vin   = Vip;
            Vinjn = Vipjp;
            Vinjp = Vipjn;
            Vinkn = Vipkp;
            Vinkp = Vipkn;
        }
        else if(obst.normal == normal_ip) //i=nx-1
        {
            Vip   = Vin;
            Vipjp = Vinjn;
            Vipjn = Vinjp;
            Vipkp = Vinkn;
            Vipkn = Vinkp;
        }
        else if(obst.normal == normal_jn) //j=0
        {
            Vjn  = Vjp;
            Vinjn = Vipjp;
            Vipjn = Vinjp;
            Vjnkn = Vjpkp;
            Vjnkp = Vjpkn;
        }
        else if(obst.normal == normal_jp) //j=ny-1
        {
            Vjp  = Vjn;
            Vinjp = Vipjn;
            Vipjp = Vinjn;
            Vjpkn = Vjnkp;
            Vjpkp = Vjnkn;
        }
        else if(obst.normal == normal_kn) //k=0
        {
            Vkn   = Vkp;
            Vinkn = Vipkp;
            Vipkn = Vinkp;
            Vjnkn = Vjpkp;
            Vjpkn = Vjnkp;
        }
        else if(obst.normal == normal_kp) //k=nz-1
        {
            Vkp   = Vkn;
            Vinkp = Vipkn;
            Vipkp = Vinkn;
            Vjnkp = Vjpkn;
            Vjpkp = Vjnkn;
        }
        else if(obst.normal == normal_injn) //i=0 && j=0
        {
            Vin   = Vip;
            Vinjn = Vipjp;
            // Vinjp = Vipjn;
            Vinkn = Vipkp;
            Vinkp = Vipkn;

            Vjn  = Vjp;
            Vinjn = Vipjp;
            // Vipjn = Vinjp;
            Vjnkn = Vjpkp;
            Vjnkp = Vjpkn;

            Vinjp = 1.0/36.0;
            Vipjn = 1.0/36.0;
        }
        else if(obst.normal == normal_injp) //i=0 && j=ny-1
        {
            Vin   = Vip;
            // Vinjn = Vipjp;
            Vinjp = Vipjn;
            Vinkn = Vipkp;
            Vinkp = Vipkn;

            Vjp  = Vjn;
            Vinjp = Vipjn;
            // Vipjp = Vinjn;
            Vjpkn = Vjnkp;
            Vjpkp = Vjnkn;

            Vinjn = 1.0/36.0;
            Vipjp = 1.0/36.0;
        }
        else if(obst.normal == normal_ipjn) //i=nx-1 && j=0
        {
            Vip   = Vin;
            // Vipjp = Vinjn;
            Vipjn = Vinjp;
            Vipkp = Vinkn;
            Vipkn = Vinkp;

            Vjn  = Vjp;
            // Vinjn = Vipjp;
            Vipjn = Vinjp;
            Vjnkn = Vjpkp;
            Vjnkp = Vjpkn;

            Vipjp = 1.0/36.0;
            Vinjn = 1.0/36.0;
        }
        else if(obst.normal == normal_ipjp) //i=nx-1 && j=ny-1
        {
            Vip   = Vin;
            Vipjp = Vinjn;
            // Vipjn = Vinjp;
            Vipkp = Vinkn;
            Vipkn = Vinkp;

            Vjp  = Vjn;
            // Vinjp = Vipjn;
            Vipjp = Vinjn;
            Vjpkn = Vjnkp;
            Vjpkp = Vjnkn;

            Vipjn = 1.0/36.0;
            Vinjp = 1.0/36.0;
        }
        else if(obst.normal == normal_inkn) //i=0 && k=0
        {
            Vin   = Vip;
            Vinjn = Vipjp;
            Vinjp = Vipjn;
            Vinkn = Vipkp;
            // Vinkp = Vipkn;

            Vkn   = Vkp;
            Vinkn = Vipkp;
            // Vipkn = Vinkp;
            Vjnkn = Vjpkp;
            Vjpkn = Vjnkp;

            Vinkp = 1.0/36.0;
            Vipkn = 1.0/36.0;
        }
        else if(obst.normal == normal_inkp) //i=0 && k=nz-1
        {
            Vin   = Vip;
            Vinjn = Vipjp;
            Vinjp = Vipjn;
            // Vinkn = Vipkp;
            Vinkp = Vipkn;

            Vkp   = Vkn;
            Vinkp = Vipkn;
            // Vipkp = Vinkn;
            Vjnkp = Vjpkn;
            Vjpkp = Vjnkn;

            Vinkn = 1.0/36.0;
            Vipkp = 1.0/36.0;
        }
        else if(obst.normal == normal_ipkn) //i=nx-1 && k=0
        {
            Vip   = Vin;
            Vipjp = Vinjn;
            Vipjn = Vinjp;
            // Vipkp = Vinkn;
            Vipkn = Vinkp;

            Vkn   = Vkp;
            // Vinkn = Vipkp;
            Vipkn = Vinkp;
            Vjnkn = Vjpkp;
            Vjpkn = Vjnkp;

            Vipkp = 1.0/36.0;
            Vinkn = 1.0/36.0;
        }
        else if(obst.normal == normal_ipkp) //i=nx-1 && k=nz-1
        {
            Vip   = Vin;
            Vipjp = Vinjn;
            Vipjn = Vinjp;
            Vipkp = Vinkn;
            // Vipkn = Vinkp;

            Vkp   = Vkn;
            // Vinkp = Vipkn;
            Vipkp = Vinkn;
            Vjnkp = Vjpkn;
            Vjpkp = Vjnkn;

            Vipkn = 1.0/36.0;
            Vinkp = 1.0/36.0;
        }
        else if(obst.normal == normal_jnkn) //j=0 && k=0
        {
            Vjn  = Vjp;
            Vinjn = Vipjp;
            Vipjn = Vinjp;
            Vjnkn = Vjpkp;
            // Vjnkp = Vjpkn;

            Vkn   = Vkp;
            Vinkn = Vipkp;
            Vipkn = Vinkp;
            Vjnkn = Vjpkp;
            // Vjpkn = Vjnkp;

            Vjnkp = 1.0/36.0;
            Vjpkn = 1.0/36.0;
        }
        else if(obst.normal == normal_jnkp) //j=0 && k=nz-1
        {
            Vjn  = Vjp;
            Vinjn = Vipjp;
            Vipjn = Vinjp;
            // Vjnkn = Vjpkp;
            Vjnkp = Vjpkn;

            Vkp   = Vkn;
            Vinkp = Vipkn;
            Vipkp = Vinkn;
            Vjnkp = Vjpkn;
            // Vjpkp = Vjnkn;

            Vjnkn = 1.0/36.0;
            Vjpkp = 1.0/36.0;
        }
        else if(obst.normal == normal_jpkn) //j=ny-1 && k=0
        {
            Vjp  = Vjn;
            Vinjp = Vipjn;
            Vipjp = Vinjn;
            Vjpkn = Vjnkp;
            // Vjpkp = Vjnkn;

            Vkn   = Vkp;
            Vinkn = Vipkp;
            Vipkn = Vinkp;
            // Vjnkn = Vjpkp;
            Vjpkn = Vjnkp;

            Vjpkp = 1.0/36.0;
            Vjnkn = 1.0/36.0;
        }
        else if(obst.normal == normal_jpkp) //j=ny-1 && k=nz-1
        {
            Vjp  = Vjn;
            Vinjp = Vipjn;
            Vipjp = Vinjn;
            // Vjpkn = Vjnkp;
            Vjpkp = Vjnkn;

            Vkp   = Vkn;
            Vinkp = Vipkn;
            Vipkp = Vinkn;
            // Vjnkp = Vjpkn;
            Vjpkp = Vjnkn;

            Vjpkn = 1.0/36.0;
            Vjnkp = 1.0/36.0;
        }
        else if(obst.normal == normal_injnkn) //i=0 && j=0 && k=0
        {
            Vin   = Vip;
            Vinjn = Vipjp;
            Vinjp = Vipjn;
            Vinkn = Vipkp;
            Vinkp = Vipkn;

            Vjn  = Vjp;
            Vinjn = Vipjp;
            Vipjn = Vinjp;
            Vjnkn = Vjpkp;
            Vjnkp = Vjpkn;

            Vkn   = Vkp;
            Vinkn = Vipkp;
            Vipkn = Vinkp;
            Vjnkn = Vjpkp;
            Vjpkn = Vjnkp;
        }
        else if(obst.normal == normal_injpkn) //i=0 && j=ny-1 && k=0
        {
            Vin   = Vip;
            // Vinjn = Vipjp;
            Vinjp = Vipjn;
            Vinkn = Vipkp;
            // Vinkp = Vipkn;

            Vjp  = Vjn;
            Vinjp = Vipjn;
            // Vipjp = Vinjn;
            Vjpkn = Vjnkp;
            Vjpkp = Vjnkn;

            Vkn   = Vkp;
            Vinkn = Vipkp;
            // Vipkn = Vinkp;
            Vjnkn = Vjpkp;
            Vjpkn = Vjnkp;

            Vinjn = 1.0/36.0;
            Vipjp = 1.0/36.0;
            Vinkp = 1.0/36.0;
            Vipkn = 1.0/36.0;
        }
        else if(obst.normal == normal_injnkp) //i=0 && j=0 && k=nz-1
        {
            Vin   = Vip;
            Vinjn = Vipjp;
            // Vinjp = Vipjn;
            // Vinkn = Vipkp;
            Vinkp = Vipkn;

            Vjn  = Vjp;
            Vinjn = Vipjp;
            // Vipjn = Vinjp;
            Vjnkn = Vjpkp;
            Vjnkp = Vjpkn;

            Vkp   = Vkn;
            Vinkp = Vipkn;
            // Vipkp = Vinkn;
            Vjnkp = Vjpkn;
            Vjpkp = Vjnkn;

            Vinjp = 1.0/36.0;
            Vipjn = 1.0/36.0;
            Vinkn = 1.0/36.0;
            Vipkp = 1.0/36.0;
        }
        else if(obst.normal == normal_injnkp) //i=0 && j=ny-1 && k=nz-1
        {
            Vin   = Vip;
            // Vinjn = Vipjp;
            Vinjp = Vipjn;
            // Vinkn = Vipkp;
            Vinkp = Vipkn;

            Vjp  = Vjn;
            Vinjp = Vipjn;
            // Vipjp = Vinjn;
            Vjpkn = Vjnkp;
            Vjpkp = Vjnkn;

            Vkp   = Vkn;
            Vinkp = Vipkn;
            // Vipkp = Vinkn;
            Vjnkp = Vjpkn;
            Vjpkp = Vjnkn;

            Vinjn = 1.0/36.0;
            Vipjp = 1.0/36.0;
            Vinkn = 1.0/36.0;
            Vipkp = 1.0/36.0;
        }
        else if(obst.normal == normal_ipjnkn) //i=nx-1 && j=0 && k=0
        {
            Vip   = Vin;
            // Vipjp = Vinjn;
            Vipjn = Vinjp;
            // Vipkp = Vinkn;
            Vipkn = Vinkp;

            Vjn  = Vjp;
            // Vinjn = Vipjp;
            Vipjn = Vinjp;
            Vjnkn = Vjpkp;
            Vjnkp = Vjpkn;

            Vkn   = Vkp;
            // Vinkn = Vipkp;
            Vipkn = Vinkp;
            Vjnkn = Vjpkp;
            Vjpkn = Vjnkp;

            Vipjp = 1.0/36.0;
            Vinjn = 1.0/36.0;
            Vipkp = 1.0/36.0;
            Vinkn = 1.0/36.0;
        }
        else if(obst.normal == normal_ipjpkn) //i=nx-1 && j=ny-1 && k=0
        {
            Vip   = Vin;
            Vipjp = Vinjn;
            // Vipjn = Vinjp;
            // Vipkp = Vinkn;
            Vipkn = Vinkp;

            Vjp  = Vjn;
            // Vinjp = Vipjn;
            Vipjp = Vinjn;
            Vjpkn = Vjnkp;
            Vjpkp = Vjnkn;

            Vkn   = Vkp;
            // Vinkn = Vipkp;
            Vipkn = Vinkp;
            Vjnkn = Vjpkp;
            Vjpkn = Vjnkp;

            Vipjn = 1.0/36.0;
            Vinjp = 1.0/36.0;
            Vipkp = 1.0/36.0;
            Vinkn = 1.0/36.0;
        }
        else if(obst.normal == normal_ipjpkn) //i=nx-1 && j=0 && k=nz-1
        {
            Vip   = Vin;
            // Vipjp = Vinjn;
            Vipjn = Vinjp;
            Vipkp = Vinkn;
            // Vipkn = Vinkp;

            Vjn  = Vjp;
            // Vinjn = Vipjp;
            Vipjn = Vinjp;
            Vjnkn = Vjpkp;
            Vjnkp = Vjpkn;

            Vkp   = Vkn;
            // Vinkp = Vipkn;
            Vipkp = Vinkn;
            Vjnkp = Vjpkn;
            Vjpkp = Vjnkn;

            Vipjp = 1.0/36.0;
            Vinjn = 1.0/36.0;
            Vipkn = 1.0/36.0;
            Vinkp = 1.0/36.0;
        }
        else if(obst.normal == normal_ipjpkp) //i=nx-1 && j=ny-1 && k=nz-1
        {
            Vip   = Vin;
            Vipjp = Vinjn;
            // Vipjn = Vinjp;
            Vipkp = Vinkn;
            // Vipkn = Vinkp;

            Vjp  = Vjn;
            // Vinjp = Vipjn;
            Vipjp = Vinjn;
            Vjpkn = Vjnkp;
            Vjpkp = Vjnkn;

            Vkp   = Vkn;
            // Vinkp = Vipkn;
            Vipkp = Vinkn;
            Vjnkp = Vjpkn;
            Vjpkp = Vjnkn;

            Vipjn = 1.0/36.0;
            Vinjp = 1.0/36.0;
            Vipkn = 1.0/36.0;
            Vinkp = 1.0/36.0;
        }
    }
    else if(obst.boundary == 2) //wall velocity by bounce back
    {
        double u = obst.u0;
        double v = obst.v0;
        double w = obst.w0;

        if(obst.normal == normal_in)
        {
            double rho =(Vc+Vjn+Vjp+Vkn+Vkp+Vjnkn+Vjpkp+Vjnkp+Vjpkn +2.0*(Vip+Vipjn+Vipjp+Vipkn+Vipkp))/(1.0-u);
            Vin = Vip +rho*u/3.0;
            double Nyx = -rho*v/3.0 +0.5*(Vjn+Vjnkn+Vjnkp-(Vjp+Vjpkn+Vjpkp));
            double Nzx = -rho*w/3.0 +0.5*(Vkn+Vjnkn+Vjpkn-(Vkp+Vjnkp+Vjpkp));
            Vinjn = Vipjp +rho*(u+v)/6.0 -Nyx;
            Vinjp = Vipjn +rho*(u-v)/6.0 +Nyx;
            Vinkn = Vipkp +rho*(u+w)/6.0 -Nzx;
            Vinkp = Vipkn +rho*(u-w)/6.0 +Nzx;
        }
        else if(obst.normal == normal_ip)
        {
            double rho =(Vc+Vjn+Vjp+Vkn+Vkp+Vjnkn+Vjpkp+Vjnkp+Vjpkn +2.0*(Vin+Vinjn+Vinjp+Vinkn+Vinkp))/(1.0+u);
            Vip = Vin -rho*u/3.0;
            double Nyx = -rho*v/3.0 +0.5*(Vjn+Vjnkn+Vjnkp-(Vjp+Vjpkn+Vjpkp));
            double Nzx = -rho*w/3.0 +0.5*(Vkn+Vjnkn+Vjpkn-(Vkp+Vjnkp+Vjpkp));
            Vipjp = Vinjn -rho*(u+v)/6.0 +Nyx;
            Vipjn = Vinjp -rho*(u-v)/6.0 -Nyx;
            Vipkp = Vinkn -rho*(u+w)/6.0 +Nzx;
            Vipkn = Vinkp -rho*(u-w)/6.0 -Nzx;
        }
        else if(obst.normal == normal_jn)
        {
            double rho = (Vc+Vin+Vip+Vkn+Vkp+Vinkn+Vipkp+Vipkn+Vinkp +2.0*(Vjp+Vinjp+Vipjp+Vjpkn+Vjpkp))/(1.0-v);
            Vjn = Vjp +rho*v/3.0;
            double Nxy = -rho*u/3.0 +0.5*(Vin+Vinkn+Vinkp-(Vip+Vipkn+Vipkp));
            double Nzy = -rho*w/3.0 +0.5*(Vkn+Vinkn+Vipkn-(Vkp+Vinkp+Vipkp));
            Vinjn = Vipjp +rho*(v+u)/6.0 -Nxy;
            Vipjn = Vinjp +rho*(v-u)/6.0 +Nxy;
            Vjnkn = Vjpkp +rho*(v+w)/6.0 -Nzy;
            Vjnkp = Vjpkn +rho*(v-w)/6.0 +Nzy;
        }
        else if(obst.normal == normal_jp)
        {
            double rho =(Vc+Vin+Vip+Vkn+Vkp+Vinkn+Vipkp+Vipkn+Vinkp +2.0*(Vjn+Vinjn+Vipjn+Vjnkn+Vjnkp))/(1.0+v);
            Vjp = Vjn -rho*v/3.0;
            double Nxy = -rho*u/3.0 +0.5*(Vin+Vinkn+Vinkp-(Vip+Vipkn+Vipkp));
            double Nzy = -rho*w/3.0 +0.5*(Vkn+Vinkn+Vipkn-(Vkp+Vinkp+Vipkp));
            Vipjp = Vinjn  -rho*(v+u)/6.0 +Nxy;
            Vinjp = Vipjn  -rho*(v-u)/6.0 -Nxy;
            Vjpkp = Vjnkn  -rho*(v+w)/6.0 +Nzy;
            Vjpkn = Vjnkp  -rho*(v-w)/6.0 -Nzy;
        }
        else if(obst.normal == normal_kn)
        {
            double rho = (Vc+Vin+Vip+Vjn+Vjp+Vinjn+Vipjp+Vinjp+Vipjn+2.0*(Vkp+Vipkp+Vinkp+Vjpkp+Vjnkp))/(1.0-w);
            Vkn = Vkp +rho*w/3.0;
            double Nxz = -rho*u/3.0 +0.5*(Vin+Vinjn+Vinjp-(Vip+Vipjn+Vipjp));
            double Nyz = -rho*v/3.0 +0.5*(Vjn+Vinjn+Vipjp-(Vjp+Vinjp+Vipjp));
            Vinkn = Vipkp +rho*(w+u)/6.0 -Nxz;
            Vipkn = Vinkp +rho*(w-u)/6.0 +Nxz;
            Vjnkn = Vjpkp +rho*(w+v)/6.0 -Nyz;
            Vjpkn = Vjnkn +rho*(w-v)/6.0 +Nyz;
        }
        else if(obst.normal == normal_kp)
        {
            double rho = (Vc+Vin+Vip+Vjn+Vjp+Vinjn+Vipjp+Vinjp+Vipjn+2.0*(Vkn+Vipkn+Vinkn+Vjpkn+Vjnkn))/(1.0+w);
            Vkp = Vkn -rho*w/3.0;
            double Nxz = -rho*u/3.0 +0.5*(Vin+Vinjn+Vinjp-(Vip+Vipjn+Vipjp));
            double Nyz = -rho*v/3.0 +0.5*(Vjn+Vinjn+Vipjp-(Vjp+Vinjp+Vipjp));
            Vipkp = Vinkn -rho*(w+u)/6.0 +Nxz;
            Vinkp = Vipkn -rho*(w-u)/6.0 -Nxz;
            Vjpkp = Vjnkn -rho*(w+v)/6.0 +Nyz;
            Vjnkp = Vjpkn -rho*(w-v)/6.0 -Nyz;
        }
    }
}

void collision(const double omega, const int ic, const int nx, const int ny, const int nz, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, std::vector<double>& f)
{
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double w = 0.0;
    for(int q = 0; q < 19; q++)
    {
        int qic = idf(q,ic,nx,ny,nz);

        rho += f[qic];

        u += f[qic]*cx[q];
        v += f[qic]*cy[q];
        w += f[qic]*cz[q];
    }
    u /= rho;
    v /= rho;
    w /= rho;

    for(int q = 0; q < 19; q++)
    {
        double uSqr =u*u+v*v+w*w;
        double uDotC = u*cx[q]+v*cy[q]+w*cz[q];
        double feq = (1.0+3.0*uDotC +4.5*uDotC*uDotC -1.5*uSqr)*wt[q]*rho;

        int qic = idf(q,ic,nx,ny,nz);

        f[qic] += -omega *(f[qic] -feq);
    }
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

void cal_rhoUVW(int ic, int nx, int ny, int nz, const std::vector<double>& f, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, double& rho, double& u, double& v, double& w)
{
    rho = 0.0;
    u = 0.0;
    v = 0.0;
    w = 0.0;

    for(int q = 0; q < 19; q++)
    {
        int qic = idf(q,ic,nx,ny,nz);
        rho += f[qic];
        u += f[qic]*cx[q];
        v += f[qic]*cy[q];
        w += f[qic]*cz[q];
    }
    u /= rho;
    v /= rho;
    w /= rho;                    
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

    std::string uMaxStr("uMax");//(m/s)
    const double uMax = lookup<double>(lines, uMaxStr);

    std::string rho0Str("rho0");
    const double rho0 = lookup<double>(lines, rho0Str);

    std::string ReStr("Re");
    const double Re = lookup<double>(lines, ReStr);
    std::cout << std::endl;

    inputFile.close();

    // Single Relaxation Time model
    // // -- For cavity flow --
    const double a = 1.0; //(m)
    const double L = a/double(nx-1);

    const double u0 = 0.1;
    const double c = uMax/u0;

    const double deltaT = L/c;

    double nu = uMax*a/Re;
    nu = nu/(L*c);
    double dpdx = 0.0;


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

    double omega = 1.0/(3.0*nu +0.5);
    // double omega = 1.0/0.56;

    std::cout << "tau = " << 1.0/omega << std::endl;

    const std::vector<double> wt = {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

    // D3Q19 model
    const std::vector<double> cx = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};
    const std::vector<double> cy = {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0};
    const std::vector<double> cz = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0};

    std::vector<double> f(19*nx*ny*nz);
    std::vector<double> ftmp(19*nx*ny*nz);
    
    std::vector<double> rho(nx*ny*nz);
    std::vector<double> u(nx*ny*nz);
    std::vector<double> v(nx*ny*nz);
    std::vector<double> w(nx*ny*nz);

    std::vector<obstructure> obst(nx*ny*nz);
    

    // Setting conditions
    #include "boundaryConditionCavityFlow2D.hpp"
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
    std::cout << "Speed: " << double(endTimeStep)*double(nx*ny*nz)/time << " (LUPS)" << std::endl;

    return EXIT_SUCCESS;
}
