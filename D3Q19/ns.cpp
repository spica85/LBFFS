#include <boost/algorithm/string.hpp>
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
    double u0;
    double v0;
    double w0;
    double rho0;
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

void streaming(const int ic, const int i, const int j, const int k,const int nx, const int ny, const int nz, std::vector<double>* pf, std::vector<double>* pftmp)
{
    double& Vc    = (*pftmp)[idf(0, ic,nx,ny,nz)]; //(0)
    double& Vin   = (*pftmp)[idf(1, ic,nx,ny,nz)]; //(+x)
    double& Vip   = (*pftmp)[idf(2, ic,nx,ny,nz)]; //(-x)
    double& Vjn   = (*pftmp)[idf(3, ic,nx,ny,nz)]; //(+y)
    double& Vjp   = (*pftmp)[idf(4, ic,nx,ny,nz)]; //(-y)
    double& Vkn   = (*pftmp)[idf(5, ic,nx,ny,nz)]; //(+z)
    double& Vkp   = (*pftmp)[idf(6, ic,nx,ny,nz)]; //(-z)
    double& Vinjn = (*pftmp)[idf(7, ic,nx,ny,nz)]; //(+x,+y)
    double& Vipjp = (*pftmp)[idf(8, ic,nx,ny,nz)]; //(-x,-y)
    double& Vinjp = (*pftmp)[idf(9, ic,nx,ny,nz)]; //(+x,-y)
    double& Vipjn = (*pftmp)[idf(10,ic,nx,ny,nz)];//(-x,+y)
    double& Vinkn = (*pftmp)[idf(11,ic,nx,ny,nz)];//(+x,+z)
    double& Vipkp = (*pftmp)[idf(12,ic,nx,ny,nz)];//(-x,-z)
    double& Vinkp = (*pftmp)[idf(13,ic,nx,ny,nz)];//(+x,-z)
    double& Vipkn = (*pftmp)[idf(14,ic,nx,ny,nz)];//(-x,+z)
    double& Vjnkn = (*pftmp)[idf(15,ic,nx,ny,nz)];//(+y,+z)
    double& Vjpkp = (*pftmp)[idf(16,ic,nx,ny,nz)];//(-y,-z)
    double& Vjnkp = (*pftmp)[idf(17,ic,nx,ny,nz)];//(+y,-z)
    double& Vjpkn = (*pftmp)[idf(18,ic,nx,ny,nz)];//(-y,+z)

    int in = (i != 0)    ? index1d(i-1, j, k, nx, ny) : index1d(nx-1, j, k, nx, ny);
    int ip = (i != nx-1) ? index1d(i+1, j, k, nx, ny) : index1d(0, j, k, nx, ny);
    int jn = (j != 0 )   ? index1d(i, j-1, k, nx, ny) : index1d(i, ny-1, k, nx, ny);
    int jp = (j != ny-1) ? index1d(i, j+1, k, nx, ny) : index1d(i, 0, k, nx, ny);
    int kn = (k != 0 )   ? index1d(i, j, k-1, nx, ny) : index1d(i, j, nz-1, nx, ny);
    int kp = (k != nz-1) ? index1d(i, j, k+1, nx, ny) : index1d(i, j, 0, nx, ny);

    int injn = (i != 0 && j != 0) ? index1d(i-1, j-1, k, nx, ny) :
               (i == 0 && j != 0) ? index1d(nx-1, j-1, k, nx, ny) :
               (i != 0 && j == 0) ? index1d(i-1, ny-1, k, nx, ny) :
               index1d(nx-1, ny-1, k, nx, ny);
    int injp = (i != 0 && j != ny-1) ? index1d(i-1, j+1, k, nx, ny) :
               (i == 0 && j != ny-1) ? index1d(nx-1, j+1, k, nx, ny) :
               (i != 0 && j == ny-1) ? index1d(i-1, 0, k, nx, ny) :
               index1d(nx-1, 0, k, nx, ny);
    int inkn = (i != 0 && k != 0) ? index1d(i-1, j, k-1, nx, ny) :
               (i == 0 && k != 0) ? index1d(nx-1, j, k-1, nx, ny) :
               (i != 0 && k == 0) ? index1d(i-1, j, nz-1, nx, ny) :
               index1d(nx-1, j, nz-1, nx, ny);
    int inkp = (i != 0 && k != nz-1) ? index1d(i-1, j, k+1, nx, ny) :
               (i == 0 && k != nz-1) ? index1d(nx-1, j, k+1, nx, ny) :
               (i != 0 && k == nz-1) ? index1d(i-1, j, 0, nx, ny) :
               index1d(nx-1, j, 0, nx, ny);

    int ipjn = (i != nx-1 && j != 0) ? index1d(i+1, j-1, k, nx, ny) :
               (i == nx-1 && j != 0) ? index1d(0, j-1, k, nx, ny) :
               (i != nx-1 && j == 0) ? index1d(i+1, ny-1, k, nx, ny) :
               index1d(0, ny-1, k, nx, ny);
    int ipjp = (i != nx-1 && j != ny-1) ? index1d(i+1, j+1, k, nx, ny) :
               (i == nx-1 && j != ny-1) ? index1d(0, j+1, k, nx, ny) :
               (i != nx-1 && j == ny-1) ? index1d(i+1, 0, k, nx, ny) :
               index1d(0, 0, k, nx, ny);
    int ipkn = (i != nx-1 && k != 0) ? index1d(i+1, j, k-1, nx, ny) :
               (i == nx-1 && k != 0) ? index1d(0, j, k-1, nx, ny) :
               (i != nx-1 && k == 0) ? index1d(i+1, j, nz-1, nx, ny) :
               index1d(0, j, nz-1, nx, ny);
    int ipkp = (i != nx-1 && k != nz-1) ? index1d(i+1, j, k+1, nx, ny) :
               (i == nx-1 && k != nz-1) ? index1d(0, j, k+1, nx, ny) :
               (i != nx-1 && k == nz-1) ? index1d(i+1, j, 0, nx, ny) :
               index1d(0, j, 0, nx, ny);

    int jnkn = (j != 0 && k != 0) ? index1d(i, j-1, k-1, nx, ny) :
               (j == 0 && k != 0) ? index1d(i, ny-1, k-1, nx, ny) :
               (j != 0 && k == 0) ? index1d(i, j-1, nz-1, nx, ny) :
               index1d(i, ny-1, nz-1, nx, ny);
    int jnkp = (j != 0 && k != nz-1) ? index1d(i, j-1, k+1, nx, ny) :
               (j == 0 && k != nz-1) ? index1d(i, ny-1, k+1, nx, ny) :
               (j != 0 && k == nz-1) ? index1d(i, j-1, 0, nx, ny) :
               index1d(i, ny-1, 0, nx, ny);
    int jpkn = (j != ny-1 && k != 0) ? index1d(i, j+1, k-1, nx, ny) :
               (j == ny-1 && k != 0) ? index1d(i, 0, k-1, nx, ny) :
               (j != ny-1 && k == 0) ? index1d(i, j+1, nz-1, nx, ny) :
               index1d(i, 0, nz-1, nx, ny);
    int jpkp = (j != ny-1 && k != nz-1) ? index1d(i, j+1, k+1, nx, ny) :
               (j == ny-1 && k != nz-1) ? index1d(i, 0, k+1, nx, ny) :
               (j != ny-1 && k == nz-1) ? index1d(i, j+1, 0, nx, ny) :
               index1d(i, 0, 0, nx, ny);

    Vc    = (*pf)[idf(0, ic  ,nx,ny,nz)];
    Vin   = (*pf)[idf(1, in  ,nx,ny,nz)];
    Vip   = (*pf)[idf(2, ip  ,nx,ny,nz)];
    Vjn   = (*pf)[idf(3, jn  ,nx,ny,nz)];
    Vjp   = (*pf)[idf(4, jp  ,nx,ny,nz)];
    Vkn   = (*pf)[idf(5, kn  ,nx,ny,nz)];
    Vkp   = (*pf)[idf(6, kp  ,nx,ny,nz)];
    Vinjn = (*pf)[idf(7, injn,nx,ny,nz)];
    Vipjp = (*pf)[idf(8, ipjp,nx,ny,nz)];
    Vinjp = (*pf)[idf(9, injp,nx,ny,nz)];
    Vipjn = (*pf)[idf(10,ipjn,nx,ny,nz)];
    Vinkn = (*pf)[idf(11,inkn,nx,ny,nz)];
    Vipkp = (*pf)[idf(12,ipkp,nx,ny,nz)];
    Vinkp = (*pf)[idf(13,inkp,nx,ny,nz)];
    Vipkn = (*pf)[idf(14,ipkn,nx,ny,nz)];
    Vjnkn = (*pf)[idf(15,jnkn,nx,ny,nz)];
    Vjpkp = (*pf)[idf(16,jpkp,nx,ny,nz)];
    Vjnkp = (*pf)[idf(17,jnkp,nx,ny,nz)];
    Vjpkn = (*pf)[idf(18,jpkn,nx,ny,nz)];
}

void boundaryConditions(obstructure& obst, std::vector<double>* pftmp, int ic, int nx, int ny, int nz)
{
    double& Vc    = (*pftmp)[idf(0, ic,nx,ny,nz)];//(0)
    double& Vin   = (*pftmp)[idf(1, ic,nx,ny,nz)];//(+x)
    double& Vip   = (*pftmp)[idf(2, ic,nx,ny,nz)];//(-x)
    double& Vjn   = (*pftmp)[idf(3, ic,nx,ny,nz)];//(+y)
    double& Vjp   = (*pftmp)[idf(4, ic,nx,ny,nz)];//(-y)
    double& Vkn   = (*pftmp)[idf(5, ic,nx,ny,nz)];//(+z)
    double& Vkp   = (*pftmp)[idf(6, ic,nx,ny,nz)];//(-z)
    double& Vinjn = (*pftmp)[idf(7, ic,nx,ny,nz)];//(+x,+y)
    double& Vipjp = (*pftmp)[idf(8, ic,nx,ny,nz)];//(-x,-y)
    double& Vinjp = (*pftmp)[idf(9, ic,nx,ny,nz)];//(+x,-y)
    double& Vipjn = (*pftmp)[idf(10,ic,nx,ny,nz)];//(-x,+y)
    double& Vinkn = (*pftmp)[idf(11,ic,nx,ny,nz)];//(+x,+z)
    double& Vipkp = (*pftmp)[idf(12,ic,nx,ny,nz)];//(-x,-z)
    double& Vinkp = (*pftmp)[idf(13,ic,nx,ny,nz)];//(+x,-z)
    double& Vipkn = (*pftmp)[idf(14,ic,nx,ny,nz)];//(-x,+z)
    double& Vjnkn = (*pftmp)[idf(15,ic,nx,ny,nz)];//(+y,+z)
    double& Vjpkp = (*pftmp)[idf(16,ic,nx,ny,nz)];//(-y,-z)
    double& Vjnkp = (*pftmp)[idf(17,ic,nx,ny,nz)];//(+y,-z)
    double& Vjpkn = (*pftmp)[idf(18,ic,nx,ny,nz)];//(-y,+z)

    if(obst.boundary == 1) //bounce back
    {
        if(obst.normal == -1) //i=0
        {
            Vin   = Vip;
            Vinjn = Vipjp;
            Vinjp = Vipjn;
            Vinkn = Vipkp;
            Vinkp = Vipkn;
        }
        else if(obst.normal == 1) //i=nx-1
        {
            Vip   = Vin;
            Vipjp = Vinjn;
            Vipjn = Vinjp;
            Vipkp = Vinkn;
            Vipkn = Vinkp;
        }
        else if(obst.normal == -2) //j=0
        {
            Vjn  = Vjp;
            Vinjn = Vipjp;
            Vipjn = Vinjp;
            Vjnkn = Vjpkp;
            Vjnkp = Vjpkn;
        }
        else if(obst.normal == 2) //j=ny-1
        {
            Vjp  = Vjn;
            Vinjp = Vipjn;
            Vipjp = Vinjn;
            Vjpkn = Vjnkp;
            Vjpkp = Vjnkn;
        }
        else if(obst.normal == -3) //k=0
        {
            Vkn   = Vkp;
            Vinkn = Vipkp;
            Vipkn = Vinkp;
            Vjnkn = Vjpkp;
            Vjpkn = Vjnkp;
        }
        else if(obst.normal == 3) //k=nz-1
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
        double u = obst.u0;
        double v = obst.v0;
        double w = obst.w0;

        // if(obst.normal == -1)
        // {
        //     double rho =(V[0]+V[3]+V[4]+V[5]+V[6]+V[13]+V[14]+V[17]+V[18]+2.0*(V[2]+V[8]+V[9]+Vipkp+V[15]))/(1.0+u);
        //     V[1]  = V[2]  +rho*u/3.0;
        //     V[10] = V[9]  +rho*u/6.0 +0.5*(V[3]+V[13]+V[18] -(V[4]+V[14]+V[17])) -rho*v/2.0;
        //     V[7]  = V[8]  +rho*u/6.0 -0.5*(V[3]+V[13]+V[18] -(V[4]+V[14]+V[17])) +rho*v/2.0;
        //     V[11] = V[12] +rho*u/6.0 -0.5*(V[5]+V[9]+V[13]  -(V[6]+V[14]+V[18])) +rho*w/2.0;
        //     V[16] = V[15] +rho*u/6.0 +0.5*(V[5]+V[9]+V[13]  -(V[6]+V[14]+V[18])) -rho*w/2.0;
        // }
        // else if(obst.normal == 1)
        // {
        //     double rho =(V[0]+V[3]+V[4]+V[5]+V[6]+V[13]+V[14]+V[17]+V[18]+2.0*(V[1]+V[7]+V[10]+V[11]+V[16]))/(1.0-u);
        //     V[2]  = V[1]  -rho*u/3.0;
        //     V[9]  = V[10] -rho*u/6.0 -0.5*(V[3]+V[13]+V[18] -(V[4]+V[14]+V[17])) +rho*v/2.0;
        //     V[8]  = V[7]  -rho*u/6.0 +0.5*(V[3]+V[13]+V[18] -(V[4]+V[14]+V[17])) -rho*v/2.0;
        //     V[12] = V[11] -rho*u/6.0 +0.5*(V[5]+V[9]+V[13]  -(V[6]+V[14]+V[18])) -rho*w/2.0;
        //     V[15] = V[16] -rho*u/6.0 -0.5*(V[5]+V[9]+V[13]  -(V[6]+V[14]+V[18])) +rho*w/2.0;
        // }
        // else if(obst.normal == -2)
        // {
        //     double rho =(V[0]+V[1]+V[2]+V[5]+V[6]+V[11]+V[12]+V[15]+V[16] +2.0*(V[4]+V[8]+V[10]+V[14]+V[17]))/(1.0-v);
        //     V[3]  = V[4]  +rho*v/3.0;
        //     V[7]  = V[8]  +rho*v/6.0 -0.5*(V[1]+V[11]+V[16] -(V[2]+V[12]+V[15])) +rho*u/2.0;
        //     V[9]  = V[10] +rho*v/6.0 +0.5*(V[1]+V[11]+V[16] -(V[2]+V[12]+V[15])) -rho*u/2.0;
        //     V[13] = V[14] +rho*v/6.0 -0.5*(V[5]+V[11]+V[15] -(V[6]+V[12]+V[16])) +rho*w/2.0;
        //     V[18] = V[17] +rho*v/6.0 +0.5*(V[5]+V[11]+V[15] -(V[6]+V[12]+V[16])) -rho*w/2.0;
        // }
        // else if(obst.normal == 2)
        // {
        //     double rho =(V[0]+V[1]+V[2]+V[5]+V[6]+V[11]+V[12]+V[15]+V[16] +2.0*(V[3]+V[7]+V[9]+V[13]+V[18]))/(1.0+v);
        //     V[4]      = V[3]  -rho*v/3.0;
        //     V[8]  = V[7]  -rho*v/6.0 +0.5*(V[1]+V[11]+V[16] -(V[2]+V[12]+V[15])) -rho*u/2.0;
        //     V[10] = V[9]  -rho*v/6.0 -0.5*(V[1]+V[11]+V[16] -(V[2]+V[12]+V[15])) +rho*u/2.0;
        //     V[14] = V[13] -rho*v/6.0 +0.5*(V[5]+V[11]+V[15] -(V[6]+V[12]+V[16])) -rho*w/2.0;
        //     V[17] = V[18] -rho*v/6.0 -0.5*(V[5]+V[11]+V[15] -(V[6]+V[12]+V[16])) +rho*w/2.0;
        // }
        // else if(obst.normal == -3)
        // {
        //     double rho = (Vc_+Vin_+Vip_+Vjn_+Vjp_+Vinjn_+Vipjp_+Vinjp_+Vipjn_+2.0*(Vkp_+Vipkp_+Vipkn_+Vjpkp_+Vjpkn_))/(1.0-w);
        //     Vkn  = Vkp_  +rho*w/3.0;
        //     Vinkn = Vipkp_ +rho*w/6.0 -0.5*(Vin_+Vinjn_+Vipjn_ -(Vip_+Vinjp_+Vipjp_))  +rho*u/2.0;
        //     Vjnkn = Vjpkp_ +rho*w/6.0 +0.5*(Vin_+Vinjn_+Vipjn_ -(Vip_+Vinjp_+Vipjp_))  -rho*u/2.0;
        //     Vinkp = Vipkn_ +rho*w/6.0 -0.5*(Vjn_+Vinjn_+Vinjp_  -(Vjp_+Vipjp_+Vipjn_)) +rho*v/2.0;
        //     Vjnkp = Vjpkn_ +rho*w/6.0 +0.5*(Vjn_+Vinjn_+Vinjp_  -(Vjp_+Vipjp_+Vipjn_)) -rho*v/2.0;
        // }
        // else 
        if(obst.normal == 3)
        {
            double rho = (Vc+Vin+Vip+Vjn+Vjp+Vinjn+Vipjp+Vinjp+Vipjn+2.0*(Vkp+Vipkp+Vinkp+Vjpkp+Vjnkp))/(1.0-w);
            Vkp   = Vkn   -rho*w/3.0;
            Vipkp = Vinkn -rho*w/6.0 +0.5*(Vin+Vinjn+Vinjp -(Vip+Vipjp+Vipjn))  -rho*u/2.0;
            Vinkp = Vipkn -rho*w/6.0 -0.5*(Vin+Vinjn+Vinjp  -(Vip+Vipjp+Vipjn)) +rho*u/2.0;
            Vjpkp = Vjnkn -rho*w/6.0 +0.5*(Vjn+Vinjn+Vipjn -(Vjp+Vipjp+Vinjp))  -rho*v/2.0;
            Vjnkp = Vjpkn -rho*w/6.0 -0.5*(Vjn+Vinjn+Vipjn  -(Vjp+Vipjp+Vinjp)) +rho*v/2.0;
        }
    }
}

void collision(const double omega, const int ic, const int nx, const int ny, const int nz, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, std::vector<double>* pftmp)
{
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double w = 0.0;
    for(int q = 0; q < 19; q++)
    {
        int idc = idf(q,ic,nx,ny,nz);

        rho += (*pftmp)[idc];

        u += (*pftmp)[idc]*cx[q];
        v += (*pftmp)[idc]*cy[q];
        w += (*pftmp)[idc]*cz[q];
    }
    u /= rho;
    v /= rho;
    w /= rho;

    for(int q = 0; q < 19; q++)
    {
        double uSqr =u*u+v*v+w*w;
        double uDotC = u*cx[q]+v*cy[q]+w*cz[q];
        double feq = (1.0+3.0*uDotC +4.5*uDotC*uDotC -1.5*uSqr)*wt[q]*rho;

        int idc = idf(q,ic,nx,ny,nz);

        (*pftmp)[idc] += -omega *((*pftmp)[idc] -feq);
    }
}

void externalForce(const double dpdx, const int ic, const int nx, const int ny, const int nz, const std::vector<double>& cx, const std::vector<double>& cy, const std::vector<double>& cz, const std::vector<double>& wt, std::vector<double>* pftmp)
{
    double rho = 0.0;
    for(int q = 0; q < 19; q++)
    {
        int idc = idf(q,ic,nx,ny,nz);
        rho += (*pftmp)[idc];
    }
    for(int q = 0; q < 19; q++)
    {
        int idc = idf(q,ic,nx,ny,nz);
        (*pftmp)[idc] += rho*wt[q]*3.0*dpdx*cx[q];
    }
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

    std::string u0Str("u0");
    const double u0 = lookup<double>(lines, u0Str);

    std::string rho0Str("rho0");
    const double rho0 = lookup<double>(lines, rho0Str);

    std::string ReStr("Re");
    const double Re = lookup<double>(lines, ReStr);
    std::cout << std::endl;

    inputFile.close();

    // Single Relaxation Time model

    // For cavity flow
    // double nu = std::abs(u0)*double(nx)/Re;

    //For channel flow
    // double Retau = 10;
    // double Retau = 1;
    // double utau = 0.005;
    // double nu = utau*0.5*ny/Retau;
    // double dpdx = utau*utau/(0.5*ny);

    //For Poiseuille flow
    double umax = 0.1;
    double h = 1.0;
    double nu = umax*h/Re;
    // double dpdx = umax/(h*h)*8.0*nu/(ny-1);
    double dpdx = 0.00001;

    std::cout << "dpdx = " << dpdx << std::endl;

    std::cout << "nu = " << nu << std::endl;
    
    nu = nu*(ny-1);

    // double omega = 1.0/(3.0*nu +0.5);
    double omega = 1.0/0.56;

    std::cout << "tau = " << 1.0/omega << std::endl;

    const std::vector<double> wt = {1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

    // D3Q19 model
    const std::vector<double> cx = {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0};
    const std::vector<double> cy = {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0};
    const std::vector<double> cz = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0};

    // for(int i = 0; i < 19; i++)
    // {
    //     std::cout << cx[i] << ", "
    //           << cy[i] << ", "
    //           << cz[i] << ", "
    //           << std::endl;
    // }

    std::vector<double> f(19*nx*ny*nz);
    std::vector<double> ftmp(19*nx*ny*nz);

    std::vector<double>* pf = &f;
    std::vector<double>* pftmp = &ftmp;
    std::vector<double>* ptmp;
    
    std::vector<double> rho(nx*ny*nz);
    std::vector<double> u(nx*ny*nz);
    std::vector<double> v(nx*ny*nz);
    std::vector<double> w(nx*ny*nz);

    std::vector<obstructure> obst(nx*ny*nz);
    

    // Setting conditions
    // #include "boundaryCondition.hpp"
    #include "boundaryCondition_channelFlow.hpp"
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

                    streaming(ic,i,j,k,nx,ny,nz,pf,pftmp);
                    collision(omega,ic,nx,ny,nz,cx,cy,cz,wt,pftmp);
                    externalForce(dpdx,ic,nx,ny,nz,cx,cy,cz,wt,pftmp);
                    boundaryConditions(obst[ic],pftmp,ic,nx,ny,nz);
                }
            }
        }

        ptmp = pftmp;
        pftmp = pf;
        pf = ptmp;
    }

    return EXIT_SUCCESS;
}
