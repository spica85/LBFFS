#ifndef walls_H
#define walls_H

#include <vector>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <omp.h>
#include "input.hpp"
#include "stl.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

float det(const std::vector<float>& a, const std::vector<float>& b, const std::vector<float>& c)
{
    return (a[0]*b[1]*c[2]) +(a[1]*b[2]*c[0]) +(a[2]*b[0]*c[1])
            -(a[0]*b[2]*c[1]) -(a[1]*b[0]*c[2]) -(a[2]*b[1]*c[0]);
}

bool isCross(const std::vector<float>& v0, const std::vector<float>& v1, const std::vector<float>& v2, const int i, const int j, const int k)
{
    const std::vector<float> ray = {1.f, 0.f, 0.f};
    const std::vector<float> invRay = {-1.f, 0.f, 0.f};

    std::vector<float> a(3);
    a[0] = v1[0] -v0[0];
    a[1] = v1[1] -v0[1];
    a[2] = v1[2] -v0[2];
    std::vector<float> b(3);
    b[0] = v2[0] -v0[0];
    b[1] = v2[1] -v0[1];
    b[2] = v2[2] -v0[2];
    const std::vector<float> c = invRay;

    std::vector<float> d(3);
    d[0] = float(i) -v0[0];
    d[1] = float(j) -v0[1];
    d[2] = float(k) -v0[2];

    float x = det(d,b,c)/det(a,b,c);
    float y = det(a,d,c)/det(a,b,c);
    float z = det(a,b,d)/det(a,b,c);

    if(x < 0.f || x > 1.f)
    {
        return false;
    }
    else if(y < 0.f || y > 1.f)
    {
        return false;
    }
    else if(x+y < 0.f || x+y > 1.f)
    {
        return false;
    }
    else if(z < 0.f)
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool isInner(const STL& stl, const int i, const int j, const int k)
{
    if(float(i) < stl.xMin)
    {
        return false;
    }
    else if(float(i) > stl.xMax)
    {
        return false;
    }
    else if(float(j) < stl.yMin)
    {
        return false;
    }
    else if(float(j) > stl.yMax)
    {
        return false;
    }
    else if(float(k) < stl.zMin)
    {
        return false;
    }
    else if(float(k) > stl.zMax)
    {
        return false;
    }
    else
    {
        int nCross = 0;
        int nCount = 0;
        for(int iPatch = 0; iPatch < stl.patch.size(); iPatch++)
        {
            const std::unique_ptr<STLpatch>& patch = stl.patch[iPatch];
            
            std::vector<std::vector<float> >& STLv0 = patch->STLv0;
            std::vector<std::vector<float> >& STLv1 = patch->STLv1;
            std::vector<std::vector<float> >& STLv2 = patch->STLv2;
            const int nSTL = STLv0[0].size();
            
            for(int iSTL = 0; iSTL < nSTL; iSTL++)
            {
                float xMin = patch->eleXMin[iSTL];
                float xMax = patch->eleXMax[iSTL];
                float yMin = patch->eleYMin[iSTL];
                float yMax = patch->eleYMax[iSTL];
                float zMin = patch->eleZMin[iSTL];
                float zMax = patch->eleZMax[iSTL];
                
                if
                (
                    float(j) >= yMin &&
                    float(j) <= yMax &&
                    float(k) >= zMin &&
                    float(k) <= zMax
                ) 
                {
                    std::vector<float> v0(3);
                    v0[0] = STLv0[0][iSTL];
                    v0[1] = STLv0[1][iSTL];
                    v0[2] = STLv0[2][iSTL];
                    std::vector<float> v1(3);
                    v1[0] = STLv1[0][iSTL];
                    v1[1] = STLv1[1][iSTL];
                    v1[2] = STLv1[2][iSTL];
                    std::vector<float> v2(3);
                    v2[0] = STLv2[0][iSTL];
                    v2[1] = STLv2[1][iSTL];
                    v2[2] = STLv2[2][iSTL];

                    nCount++;
                    if(isCross(v0,v1,v2,i,j,k))
                    {
                        nCross++;
                    }   
                }
            }
        }
        
        if(nCross%2 == 0)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
}

void setSolid(std::vector<unsigned char>& solid, const STL& stl, const int nx, const int ny, const int nz, const bool fluid = true)
{
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                int ic = index1d(i,j,k,nx,ny);
                if(isInner(stl,i,j,k))
                {
                    solid[ic] = 1;
                }
            }
        }
    }

    const int elements = nx*ny*nz;
    if(!fluid)
    {
        for(int ic = 0; ic < elements; ic++)
        {            
            solid[ic] = (solid[ic] == 1) ? 0 : 1;
        }        
    }
}


void setNeiSolid(std::vector<unsigned char>& neiSolid, const std::vector<unsigned char>& solid, const int nx, const int ny, const int nz)
{
    const int elements = nx*ny*nz;
    for(int ic = 0; ic < elements; ic++)
    {
        int i = ic2i(ic,nx,ny);
        int j = ic2j(ic,nx,ny);
        int k = ic2k(ic,nx,ny);
        
        if
        (
            (nx == 1 && (j != 0 && j != ny-1 && k != 0 && k != nz-1)) ||
            (ny == 1 && (i != 0 && i != nx-1 && k != 0 && k != nz-1)) || 
            (nz == 1 && (i != 0 && i != nx-1 && j != 0 && j != ny-1)) ||
            (i != 0 && i != nx-1 && j != 0 && j != ny-1 && k != 0 && k != nz-1)
        )
        {
            for(int q = 0; q < 19; q++)
            {
                int qic = q*elements +ic;
                const int upID = upwindID(q,i,j,k,nx,ny,nz);
                
                if(solid[ic] == 0 && solid[upID] == 1)
                {
                    neiSolid[ic] = 1;
                }
            }
        }
    }
}

void setSDF(std::vector<float>& sdf, const float sdfIni, const float dr, const float p, const std::unique_ptr<STLpatch>& patch, const int nx, const int ny, const int nz, const bool boundary, std::vector<unsigned char>& neiSolid)
{
    const int elements = nx*ny*nz;
    const int drn = int(dr);
    
    std::vector<int> nearSTL;
    for(int ic = 0; ic < elements; ic++)
    {
        if(neiSolid[ic] == 1)
        {
            int i = ic2i(ic,nx,ny);
            int j = ic2j(ic,nx,ny);
            int k = ic2k(ic,nx,ny);

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
    std::cout << "Number of near walls: " << nearSTL.size() << std::endl;

    std::vector<std::vector<float> >& STLc = patch->STLc;
    std::vector<std::vector<float> >& STLv0 = patch->STLv0;
    std::vector<std::vector<float> >& STLv1 = patch->STLv1;
    std::vector<std::vector<float> >& STLv2 = patch->STLv2;
    std::vector<std::vector<float> >& STLnormal = patch->STLnormal;
    const int nSTL = STLc[0].size();
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
    }
}

void setSDF(std::vector<float>& sdf, const float sdfIni, const float dr, const float p, std::vector<std::vector<float> >& STLc, std::vector<std::vector<float> >& STLnormal, const int nx, const int ny, const int nz, const bool boundary)
{
    const int drn = int(dr);
    const int nSTL = STLc[0].size();

    std::vector<int> nearSTL;
    for(int iSTL = 0; iSTL < nSTL; iSTL++)
    {
        int i = int(STLc[0][iSTL]);
        int j = int(STLc[1][iSTL]);
        int k = int(STLc[2][iSTL]);
        // if(boundary || (0 < i && i < nx-1 && 0 < j && j < ny-1 && 0 < k && k < nz-1))
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

    std::cout << "Number of near walls: " << nearSTL.size() << std::endl;
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
    }
}

void setBwalls(std::vector<std::vector<float> >& bWallsC, std::vector<std::vector<float> >& bWallsNormal, const std::vector<int>& boundary1, const std::vector<int>& boundary2, const std::vector<int>& boundary3, const int nx, const int ny, const int nz)
{
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
}

void setQf(std::vector<float>& qf, const std::vector<unsigned char>& neiSolid, std::vector<float>& sdf, const int nx, const int ny, const int nz)
{
    const int elements = nx*ny*nz;
    for(int ic = 0; ic < elements; ic++)
    {
        int i = ic2i(ic,nx,ny);
        int j = ic2j(ic,nx,ny);
        int k = ic2k(ic,nx,ny);
        
        if
        (
            (nx == 1 && (j != 0 && j != ny-1 && k != 0 && k != nz-1)) ||
            (ny == 1 && (i != 0 && i != nx-1 && k != 0 && k != nz-1)) || 
            (nz == 1 && (i != 0 && i != nx-1 && j != 0 && j != ny-1)) ||
            (i != 0 && i != nx-1 && j != 0 && j != ny-1 && k != 0 && k != nz-1)
        )
        {
            for(int q = 0; q < 19; q++)
            {
                int qic = q*elements +ic;
                const float sdf0 = sdf[ic];
                const int upID = upwindID(q,i,j,k,nx,ny,nz);
                
                const float sdf1 = sdf[upID];
                if(neiSolid[ic] == 1)
                {
                    qf[qic] = abs(sdf0)/(abs(sdf0)+abs(sdf1));
                }
            }
        }
    }
}

void readMotions
(
    float& uMovingTrans, float& vMovingTrans, float& wMovingTrans,
    float& rotOmega,
    float& rotX, float& rotY, float& rotZ,
    float& rotAxisX, float& rotAxisY, float& rotAxisZ,
    const float c, const float L
)
{
    std::string inputFileName("input.txt");
    std::vector<std::string> lines;
    std::ifstream inputFile(inputFileName);
    if(!inputFile)
    {
        std::cerr << "Could not open input.txt" << std::endl;
        exit(EXIT_FAILURE);
    }

    readToLines(inputFile, lines);

    std::string uMovingTransStr("uMovingTrans");
    uMovingTrans = lookup<float>(lines, uMovingTransStr)/c;
    std::string vMovingTransStr("vMovingTrans");
    vMovingTrans = lookup<float>(lines, vMovingTransStr)/c;
    std::string wMovingTransStr("wMovingTrans");
    wMovingTrans = lookup<float>(lines, wMovingTransStr)/c;

    std::string rotOmegaStr("rotOmega");// rpm
    rotOmega = lookup<float>(lines, rotOmegaStr)*2.f*M_PI/60.f/(c/L);// rad/s
    std::string rotXStr("rotX");
    rotX = lookup<float>(lines, rotXStr)/L;
    std::string rotYStr("rotY");
    rotY = lookup<float>(lines, rotYStr)/L;
    std::string rotZStr("rotZ");
    rotZ = lookup<float>(lines, rotZStr)/L;
    std::string rotAxisXStr("rotAxisX");
    rotAxisX = lookup<float>(lines, rotAxisXStr);
    std::string rotAxisYStr("rotAxisY");
    rotAxisY = lookup<float>(lines, rotAxisYStr);
    std::string rotAxisZStr("rotAxisZ");
    rotAxisZ = lookup<float>(lines, rotAxisZStr);

    inputFile.close();
}



#endif