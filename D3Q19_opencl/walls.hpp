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

void readSTL(const std::string STLname, std::vector<std::vector<float> >& STLnormal, std::vector<std::vector<float> >& STLv0, std::vector<std::vector<float> >& STLv1, std::vector<std::vector<float> >& STLv2, int& nSTL, std::vector<std::vector<float> >& STLc, float L)
{
    std::ifstream STLfile(STLname);
    if(!STLfile)
    {
        std::cout << "\nSTL (" << STLname << ") was not read\n" << std::endl;
        nSTL = 0;
        return;
    }
    else
    {
        std::cout << "\nSTL (" << STLname << ") was read" << std::endl;
    }

    std::vector<std::string> lines;
    readToLines(STLfile, lines);

    for(int i = 0; i < lines.size(); i++)
    {
        // std::cout << i << " " << lines[i] << std::endl;
        if(lines[i] == "facet")
        {
            STLnormal[0].push_back(std::stof(lines[i+2]));
            STLnormal[1].push_back(std::stof(lines[i+3]));
            STLnormal[2].push_back(std::stof(lines[i+4]));
            // std::cout << i << " " 
            //     << lines[i+2] << " " 
            //     << lines[i+3] << " " 
            //     << lines[i+4] << " " 
            //     << std::endl;
            i += 4;
        }
        else if(lines[i] == "vertex")
        {
            // std::cout << i << " " 
            //     << lines[i+1] << " " 
            //     << lines[i+2] << " " 
            //     << lines[i+3] << " " 
            //     << std::endl;
            STLv0[0].push_back(std::stof(lines[i+1]));
            STLv0[1].push_back(std::stof(lines[i+2]));
            STLv0[2].push_back(std::stof(lines[i+3]));
            i += 5;

            // std::cout << i << " " 
            //     << lines[i+1] << " " 
            //     << lines[i+2] << " " 
            //     << lines[i+3] << " " 
            //     << std::endl;
            STLv1[0].push_back(std::stof(lines[i+1]));
            STLv1[1].push_back(std::stof(lines[i+2]));
            STLv1[2].push_back(std::stof(lines[i+3]));
            i += 5;

            // std::cout << i << " " 
            //     << lines[i+1] << " " 
            //     << lines[i+2] << " " 
            //     << lines[i+3] << " " 
            //     << std::endl;
            STLv2[0].push_back(std::stof(lines[i+1]));
            STLv2[1].push_back(std::stof(lines[i+2]));
            STLv2[2].push_back(std::stof(lines[i+3]));
            i += 3;
        }
    }

    nSTL = STLnormal[0].size();
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

        STLc[0].push_back((STLv0[0][i]+STLv1[0][i]+STLv2[0][i])/3.f);
        STLc[1].push_back((STLv0[1][i]+STLv1[1][i]+STLv2[1][i])/3.f);
        STLc[2].push_back((STLv0[2][i]+STLv1[2][i]+STLv2[2][i])/3.f);
    }

    std::cout << "Number of elements of " << STLname << ": " << nSTL << std::endl;
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

void setSolid(std::vector<unsigned char>& solid, std::vector<float>& sdf, const float sdfIni, const int nx, const int ny, const int nz)
{
    const int elements = sdf.size();

    for(int ic = 0; ic < elements; ic++)
    {
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

void setQf(std::vector<float>& qf, std::vector<unsigned char>& neiSolid, std::vector<float>& sdf, std::vector<unsigned char>& solid, const int nx, const int ny, const int nz)
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
                if(solid[ic] == 0 && solid[upID] == 1)
                {
                    neiSolid[ic] = 1;
                    qf[qic] = abs(sdf0)/(abs(sdf0)+abs(sdf1));
                }
            }
        }
    }
}



#endif