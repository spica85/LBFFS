#ifndef stl_H
#define stl_H

#include <vector>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <omp.h>
#include "input.hpp"
#include <memory>

#define _USE_MATH_DEFINES
#include <cmath>

int readSTL(const std::string STLname, std::vector<std::vector<float> >& STLnormal, std::vector<std::vector<float> >& STLv0, std::vector<std::vector<float> >& STLv1, std::vector<std::vector<float> >& STLv2, int& nSTL, std::vector<std::vector<float> >& STLc, float L)
{
    std::ifstream STLfile(STLname);
    if(!STLfile)
    {
        std::cout << "\nSTL (" << STLname << ") was not read\n" << std::endl;
        nSTL = 0;
        return 0;
    }
    else
    {
        std::cout << "\nSTL (" << STLname << ") was read" << std::endl;
    }

    std::vector<std::string> lines;
    readToLines(STLfile, lines);

    for(int i = 0; i < lines.size(); i++)
    {
        if(lines[i] == "facet")
        {
            STLnormal[0].push_back(std::stof(lines[i+2]));
            STLnormal[1].push_back(std::stof(lines[i+3]));
            STLnormal[2].push_back(std::stof(lines[i+4]));
            i += 4;
        }
        else if(lines[i] == "vertex")
        {
            STLv0[0].push_back(std::stof(lines[i+1]));
            STLv0[1].push_back(std::stof(lines[i+2]));
            STLv0[2].push_back(std::stof(lines[i+3]));
            i += 5;

            STLv1[0].push_back(std::stof(lines[i+1]));
            STLv1[1].push_back(std::stof(lines[i+2]));
            STLv1[2].push_back(std::stof(lines[i+3]));
            i += 5;

            STLv2[0].push_back(std::stof(lines[i+1]));
            STLv2[1].push_back(std::stof(lines[i+2]));
            STLv2[2].push_back(std::stof(lines[i+3]));
            i += 3;
        }
    }

    nSTL = STLnormal[0].size();
    for(int i = 0; i < nSTL; i++)
    {
        STLv0[0][i] = STLv0[0][i]/L;
        STLv0[1][i] = STLv0[1][i]/L;
        STLv0[2][i] = STLv0[2][i]/L;
        STLv1[0][i] = STLv1[0][i]/L;
        STLv1[1][i] = STLv1[1][i]/L;
        STLv1[2][i] = STLv1[2][i]/L;
        STLv2[0][i] = STLv2[0][i]/L;
        STLv2[1][i] = STLv2[1][i]/L;
        STLv2[2][i] = STLv2[2][i]/L;

        STLc[0].push_back((STLv0[0][i]+STLv1[0][i]+STLv2[0][i])/3.f);
        STLc[1].push_back((STLv0[1][i]+STLv1[1][i]+STLv2[1][i])/3.f);
        STLc[2].push_back((STLv0[2][i]+STLv1[2][i]+STLv2[2][i])/3.f);
    }

    std::cout << "Number of elements of " << STLname << ": " << nSTL << std::endl;
    return 1;
}


class STLpatch
{
    void setPoint(std::vector<std::vector<float> >& point, float x, float y, float z)
    {
        point[0].push_back(x);
        point[1].push_back(y);
        point[2].push_back(z);
    }

    bool pointInElement(const float x, const float y, const float x1, const float y1, const float x2, const float y2)
    {
        float s = (x*y2 -x2*y)/(x1*y2 -x2*y1);
        float t = (x*y1 -x1*y)/(x2*y1 -x1*y2);

        if(s+t >= 0.f && s+t <= 1.f)
        {
            if(s >= 0.f && s <= 1.f)
            {
                if(t >= 0.f && t <= 1.f)
                {
                    return true;
                }
            }
        }
        return false;
    }

    public:
    const std::string patchName;
    std::vector<std::vector<float> > STLv0;
    std::vector<std::vector<float> > STLv1;
    std::vector<std::vector<float> > STLv2;
    std::vector<std::vector<float> > STLnormal;
    std::vector<std::vector<float> > STLc;
    std::vector<float> eleXMin;
    std::vector<float> eleXMax;
    std::vector<float> eleYMin;
    std::vector<float> eleYMax;
    std::vector<float> eleZMin;
    std::vector<float> eleZMax;
    
    float xMin;
    float xMax;
    float yMin;
    float yMax;
    float zMin;
    float zMax;

    STLpatch(const std::string patchName): 
        patchName(patchName), STLnormal(3),
        STLv0(3), STLv1(3), STLv2(3),
        STLc(3)
    { }

    
    void setSTLnormal(float x, float y, float z)
    {
        setPoint(STLnormal, x, y, z);
    }

    void setSTLv0(float x, float y, float z)
    {
        setPoint(STLv0, x, y, z);
    }

    void setSTLv1(float x, float y, float z)
    {
        setPoint(STLv1, x, y, z);
    }

    void setSTLv2(float x, float y, float z)
    {
        setPoint(STLv2, x, y, z);
    }

    void normalizeSTLv(const float L)
    {
        for(int i = 0; i < STLv0[0].size(); i++)
        {
            STLv0[0][i] /= L;
            STLv0[1][i] /= L;
            STLv0[2][i] /= L;
            STLv1[0][i] /= L;
            STLv1[1][i] /= L;
            STLv1[2][i] /= L;
            STLv2[0][i] /= L;
            STLv2[1][i] /= L;
            STLv2[2][i] /= L;
        }
    }

    void setSTLc()
    {
        for(int i = 0; i < STLv0[0].size(); i++)
        {
            float x = (STLv0[0][i] +STLv1[0][i] +STLv2[0][i])/3.f;
            float y = (STLv0[1][i] +STLv1[1][i] +STLv2[1][i])/3.f;
            float z = (STLv0[2][i] +STLv1[2][i] +STLv2[2][i])/3.f;
            setPoint(STLc, x, y, z);
        }
    }

    void setBoundingBox()
    {
        float v0xMin = *min_element(STLv0[0].begin(),STLv0[0].end());
        float v0yMin = *min_element(STLv0[1].begin(),STLv0[1].end());
        float v0zMin = *min_element(STLv0[2].begin(),STLv0[2].end());

        float v1xMin = *min_element(STLv1[0].begin(),STLv1[0].end());
        float v1yMin = *min_element(STLv1[1].begin(),STLv1[1].end());
        float v1zMin = *min_element(STLv1[2].begin(),STLv1[2].end());

        float v2xMin = *min_element(STLv2[0].begin(),STLv2[0].end());
        float v2yMin = *min_element(STLv2[1].begin(),STLv2[1].end());
        float v2zMin = *min_element(STLv2[2].begin(),STLv2[2].end());

        xMin = fmin(v0xMin, fmin(v1xMin, v2xMin));
        yMin = fmin(v0yMin, fmin(v1yMin, v2yMin));
        zMin = fmin(v0zMin, fmin(v1zMin, v2zMin));


        float v0xMax = *max_element(STLv0[0].begin(),STLv0[0].end());
        float v0yMax = *max_element(STLv0[1].begin(),STLv0[1].end());
        float v0zMax = *max_element(STLv0[2].begin(),STLv0[2].end());

        float v1xMax = *max_element(STLv1[0].begin(),STLv1[0].end());
        float v1yMax = *max_element(STLv1[1].begin(),STLv1[1].end());
        float v1zMax = *max_element(STLv1[2].begin(),STLv1[2].end());

        float v2xMax = *max_element(STLv2[0].begin(),STLv2[0].end());
        float v2yMax = *max_element(STLv2[1].begin(),STLv2[1].end());
        float v2zMax = *max_element(STLv2[2].begin(),STLv2[2].end());

        xMax = fmax(v0xMax, fmax(v1xMax, v2xMax));
        yMax = fmax(v0yMax, fmax(v1yMax, v2yMax));
        zMax = fmax(v0zMax, fmax(v1zMax, v2zMax));
    }

    void eleBoundingBox(const int iSTL, float& xMin, float& xMax, float& yMin, float& yMax, float& zMin, float& zMax)
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

        xMin = fmin(v0[0], fmin(v1[0],v2[0]));
        yMin = fmin(v0[1], fmin(v1[1],v2[1]));
        zMin = fmin(v0[2], fmin(v1[2],v2[2]));
        xMax = fmax(v0[0], fmax(v1[0],v2[0]));
        yMax = fmax(v0[1], fmax(v1[1],v2[1]));
        zMax = fmax(v0[2], fmax(v1[2],v2[2]));
    }

    void setEleBoundingBox()
    {
        for(int i = 0; i < STLnormal[0].size(); i++)
        {
            float xMin;
            float xMax;
            float yMin;
            float yMax;
            float zMin;
            float zMax;
            eleBoundingBox(i,xMin,xMax,yMin,yMax,zMin,zMax);
            eleXMin.push_back(xMin);
            eleXMax.push_back(xMax);
            eleYMin.push_back(yMin);
            eleYMax.push_back(yMax);
            eleZMin.push_back(zMin);
            eleZMax.push_back(zMax);
        }
    }


    void printSTLnormal()
    {
        for(int i = 0; i < STLnormal[0].size(); i++)
        {
            std::cout << "("
                << STLnormal[0][i] << ", "
                << STLnormal[1][i] << ", "
                << STLnormal[2][i] << ")"
                << std::endl;
        }
    }

    void printSTLv()
    {
        for(int i = 0; i < STLv0[0].size(); i++)
        {
            std::cout << "("
                << STLv0[0][i] << ", "
                << STLv0[1][i] << ", "
                << STLv0[2][i] << "), "
                << "("
                << STLv1[0][i] << ", "
                << STLv1[1][i] << ", "
                << STLv1[2][i] << "), "
                << "("
                << STLv2[0][i] << ", "
                << STLv2[1][i] << ", "
                << STLv2[2][i] << ")"
                << std::endl;
        }
    }

    void printSTLc()
    {
        for(int i = 0; i < STLc[0].size(); i++)
        {
            std::cout << "("
                << STLc[0][i] << ", "
                << STLc[1][i] << ", "
                << STLc[2][i] << ")"
                << std::endl;
        }
    }

    int nSize()
    {
        return STLv0[0].size();
    }

    void findPoints(std::vector<int>& points, const int nx, const int ny, const int nz)
    {
        for(int id = 0; id < nSize(); id++)
        {
            float maxNormal = fmax(abs(STLnormal[0][id]), fmax(abs(STLnormal[1][id]), abs(STLnormal[2][id])));

            if(abs(STLnormal[0][id]) == maxNormal) // xPlane
            {
                if(int(STLc[0][id]) == 0 || int(STLc[0][id])+1 == 0) // i = 0
                {
                    for(int j = 0; j < ny; j++)
                    {
                        for(int k = 0; k < nz; k++)
                        {
                            float x = float(j) -STLv0[1][id];
                            float y = float(k) -STLv0[2][id];
                            float x1 = STLv1[1][id] -STLv0[1][id];
                            float y1 = STLv1[2][id] -STLv0[2][id];
                            float x2 = STLv2[1][id] -STLv0[1][id];
                            float y2 = STLv2[2][id] -STLv0[2][id];

                            if(pointInElement(x, y, x1, y1, x2, y2))
                            {
                                points.push_back(index1d(0,j,k,nx,ny));
                            }
                        }
                    }
                }
                else if(int(STLc[0][id]) == nx-1 || int(STLc[0][id])+1 == nx-1) // i = nx-1
                {
                    for(int j = 0; j < ny; j++)
                    {
                        for(int k = 0; k < nz; k++)
                        {
                            float x = float(j) -STLv0[1][id];
                            float y = float(k) -STLv0[2][id];
                            float x1 = STLv1[1][id] -STLv0[1][id];
                            float y1 = STLv1[2][id] -STLv0[2][id];
                            float x2 = STLv2[1][id] -STLv0[1][id];
                            float y2 = STLv2[2][id] -STLv0[2][id];

                            if(pointInElement(x, y, x1, y1, x2, y2))
                            {
                                points.push_back(index1d(nx-1,j,k,nx,ny));
                            }
                        }
                    }
                }
            }
            else if(abs(STLnormal[1][id]) == maxNormal) // yPlane
            {
                if(int(STLc[1][id]) == 0 || int(STLc[1][id])+1 == 0) // j = 0
                {
                    for(int i = 0; i < nx; i++)
                    {
                        for(int k = 0; k < nz; k++)
                        {
                            float x = float(i) -STLv0[0][id];
                            float y = float(k) -STLv0[2][id];
                            float x1 = STLv1[0][id] -STLv0[0][id];
                            float y1 = STLv1[2][id] -STLv0[2][id];
                            float x2 = STLv2[0][id] -STLv0[0][id];
                            float y2 = STLv2[2][id] -STLv0[2][id];

                            if(pointInElement(x, y, x1, y1, x2, y2))
                            {
                                points.push_back(index1d(i,0,k,nx,ny));
                            }
                        }
                    }
                }
                else if(int(STLc[1][id]) == ny-1 || int(STLc[1][id])+1 == ny-1) // j = ny-1
                {
                    for(int i = 0; i < nx; i++)
                    {
                        for(int k = 0; k < nz; k++)
                        {
                            float x = float(i) -STLv0[0][id];
                            float y = float(k) -STLv0[2][id];
                            float x1 = STLv1[0][id] -STLv0[0][id];
                            float y1 = STLv1[2][id] -STLv0[2][id];
                            float x2 = STLv2[0][id] -STLv0[0][id];
                            float y2 = STLv2[2][id] -STLv0[2][id];

                            if(pointInElement(x, y, x1, y1, x2, y2))
                            {
                                points.push_back(index1d(i,ny-1,k,nx,ny));
                            }
                        }
                    }
                }
            }
            else // zPlane
            {
                if(int(STLc[2][id]) == 0 || int(STLc[2][id])+1 == 0) // k = 0
                {
                    for(int i = 0; i < nx; i++)
                    {
                        for(int j = 0; j < ny; j++)
                        {
                            float x = float(i) -STLv0[0][id];
                            float y = float(j) -STLv0[1][id];
                            float x1 = STLv1[0][id] -STLv0[0][id];
                            float y1 = STLv1[1][id] -STLv0[1][id];
                            float x2 = STLv2[0][id] -STLv0[0][id];
                            float y2 = STLv2[1][id] -STLv0[1][id];

                            if(pointInElement(x, y, x1, y1, x2, y2))
                            {
                                points.push_back(index1d(i,j,0,nx,ny));
                            }
                        }
                    }
                }
                else if(int(STLc[2][id]) == nz-1 || int(STLc[2][id])+1 == nz-1) // k = nz-1
                {
                    for(int i = 0; i < nx; i++)
                    {
                        for(int j = 0; j < ny; j++)
                        {
                            float x = float(i) -STLv0[0][id];
                            float y = float(j) -STLv0[1][id];
                            float x1 = STLv1[0][id] -STLv0[0][id];
                            float y1 = STLv1[1][id] -STLv0[1][id];
                            float x2 = STLv2[0][id] -STLv0[0][id];
                            float y2 = STLv2[1][id] -STLv0[1][id];

                            if(pointInElement(x, y, x1, y1, x2, y2))
                            {
                                points.push_back(index1d(i,j,nz-1,nx,ny));
                            }
                        }
                    }
                }
            }
        }
    }
};

class STL
{
    const std::string fileName;
    const float L;   
    
    int readSTL()
    {
        std::ifstream STLfile(fileName);
        if(!STLfile)
        {
            std::cout << "\nSTL (" << fileName << ") was not read\n" << std::endl;
            return 0;
        }
        else
        {
            std::cout << "\nSTL (" << fileName << ") was read" << std::endl;
        }      

        std::vector<std::string> lines;
        readToLines(STLfile, lines);

        for(int i = 0; i < lines.size(); i++)
        {
            if(lines[i] == "solid")
            {
                patch.push_back(std::unique_ptr<STLpatch>(new STLpatch(lines[i+1])));
            }
        }


        int iPatch = 0;
        for(int i = 0; i < lines.size(); i++)
        {
            if(lines[i] == "solid")
            {
                if(lines[i+1] != patch[iPatch]->patchName)
                {
                    iPatch++;
                    i++;
                }
            }

            if(lines[i] == "facet")
            {
                patch[iPatch]->setSTLnormal(std::stof(lines[i+2]), std::stof(lines[i+3]), std::stof(lines[i+4]));
                i += 4;
            }
            else if(lines[i] == "vertex")
            {
                patch[iPatch]->setSTLv0(std::stof(lines[i+1]), std::stof(lines[i+2]), std::stof(lines[i+3]));
                i += 5;

                patch[iPatch]->setSTLv1(std::stof(lines[i+1]), std::stof(lines[i+2]), std::stof(lines[i+3]));
                i += 5;

                patch[iPatch]->setSTLv2(std::stof(lines[i+1]), std::stof(lines[i+2]), std::stof(lines[i+3]));
                i += 3;
            }
        }

        for(int iPatch = 0; iPatch < patch.size(); iPatch++)
        {
            patch[iPatch]->normalizeSTLv(L);
            patch[iPatch]->setSTLc();
            patch[iPatch]->setBoundingBox();
        }

        return 1;
    }

    public:
    std::vector<std::unique_ptr<STLpatch> > patch;
    float xMin;
    float xMax;
    float yMin;
    float yMax;
    float zMin;
    float zMax;
    
    STL(const std::string fileName, const float L): fileName(fileName), L(L)
    {
        readSTL();
        xMin = patch[0]->xMin;
        xMax = patch[0]->xMax;
        yMin = patch[0]->yMin;
        yMax = patch[0]->yMax;
        zMin = patch[0]->zMin;
        zMax = patch[0]->zMax;
        for(int iPatch = 0; iPatch < patch.size(); iPatch++)
        {
            std::cout << "Number of elements of " << patch[iPatch]->patchName << ": " << patch[iPatch]->nSize() << std::endl;
            xMin = fmin(xMin, patch[iPatch]->xMin);
            xMax = fmax(xMax, patch[iPatch]->xMax);
            yMin = fmin(yMin, patch[iPatch]->yMin);
            yMax = fmax(yMax, patch[iPatch]->yMax);
            zMin = fmin(zMin, patch[iPatch]->zMin);
            zMax = fmax(zMax, patch[iPatch]->zMax);
            patch[iPatch]->setEleBoundingBox();
        }
    }

};


#endif