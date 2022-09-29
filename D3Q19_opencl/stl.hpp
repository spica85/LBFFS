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
    std::vector<std::vector<float> > STLv0;
    std::vector<std::vector<float> > STLv1;
    std::vector<std::vector<float> > STLv2;

    void setPoint(std::vector<std::vector<float> >& point, float x, float y, float z)
    {
        point[0].push_back(x);
        point[1].push_back(y);
        point[2].push_back(z);
    }

    public:
    const std::string patchName;
    std::vector<std::vector<float> > STLnormal;
    std::vector<std::vector<float> > STLc;

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
        }

        return 1;
    }

    public:
    std::vector<std::unique_ptr<STLpatch> > patch;
    
    STL(const std::string fileName, const float L): fileName(fileName), L(L)
    {
        readSTL();
        for(int iPatch = 0; iPatch < patch.size(); iPatch++)
        {
            std::cout << "Number of elements of " << patch[iPatch]->patchName << ": " << patch[iPatch]->nSize() << std::endl;
        }
    }

};


#endif