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


#endif