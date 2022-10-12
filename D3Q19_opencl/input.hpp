#ifndef input_H
#define input_H

#include <vector>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>

template <typename Type>
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

template bool lookup<bool>(std::vector<std::string>& lines, std::string& str);
template int lookup<int>(std::vector<std::string>& lines, std::string& str);
template double lookup<double>(std::vector<std::string>& lines, std::string& str);
template std::string lookup<std::string>(std::vector<std::string>& lines, std::string& str);

template<typename Type> Type lookupOrDefault(std::vector<std::string>& lines, std::string& str, Type defaultValue)
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
    return defaultValue;
}

template bool lookupOrDefault<bool>(std::vector<std::string>& lines, std::string& str, bool defaultValue);
template int lookupOrDefault<int>(std::vector<std::string>& lines, std::string& str, int defaultValue);
template double lookupOrDefault<double>(std::vector<std::string>& lines, std::string& str, double defaultValue);
template std::string lookupOrDefault<std::string>(std::vector<std::string>& lines, std::string& str, std::string defaultValue);


void readToLines(std::ifstream& inputFile, std::vector<std::string>& lines)
{
    std::string line;

    // if(!inputFile.is_open())
    // {
    //     std::cerr << "Could not open the file" << std::endl;
    //     exit(EXIT_FAILURE);
    // }


    while(std::getline(inputFile, line))
    {
        std::vector<std::string> list_string;
        boost::split(list_string,line, boost::is_space(), boost::algorithm::token_compress_on);
        for(auto& str: list_string)
        {
            lines.push_back(str);
        }
    }
}




void input(bool& restart, bool& Fwrite, bool& writeBinary, int& startTimeStep, int& endTimeStep, int& nextOutTime, int& outInterval, int& nx, int& ny, int& nz, float& Lx, float& uMax, float& rho0, float& U0, float& nu, float& dpdx, float& LES, bool& forceCoeffs, float& Dref, float& omegaB, float& spzWidth)
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

    std::string restartStr("restart");
    restart = lookup<bool>(lines, restartStr);

    std::string FwriteStr("Fwrite");
    Fwrite = lookup<bool>(lines, FwriteStr);

    std::string writeBinaryStr("writeBinary");
    writeBinary = lookup<bool>(lines, writeBinaryStr);
    
    std::string startTimeStepStr("startTimeStep");
    startTimeStep = lookup<int>(lines, startTimeStepStr);
    
    std::string endTimeStepStr("endTimeStep");
    endTimeStep = lookup<int>(lines, endTimeStepStr);

    nextOutTime = startTimeStep;

    std::string outIntervalStr("outInterval");
    outInterval = lookup<int>(lines, outIntervalStr);
    std::cout << std::endl;

    std::string nxStr("nx");
    nx = lookup<int>(lines, nxStr);

    std::string nyStr("ny");
    ny = lookup<int>(lines, nyStr);

    std::string nzStr("nz");
    nz = lookup<int>(lines, nzStr);
    std::cout << std::endl;

    std::string LxStr("Lx");
    Lx = lookup<float>(lines, LxStr);
    std::cout << std::endl;

    std::string uMaxStr("uMax");//(m/s)
    uMax = lookup<float>(lines, uMaxStr);

    std::string rho0Str("rho0");
    rho0 = lookup<float>(lines, rho0Str);

    std::string U0Str("U0"); //Dimensionless maximum velocity
    U0 = lookup<float>(lines, U0Str);
    std::cout << std::endl;

    std::string nuStr("nu");
    nu = lookup<float>(lines, nuStr);
    std::cout << std::endl;

    std::string dpdxStr("dpdx");
    dpdx = lookup<float>(lines, dpdxStr);
    std::cout << std::endl;

    std::string LESStr("LES");
    LES = lookup<bool>(lines, LESStr) ? 1.f : 0.f;

    std::string forceCoeffsStr("forceCoeffs");
    forceCoeffs = lookup<bool>(lines, forceCoeffsStr);

    if(forceCoeffs)
    {
        std::string DrefStr("Dref");
        Dref = lookup<float>(lines, DrefStr);
    }
    else
    {
        Dref = 1.f;
    }

    std::string omegaBStr("omegaB");
    omegaB = lookupOrDefault<float>(lines, omegaBStr, 0.1f);

    std::string spzWidthStr("spzWidth");
    spzWidth = lookupOrDefault<float>(lines, spzWidthStr, 0.1f);

    inputFile.close();
}

int BCnameToNum(std::string name)
{
    if(name == "Cyclic")
    {
        return 0;
    }
    else if(name == "BounceBack")
    {
        return 1;
    }
    else if(name == "Outlet")
    {
        return 3;
    }
    else if(name == "Symmetry")
    {
        return 4;
    }
    else if(name == "FixedVelocity")
    {
        return 5;
    }
    else if(name == "FixedDensity")
    {
        return 6;
    }
    else
    {
        std::cerr << "Please select from Cyclic, BounceBack, Equilibrium, Symmetry, and Outlet" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void inputBoundaryCondition(const std::string bName, int& boundary1, int& boundary2, int& boundary3, float& U0, float& V0, float& W0, const float c)
{
    std::string inputFileName("boundaryConditions.txt");
    std::vector<std::string> lines;
    std::ifstream inputFile(inputFileName);
    if(!inputFile)
    {
        std::cerr << "Could not open input.txt" << std::endl;
        exit(EXIT_FAILURE);
    }

    readToLines(inputFile, lines);

    std::string boundary1Str(bName+"Boundary1");
    std::string boundary1name = lookup<std::string>(lines, boundary1Str);
    boundary1 = BCnameToNum(boundary1name);
    std::string boundary2Str(bName+"Boundary2");
    std::string boundary2name = lookup<std::string>(lines, boundary2Str);
    boundary2 = BCnameToNum(boundary2name);
    std::string boundary3Str(bName+"Boundary3");
    std::string boundary3name = lookup<std::string>(lines, boundary3Str);
    boundary3 = BCnameToNum(boundary3name);
    std::string U0Str(bName+"U0");
    U0 = lookup<float>(lines, U0Str)/c;
    std::string V0Str(bName+"V0");
    V0 = lookup<float>(lines, V0Str)/c;
    std::string W0Str(bName+"W0");
    W0 = lookup<float>(lines, W0Str)/c;
}

int inputSpzWidth(const std::string bName, const float L)
{
    std::string inputFileName("boundaryConditions.txt");
    std::vector<std::string> lines;
    std::ifstream inputFile(inputFileName);
    if(!inputFile)
    {
        std::cerr << "Could not open input.txt" << std::endl;
        exit(EXIT_FAILURE);
    }

    readToLines(inputFile, lines);

    std::string spzWidthStr(bName+"SpzWidth");
    float spzWidth = lookup<float>(lines, spzWidthStr);
    return int(spzWidth/L);
}

void inputBoundaryConditions(
    int& nxMinBoundary1, int& nxMinBoundary2, int& nxMinBoundary3, float& nxMinU0, float& nxMinV0, float& nxMinW0,
    int& nxMaxBoundary1, int& nxMaxBoundary2, int& nxMaxBoundary3, float& nxMaxU0, float& nxMaxV0, float& nxMaxW0,
    int& nyMinBoundary1, int& nyMinBoundary2, int& nyMinBoundary3, float& nyMinU0, float& nyMinV0, float& nyMinW0,
    int& nyMaxBoundary1, int& nyMaxBoundary2, int& nyMaxBoundary3, float& nyMaxU0, float& nyMaxV0, float& nyMaxW0,
    int& nzMinBoundary1, int& nzMinBoundary2, int& nzMinBoundary3, float& nzMinU0, float& nzMinV0, float& nzMinW0,
    int& nzMaxBoundary1, int& nzMaxBoundary2, int& nzMaxBoundary3, float& nzMaxU0, float& nzMaxV0, float& nzMaxW0,
    const float c
    )
{
    inputBoundaryCondition("nxMin",nxMinBoundary1,nxMinBoundary2,nxMinBoundary3,nxMinU0,nxMinV0,nxMinW0,c);
    inputBoundaryCondition("nxMax",nxMaxBoundary1,nxMaxBoundary2,nxMaxBoundary3,nxMaxU0,nxMaxV0,nxMaxW0,c);
    inputBoundaryCondition("nyMin",nyMinBoundary1,nyMinBoundary2,nyMinBoundary3,nyMinU0,nyMinV0,nyMinW0,c);
    inputBoundaryCondition("nyMax",nyMaxBoundary1,nyMaxBoundary2,nyMaxBoundary3,nyMaxU0,nyMaxV0,nyMaxW0,c);
    inputBoundaryCondition("nzMin",nzMinBoundary1,nzMinBoundary2,nzMinBoundary3,nzMinU0,nzMinV0,nzMinW0,c);
    inputBoundaryCondition("nzMax",nzMaxBoundary1,nzMaxBoundary2,nzMaxBoundary3,nzMaxU0,nzMaxV0,nzMaxW0,c);   
}

#endif
