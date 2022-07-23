#ifndef input_H
#define input_H

#include <vector>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <omp.h>

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




void input(bool& restart, bool& Fwrite, bool& writeBinary, int& startTimeStep, int& endTimeStep, int& nextOutTime, int& outInterval, int& nx, int& ny, int& nz, float& Lx, float& uMax, float& rho0, float& U0, float& nu, float& dpdx, float& LES, bool& forceCoeffs, float& Dref)
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

    std::string nThreadsStr("nThreads");
    omp_set_num_threads(lookup<int>(lines, nThreadsStr));

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
    else if(name == "Equilibrium")
    {
        return 2;
    }
    else if(name == "Outlet")
    {
        return 3;
    }
    else
    {
        std::cerr << "Please select from Cyclic, BounceBack, Equilibrium, and Outlet" << std::endl;
        exit(EXIT_FAILURE);
    }
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
    std::string inputFileName("boundaryConditions.txt");
    std::vector<std::string> lines;
    std::ifstream inputFile(inputFileName);
    if(!inputFile)
    {
        std::cerr << "Could not open input.txt" << std::endl;
        exit(EXIT_FAILURE);
    }

    readToLines(inputFile, lines);

    std::string nxMinBoundary1Str("nxMinBoundary1");
    std::string nxMinBoundary1name = lookup<std::string>(lines, nxMinBoundary1Str);
    nxMinBoundary1 = BCnameToNum(nxMinBoundary1name);
    std::string nxMinBoundary2Str("nxMinBoundary2");
    std::string nxMinBoundary2name = lookup<std::string>(lines, nxMinBoundary2Str);
    nxMinBoundary2 = BCnameToNum(nxMinBoundary2name);
    std::string nxMinBoundary3Str("nxMinBoundary3");
    std::string nxMinBoundary3name = lookup<std::string>(lines, nxMinBoundary3Str);
    nxMinBoundary3 = BCnameToNum(nxMinBoundary3name);
    std::string nxMinU0Str("nxMinU0");
    nxMinU0 = lookup<float>(lines, nxMinU0Str)/c;
    std::string nxMinV0Str("nxMinV0");
    nxMinV0 = lookup<float>(lines, nxMinV0Str)/c;
    std::string nxMinW0Str("nxMinW0");
    nxMinW0 = lookup<float>(lines, nxMinW0Str)/c;

    std::string nxMaxBoundary1Str("nxMaxBoundary1");
    std::string nxMaxBoundary1name = lookup<std::string>(lines, nxMaxBoundary1Str);
    nxMaxBoundary1 = BCnameToNum(nxMaxBoundary1name);
    std::string nxMaxBoundary2Str("nxMaxBoundary2");
    std::string nxMaxBoundary2name = lookup<std::string>(lines, nxMaxBoundary2Str);
    nxMaxBoundary2 = BCnameToNum(nxMaxBoundary2name);
    std::string nxMaxBoundary3Str("nxMaxBoundary3");
    std::string nxMaxBoundary3name = lookup<std::string>(lines, nxMaxBoundary3Str);
    nxMaxBoundary3 = BCnameToNum(nxMaxBoundary3name);
    std::string nxMaxU0Str("nxMaxU0");
    nxMaxU0 = lookup<float>(lines, nxMaxU0Str)/c;
    std::string nxMaxV0Str("nxMaxV0");
    nxMaxV0 = lookup<float>(lines, nxMaxV0Str)/c;
    std::string nxMaxW0Str("nxMaxW0");
    nxMaxW0 = lookup<float>(lines, nxMaxW0Str)/c;


    std::string nyMinBoundary1Str("nyMinBoundary1");
    std::string nyMinBoundary1name = lookup<std::string>(lines, nyMinBoundary1Str);
    nyMinBoundary1 = BCnameToNum(nyMinBoundary1name);
    std::string nyMinBoundary2Str("nyMinBoundary2");
    std::string nyMinBoundary2name = lookup<std::string>(lines, nyMinBoundary2Str);
    nyMinBoundary2 = BCnameToNum(nyMinBoundary2name);
    std::string nyMinBoundary3Str("nyMinBoundary3");
    std::string nyMinBoundary3name = lookup<std::string>(lines, nyMinBoundary3Str);
    nyMinBoundary3 = BCnameToNum(nyMinBoundary3name);
    std::string nyMinU0Str("nyMinU0");
    nyMinU0 = lookup<float>(lines, nyMinU0Str)/c;
    std::string nyMinV0Str("nyMinV0");
    nyMinV0 = lookup<float>(lines, nyMinV0Str)/c;
    std::string nyMinW0Str("nyMinW0");
    nyMinW0 = lookup<float>(lines, nyMinW0Str)/c;

    std::string nyMaxBoundary1Str("nyMaxBoundary1");
    std::string nyMaxBoundary1name = lookup<std::string>(lines, nyMaxBoundary1Str);
    nyMaxBoundary1 = BCnameToNum(nyMaxBoundary1name);
    std::string nyMaxBoundary2Str("nyMaxBoundary2");
    std::string nyMaxBoundary2name = lookup<std::string>(lines, nyMaxBoundary2Str);
    nyMaxBoundary2 = BCnameToNum(nyMaxBoundary2name);
    std::string nyMaxBoundary3Str("nyMaxBoundary3");
    std::string nyMaxBoundary3name = lookup<std::string>(lines, nyMaxBoundary3Str);
    nyMaxBoundary3 = BCnameToNum(nyMaxBoundary3name);
    std::string nyMaxU0Str("nyMaxU0");
    nyMaxU0 = lookup<float>(lines, nyMaxU0Str)/c;
    std::string nyMaxV0Str("nyMaxV0");
    nyMaxV0 = lookup<float>(lines, nyMaxV0Str)/c;
    std::string nyMaxW0Str("nyMaxW0");
    nyMaxW0 = lookup<float>(lines, nyMaxW0Str)/c;


    std::string nzMinBoundary1Str("nzMinBoundary1");
    std::string nzMinBoundary1name = lookup<std::string>(lines, nzMinBoundary1Str);
    nzMinBoundary1 = BCnameToNum(nzMinBoundary1name);
    std::string nzMinBoundary2Str("nzMinBoundary2");
    std::string nzMinBoundary2name = lookup<std::string>(lines, nzMinBoundary2Str);
    nzMinBoundary2 = BCnameToNum(nzMinBoundary2name);
    std::string nzMinBoundary3Str("nzMinBoundary3");
    std::string nzMinBoundary3name = lookup<std::string>(lines, nzMinBoundary3Str);
    nzMinBoundary3 = BCnameToNum(nzMinBoundary3name);
    std::string nzMinU0Str("nzMinU0");
    nzMinU0 = lookup<float>(lines, nzMinU0Str)/c;
    std::string nzMinV0Str("nzMinV0");
    nzMinV0 = lookup<float>(lines, nzMinV0Str)/c;
    std::string nzMinW0Str("nzMinW0");
    nzMinW0 = lookup<float>(lines, nzMinW0Str)/c;

    std::string nzMaxBoundary1Str("nzMaxBoundary1");
    std::string nzMaxBoundary1name = lookup<std::string>(lines, nzMaxBoundary1Str);
    nzMaxBoundary1 = BCnameToNum(nzMaxBoundary1name);
    std::string nzMaxBoundary2Str("nzMaxBoundary2");
    std::string nzMaxBoundary2name = lookup<std::string>(lines, nzMaxBoundary2Str);
    nzMaxBoundary2 = BCnameToNum(nzMaxBoundary2name);
    std::string nzMaxBoundary3Str("nzMaxBoundary3");
    std::string nzMaxBoundary3name = lookup<std::string>(lines, nzMaxBoundary3Str);
    nzMaxBoundary3 = BCnameToNum(nzMaxBoundary3name);
    std::string nzMaxU0Str("nzMaxU0");
    nzMaxU0 = lookup<float>(lines, nzMaxU0Str)/c;
    std::string nzMaxV0Str("nzMaxV0");
    nzMaxV0 = lookup<float>(lines, nzMaxV0Str)/c;
    std::string nzMaxW0Str("nzMaxW0");
    nzMaxW0 = lookup<float>(lines, nzMaxW0Str)/c;
    
}

#endif
