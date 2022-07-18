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




void input(bool& restart, bool& Fwrite, bool& writeBinary, int& startTimeStep, int& endTimeStep, int& nextOutTime, int& outInterval, int& nx, int& ny, int& nz, float& Lx, float& uMax, float& rho0, float& U0, float& nu, float& dpdx)
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

    inputFile.close();
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
    nxMinBoundary1 = lookup<int>(lines, nxMinBoundary1Str);
    std::string nxMinBoundary2Str("nxMinBoundary2");
    nxMinBoundary2 = lookup<int>(lines, nxMinBoundary2Str);
    std::string nxMinBoundary3Str("nxMinBoundary3");
    nxMinBoundary3 = lookup<int>(lines, nxMinBoundary3Str);
    std::string nxMinU0Str("nxMinU0");
    nxMinU0 = lookup<float>(lines, nxMinU0Str)/c;
    std::string nxMinV0Str("nxMinV0");
    nxMinV0 = lookup<float>(lines, nxMinV0Str)/c;
    std::string nxMinW0Str("nxMinW0");
    nxMinW0 = lookup<float>(lines, nxMinW0Str)/c;

    std::string nxMaxBoundary1Str("nxMaxBoundary1");
    nxMaxBoundary1 = lookup<int>(lines, nxMaxBoundary1Str);
    std::string nxMaxBoundary2Str("nxMaxBoundary2");
    nxMaxBoundary2 = lookup<int>(lines, nxMaxBoundary2Str);
    std::string nxMaxBoundary3Str("nxMaxBoundary3");
    nxMaxBoundary3 = lookup<int>(lines, nxMaxBoundary3Str);
    std::string nxMaxU0Str("nxMaxU0");
    nxMaxU0 = lookup<float>(lines, nxMaxU0Str)/c;
    std::string nxMaxV0Str("nxMaxV0");
    nxMaxV0 = lookup<float>(lines, nxMaxV0Str)/c;
    std::string nxMaxW0Str("nxMaxW0");
    nxMaxW0 = lookup<float>(lines, nxMaxW0Str)/c;


    std::string nyMinBoundary1Str("nyMinBoundary1");
    nyMinBoundary1 = lookup<int>(lines, nyMinBoundary1Str);
    std::string nyMinBoundary2Str("nyMinBoundary2");
    nyMinBoundary2 = lookup<int>(lines, nyMinBoundary2Str);
    std::string nyMinBoundary3Str("nyMinBoundary3");
    nyMinBoundary3 = lookup<int>(lines, nyMinBoundary3Str);
    std::string nyMinU0Str("nyMinU0");
    nyMinU0 = lookup<float>(lines, nyMinU0Str)/c;
    std::string nyMinV0Str("nyMinV0");
    nyMinV0 = lookup<float>(lines, nyMinV0Str)/c;
    std::string nyMinW0Str("nyMinW0");
    nyMinW0 = lookup<float>(lines, nyMinW0Str)/c;

    std::string nyMaxBoundary1Str("nyMaxBoundary1");
    nyMaxBoundary1 = lookup<int>(lines, nyMaxBoundary1Str);
    std::string nyMaxBoundary2Str("nyMaxBoundary2");
    nyMaxBoundary2 = lookup<int>(lines, nyMaxBoundary2Str);
    std::string nyMaxBoundary3Str("nyMaxBoundary3");
    nyMaxBoundary3 = lookup<int>(lines, nyMaxBoundary3Str);
    std::string nyMaxU0Str("nyMaxU0");
    nyMaxU0 = lookup<float>(lines, nyMaxU0Str)/c;
    std::string nyMaxV0Str("nyMaxV0");
    nyMaxV0 = lookup<float>(lines, nyMaxV0Str)/c;
    std::string nyMaxW0Str("nyMaxW0");
    nyMaxW0 = lookup<float>(lines, nyMaxW0Str)/c;


    std::string nzMinBoundary1Str("nzMinBoundary1");
    nzMinBoundary1 = lookup<int>(lines, nzMinBoundary1Str);
    std::string nzMinBoundary2Str("nzMinBoundary2");
    nzMinBoundary2 = lookup<int>(lines, nzMinBoundary2Str);
    std::string nzMinBoundary3Str("nzMinBoundary3");
    nzMinBoundary3 = lookup<int>(lines, nzMinBoundary3Str);
    std::string nzMinU0Str("nzMinU0");
    nzMinU0 = lookup<float>(lines, nzMinU0Str)/c;
    std::string nzMinV0Str("nzMinV0");
    nzMinV0 = lookup<float>(lines, nzMinV0Str)/c;
    std::string nzMinW0Str("nzMinW0");
    nzMinW0 = lookup<float>(lines, nzMinW0Str)/c;

    std::string nzMaxBoundary1Str("nzMaxBoundary1");
    nzMaxBoundary1 = lookup<int>(lines, nzMaxBoundary1Str);
    std::string nzMaxBoundary2Str("nzMaxBoundary2");
    nzMaxBoundary2 = lookup<int>(lines, nzMaxBoundary2Str);
    std::string nzMaxBoundary3Str("nzMaxBoundary3");
    nzMaxBoundary3 = lookup<int>(lines, nzMaxBoundary3Str);
    std::string nzMaxU0Str("nzMaxU0");
    nzMaxU0 = lookup<float>(lines, nzMaxU0Str)/c;
    std::string nzMaxV0Str("nzMaxV0");
    nzMaxV0 = lookup<float>(lines, nzMaxV0Str)/c;
    std::string nzMaxW0Str("nzMaxW0");
    nzMaxW0 = lookup<float>(lines, nzMaxW0Str)/c;
}

#endif
