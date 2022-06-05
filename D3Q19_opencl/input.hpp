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

    if(!inputFile.is_open())
    {
        std::cerr << "Could not open the file" << std::endl;
        exit(EXIT_FAILURE);
    }


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




void input(bool& restart, bool& Fwrite, bool& writeBinary, int& startTimeStep, int& endTimeStep, int& nextOutTime, int& outInterval, int& nx, int& ny, int& nz, float& uMax, float& rho0, float& Re, float& U0)
{
    std::string inputFileName("input.txt");
    std::vector<std::string> lines;
    std::ifstream inputFile(inputFileName);

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

    std::string uMaxStr("uMax");//(m/s)
    uMax = lookup<float>(lines, uMaxStr);

    std::string rho0Str("rho0");
    rho0 = lookup<float>(lines, rho0Str);

    std::string ReStr("Re");
    Re = lookup<float>(lines, ReStr);
    std::cout << std::endl;

    std::string U0Str("U0"); //Dimensionless maximum velocity
    U0 = lookup<float>(lines, U0Str);
    std::cout << std::endl;

    inputFile.close();
}

void readSTL(const std::string STLname, std::vector<std::vector<float> >& STLnormal, std::vector<std::vector<float> >& STLv0, std::vector<std::vector<float> >& STLv1, std::vector<std::vector<float> >& STLv2)
{
    std::ifstream STLfile(STLname);
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

}
#endif
