#ifndef input_H
#define input_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <omp.h>

template <typename Type = std::string>
Type returnWrapper(std::string arg);

template <>
int returnWrapper<int>(std::string arg);

template <>
bool returnWrapper<bool>(std::string arg);

template <>
float returnWrapper<float>(std::string arg);

template <>
double returnWrapper<double>(std::string arg);

template<typename Type> Type lookup(std::vector<std::string>& lines, std::string& str);

void readToLines(std::ifstream& inputFile, std::vector<std::string>& lines);

void input(bool& restart, bool& Fwrite, bool& writeBinary, int& startTimeStep, int& endTimeStep, int& nextOutTime, int& outInterval, int& nx, int& ny, int& nz);

#endif