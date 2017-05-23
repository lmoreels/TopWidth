#ifndef INTERFACE_HELPERTOOLS_H
#define INTERFACE_HELPERTOOLS_H

#include <stdio.h>
#include <sys/stat.h>
#include <errno.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <map>

class HelperTools{
  public:
    HelperTools();
    ~HelperTools();
    bool fexists(const char *filename);
    std::string ConvertDoubleToString(double Number);
    std::string DotReplace(double var);
    std::string ConvertIntToString(int Number);
    std::string ConvertIntToString(int Number, int pad);
    std::string MakeTimeStamp();
};


#endif
