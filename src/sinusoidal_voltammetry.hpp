#ifndef SINUSOIDAL_VOLTAMMETRY_HPP
#define SINUSOIDAL_VOLTAMMETRY_HPP

#include <vector>
#include <map>
#include <string>


typedef std::vector<double> vector; 
typedef std::map<std::string,double> map; 

#include "utilities.hpp"

void e_surface(map& params, vector& Itot, vector& t, const vector* Edata=NULL);

#endif //SINUSOIDAL_VOLTAMMETRY_HPP
