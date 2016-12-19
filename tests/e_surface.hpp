#include <cxxtest/TestSuite.h>
#include "sinusoidal_voltammetry.hpp"
#include <math.h>
#include <fstream>

class Test_e_surface: public CxxTest::TestSuite {
public:
    void testDefault(void) {
        map params;

        vector Itot;
        vector t;

        e_surface(params,Itot,t);

    }
};
