#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include "sinusoidal_voltammetry.hpp"

//BOOST_PYTHON_FUNCTION_OVERLOADS(e_surface_overloads, e_surface, 3, 4)
void e_surface1(map& params, vector& Itot, vector& t) { 
    return e_surface(params,Itot,t);
}
void e_surface2(map& params, vector& Itot, vector& t, const vector* Edata) {
    return e_surface(params,Itot,t,Edata);
}

BOOST_PYTHON_MODULE(exp_design)
{
        using namespace boost::python;

        class_<vector>("vector")
            .def(vector_indexing_suite<vector>())
            ;

        class_<map>("map")
            .def(map_indexing_suite<map>())
            ;

        def("e_surface", e_surface1, (arg("params"), arg("Itot"), arg("t")), "Solves E problem for a surface confined reactant");
        def("e_surface", e_surface2,(arg("params"), arg("Itot"), arg("t"),arg("Edata")), "Solves E problem for a surface confined reactant");
}
