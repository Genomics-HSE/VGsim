#include <boost/python.hpp>

#include "simulator.cpp"

using namespace boost::python;

BOOST_PYTHON_MODULE( VGsim_test ) {
    class_<Simulator>( "Simulator" )
        .def( init<uint64_t,uint64_t,uint64_t,uint64_t,uint64_t>( args( "number_of_sites", "number_of_populations", "number_of_susceptible_groups", "seed", "number_attempts" ) ) )
        .def( "DirectSimulate", static_cast< void (Simulator::*)(uint64_t) >( &Simulator::DirectSimulate ), args( "iterations" ) )
    ;
}