#include <boost/python.hpp>

#include "simulator.cpp"

using namespace boost::python;

BOOST_PYTHON_MODULE( VGsim ) {
    class_<Simulator>( "Simulator" )
        .def( init<uint64_t,uint64_t,uint64_t,uint64_t>( args( "number_of_sites", "number_of_populations", "number_of_susceptible_groups", "seed" ) ) )
        .def( "Simulate", static_cast< void (Simulator::*)(uint64_t, std::string, uint64_t) >( &Simulator::Simulate ), args( "iterations", "type", "number_attempts" ) )
        .def( "Genealogy", static_cast< void (Simulator::*)() > ( &Simulator::Genealogy ))
        .def( "Debug", static_cast< void (Simulator::*)() > ( &Simulator::Debug ))
    ;
}