#include <boost/python.hpp>

#include "simulator.cpp"

using namespace boost::python;

BOOST_PYTHON_MODULE( VGsim_new ) {
    class_<Simulator>( "Simulator" )
        .def( init<uint64_t,uint64_t,uint64_t,uint64_t>( args( "number_of_sites", "number_of_populations", "number_of_susceptible_groups", "seed" ) ) )
        // .def( init<uint64_t,uint64_t,uint64_t,uint64_t>( (arg("number_of_sites")=0, arg("number_of_populations")=1, arg("number_of_susceptible_groups")=1, arg("seed")=1234) ) )
        // .def( "Simulate", static_cast< void (Simulator::*)(uint64_t, std::string, uint64_t) >( &Simulator::Simulate ), args( "iterations", "type", "number_attempts" ) )
        .def( "Simulate", static_cast< void (Simulator::*)(uint64_t, std::string, uint64_t) >( &Simulator::Simulate ), (arg("iterations")=1000000, arg("type")="direct", arg("number_attempts")=100) )
        .def( "Genealogy", static_cast< void (Simulator::*)() > ( &Simulator::Genealogy ))
        .def( "Debug", static_cast< void (Simulator::*)() > ( &Simulator::Debug ))
    ;
}