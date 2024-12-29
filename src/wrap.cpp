#include <boost/python.hpp>

#include "simulator.cpp"

using namespace boost::python;

BOOST_PYTHON_MODULE( source_VGsim ) {
    class_<Simulator>( "Simulator" )
        .def( init<uint64_t,uint64_t,uint64_t,uint64_t>( args( "number_of_sites", "number_of_populations", "number_of_susceptible_groups", "seed" ) ) )
        // .def( init<uint64_t,uint64_t,uint64_t,uint64_t>( (arg("number_of_sites")=0, arg("number_of_populations")=1, arg("number_of_susceptible_groups")=1, arg("seed")=1234) ) )
        // .def( "Simulate", static_cast< void (Simulator::*)(uint64_t, std::string, uint64_t) >( &Simulator::Simulate ), args( "iterations", "type", "number_attempts" ) )
        .def( "simulate", static_cast< void (Simulator::*)(uint64_t, std::string, uint64_t) >( &Simulator::simulate ), (arg("iterations")=1000000, arg("type")="direct", arg("number_attempts")=100) )
        .def( "Genealogy", static_cast< void (Simulator::*)() > ( &Simulator::Genealogy ))
        .def( "Debug", static_cast< void (Simulator::*)() > ( &Simulator::Debug ))

        .def( "get_flat_chain", static_cast< PyObject* (Simulator::*)() > ( &Simulator::get_flat_chain ))

        .def( "set_transmission_rate", static_cast< void (Simulator::*)(double, uint64_t) > ( &Simulator::set_transmission_rate ), args("rate", "haplotype"))
        .def( "get_transmission_rate", static_cast< PyObject* (Simulator::*)() > ( &Simulator::get_transmission_rate ))
    ;
}