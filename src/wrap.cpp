#include <boost/python.hpp>

#include "simulator.cpp"

using namespace boost::python;

BOOST_PYTHON_MODULE( source_VGsim ) {
    class_<Simulator>( "Simulator" )
        .def( init<uint64_t,uint64_t,uint64_t,uint64_t>( args( "number_of_sites", "number_of_populations", "number_of_susceptible_groups", "seed" ) ) )
        .def( "simulate", static_cast< void (Simulator::*)(uint64_t, uint64_t, double, std::string, uint64_t) >( &Simulator::simulate ), args("iterations", "sampling", "epidemic_time", "type", "number_attempts"))
        .def( "Genealogy", static_cast< void (Simulator::*)() > ( &Simulator::Genealogy ))
        .def( "Debug", static_cast< void (Simulator::*)() > ( &Simulator::Debug ))

        .def( "get_flat_chain", static_cast< PyObject* (Simulator::*)() > ( &Simulator::get_flat_chain ))

        // Infectious
        .def( "set_susceptibility_group", static_cast< void (Simulator::*)(uint64_t, uint64_t) > ( &Simulator::set_susceptibility_group ), args("group", "haplotype"))
        .def( "get_susceptibility_group", static_cast< PyObject* (Simulator::*)() > ( &Simulator::get_susceptibility_group ))
        .def( "set_transmission_rate", static_cast< void (Simulator::*)(double, uint64_t) > ( &Simulator::set_transmission_rate ), args("rate", "haplotype"))
        .def( "get_transmission_rate", static_cast< PyObject* (Simulator::*)() > ( &Simulator::get_transmission_rate ))
        .def( "set_recovery_rate", static_cast< void (Simulator::*)(double, uint64_t) > ( &Simulator::set_recovery_rate ), args("rate", "haplotype"))
        .def( "get_recovery_rate", static_cast< PyObject* (Simulator::*)() > ( &Simulator::get_recovery_rate ))
        .def( "set_sampling_rate", static_cast< void (Simulator::*)(double, uint64_t) > ( &Simulator::set_sampling_rate ), args("rate", "haplotype"))
        .def( "get_sampling_rate", static_cast< PyObject* (Simulator::*)() > ( &Simulator::get_sampling_rate ))
        .def( "set_mutation_rate", static_cast< void (Simulator::*)(double, uint64_t, uint64_t) > ( &Simulator::set_mutation_rate ), args("rate", "haplotype", "mutation"))
        .def( "get_mutation_rate", static_cast< PyObject* (Simulator::*)() > ( &Simulator::get_mutation_rate ))
        .def( "set_mutation_probabilities", static_cast< void (Simulator::*)(double, uint64_t, uint64_t, uint64_t) > ( &Simulator::set_mutation_probabilities ), args("rate", "haplotype", "mutation", "index"))
        .def( "get_mutation_probabilities", static_cast< PyObject* (Simulator::*)() > ( &Simulator::get_mutation_probabilities ))
        .def( "set_susceptibility", static_cast< void (Simulator::*)(double, uint64_t, uint64_t) > ( &Simulator::set_susceptibility ), args("rate", "haplotype", "group"))
        .def( "get_susceptibility", static_cast< PyObject* (Simulator::*)() > ( &Simulator::get_susceptibility ))
    ;
}