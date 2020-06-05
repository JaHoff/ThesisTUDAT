/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


/* Optimize the configuration of a small satellite swarm in orbit around the L4 point
Base adapted from the propagationTargetingExample included by TUDAT*/

#include <pagmo/problem.hpp>
#include <pagmo/algorithms/sade.hpp>
#include <pagmo/algorithms/de1220.hpp>
#include <pagmo/algorithms/de.hpp>
#include <pagmo/algorithms/simulated_annealing.hpp>
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>

#include "Problems/SwarmOptimization.h"
#include "Problems/applicationOutput.h"
#include "Problems/saveOptimizationResults.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"

using namespace pagmo;
using namespace swarm_optimization;
using namespace tudat;

int main( )
{
    std::cout << "General optimization loop start!" << std::endl;
    //Set seed for reproducible results
    pagmo::random_device::set_seed(42);

    //Load spice kernels
    tudat::spice_interface::loadStandardSpiceKernels( );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
//    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
//                                          propagators::body_fixed_relative_cartesian_position, "Moon", "Earth" ) );
//    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
//                                          propagators::relative_velocity_dependent_variable, "Moon", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared < propagators::SingleDependentVariableSaveSettings> (
                                          propagators::keplerian_state_dependent_variable, "Moon", "Earth"));

    // Create object with list of dependent variables.
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    std::cout << "Start defining the general optimization problem" << std::endl;
    int numberOfSatellites = 5;
    double missionLength = 160*tudat::physical_constants::JULIAN_DAY;


    SwarmOptimization swarmProblem( numberOfSatellites, dependentVariablesToSave, missionLength );
    std::cout << "Problemize the problem" << std::endl;
    problem prob{ swarmProblem };

    // Instantiate a pagmo algorithm
    algorithm algo{de1220( )};

    // Create an island with 128 individuals
    pagmo::population::size_type populationSize = 128;
    island isl{algo, prob, populationSize};

    std::map <int, std::vector< double >> fitnessmap;
    std::cout << "Starting evolving optimization problem for 25 generations!" << std::endl;
    // Evolve for 25 generations
    for( int i = 0; i < 75; i++ )
    {
        std::cout << "Starting generation " + std::to_string(i) << std::endl;
        isl.evolve( );
        while( isl.status( ) != pagmo::evolve_status::idle &&
               isl.status( ) != pagmo::evolve_status::idle_error )
        {
            isl.wait( );
        }
        isl.wait_check( ); // Raises errors

        fitnessmap.insert( std::pair<double, std::vector< double > >( i, isl.get_population().get_f().at(0) ) );
        // Write current iteration results to file
        printPopulationToFile( isl.get_population( ).get_x( ), "swarmPropagation_" + std::to_string( i ) + "_" + std::to_string( i ) , false );
        printPopulationToFile( isl.get_population( ).get_f( ), "swarmPropagation_" + std::to_string( i ) + "_" + std::to_string( i ) , true );
    }

    std::cout << "Best solution found was for cost: " << swarmProblem.getBestCost() << std::endl;
    // Retrieve final Cartesian states for population in last generation, and save final states to a file.
    std::vector<std::vector< double > > decisionVariables = isl.get_population( ).get_x( );
    std::map< int, Eigen::VectorXd > finalStates;
    for( unsigned int i = 0; i < decisionVariables.size( ); i++ )
    {
        swarmProblem.fitness( decisionVariables.at( i ) );
        finalStates[ i ] = swarmProblem.getPreviousFinalState( );
    }
    tudat::input_output::writeDataMapToTextFile(
                finalStates, "swarmFinalStates.dat", swarm_optimization::getOutputPath( ) );

    // Write lunar state history to file.
    input_output::writeDataMapToTextFile( swarmProblem.ComputeLunarOrbit(),
                                          "propagationHistory_moon.dat",
                                          swarm_optimization::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write core position to file.
    input_output::writeMatrixToFile(swarmProblem.getCorePosition(),
                                    "corePosition_best.dat",10,
                                    swarm_optimization::getOutputPath( ) ) ;

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( swarmProblem.getBestStateHistory(),
                                          "propagationHistory_best.dat",
                                          swarm_optimization::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );


//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    //////////////// PERTURBED PROBLEM OPTIMIZATION ///////////////////////////////////////////////////////////////////////////////
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    // Create object to compute the problem fitness; with perturbations
//    problem prob_pert{SwarmOptimization( numberOfSatellites, dependentVariablesToSave, missionLength) };


//    // Instantiate a pagmo algorithm for the new problem.
//    algorithm algo_pert{de1220( )};

//    // Create an empty population for perturbed problem
//    population population_pert = population( prob_pert, 0 );

//    // Retrieve population of unperturbed problem, and instantiate population of perturbed problem
//    std::vector<vector_double> original_population = isl.get_population( ).get_x( );
//    for( unsigned int k = 0; k < populationSize; k++ )
//    {
//        population_pert.push_back( original_population.at( k ) );
//    }

//    // Create island for perturbed problem
//    island isl_pert{algo_pert, population_pert};

//    // Write original (unevolved) population to file
//    printPopulationToFile( isl_pert.get_population( ).get_x( ), "swarmPropagation_pert_orig" , false );
//    printPopulationToFile( isl_pert.get_population( ).get_f( ), "swarmPropagation_pert_orig" , true );


//    // Evolve for 4 generations
//    for( int i = 0; i < 4; i++ )
//    {
//        isl_pert.evolve( );
//        while( isl_pert.status( ) != pagmo::evolve_status::idle &&
//               isl_pert.status( ) != pagmo::evolve_status::idle_error )
//        {
//            isl_pert.wait( );
//        }
//        isl_pert.wait_check( ); // Raises errors

//        // Write current iteration results to file
//        printPopulationToFile( isl_pert.get_population( ).get_x( ), "swarmPropagation_pert_" + std::to_string( i ) + "_" + std::to_string( i ) , false );
//        printPopulationToFile( isl_pert.get_population( ).get_f( ), "swarmPropagation_pert_" + std::to_string( i ) + "_" + std::to_string( i ) , true );

//        std::cout<<i<<std::endl;
//    }
}

