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

#include <chrono>

#include <pagmo/problem.hpp>
#include <pagmo/algorithms/sade.hpp>
#include <pagmo/algorithms/de1220.hpp>
#include <pagmo/algorithms/de.hpp>
#include <pagmo/algorithms/simulated_annealing.hpp>
#include <pagmo/algorithms/gaco.hpp>
#include <pagmo/algorithms/pso.hpp>
#include <pagmo/algorithms/pso_gen.hpp>
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>

#include <pagmo/bfe.hpp>
#include <pagmo/batch_evaluators/default_bfe.hpp>

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
    int n_generations = 75;
    int n_pops = 128;
    int r_seed = 42;
    string subfolder = "/AlgComparison/";
    std::cout << "General optimization start!" << std::endl;

    std::vector<std::string> algo_list_names{"Differential Evolution", "Self Adjusting Differential Evolution",
                                            "Particle Swarm Optimization",
                                            "Particle Swarm Optimization Generational",
                                             "Generational Ant Colony"};
    std::vector<std::string> algo_names_shorthand{"de1220", "sade", "pso", "pso_gen", "gaco"};


    //Set seed for reproducible results
    pagmo::random_device::set_seed(r_seed);

    //Load spice kernels
    tudat::spice_interface::loadStandardSpiceKernels( );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
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

    for (int g = 0; g < 5; g++){
        auto t1 = std::chrono::high_resolution_clock::now();
        string namesnip = algo_names_shorthand[g] + "_sd" + std::to_string(r_seed);

        // reset the internal tracker for the best cost
        swarmProblem.resetBestCost();

        // Instantiate a pagmo algorithm
        algorithm algo;
        switch(g){
        case 0:{
            algo = de1220();
            break;
        }
        case 1:{
            algo = sade();
            break;
        }
        case 2:{
            algo = pso();
            break;
        }
        case 3:{
            algo = pso_gen();
            break;
        }
        case 4:{
            algo = gaco();
            break;
        }
        }


        // Create an island with 128 individuals
        pagmo::population::size_type populationSize = n_pops;


        island isl{algo, prob, populationSize};

        std::cout << "Starting evolving optimization problem for "<< n_generations << " generations!" << std::endl;
        int i = 0;
        bool iterate = true;
        while( iterate)
        {
            std::cout << "Now at generation " + std::to_string(i) << std::endl;
            i++;
            isl.evolve( );
            while( isl.status( ) != pagmo::evolve_status::idle &&
                   isl.status( ) != pagmo::evolve_status::idle_error )
            {
                isl.wait( );
            }
            isl.wait_check( ); // Raises errors


            std::vector<vector_double> popsf = isl.get_population().get_f();
            printPopulationToFile( isl.get_population( ).get_x( ), namesnip+"_" + std::to_string( i ) + "_" + std::to_string( i ) , false, subfolder );
            printPopulationToFile( popsf, namesnip+"_" + std::to_string( i ) + "_" + std::to_string( i ) , true, subfolder );

            std::cout << "printed file: " << namesnip+"_" + std::to_string( i ) + "_" + std::to_string( i ) << std::endl;
            population pops = isl.get_population();
            double bestcost = pops.get_f()[pops.best_idx()][0];

            if (i > n_generations || bestcost == 0){ iterate = false;}
        }


        auto t2 = std::chrono::high_resolution_clock::now();
        std::cout << "Full optimization of " << n_generations << " generations took " <<
                     std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count()<< " seconds using " <<
                     algo_list_names[g] << std::endl;


        // Retrieve final Cartesian states for population in last generation, and save final states to a file.
        std::vector<std::vector< double > > decisionVariables = isl.get_population( ).get_x( );
        std::map< int, Eigen::VectorXd > finalStates;
        for( unsigned int i = 0; i < decisionVariables.size( ); i++ )
        {
            swarmProblem.fitness( decisionVariables.at( i ) );
            finalStates[ i ] = swarmProblem.getPreviousFinalState( );
        }
        tudat::input_output::writeDataMapToTextFile(
                    finalStates, "swarmFinalStates_"+namesnip+".dat", swarm_optimization::getOutputPath( ) + subfolder );

        // Write lunar state history to file.
        input_output::writeDataMapToTextFile( swarmProblem.ComputeLunarOrbit(),
                                              "propagationHistory_moon.dat",
                                              swarm_optimization::getOutputPath( ) + subfolder,
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        // Write core position to file.
        input_output::writeMatrixToFile(swarmProblem.getCorePosition(),
                                        "corePosition_"+namesnip+"_best.dat",10,
                                        swarm_optimization::getOutputPath( ) + subfolder ) ;

        // Write perturbed satellite propagation history to file.
        input_output::writeDataMapToTextFile( swarmProblem.getBestStateHistory(),
                                              "propagationHistory_"+namesnip+"_best.dat",
                                              swarm_optimization::getOutputPath( ) + subfolder,
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        input_output::writeDataMapToTextFile(swarmProblem.getPenalizedBaselineHistoryMap(),
                                              "penaltyHistory_"+namesnip+".dat",
                                              swarm_optimization::getOutputPath( ) + subfolder,
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              ",");
    }


}

