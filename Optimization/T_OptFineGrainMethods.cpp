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
Base adapted from the propagationTargetingExample included by TUDAT

Compare the performance of fine-grain optimisation methods using batch fitness evaluation.

DEPRECATED - Not supported in TUDAT*/

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

//#include <pagmo/bfe.hpp>
//#include <pagmo/batch_evaluators/default_bfe.hpp>

#include "Problems/SwarmOptimization.h"
#include "Problems/applicationOutput.h"
#include "Problems/saveOptimizationResults.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"

using namespace pagmo;
using namespace swarm_optimization;
using namespace tudat;

#define n_islands 2;

int main( )
{
    int n_generations = 5;
    int r_seed = 42;
    std::cout << "General optimization start!" << std::endl;

    std::vector<std::string> algo_list_names{"Particle Swarm Optimization Generational BFE",
                                             "Generational Ant Colony BFE"};
    std::vector<std::string> algo_names_shorthand{ "pso_genbfe","gacobfe"};


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

    for (int g = 0; g < 1; g++){
        auto t1 = std::chrono::high_resolution_clock::now();
        string namesnip = algo_names_shorthand[g] + "_sd" + std::to_string(r_seed);

        // reset the internal tracker for the best cost
        swarmProblem.resetBestCost();

        // Instantiate a pagmo algorithm
        algorithm algo;
        switch(g){
        case 0:{
            pso_gen psog = pso_gen();
            psog.set_bfe(bfe());
            algo = psog;
            break;
        }
        case 1:{
            gaco genalg = gaco();
            genalg.set_bfe(bfe());
            algo = genalg;
            break;
        }
        default:{
            std::cout << "invalid case!"<< std::endl;
        }
        }

        // Create an island with 128 individuals
        pagmo::population::size_type populationSize = 32;

        population pops = population(swarmProblem, populationSize);


        std::cout << "Starting evolving optimization problem for "<< n_generations << " generations!" << std::endl;
        int i = 0;
        bool iterate = true;
        while( iterate)
        {
            std::cout << "Now at generation " + std::to_string(i) << std::endl;
            i++;

            algo.evolve(pops);


            std::vector<vector_double> popsf = pops.get_f();
            printPopulationToFile( pops.get_x( ), "swarmPropagation_"+namesnip+"_" + std::to_string( i ) + "_" + std::to_string( i ) , false );
            printPopulationToFile( popsf, "swarmPropagation_"+namesnip+"_" + std::to_string( i ) + "_" + std::to_string( i ) , true );

            double bestcost = pops.get_f()[pops.best_idx()][0];

            if (i > n_generations || bestcost == 0){ iterate = false;}
        }


        auto t2 = std::chrono::high_resolution_clock::now();
        std::cout << "Full optimization of " << n_generations << " generations took " <<
                     std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count()<< " seconds using " <<
                     algo_list_names[g] << std::endl;

        std::cout << "Best solution found was for cost: " << swarmProblem.getBestCost() << std::endl;
        // Retrieve final Cartesian states for population in last generation, and save final states to a file.
        std::vector<std::vector< double > > decisionVariables = pops.get_x( );
        std::map< int, Eigen::VectorXd > finalStates;
        for( unsigned int i = 0; i < decisionVariables.size( ); i++ )
        {
            swarmProblem.fitness( decisionVariables.at( i ) );
            finalStates[ i ] = swarmProblem.getPreviousFinalState( );
        }
        tudat::input_output::writeDataMapToTextFile(
                    finalStates, "swarmFinalStates_"+namesnip+".dat", swarm_optimization::getOutputPath( ) + "/FineGrain/" );

        // Write lunar state history to file.
        input_output::writeDataMapToTextFile( swarmProblem.ComputeLunarOrbit(),
                                              "propagationHistory_moon.dat",
                                              swarm_optimization::getOutputPath( ) + "/FineGrain/",
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        // Write core position to file.
        input_output::writeMatrixToFile(swarmProblem.getCorePosition(),
                                        "corePosition_"+namesnip+"_best.dat",10,
                                        swarm_optimization::getOutputPath( ) + "/FineGrain/" ) ;

        // Write perturbed satellite propagation history to file.
        input_output::writeDataMapToTextFile( swarmProblem.getBestStateHistory(),
                                              "propagationHistory_"+namesnip+"_best.dat",
                                              swarm_optimization::getOutputPath( ) + "/FineGrain/",
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        input_output::writeDataMapToTextFile(swarmProblem.getPenalizedBaselineHistoryMap(),
                                              "penaltyHistory_"+namesnip+".dat",
                                              swarm_optimization::getOutputPath( ) + "/FineGrain/",
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              ",");
    }



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

