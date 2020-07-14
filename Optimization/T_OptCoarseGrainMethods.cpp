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
    int n_generations = 50;
    int n_islands = 5;
    int n_pops = 32;
    int r_seed = 72;
    int n_sats = 20;
    int n_days = 365;

    double missionLength = n_days*tudat::physical_constants::JULIAN_DAY;


    // The number of internal iterations a island goes through before the next global generation, yields more efficient progress per iteration, but slower generation computations
    int internalIterations = 5;

    string subfolder = "/Coarse/";
    std::cout << "General optimization start!" << std::endl;

    std::vector<std::string> algo_list_names{"Differential Evolution", "Self Adjusting Differential Evolution",
                                            "Particle Swarm Optimization",
                                            "Particle Swarm Optimization Generational",
                                             "Generational Ant Colony"};
    std::vector<std::string> algo_names_shorthand{"de1220", "sade", "pso", "pso_gen", "gaco"};

    int algochoice = 0;
    string namesnip = algo_names_shorthand.at(algochoice) + "_sd" + std::to_string(r_seed) +
            "_sats" + std::to_string(n_sats) + "_nisl" + std::to_string(n_islands) + "_npop" + std::to_string(n_pops) +
            "_int" + std::to_string(internalIterations) + "_d" + std::to_string(n_days);

    //Set seed for reproducible results
    pagmo::random_device::set_seed(r_seed);


    // Instantiate a pagmo algorithm
     // was 20 gens internally the magic nr?
    algorithm algo;

    switch (algochoice){
    case 0:
        algo = de1220(internalIterations);
        break;
    default:
        algo = de1220();
        break;
    }



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

    auto t1 = std::chrono::high_resolution_clock::now();



    // Convert the size settings to the appropiate types
    pagmo::population::size_type populationSize = n_pops;


//     swarmProblem;
    std::vector< SwarmOptimization> swarmProblems;
    archipelago arch;
    for (int i = 0; i < n_islands; i++){

        SwarmOptimization swarmProblem = SwarmOptimization( n_sats, dependentVariablesToSave, missionLength );\
        swarmProblems.push_back(swarmProblem);
        std::cout << "Problemize the problem" << std::endl;
        problem prob{ swarmProblem };

        arch.push_back(algo,prob,populationSize);
        std::cout << "added island nr " << i << std::endl;
    }


    std::map <int, std::vector< double >> fitnessmap;
    std::cout << "Starting evolving optimization problem for "<< n_generations << " generations!" << std::endl;
    int i = 0;
    bool iterate = true;
    while( iterate)
    {
        std::cout << "Now at generation " + std::to_string(i) << std::endl;
        i++;
        arch.evolve(2 );
        while( arch.status( ) != pagmo::evolve_status::idle &&
               arch.status( ) != pagmo::evolve_status::idle_error )
        {
            arch.wait( );
        }
        arch.wait_check( ); // Raises errors
        std::cout << "evolution done!" << std::endl;

        fitnessmap.insert( std::pair<double, std::vector< double > >( i, arch.begin()->get_population().get_f().at(0) ) );
        // Write current iteration results to file


        double bestcost = 1e8;
        int islandcount = 0;

        // iterate through the islands and store their generational data in files
        for (auto it = arch.begin(); it != arch.end(); it++){
            std::vector<vector_double> popsf = it->get_population().get_f();
            auto popsx = it->get_population().get_x();
            printPopulationToFile( popsx, namesnip+"_g" + std::to_string( i ) + "_i" + std::to_string(islandcount) , false, subfolder);
            printPopulationToFile( popsf, namesnip+"_g" + std::to_string( i ) + "_i" + std::to_string(islandcount) , true ,subfolder);

            population pops = it->get_population();

            // Store the best fitness of this generation
            double bestgencost = pops.get_f()[pops.best_idx()][0];
            if (bestgencost < bestcost) {bestcost = bestgencost;}
            islandcount++;
        }

        // If an optima is found, stop iterating
        // Else continue until the given maximum
        if (i > n_generations || bestcost == 0){ iterate = false;}
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Full optimization of " << n_generations << " generations took " <<
                 std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count()<< " seconds" <<  std::endl;

    int c = 0;
    // Iterate through the islands to write the relevant data to files
    for (auto it = arch.begin(); it != arch.end(); it++){

        auto SP = swarmProblems.at( c);
        // Retrieve final Cartesian states for population in last generation, and save final states to a file.
        std::vector<std::vector< double > > decisionVariables = it->get_population( ).get_x( );
        std::map< int, Eigen::VectorXd > finalStates;
        for( unsigned int i = 0; i < decisionVariables.size( ); i++ )
        {
            SP.fitness( decisionVariables.at( i ) );
            finalStates[ i ] = SP.getPreviousFinalState( );
        }
        tudat::input_output::writeDataMapToTextFile(
                    finalStates, "swarmFinalStates_"+namesnip+".dat", swarm_optimization::getOutputPath( ) + subfolder );

        // Write lunar state history to file.
        input_output::writeDataMapToTextFile( SP.ComputeLunarOrbit(),
                                              "propagationHistory_moon.dat",
                                              swarm_optimization::getOutputPath( ) + subfolder,
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        // Write core position to file.
        input_output::writeMatrixToFile(SP.getCorePosition(),
                                        "corePosition_"+namesnip+"_best.dat",10,
                                        swarm_optimization::getOutputPath( ) + subfolder ) ;

        // Write the champion orbit data to a file.
        input_output::writeDataMapToTextFile( SP.getBestStateHistory(),
                                              "propagationHistory_"+namesnip+"_best.dat",
                                              swarm_optimization::getOutputPath( ) + subfolder,
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        // Write down files that show the baselines which triggered the cost function
        input_output::writeDataMapToTextFile(SP.getPenalizedBaselineHistoryMap(),
                                              "penaltyHistory_"+namesnip+".dat",
                                              swarm_optimization::getOutputPath( ) + subfolder,
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              ",");

    }




}

