#include <iostream>

#include <boost/filesystem.hpp>

#include "applicationOutput.h"

#include "Tudat/InputOutput/basicInputOutput.h"

void printPopulationToFile( const std::vector< std::vector< double > >& population,
                            const std::string fileSuffix,
                            const bool isFitness )
{

    Eigen::MatrixXd matrixToPrint( population.size( ), population.at( 0 ).size( ) );
    for( unsigned int i = 0; i < population.size( ); i++ )
    {
        for( unsigned int j = 0; j < population.at( 0 ).size( ); j++ )
        {
            matrixToPrint( i, j ) = population.at( i ).at( j );
        }
    }

    if( !isFitness )
    {
        tudat::input_output::writeMatrixToFile( matrixToPrint, "population_" + fileSuffix + ".dat", 16,
                                                swarm_optimization::getOutputPath( ) );
    }
    else
    {
        tudat::input_output::writeMatrixToFile( matrixToPrint, "fitness_" + fileSuffix + ".dat", 16,
                                                swarm_optimization::getOutputPath( )  );
    }
}

void printPopulationMapToFile( std::map <int, const std::vector< std::vector< double > >&> populationMap,
                            const std::string fileSuffix,
                            const bool isFitness )
{

    auto pop = populationMap.begin()->second;

    int sizemap = populationMap.size();
    int size1 = pop.size();
    int size2 = pop.at(0).size();

    std::cout << "Full map size: " << sizemap << " pop lvl 1 size: " << size1 << " second level size: " << size2 << std::endl;

    Eigen::MatrixXd matrixToPrint( populationMap.size( ), populationMap.at( 0 ).size( )*size2 );
    std::map <int, const std::vector< std::vector< double > >&>::iterator it;

    for( it = populationMap.begin(); it != populationMap.end(); it++){
        const std::vector< std::vector< double > >& population = it->second;
        int i = it->first;
        matrixToPrint(i,0) = it->first;

        for (unsigned int k = 0; k < population.size(); k++){
            for( unsigned int j = 0; j < size2; j++ )
            {

                matrixToPrint( i, j + k*population.at(0).size() ) = population.at( k ).at( j );
            }

        }


        if( !isFitness )
        {
            tudat::input_output::writeMatrixToFile( matrixToPrint, "population_" + fileSuffix + ".dat", 16,
                                                    swarm_optimization::getOutputPath( ) );
        }
        else
        {
            tudat::input_output::writeMatrixToFile( matrixToPrint, "fitness_" + fileSuffix + ".dat", 16,
                                                    swarm_optimization::getOutputPath( )  );
        }
    }


}
