/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_SWARM_OPTIMIZATION_HPP
#define TUDAT_SWARM_OPTIMIZATION_HPP

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include <utility>
#include <vector>

//#include <boost/python.hpp>
//#include <cstdlib>

using namespace tudat;

// Define the problem PaGMO-style
struct SwarmOptimization {

    // Empty constructor
    SwarmOptimization( ){ }

    SwarmOptimization( const int swarmSize,
                                 const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave ,
                       const double missionLength);

    // interpolating data
    std::map< double, Eigen::VectorXd> InterpolateData(std::map< double, Eigen::VectorXd > integrationResult, double stepsize ) const;


    // Fitness: takes the value of the RAAN and returns the value of the closest distance from target
    std::vector<double> fitness(const std::vector<double> &x) const;

    // Boundaries of the problem set between 0 and (360) degrees
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const;

    Eigen::VectorXd getPreviousFinalState( )
    {
        return previousFinalState_;
    }

    std::map< double, Eigen::VectorXd > getPreviousStateHistory( )
    {
        return previousStateHistory_;
    }

    std::map< double, Eigen::VectorXd > getBestStateHistory( )
    {
        return bestStateHistory_;
    }

    std::map< double, Eigen::VectorXd > getPreviousDependentVariablesHistory()
    {
        return previousDependentVariablesHistory_;
    }

    Eigen::VectorXd getPreviousDependentVariablesFinalValues()
    {
        return previousDependentVariablesFinalValues_;
    }

    Eigen::Vector3d getCorePosition(){
        return corePosition_;
    }

    double getBestCost(){
        return bestCost_;
    }

    std::map< double, Eigen::VectorXd > ComputeLunarOrbit()
    {
        using namespace tudat::orbital_element_conversions;
        double mu = bodyMap_["Earth"]->getGravityFieldModel()->getGravitationalParameter();
        std::map< double, Eigen::VectorXd > newmap;
        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = lunarkeplerMap_.begin( );
                 stateIterator != lunarkeplerMap_.end( ); stateIterator++ ){
            Eigen::Vector6d grabber = stateIterator->second;
            newmap.insert(std::pair<double,Eigen::VectorXd>(stateIterator->first , convertKeplerianToCartesianElements(grabber,mu) ));
        }
        return InterpolateData(newmap,interpolationTime_);
    }

    std::vector<double> getPenalizedBaselineHistory(){

        return penalizedBaselineHistory_;
    }
    std::map<int, double> getPenalizedBaselineHistoryMap(){
        std::map<int,double> map;
        for (auto const& value: penalizedBaselineHistory_){
            map.insert(std::pair<int,double>(1,value));
        }
        return map;
    }

private:
    // used
    int swarmSize_;
    double simulationStartEpoch_;
    double simulationEndEpoch_;
    double timestep_;
    double interpolationTime_ = 1*3600;
    mutable double bestCost_;

    mutable Eigen::Vector3d corePosition_;

    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave_;

    mutable std::map< double, Eigen::VectorXd > previousStateHistory_;
    mutable std::map< double, Eigen::VectorXd > bestStateHistory_;
    mutable Eigen::VectorXd previousFinalState_;
    mutable tudat::simulation_setup::NamedBodyMap bodyMap_;
    mutable std::map< double, Eigen::VectorXd > previousDependentVariablesHistory_;
    mutable Eigen::VectorXd previousDependentVariablesFinalValues_;

    mutable std::map< double, Eigen::VectorXd> interpolatedMap_;
    mutable std::map< double, Eigen::VectorXd> lunarkeplerMap_;

    mutable std::vector<double> penalizedBaselineHistory_;

};

#endif // TUDAT_SWARM_OPTIMIZATION_HPP