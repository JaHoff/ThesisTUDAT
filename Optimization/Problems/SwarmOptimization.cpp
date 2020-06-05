/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Mathematics/RootFinders/secantRootFinder.h>

#include "SwarmOptimization.h"

using namespace tudat;


std::map< double, Eigen::VectorXd> SwarmOptimization::InterpolateData(std::map< double, Eigen::VectorXd > integrationResult, double stepsize ) const
{
    /* Interpolate a integrationresult dataset to a fixed time step */

    //std::cout << "Start interpolating the data to a fixed step size" << std::endl;
    using namespace tudat;

//    std::map<double, Eigen::VectorXd> statemap = integrationResult;

    double startTime = simulationStartEpoch_;
    double endTime = simulationEndEpoch_;

    //std::cout << "startTime: " + std::to_string(startTime) + "   endTime: " + std::to_string(endTime) << std::endl;
    int N = std::floor((endTime-startTime)/stepsize);
    //std::cout << "Going to interpolate to " + std::to_string(N) + " data points" << std::endl;
    // create interpolator
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 );

    std::shared_ptr< interpolators::OneDimensionalInterpolator < double, Eigen::VectorXd> > interpolator =
            interpolators::createOneDimensionalInterpolator(
                integrationResult, interpolatorSettings);


    std::vector<double> time(N);
    std::iota(std::begin(time), std::end(time),0);
    // convert the set of points to time series
    std::transform(time.begin(),time.end(), time.begin(),
                   std::bind(std::multiplies<double>(),std::placeholders::_1, stepsize)); //[&day](auto& c){return c*day;}
    // offset these by T0:
    std::transform(time.begin(),time.end(), time.begin(),
                    std::bind2nd(std::plus<double>(), startTime));

    std::map< double, Eigen::VectorXd> interpolatedMap;

    for (double t : time){
        // interpolate
        Eigen::VectorXd interpolated = interpolator->interpolate(t);
        std::pair<double, Eigen::VectorXd> insert(t,interpolated);
        interpolatedMap.insert(insert);
    }

    //std::cout << "Finished interpolating" << std::endl;
    return interpolatedMap;
}


/*Initialisation of the problem should include most of the general setup, that would not vary with iterations*/
SwarmOptimization::SwarmOptimization(const int swarmSize,
        const std::shared_ptr< propagators::DependentVariableSaveSettings> dependentVariablesToSave,
                                     const double missionLength) :
    swarmSize_( swarmSize ),
    dependentVariablesToSave_( dependentVariablesToSave )
{
    //std::cout << "Started constructing the base problem" << std::endl;
//    using namespace tudat;
//    using namespace tudat::simulation_setup;
//    using namespace tudat::propagators;
//    using namespace tudat::numerical_integrators;
//    using namespace tudat::basic_astrodynamics;
//    using namespace tudat::basic_mathematics;
//    using namespace tudat::orbital_element_conversions;
//    using namespace tudat::unit_conversions;
//    using namespace tudat::input_output;

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;

    bestCost_ = 1e8;

    // Load Spice kernels.
    //const string kernelfile = input_output::getSpiceKernelPath() + "tudat_merged_edit.bsp";
    std::vector< std::string > kernelfile;
    kernelfile.push_back(input_output::getSpiceKernelPath() + "tudat_merged_edit.bsp");
    spice_interface::loadStandardSpiceKernels(kernelfile);

    // Set simulation time settings.
    simulationStartEpoch_ = 30*tudat::physical_constants::JULIAN_YEAR; // 2030 start date
    simulationEndEpoch_ = simulationStartEpoch_ +  missionLength; //
    const double timestep_ = 30*60;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mercury" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Venus" );
    bodiesToCreate.push_back("Jupiter");
    bodiesToCreate.push_back("Saturn");

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch_ - 3*timestep_, simulationEndEpoch_ + 3*timestep_ );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    bodyMap_ = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLES           /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create radiation pressure settings
    double referenceAreaRadiation = 0.3416;
    double radiationPressureCoefficient = 1.1421;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > RadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    for (int i=0; i < swarmSize_; i++){
        // Create spacecraft object.
        std::string name = std::to_string(i);
        bodyMap_[ name ] = std::make_shared< simulation_setup::Body >( );
        bodyMap_[ name ]->setConstantBodyMass( 5.0 );

        // Create and set radiation pressure settings
        bodyMap_[ name ]->setRadiationPressureInterface(
                    "Sun", createRadiationPressureInterface(
                        RadiationPressureSettings, name, bodyMap_ ));
    }

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap_, "SSB", "J2000" );

    //std::cout << "Finalized constructing the base problem" << std::endl;
}

/*Fitness function, convert the variables to a swarm constellation in orbit and propagate,
Return the evaluated cost. This will need the inclusion of a Python hook*/
std::vector<double> SwarmOptimization::fitness(const std::vector<double> &x) const
{
//    using namespace tudat;
//    using namespace tudat::simulation_setup;
//    using namespace tudat::propagators;
//    using namespace tudat::numerical_integrators;
//    using namespace tudat::basic_astrodynamics;
//    using namespace tudat::basic_mathematics;
//    using namespace tudat::orbital_element_conversions;
//    using namespace tudat::unit_conversions;
//    using namespace tudat::input_output;

    //std::cout << "Start evaluating the cost function" << std::endl;

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Starting fitness function for a problem with " + std::to_string(x.size()) + " variables" << std::endl;
    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSats;
    accelerationsOfSats[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 6, 6 ) );
    accelerationsOfSats[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfSats[ "Moon" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 6, 6 )  );
    accelerationsOfSats[ "Mercury" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfSats[ "Venus" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfSats[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfSats["Saturn"].push_back( std::make_shared< AccelerationSettings> (
                                                 basic_astrodynamics::central_gravity));
    accelerationsOfSats[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::cannon_ball_radiation_pressure ) );

    for ( int i =0; i < swarmSize_; i++){
        std::string name = std::to_string(i);
        bodiesToPropagate.push_back( name );
        centralBodies.push_back( "Earth" );

        accelerationMap[ name ] = accelerationsOfSats;
    }

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, bodiesToPropagate, centralBodies );
    //std::cout << "Acceleration map succesfully defined" << std::endl;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Could be moved to setup phase to speed up code
    Eigen::Vector6d moonInitialState = getInitialStateOfBody("Moon", "Earth", bodyMap_, simulationStartEpoch_);
    Eigen::Vector3d earthMoonVector = moonInitialState.segment(0,3);
    Eigen::Vector3d moonVelocity = moonInitialState.segment(3,3);

    Eigen::Vector3d moonMomentum = earthMoonVector.cross(moonVelocity);
    Eigen::Vector3d L4dir = -earthMoonVector.cross(moonMomentum);
    Eigen::Vector3d L4Cart = 0.5*earthMoonVector + 0.5*sqrt(3)*earthMoonVector.norm()*L4dir.normalized();

//    std::cout << "L4 position: " << L4Cart << std::endl;
    // Displace the core of the swarm based off a variable
    Eigen::Vector3d coreDisplacement = Eigen::Vector3d();
    coreDisplacement(0) = x[0];
    coreDisplacement(1) = x[1];
    coreDisplacement(2) = x[2];
    Eigen::Vector3d corePosition = L4Cart + coreDisplacement;
    // Set the swarm velocity based off the core position plus a additionalVelocity variable
    Eigen::Vector3d stableVelocity = -L4Cart.cross(moonMomentum).normalized() * moonVelocity.norm();
    Eigen::Vector3d additionalVelocity = {x[3],x[4],x[5]};

    // Distribute the initial states of swarm elements around the core:

    // update position of the core
    corePosition += coreDisplacement;

    // Set initial state
    Eigen::VectorXd systemInitialState = Eigen::VectorXd( 6*swarmSize_);

    //std::cout << "Start displacing satellites" << std::endl;
    for (int i = 0; i < swarmSize_; i++){
        Eigen::Vector3d Displacement = Eigen::Vector3d();
        Displacement(0) = x[6+3*i]; // 1: 100
        Displacement(1) = x[7+3*i];
        Displacement(2) = x[8+3*i];

        Eigen::Vector6d InitialState = Eigen::Vector6d();
        InitialState.segment(0,3) = corePosition + Displacement;
        InitialState.segment(3,3) = stableVelocity + additionalVelocity;

        systemInitialState.segment( 6*i, 6 ) = InitialState;
    }

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch_,propagators::cowell, dependentVariablesToSave_ );

    const double minimumStepSize = 0.01;
    const double maximumStepSize = 24*3600;
    const double relativeErrorTolerance = 1E-11;
    const double absoluteErrorTolerance = 1E-11;
    std::shared_ptr< RungeKuttaVariableStepSizeSettings<> > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings <> >
            ( simulationStartEpoch_, 10, RungeKuttaCoefficients::rungeKuttaFehlberg45,minimumStepSize,maximumStepSize,
              relativeErrorTolerance, absoluteErrorTolerance);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Starting propagation from epoch " + std::to_string(simulationStartEpoch_) + " to " + std::to_string(simulationEndEpoch_) << std::endl;
    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettings, propagatorSettings );

    //std::cout << "Fully simulated the satellite motion" << std::endl;
    //Retrieve results
    std::map< double, Eigen::VectorXd > integrationResult = InterpolateData(dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),interpolationTime_);
    previousStateHistory_ = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    previousFinalState_ = previousStateHistory_.rbegin( )->second;

//    previousDependentVariablesFinalValues_ = previousDependentVariablesHistory_.rbegin()->second;

    /* COST FUNCTION */

    // Simple cost function; penalize too large or too little baselines
    //std::cout << "Starting cost function" << std::endl;
    double cost = 0.;
    int count = 0;
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
             stateIterator != integrationResult.end( ); stateIterator++ ){
        for ( int i = 0; i < swarmSize_; i++){
            for ( int j = 0; j <=  swarmSize_ - i; j++){
                count++;
                double baseline = (stateIterator->second.segment(6*i,3) - stateIterator->second.segment(6*(swarmSize_ -1- j),3)).norm();
                if (baseline > 0 && (baseline > 100e3 || baseline< 500)){ cost++;
                //std::cout << "Penalized a baseline with size: " + std::to_string(baseline.norm()) << std::endl;
                }
            }
        }

    }

    //std::cout << "cost function evaluated a total of " << count << " baselines" << std::endl;
    if (cost < bestCost_){
        std::cout << "Updated bestCost_ to:" << cost << std::endl;
        bestCost_= cost;
        bestStateHistory_ = integrationResult;

        lunarkeplerMap_ = dynamicsSimulator.getDependentVariableHistory();
        corePosition_ = corePosition;
    }

    std::vector<double> output = {cost};

    //std::cout << "Succesfully computed cost function with cost: " + std::to_string(cost) << std::endl;
    return output;


}


std::pair<std::vector<double>, std::vector<double>> SwarmOptimization::get_bounds() const
{
    //std::cout << "I just tried to access the problem bounds!" << std::endl;
    std::vector<double> lowerbounds = {-50.e3, -50.e3, -50.e3,-20, -20., -20.};
    std::vector<double> upperbounds = {50.e3, 50.e3, 50.e3 , 20., 20., 20.};
    std::vector<double> lowerelementbounds = {-50.e3, -50.e3, -50.e3};
    std::vector<double> upperelementbounds = {50.e3, 50.e3, 50.e3};

    for (int i =0; i < swarmSize_; i++){
        lowerbounds.insert(lowerbounds.end(),lowerelementbounds.begin(), lowerelementbounds.end());
        upperbounds.insert(upperbounds.end(), upperelementbounds.begin(), upperelementbounds.end());
    }

    //std::cout << "No problems in defining the bounds" << std::endl;
    return {lowerbounds,
        upperbounds};
}


