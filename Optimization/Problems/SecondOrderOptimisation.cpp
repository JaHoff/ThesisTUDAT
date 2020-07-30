/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *
 * Implementation by Jurriaan van 't Hoff
 * Code contains a lot of console output commands that track progress, which are commented out by default.
 * Recommend decommenting these using find and replace for debugging crashes
 */



#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Mathematics/RootFinders/secantRootFinder.h>

#include "SecondOrderOptimisation.h"

#include <thread>
#include <future>

using namespace tudat;

std::map< double, Eigen::VectorXd> SecondOrderOptimisation::InterpolateData(std::map< double, Eigen::VectorXd > integrationResult, double stepsize ) const
{
    /* Interpolate a integrationresult dataset to a fixed time step */

    //std::cout << "Start interpolating the data to a fixed step size" << std::endl;
    using namespace tudat;

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
SecondOrderOptimisation::SecondOrderOptimisation(const int swarmSize,
        const std::shared_ptr< propagators::DependentVariableSaveSettings> dependentVariablesToSave,
                                     const double missionLength) :
    swarmSize_( swarmSize ),
    dependentVariablesToSave_( dependentVariablesToSave )
{
    //std::cout << "Started constructing the base problem" << std::endl;
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;

    // Load Spice kernels.
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
        bodySettings[ bodiesToCreate.at( i )]->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            simulationStartEpoch_, simulationEndEpoch_, 30*60, "SSB", "J2000" );
        bodySettings[ bodiesToCreate.at( i )] -> rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                    "J2000", "J2000", spice_interface::computeRotationQuaternionBetweenFrames(
                        "J2000", "J2000", simulationStartEpoch_ ),
                    simulationStartEpoch_, 2.0 * mathematical_constants::PI / physical_constants::SIDEREAL_DAY );
    }

    bodySettings[ "Earth" ] -> rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "J2000", "IAU_Earth", spice_interface::computeRotationQuaternionBetweenFrames(
                    "J2000", "IAU_Earth", simulationStartEpoch_ ),
                simulationStartEpoch_, 2.0 * mathematical_constants::PI / physical_constants::SIDEREAL_DAY );
    bodySettings[ "Moon" ] -> rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "J2000", "IAU_Moon", spice_interface::computeRotationQuaternionBetweenFrames(
                    "J2000", "IAU_Moon", simulationStartEpoch_ ),
                simulationStartEpoch_, 2.0 * mathematical_constants::PI / physical_constants::SIDEREAL_DAY );

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

    // Could be moved to setup phase to speed up code
    Eigen::Vector6d moonInitialState = getInitialStateOfBody("Moon", "Earth", bodyMap_, simulationStartEpoch_);
    Eigen::Vector3d earthMoonVector = moonInitialState.segment(0,3);
    moonVelocity_ = moonInitialState.segment(3,3);

    moonMomentum_ = earthMoonVector.cross(moonVelocity_);
    Eigen::Vector3d L4dir = -earthMoonVector.cross(moonMomentum_);
    L4Cart_ = 0.5*earthMoonVector + 0.5*sqrt(3)*earthMoonVector.norm()*L4dir.normalized();
    //std::cout << "Finalized constructing the base problem" << std::endl;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;

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
        bodiesToPropagate_.push_back( name );
        centralBodies_.push_back( "Earth" );

        accelerationMap[ name ] = accelerationsOfSats;
    }

    accelerationModelMap_ = createAccelerationModelsMap(
                bodyMap_, accelerationMap, bodiesToPropagate_, centralBodies_ );

    //std::cout << "Acceleration model map succesfully defined" << std::endl;

    const double minimumStepSize = 0.01;
    const double maximumStepSize = 24*3600;
    const double relativeErrorTolerance = 1E-11;
    const double absoluteErrorTolerance = 1E-11;
    integratorSettings_ = std::make_shared< RungeKuttaVariableStepSizeSettings <> >
                            ( simulationStartEpoch_, 10, RungeKuttaCoefficients::rungeKuttaFehlberg45,minimumStepSize,maximumStepSize,
                              relativeErrorTolerance, absoluteErrorTolerance);
}

/*Fitness function, convert the variables to a swarm constellation in orbit and propagate,
Return the evaluated cost. This will need the inclusion of a Python hook*/
std::vector<double> SecondOrderOptimisation::fitness(const std::vector<double> &x) const
{

    //std::cout << "Starting fitness function for a problem with " + std::to_string(x.size()) + " variables" << std::endl;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///

//    //std::cout << "L4 position: " << L4Cart << std::endl;
    // Displace the core of the swarm based off a variable
    Eigen::Vector3d coreDisplacement = Eigen::Vector3d();
    coreDisplacement(0) = x[0];
    coreDisplacement(1) = x[1];
    coreDisplacement(2) = x[2];
    Eigen::Vector3d corePosition = L4Cart_ + coreDisplacement;
    // Set the swarm velocity based off the core position plus a additionalVelocity variable
    Eigen::Vector3d stableVelocity = -L4Cart_.cross(moonMomentum_).normalized() * moonVelocity_.norm();
    Eigen::Vector3d additionalVelocity = {x[3],x[4],x[5]};

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
            ( centralBodies_, accelerationModelMap_, bodiesToPropagate_, systemInitialState,
              simulationEndEpoch_,propagators::cowell, dependentVariablesToSave_ );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Starting propagation from epoch " + std::to_string(simulationStartEpoch_) + " to " + std::to_string(simulationEndEpoch_) << std::endl;
    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettings_, propagatorSettings );

    //std::cout << "Fully simulated the satellite motion" << std::endl;

    //Retrieve results
    std::map< double, Eigen::VectorXd > integrationResult = InterpolateData(dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),interpolationTime_);
    previousStateHistory_ = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    previousFinalState_ = previousStateHistory_.rbegin( )->second;


    /* COST FUNCTION */

    // Simple cost function; penalize too large or too little baselines


    //std::cout << "Starting cost function" << std::endl;
    std::vector<double> penalizedBaselineHistory;
    double cost = 0.;
    int count = 0;
    double baseline, baselinerate;

    int bsit1, bsit2, bsrit1, bsrit2;
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
             stateIterator != integrationResult.end( ); stateIterator++ ){
        for ( int i = 0; i < swarmSize_; i++){
            for ( int j = 0; j <  swarmSize_ -1 -i; j++){
                count++;

                bsit1 = 6*i;
                bsit2 = 6*(swarmSize_- j -1);
                bsrit1 = bsit1 + 3;
                bsrit2 = bsit2 + 3;

                ////std::cout << "baseline comp will try to access values " << bsit1 << " to " << bsit2 << std::endl;
                ////std::cout << "baselinerate comp will try to access values " << bsrit1 << " to " << bsrit2 << std::endl;
                baseline = (stateIterator->second.segment(bsit1,3) - stateIterator->second.segment(bsit2,3)).norm();
                baselinerate = (stateIterator->second.segment(bsrit1,3) - stateIterator->second.segment(bsrit2,3)).norm();
                if (baselinerate > 1 || (baseline > 100e3 || baseline< 500)){ cost++;
                    penalizedBaselineHistory.push_back(baseline);
                    penalizedBaselineHistory.push_back(baselinerate);
                //std::cout << "Penalized a baseline with size: " + std::to_string(baseline) << std::endl;
                }
            }
        }

    }

    //std::cout << "cost function evaluated a total of " << count << " baselines" << std::endl;
    if (cost < bestCost_){

        bestCost_= cost;
        std::cout << "Updated bestCost_ to:" << bestCost_ << std::endl;
        bestStateHistory_ = integrationResult;
        penalizedBaselineHistory_ = penalizedBaselineHistory;
        lunarkeplerMap_ = dynamicsSimulator.getDependentVariableHistory();
        corePosition_ = corePosition;
        coreVelocity_ = additionalVelocity;
        bestPopulationData_ = x;
    }

    std::vector<double> output = {cost};

    //std::cout << "Succesfully computed cost function with cost: " + std::to_string(cost) << std::endl;
    return output;


}


std::pair<std::vector<double>, std::vector<double>> SecondOrderOptimisation::get_bounds() const
{
    std::vector<double> lowerbounds;
    std::vector<double> upperbounds;
    std::vector<double> lowerpositionbounds = {-50.e3, -50.e3, -50.e3};
    std::vector<double> upperpositionbounds = {50.e3, 50.e3, 50.e3};
    std::vector<double> lowervelocitybounds = {-50, -50, -50};
    std::vector<double> uppervelocitybounds = {50, 50, 50};

    // Insert a position, velocity pair
    for (int i =0; i < swarmSize_; i++){
        lowerbounds.insert(lowerbounds.end(),lowerpositionbounds.begin(),lowerpositionbounds.end());
        upperbounds.insert(upperbounds.end(),upperpositionbounds.begin(),upperpositionbounds.end());

        lowerbounds.insert(lowerbounds.end(),lowervelocitybounds.begin(), lowervelocitybounds.end());
        upperbounds.insert(upperbounds.end(), uppervelocitybounds.begin(), uppervelocitybounds.end());
    }

    //std::cout << "No problems in defining the bounds" << std::endl;
    return {lowerbounds,
        upperbounds};
}






