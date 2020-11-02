/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


/* Propagate the orbit of a optimisation result / swarm constellation for a set timeframe */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Problems/applicationOutput.h"
//#include "Problems/saveOptimizationResults.h"


std::map< double, Eigen::VectorXd> InterpolateData(std::map< double, Eigen::VectorXd > integrationResult, double stepsize , double startTime, double endTime)
{
    /* Interpolate a integrationresult dataset to a fixed time step
       Adjusted from the swarmoptimisation interpolation variant*/

    //std::cout << "Start interpolating the data to a fixed step size" << std::endl;
    using namespace tudat;

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

//! Place two satellites in a halo orbit around L4, using a variable time-step propagator
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            SETTINGS                      //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2nd order optimisation result?
    bool SecondOrder = true;
    // Filename to look for
    std::string filename = "population_35_2ndorder_25.dat";
    // Relative folder of file (RELATIVE TO THE .EXE - /bin/applications/xx)
    std::string relFolder = "Propagate/";
    //Relative subfolder to place the output
    std::string outputSubFolder = "champions_propagated/35sats2ndorder25/";

    // Name attachment to differentiate result files
    std::string attachment = "35sat_champ";

    // Integer name of the population member to propagate for
    int popmember = 0; // Remember - count starts at 0

    // for how long to propagate the orbit
    double daysToPropagate = 365*5;

    // Interpolation timesteps
    double interpolationTime = 4*3600;

    /// Dependent variables
    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared < propagators::SingleDependentVariableSaveSettings> (
                                          propagators::keplerian_state_dependent_variable, "Moon", "Earth"));

    // Create object with list of dependent variables.
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            LOADING FILE                  //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //auto x = tudat::input_output::readEigenVectorMapFromFile( relFolder + filename);

    auto x = tudat::input_output::readMatrixFromFile(relFolder + filename);

    std::cout << "found data!" << std::endl;

    Eigen::VectorXd popdata = x.row(popmember);
    int n_satellites;
    if (SecondOrder){
       n_satellites  = (popdata.size()) / 6;
    }
    else{
       n_satellites = (popdata.size() - 6) / 3;
    }


    std::cout << "Grabbed population number " << popmember << ", retrieved following vector:" << std::endl;
    std::cout << popdata << std::endl;
    std::cout << "Identified " << n_satellites << " satellites from this data file" << std::endl;



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            SETUP                         //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Load Spice kernels.
    std::vector< std::string > kernelfile;
    kernelfile.push_back(input_output::getSpiceKernelPath() + "tudat_merged_edit.bsp");
    spice_interface::loadStandardSpiceKernels(kernelfile);

    // Set simulation time settings.
    double simulationStartEpoch = 30*tudat::physical_constants::JULIAN_YEAR; // 2030 start date
    double simulationEndEpoch = simulationStartEpoch +  daysToPropagate*tudat::physical_constants::JULIAN_DAY; //
    const double timestep_ = 30*60;

    std::cout << "Will propagate from epoch" << simulationStartEpoch << " to " <<simulationEndEpoch << std::endl;
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
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 3*timestep_, simulationEndEpoch + 3*timestep_ );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
        bodySettings[ bodiesToCreate.at( i )]->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            simulationStartEpoch, simulationEndEpoch, 30*60, "SSB", "J2000" );
        bodySettings[ bodiesToCreate.at( i )] -> rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                    "J2000", "J2000", spice_interface::computeRotationQuaternionBetweenFrames(
                        "J2000", "J2000", simulationStartEpoch ),
                    simulationStartEpoch, 2.0 * mathematical_constants::PI / physical_constants::SIDEREAL_DAY );
    }

    bodySettings[ "Earth" ] -> rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "J2000", "IAU_Earth", spice_interface::computeRotationQuaternionBetweenFrames(
                    "J2000", "IAU_Earth", simulationStartEpoch ),
                simulationStartEpoch, 2.0 * mathematical_constants::PI / physical_constants::SIDEREAL_DAY );
    bodySettings[ "Moon" ] -> rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "J2000", "IAU_Moon", spice_interface::computeRotationQuaternionBetweenFrames(
                    "J2000", "IAU_Moon", simulationStartEpoch ),
                simulationStartEpoch, 2.0 * mathematical_constants::PI / physical_constants::SIDEREAL_DAY );

    auto bodyMap = createBodies( bodySettings );

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

    for (int i=0; i < n_satellites; i++){
        // Create spacecraft object.
        std::string name = std::to_string(i);
        bodyMap[ name ] = std::make_shared< simulation_setup::Body >( );
        bodyMap[ name ]->setConstantBodyMass( 5.0 );

        // Create and set radiation pressure settings
        bodyMap[ name ]->setRadiationPressureInterface(
                    "Sun", createRadiationPressureInterface(
                        RadiationPressureSettings, name, bodyMap ));
    }

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    // Could be moved to setup phase to speed up code
    Eigen::Vector6d moonInitialState = getInitialStateOfBody("Moon", "Earth", bodyMap, simulationStartEpoch);
    Eigen::Vector3d earthMoonVector = moonInitialState.segment(0,3);
    Eigen::Vector3d moonVelocity = moonInitialState.segment(3,3);

    Eigen::Vector3d moonMomentum = earthMoonVector.cross(moonVelocity);
    Eigen::Vector3d L4dir = -earthMoonVector.cross(moonMomentum);
    Eigen::Vector3d L4Cart = 0.5*earthMoonVector + 0.5*sqrt(3)*earthMoonVector.norm()*L4dir.normalized();
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

    // Create and fill bodymaps
    std::vector< std::string > bodiesToPropagate,centralBodies;

    for ( int i =0; i < n_satellites; i++){
        std::string name = std::to_string(i);
        bodiesToPropagate.push_back( name );
        centralBodies.push_back( "Earth" );

        accelerationMap[ name ] = accelerationsOfSats;
    }


    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    //std::cout << "Acceleration model map succesfully defined" << std::endl;

    const double minimumStepSize = 1e-6;
    const double maximumStepSize = 24*3600;
    const double relativeErrorTolerance = 1E-11;
    const double absoluteErrorTolerance = 1E-11;
    std::shared_ptr< RungeKuttaVariableStepSizeSettings<> > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings <> >
                            ( simulationStartEpoch, 10, RungeKuttaCoefficients::rungeKuttaFehlberg45,minimumStepSize,maximumStepSize,
                              relativeErrorTolerance, absoluteErrorTolerance);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///

//std::cout << "L4 position: " << L4Cart << std::endl;
    // Displace the core of the swarm based off a variable
    Eigen::Vector3d coreDisplacement = Eigen::Vector3d();
    Eigen::Vector3d additionalVelocity = Eigen::Vector3d();
    if (SecondOrder){
        coreDisplacement = {4451595.805418905,
                            45662.36398591232,
                            -3765.688752660145};
        additionalVelocity = {2.181963044975345,
                              -13.29823474082929,
                              -12.3799288546268};
    }
    else{
        coreDisplacement(0) = popdata[0];
        coreDisplacement(1) = popdata[1];
        coreDisplacement(2) = popdata[2];
        additionalVelocity = {popdata[3],popdata[4],popdata[5]};
    }


    Eigen::Vector3d corePosition = L4Cart + coreDisplacement;
    // Set the swarm velocity based off the core position plus a additionalVelocity variable
    Eigen::Vector3d stableVelocity = -L4Cart.cross(moonMomentum).normalized() * moonVelocity.norm();


    // Set initial state
    Eigen::VectorXd systemInitialState = Eigen::VectorXd( 6*n_satellites);

    //std::cout << "Start displacing satellites" << std::endl;
    for (int i = 0; i < n_satellites; i++){
        Eigen::Vector3d Displacement = Eigen::Vector3d();

        if (SecondOrder){
            Displacement(0) = popdata[0 + 6*i]; // 1: 100
            Displacement(1) = popdata[1 + 6*i];
            Displacement(2) = popdata[2 + 6*i];
        }
        else{
            Displacement(0) = popdata[6+3*i]; // 1: 100
            Displacement(1) = popdata[7+3*i];
            Displacement(2) = popdata[8+3*i];
        }

        Eigen::Vector6d InitialState = Eigen::Vector6d();
        InitialState.segment(0,3) = corePosition + Displacement;
        InitialState.segment(3,3) = stableVelocity + additionalVelocity;

        systemInitialState.segment( 6*i, 6 ) = InitialState;
    }

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              simulationEndEpoch,propagators::cowell, dependentVariablesToSave );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Starting propagation from epoch " + std::to_string(simulationStartEpoch) + " to " + std::to_string(simulationEndEpoch_) << std::endl;
    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );

    //std::cout << "Fully simulated the satellite motion" << std::endl;

    //Retrieve results
    std::map< double, Eigen::VectorXd > integrationResult =
            InterpolateData(dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),interpolationTime, simulationStartEpoch, simulationEndEpoch);

    //Interpolate moon position

    std::map<double, Eigen::VectorXd> lunarkeplerMap = dynamicsSimulator.getDependentVariableHistory();
    double mu = bodyMap["Earth"]->getGravityFieldModel()->getGravitationalParameter();
    std::map< double, Eigen::VectorXd > newmap;
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = lunarkeplerMap.begin( );
             stateIterator != lunarkeplerMap.end( ); stateIterator++ ){
        Eigen::Vector6d grabber = stateIterator->second;
        newmap.insert(std::pair<double,Eigen::VectorXd>(stateIterator->first , convertKeplerianToCartesianElements(grabber,mu) ));
    }
    std::map< double, Eigen::VectorXd > moonState = InterpolateData(newmap,interpolationTime,simulationStartEpoch,simulationEndEpoch);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( integrationResult,
                                          "propagationHistory_" + attachment +".dat",
                                          swarm_optimization::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write lunar state history to file.
    input_output::writeDataMapToTextFile( moonState,
                                          "propagationHistory_"+attachment+"_moon.dat",
                                          swarm_optimization::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write core position to file.
    input_output::writeMatrixToFile(corePosition,
                                    "corePosition_"+attachment+".dat",10,
                                    swarm_optimization::getOutputPath( ) + outputSubFolder ) ;

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
