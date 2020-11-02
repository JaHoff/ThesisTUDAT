/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

/* Generate a set of orbits for relative placement along a grid, used to map the cost function along a 2-d plane */
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "TestApp/applicationOutput.h"

void InterpolateAndSave(std::map< double, Eigen::VectorXd > integrationResult, std::string name, double startTime,
                        int N = 365, double stepsize = tudat::physical_constants::JULIAN_DAY){
    using namespace tudat;
    std::string outputSubFolder = "costMapping/";

    std::map<double, Eigen::VectorXd> statemap = integrationResult;



    // create interpolator
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 );

    std::shared_ptr< interpolators::OneDimensionalInterpolator < double, Eigen::VectorXd> > interpolator =
            interpolators::createOneDimensionalInterpolator(
                statemap, interpolatorSettings);


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


    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( interpolatedMap,
                                          name + ".dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );
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
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Save name settings

    string attachment;

    // Load Spice kernels.
    //const string kernelfile = input_output::getSpiceKernelPath() + "tudat_merged_edit.bsp";
    std::vector< std::string > kernelfile;
    kernelfile.push_back(input_output::getSpiceKernelPath() + "tudat_merged_edit.bsp");
    spice_interface::loadStandardSpiceKernels(kernelfile);

    // Set simulation time settings.
    const double simulationStartEpoch = 30*tudat::physical_constants::JULIAN_YEAR; // 2030 start date
    const double simulationEndEpoch = simulationStartEpoch +  1*tudat::physical_constants::JULIAN_YEAR; //
    const double timestep = 30*60;

    const double interpolationStep  = 4*3600;
    const int interpolationN = int(1*tudat::physical_constants::JULIAN_YEAR/interpolationStep);

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
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 3*timestep, simulationEndEpoch + 3*timestep );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLES           /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Asterix" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Asterix" ]->setConstantBodyMass( 5.0 );
    bodyMap[ "Obelix" ] = std::make_shared< simulation_setup::Body >();
    bodyMap[ "Obelix" ]->setConstantBodyMass(5.0);
    bodyMap[ "Tetris" ] = std::make_shared< simulation_setup::Body >();
    bodyMap[ "Tetris" ]->setConstantBodyMass(5.0);

    // Create radiation pressure settings
    double referenceAreaRadiation = 0.3416;
    double radiationPressureCoefficient = 1.1421;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Asterix" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Asterix", bodyMap ));
    bodyMap[ "Obelix" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Obelix", bodyMap));
    bodyMap[ "Tetris" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Obelix", bodyMap));
    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );
    bodiesToPropagate.push_back("Obelix");
    centralBodies.push_back("Earth");
    bodiesToPropagate.push_back("Tetris");
    centralBodies.push_back("Earth");

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

    accelerationMap[ "Asterix" ] = accelerationsOfSats;
    accelerationMap[ "Obelix" ] = accelerationsOfSats;
    accelerationMap[ "Tetris" ] = accelerationsOfSats;

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::Vector6d moonInitialState = getInitialStateOfBody("Moon", "Earth", bodyMap, simulationStartEpoch);
    Eigen::Vector3d earthMoonVector = moonInitialState.segment(0,3);
    Eigen::Vector3d moonVelocity = moonInitialState.segment(3,3);

    Eigen::Vector3d moonMomentum = earthMoonVector.cross(moonVelocity);
    Eigen::Vector3d L4dir = -earthMoonVector.cross(moonMomentum);
    Eigen::Vector3d L4Cart = 0.5*earthMoonVector + 0.5*sqrt(3)*earthMoonVector.norm()*L4dir.normalized();


    Eigen::Vector3d L4initVelocity = -L4Cart.cross(moonMomentum).normalized() * moonVelocity.norm();
    Eigen::Vector3d displaceVelocity = moonMomentum.normalized() * 11;

    Eigen::Vector6d obelixInitialState = Eigen::Vector6d();
    obelixInitialState.segment(0,3) = L4Cart;
    obelixInitialState.segment(3,3) = L4initVelocity + displaceVelocity;

    Eigen::Vector6d asterixDisplacement = Eigen::Vector6d();
    asterixDisplacement(0) = 10e3; // 1: 100
    asterixDisplacement(1) = 1e3; // 1: 1e3

    Eigen::Vector6d asterixInitialState = Eigen::Vector6d();
    asterixInitialState.segment(0,3) = L4Cart + asterixDisplacement.segment(0,3);
    asterixInitialState.segment(3,3) = L4initVelocity + displaceVelocity;

    int it = 0;
    int N = 20;
    double range = 100e3;
    double step = 2*range/N;
    for (int x =0; x < N; x++){
        for (int y = 0; y < N; y++){

            std::cout << "Starting iteration number " << it << std::endl;
            it++;
            Eigen::Vector6d tetrisDisplacement = Eigen::Vector6d();
            tetrisDisplacement(0) = -range + x*step;
            tetrisDisplacement(1) = - range + y*step;

            std::cout << "Position tetris" << std::endl;
            Eigen::Vector6d tetrisInitialState = Eigen::Vector6d();
            tetrisInitialState.segment(0,3) = L4Cart + tetrisDisplacement.segment(0,3);
            tetrisInitialState.segment(3,3) = L4initVelocity + displaceVelocity;


            // Set initial state
            Eigen::VectorXd systemInitialState = Eigen::VectorXd( 18 );
            systemInitialState.segment( 0, 6 ) = asterixInitialState;
            systemInitialState.segment( 6, 6 ) = obelixInitialState;
            systemInitialState.segment( 12, 6 ) = tetrisInitialState;


            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch );

            const double minimumStepSize = 0.1;
            const double maximumStepSize = 24*3600;
            const double relativeErrorTolerance = 1E-11;
            const double absoluteErrorTolerance = 1E-11;
            std::shared_ptr< RungeKuttaVariableStepSizeSettings<> > integratorSettings =
                    std::make_shared< RungeKuttaVariableStepSizeSettings <> >
                    ( simulationStartEpoch, 10, RungeKuttaCoefficients::rungeKuttaFehlberg45,minimumStepSize,maximumStepSize,
                      relativeErrorTolerance, absoluteErrorTolerance);


            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            std::cout << "Start propagating" << std::endl;
            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution();


            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            std::string outputSubFolder = "costMapping/";

            attachment = "_" + std::to_string(x) + "_" + std::to_string(y);
            InterpolateAndSave(integrationResult,
                               "propagationHistory" + attachment +".dat",
                               simulationStartEpoch,
                               interpolationN,
                               interpolationStep);

        }
    }


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}

