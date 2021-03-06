/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

/* Used to study the effect of changes in body pointing on solar radiation pressure and orbit evolution */
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "TestApp/applicationOutput.h"

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

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 365*tudat::physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Venus" );

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared < propagators::SingleDependentVariableSaveSettings> (
                                          propagators::keplerian_state_dependent_variable, "Moon", "Earth"));
    // Create object with list of dependent variables.
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    for (int i = 0; i<=1; i++){

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create spacecraft object.
        bodyMap[ "Asterix" ] = std::make_shared< simulation_setup::Body >( );
        bodyMap[ "Asterix" ]->setConstantBodyMass( 5.0 );
        bodyMap[ "Obelix" ] = std::make_shared< simulation_setup::Body >();
        bodyMap[ "Obelix" ]->setConstantBodyMass(5.0);

//        // Create aerodynamic coefficient interface settings.
//        double referenceArea = 0.32;
//        double aerodynamicCoefficient = 1.2;
//        std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
//                std::make_shared< ConstantAerodynamicCoefficientSettings >(
//                    referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

//        // Create and set aerodynamic coefficients object
//        bodyMap[ "Asterix" ]->setAerodynamicCoefficientInterface(
//                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Asterix" ) );

//        bodyMap[ "Obelix"]->setAerodynamicCoefficientInterface(
//                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Obelix"));

        // Create radiation pressure settings
        double referenceAreaRadiation = 0.32;
        double radiationPressureCoefficient = 1.300;
        std::vector< std::string > occultingBodies;
        occultingBodies.push_back( "Earth" );
        std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
                std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

        // Create radiation pressure settings
        referenceAreaRadiation = 0.3416;
        radiationPressureCoefficient = 1.1421;
        std::shared_ptr< RadiationPressureInterfaceSettings > obelixRadiationPressureSettings =
                std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

        // Create and set radiation pressure settings
        bodyMap[ "Asterix" ]->setRadiationPressureInterface(
                    "Sun", createRadiationPressureInterface(
                        asterixRadiationPressureSettings, "Asterix", bodyMap ));
        if (i==0){
            bodyMap["Obelix"]->setRadiationPressureInterface(
                        "Sun", createRadiationPressureInterface(
                            asterixRadiationPressureSettings,"Obelix", bodyMap));
        }
        else{
            bodyMap[ "Obelix" ]->setRadiationPressureInterface(
                        "Sun", createRadiationPressureInterface(
                            obelixRadiationPressureSettings, "Obelix", bodyMap));
        }


        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define propagation settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;

        accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

        accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::central_gravity ) );
        accelerationsOfAsterix[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
        accelerationsOfAsterix[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
        accelerationsOfAsterix[ "Venus" ].push_back( std::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );

        accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::cannon_ball_radiation_pressure ) );

        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
        accelerationMap[ "Obelix" ] = accelerationsOfAsterix;
        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth" );
        bodiesToPropagate.push_back("Obelix");
        centralBodies.push_back("Earth");

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
        asterixDisplacement(0) = 0;
        asterixDisplacement(1) = 1e3;

        Eigen::Vector6d asterixInitialState = Eigen::Vector6d();
        asterixInitialState.segment(0,3) = L4Cart + asterixDisplacement.segment(0,3);
        asterixInitialState.segment(3,3) = L4initVelocity + displaceVelocity;


        // Set initial state
        Eigen::VectorXd systemInitialState = Eigen::VectorXd( 12 );
        systemInitialState.segment( 0, 6 ) = asterixInitialState;
        systemInitialState.segment( 6, 6 ) = obelixInitialState;

        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch,
                  propagators::cowell, dependentVariablesToSave  );

        const double fixedStepSize = 30.0*60.0;
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, fixedStepSize );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::string outputSubFolder = "PointingAnalysis/";

        Eigen::VectorXd finalIntegratedState = (--integrationResult.end( ) )->second;
        // Print the position (in km) and the velocity (in km/s) at t = 0.
        std::cout << "Single Earth-Orbiting Satellite Example." << std::endl <<
                     "The initial position vector of Asterix is [km]:" << std::endl <<
                     asterixInitialState.segment( 0, 3 ) / 1E3 << std::endl <<
                     "The initial velocity vector of Asterix is [km/s]:" << std::endl <<
                     asterixInitialState.segment( 3, 3 ) / 1E3 << std::endl;

        // Print the position (in km) and the velocity (in km/s) at t = 86400.
        std::cout << "After " << simulationEndEpoch <<
                     " seconds, the position vector of Asterix is [km]:" << std::endl <<
                     finalIntegratedState.segment( 0, 3 ) / 1E3 << std::endl <<
                     "And the velocity vector of Asterix is [km/s]:" << std::endl <<
                     finalIntegratedState.segment( 3, 3 ) / 1E3 << std::endl;

        const string nameAttach = "_"+ std::to_string(i);
        // Write perturbed satellite propagation history to file.
        input_output::writeDataMapToTextFile( integrationResult,
                                              "propagationHistory" + nameAttach + ".dat",
                                              tudat_applications::getOutputPath( ) + outputSubFolder,
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        std::map<double, Eigen::VectorXd> lunarkeplerMap = dynamicsSimulator.getDependentVariableHistory();
        double mu = bodyMap["Earth"]->getGravityFieldModel()->getGravitationalParameter();
        std::map< double, Eigen::VectorXd > moonstate;
        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = lunarkeplerMap.begin( );
                 stateIterator != lunarkeplerMap.end( ); stateIterator++ ){
            Eigen::Vector6d grabber = stateIterator->second;
            moonstate.insert(std::pair<double,Eigen::VectorXd>(stateIterator->first , convertKeplerianToCartesianElements(grabber,mu) ));
        }

        // Write lunar state history to file.
        input_output::writeDataMapToTextFile( moonstate,
                                              "propagationHistory_moon.dat",
                                              tudat_applications::getOutputPath( ) + outputSubFolder,
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );
    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
