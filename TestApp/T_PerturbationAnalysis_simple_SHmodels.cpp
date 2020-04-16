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
    //const string kernelfile = input_output::getSpiceKernelPath() + "tudat_merged_edit.bsp";
    std::vector< std::string > kernelfile;
    kernelfile.push_back(input_output::getSpiceKernelPath() + "tudat_merged_edit.bsp");
    spice_interface::loadStandardSpiceKernels(kernelfile);

    // Set simulation time settings.
    const double simulationStartEpoch = 30*tudat::physical_constants::JULIAN_YEAR; // 2030 start date
    const double simulationEndEpoch = simulationStartEpoch +  1*tudat::physical_constants::JULIAN_YEAR; //
    const double timestep = 30*60;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Venus" );
    bodiesToCreate.push_back("Jupiter");
    bodiesToCreate.push_back("Saturn");
    bodiesToCreate.push_back("Neptune");
    bodiesToCreate.push_back("Uranus");
    bodiesToCreate.push_back("Mercury");


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

    // Create radiation pressure settings
    referenceAreaRadiation = tudat::mathematical_constants::PI * (1734E3)*(1734E3);
    radiationPressureCoefficient = 1.14; // Moon albedo is about 0.14 - wiki
    std::shared_ptr< RadiationPressureInterfaceSettings > moonRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Moon" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    moonRadiationPressureSettings, "Moon", bodyMap ));

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

    const int numcases = 10;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////                   CASES                  //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// 0: Default full model
    /// 1: Remove Mercury
    /// 2: Remove Mars
    /// 3: Remove Venus
    /// 4: Remove Sun
    /// 5: Remove Radiation Pressure
    /// 6: Remove Jupiter
    /// 7: Remove Saturn
    /// 8: Remove Neptune
    /// 9: Remove Uranus
    /// 10: Remove Spherical Harmonics
    /// 11: SH: 10,10
    /// 12: SH: 20,20
    /// 13: SHM: 5,5
    /// 14: SHM: 10,10
    /// 15: SHM: 20,20
    ///
    ///

    for (int j = 0; j <=2; j++){
        for (int i = 1; i <= numcases; i++){
            const string nameAttach = "_"+ std::to_string(j) + "_" + std::to_string(i);

            int D = 2*i;
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSats;
            switch(j){
            case 0:{
                // Define propagation settings.
                accelerationsOfSats[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                              basic_astrodynamics::central_gravity )  );
                accelerationsOfSats[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                               basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Mercury" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Venus" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats["Saturn"].push_back( std::make_shared< AccelerationSettings> (
                                                             basic_astrodynamics::central_gravity));
                accelerationsOfSats["Neptune"].push_back( std::make_shared< AccelerationSettings> (
                                                              basic_astrodynamics::central_gravity));
                accelerationsOfSats["Uranus"].push_back( std::make_shared< AccelerationSettings >(
                                                             basic_astrodynamics::central_gravity));
                accelerationsOfSats[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::cannon_ball_radiation_pressure ) );

                break;
            }
            case 1:{
                // Define propagation settings.
                accelerationsOfSats[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( D, D ) );
                accelerationsOfSats[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                               basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Mercury" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Venus" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats["Saturn"].push_back( std::make_shared< AccelerationSettings> (
                                                             basic_astrodynamics::central_gravity));
                accelerationsOfSats["Neptune"].push_back( std::make_shared< AccelerationSettings> (
                                                              basic_astrodynamics::central_gravity));
                accelerationsOfSats["Uranus"].push_back( std::make_shared< AccelerationSettings >(
                                                             basic_astrodynamics::central_gravity));
                accelerationsOfSats[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::cannon_ball_radiation_pressure ) );

                break;
            }
            case 2:{
                // Define propagation settings.
                accelerationsOfSats[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );
                accelerationsOfSats[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                               basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Moon" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( D, D ) );
                accelerationsOfSats[ "Mercury" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Venus" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::central_gravity ) );
                accelerationsOfSats["Saturn"].push_back( std::make_shared< AccelerationSettings> (
                                                             basic_astrodynamics::central_gravity));
                accelerationsOfSats["Neptune"].push_back( std::make_shared< AccelerationSettings> (
                                                              basic_astrodynamics::central_gravity));
                accelerationsOfSats["Uranus"].push_back( std::make_shared< AccelerationSettings >(
                                                             basic_astrodynamics::central_gravity));
                accelerationsOfSats[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 basic_astrodynamics::cannon_ball_radiation_pressure ) );

                break;
            }

            }


            accelerationMap[ "Asterix" ] = accelerationsOfSats;
            accelerationMap[ "Obelix" ] = accelerationsOfSats;

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


            // Save name settings
            const string attachment = "_100km_SH";

            Eigen::Vector6d asterixDisplacement = Eigen::Vector6d();
            asterixDisplacement(0) = 100e3; // 1: 100
            asterixDisplacement(1) = 0e3; // 1: 1e3

            Eigen::Vector6d asterixInitialState = Eigen::Vector6d();
            asterixInitialState.segment(0,3) = L4Cart + asterixDisplacement.segment(0,3);
            asterixInitialState.segment(3,3) = L4initVelocity + displaceVelocity;


            // Set initial state
            Eigen::VectorXd systemInitialState = Eigen::VectorXd( 12 );
            systemInitialState.segment( 0, 6 ) = asterixInitialState;
            systemInitialState.segment( 6, 6 ) = obelixInitialState;

            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch );


            std::shared_ptr< IntegratorSettings< > > integratorSettings =
                    std::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, simulationStartEpoch, timestep );
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution();


            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            std::string outputSubFolder = "PA_SH/";

            std::cout << "Finished simulation run " << std::to_string(i) << " of " << std::to_string(numcases) << std::endl;

            // Write perturbed satellite propagation history to file.
            input_output::writeDataMapToTextFile( integrationResult,
                                                  "propagationHistory"+nameAttach + attachment +".dat",
                                                  tudat_applications::getOutputPath( ) + outputSubFolder,
                                                  "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );

            if( j ==0){break;}
        }

    }





    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
