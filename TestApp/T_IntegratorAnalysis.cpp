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

void InterpolateAndSave(std::map< double, Eigen::VectorXd > integrationResult, std::string name, double startTime,
                        int N = 365, double stepsize = tudat::physical_constants::JULIAN_DAY){
    using namespace tudat;
    std::string outputSubFolder = "IntegratorAnalysis/";

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

//    // Write perturbed satellite propagation history to file.
//    input_output::writeDataMapToTextFile( statemap,
//                                          "original_"+ name + ".dat",
//                                          tudat_applications::getOutputPath( ) + outputSubFolder ,
//                                          "",
//                                          std::numeric_limits< double >::digits10,
//                                          std::numeric_limits< double >::digits10,
//                                          "," );
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

    using namespace tudat_applications;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Save name settings

    const string attachment = "_365d";
    const string basename = "interpolatedHistory";



    // Export interesting data
    // Vector of function evaluations and times
    std::vector< unsigned int > numberOfFunctionEvaluations;
    std::vector< double > totalPropagationTime;
    std::vector< string > labels;

    // Load Spice kernels.
    //const string kernelfile = input_output::getSpiceKernelPath() + "tudat_merged_edit.bsp";
    std::vector< std::string > kernelfile;
    kernelfile.push_back(input_output::getSpiceKernelPath() + "tudat_merged_edit.bsp");
    spice_interface::loadStandardSpiceKernels(kernelfile);

    // Set simulation time settings.
    const double simulationDuration = 365*tudat::physical_constants::JULIAN_DAY;
    const double simulationStartEpoch = 30*tudat::physical_constants::JULIAN_YEAR; // 2030 start date
    const double simulationEndEpoch = simulationStartEpoch + simulationDuration; //
    double timestep = 30*60;

    // interpolation settings

    const double interpolationStep  = 4*3600;
    const int interpolationN = int(simulationDuration/interpolationStep);

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
    asterixDisplacement(0) = 100e3; // 1: 100
    asterixDisplacement(1) = 0; // 1: 1e3

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

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// TOLERANCE - AFFIX
    /// E6: 1E-6
    /// E7: 1E-7
    /// E8: 1E-8
    /// E9: 1E-9
    /// E10: 1E-10
    /// E11: 1E-11
    /// E12: 1E-12

    std::vector<std::string> MethodNames = {"RK4 10s", "RK4 1M", "RK4 5M", "RK4 15M", "RKF45", "RKF56", "RKF78", "DOPRI87", "ABM", "BS"};
    std::vector<std::string> ToleranceNames = {"E6", "E7", "E8", "E9", "E10", "E11", "E12"};
    std::vector<float> Tolerance = {1E-6, 1E-7, 1E-8, 1E-9, 1E-10, 1E-11, 1E-12 };

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// PREFIX
    /// 0: default / control
    /// 1: RK4 1M
    /// 2: RK4 5M
    /// 3: RK4 15M
    /// 4: RKF45
    /// 5: RKF56
    /// 6: RKF78
    /// 7: DOPRI 87
    /// 8: ABM
    /// 9: BS
    for (int i =0; i <= 9; i++){

        std::cout << "Starting run " + std::to_string(i) + " out of 9" << std::endl;
        double minimumStepSize = 1; // 0.1? Caused ABM to hang. Alternatively use a terminationCondition on CPU time (300s?)
        double maximumStepSize = 24*3600;
        double relativeErrorTolerance = 1E-10;
        double absoluteErrorTolerance = 1E-10;
        double initialTimestep = 10;
        timestep = 10;

        if (i <= 3){
            std::shared_ptr< IntegratorSettings< > > integratorSettings;
            switch(i){
            case 0:{
                integratorSettings = std::make_shared< IntegratorSettings< > >
                                    ( rungeKutta4, simulationStartEpoch, timestep );
                break;
            }
            case 1:{
                integratorSettings = std::make_shared< IntegratorSettings< > >
                                    ( rungeKutta4, simulationStartEpoch, 60 );
                break;
            }
            case 2:{
                integratorSettings = std::make_shared< IntegratorSettings< > >
                                    ( rungeKutta4, simulationStartEpoch, 5*60 );
                break;
            }
            case 3:{
                integratorSettings = std::make_shared< IntegratorSettings< > >
                                    ( rungeKutta4, simulationStartEpoch, 15*60 );
                break;
            }
            default:{
                integratorSettings = std::make_shared< IntegratorSettings< > >
                                    ( rungeKutta4, simulationStartEpoch, timestep );
                std::cout << "THIS SHOULD NOT BE VISIBLE - LOOP 1 DEFAULT" << std::endl;
                break;
            }
            }


            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
            numberOfFunctionEvaluations.push_back( dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second );
            totalPropagationTime.push_back( dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second );
            labels.push_back( MethodNames[i] );
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//            std::string outputSubFolder = "IntegratorAnalysis/";

            std::string name = basename +"_" + std::to_string(i) + attachment;
            InterpolateAndSave(integrationResult, name, simulationStartEpoch,interpolationN,interpolationStep);

        }
        else{
            for (int j = 0; j < Tolerance.size(); j++){

                std::cout << "Starting subroutine " + std::to_string(j+1) + " out of " + std::to_string(Tolerance.size()) << std::endl;
                relativeErrorTolerance = Tolerance[j];
                absoluteErrorTolerance = Tolerance[j];

                std::map< double, Eigen::VectorXd > integrationResult;
                switch(i){
                case 4:{
                    std::shared_ptr< RungeKuttaVariableStepSizeSettings<> > integratorSettings =
                        std::make_shared< RungeKuttaVariableStepSizeSettings <> >
                        ( simulationStartEpoch, initialTimestep, RungeKuttaCoefficients::rungeKuttaFehlberg45,minimumStepSize,maximumStepSize,
                          relativeErrorTolerance, absoluteErrorTolerance);

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    // Create simulation object and propagate dynamics.
                    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
                    integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution();

                    numberOfFunctionEvaluations.push_back( dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second );
                    totalPropagationTime.push_back( dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second );
                    labels.push_back( MethodNames[i] + ToleranceNames[j] );
                    break;
                }
                case 5:{
                    std::shared_ptr< RungeKuttaVariableStepSizeSettings<> > integratorSettings =
                        std::make_shared< RungeKuttaVariableStepSizeSettings <> >
                        ( simulationStartEpoch, initialTimestep, RungeKuttaCoefficients::rungeKuttaFehlberg56,minimumStepSize,maximumStepSize,
                          relativeErrorTolerance, absoluteErrorTolerance);

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    // Create simulation object and propagate dynamics.
                    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
                    integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
                    numberOfFunctionEvaluations.push_back( dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second );
                    totalPropagationTime.push_back( dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second );
                    labels.push_back( MethodNames[i] + ToleranceNames[j] );
                    break;
                }
                case 6:{
                    std::shared_ptr< RungeKuttaVariableStepSizeSettings<> > integratorSettings =
                        std::make_shared< RungeKuttaVariableStepSizeSettings <> >
                        ( simulationStartEpoch, initialTimestep, RungeKuttaCoefficients::rungeKuttaFehlberg78,minimumStepSize,maximumStepSize,
                          relativeErrorTolerance, absoluteErrorTolerance);

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    // Create simulation object and propagate dynamics.
                    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
                    integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
                    numberOfFunctionEvaluations.push_back( dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second );
                    totalPropagationTime.push_back( dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second );
                    labels.push_back( MethodNames[i] + ToleranceNames[j] );
                    break;
                }
                case 7:{
                    std::shared_ptr< RungeKuttaVariableStepSizeSettings<> > integratorSettings =
                        std::make_shared< RungeKuttaVariableStepSizeSettings <> >
                        ( simulationStartEpoch, initialTimestep, RungeKuttaCoefficients::rungeKutta87DormandPrince,minimumStepSize,maximumStepSize,
                          relativeErrorTolerance, absoluteErrorTolerance);

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    // Create simulation object and propagate dynamics.
                    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
                    integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
                    numberOfFunctionEvaluations.push_back( dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second );
                    totalPropagationTime.push_back( dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second );
                    labels.push_back( MethodNames[i] + ToleranceNames[j] );
                    break;
                }
                case 8:{
                    std::shared_ptr< AdamsBashforthMoultonSettings<> > integratorSettings =
                        std::make_shared< AdamsBashforthMoultonSettings <> >
                        ( simulationStartEpoch, initialTimestep ,minimumStepSize,maximumStepSize,
                          relativeErrorTolerance, absoluteErrorTolerance, 6, 11);

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    // Create simulation object and propagate dynamics.
                    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
                    integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
                    numberOfFunctionEvaluations.push_back( dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second );
                    totalPropagationTime.push_back( dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second );
                    labels.push_back( MethodNames[i] + ToleranceNames[j] );
                    break;
                }
                case 9:{
                    std::shared_ptr< BulirschStoerIntegratorSettings<> > integratorSettings =
                        std::make_shared< BulirschStoerIntegratorSettings <> >
                        ( simulationStartEpoch, initialTimestep, ExtrapolationMethodStepSequences::bulirsch_stoer_sequence,8,minimumStepSize,maximumStepSize,
                          relativeErrorTolerance, absoluteErrorTolerance);

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    // Create simulation object and propagate dynamics.
                    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
                    integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
                    numberOfFunctionEvaluations.push_back( dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second );
                    totalPropagationTime.push_back( dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second );
                    labels.push_back( MethodNames[i] + ToleranceNames[j] );
                    break;
                }
                default:{
                    std::cout << "THIS SHOULD NOT BE HERE!" << std::endl;
                }
                }







                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                std::string name = basename +"_" + std::to_string(i) + attachment + "_" + ToleranceNames[j];

                InterpolateAndSave(integrationResult, name, simulationStartEpoch,interpolationN,interpolationStep);



        }
    }

    }
    // Write function evaluations and times to file
    tudat::input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( numberOfFunctionEvaluations ),
                       "functionEvaluations"+attachment+".dat", 16, getOutputPath("IntegratorAnalysis/") );
    tudat::input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( totalPropagationTime ),
                       "propagationTime"+attachment+".dat", 16, getOutputPath( "IntegratorAnalysis/" ) );

    //tudat::input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( labels ),
    //                   "propagationTime"+attachment+".dat", 16, getOutputPath( "IntegratorAnalysis/" ) );
    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}


