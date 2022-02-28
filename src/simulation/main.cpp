#include "mpi.h"
#include <Eigen/Geometry>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "General_Utilities.hpp"
#include "IC_Generator.hpp"
#include "PREM.hpp"
#include "Physical_Parameters.hpp"
#include "RN_Generators.hpp"
#include "Trajectory_Simulation.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[])
{

	// INITIALIZATION
	////////////////////////////////////////////////////////////
	// MPI Enviroment
	//  Initialize the MPI environment
	MPI_Init(NULL, NULL);
	// Get the number of processes
	int numprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	// Get the ID number of the process
	int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	// Starting time
	high_resolution_clock::time_point tStart, t1, t2, tEnd;
	double durationIni	= 0.0;
	double durationMain = 0.0;
	if(myRank == 0)
	{
		tStart = high_resolution_clock::now();
	}
	// Read in configuration file, set up the detector and define the earth velocity
	Read_Config_File(argv[1]);
	// Initialize Logfile and create folders for inputfiles
	if(myRank == 0)
	{
		Copy_Config_File(argv[1]);
		std::chrono::time_point<std::chrono::system_clock> start;
		start				   = std::chrono::system_clock::now();
		std::time_t start_time = std::chrono::system_clock::to_time_t(start);
		cout << "\n##############################" << endl
			 << "DaMaSCUS" << version << " - Simulation" << endl
			 << endl
			 << "Starting Time: " << std::ctime(&start_time)
			 << "Simulation ID: " << SimID << endl
			 << endl
			 << "Creating logfile." << endl;
		LogFile_Initialize(numprocs);
	}
	// Initialize Random Number Generator:
	std::random_device rd;
	std::mt19937 PRNG(rd());
	// Synchronize processes:
	MPI_Barrier(MPI_COMM_WORLD);

	////////////////////////////////////////////////////////////

	// 1. Initial MC run without scatterings:
	// Initialize the earth model for initial run
	double sigmaInitial = 1.0e-60 * cm * cm;
	Initialize_PREM(mChi, sigmaInitial);
	// Deactivate formfactor, if it is used:
	string formfactor0 = FormFactor;
	FormFactor		   = "None";
	if(myRank == 0)
		cout << "Start initial MC simulation run without DM scatterings." << endl;
	// Desired Data Sample for inital and main MC simulation
	unsigned long long int Local_SampleSize_Initial = ceil((double) Global_SampleSize_Initial / numprocs);
	Global_SampleSize_Initial						= Local_SampleSize_Initial * numprocs;
	// Particle counter per isodetection ring
	unsigned long long int Global_N0[Isodetection_Rings];
	unsigned long long int Local_N0[Isodetection_Rings];
	for(int i = 0; i < Isodetection_Rings; i++)
		Local_N0[i] = 0;
	long double Global_W0[Isodetection_Rings];
	long double Local_W0[Isodetection_Rings];
	for(int i = 0; i < Isodetection_Rings; i++)
		Local_W0[i] = 0.0;
	long double Global_W0sq[Isodetection_Rings];
	long double Local_W0sq[Isodetection_Rings];
	for(int i = 0; i < Isodetection_Rings; i++)
		Local_W0sq[i] = 0.0;

	// Depth crossing counter
	unsigned long long int Global_Counter_Crossings0 = 0;
	unsigned long long int Local_Counter_Crossings0	 = 0;

	// Simulation of trajectories without scatterings and without formfactor implementation
	for(unsigned int i = 0; i < Local_SampleSize_Initial; i++)
	{
		Event IC						  = InitialCondition(0, vEarth, PRNG);
		Trajectory trajectory			  = ParticleTrack(mChi, sigmaInitial, IC, vcut, PRNG);
		std::vector<Event> CrossingEvents = trajectory.DepthCrossing(Detector_Depth);
		Local_Counter_Crossings0 += CrossingEvents.size();
		// Increment the isodetection counters and record the weighted velocity norm.
		for(unsigned int j = 0; j < CrossingEvents.size(); j++)
		{
			int IsoRing = CrossingEvents[j].IsodetectionRing(Isodetection_Rings);
			Local_N0[IsoRing]++;

			double speed_ratio = IC.NormVelocity() / CrossingEvents[j].NormVelocity();
			double weight	   = CrossingEvents[j].DataWeight(speed_ratio);
			Local_W0[IsoRing] += weight;
			Local_W0sq[IsoRing] += weight * weight;
		}
	}
	// Reactivate form factor
	FormFactor = formfactor0;
	// Reduce the local isodetection counter
	MPI_Reduce(&Local_N0, &Global_N0, Isodetection_Rings, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&Local_Counter_Crossings0, &Global_Counter_Crossings0, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&Local_W0, &Global_W0, Isodetection_Rings, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&Local_W0sq, &Global_W0sq, Isodetection_Rings, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	// Status update on the console
	if(myRank == 0)
	{
		// computing time
		t1			= high_resolution_clock::now();
		durationIni = 1e-6 * duration_cast<microseconds>(t1 - tStart).count();
		cout << "\tInitial run finished\t(" << floor(durationIni) << " s)." << endl
			 << endl;
	}
	// 2. MC run with scatterings and velocity data collection.
	// Initialize the earth model for the main run
	Initialize_PREM(mChi, sigma0);
	if(myRank == 0)
		cout << "Start main MC simulation run with scatterings." << endl
			 << "\tDM mass [MeV]:\t" << mChi / MeV << endl
			 << "\tSigma [cm^2]:\t" << sigma0 / (cm * cm) << endl
			 << "\tForm factor:\t" << FormFactor << endl
			 << "\tMean free path:" << std::endl
			 << "\t\tCore [rEarth]:\t\t" << Mean_Free_Path(0.0, mChi, sigma0, 300.0 * km / sec) / rEarth << std::endl
			 << "\t\tMantle [rEarth]:\t" << Mean_Free_Path(0.8 * rEarth, mChi, sigma0, 300.0 * km / sec) / rEarth << std::endl
			 << std::endl;
	// Particle counter per isodetection ring
	unsigned long long int Global_N[Isodetection_Rings];
	unsigned long long int Local_N[Isodetection_Rings];
	for(int i = 0; i < Isodetection_Rings; i++)
		Local_N[i] = 0;
	// Sum of data weights, necessary for determination of local DM number density
	long double Global_W[Isodetection_Rings];
	long double Local_W[Isodetection_Rings];
	for(int i = 0; i < Isodetection_Rings; i++)
		Local_W[i] = 0.0;
	// Sum of squared data weights, necessary for determination of the local DM number density uncertainty.
	long double Global_Wsq[Isodetection_Rings];
	long double Local_Wsq[Isodetection_Rings];
	for(int i = 0; i < Isodetection_Rings; i++)
		Local_Wsq[i] = 0.0;
	// Work Load Distribution: Desired Velocity Data Sample Size for each isodetection ring for each MPI process.
	unsigned int Local_SampleSize_Velocity = ceil((double) Global_SampleSize_Velocity / numprocs);
	Global_SampleSize_Velocity			   = Local_SampleSize_Velocity * numprocs;
	// Counter for simulated trajectories
	unsigned long long int Global_Counter_Tracks = 0;
	unsigned long long int Local_Counter_Tracks	 = 0;
	// Counter for scattering events
	unsigned long long int Global_Counter_Scatterings = 0;
	unsigned long long int Local_Counter_Scatterings  = 0;
	// Counter for free particles:
	unsigned long long int Global_Counter_Free = 0;
	unsigned long long int Local_Counter_Free  = 0;
	// Counter for depth crossing
	unsigned long long int Global_Counter_Crossings = 0;
	unsigned long long int Local_Counter_Crossings	= 0;
	// Counter for velocity cutoffs
	unsigned long long int Global_Counter_vCutoff = 0;
	unsigned long long int Local_Counter_vCutoff  = 0;
	// Local Counter for velocity data
	unsigned int Local_Counter_DatapointsTotal = 0;
	unsigned int Local_Counter_Datapoints[Isodetection_Rings];
	for(int i = 0; i < Isodetection_Rings; i++)
		Local_Counter_Datapoints[i] = 0;

	// MPI Output File and Offset
	MPI_Offset offset;
	MPI_File file_Velocity[Isodetection_Rings];
	MPI_File file_Weights[Isodetection_Rings];
	MPI_Status status;
	// Velocity Output files including process dependent offset
	for(int i = 0; i < Isodetection_Rings; i++)
	{
		// Velocity output files
		string filename = "../data/" + SimID + "_data/velocity." + std::to_string(i);
		int test		= MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(filename.c_str()), MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_Velocity[i]);
		// if it already exists, it will be overwritten!
		if(test != MPI_SUCCESS)
		{
			if(myRank == 0)
				MPI_File_delete(const_cast<char*>(filename.c_str()), MPI_INFO_NULL);
			// MPI_Barrier(MPI_COMM_WORLD);
			test = MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(filename.c_str()), MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_Velocity[i]);
		}
		// Offset
		offset = myRank * Local_SampleSize_Velocity * sizeof(Eigen::Vector3d);
		MPI_File_seek(file_Velocity[i], offset, MPI_SEEK_SET);
		// Weights output files
		filename = "../data/" + SimID + "_data/weights." + std::to_string(i);
		test	 = MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(filename.c_str()), MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_Weights[i]);
		// if it already exists, it will be overwritten!
		if(test != MPI_SUCCESS)
		{
			if(myRank == 0)
				MPI_File_delete(const_cast<char*>(filename.c_str()), MPI_INFO_NULL);
			// MPI_Barrier(MPI_COMM_WORLD);
			test = MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(filename.c_str()), MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_Weights[i]);
		}
		// Offset
		offset = myRank * Local_SampleSize_Velocity * sizeof(double);
		MPI_File_seek(file_Weights[i], offset, MPI_SEEK_SET);
	}
	// Simulation of Tracks
	while(Local_Counter_DatapointsTotal < Isodetection_Rings * Local_SampleSize_Velocity)
	{
		// Simulate track.
		Event IC			  = InitialCondition(0, vEarth, PRNG);
		Trajectory trajectory = ParticleTrack(mChi, sigma0, IC, vcut, PRNG);
		// Increase counter of tracks and scatterings and vcutoff reachers
		Local_Counter_Tracks++;
		int scatterings = trajectory.NoOfScatterings();
		if(scatterings == 0)
			Local_Counter_Free++;
		else
			Local_Counter_Scatterings += scatterings;
		if(trajectory.Trajectory_Type() == 2)
			Local_Counter_vCutoff++;
		// Where did the particle pass the isodetection rings.
		std::vector<Event> CrossingEvents = trajectory.DepthCrossing(Detector_Depth);
		Local_Counter_Crossings += CrossingEvents.size();
		for(unsigned int j = 0; j < CrossingEvents.size(); j++)
		{
			// Which ring?
			int IsoRing = CrossingEvents[j].IsodetectionRing(Isodetection_Rings);
			// Increase the weight sums
			double speed_ratio = IC.NormVelocity() / CrossingEvents[j].NormVelocity();
			double weight	   = CrossingEvents[j].DataWeight(speed_ratio);
			Local_W[IsoRing] += weight;
			Local_Wsq[IsoRing] += weight * weight;
			// Increase the particle counter.
			Local_N[IsoRing]++;
			// If this isodetection ring is not finished collecting data, save the velocity and weight
			if(Local_Counter_Datapoints[IsoRing] < Local_SampleSize_Velocity)
			{
				// Save velocity.
				Eigen::Vector3d vel = CrossingEvents[j].Velocity();
				// Increase data counters
				Local_Counter_Datapoints[IsoRing]++;
				Local_Counter_DatapointsTotal++;
				// Save velocity and weights to file
				MPI_File_write(file_Velocity[IsoRing], vel.data(), vel.size(), MPI_DOUBLE, &status);
				MPI_File_write(file_Weights[IsoRing], &weight, 1, MPI_DOUBLE, &status);
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// MPI Reductions
	MPI_Reduce(&Local_Counter_Tracks, &Global_Counter_Tracks, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&Local_Counter_Scatterings, &Global_Counter_Scatterings, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&Local_Counter_vCutoff, &Global_Counter_vCutoff, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&Local_Counter_Free, &Global_Counter_Free, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&Local_Counter_Crossings, &Global_Counter_Crossings, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&Local_N, &Global_N, Isodetection_Rings, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&Local_W, &Global_W, Isodetection_Rings, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&Local_Wsq, &Global_Wsq, Isodetection_Rings, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	// Calculate DM Densities
	vector<vector<double>> Edensity = DM_EnergyDensity(Global_W0, Global_SampleSize_Initial, Global_W, Global_Counter_Tracks, Global_W0sq, Global_Wsq, Isodetection_Rings);
	vector<vector<double>> Ndensity = DM_NumberDensity(mChi, Edensity);
	double averageNdensity			= DM_AverageDensity(Ndensity, Isodetection_Rings, Detector_Depth);
	// Status update and save N(theta) and rho(theta)
	if(myRank == 0)
	{
		// Save density
		ofstream f;
		f.open("../data/" + SimID + ".rho");
		for(int i = 0; i < Isodetection_Rings; i++)
			f << 180.0 / Isodetection_Rings * i << "\t" << InUnits(Edensity[i][0], GeV / cm / cm / cm) << "\t" << InUnits(Edensity[i][1], GeV / cm / cm / cm) << endl;	 //<<"\t" <<InUnits(AverageVelocity_0[i],km/sec)<<"\t" <<InUnits(AverageVelocity[i],km/sec)<<endl;// <<"\t" <<Global_N0[i]<<"\t" <<Global_N[i] <<endl;
		f.close();
		// computing time
		t2			 = high_resolution_clock::now();
		durationMain = 1e-6 * duration_cast<microseconds>(t2 - t1).count();
		cout << "Main MC run finished\t"
			 << "(" << floor(durationMain) << " s)." << endl;
	}
	// Close all output files
	for(int i = 0; i < Isodetection_Rings; i++)
	{
		MPI_File_close(&file_Velocity[i]);
		MPI_File_close(&file_Weights[i]);
	}
	////////////////////////////////////////////////////////////

	// Finalize simulation
	////////////////////////////////////////////////////////////
	// Master process finishes logfile
	if(myRank == 0)
	{
		// Ending time and computing time
		tEnd				 = high_resolution_clock::now();
		double durationTotal = 1e-6 * duration_cast<microseconds>(tEnd - tStart).count();
		cout << "\nProcessing Time:\t" << durationTotal << "s (" << floor(durationTotal / 3600.0) << ":" << floor(fmod(durationTotal / 60.0, 60.0)) << ":" << floor(fmod(durationTotal, 60.0)) << ":" << floor(fmod(1000 * durationTotal, 1000.0)) << ")." << endl
			 << "##############################" << endl;
		// Finalize Logfile
		LogFile_Finalize(Global_Counter_Tracks, Global_Counter_Scatterings, Global_Counter_Free, Global_Counter_vCutoff, Global_SampleSize_Initial, Global_SampleSize_Velocity, Global_Counter_Crossings0, Global_Counter_Crossings, averageNdensity, durationIni, durationMain, durationTotal);
	}
	// Finalize the MPI environment.
	MPI_Finalize();
	return 0;
	////////////////////////////////////////////////////////////
}
