#include "General_Utilities.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include <libconfig.h++>

#include "Physical_Parameters.hpp"

using namespace libconfig;
	

//Config file
	string FormFactor="None";
	string SimID="default";
	string experiment="default";
	string version="v1.0";
	double mChi=0,sigma=0;
	double rhoDM=0.0;
	int SampleSize=0;
	int SimDate [3] = 	{0,0,0};
	int SimTime [3] = 	{0,0,0};
	double nJ2000 = 0.0;
	Eigen::Vector3d vEarth(0,0,0);
	double vcut=0.0;
	void Read_Config_File(char *inputfile)
	{
		Config cfg;
		std::string id(inputfile);
		std::string path = "../data/"+id+".cfg";
		const char *cstr = path.c_str();
		// Read the file. If there is an error, report it and exit.
			try
			{
			  cfg.readFile(cstr);
			}
			catch(const FileIOException &fioex)
			{
				std::cerr << "I/O error while reading configuration file." << std::endl;
				exit(EXIT_FAILURE);
			}
			catch(const ParseException &pex)
			{
				std::cerr << "Configurate file parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
				exit(EXIT_FAILURE);
			}
		//Simulation ID
			try
			{
				SimID = cfg.lookup("simID").c_str();
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'simID' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Experiment
			try
			{
				experiment = cfg.lookup("experiment").c_str();
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'simID' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//sample size
			try
			{
				SampleSize = cfg.lookup("samplesize");
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'samplesize' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Mass
			try
			{
				mChi = cfg.lookup("mass");
				mChi*=MeV;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'mass' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Cross section
			try
			{
				sigma = cfg.lookup("sigma");
				sigma*=pb;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "Error: While reading 'sigma' setting in configuration file. Computation cancelled." << endl;
				exit(EXIT_FAILURE);
			}
		//DM energy density
			try
			{
				rhoDM = cfg.lookup("rho");
				rhoDM*=GeV/cm/cm/cm;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "Error: While reading 'rho' setting in configuration file. Computation cancelled." << endl;
				exit(EXIT_FAILURE);
			}
		//Date
			try
			{
				SimDate[0] = cfg.lookup("date")[0];
				SimDate[1] = cfg.lookup("date")[1];
				SimDate[2] = cfg.lookup("date")[2];
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "Error: While reading 'date' setting in configuration file. Computation cancelled." << endl;
				exit(EXIT_FAILURE);
			}
		//Time
			try
			{
				SimTime[0] = cfg.lookup("time")[0];
				SimTime[1] = cfg.lookup("time")[1];
				SimTime[2] = cfg.lookup("time")[2];
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "Error: While reading 'time' setting in configuration file. Computation cancelled." << endl;
				exit(EXIT_FAILURE);
			}
		//Time and Detector
			nJ2000= FractionalDays(SimDate,SimTime);
		//Earth Velocity
			vEarth=EarthVelocity(nJ2000);
	}

	void LogAnalysis(double duration,int worldsize)
	{
			ofstream f;
			f.open("../results/"+SimID+".log",std::ofstream::app);
			std::chrono::time_point<std::chrono::system_clock> end;
			end = std::chrono::system_clock::now();
			std::time_t end_time = std::chrono::system_clock::to_time_t(end);
			f	<<"\n////////////////////////////////////////////////////\n\n"
			<<"//Data analysis performed:\t" <<std::ctime(&end_time)<<endl
			<<"\tMPI processes:\t\t" <<worldsize<<endl
			<<"\tExperiment:\t\t" <<experiment<<endl
			<<"\tComputation time:\t" <<duration <<"\t("<<floor(duration/3600.0)<<":"<<floor(fmod(duration/60.0,60.0))<<":"<<floor(fmod(duration,60.0))<<":"<<floor(fmod(1000*duration,1000.0)) <<")"<<endl;
			
	}

	std::vector<int> WorkDivision(int WorldSize)
	{
		std::vector<int> output;
		int overlap = WorldSize*ceil(180.0/WorldSize)-180;
		for(int i =0 ; i<=WorldSize; i++)
		{
			if(i==0) 			output.push_back(0);
			else if(i<=overlap)	output.push_back(i*floor(180.0/WorldSize));
			else				output.push_back(output[i-1]+ceil(180.0/WorldSize));
		}
		return output;
	}