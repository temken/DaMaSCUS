#ifndef __General_Utilities_hpp_
#define __General_Utilities_hpp_

#include <fstream>
#include <Eigen/Geometry>
#include <string>

//Input Parameter
	extern std::string SimID,FormFactor,version;
	extern double vcut;
	extern double mChi,sigma0,nJ2000,rhoDM;
	extern double Detector_Depth;
	extern unsigned long long int Global_SampleSize_Initial;		
	extern unsigned int Global_SampleSize_Velocity;
	extern int SimDate [3];
	extern int SimTime [3];

	extern Eigen::Vector3d vEarth;


	//Read config file and define Input Parameters
	 	extern void Read_Config_File(char *inputfile);

//LogFile
	//Global Ofstream for logfile
		extern std::ofstream LogStream;
	//Create and write stuff into the logfile
		extern void LogFile_Initialize(int numprocs, std::ofstream& f=LogStream);
		extern void LogFile_Detector(std::ofstream& f=LogStream);
		extern void LogFile_Finalize(unsigned long long int nParticles,unsigned long long int nScatterings,unsigned long long int nFree,unsigned long long int vcutoff,unsigned long long int nIni,unsigned int nData,unsigned long long int nCrossing0,unsigned long long int nCrossing,double AvDensity,double duration1,double duration2,double duration3,std::ofstream& f=LogStream);


#endif