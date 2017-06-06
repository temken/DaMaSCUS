#ifndef __General_Utilities_hpp_
#define __General_Utilities_hpp_

#include <Eigen/Geometry>
#include <string>
#include <vector>

//Config File
	//Read config file and Input Parameters
		extern std::string FormFactor,version;
		extern std::string SimID,experiment;
		extern double mChi,sigma,nJ2000,rhoDM;
		extern int SimDate [3];
		extern int SimTime [3];
		extern Eigen::Vector3d vEarth;
		extern int SampleSize;
		
	//Read config file and define Input Parameters. Also calculate fractional days, earth velocity.
		extern void Read_Config_File(char *inputfile);

//Create an entry in the log file after successful data analysis
	extern void LogAnalysis(double duration,int worldsize);

//Divide up the isodetection rings between the MPI processes as evenly as possible
	std::vector<int> WorkDivision(int WorldSize);
#endif