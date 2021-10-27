#ifndef __General_Utilities_hpp_
#define __General_Utilities_hpp_

#include <Eigen/Geometry>
#include <fstream>
#include <string>
#include <vector>

// Input Parameter from Config file
extern std::string SimID, FormFactor, version, experiment;
extern double vcut;
extern double mChi, sigma0, nJ2000, rhoDM;
extern double Detector_Depth;
extern unsigned long long int Global_SampleSize_Initial;
extern unsigned int Global_SampleSize_Velocity;
extern int Isodetection_Rings;
extern int SimDate[3];
extern int SimTime[3];
extern Eigen::Vector3d vEarth;
extern int SampleSize;

// Read config file and define Input Parameters
extern void Read_Config_File(const char* inputfile);
extern void Copy_Config_File(const char* inputfile);
// For the simulator
// LogFile
// Global Ofstream for logfile
extern std::ofstream LogStream;
// Create and write stuff into the logfile
extern void LogFile_Initialize(int numprocs, std::ofstream& f = LogStream);
extern void LogFile_Detector(std::ofstream& f = LogStream);
extern void LogFile_Finalize(unsigned long long int nParticles, unsigned long long int nScatterings, unsigned long long int nFree, unsigned long long int vcutoff, unsigned long long int nIni, unsigned int nData, unsigned long long int nCrossing0, unsigned long long int nCrossing, double AvDensity, double duration1, double duration2, double duration3, std::ofstream& f = LogStream);

// DM Density
extern double IsoDetectionRing_Area(int theta, int rings, double depth = Detector_Depth);
extern std::vector<std::vector<double>> DM_EnergyDensity(long double w0[], int long long unsigned n0total, long double w[], int unsigned long long ntotal, long double w0sq[], long double wsq[], int rings);
extern std::vector<std::vector<double>> DM_NumberDensity(double mDM, std::vector<std::vector<double>> density);
extern double DM_AverageDensity(std::vector<std::vector<double>> density, int rings, double depth = Detector_Depth);

// For the analyzer
// Create an entry in the log file after successful data analysis
extern void LogFile_Analysis(double duration, int worldsize);

// Divide up the isodetection rings between the MPI processes as evenly as possible
std::vector<int> WorkDivision(int WorldSize, int rings);
#endif