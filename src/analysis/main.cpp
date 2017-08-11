#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <Eigen/Geometry>
#include "mpi.h"

//for creating folders
#include <sys/types.h> // required for stat.h
#include <sys/stat.h>	//required to create a folder

#include "General_Utilities.hpp"
#include "Physical_Parameters.hpp"
#include "ER_General.hpp"
#include "ER_CRESSTII.hpp"
#include "ER_LUX.hpp"

using namespace std;
using namespace std::chrono;



int main(int argc, char *argv[])
{	
	
//INITIALIZATION 	
////////////////////////////////////////////////////////////
	//Read config file created by the simulation code. argv is the SimID
		std::string id(argv[1]);
		std::string path = "../data/"+id+".cfg";
		Read_Config_File(path.c_str());
	//MPI Enviroment
		// Initialize the MPI environment
	    	MPI_Init(NULL, NULL);
	    // Get the number of processes
	    	int numprocs;
	    	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	    // Get the ID number of the process
	    	int myRank;
	    	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	//Starting time
	    high_resolution_clock::time_point tStart,tEnd;
	    if(myRank==0)
	    {
	    	tStart = high_resolution_clock::now();
	    }

	//Divide the isodetection rings between the processes.
    	std::vector<int> iList = WorkDivision(numprocs);
	 	MPI_Barrier(MPI_COMM_WORLD);
	 	if(myRank==0&&numprocs>180)
	 	{
	 		cout <<"Warning: More MPI processes than isodetection rings. " <<numprocs-180 <<" have nothing to do." <<endl;
	 	}
	
	//Create folder for histograms DaMaSCUS/results/SimID_histograms
		if(myRank==0)
		{
			std::chrono::time_point<std::chrono::system_clock> start;
			start = std::chrono::system_clock::now();
			std::time_t start_time = std::chrono::system_clock::to_time_t(start);
			cout <<"\n##############################"<<endl
			<<"DaMaSCUS"<<version<<" - Data Analysis"<<endl<<endl
			<<"Starting Time:\t" <<std::ctime(&start_time)
			<<"Simulation ID:\t" <<SimID <<endl
			<<"Experiment:\t" <<experiment <<endl<<endl;
			cout <<"Creating folder for histograms."<<endl;
			std::string sPath = "../results/"+SimID+"_histograms";
			mode_t nMode = 0733; // UNIX style permissions
			int nError = 0;
			#if defined(_WIN32)
				nError = _mkdir(sPath.c_str()); // can be used on Windows
			#else 
			  nError = mkdir(sPath.c_str(),nMode); // can be used on non-Windows
			#endif
			if (nError != 0) 
			{
			  cout <<"Folder already exists, existing data will be overwritten." <<endl<<endl;
			}
			else
			{
				cout <<"Done." <<endl<<endl;
			}
		}
	//Temporary binary output files, which will be translated to ascii at the end.
		if(myRank==0)cout <<"Creating temporary files."<<endl;
		MPI_Offset offset;
 		MPI_File   file_speed;
 		MPI_File   file_rate;
 	  	MPI_Status status;
 	  	//Speed output
		 	string SpeedFilename="../results/"+SimID+".speedTemp";
		 	int test=MPI_File_open(MPI_COMM_WORLD,SpeedFilename.c_str(),MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,MPI_INFO_NULL, &file_speed);
		    //if it already exists, it will be overwritten!
			if(test != MPI_SUCCESS)
			{
	  			if (myRank == 0) MPI_File_delete(SpeedFilename.c_str(),MPI_INFO_NULL);
	  			//MPI_Barrier(MPI_COMM_WORLD);
		  		test=MPI_File_open(MPI_COMM_WORLD,SpeedFilename.c_str(),MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,MPI_INFO_NULL, &file_speed);
		  	}
		  	//Offset for (ring,vMean,vError)
	   		offset = iList[myRank] * 3 * sizeof(double);
	  		MPI_File_seek(file_speed, offset,MPI_SEEK_SET);
	  	//Event rate output
	  		string RateFilename="../results/"+SimID+".rateTemp";
	 		test=MPI_File_open(MPI_COMM_WORLD,RateFilename.c_str(),MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,MPI_INFO_NULL, &file_rate);
		   	//if it already exists, it will be overwritten!
				if(test != MPI_SUCCESS)
				{
		 				if (myRank == 0) MPI_File_delete(RateFilename.c_str(),MPI_INFO_NULL);
		  			test=MPI_File_open(MPI_COMM_WORLD,RateFilename.c_str(),MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,MPI_INFO_NULL, &file_rate);
		  		}
		 		//Offset for (ring,rate,error,analytic rate)
	   			offset = iList[myRank] * 4 * sizeof(double);
	  			MPI_File_seek(file_rate, offset,MPI_SEEK_SET);  	
	
	//Synchronize processes:
		MPI_Barrier(MPI_COMM_WORLD);

////////////////////////////////////////////////////////////
		//DM Density	
 		//The master reads and copies the rho file from the /data directory. Afterwards he shares the vector with everyone.
			std::vector<std::vector<double>> rho(180, vector<double>(2));
			if(myRank==0)
			{
				cout <<"Reading in local DM densities."<<endl;
				ifstream f2;
				ofstream copy;
				f2.open("../data/"+SimID+".rho");
				copy.open("../results/"+SimID+".rho");
				if(f2.is_open())
				{
					for(int i =0;i<180;i++)
					{
						int ring;
						double b1,b2;
						f2 >> ring;
						f2 >> b1;
						f2 >> b2;
						copy <<ring <<"\t" <<b1 <<"\t" <<b2 <<endl;
						rho[i][0]=b1*GeV/cm/cm/cm;
						rho[i][1]=b2*GeV/cm/cm/cm;
					}
				}
				else
				{
					cout <<SimID+".rho does not exist."<<endl;
				}
				f2.close();
				copy.close();
			}
			//Broadcast the rho vector
			if(myRank==0) cout <<"Broadcast local DM densities to all MPI processes." <<endl;
				for(int k=0;k<180;k++)
				{
					MPI_Bcast(&rho[k].front(),2,MPI_DOUBLE,0,MPI_COMM_WORLD);
				}
			


	//Go through the files and analyze isodetection ring by ring.
		if(myRank==0) 
		{
			cout <<"Start data analysis..." <<endl<<endl
			<<"MPI rank\tIsodetection ring\tLocal Progress\tComputing time[s]\tResidual time estimate[s]"<<endl;
		}
			
		double R_A; //the analytic result for the event rate must only be computed once.
		double durationMean=0.0; //for the estimation of the residual time
		for(int i=iList[myRank];i<iList[myRank+1];i++)
		{
			//Beginning time of this isodetection ring
				high_resolution_clock::time_point tIRstart;
				tIRstart = high_resolution_clock::now();
			
			//Open Files
				MPI_File VelocityFile,WeightFile;
				MPI_Status status;
				string VelocityFilename="../data/"+SimID+"_data/velocity."+std::to_string(i);
				string WeightFilename="../data/"+SimID+"_data/weights."+std::to_string(i);
				int rc = MPI_File_open(MPI_COMM_SELF, VelocityFilename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &VelocityFile );
				if (rc) cout <<"Unable to read file " <<VelocityFilename <<"."<<endl;
				rc = MPI_File_open(MPI_COMM_SELF, WeightFilename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &WeightFile );
				if (rc) cout <<"Unable to read file " <<WeightFilename <<"."<<endl;

			//First we calculate the histogram domain, the weighted speed average, variance and standard error. 
				std::vector<double> AverageSpeed;
				//1. Compute the weighted average speed and find the maximum speed in the sample, which defines our histogram's domain. 
					double vMax = 0.0;
					double vMin=0.0;
					double WeightSum=0.0;
					double SpeedSum=0.0;
					for(int j=0;j<SampleSize;j++)
					{
						Eigen::Vector3d v;
						MPI_File_read(VelocityFile,v.data(),v.size(),MPI_DOUBLE,&status);
						double speed=v.norm();
						if(speed>vMax)			vMax=speed;
						//Compute weighted mean speed;
							double weight;
							MPI_File_read(WeightFile,&weight,1,MPI_DOUBLE,&status);
							WeightSum+=weight;
							SpeedSum+=speed*weight;
					}
					double MeanSpeed=SpeedSum/WeightSum;
				
				//2. Now we can compute the variance for Scott's rule, as well as the standard error (not deviation!) of the Mean Velocity with Cochran formula
					double SpeedVariance=0.0;
					double MeanWeight=WeightSum/SampleSize;
					double sum1=0.0;
					double sum2=0.0;
					double sum3=0.0;
					//Jump back to the top of file 
						MPI_File_seek(VelocityFile, 0,MPI_SEEK_SET);
						MPI_File_seek(WeightFile,0,MPI_SEEK_SET);
					for(int j=0;j<SampleSize;j++)
					{
						//Read in velocity and data weight
							Eigen::Vector3d v;
							double weight;
							MPI_File_read(VelocityFile,v.data(),v.size(),MPI_DOUBLE,&status);
							MPI_File_read(WeightFile,&weight,1,MPI_DOUBLE,&status);
							double speed=v.norm();					
						//Compute variance of the sample;
							SpeedVariance+=weight*pow(speed-MeanSpeed,2.0)/WeightSum;
						//Standard Error
							sum1+=pow(weight*speed-MeanWeight*MeanSpeed,2.0);
							sum2+=(weight-MeanWeight)*(weight*speed-MeanWeight*MeanSpeed);
							sum3+=pow(weight-MeanWeight,2.0);
					}
					double StandardError = sqrt(SampleSize/(SampleSize-1)/pow(WeightSum,2.0)*(sum1-2.0*MeanSpeed*sum2+pow(MeanSpeed,2.0)*sum3));
					//Save the average speed and standard error into a file
						AverageSpeed.push_back(i);
						AverageSpeed.push_back(MeanSpeed);
						AverageSpeed.push_back(StandardError);
						MPI_File_write(file_speed,AverageSpeed.data(),AverageSpeed.size(),MPI_DOUBLE,&status);
				//3. Use Scott's rule for the bin width
					double h = 3.5*sqrt(SpeedVariance)/pow(SampleSize,1.0/3.0);
					double bins = ceil((vMax-vMin)/h);
			
			//Create Histogram including errors
				vector< vector<double> > VelocityHistogram (bins, vector<double>(4));
				for(int j =0;j<bins;j++)
				{
					VelocityHistogram[j][0]=vMin+j*h+h/2.0;	//x coordinate / bin position
					VelocityHistogram[j][1]=0.0;			//y coordinate / bin height
					VelocityHistogram[j][2]=h/2.0;			//error on x 
					VelocityHistogram[j][3]=0.0;			//error on y
				} 
				//Jump back to the top of file 
					MPI_File_seek(VelocityFile, 0,MPI_SEEK_SET);
					MPI_File_seek(WeightFile,0,MPI_SEEK_SET);
				for(int j=0;j<SampleSize;j++)
				{
					Eigen::Vector3d v;
					double weight;
					int bin;
					MPI_File_read(VelocityFile,v.data(),v.size(),MPI_DOUBLE,&status);
					MPI_File_read(WeightFile,&weight,1,MPI_DOUBLE,&status);
					bin=(v.norm()-vMin)/h;
					VelocityHistogram[bin][1]+=weight;
					VelocityHistogram[bin][3]+=weight*weight;
					// WeightSum+=weight;
					// VelocitySum+=v.norm()*weight;
				}
				//Normalize the histogram and the average, and the errors
					for (int j=0;j<bins;j++)
					{
						VelocityHistogram[j][1]/=h*WeightSum;
						VelocityHistogram[j][3]/=h*h*WeightSum*WeightSum;
						VelocityHistogram[j][3]=sqrt(VelocityHistogram[j][3]);
					} 			

			//Close Files
				MPI_File_close( &VelocityFile );
				MPI_File_close( &WeightFile );
			//Save velocity histogram
				ofstream f;
				f.open("../results/"+SimID+"_histograms/speed."+std::to_string(i));
				for(int j=0;j<bins;j++) f<<VelocityHistogram[j][0] <<"\t" <<VelocityHistogram[j][1] <<"\t" <<VelocityHistogram[j][2]<<"\t" <<VelocityHistogram[j][3] <<endl;
				f.close();
			//Create and save eta-Histogram
				std::vector<std::vector<double>> etaH = EtaHistogram(VelocityHistogram);
				f.open("../results/"+SimID+"_histograms/eta."+std::to_string(i));
				for(unsigned int j=0;j<etaH.size();j++) f<<etaH[j][0] <<"\t"<<etaH[j][1] <<"\t"<<etaH[j][2] <<"\t"<<etaH[j][3]<<endl;
				f.close();
			//Calculate the event rate.
				std::vector<double> TotalRate;
				if(experiment=="CRESST-II")
				{
					//Create Recoil Spectrum Histograms
						std::vector<std::vector<double>> dRdEH_O = dRdEHistogram(0.22,16,rho[i],mChi,sigma0,etaH);
						std::vector<std::vector<double>> dRdEH_Ca = dRdEHistogram(0.14,40,rho[i],mChi,sigma0,etaH);
						std::vector<std::vector<double>> dRdEH_W = dRdEHistogram(0.64,184,rho[i],mChi,sigma0,etaH);
					 //Save dRdE histograms for CRESST-II	
						f.open("../results/"+SimID+"_histograms/dRdE."+std::to_string(i));
						//Emin and Emax for the histogram
							double Emin=0.3*keV;
							std::vector<double> EmaxV;
							EmaxV.push_back(dRdEH_O.back()[0]+3.0*62*eV);
							EmaxV.push_back(dRdEH_Ca.back()[0]+3.0*62*eV);
							EmaxV.push_back(dRdEH_W.back()[0]+3.0*62*eV);
							double Emax = *std::max_element( EmaxV.begin(), EmaxV.end() );	
							//For very low signal rates we only save the first low 0.1keV of the spectrum.
							if(Emax<0.4*keV) Emax=0.4*keV;
						double dE=(Emax-Emin)/100.;
						//Below 0.3keV the efficiency is assumed to zero. 0.3 keV is the threshold
						for(double E=0.3*keV;E<=Emax;E+=dE)
						{
							std::vector<double> dR = dRdE_CRESSTII_MC(E,dRdEH_O,dRdEH_Ca,dRdEH_W,vEarth.norm());
							f <<E <<"\t" <<dR[0] <<"\t" <<dE/2.0<<"\t" <<dR[1] <<"\t" <<dRdE_CRESSTII_A(E,rhoDM,mChi,sigma0,vEarth.norm()) <<endl; 
						}
						f.close();
					//Calculate Total event rate
						//Compute the analytic result at the first iteration
							if(i==iList[myRank])
							{
								R_A = R_CRESSTII_A(rhoDM,mChi,sigma0,vEarth.norm());
							}
						//Compute result for MC
							std::vector<double> R_MC = R_CRESSTII_MC(dRdEH_O,dRdEH_Ca,dRdEH_W,mChi,vEarth.norm());;
							// cout <<"MPI Process:\t" <<myRank <<"\t IsoRing:\t" <<i<<"\tR = "<<R_MC[0] <<"+-"<<R_MC[1] <<"\t("<<R_A<<")"<<endl;
					//Save rate including the analytic result for the transparent earth
						TotalRate.push_back(i);
						TotalRate.push_back(R_MC[0]);
						TotalRate.push_back(R_MC[1]);
						TotalRate.push_back(R_A);
						MPI_File_write(file_rate,TotalRate.data(),TotalRate.size(),MPI_DOUBLE,&status);

				}
				else if(experiment=="LUX")
				{
					//Create and save dRdE-Histogram
						std::vector<std::vector<double>> dRdEH = dRdEHistogram(1,131,rho[i],mChi,sigma0,etaH);
						f.open("../results/"+SimID+"_histograms/dRdE."+std::to_string(i));
						for(unsigned int j=0;j<dRdEH.size();j++) f<<dRdEH[j][0] <<"\t"<<dRdEH[j][1] <<"\t"<<dRdEH[j][2] <<"\t"<<dRdEH[j][3]<<"\t"<<dRdErA(dRdEH[j][0],1,131,rhoDM,mChi,sigma0,vEarth.norm())<<endl;
						f.close();
					//Calculate Total event rate
						//Compute the analytic result at the first iteration
							if(i==iList[myRank]) R_A = R_LUX_A(rhoDM,mChi,sigma0,vEarth.norm());
								
						//Compute the MC result
							std::vector<double> R_MC = R_LUX_MC(dRdEH);

					//Save rate including the analytic result for the transparent earth
						TotalRate.push_back(i);
						TotalRate.push_back(R_MC[0]);
						TotalRate.push_back(R_MC[1]);
						TotalRate.push_back(R_A);
						MPI_File_write(file_rate,TotalRate.data(),TotalRate.size(),MPI_DOUBLE,&status);
				}
				else if(experiment!="None")
				{
					if(myRank==0) cout <<"Error: Experiment not recognized. No detection rate will be computed." <<endl;
					experiment="None";
				}
				//Console output to give an indication of the progress
					//computing time of this isodetection ring
						high_resolution_clock::time_point tIRend;
						tIRend = high_resolution_clock::now();
						double durationIR =1e-6*duration_cast<microseconds>( tIRend - tIRstart ).count();
					//Mean computing time of all finished isodetection rings so far. This will be used to make the prognosis.
						durationMean=1.0*(i-iList[myRank])/(i-iList[myRank]+1)*durationMean+1.0*durationIR/(i-iList[myRank]+1);
					
					cout <<myRank<<"\t\t" <<i<<"\t\t\t" <<floor(100.0*(i-iList[myRank]+1)/(iList[myRank+1]-iList[myRank]))<<"\%\t\t"<<ceil(durationIR)<<"\t\t\t"<<ceil((iList[myRank+1]-i-1)*durationMean)<<endl;	
				

		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(myRank==0) cout <<"\nData analysis complete." <<endl;
		MPI_File_close(&file_speed);
		MPI_File_close(&file_rate);
		
		

	//Create ASCII Output
		if(myRank==0)
		{
			cout <<"Creating ASCII output." <<endl;
			//Binary Input files
			//Speed
				MPI_File speedfile;
				int rc = MPI_File_open(MPI_COMM_SELF, SpeedFilename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &speedfile );
				if (rc) cout <<"Unable to read file " <<SpeedFilename <<"."<<endl;
				//ASCII Output files
					ofstream f;
					f.open("../results/"+SimID+".speed");
				//Write results in ascii files.
					for(int ring=0;ring<180;ring++)
					{
						std::vector<double> speedbuffer (3);
						MPI_File_read(speedfile,speedbuffer.data(),speedbuffer.size(),MPI_DOUBLE,&status);
						f <<ring <<"\t" <<speedbuffer[1]<<"\t"<<speedbuffer[2]<<endl;
					}
				//Close file
					f.close();
					MPI_File_close(&speedfile);
			//Rate
				if(experiment!="None")
				{
					MPI_File ratefile;
					rc = MPI_File_open(MPI_COMM_SELF, RateFilename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &ratefile );
					if (rc) cout <<"Unable to read file " <<RateFilename <<"."<<endl;
					//ASCII Output files
						ofstream g;
						g.open("../results/"+SimID+"."+experiment);
					//Write results in ascii files.
						for(int ring=0;ring<180;ring++)
						{
							std::vector<double> ratebuffer (4);
							MPI_File_read(ratefile,ratebuffer.data(),ratebuffer.size(),MPI_DOUBLE,&status);
							g <<ring <<"\t" <<ratebuffer[1]<<"\t"<<ratebuffer[2]<<"\t"<<ratebuffer[3]<<endl;
						}
					//Close files
						g.close();
						MPI_File_close(&ratefile);
				}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		
	
//FINALIZE
	//Master process finishes logfile
		if(myRank==0)
		{
			//Delete temporary files
				cout <<"Delete temporary files and finish." <<endl;
				MPI_File_delete(SpeedFilename.c_str(),MPI_INFO_NULL);
				MPI_File_delete(RateFilename.c_str(),MPI_INFO_NULL);
			//Ending time and computing time
				tEnd = high_resolution_clock::now();
				double durationTotal =1e-6*duration_cast<microseconds>( tEnd - tStart ).count();
				cout <<"\nProcessing Time:\t"<< durationTotal<<"s ("<< floor(durationTotal/3600.0)<<":"<<floor(fmod(durationTotal/60.0,60.0))<<":"<<floor(fmod(durationTotal,60.0))<<":"<<floor(fmod(1000*durationTotal,1000.0))<<")."<<endl
				<<"##############################"<<endl;
			//LogFile
				LogFile_Analysis(durationTotal,numprocs);
				
		}	
	// Finalize the MPI environment.
	    MPI_Finalize();
   
	return 0;
}


