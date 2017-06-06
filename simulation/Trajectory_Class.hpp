//The simulation algorithm will return a list of "events", i.e. a time, a position vector and a velocity vector. For this we define the class "Event".
#ifndef __Trajectory_Class_hpp_
#define __Trajectory_Class_hpp_

#include <Eigen/Geometry>
#include <vector>

//Format of the vector output of eigen: {x,y,z}
	const Eigen::IOFormat VectorFormat(Eigen::StreamPrecision, Eigen::DontAlignCols,",",",","","","{","}");

class Event {
	double time;
	Eigen::Vector3d position;
	Eigen::Vector3d velocity;
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		//Constructors
		Event();
		Event(double t, Eigen::Vector3d& x, Eigen::Vector3d& v);
		Event(double t, double x, double y, double z, double vx, double vy, double vz);
		//Set Values
		void SetTime(double t);
		void SetPosition(Eigen::Vector3d& newposition);
		void SetPosition(double x,double y,double z);
		void SetVelocity(Eigen::Vector3d& newvelocity);
		void SetVelocity(double vx,double vy,double vz);
		//Return Values
		double Time();
		Eigen::Vector3d Position();
		Eigen::Vector3d Velocity();
		Event kmsec();
		//Vector Norms
		double NormPosition();
		double NormVelocity();
		int IsodetectionRing();
		double DataWeight();
		//Overloading the output operator <<
		friend std::ostream& operator<<(std::ostream &output,const Event &event);
};


class Trajectory {
	std::vector<Event> events;
	public:
		//Constructors:
			Trajectory();
			Trajectory(std::vector<Event> input);
		//Extract Info from a single Track
			int Length();
			double tStart();
			double tEarthEntry();
			double tEnd();
			Event TrackInterpolation(double t);
			int Trajectory_Type();
			int NoOfScatterings();
			double TrajectoryLength();
			std::vector<Event> DepthCrossing(double d);
			// std::vector<Detection> DetectionCrossing(double d,double L);

		//Overloading the output operator <<
			friend std::ostream& operator<<(std::ostream &output,const Trajectory &trajectory);
		//Overloading the brackets[]
			Event& operator[](unsigned int i){
				return events[i];
			}

};


#endif