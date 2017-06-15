#include "Trajectory_Class.hpp"
#include <iostream>
#include "Physical_Parameters.hpp"
#include "General_Utilities.hpp"

//Event Class
	//Constructors
		Event::Event(){
			time=0.0;
			position=Eigen::Vector3d (0.0,0.0,0.0);
			velocity=Eigen::Vector3d (0.0,0.0,0.0);
		}
		Event::Event(double t, Eigen::Vector3d& x, Eigen::Vector3d& v){
			time=t;
			position=x;
			velocity=v;
		}
		Event::Event(double t, double x, double y, double z, double vx, double vy, double vz){
			time=t;
			position=Eigen::Vector3d (x,y,z);
			velocity=Eigen::Vector3d (vx,vy,vz);
		}
	//Set Values
		void Event::SetTime(double t){
			time=t;
		}
		void Event::SetPosition(Eigen::Vector3d& newposition){
			position=newposition;
		}
		void Event::SetVelocity(Eigen::Vector3d& newvelocity){
			velocity=newvelocity;
		}
		void Event::SetPosition(double x,double y,double z){
			position=Eigen::Vector3d (x,y,z);
		}
		void Event::SetVelocity(double vx,double vy,double vz){
			velocity=Eigen::Vector3d (vx,vy,vz);
		}
	//Return Values
		double Event::Time()
		{
			return time;
		}
		Eigen::Vector3d Event::Position(){
			return position;
		}
		Eigen::Vector3d Event::Velocity()
		{
			return velocity;
		}
		//Returns the event in sec, km, km/sec.
		Event Event::kmsec()
		{
			double t=time/sec;
			Eigen::Vector3d x=1/km*position;
			Eigen::Vector3d v=sec/km*velocity;
			return Event(t,x,v);
		}
	//Vector Norms
		double Event::NormPosition(){
			return position.norm();
		}
		double Event::NormVelocity(){
			return velocity.norm();
		}
	//Returns the number of the iso ring 1-180;
		int Event::IsodetectionRing()
		{
			return acos(position.normalized().dot(vEarth.normalized()))/M_PI*180.0;
		}
	//Returns the weight of a velocity data point. Should only be used for events returned by Trajectory::DepthCrossing
		double Event::DataWeight()
		{
			//1/cos gamma, where gamma is the angle between the velocity and the normal vector of the isoring, where the particle crossed the isoring
			return 1.0/abs(position.normalized().dot(velocity.normalized()));
		}
	//Overload <<
		std::ostream& operator<<(std::ostream &output,const Event &event){
			return output 	<<"{"
							<<event.time
							<<","
							<<event.position.format(VectorFormat)
							<<","
							<<event.velocity.format(VectorFormat)
							<<"}";
		}



//Trajectory Class
	//Constructors:
		Trajectory::Trajectory(){

		}
		Trajectory::Trajectory(std::vector<Event> input){
			events=input;
		}
	//Extract Info from a single Track
		int Trajectory::Length()
		{
			return events.size();
		}

		double Trajectory::tStart(){
			return events.front().Time();
		}
		double Trajectory::tEnd(){
			return events.back().Time();
		}
		double Trajectory::tEarthEntry(){
			return events[1].Time();
		}
		Event Trajectory::TrackInterpolation(double t){
			if(t<tStart()||t>tEnd())
			{
				std::cout <<"Error: TrackInterpolation was given an argument outside the Trajectorys range."<<std::endl;
				printf("Range is [%f,%f], but t = %f\n",tStart()/sec,tEnd()/sec,t);
				return Event();
			}
			else
			{
				for(unsigned int i=1;i<events.size();i++)
				{
					if(t<=events[i].Time())
					{
						Eigen::Vector3d x=events[i-1].Position()+(t-events[i-1].Time())*events[i-1].Velocity();
						Eigen::Vector3d v=events[i-1].Velocity();
						return Event(t,x,v);
					}
				}
			}
			cout <<"Error in TrackInterpolation(): This shouldn't happen..."<<endl;
			return Event();	
		}
		Event Trajectory::TrajectoryStart()
		{
			return events.front();
		}
		Event Trajectory::TrajectoryEnd()
		{
			return events.back();
		}

		//We differentiate three types of Trajectory.
		// 1 - particle enters the earth and leaves
		// 2 - particle enters the earth and does not leave (typically, the velocity cutoff is reached)
		// 3 - particle misses the earth, this should not be passible since the IC generator makes sure that the earth is hit. But with manually set IC, it becomes possible. 
		int Trajectory::Trajectory_Type()
		{
			int Wcase;
			if(events.back().Position().norm()<rEarth) Wcase=2;
			else if(events.size()>2) Wcase =1;
			else Wcase=3;
			return Wcase;
		}
		int Trajectory::NoOfScatterings()
		{
			int nScatterings;
			int Wcase =  Trajectory_Type();
			if (Wcase==1) nScatterings=events.size()-4;
			else if (Wcase==2) nScatterings=events.size()-2;
			else nScatterings=0;
			return nScatterings;
		}
		double Trajectory::TrajectoryLength()
		{
			double length=0.0;
			for(unsigned int i=1;i<events.size();i++)
			{
				length+= (events[i].Position()-events[i-1].Position()).norm();
			}
			return length;
		}

		std::vector<Event> Trajectory::DepthCrossing(double d)
		{
			vector<Event> EventList;
			double a,b,c;
			for(int i=1;i<Length()-2;i++)
			{
				Eigen::Vector3d r1=events[i].Position();
				Eigen::Vector3d r2=events[i+1].Position();
				//Both points inside the detection sphere:
				if(r1.norm()<rEarth-d&&r2.norm()<rEarth-d) continue;
				//Both points outside the detection sphere
				else if(r1.norm()>rEarth-d&&r2.norm()>rEarth-d)
				{
					a=pow(r2.norm(),2.0)+pow(r1.norm(),2.0)-2.0*r1.dot(r2);
					b=2.0*(r1.dot(r2)-pow(r1.norm(),2.0));
					c=pow(r1.norm(),2.0)-pow(rEarth-d,2.0);
					//Real solutions
						//Two solutions, but we have to check if they are between 0 and 1
							if(b*b-4*a*c>0)
							{
								double l1= (-b-sqrt(b*b-4*a*c))/(2*a);
								double l2= (-b+sqrt(b*b-4*a*c))/(2*a);
								if(l1>=0.0&&l1<=1.0) 
								{
									double t1=events[i].Time();
									double t2=events[i+1].Time();
									EventList.push_back(TrackInterpolation(t1+l1*(t2-t1)));
								}
								if(l2>=0.0&&l2<=1.0) 
								{
									double t1=events[i].Time();
									double t2=events[i+1].Time();
									EventList.push_back(TrackInterpolation(t1+l2*(t2-t1)));
								}

							}
						//One Solution (unlikely...)
							else if(b*b-4*a*c==0.0)//will never happen...
							{
								double l=-b/(2*a);
								if(l>=0.0&&l<=1.0) 
								{
									double t1=events[i].Time();
									double t2=events[i+1].Time();
									EventList.push_back(TrackInterpolation(t1+l*(t2-t1)));
								}
							}
				}
				//One point inside the detection sphere, one outside
				else
				{
					//cout <<"|r1|=" << r1.norm()/km <<" |r2|=" << r2.norm()/km <<endl;
					a=pow(r1.norm(),2)+pow(r2.norm(),2)-2*r1.dot(r2);
					b=2*(r1.dot(r2)-pow(r1.norm(),2));
					c=pow(r1.norm(),2)-pow(rEarth-d,2);
					double l;
					if(r1.norm()<rEarth-d)
					{
						l=(-b+sqrt(b*b-4*a*c))/(2*a);
						//cout <<"Outgoing!!! l= " <<l<<endl;

					} 
					else 
					{
						l=(-b-sqrt(b*b-4*a*c))/(2*a);
						//cout <<"Incoming!!! l=" <<l<<endl;
					}
				
					if(l>=0.0&&l<=1.0) 
					{

						double t1=events[i].Time();
						double t2=events[i+1].Time();
						EventList.push_back(TrackInterpolation(t1+l*(t2-t1)));
					}
				}
			}
			return EventList;
		}

		


	//Overload <<
		std::ostream& operator<<(std::ostream &output,const Trajectory &trajectory){
			output <<"{";
			for(unsigned int i=0;i<trajectory.events.size();i++)
			{
				output << trajectory.events[i];
				if(i<(trajectory.events.size()-1))
				{
					output <<",";
				}
			}
			output <<"}";
			return output;
		}
		






