/*
 *  WormAgent.h
 *
 *  Created by EJI on 17/08/2008.
 *  Copyright 2008 Centre for Computational Neuroscience and Robotics. All rights reserved.
 *
 */

#pragma once
#include "CTRNN.h"


//#define GRAD_CONE
//#define GRAD_INVCONE
//#define GRAD_ISO
#define GRAD_GAUS				

// Global constants
const double	StepSize		=	0.01;			// Fastest time-constant is now 0.1 (10ms)
const double	Pi				=	3.1415926;
const double	MaxDist			=	4.5; // Half the radius of the big petri dish (in cm)
const double	MaxVel			=	0.022;			// Forward velocity (in cm/s) CHECK WHAT THE REAL VELOCITY SHOULD BE
const double	ChemDiffConst	=	2*pow(1.61,2);	// Simulated chemical environment according to Ward, 1973 as described in Ferree and Lockery 1999 equation 14.
const double	HST				=	4.2;			// Head Sweep Time, T=4.2sec, According to Ferree, Marcotte, Lockery, 1997.
const double	HSP				=	(2*Pi)/HST;		//Head-sweep period 2*Pi/T, According to Ferree, Marcotte, Lockery, 1997.
const int		VelDelta		=	(int) (HST/StepSize);
const double	SensorNoiseStd	=	0.0; //0.0000001;		// Standard deviation of the Gaussian noise to the sensor
const double	MotorNoiseStd	=	0.05;			// Standard deviation of the Gaussian noise to the sensor
//const double	NeckLimit		=	(Pi/4);			// Lv is 2Ld. Lo is 1. So, ArcTan[2-1/1]

// The WormAgent class declaration
class WormAgent {
public:
	// The constructor
	WormAgent(int newsize = 0);	// Construct from evolutionary algorithm
	WormAgent(char* fnm);	// Construct from file
	// The destructor
	~WormAgent();
	
	// Accessors
	int CircuitSize(void) {return size;};
	void SetCircuitSize(int news) {size = news;};
	void InitialiseCircuit(int newcs);
	void InitialiseSensors();
	void InitialiseMotors();
	void SetWormParametersFromFile(char* fnm);
	double PositionX(void) {return px;};
	void SetPositionX(double newx) {px = newx;};
	double PositionY(void) {return py;};
	void SetPositionY(double newy) {py = newy;};
	double VelocityX(void) {return vx;};
	void SetVelocityX(double newx) {vx = newx;};
	double VelocityY(void) {return vy;};
	void SetVelocityY(double newy) {vy = newy;};
	double Orientation(void) {return orient;};
	void SetOrientation(double newo) {orient = newo;};
	double OffsetCPG(void) {return CPGoffset;};
	void SetOffsetCPG(double newo) {CPGoffset = newo;};
	double ChemConDelayed(void) {return chemConDelayed;};
	double ChemCon(void) {return chemCon;};
	void SetChemCon(double newc) {chemCon = newc;};
	double SensorN(void) {return sensorN;};
	void SetSensorN(double newn) {sensorN = newn;};
	double SensorM(void) {return sensorM;};
	void SetSensorM(double newn) {sensorM = newn;};
	double SensorD(void) {return sensorD;};
	void SetSensorD(double newn) {sensorD = newn;};
	double OutputGain(void) {return outputGain;};
	void SetOutputGain(double newdc) {outputGain = newdc;};
	int EnvironmentType(void) {return envtype;};
	void SetEnvironmentType(int i) {envtype = i;};
	double DistanceToCentre(void) {return distanceToCentre;};
	double OnWeight(int i) {return onWeight[i];};
	void SetOnWeight(int i, double news) {onWeight[i] = news;};
	double OffWeight(int i) {return offWeight[i];};
	void SetOffWeight(int i, double news) {offWeight[i] = news;};
	double CpgWeight(int i) {return cpgWeight[i];};
	void SetCpgWeight(int i, double news) {cpgWeight[i] = news;};
	double OnCell() {return V_on;};
	void SetOnCell(double newc) {V_on = newc;};
	double OffCell() {return V_off;};
	void SetOffCell(double newc) {V_off = newc;};
	double GradSteepness() {return gradSteepness;};
	void SetGradientSteepness(double newc) {gradSteepness = newc;};
	double Theta() {return theta;};
	
	// Control
	void ResetChemCon( RandomState &rs );
	void UpdateChemCon( RandomState &rs );
	void UpdateChemConHistory();
	void UpdateSensors();
    void ResetAgentsBody(RandomState &rs);
    void ResetAgentIntState(RandomState &rs);
	void ResetEnvType(RandomState &rs);
    void ResetChemical(RandomState &rs);
	void SetParameters(TVector<double> &v);
	void InitialiseAgent(double stepsize);
	void InitialiseSensorHistory();	
	void PrintDetail(ofstream &file);
	void PrintPath(ofstream &file);
	void Step(double StepSize, RandomState &rs, double timestep);
	
	double avgvel;	
	double pastTheta,DthetaDt;
	int timer;
	double NMdiff;
	double concave, convex, pathlength;
	double pushConcave, pushConvex;
	double gradSteepness;
	TVector<double> histConcave,histConvex;
	double px, py, vx, vy, orient, distanceToCentre, theta ;
	double CPGoffset, chemCon, chemConDelayed, V_on, V_off, outputGain, whiteNoise;
	int envtype,size,dorsalMotorNeuron,ventralMotorNeuron,delta,tau;
	TVector<double> onWeight,offWeight,cpgWeight;
	double chemConHistory[50300];	// RunDuration / StepSize
	double currentConc, pastConc;
	double sensorD, sensorN, sensorM;
	double dSensorD, dSensorN, dSensorM;
	int iSensorD, iSensorN, iSensorM;
	double tempDiff;
	CTRNN NervousSystem;
};

