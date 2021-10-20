// ***************************************
// Worm Chemotaxis
// ***************************************

#include "CTRNN.h"
#include "WormAgent.h"
#include "TSearch.h"
#include <sstream>

//#define EVOLVE
//#define PRINTTOFILE

//#define PSON
//#define PSOFF
//#define PSBOTH

//#define ONKILL
//#define OFFKILL

// Global constants
const int		CircuitSize				= 2;
const double	RunDuration				= 500.0;
const int		MaxRepetitions			= 5;
const int		MaxReplacements			= 10;
const double	WeightRange				= 15.0;
const double	SensoryWeightRange		= 1500.0;
const double	BiasRange				= 15.0;
const double	CloseEnoughRadius		= 0.5;
const double	MinSensorN				= 10*StepSize;
const double	MaxSensorN				= HST;
const double	MinSensorM				= 10*StepSize;
const double	MaxSensorM				= HST;
const double	MinSensorD				= 10*StepSize;
const double	MaxSensorD				= HST;
const double	MinNeckTurnGain			= 1.0;
const double	MaxNeckTurnGain			= 3.0;
//const double	PirouetteProb			= 0.000333333; //(2/60)*0.01; // They reverse about twice a minute (according to Shawn and Serge). In 500 secs, that's 16.66 on avg.
const double	PirouetteProb			= 0.00004; // Smaller. So that it doesn't evolve to depend on it. On avg, 2 per trial (500secs) [(2/500) * 0.01]
const double	SizeOfStep				= 0.002;	// Size of steps for the riverchip experiment

const double	RiverChipSteepness		= 50.0;

// Global variables
int	VectSize = 9;

// ------------------------------------
// Genotype-Phenotype Mapping Functions
// ------------------------------------
void GenPhenMapping(TVector<double> &gen, TVector<double> &phen)
{
	phen(1) = MapSearchParameter( gen(1), -SensoryWeightRange, SensoryWeightRange); // w_on
	phen(2) = MapSearchParameter( gen(2), -SensoryWeightRange, SensoryWeightRange);	// w_off
	phen(3) = MapSearchParameter( gen(3), 0.0, WeightRange);						// w_osc
	phen(4) = MapSearchParameter( gen(4), -WeightRange, WeightRange);				// w_osc
	phen(5) = MapSearchParameter( gen(5), -BiasRange, BiasRange);					// bias or threshold
	phen(6) = MapSearchParameter( gen(6), MinSensorN, MaxSensorN);					// Sensory cell integration time, current
	phen(7) = MapSearchParameter( gen(7),  MinSensorM, MaxSensorM);					// Sensory cell integration time, past
	phen(8) = MapSearchParameter( gen(8),  MinSensorD, MaxSensorD);					// Sensory cell time constant delay
	phen(9) = MapSearchParameter( gen(9), MinNeckTurnGain, MaxNeckTurnGain);		// Neck turning gain
}

// ------------------------------------
// Fitness function
// ------------------------------------
double EvaluationFunction(TVector<double> &v, RandomState &rs)
{
	double t,fitness=0.0;
	double fA,accdist,totaldist;
	int k,condition,repetitions,replacements;
	TVector<double> phenotype;

	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(v, phenotype);
	WormAgent Worm(CircuitSize);
	Worm.SetParameters(phenotype);
	Worm.InitialiseAgent(StepSize);

	k=0;
	for (condition = 4; condition <= 4; condition++){
		for (repetitions = 1; repetitions <= MaxRepetitions; repetitions++)
		{
			Worm.ResetAgentIntState(rs);
			for (replacements = 1; replacements <= MaxReplacements; replacements++)
			{
				Worm.ResetAgentsBody(rs);
				Worm.ResetChemCon(rs);
				Worm.UpdateChemCon(rs);
				Worm.InitialiseSensorHistory();
				accdist = 0.0;
				for (t = StepSize; t <= RunDuration; t += StepSize)
				{
					Worm.Step(StepSize,rs,t);
					Worm.UpdateChemCon(rs);
					Worm.UpdateChemConHistory();
					Worm.UpdateSensors();

					if (rs.UniformRandom(0.0, 1.0) < PirouetteProb)
						Worm.SetOrientation(rs.UniformRandom(0, 2*Pi));

					if (condition == 1)
						Worm.SetOnCell(0.0);
					if (condition == 2)
						Worm.SetOffCell(0.0);

					accdist += Worm.DistanceToCentre();
				}

				totaldist = (accdist/(RunDuration/StepSize));
				fA = (MaxDist - totaldist)/MaxDist;
				fA = fA < 0 ? 0.0 : fA;
				fitness += fA;
			}
		}
		k++;
	}
	return fitness/(k*MaxRepetitions*MaxReplacements);
}

// ------------------------------------
// Behavioral Analysis
// ------------------------------------
void ExampleRun()
{
	ofstream ExampleFile,PirFile;
	double t;
	long IDUM=time(0);
	int repetitions,replacements;
	RandomState rs;
	rs.SetRandomSeed(IDUM);
	ExampleFile.open("/Users/carlos/Downloads/Klino/examples.dat");
	PirFile.open("/Users/carlos/Downloads/Klino/pirxy.dat");
	WormAgent Worm("/Users/carlos/Downloads/Klino/best.ns.dat");
	Worm.InitialiseAgent(StepSize);

	for (repetitions = 1; repetitions <= 1; repetitions++)
	{
		Worm.ResetAgentIntState(rs);
		for (replacements = 1; replacements <= 1; replacements++)
		{
			Worm.ResetAgentsBody(rs);
			Worm.ResetChemCon(rs);
			Worm.UpdateChemCon(rs);
			Worm.InitialiseSensorHistory();
			Worm.SetOffsetCPG(Pi/2);
			Worm.SetOrientation(-Pi/2);
			Worm.PrintDetail(ExampleFile);
			//Worm.PrintPath(ExampleFile);
			for (t = StepSize; t <= RunDuration; t += StepSize)
			{
				Worm.Step(StepSize, rs,t);
				Worm.UpdateChemCon(rs);
				Worm.UpdateChemConHistory();
				Worm.UpdateSensors();
				if (rs.UniformRandom(0.0,1.0) < PirouetteProb){
					PirFile << Worm.PositionX() << " " << Worm.PositionY() << endl;
					Worm.SetOrientation(rs.UniformRandom(0,2*Pi));
				}

#ifdef ONKILL
				Worm.SetOnCell(0.0);
#endif
#ifdef OFFKILL
				Worm.SetOffCell(0.0);
#endif
				Worm.PrintDetail(ExampleFile);
				//Worm.PrintPath(ExampleFile);
			}
		}
	}

	ExampleFile.close();
	PirFile.close();
}

void BehavioralOverview()
{
	ofstream ExampleFile;
	double t;
	long IDUM=time(0);
	int condition,repetitions;
	RandomState rs;
	rs.SetRandomSeed(IDUM);
	WormAgent Worm("best.ns.dat");
	Worm.InitialiseAgent(StepSize);
	std::stringstream out;
	std::string s;
	char str[80];
	for (repetitions = 1; repetitions <= 1; repetitions++)
	{
		out.clear();
		out << repetitions;
		s = out.str();
		for (condition = 1; condition <= 4; condition++){
			switch (condition){
				case 2:
#ifdef GRAD_CONE
					strcpy(str,"r.onkill.nor.");
#endif
#ifdef GRAD_INVCONE
					strcpy(str,"r.onkill.inv.");
#endif
#ifdef GRAD_GAUS
					strcpy(str,"r.onkill.gau.");
#endif
#ifdef GRAD_ISO
					strcpy(str,"r.onkill.iso.");
#endif
					strcat(str,s.c_str());
					strcat(str,".dat");
					cout << str << endl;
					ExampleFile.open(str);
					break;
				case 3:
#ifdef GRAD_CONE
					strcpy(str,"r.offkill.nor.");
#endif
#ifdef GRAD_INVCONE
					strcpy(str,"r.offkill.inv.");
#endif
#ifdef GRAD_GAUS
					strcpy(str,"r.offkill.gau.");
#endif
#ifdef GRAD_ISO
					strcpy(str,"r.offkill.iso.");
#endif
					strcat(str,s.c_str());
					strcat(str,".dat");
					cout << str << endl;
					ExampleFile.open(str);
					break;
				case 4:
#ifdef GRAD_CONE
					strcpy(str,"r.bothkill.nor.");
#endif
#ifdef GRAD_INVCONE
					strcpy(str,"r.bothkill.inv.");
#endif
#ifdef GRAD_GAUS
					strcpy(str,"r.bothkill.gau.");
#endif
#ifdef GRAD_ISO
					strcpy(str,"r.bothkill.iso.");
#endif
					strcat(str,s.c_str());
					strcat(str,".dat");
					cout << str << endl;
					ExampleFile.open(str);
					break;
				default:
#ifdef GRAD_CONE
					strcpy(str,"r.normal.nor.");
#endif
#ifdef GRAD_INVCONE
					strcpy(str,"r.normal.inv.");
#endif
#ifdef GRAD_GAUS
					strcpy(str,"r.normal.gau.");
#endif
#ifdef GRAD_ISO
					strcpy(str,"r.normal.iso.");
#endif
					strcat(str,s.c_str());
					strcat(str,".dat");
					cout << str << endl;
					ExampleFile.open(str);
					break;
			}
			Worm.ResetAgentIntState(rs);
			for (double orient = 0.0;  orient < 2*Pi; orient += Pi/6)
				//for (replacements = 1; replacements <= 1; replacements++)
			{
				Worm.ResetAgentsBody(rs);
				Worm.SetOrientation(orient);
				Worm.ResetChemCon(rs);
				Worm.UpdateChemCon(rs);
				Worm.InitialiseSensorHistory();
				Worm.PrintPath(ExampleFile);
				for (t = StepSize; t <= 1.5*RunDuration; t += StepSize)
				{
					Worm.Step(StepSize, rs,t);
					Worm.UpdateChemCon(rs);
					Worm.UpdateChemConHistory();
					Worm.UpdateSensors();
					switch (condition){
						case 2:
							Worm.SetOnCell(0.0);
							break;
						case 3:
							Worm.SetOffCell(0.0);
							break;
						case 4:
							Worm.SetOnCell(0.0);
							Worm.SetOffCell(0.0);
							break;
						default:
							break;
					}
					Worm.PrintPath(ExampleFile);
				}
			}
			ExampleFile.close();
		}
	}
}

void BehavioralAnalysis()
{
	ofstream ExampleFile;
	int condition,reached;
	double timeAfterReached;
	double accdist,dist,timespeaked,totalrepeats;
	double t,timetopeak,avg;
	long IDUM=time(0);
	RandomState rs;
	rs.SetRandomSeed(IDUM);
#ifdef GRAD_CONE
        ExampleFile.open("ba.con.dat");
#endif
#ifdef GRAD_INVCONE
        ExampleFile.open("ba.inv.dat");
#endif
#ifdef GRAD_GAUS
        ExampleFile.open("ba.gau.dat");
#endif
#ifdef GRAD_ISO
        ExampleFile.open("ba.iso.dat");
#endif
	WormAgent Worm("best.ns.dat");
	Worm.InitialiseAgent(StepSize);
	double 	TempOnWeight = Worm.OnWeight(0);
	double TempOffWeight = Worm.OffWeight(0);

	//cout << "Overall behavioral analysis:" << endl;

	for (condition = 1; condition <= 4; condition++){
		//cout << "	" << condition << endl;
		avg = 0.0;
		dist = 0.0;
		timespeaked = 0.0;
		totalrepeats = 0.0;

		if (condition == 7){
			Worm.SetOffWeight(0, TempOnWeight);
			Worm.SetOffWeight(1, TempOnWeight);
		}
		if (condition == 8){
			// Reset the OFF weights to their original value
			Worm.SetOffWeight(0, TempOffWeight);
			Worm.SetOffWeight(1, TempOffWeight);
			// Set the ON weights to the off
			Worm.SetOnWeight(0, TempOffWeight);
			Worm.SetOnWeight(1, TempOffWeight);
		}

		for (int repetitions = 1; repetitions <= 1000; repetitions++)
		{
			Worm.ResetAgentIntState(rs);
			Worm.ResetAgentsBody(rs);
			Worm.ResetChemCon(rs);
			Worm.UpdateChemCon(rs);
			Worm.InitialiseSensorHistory();

			reached = 0;
			timeAfterReached = 0.0;
			accdist = 0.0;
			for (t = StepSize; t <= 2*RunDuration; t += StepSize)
			{
				Worm.Step(StepSize, rs,t);
				Worm.UpdateChemCon(rs);
				Worm.UpdateChemConHistory();
				Worm.UpdateSensors();
				switch (condition){
					case 2:
						Worm.SetOnCell(0.0);
						break;
					case 3:
						Worm.SetOffCell(0.0);
						break;
					case 4:
						Worm.SetOnCell(0.0);
						Worm.SetOffCell(0.0);
						break;
					case 5:
						Worm.SetOffCell(Worm.OnCell());
						break;
					case 6:
						Worm.SetOnCell(Worm.OffCell());
						break;
					case 7:
						Worm.SetOffCell(Worm.OnCell());
						break;
					case 8:
						Worm.SetOnCell(Worm.OffCell());
						break;
					default:
						break;
				}
				if ((reached == 0) and (Worm.DistanceToCentre() <= CloseEnoughRadius)){
					reached = 1;
					timetopeak = t;
				}
				if (reached == 1){
					accdist += Worm.DistanceToCentre();
					timeAfterReached += 1;
				}
			}
			if (reached == 1){
				avg += timetopeak;
				timespeaked++;
				dist += accdist/timeAfterReached;
			}
			totalrepeats++;
		}
		ExampleFile << avg/timespeaked << " " << timespeaked/totalrepeats << " " << dist/timespeaked << endl;
	}
	ExampleFile.close();
}

void AccConcentration()
{
	ofstream ExampleFile;
	long IDUM=time(0);
	int condition;
	double t;
	double meanAccConc, accconc, reps;
	RandomState rs;
	double maxrep = 1000;
	rs.SetRandomSeed(IDUM);
#ifdef GRAD_CONE
	ExampleFile.open("accconc.con.dat");
#endif
#ifdef GRAD_INVCONE
	ExampleFile.open("accconc.inv.dat");
#endif
#ifdef GRAD_GAUS
	ExampleFile.open("accconc.gau.dat");
#endif
#ifdef GRAD_ISO
	ExampleFile.open("accconc.iso.dat");
#endif

	WormAgent Worm("best.ns.dat");
	Worm.InitialiseAgent(StepSize);

	for (condition = 1; condition <= 4; condition++){
		meanAccConc = 0.0;
		for (int repetitions = 1; repetitions <= maxrep; repetitions++)
		{
			Worm.ResetAgentIntState(rs);
			Worm.ResetAgentsBody(rs);
			Worm.ResetChemCon(rs);
			Worm.UpdateChemCon(rs);
			Worm.InitialiseSensorHistory();

			accconc = 0.0;
			reps = 0.0;
			for (t = StepSize; t <= 2*RunDuration; t += StepSize)
			{
				Worm.Step(StepSize, rs,t);
				Worm.UpdateChemCon(rs);
				Worm.UpdateChemConHistory();
				Worm.UpdateSensors();
				switch (condition){
					case 2:
						Worm.SetOnCell(0.0);
						break;
					case 3:
						Worm.SetOffCell(0.0);
						break;
					case 4:
						Worm.SetOnCell(0.0);
						Worm.SetOffCell(0.0);
						break;
					default:
						break;
				}
				accconc += Worm.ChemCon();
				reps++;
			}
			meanAccConc += accconc/reps;
		}
		ExampleFile << meanAccConc/maxrep << endl;
	}
	ExampleFile.close();
}

void OrientationAnalysis()
{
	ofstream ExampleFile;
	int reached,condition;
	double timespeaked,totalrepeats;
	double t,orient,timetopeak,avg,dist,timeAfterReached,accdist;
	long IDUM=time(0);
	RandomState rs;
	rs.SetRandomSeed(IDUM);
	WormAgent Worm("best.ns.dat");
	Worm.InitialiseAgent(StepSize);
	double 	TempOnWeight = Worm.OnWeight(0);
	double TempOffWeight = Worm.OffWeight(0);

	//cout << "Orientation analysis:" << endl;

	for (condition = 1; condition <= 4; condition++){
		//cout << "	" << condition << endl;
		switch (condition){
			case 2:
				ExampleFile.open("orient.perf.onkill.dat");
				break;
			case 3:
				ExampleFile.open("orient.perf.offkill.dat");
				break;
			case 4:
				ExampleFile.open("orient.perf.bothkill.dat");
				break;
			case 5:
				ExampleFile.open("orient.perf.2ons.dat");
				break;
			case 6:
				ExampleFile.open("orient.perf.2offs.dat");
				break;
			case 7:
				Worm.SetOffWeight(0, TempOnWeight);
				Worm.SetOffWeight(1, TempOnWeight);
				ExampleFile.open("orient.perf.2ons.rewired.dat");
				break;
			case 8:
				// Reset the OFF weights to their original value
				Worm.SetOffWeight(0, TempOffWeight);
				Worm.SetOffWeight(1, TempOffWeight);
				// Set the ON weights to the off
				Worm.SetOnWeight(0, TempOffWeight);
				Worm.SetOnWeight(1, TempOffWeight);
				ExampleFile.open("orient.perf.2offs.rewired.dat");
				break;
			default:
				ExampleFile.open("orient.perf.normal.dat");
				break;
		}

		for (orient = 0.0;  orient <= 2*Pi; orient += Pi/90)
		{
			avg = 0.0;
			dist = 0.0;
			timespeaked = 0.0;
			totalrepeats = 0.0;
			for (int repetitions = 1; repetitions <= 10; repetitions++)
			{
				Worm.ResetAgentIntState(rs);
				Worm.ResetAgentsBody(rs);
				Worm.SetOrientation(orient);
				Worm.ResetChemCon(rs);
				Worm.UpdateChemCon(rs);
				Worm.InitialiseSensorHistory();
				reached = 0;
				timeAfterReached = 0.0;
				accdist = 0.0;
				for (t = StepSize; t <= RunDuration; t += StepSize)
				{
					Worm.Step(StepSize, rs,t);
					Worm.UpdateChemCon(rs);
					Worm.UpdateChemConHistory();
					Worm.UpdateSensors();

					if (rs.UniformRandom(0.0, 1.0) < PirouetteProb)
						Worm.SetOrientation(rs.UniformRandom(0, 2*Pi));

					switch (condition){
						case 2:
							Worm.SetOnCell(0.0);
							break;
						case 3:
							Worm.SetOffCell(0.0);
							break;
						case 4:
							Worm.SetOnCell(0.0);
							Worm.SetOffCell(0.0);
							break;
						case 5:
							Worm.SetOffCell(Worm.OnCell());
							break;
						case 6:
							Worm.SetOnCell(Worm.OffCell());
							break;
						case 7:
							Worm.SetOffCell(Worm.OnCell());
							break;
						case 8:
							Worm.SetOnCell(Worm.OffCell());
							break;
						default:
							break;
					}

					if ((reached == 0) and (Worm.DistanceToCentre() <= CloseEnoughRadius)){
						reached = 1;
						timetopeak = t;
					}
					if (reached == 1){
						accdist += Worm.DistanceToCentre();
						timeAfterReached += 1;
					}
				}
				if (reached == 1){
					avg += timetopeak;
					timespeaked++;
					dist += accdist/timeAfterReached;
				}
				totalrepeats++;
			}
			ExampleFile << orient << " " << avg/timespeaked << " " << timespeaked/totalrepeats << " " << dist/timespeaked << endl;
		}
		ExampleFile.close();
	}
}

void BearingCurvature()
{
	double peakX = 0.0, peakY = 0.0;
	int maxpoints = (int) (RunDuration/StepSize + 10);	//
	int k,arrived;
	int InflexionPoints[maxpoints];
	double aX,aY,bX,bY,cX,cY,transOrient,peakOrient,tempTheta,BearingDegree,ThetaDegree;
	double Xpoints[maxpoints],Ypoints[maxpoints],ThetaPoints[maxpoints],Bearing[maxpoints];
	ofstream ExampleFile;
	double tempBearing;
	long IDUM=time(0);
	RandomState rs;
	rs.SetRandomSeed(IDUM);
	WormAgent Worm("best.ns.dat");
	Worm.InitialiseAgent(StepSize);
	int bearingIndex;
	int BinSize = 10;
	int NumBins = (int) ((360/BinSize)+1);
	//int MaxPointsPerBin = 20000;
	double AccBearing[NumBins],BearingBins[NumBins]; //BearingBins[NumBins][MaxPointsPerBin];
	int BearingBinCounter[NumBins];
	double 	TempOnWeight = Worm.OnWeight(0);
	double TempOffWeight = Worm.OffWeight(0);
	double avgBearing; //,stddev;
	//cout << "Klinotaxis curvature versus bearing analysis:" << endl;

	for (int condition = 1; condition <= 4; condition++){
		//cout << "	" << condition << endl;
		switch (condition){
			case 2:
				ExampleFile.open("bearing.curvature.onkill.dat");
				break;
			case 3:
				ExampleFile.open("bearing.curvature.offkill.dat");
				break;
			case 4:
				ExampleFile.open("bearing.curvature.bothkill.dat");
				break;
			case 5:
				ExampleFile.open("bearing.curvature.2ons.dat");
				break;
			case 6:
				ExampleFile.open("bearing.curvature.2offs.dat");
				break;
			case 7:
				Worm.SetOffWeight(0, TempOnWeight);
				Worm.SetOffWeight(1, TempOnWeight);
				ExampleFile.open("bearing.curvature.2ons.rewired.dat");
				break;
			case 8:
				// Reset the OFF weights to their original value
				Worm.SetOffWeight(0, TempOffWeight);
				Worm.SetOffWeight(1, TempOffWeight);
				// Set the ON weights to the off
				Worm.SetOnWeight(0, TempOffWeight);
				Worm.SetOnWeight(1, TempOffWeight);
				ExampleFile.open("bearing.curvature.2offs.rewired.dat");
				break;
			default:
				ExampleFile.open("bearing.curvature.normal.dat");
				break;
		}

		for (int i = 0; i < NumBins; i++){
			AccBearing[i]= 0.0;
			BearingBinCounter[i] = 0;
			//for (int j = 0; j < MaxPointsPerBin; j++){
			//	BearingBins[i][j] = 0.0;
			//}
			BearingBins[i] = 0.0;
		}
		for (int repetitions = 1; repetitions <= 1000; repetitions++)
		{
			// Run a trial
			k = 0;
			arrived = 0;
			tempTheta = 0.0;
			Worm.ResetAgentIntState(rs);
			Worm.ResetAgentsBody(rs);
			Worm.ResetChemCon(rs);
			Worm.UpdateChemCon(rs);
			Worm.InitialiseSensorHistory();
			for (double t = StepSize; t <= RunDuration; t += StepSize)
			{
				Worm.Step(StepSize, rs,t);
				Worm.UpdateChemCon(rs);
				Worm.UpdateChemConHistory();
				Worm.UpdateSensors();

				Xpoints[k] = Worm.PositionX();
				Ypoints[k] = Worm.PositionY();
				ThetaPoints[k] = Worm.Theta();
				k++;

				switch (condition){
					case 2:
						Worm.SetOnCell(0.0);
						break;
					case 3:
						Worm.SetOffCell(0.0);
						break;
					case 4:
						Worm.SetOnCell(0.0);
						Worm.SetOffCell(0.0);
						break;
					case 5:
						Worm.SetOffCell(Worm.OnCell());
						break;
					case 6:
						Worm.SetOnCell(Worm.OffCell());
						break;
					default:
						break;
				}

			}
			// Analyse the angles from that trial (if it arrived to the peak)
			// First calculate the inflexion points from the turning angle theta
			int j = 0;
			for (int i = 0; i < k; i++)
			{
				if (((ThetaPoints[i] < 0.0) and (ThetaPoints[i+1] > 0.0)) or
					((ThetaPoints[i] > 0.0) and (ThetaPoints[i+1] < 0.0)) ) {
					InflexionPoints[j] = i;
					j++;
				}
			}

			for (int l = 0; l < j; l+=2)
			{
				// Obtain the X and Y coordinates at the points of inflexion
				aX = Xpoints[InflexionPoints[l]];
				aY = Ypoints[InflexionPoints[l]];
				bX = Xpoints[InflexionPoints[l+1]];
				bY = Ypoints[InflexionPoints[l+1]];
				cX = Xpoints[InflexionPoints[l+2]];
				cY = Ypoints[InflexionPoints[l+2]];

				// Calculate the translational orientation of the worm
				if (cX - aX == 0.0){
					if (cY > aY)
						transOrient = Pi/2;
					else
						transOrient	= -Pi/2;
				}
				else
					transOrient = atan((cY - aY)/(cX - aX));
				transOrient = cX < aX ? transOrient + Pi: transOrient;
				transOrient = transOrient > Pi ? transOrient - 2*Pi: transOrient;

				// Calculate the peak orientation
				if (peakX - bX == 0.0)
					peakOrient = Pi/2;
				else
					peakOrient = atan((peakY - bY)/(peakX - bX));
				peakOrient = peakX < bX ? peakOrient + Pi: peakOrient;
				peakOrient = peakOrient > Pi ? peakOrient - 2*Pi: peakOrient;

				// Calculate the bearing (i.e., difference in angle between the translational orientation and the orientation of the peak).
				Bearing[l] = transOrient - peakOrient;
				if (Bearing[l] > Pi){
					Bearing[l] = Bearing[l] - 2*Pi;
				}
				if (Bearing[l] < -Pi) {
					Bearing[l] = Bearing[l]	+ 2*Pi;
				}
			}

			// Write on file the 2-tuple {Theta, Bearing}
			for (int l = 0; l < j; l+=2)
			{
				for (int m = InflexionPoints[l]; m < InflexionPoints[l+2]; m++){
					BearingDegree = Bearing[l]*(180/Pi);
					ThetaDegree = ThetaPoints[m]*(180/Pi);
					if ((-180.0 <= BearingDegree) and (BearingDegree < -175.0))
						bearingIndex = 0;
					if ((-175.0 <= BearingDegree) and (BearingDegree < -165.0))
						bearingIndex = 1;
					if ((-165.0 <= BearingDegree) and (BearingDegree < -155.0))
						bearingIndex = 2;
					if ((-155.0 <= BearingDegree) and (BearingDegree < -145.0))
						bearingIndex = 3;
					if ((-145.0 <= BearingDegree) and (BearingDegree < -135.0))
						bearingIndex = 4;
					if ((-135.0 <= BearingDegree) and (BearingDegree < -125.0))
						bearingIndex = 5;
					if ((-125.0 <= BearingDegree) and (BearingDegree < -115.0))
						bearingIndex = 6;
					if ((-115.0 <= BearingDegree) and (BearingDegree < -105.0))
						bearingIndex = 7;
					if ((-105.0 <= BearingDegree) and (BearingDegree < -95.0))
						bearingIndex = 8;
					if ((-95.0 <= BearingDegree) and (BearingDegree < -85.0))
						bearingIndex = 9;
					if ((-85.0 <= BearingDegree) and (BearingDegree < -75.0))
						bearingIndex = 10;
					if ((-75.0 <= BearingDegree) and (BearingDegree < -65.0))
						bearingIndex = 11;
					if ((-65.0 <= BearingDegree) and (BearingDegree < -55.0))
						bearingIndex = 12;
					if ((-55.0 <= BearingDegree) and (BearingDegree < -45.0))
						bearingIndex = 13;
					if ((-45.0 <= BearingDegree) and (BearingDegree < -35.0))
						bearingIndex = 14;
					if ((-35.0 <= BearingDegree) and (BearingDegree < -25.0))
						bearingIndex = 15;
					if ((-25.0 <= BearingDegree) and (BearingDegree < -15.0))
						bearingIndex = 16;
					if ((-15.0 <= BearingDegree) and (BearingDegree < -5.0))
						bearingIndex = 17;
					if ((-5.0 <= BearingDegree) and (BearingDegree < 5.0))
						bearingIndex = 18;
					if ((5.0 <= BearingDegree) and (BearingDegree < 15.0))
						bearingIndex = 19;
					if ((15.0 <= BearingDegree) and (BearingDegree < 25.0))
						bearingIndex = 20;
					if ((25.0 <= BearingDegree) and (BearingDegree < 35.0))
						bearingIndex = 21;
					if ((35.0 <= BearingDegree) and (BearingDegree < 45.0))
						bearingIndex = 22;
					if ((45.0 <= BearingDegree) and (BearingDegree < 55.0))
						bearingIndex = 23;
					if ((55.0 <= BearingDegree) and (BearingDegree < 65.0))
						bearingIndex = 24;
					if ((65.0 <= BearingDegree) and (BearingDegree < 75.0))
						bearingIndex = 25;
					if ((75.0 <= BearingDegree) and (BearingDegree < 85.0))
						bearingIndex = 26;
					if ((85.0 <= BearingDegree) and (BearingDegree < 95.0))
						bearingIndex = 27;
					if ((95.0 <= BearingDegree) and (BearingDegree < 105.0))
						bearingIndex = 28;
					if ((105.0 <= BearingDegree) and (BearingDegree < 115.0))
						bearingIndex = 29;
					if ((115.0 <= BearingDegree) and (BearingDegree < 125.0))
						bearingIndex = 30;
					if ((125.0 <= BearingDegree) and (BearingDegree < 135.0))
						bearingIndex = 31;
					if ((135.0 <= BearingDegree) and (BearingDegree < 145.0))
						bearingIndex = 32;
					if ((145.0 <= BearingDegree) and (BearingDegree < 155.0))
						bearingIndex = 33;
					if ((155.0 <= BearingDegree) and (BearingDegree < 165.0))
						bearingIndex = 34;
					if ((165.0 <= BearingDegree) and (BearingDegree < 175.0))
						bearingIndex = 35;
					if ((175.0 <= BearingDegree) and (BearingDegree < 180.0))
						bearingIndex = 36;

					//if (BearingBinCounter[bearingIndex] < MaxPointsPerBin){
					AccBearing[bearingIndex] += ThetaDegree;
					//	BearingBins[bearingIndex][BearingBinCounter[bearingIndex]] = ThetaDegree;
					BearingBinCounter[bearingIndex] += 1;
					//}
				}
			}
		}
		for (int i = 0; i < NumBins; i++){
			if (BearingBinCounter[i] > 0){
				tempBearing = i*BinSize - 180.0;
				avgBearing = AccBearing[i]/BearingBinCounter[i];
				//stddev = 0.0;
				//for (int p = 0; p < BearingBinCounter[i]; p++){
				//	stddev = pow(BearingBins[i][p] - avgBearing, 2);
				//}
				//stddev = sqrt(stddev / (BearingBinCounter[i] - 1));
				//ExampleFile << tempBearing << " " << avgBearing << " " << stddev << endl;
				ExampleFile << tempBearing << " " << avgBearing << endl;
			}
		}
		ExampleFile.close();
	}
}

// ------------------------------------
// Constrained body analysis
// ------------------------------------
void Riverchip()
{
	ofstream ExampleFile,IndexFile;
	double t,c,fitnessB;
	double turningBias_N = 0,turningBias_DP = 0.0,turningBias_DN = 0.0,turningBias_VP = 0.0,turningBias_VN = 0.0;
	int repetitions,riverReps=1,k,condition,length;
	RandomState rs;
	long IDUM=time(0);
	rs.SetRandomSeed(IDUM);
	WormAgent Worm("best.ns.dat");
	ExampleFile.open("micro.dat");
	IndexFile.open("riverchipindex.dat");
	Worm.ResetAgentIntState(rs);
	Worm.InitialiseAgent(StepSize);
	Worm.InitialiseSensorHistory();
	Worm.ResetAgentsBody(rs);
	Worm.SetOffsetCPG(0.0);
	Worm.SetOrientation(0.0);
	Worm.ResetChemCon(rs);
	length = (int) ((7*HST)/StepSize);
	t = 0.0;

	for (k = 0; k < length; k++)
	{
		Worm.Step(StepSize, rs, t);
		Worm.SetChemCon(0.0);
		Worm.UpdateChemConHistory();
		Worm.UpdateSensors();
		t += StepSize;
	}
	for (repetitions = 1; repetitions <= riverReps; repetitions++){
		for (condition = 0; condition < 8; condition++){
			for (k = 0; k < length; k++)
			{
#ifdef ONKILL
				Worm.SetOnCell(0.0);
#endif
#ifdef OFFKILL
				Worm.SetOffCell(0.0);
#endif
				Worm.Step(StepSize, rs, t);
				switch (condition) {
					case 1:
						turningBias_DP += Worm.Theta();
						c = SizeOfStep*(1/(1+exp(-RiverChipSteepness*Worm.Theta())));
						break;
					case 3:
						turningBias_DN += Worm.Theta();
						c = - SizeOfStep*(1/(1+exp(-RiverChipSteepness*Worm.Theta())));
						break;
					case 5:
						turningBias_VP += Worm.Theta();
						c = SizeOfStep*(1/(1+exp(RiverChipSteepness*Worm.Theta())));
						break;
					case 7:
						turningBias_VN += Worm.Theta();
						c = - SizeOfStep*(1/(1+exp(RiverChipSteepness*Worm.Theta())));
						break;
					default:
						turningBias_N += Worm.Theta();
						c = 0.0;
						break;
				}
				Worm.SetChemCon(c);
				Worm.UpdateChemConHistory();
				Worm.UpdateSensors();
				t += StepSize;

				ExampleFile << Worm.ChemCon() << " " << Worm.OnCell() << " " << Worm.OffCell() << " " << Worm.PositionX() << " " << Worm.PositionY() << " " << Worm.Theta() << " " << Worm.avgvel << " " << Worm.NervousSystem.NeuronOutput(1) << " " << Worm.NervousSystem.NeuronOutput(2) << endl;
			}
		}
	}

	turningBias_N = turningBias_N/(4*length*riverReps);
	turningBias_DP = turningBias_DP/(length*riverReps);
	turningBias_DN = turningBias_DN/(length*riverReps);
	turningBias_VP = turningBias_VP/(length*riverReps);
	turningBias_VN = turningBias_VN/(length*riverReps);
	fitnessB = (turningBias_DP + turningBias_VN) - (turningBias_DN + turningBias_VP) - fabs(turningBias_N);

	IndexFile << turningBias_N << " " << turningBias_DP << " " << turningBias_DN << " " << turningBias_VP << " " << turningBias_VN << " " << fitnessB << endl;
	ExampleFile.close();
	IndexFile.close();
}

// ------------------------------------
// Sensory cell dynamics
// ------------------------------------
void ConcStepSensoryDyn()
{
	ofstream ExampleFile;
	long IDUM=time(0);
	RandomState rs;
	rs.SetRandomSeed(IDUM);
	ExampleFile.open("downstep.0.1.dat");
	WormAgent Worm("best.ns.dat");
	Worm.SetSensorN(1.0);
	Worm.SetSensorM(2.0);
	Worm.SetSensorD(0.0);
	Worm.InitialiseAgent(StepSize);
	Worm.InitialiseSensorHistory();
	Worm.ResetAgentIntState(rs);
	Worm.ResetAgentsBody(rs);
	Worm.ResetChemCon(rs);
	Worm.SetChemCon(0.0);
	Worm.UpdateChemConHistory();
	// Initial phase
	for (double t = StepSize; t <= 10.0; t += StepSize){
		Worm.Step(StepSize, rs,t);
		Worm.SetChemCon(0.0);
		Worm.UpdateChemConHistory();
		Worm.UpdateSensors();
		ExampleFile << Worm.ChemCon() << " " << Worm.OffCell() << endl;
	}

	double c = 0.0;
	for (int i=0; i<1; i++){
		// Initial phase
		c -= 1.0;
		for (double t = StepSize; t <= 10.0; t += StepSize){
			Worm.Step(StepSize, rs,t);
			Worm.SetChemCon(c);
			Worm.UpdateChemConHistory();
			Worm.UpdateSensors();
			ExampleFile << Worm.ChemCon() << " " << Worm.OffCell() << endl;
		}
		// Stay up
		for (double t = StepSize; t <= 0.0; t += StepSize){
			Worm.Step(StepSize, rs,t);
			Worm.SetChemCon(0.0);
			Worm.UpdateChemConHistory();
			Worm.UpdateSensors();
			ExampleFile << Worm.ChemCon() << " " << Worm.OffCell() << endl;
		}
	}

	/*	for (double t = StepSize; t <= 10*Pi; t += StepSize){
	 Worm.Step(StepSize, rs,t);
	 Worm.SetChemCon(sin(t));
	 Worm.UpdateChemConHistory();
	 Worm.UpdateSensors();
	 ExampleFile << Worm.ChemCon() << " " << Worm.OnCell() << endl;
	 }
	 */
	for (double t = StepSize; t <= 0.0; t += StepSize){
		Worm.Step(StepSize, rs,t);
		Worm.SetChemCon(0.0);
		Worm.UpdateChemConHistory();
		Worm.UpdateSensors();
		ExampleFile << Worm.ChemCon() << " " << Worm.OffCell() << endl;
	}

	ExampleFile.close();
}

// ------------------------------------
// Dynamical Sytems analysis
// ------------------------------------
void AutonomousDynamics()
{
	ofstream OutputFile;
	double stepinc=0.01;
	double t,onWeight,offWeight,oscInput,selfWeight,oscWeight,bias,initActivation,outputGain,sensorN,sensorM,sensorD;
	double EPList[2],NumEPs; // Only two because that's the most one neuron with a self-connection can have
	double NMdiff,theta;

	// Read parameters from file
	ifstream InputFile;
	InputFile.open("quickgen.dat");

	InputFile >> onWeight;
	InputFile >> offWeight;
	InputFile >> oscWeight;
	InputFile >> selfWeight;
	InputFile >> sensorN;
	InputFile >> sensorM;
	InputFile >> sensorD;
	InputFile >> bias;
	InputFile >> outputGain;		// XXX

	// Initialise parameters of the neuron
	CTRNN MotorNeuron;
	MotorNeuron.SetCircuitSize(1);
	MotorNeuron.SetNeuronTimeConstant(1,0.1);
	MotorNeuron.SetConnectionWeight(1, 1, selfWeight);
	MotorNeuron.SetNeuronBias(1, bias);

	CTRNN MotorNeuronD;
	MotorNeuronD.SetCircuitSize(1);
	MotorNeuronD.SetNeuronTimeConstant(1,0.1);
	MotorNeuronD.SetConnectionWeight(1, 1, selfWeight);
	MotorNeuronD.SetNeuronBias(1, bias);

	CTRNN MotorNeuronV;
	MotorNeuronV.SetCircuitSize(1);
	MotorNeuronV.SetNeuronTimeConstant(1,0.1);
	MotorNeuronV.SetConnectionWeight(1, 1, selfWeight);
	MotorNeuronV.SetNeuronBias(1, bias);

	//---------------------------------------
	// Steady-state given different external input from the Oscillator from its full range
	OutputFile.open("autonomousdynamics.dat");
	for (oscInput = -5.0;  oscInput <= 5.0 + stepinc; oscInput += stepinc)
	{
		NumEPs = 0.0;
		MotorNeuron.SetNeuronExternalInput(1, oscWeight*oscInput);
		for (initActivation= -20.0; initActivation <= 20.0; initActivation += 5.0)
		{
			MotorNeuron.SetNeuronState(1,initActivation);
			for (t = 0.0; t <= 10.0; t+=StepSize)
				MotorNeuron.EulerStep(StepSize);

			if (NumEPs==0){
				// If no EP recorded, save first.
				EPList[0] = MotorNeuron.NeuronOutput(1);
				NumEPs+=1.0;
			}
			else{
				// Check if it is different from the one that's already recorded, if so record it.
				if (fabs(MotorNeuron.NeuronOutput(1) - EPList[0]) > 0.01){
					if ((NumEPs==2) and (fabs(MotorNeuron.NeuronOutput(1) - EPList[1]) > 0.01)){
						cout << "Error" << endl;
					}
					else{
						EPList[1] = MotorNeuron.NeuronOutput(1);
						NumEPs=2.0;
					}
				}
			}
		}
		for (int i=0; i<NumEPs; i++){
			OutputFile << oscInput << " " << EPList[i] << " " << endl;
		}
	}
	OutputFile.close();

	//---------------------------------------
	// Dynamics given the driven oscillations
	OutputFile.open("drivendynamics.dat");
	MotorNeuronD.SetNeuronState(1,initActivation);
	MotorNeuronV.SetNeuronState(1,initActivation);
	theta = 0.0;
	for (t = 0.0; t < 10.0; t+=StepSize){
		MotorNeuronD.EulerStep(StepSize);
		MotorNeuronD.SetNeuronExternalInput(1, oscWeight*sin(0 + HSP*t));
		MotorNeuronV.EulerStep(StepSize);
		MotorNeuronV.SetNeuronExternalInput(1, oscWeight*sin(Pi + HSP*t));
		NMdiff = MotorNeuronV.NeuronOutput(1) - MotorNeuronD.NeuronOutput(1);
		theta += StepSize * (outputGain * NMdiff);
	}
	for (t = 10.0; t <= 20.0; t+=StepSize){
		MotorNeuronD.EulerStep(StepSize);
		MotorNeuronD.SetNeuronExternalInput(1, oscWeight*sin(0 + HSP*t));
		MotorNeuronV.EulerStep(StepSize);
		MotorNeuronV.SetNeuronExternalInput(1, oscWeight*sin(Pi + HSP*t));
		NMdiff = MotorNeuronV.NeuronOutput(1) - MotorNeuronD.NeuronOutput(1);
		theta += StepSize * (outputGain * NMdiff);
		OutputFile << sin(HSP*t) << " " << MotorNeuronD.NeuronOutput(1) << " " << MotorNeuronV.NeuronOutput(1) << " " << NMdiff << " " << theta << " " << endl;
	}
	OutputFile.close();

	//---------------------------------------
	// Dynamics given the driven oscillations and an ON signal
	// Von max 0.0277096
	/*
	 double currentConc=0,pastConc=0,tempDiff,V_on,V_off;
	 double chemConHistory[];
	 int timer=0;
	 OutputFile.open("drivendynamicsON.dat");
	 MotorNeuronD.SetNeuronState(1,initActivation);
	 MotorNeuronV.SetNeuronState(1,initActivation);
	 theta = 0.0;
	 for (t = 0.0; t < 10.0; t+=StepSize){
	 // Update chemical concentration

	 // Update sensory cells
	 currentConc += chemConHistory[timer - iSensorD - 1] - chemConHistory[timer - iSensorD - iSensorN - 1];
	 pastConc += chemConHistory[timer - iSensorD - iSensorN - 1] - chemConHistory[timer - iSensorD - iSensorN - iSensorM - 1];
	 tempDiff = (currentConc/dSensorN) - (pastConc/dSensorM);
	 V_on = tempDiff > 0.0 ? tempDiff: 0.0;
	 V_off = tempDiff < 0.0 ? fabs(tempDiff): 0.0;

	 // Update motor neurons
	 MotorNeuronD.EulerStep(StepSize);
	 MotorNeuronD.SetNeuronExternalInput(1, onWeight*V_on + oscWeight*sin(0 + HSP*t));
	 MotorNeuronV.EulerStep(StepSize);
	 MotorNeuronV.SetNeuronExternalInput(1, onWeight*V_on + oscWeight*sin(Pi + HSP*t));
	 NMdiff = MotorNeuronV.NeuronOutput(1) - MotorNeuronD.NeuronOutput(1);
	 theta += StepSize * (outputGain * NMdiff);
	 }
	 for (t = 10.0; t <= 20.0; t+=StepSize){
	 MotorNeuronD.EulerStep(StepSize);
	 MotorNeuronD.SetNeuronExternalInput(1, onWeight*V_on + oscWeight*sin(0 + HSP*t));
	 MotorNeuronV.EulerStep(StepSize);
	 MotorNeuronV.SetNeuronExternalInput(1, onWeight*V_on + oscWeight*sin(Pi + HSP*t));
	 NMdiff = MotorNeuronV.NeuronOutput(1) - MotorNeuronD.NeuronOutput(1);
	 theta += StepSize * (outputGain * NMdiff);
	 OutputFile << sin(HSP*t) << " " << MotorNeuronD.NeuronOutput(1) << " " << MotorNeuronV.NeuronOutput(1) << " " << NMdiff << " " << theta << " " << endl;
	 }
	 OutputFile.close();

	 //---------------------------------------
	 // Dynamics given the driven oscillations and an ON signal
	 // Voff max 0.0277941
	 OutputFile.open("drivendynamicsOFF.dat");
	 MotorNeuronD.SetNeuronState(1,initActivation);
	 MotorNeuronV.SetNeuronState(1,initActivation);
	 theta = 0.0;
	 for (t = 0.0; t < 10.0; t+=StepSize){
	 MotorNeuronD.EulerStep(StepSize);
	 MotorNeuronD.SetNeuronExternalInput(1, offWeight*0.0277941 + oscWeight*sin(0 + HSP*t));
	 MotorNeuronV.EulerStep(StepSize);
	 MotorNeuronV.SetNeuronExternalInput(1, offWeight*0.0277941 + oscWeight*sin(Pi + HSP*t));
	 NMdiff = MotorNeuronV.NeuronOutput(1) - MotorNeuronD.NeuronOutput(1);
	 theta += StepSize * (outputGain * NMdiff);
	 }
	 for (t = 10.0; t <= 20.0; t+=StepSize){
	 MotorNeuronD.EulerStep(StepSize);
	 MotorNeuronD.SetNeuronExternalInput(1, offWeight*0.0277941 + oscWeight*sin(0 + HSP*t));
	 MotorNeuronV.EulerStep(StepSize);
	 MotorNeuronV.SetNeuronExternalInput(1, offWeight*0.0277941 + oscWeight*sin(Pi + HSP*t));
	 NMdiff = MotorNeuronV.NeuronOutput(1) - MotorNeuronD.NeuronOutput(1);
	 theta += StepSize * (outputGain * NMdiff);
	 OutputFile << sin(HSP*t) << " " << MotorNeuronD.NeuronOutput(1) << " " << MotorNeuronV.NeuronOutput(1) << " " << NMdiff << " " << theta << " " << endl;
	 }
	 OutputFile.close();
	 */
}
// ------------------------------------
// Phase-sensititivy analysis
// ------------------------------------
void PhaseSensitivityExample()
{
	ofstream ExampleFile;
	double t,conc,concStep;
	long IDUM=time(0);
	double avgtheta,k,phasedelay;
	double duration = 5*HST;

	RandomState rs;
	rs.SetRandomSeed(IDUM);
	WormAgent Worm("best.ns.dat");
	Worm.InitialiseAgent(StepSize);

	// Get it to initialise
	Worm.ResetAgentIntState(rs);
	Worm.ResetAgentsBody(rs);
	Worm.SetChemCon(conc);
	Worm.InitialiseSensorHistory();
	Worm.SetOffsetCPG(0.0);
	Worm.SetOrientation(0.0);
	for (t = StepSize; t <= duration; t += StepSize)
	{
		Worm.Step(StepSize, rs,t);
		Worm.UpdateChemConHistory();
		Worm.UpdateSensors();
	}

	// No sensory input
	ExampleFile.open("exampleUpstepB.dat");
	concStep = 0.010;
	phasedelay = (3*HST)/4;
	Worm.ResetAgentIntState(rs);
	Worm.ResetAgentsBody(rs);
	Worm.SetChemCon(conc);
	Worm.InitialiseSensorHistory();
	Worm.SetOffsetCPG(0.0);
	Worm.SetOrientation(0.0);
	Worm.PrintDetail(ExampleFile);
	conc = 0.0;
	for (t = StepSize; t <= duration; t += StepSize)
	{
		Worm.Step(StepSize, rs,t);
		Worm.SetChemCon(conc); //Worm.UpdateChemCon(rs);
		if (fmod(t+phasedelay,HST) < StepSize){
			conc += concStep;
		}
		Worm.UpdateChemConHistory();
		Worm.UpdateSensors();
		Worm.PrintDetail(ExampleFile);
	}
	avgtheta = avgtheta/k;
	ExampleFile.close();
}

// ------------------------------------
// Phase-sensititivy analysis
// ------------------------------------
void PhaseSensitivity()
{
	ofstream ExampleFile;
	double t,conc,concStep;
	long IDUM=time(0);
	double avgtheta,k,phasedelay,actualtheta,actualorient;
	double iX,iY,fX,fY,transOrient,transOrientZ,diffOrient;
	double phasedelaystepsize=0.01;
	double duration = 10*HST;
	double startrecordingtime=5*HST;

	RandomState rs;
	rs.SetRandomSeed(IDUM);
	WormAgent Worm("best.ns.dat");
	Worm.InitialiseAgent(StepSize);

	// Get it to initialise
	Worm.ResetAgentIntState(rs);
	Worm.ResetAgentsBody(rs);
	Worm.SetChemCon(conc);
	Worm.InitialiseSensorHistory();
	Worm.SetOffsetCPG(0.0);
	Worm.SetOrientation(0.0);
	for (t = StepSize; t <= duration; t += StepSize)
	{
		Worm.Step(StepSize, rs,t);
		Worm.UpdateChemConHistory();
		Worm.UpdateSensors();
	}

	// No sensory input
	ExampleFile.open("phasesensNormal.dat");
	concStep = 0.0;
	for (phasedelay = 0.0; phasedelay <= HST; phasedelay+=phasedelaystepsize)
	{
		Worm.ResetAgentIntState(rs);
		Worm.ResetAgentsBody(rs);
		Worm.SetChemCon(conc);
		Worm.InitialiseSensorHistory();
		Worm.SetOffsetCPG(0.0);
		Worm.SetOrientation(0.0);
		conc = 0.0;
		avgtheta = 0.0;
		k=0;
		iX = Worm.PositionX();
		iY = Worm.PositionY();
		for (t = StepSize; t <= duration; t += StepSize)
		{
			Worm.Step(StepSize, rs,t);
			Worm.SetChemCon(conc); //Worm.UpdateChemCon(rs);
			if (fmod(t+phasedelay,HST) < StepSize){
				//cout << "	" << t << " " << phasedelay << endl;
				conc += concStep;
				actualtheta = Worm.Theta();
				actualorient = Worm.Orientation();
			}
			Worm.UpdateChemConHistory();
			Worm.UpdateSensors();
			avgtheta += Worm.Theta();
			k++;
		}
		avgtheta = avgtheta/k;
		fX = Worm.PositionX();
		fY = Worm.PositionY();
		// Calculate the translational orientation of the worm
		if (fX - iX == 0.0){
			if (fY > iY)
				transOrient = Pi/2;
			else
				transOrient	= -Pi/2;
		}
		else
			transOrient = atan((fY - iY)/(fX - iX));
		transOrient = fX < iX ? transOrient + Pi: transOrient;
		transOrient = transOrient > Pi ? transOrient - 2*Pi: transOrient;
		//cout << iX << " " << iY << " " << fX << " " << fY << " " << transOrient << endl;
		ExampleFile << phasedelay*HSP << " " << avgtheta << " " << actualtheta << " " << actualorient << " " << transOrient << endl;
	}
	ExampleFile.close();

	// Up steps
	ExampleFile.open("phasesensDownstep.0005.dat");
	concStep = -0.0005;
	for (phasedelay = 0.0; phasedelay <= HST; phasedelay+=0.1)
	{
		Worm.ResetAgentIntState(rs);
		Worm.ResetAgentsBody(rs);
		conc = 0.0;
		Worm.SetChemCon(conc);
		Worm.InitialiseSensorHistory();
		Worm.SetOffsetCPG(0.0);
		Worm.SetOrientation(0.0);
		avgtheta = 0.0;
		k=0;
		iX = Worm.PositionX();
		iY = Worm.PositionY();
		for (t = StepSize; t <= duration; t += StepSize)
		{
			Worm.Step(StepSize, rs,t);
			Worm.SetChemCon(conc); //Worm.UpdateChemCon(rs);
			if (fmod(t+phasedelay,HST) < StepSize){
				conc += concStep;
				actualtheta = Worm.Theta();
				actualorient = Worm.Orientation();
			}
			Worm.UpdateChemConHistory();
			Worm.UpdateSensors();
			if (t > startrecordingtime){
				avgtheta += Worm.Theta();
				k++;
			}
		}
		avgtheta = avgtheta/k;
		fX = Worm.PositionX();
		fY = Worm.PositionY();

		// Calculate the translational orientation of the worm
		if (fX - iX == 0.0){
			if (fY > iY)
				transOrientZ = Pi/2;
			else
				transOrientZ	= -Pi/2;
		}
		else
			transOrientZ = atan((fY - iY)/(fX - iX));
		transOrientZ = fX < iX ? transOrientZ + Pi: transOrientZ;
		transOrientZ = transOrientZ > Pi ? transOrientZ - 2*Pi: transOrientZ;

		// Calculate the difference (i.e., difference in angle between the translational orientation and the orientation of the peak).
		diffOrient = transOrientZ - transOrient;
		if (diffOrient > Pi){
			diffOrient = diffOrient - 2*Pi;
		}
		if (diffOrient < -Pi) {
			diffOrient = diffOrient	+ 2*Pi;
		}
		ExampleFile << phasedelay*HSP << " " << avgtheta << " " << actualtheta << " " << actualorient << " " << diffOrient << endl;
	}
	ExampleFile.close();
	// Up steps
	ExampleFile.open("phasesensDownstep.001.dat");
	concStep = -0.001;
	for (phasedelay = 0.0; phasedelay <= HST; phasedelay+=0.1)
	{
		Worm.ResetAgentIntState(rs);
		Worm.ResetAgentsBody(rs);
		conc = 0.0;
		Worm.SetChemCon(conc);
		Worm.InitialiseSensorHistory();
		Worm.SetOffsetCPG(0.0);
		Worm.SetOrientation(0.0);
		avgtheta = 0.0;
		k=0;
		iX = Worm.PositionX();
		iY = Worm.PositionY();
		for (t = StepSize; t <= duration; t += StepSize)
		{
			Worm.Step(StepSize, rs,t);
			Worm.SetChemCon(conc); //Worm.UpdateChemCon(rs);
			if (fmod(t+phasedelay,HST) < StepSize){
				conc += concStep;
				actualtheta = Worm.Theta();
				actualorient = Worm.Orientation();
			}
			Worm.UpdateChemConHistory();
			Worm.UpdateSensors();
			if (t > startrecordingtime){
				avgtheta += Worm.Theta();
				k++;
			}
		}
		avgtheta = avgtheta/k;
		fX = Worm.PositionX();
		fY = Worm.PositionY();

		// Calculate the translational orientation of the worm
		if (fX - iX == 0.0){
			if (fY > iY)
				transOrientZ = Pi/2;
			else
				transOrientZ	= -Pi/2;
		}
		else
			transOrientZ = atan((fY - iY)/(fX - iX));
		transOrientZ = fX < iX ? transOrientZ + Pi: transOrientZ;
		transOrientZ = transOrientZ > Pi ? transOrientZ - 2*Pi: transOrientZ;

		// Calculate the difference (i.e., difference in angle between the translational orientation and the orientation of the peak).
		diffOrient = transOrientZ - transOrient;
		if (diffOrient > Pi){
			diffOrient = diffOrient - 2*Pi;
		}
		if (diffOrient < -Pi) {
			diffOrient = diffOrient	+ 2*Pi;
		}
		ExampleFile << phasedelay*HSP << " " << avgtheta << " " << actualtheta << " " << actualorient << " " << diffOrient << endl;
	}
	ExampleFile.close();

	// Up steps
	ExampleFile.open("phasesensDownstep.002.dat");
	concStep = -0.002;
	for (phasedelay = 0.0; phasedelay <= HST; phasedelay+=0.1)
	{
		Worm.ResetAgentIntState(rs);
		Worm.ResetAgentsBody(rs);
		conc = 0.0;
		Worm.SetChemCon(conc);
		Worm.InitialiseSensorHistory();
		Worm.SetOffsetCPG(0.0);
		Worm.SetOrientation(0.0);
		avgtheta = 0.0;
		k=0;
		iX = Worm.PositionX();
		iY = Worm.PositionY();
		for (t = StepSize; t <= duration; t += StepSize)
		{
			Worm.Step(StepSize, rs,t);
			Worm.SetChemCon(conc); //Worm.UpdateChemCon(rs);
			if (fmod(t+phasedelay,HST) < StepSize){
				conc += concStep;
				actualtheta = Worm.Theta();
				actualorient = Worm.Orientation();
			}
			Worm.UpdateChemConHistory();
			Worm.UpdateSensors();
			if (t > startrecordingtime){
				avgtheta += Worm.Theta();
				k++;
			}
		}
		avgtheta = avgtheta/k;
		fX = Worm.PositionX();
		fY = Worm.PositionY();

		// Calculate the translational orientation of the worm
		if (fX - iX == 0.0){
			if (fY > iY)
				transOrientZ = Pi/2;
			else
				transOrientZ	= -Pi/2;
		}
		else
			transOrientZ = atan((fY - iY)/(fX - iX));
		transOrientZ = fX < iX ? transOrientZ + Pi: transOrientZ;
		transOrientZ = transOrientZ > Pi ? transOrientZ - 2*Pi: transOrientZ;

		// Calculate the difference (i.e., difference in angle between the translational orientation and the orientation of the peak).
		diffOrient = transOrientZ - transOrient;
		if (diffOrient > Pi){
			diffOrient = diffOrient - 2*Pi;
		}
		if (diffOrient < -Pi) {
			diffOrient = diffOrient	+ 2*Pi;
		}
		ExampleFile << phasedelay*HSP << " " << avgtheta << " " << actualtheta << " " << actualorient << " " << diffOrient << endl;
	}
	ExampleFile.close();

	// Up steps
	ExampleFile.open("phasesensDownstep.003.dat");
	concStep = -0.003;
	for (phasedelay = 0.0; phasedelay <= HST; phasedelay+=0.1)
	{
		Worm.ResetAgentIntState(rs);
		Worm.ResetAgentsBody(rs);
		conc = 0.0;
		Worm.SetChemCon(conc);
		Worm.InitialiseSensorHistory();
		Worm.SetOffsetCPG(0.0);
		Worm.SetOrientation(0.0);
		avgtheta = 0.0;
		k=0;
		iX = Worm.PositionX();
		iY = Worm.PositionY();
		for (t = StepSize; t <= duration; t += StepSize)
		{
			Worm.Step(StepSize, rs,t);
			Worm.SetChemCon(conc); //Worm.UpdateChemCon(rs);
			if (fmod(t+phasedelay,HST) < StepSize){
				conc += concStep;
				actualtheta = Worm.Theta();
				actualorient = Worm.Orientation();
			}
			Worm.UpdateChemConHistory();
			Worm.UpdateSensors();
			if (t > startrecordingtime){
				avgtheta += Worm.Theta();
				k++;
			}
		}
		avgtheta = avgtheta/k;
		fX = Worm.PositionX();
		fY = Worm.PositionY();

		// Calculate the translational orientation of the worm
		if (fX - iX == 0.0){
			if (fY > iY)
				transOrientZ = Pi/2;
			else
				transOrientZ	= -Pi/2;
		}
		else
			transOrientZ = atan((fY - iY)/(fX - iX));
		transOrientZ = fX < iX ? transOrientZ + Pi: transOrientZ;
		transOrientZ = transOrientZ > Pi ? transOrientZ - 2*Pi: transOrientZ;

		// Calculate the difference (i.e., difference in angle between the translational orientation and the orientation of the peak).
		diffOrient = transOrientZ - transOrient;
		if (diffOrient > Pi){
			diffOrient = diffOrient - 2*Pi;
		}
		if (diffOrient < -Pi) {
			diffOrient = diffOrient	+ 2*Pi;
		}
		ExampleFile << phasedelay*HSP << " " << avgtheta << " " << actualtheta << " " << actualorient << " " << diffOrient << endl;
	}
	ExampleFile.close();

}

// ------------------------------------
// Parameter space exploration
// ------------------------------------
void ExploreParameterSpace()
{
	double bias,selfconn;
	double t,avgtime,timetopeak,orient,timespeaked,totalrepeats,avg;
	double dist;
	int reached;
	ofstream TimeFile,ReachedFile,DistFile;
	long IDUM=time(0);
	RandomState rs;
	rs.SetRandomSeed(IDUM);

	WormAgent Worm(CircuitSize);

	//---------------------------- Fixed parameters
#ifdef PSON
	Worm.SetOnWeight(1, 2.5);
	Worm.SetOnWeight(2, 2.5);
	Worm.SetOffWeight(1, 0.0);
	Worm.SetOffWeight(2, 0.0);
	TimeFile.open("ps_time_on2.5_off0_osc12.5.dat");
	ReachedFile.open("ps_reached_on2.5_off0_osc12.5.dat");
	DistFile.open("ps_dist_on2.5_off0_osc12.5.dat");
#endif

#ifdef PSOFF
	Worm.SetOnWeight(1, 0.0);
	Worm.SetOnWeight(2, 0.0);
	Worm.SetOffWeight(1, -2.5);
	Worm.SetOffWeight(2, -2.5);
	TimeFile.open("ps_time_on0_off-2.5_osc12.5.dat");
	ReachedFile.open("ps_reached_on0_off-2.5_osc12.5.dat");
	DistFile.open("ps_dist_on0_off-2.5_osc12.5.dat");
#endif

#ifdef PSBOTH
	Worm.SetOnWeight(1, 2.5);
	Worm.SetOnWeight(2, 2.5);
	Worm.SetOffWeight(1, -2.5);
	Worm.SetOffWeight(2, -2.5);
	TimeFile.open("ps_time_on2.5_off-2.5_osc12.5.dat");
	ReachedFile.open("ps_reached_on2.5_off-2.5_osc12.5.dat");
	DistFile.open("ps_dist_on2.5_off-2.5_osc12.5.dat");
#endif

	//--------
	Worm.SetCpgWeight(1, 12.5);
	Worm.SetCpgWeight(2, 12.5);
	//--------
	Worm.NervousSystem.SetConnectionWeight(1, 2, 0.0);
	Worm.NervousSystem.SetConnectionWeight(2, 1, 0.0);
	Worm.NervousSystem.SetNeuronTimeConstant(1, 0.1);
	Worm.NervousSystem.SetNeuronTimeConstant(2, 0.1);
	Worm.SetSensorN(1.0);
	Worm.SetOutputGain(1.0);
	Worm.InitialiseAgent(StepSize);
	Worm.InitialiseSensorHistory();
	//-------------------------------------------

	for (bias = -15.0; bias <= 15.0; bias += 0.5)
	{
		Worm.NervousSystem.SetNeuronBias(1, bias);
		Worm.NervousSystem.SetNeuronBias(2, bias);
		for (selfconn = -15.0; selfconn <= 15.0; selfconn += 0.5)
		{
			Worm.NervousSystem.SetConnectionWeight(1, 1, selfconn);
			Worm.NervousSystem.SetConnectionWeight(2, 2, selfconn);
			timespeaked = 0.0;
			totalrepeats = 0.0;
			avg = 0.0;
			dist = 0.0;

			//-----------
			for (int repetitions = 1; repetitions <= 1; repetitions++)
			{
				for (orient = 0.0;  orient <= 2*Pi; orient += (2*Pi)/180)
				{
					Worm.ResetAgentIntState(rs);
					Worm.ResetAgentsBody(rs);
					Worm.SetOrientation(orient);
					Worm.ResetChemCon(rs);

					t = StepSize;
					reached = 0;
					while (t <= RunDuration)
					{
						Worm.Step(StepSize, rs, t);
						Worm.UpdateChemCon(rs);
						Worm.UpdateSensors();
						dist += Worm.DistanceToCentre();
						if (Worm.DistanceToCentre() <= CloseEnoughRadius){
							reached = 1;
							timetopeak = t;
						}
						t += StepSize;
					}
					if (reached){
						avg += timetopeak;
						timespeaked++;
					}
					totalrepeats++;
				}
			}
			avgtime = timespeaked==0.0?0.0:avg/timespeaked;	// Avoid NAN
			TimeFile << avgtime << " ";
			ReachedFile << timespeaked/totalrepeats << " ";
			DistFile << dist/totalrepeats << " ";
		}
		TimeFile << endl;
		ReachedFile << endl;
		DistFile << endl;
	}
	TimeFile.close();
	ReachedFile.close();
	DistFile.close();
}

// ------------------------------------
// Display functions
// ------------------------------------
void EvolutionaryRunDisplay(int Generation, double BestPerf, double AvgPerf, double PerfVar)
{
	cout << Generation << " " << BestPerf << " " << AvgPerf << " " << PerfVar << endl;
}

void ResultsDisplay(TSearch &s)
{
	double p;
	TVector<double> bestVector;
	ofstream BestIndividualFile;
	TVector<double> phenotype;
	phenotype.SetBounds(1,VectSize);
	// Save the genotype of the best individual
	bestVector = s.BestIndividual();
	BestIndividualFile.open("best.gen.dat");
	BestIndividualFile << bestVector << endl;
	BestIndividualFile.close();
	// Also show the best individual in the Circuit Model form
	BestIndividualFile.open("best.ns.dat");
	GenPhenMapping(bestVector, phenotype);
	WormAgent Worm(CircuitSize);
	Worm.SetParameters(phenotype);
	Worm.InitialiseAgent(StepSize);
	BestIndividualFile << Worm.NervousSystem;
	BestIndividualFile << endl;
	for (int i = 1; i <= CircuitSize; i++){
		p = Worm.onWeight(i);
		BestIndividualFile << p << " ";
	}
	BestIndividualFile << endl;
	for (int i = 1; i <= CircuitSize; i++){
		p = Worm.offWeight(i);
		BestIndividualFile << p << " ";
	}
	BestIndividualFile << endl;
	for (int i = 1; i <= CircuitSize; i++){
		p = Worm.cpgWeight(i);
		BestIndividualFile << p << " ";
	}
	BestIndividualFile << endl;
	BestIndividualFile << Worm.SensorN() << " " << endl;
	BestIndividualFile << Worm.SensorM() << " " << endl;
	BestIndividualFile << Worm.SensorD() << " " << endl;
	BestIndividualFile << Worm.OutputGain() << " " << endl;
	BestIndividualFile.close();
	BestIndividualFile.open("quickgen.dat");
	BestIndividualFile << Worm.onWeight(1) << " " << Worm.offWeight(1) << " " << Worm.cpgWeight(1) << " " << Worm.NervousSystem.ConnectionWeight(1,1) << " " << Worm.NervousSystem.NeuronBias(1) << " " << Worm.SensorN() << " " << Worm.SensorM() << " " << Worm.SensorD() << " " << Worm.OutputGain() << endl;
	BestIndividualFile.close();

	//BehavioralOverview();
}

// ------------------------------------
// The main program
// ------------------------------------
#ifdef EVOLVE
int main (int argc, const char* argv[])
{
	long IDUM=-time(0);
	TSearch s(VectSize);
	TVector<double> phenotype;
	phenotype.SetBounds(1,VectSize);

	// redirect standard output to a file
#ifdef PRINTTOFILE
	ofstream file;
	file.open ("evol.dat");
	cout.rdbuf(file.rdbuf());
#endif

	// Configure the search
	s.SetRandomSeed(IDUM);
	s.SetPopulationStatisticsDisplayFunction(EvolutionaryRunDisplay);
	s.SetSearchResultsDisplayFunction(ResultsDisplay);
	s.SetSelectionMode(RANK_BASED);			//{FITNESS_PROPORTIONATE,RANK_BASED}
	s.SetReproductionMode(HILL_CLIMBING);	// {HILL_CLIMBING, GENETIC_ALGORITHM}
	s.SetPopulationSize(10);
	s.SetMaxGenerations(100);
	s.SetMutationVariance(0.05);
	s.SetCrossoverProbability(0.5);
	s.SetCrossoverMode(TWO_POINT);			//{UNIFORM, TWO_POINT}
	s.SetMaxExpectedOffspring(1.1);
	s.SetElitistFraction(0.0);				//Default is 0.0.
	s.SetSearchConstraint(1);
	s.SetCheckpointInterval(0);
	s.SetReEvaluationFlag(1);				// CRUCIAL

	s.SetEvaluationFunction(EvaluationFunction);
	s.ExecuteSearch();

	return 0;
}
#else
int main (int argc, const char* argv[])
{
	ExampleRun();
	//BehavioralOverview();
	//AccConcentration();
	//BehavioralAnalysis();
	//OrientationAnalysis();
	//BearingCurvature();

	//AutonomousDynamics();
	//PhaseSensitivityExample();
	//PhaseSensitivity();
	//ExploreParameterSpace();
	//ConcStepSensoryDyn();
	//Riverchip();
	return 0;
}
#endif
