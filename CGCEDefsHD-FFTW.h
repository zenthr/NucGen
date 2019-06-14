//Array
#define ARRAY 1201 // (hCG/h)*(nCG-1)+1
#define HALF 600
#define OFF  26
#define MOLDARRAY 2001
#define MOLDHALF 1000
#define MOLDRANGE 1000
#define etamax 10
//Consts
# define PIE 3.141592653589793
# define GEVFM 0.1975
# define NUMSTEP 1025  // 2^10+1
#define NG 8  //No. Gluons
#define NC 3  //No. Colors
#define MinArray 104
//#define C 0.951237
#define XSect 7 //LHC NN Cross Section in fm^2
#define C 1.0

#define Order 2


const int N = 600; //Iterations for corelations

//CDF Limits
const double RhoMax = 5.0;

//CG Quantities
const double g = 1.0;
const double mIR = 1.01; // m, the InfraRed cutoff
const double h = 0.02;
const int FFTWGrain = 39;//78; //Determines Q = 2*Pi*FFTWGrain/(h*ARRAY)
const bool LKApprox = false;

int NucMethod = 1; //0 = IWS, 1 = Nucleon Sampler, 2 = flat


/************************/
/****(0) Woods-Saxon*****/
/************************/
#define WSa .6
#define WSr 6. //make sure units are fm
#define NuclDens 0.1693 //Nucleon Density at center of WS

/************************/
/***(1) Nuclei Sampler***/
/************************/

/*#define N1 197 //Nuclei Species 1 - Au
#define N2 197 //Nuclei Species 2 -Au
*/
#define N1 208 //Nuclei Species 1 - Pb
#define N2 208 //Nuclei Species 2 -Pb

const double impact = 6.;
double r_core=0.5; //radius of Nucleon for collisions
double gauss_n=0.5; //radius of gaussian profile of mu per nucleon
int Atom=1; // 0 == Au-Au ; 1 == Pb-Pb
double NucScale=39.3;//87;//400.;//8.//170.//120;//25.2;//18.;//1400.; //Scale of mu- probes different collision energies


/************************/
/******(2) Flat mu*******/
/************************/
const double mu = 60.; // Used if NucMethod 2
