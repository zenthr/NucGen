/****************************/
/****************************/
/*****Temporary Warnings*****/
/*    No Coarse Grain       */
/*     Larger Grid          */
/*     These were 80/841    */
/* Recall variable in DeclT */
/****************************/
/****************************/



/*ingw*/
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sstream>
#include <string>
#include <cstring>
#include <float.h>
#include <vector>
#include <time.h>

#include <stdio.h>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include <typeinfo> //Not sure what this does

#include "CGCEDefsHD-FFTW.h"
#include "ran_gen.cpp"
#define BIG  1e30

using namespace std;

    int flr;
    int clg;
    static double timerAlpha;
    static double timerA;
    static double timerRec;
    static double timer;
    double runtime;
    int GridRat; //Ratio of Coarse to Fine Grid
    int imp = impact/(2*h);
    double pedestal = 0;//pow(10.,-16.);

    double** Nucl1;//Nucleus 1
    double** Nucl2;//Nucleus 2
    
    double***** RhoCG;

    int size = 200;
    int Target = 1;
    static double CDFTab[NUMSTEP][2];
    static double MuShape[2][ARRAY][ARRAY];
//    static double*** MuShape;
//    static double Rho[2][NG][ARRAY][ARRAY]; //Fine charge density, only use for testing
//    static double RhoCor[2][ARRAY][ARRAY];
    double Hist[16][200];
    double Ecc[2][200];
    double EBHist[4][200];
    double CSCharge[200][3];
    double EBsq[100][4];
    double RhoStep;
    double***** CovPot;
//    static double CovPotChk[ARRAY][ARRAY];
    double** Mold;
    
    // Statistics for debugging
//    static double RhoCGCor[2][CGRho][CGRho];
//    static double RhoCGAv[2][CGRho][CGRho]; //0 not squar, 1 square
//    static double RhoCGVar[10][CGRho][CGRho];
      static double*** AlphaCor;
      static double** T00Av;
//    static double AlphaAv[2][ARRAY][ARRAY];
      static complex<double>*** AxCor;
      static complex<double>*** AyCor;
      static double**** MuShapeCenter;
//    static complex<double> AxAv[2][ARRAY][ARRAY];
//    static complex<double> AyAv[2][ARRAY][ARRAY];
//    static complex<double> E0Cor[2][ARRAY][ARRAY];
//    static complex<double> E0Av[2][ARRAY][ARRAY];
//    static complex<double> B0Cor[2][ARRAY][ARRAY];
//    static complex<double> B0Av[2][ARRAY][ARRAY];
 //   static complex<double> EDCor[2][ARRAY][ARRAY];
//    static complex<double> EDAv[2][ARRAY][ARRAY];
//    static complex<double> ED2Av[2][ARRAY][ARRAY];
//    static complex<double> PL2Av[2][ARRAY][ARRAY];
//    static complex<double> FAlphAv[2][ARRAY][ARRAY];
//    static complex<double> FDirAv[2][ARRAY][ARRAY];
//    static complex<double> DivAAv[ARRAY][ARRAY];
//    static complex<double> DivBAv[ARRAY][ARRAY];
//    static complex<double> deltav[ARRAY][ARRAY];
   // static double CovPotAv[ARRAY][ARRAY];
   // static double CovPotVar[ARRAY][ARRAY];
    
   // static double PotCorAvCnt[ARRAY][ARRAY];
   // static double PotCorVarCnt[ARRAY][ARRAY];

  //  static double PotCorAvOff[ARRAY][ARRAY];
  //  static double PotCorVarOff[ARRAY][ARRAY];
       static int Numb;
    complex<double>***** U;
    
    complex<double>****** A; //Fawkes(Fox?)-Schwinger Potential [Nuclei][Transverse Coord]@[x][y][row[pcolumn]
    complex<double>****** Ad;

    complex<double>***** A2;
    complex<double>****** Ai;
    complex<double>******* Fij;
    static double*** T00;
    static double*** T11;
    static double*** T22;
    static double*** T33;
    static double*** T01;
    static double*** T02;
    static double*** T03;
    static double*** T03FT;
    static double*** T12;
    static double*** T13;
    static double*** T23;
    static double**** Recall;
    static double*** EBGrid;
//    static complex<double>*** AScale;

    complex<double>**** AT;
    complex<double>**** AS;

    complex<double>***** AvU;
    
    static complex<double> Bx[3][3][3];
    
    double CollDens;
    double PartDens;
    double NColl;
    double NPart;
    int NCollGeo;
    int NPartGeo;
    ofstream GlauberGeo;

std::mt19937_64 eng(std::random_device{}());
//std::mt19937_64 eng(145709713); // run with set seed: 145709713
double ran() {std::uniform_real_distribution<double> uniran(0.0,1.0); return uniran(eng);}

   // static complex<double> AB[2][MinArray][MinArray][3][3];
   // static complex<double> A2nd[2][MinArray][MinArray][3][3];
   // static complex<double> A4th[2][MinArray][MinArray][3][3];
    //static complex<double> E0[ARRAY][ARRAY][3][3];
    //static complex<double> B0[ARRAY][ARRAY][3][3];
    //static complex<double> edens[ARRAY][ARRAY];
    //static complex<double> edens2[2][ARRAY][ARRAY];
    //static complex<double> PL2[2][ARRAY][ARRAY];
    //static complex<double> Flow[2][2][ARRAY][ARRAY];
    //static complex<double> delt[ARRAY][ARRAY];
    
// http://jean-pierre.moreau.pagesperso-orange.fr/c_bessel.html
// From TBESSK


//  Bessel Function of the 1st kind of order zero
/*
    double BESSI0(double X)  {
      double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
      P1=1.0; P2=3.5156229; P3=3.0899424; P4=1.2067429;
      P5=0.2659732; P6=0.360768e-1; P7=0.45813e-2;
      Q1=0.39894228; Q2=0.1328592e-1; Q3=0.225319e-2; Q4=-0.157565e-2;
      Q5=0.916281e-2; Q6=-0.2057706e-1; Q7=0.2635537e-1;
      Q8=-0.1647633e-1; Q9=0.392377e-2;
      if (fabs(X) < 3.75) {
        Y=pow((X/3.75),2);
        return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
      }
      else {
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return AX*BX;
      }
    }
*/

// Matthew things
typedef std::complex<double> cplx;
/*
cplx stat_alpha[2][ARRAY][ARRAY]{};
cplx stat_beta[2][ARRAY][ARRAY]{};
cplx stat_T[Order/2][ARRAY][ARRAY][4][4]{};
*/
 cplx I{0.0, 1.0};


 void clear_simple(cplx array[ARRAY][ARRAY])
 {
	 for (int x = 0; x < ARRAY; x++)
	 {
		 for (int y = 0; y < ARRAY; y++)
		 {
			 array[x][y] = cplx{ 0.0, 0.0 };
		 }

	 }

 }
 void clear(cplx array[ARRAY][ARRAY][3][3])
 {
	 for (int x = 0; x < ARRAY; x++)
	 {
		 for (int y = 0; y < ARRAY; y++)
		 {
			 for (int a = 0; a < 3; a++)
			 {
				 for (int b = 0; b < 3; b++)
					 array[x][y][a][b] = cplx{ 0.0, 0.0 };
			 }
		 }
	 }
 }


 void reg_dirX(int dir_number, cplx arr_1[ARRAY][ARRAY], cplx out_arr[ARRAY][ARRAY])
 {
	 clear_simple(out_arr);
	 int x = 2;
	 for (int matt_index = 0; matt_index < ARRAY; matt_index++)
	 {
		x = x% ARRAY;
		 if (x >= 2 * dir_number && x < (ARRAY - 2 * dir_number))
		 {
			 for (int y = 0; y < ARRAY; y++)
			 {
						 cplx di{};
						 di = 8.0*(arr_1[x+1][y] - arr_1[x-1][y])
							-1.0*(arr_1[x+2][y] - arr_1[x-2][y]);

						 di = di / (12.0*h);
						 out_arr[x][y] = di;
			 }

		 }
		 if(x<2)
		 {
			 for (int y = 0; y < ARRAY; y++)
			 {
						cplx AVG = out_arr[2][y];
						
						out_arr[x][y] = AVG;
			 }
		 }
		 if(x >=  ARRAY - 2)
		 {
			 for (int y = 0; y < ARRAY; y++)	 
				 {
						cplx AVG = out_arr[98][y];
						
						out_arr[x][y] = AVG;			 
		
				 }
		
		 }
		 x++;
	 }

 }

 void reg_dirY(int dir_number, cplx arr_1[ARRAY][ARRAY], cplx out_arr[ARRAY][ARRAY])
 {
	 clear_simple(out_arr);
	 int y = 2;
	 for (int matt_index = 0; matt_index < ARRAY;  matt_index++)
	 {
		 y  = y% ARRAY;

		 if (y >= 2 * dir_number && y < (ARRAY - 2 * dir_number))
		 {
			 for (int x = 0; x < ARRAY; x++)
			 {
				 cplx di{};
				di = 8.0*(arr_1[x][y+1] - arr_1[x][y-1])
				-1.0*(arr_1[x][y+2] - arr_1[x][y-2]);
				
				 di = di / (12.0*h);
				 out_arr[x][y] = di;
			 }

		 }
		 if(y<2)
		 {
			 for (int x = 0; x < ARRAY; x++)
			 {
						cplx AVG = out_arr[x][2];
						
						out_arr[x][y] = AVG ;
			 }
		 }
		 if(y >=  ARRAY - 2)
		 {
			 for (int x = 0; x < ARRAY; x++)
				 
				 {
						cplx AVG = out_arr[x][98];
						
						out_arr[x][y] = AVG;			 
		

				 }
		 }
		 y++;
	 }

 }
 
void arr_add(int option,cplx arr_1[ARRAY][ARRAY][3][3], cplx arr_2 [ARRAY][ARRAY][3][3], cplx out_arr[ARRAY][ARRAY][3][3])
{
	clear(out_arr);
	for (int x =0; x < ARRAY; x++)
		{
		for (int y = 0; y< ARRAY; y++ )
			{
			for (int a = 0; a <= 2; a++)
				{
					for (int b = 0; b <= 2; b++)
					{ 	if (option == 0 )
						{
						out_arr[x][y][a][b]= arr_1[x][y][a][b] + arr_2[x][y][a][b];
						}
						
						else
						{	
						out_arr[x][y][a][b]= arr_1[x][y][a][b] - arr_2[x][y][a][b];
						}
					
					
					}
				}
			}
		}
	}
	

void arr_mult(cplx arr_1[ARRAY][ARRAY][3][3], cplx arr_2 [ARRAY][ARRAY][3][3], cplx out_arr[ARRAY][ARRAY][3][3])
{
	clear(out_arr);
	for (int x =0; x < ARRAY; x++)
	{
		for (int y = 0; y< ARRAY; y++ )
			{
			for (int a = 0; a <= 2; a++)
				{
					for (int b = 0; b <= 2; b++)
					{
					std::complex<double> sum = 0;
					for (int c = 0; c <= 2; c++)
					{
						sum += arr_1[x][y][a][c]* arr_2[x][y][c][b];
					}

					out_arr[x][y][a][b] = sum;
				}
			}
		}
	}
	
}

void arr_comm(cplx arr_1[ARRAY][ARRAY][3][3], cplx arr_2[ARRAY][ARRAY][3][3], cplx out_arr[ARRAY][ARRAY][3][3])
{
	clear(out_arr);

	for (int x = 0; x < ARRAY; x++)
	{
		for (int y = 0; y < ARRAY; y++)
		{
			for (int a = 0; a <= 2; a++)
			{
				for (int b = 0; b <= 2; b++)
				{
					std::complex<double> sum = 0;
					for (int c = 0; c <= 2; c++)
					{
						sum += arr_1[x][y][a][c] * arr_2[x][y][c][b]
							- arr_1[x][y][c][b] * arr_2[x][y][a][c];

					}

					out_arr[x][y][a][b] = sum;
				}
			}
		}
	}
}


cplx trace(cplx arr_1[3][3])
{
	cplx sum = 0;
	for (int a = 0; a < 3; a++)
	{
		sum += arr_1[a][a];
	}
	return .5*sum;
}


void DirX(int dir_number, cplx arr_1[ARRAY][ARRAY][3][3], cplx arr_2[2][ARRAY][ARRAY][3][3], cplx out_arr[ARRAY][ARRAY][3][3])
{
	clear(out_arr);
	int x = 2;
	for (int matt_index = 0; matt_index < ARRAY; matt_index++)
	{
		x = x% ARRAY;
		if (x >= 2*dir_number && x < ARRAY-2*dir_number )
		{
			for (int y = 0; y < ARRAY; y++)
			{
				cplx sum[3][3]{};
				for (int a = 0; a <= 2; a++)
					{
						for (int b = 0; b <= 2; b++)
						{	
							sum[a][b]= cplx{0.0, 0.0};
							for (int c = 0; c <= 2; c++)
							{
								sum[a][b] += arr_1[x][y][a][c] * arr_2[0][x][y][c][b]
								- arr_1[x][y][c][b] * arr_2[0][x][y][a][c];

							}
						}
					}
				for (int a = 0; a < 3; a++)
				{
					for (int b = 0; b < 3; b++)
					{
						cplx di{};
						 di = 8.0*(arr_1[x+1][y][a][b] - arr_1[x-1][y][a][b])
							-1.0*(arr_1[x+2][y][a][b] - arr_1[x-2][y][a][b]);

						 di = di / (12.0*h) - I*g*sum[a][b];
						out_arr[x][y][a][b] = -di;
					}
				}
			}

		}
		if(x < 2)
		{
			for (int y = 0; y < ARRAY; y++)
			{
				for (int a = 0; a < 3; a++)
				{
					for (int b = 0; b < 3; b++)
					{
						cplx AVG = out_arr[2][y][a][b];
						
						out_arr[x][y][a][b] = AVG;
					}	
				}
				
			}

		 }
		if(x >= ARRAY -2)
		{
			for (int y = 0; y < ARRAY; y++)
			{
				for (int a = 0; a < 3; a++)
				{
					for (int b = 0; b < 3; b++)
					{
						cplx AVG = out_arr[98][y][a][b];
						
						out_arr[x][y][a][b] = AVG;
					}
				}
			}

		 }
		 x++;

		}
	}





void DirY(int dir_number, cplx arr_1[ARRAY][ARRAY][3][3], cplx arr_2[2][ARRAY][ARRAY][3][3], cplx out_arr[ARRAY][ARRAY][3][3])
{
	clear(out_arr);
	int y = 2;
	for (int matt_index = 0; matt_index < ARRAY; matt_index++)
	{
		y = y%ARRAY;
		if (y > 2 * dir_number -1 && y < ARRAY - 2 * dir_number)
		{
			for (int x = 0; x < ARRAY; x++)
			{
				cplx sum[3][3]{};
				for (int a = 0; a <= 2; a++)
				{
					for (int b = 0; b <= 2; b++)
					{
						sum[a][b] = cplx{0.0, 0.0};
						for (int c = 0; c <= 2; c++)
						{
							sum[a][b] += arr_1[x][y][a][c] * arr_2[1][x][y][c][b]
								- arr_1[x][y][c][b] * arr_2[1][x][y][a][c];

						}
					}
				}
				for (int a = 0; a < 3; a++)
				{
					for (int b = 0; b < 3; b++)
					{
						cplx di{};
						 di = 8.0*(arr_1[x][y+1][a][b] - arr_1[x][y-1][a][b])
							-1.0*(arr_1[x][y+2][a][b] - arr_1[x][y-2][a][b]);

						 di = di / (12.0*h) - I*g*sum[a][b];
						out_arr[x][y][a][b] = -di;
					}
				}
			}

		}
		 if(y<2)
		 {
			 for (int x = 0; x < ARRAY; x++)
			 {
				for (int a = 0; a < 3; a++)
				{
					for (int b = 0; b < 3; b++)
					{
						cplx AVG = out_arr[x][2][a][b];
						
						out_arr[x][y][a][b] = AVG;
					}
				}
				
			 }
		 }
		 if(y >= ARRAY -2)
		 {
			 for (int x = 0; x < ARRAY; x++)
			 {
				for (int a = 0; a < 3; a++)
				{
					for (int b = 0; b < 3; b++)
					{				 
						cplx AVG = out_arr[x][98][a][b];
		
						out_arr[x][y][a][b] = AVG;			 

					}	
				}
			
			}

		
		 }
		 
	y++;
	}

}



void print_array(std::complex<double> arr[ARRAY][ARRAY][3][3], int x, int y, std::string title)
{
	std::cout << title << std::endl;
	for (int m = 0; m < 3; m++)
	{
		for (int n = 0; n < 3; n++)
		{

			std::cout << arr[x][y][m][n]<< "  ";

		}
		std::cout << '\n';
	}
	std::cout << "\n \n";
}

void print_array_t(cplx arr[ARRAY][ARRAY][4][4], int x, int y, std::string title)
{
	std::cout << title << std::endl;
	for (int m = 0; m < 4; m++)
	{
		for (int n = 0; n < 4; n++)
		{

			std::cout << arr[x][y][m][n]<< "  ";

		}
		std::cout << '\n';
	}
	std::cout << "\n \n";
}

void real_print_array(std::complex<double> arr[ARRAY][ARRAY][3][3], int x, int y, std::string title)
{
	std::cout << title << std::endl;
	for (int m = 0; m < 3; m++)
	{
		for (int n = 0; n < 3; n++)
		{

			std::cout << real(arr[x][y][m][n])<< "  ";

		}
		std::cout << '\n';
	}
	std::cout << "\n \n";
}

void real_print_array_t(cplx arr[ARRAY][ARRAY][4][4], int x, int y, std::string title)
{
	std::cout << title << std::endl;
	for (int m = 0; m < 4; m++)
	{
		for (int n = 0; n < 4; n++)
		{

			std::cout << real(arr[x][y][m][n])<< "  ";

		}
		std::cout << '\n';
	}
	std::cout << "\n \n";
}


void array_check(int buffer, cplx arr_1[ARRAY][ARRAY][4][4], cplx arr_2[ARRAY][ARRAY][4][4])
{	double average = 0;
	int count = 0 ;
	for (int x =buffer; x < ARRAY-buffer; x++)
	{
		for (int y =buffer; y < ARRAY-buffer; y++)
		{
			for (int m = 0; m < 4; m++)
			{
				for (int n = m; n < 4; n++)
				{
					double comp = std::abs((real(arr_1[x][y][m][n]-arr_2[x][y][m][n]))/ (.5*(real(arr_1[x][y][m][n]+ arr_2[x][y][m][n]))+.000000001));
					average += comp;
					count++;
					if(.5 <= comp )
					{
						double wtf = std::abs((.5*(real(arr_1[x][y][m][n]+ arr_2[x][y][m][n]))+.000000001));
						if (1000. <= wtf)
						{
						std::cout <<"percent error  " << comp <<"  array x value " << x << "  array y value " << y << "  \n"; 
						real_print_array_t(arr_1, x, y, "Matt");
						real_print_array_t(arr_2, x, y, "Steven");
						//return;
						}
						
						else
						{}
					
					}
				else
				{}
			
			
				}
				
			}
			
		}
	}
	std::cout << "all good, Average difference is " << average/count << std::endl;
}

 double return_val(double array[Order/2][ARRAY][ARRAY][4][4][100], int order, int x, int y, int m, int n)
	 {
		 for (int i =0; i < 100; i++)
		 {
			 if ( .5*std::abs(((array[order][x][y][m][n][i] - array[order -1][x][y][m][n][i])/((array[order][x][y][m][n][i] + array[order -1][x][y][m][n][i])))) > .2)
			 {
				 return i/100.0;
			 }
		 else
		 {}	 
		 }
	 }

//Matthew things

    /**********/
    /*Beg Decl*/
    /**********/
    
    /*FFFFFFFFFFFUUUUUUUUUUUUU*/
   void DeclDelMe()
   {
        
      MuShapeCenter = new double***[2];
      for(int nuc=0; nuc<2; nuc++)
      {
         MuShapeCenter[nuc] = new double**[2];
         for(int kind=0; kind<2; kind++)
         {
            MuShapeCenter[nuc][kind] = new double*[ARRAY];
            for(int x=0; x<ARRAY; x++)
            {
               MuShapeCenter[nuc][kind][x] = new double[ARRAY];
               for(int y=0; y<ARRAY; y++)
               {
                  MuShapeCenter[nuc][kind][x][y] = 0;
               }
            }
         }
      }
   }
    
    
    
    
   //Definition for static double MuShape[2][CGRho][CGRho];
   //Hold the shape of Mu
/*
   void DeclMuShape()
   {
      MuShape = new double**[2];
      for(int nuc=0; nuc<2; nuc++)
      {
         MuShape[nuc] = new double*[CGRho];
         for(int x=0; x<ARRAY; x++)
         {
            MuShape[nuc][x] = new double[CGRho];
         }
      }
   }
   
   void RelMuShape()
   {
      for(int nuc=0; nuc<2; nuc++)
      {
         for(int x=0; x<ARRAY; x++)
         {
            delete [] MuShape[nuc][x];
         }
         delete [] CovPot[nuc];
      }
      delete [] CovPot;
   }
*/
    
   //Definition for static double CovPot[2][NG][etamax][ARRAY][ARRAY]
   //Covariant potential alpha by color component
   //[2] for Nuclei
   void DeclCovPot()
   {
      CovPot = new double****[2];
      for(int nuc=0; nuc<2; nuc++)
      {
         CovPot[nuc] = new double***[NG];
         for(int col=0; col<NG; col++)
         {
            CovPot[nuc][col] = new double**[etamax];
            for(int eta=0; eta<etamax; eta++)
            {
               CovPot[nuc][col][eta] = new double*[ARRAY];
               for(int x=0; x<ARRAY; x++)
               {
                  CovPot[nuc][col][eta][x] = new double[ARRAY];
                  for(int y=0; y<ARRAY; y++)
                  {
                     CovPot[nuc][col][eta][x][y] = 0;
                  }
               }
            }
         }
      }
   }
   
   void RelCovPot()
   {
      for(int nuc=0; nuc<2; nuc++)
      {
         for(int col=0; col<NG; col++)
         {
            for(int eta=0; eta<etamax; eta++)
            {
               for(int x=0; x<ARRAY; x++)
               {
                  delete [] CovPot[nuc][col][eta][x];
               }
               delete [] CovPot[nuc][col][eta];
            }
            delete [] CovPot[nuc][col];
         }
         delete [] CovPot[nuc];
      }
      delete [] CovPot;
   }
   
   //Definition for static double Mold[MOLDARRAY][MOLDARRAY]
   void DeclMold()
   {
      Mold = new double*[MOLDARRAY];
      for(int x=0; x<MOLDARRAY; x++)
      {
         Mold[x] = new double[MOLDARRAY];
         for(int y=0; y<MOLDARRAY; y++)
         {
            Mold[x][y] = 0;
         }
      }
   }
   
   void RelMold()
   {
      for(int x=0; x<MOLDARRAY; x++)
      {
         delete [] Mold[x];
      }
      delete [] Mold;
   }
   
   //Definition for complex<double> U[2][ARRAY][ARRAY][3][3]
   void DeclU()
   {
      U = new complex<double>****[2];
      for(int nuc=0; nuc<2; nuc++)
      {
         U[nuc] = new complex<double>***[ARRAY];
         for(int x=0; x<ARRAY; x++)
         {
            U[nuc][x] = new complex<double>**[ARRAY];
            for(int y=0; y<ARRAY; y++)
            {
               U[nuc][x][y] = new complex<double>*[3];
               for(int a=0; a<3; a++)
               {
                  U[nuc][x][y][a] = new complex<double>[3];
                  for(int b=0; b<3; b++)
                  {
                     U[nuc][x][y][a][b] = 0;
                  }
               }
            }
         }
      }
   }
   
   void RelU()
   {
      for(int nuc=0; nuc<2; nuc++)
      {
         for(int x=0; x<ARRAY; x++)
         {
            for(int y=0; y<ARRAY; y++)
            {
               for(int a=0; a<3; a++)
               {
                  delete [] U[nuc][x][y][a];
               }
               delete [] U[nuc][x][y];
            }
            delete [] U[nuc][x];
         }
         delete [] U[nuc];
      }
      delete [] U;
   }

   //Definition for static complex<double> A[2][2][ARRAY][ARRAY][3][3]
   //A perp for individual nuclei
   void DeclA()
   {
      A= new complex<double>*****[2];
      for(int i=0; i<2; i++)
      {
         A[i]= new complex<double>****[2];
         for(int j=0; j<2; j++)
         {
            A[i][j]= new complex<double>***[ARRAY];
            for(int ax=0; ax<ARRAY; ax++)
            {
               A[i][j][ax]= new complex<double>**[ARRAY];
               for(int ay=0; ay<ARRAY; ay++)
               {
                  A[i][j][ax][ay]= new complex<double>*[3];
                  for(int a=0; a<3; a++)
                  {
                     A[i][j][ax][ay][a]= new complex<double>[3];
                     for(int b=0; b<3; b++)
                     {
                        A[i][j][ax][ay][a][b] = 0;
                     }
                  }
               }
            }
         }
      }
      /*
      AScale = new complex<double>**[NG];
      for(int col=0; col<NG; col++)
      {
         AScale[col] = new complex<double>*[ARRAY];
         for(int ax=0; ax<ARRAY; ax++)
         {
            AScale[col][ax] = new complex<double>[ARRAY];
            for(int bx=0; bx<ARRAY; bx++)
            {
               AScale[col][ax][bx] = 0;
            }
         }
      }
      */
   }

   void RelA()
   {
      for(int i=0; i<2; i++)
      {
         for(int j=0; j<2; j++)
         {
            for(int ax=0; ax<ARRAY; ax++)
            {
               for(int ay=0; ay<ARRAY; ay++)
               {
                  for(int a=0; a<3; a++)
                  {
                     delete [] A[i][j][ax][ay][a];
                  }
                  delete [] A[i][j][ax][ay];
               }
               delete [] A[i][j][ax];
            }
            delete [] A[i][j];
         }
         delete [] A[i];
      }
      delete [] A;
   }
   
   //Definition for static complex<double> A2[Order/2 +1][ARRAY][ARRAY][3][3]
   //A longitudinal for combined field
   void DeclA2()
   {
      A2 = new complex<double>****[(Order/2)+1];
      for(int BigO=0; BigO<=(Order/2); BigO++)
      {
         A2[BigO] = new complex<double>***[ARRAY];
         for(int x=0; x<ARRAY; x++)
         {
            A2[BigO][x] = new complex<double>**[ARRAY];
            for(int y=0; y<ARRAY; y++)
            {
               A2[BigO][x][y] = new complex<double>*[3];
               for(int a=0; a<3; a++)
               {
                  A2[BigO][x][y][a] = new complex<double>[3];
                  for(int b=0; b<3; b++)
                  {
                     A2[BigO][x][y][a][b] = 0;
                  }                  
               }
            }
         }
      }
   }
   
   void RelA2()
   {
      for(int BigO=0; BigO<=(Order/2); BigO++)
      {
         for(int x=0; x<ARRAY; x++)
         {
            for(int y=0; y<ARRAY; y++)
            {
               for(int a=0; a<3; a++)
               {
                  delete [] A2[BigO][x][y][a];
               }
               delete [] A2[BigO][x][y];
            }
            delete [] A2[BigO][x];
         }
         delete [] A2[BigO];
      }
      delete [] A2;
   }
   
   //Definition for static complex<double> Ai[Order/2][2][ARRAY][ARRAY][3][3]
   //A perp for combined field
   void DeclAi()
   {
      Ai = new complex<double>*****[(Order/2)+1];
      
      for(int BigO=0; BigO<=(Order/2); BigO++)
      {

         Ai[BigO] = new complex<double>****[2];
         for(int i=0; i<2; i++)
         {
            Ai[BigO][i] = new complex<double>***[ARRAY];
            for(int x=0; x<ARRAY; x++)
            {
               Ai[BigO][i][x] = new complex<double>**[ARRAY];
               for(int y=0; y<ARRAY; y++)
               {
                  Ai[BigO][i][x][y] = new complex<double>*[3];
                  for(int a=0; a<3; a++)
                  {
                     Ai[BigO][i][x][y][a] = new complex<double>[3];
                     for(int b=0; b<3; b++)
                     {
                        Ai[BigO][i][x][y][a][b] = 0;
                     }
                  }
               }
            }
         }
      }
   }

   void RelAi()
   {
      for(int BigO=0; BigO<=(Order/2); BigO++)
      {
         for(int i=0; i<2; i++)
         {
            for(int x=0; x<ARRAY; x++)
            {
               for(int y=0; y<ARRAY; y++)
               {
                  for(int a=0; a<3; a++)
                  {
                     delete [] Ai[BigO][i][x][y][a];
                  }
                  delete [] Ai[BigO][i][x][y];
               }
               delete [] Ai[BigO][i][x];
            }
            delete [] Ai[BigO][i];
         }
         delete [] Ai[BigO];
      }
      delete [] Ai;
   }

   //Definition for static complex<double> Fij[Order/2][2][3][ARRAY][ARRAY][3][3]
   //E and B information
   //[2] 0 is E; 1 is B
   //[3] 0 is Long.; 1/2 is x/y
   //For long., order is 2*Order; for perp., order is 2*order+1
   void DeclFij()
   {
      Fij = new complex<double>******[(Order/2)+1];
      for(int BigO=0; BigO<=(Order/2); BigO++)
      {
         Fij[BigO] = new complex<double>*****[2];
         for(int em=0; em<2; em++)
         {
            Fij[BigO][em] = new complex<double>****[3];
            for(int i=0; i<3; i++)
            {
               Fij[BigO][em][i] = new complex<double>***[ARRAY];
               for(int x=0; x<ARRAY; x++)
               {
                  Fij[BigO][em][i][x] = new complex<double>**[ARRAY];
                  for(int y=0; y<ARRAY; y++)
                  {
                     Fij[BigO][em][i][x][y] = new complex<double>*[3];
                     for(int a=0; a<3; a++)
                     {
                        Fij[BigO][em][i][x][y][a] = new complex<double>[3];
                        for(int b=0; b<3; b++)
                        {
                            Fij[BigO][em][i][x][y][a][b] = 0;
                        }
                     }
                  }
               }
            }
         }
      }
   }

   void RelFij()
   {
      for(int BigO=0; BigO<=(Order/2); BigO++)
      {
         for(int em=0; em<2; em++)
         {
            for(int i=0; i<3; i++)
            {
               for(int x=0; x<ARRAY; x++)
               {
                  for(int y=0; y<ARRAY; y++)
                  {
                     for(int a=0; a<3; a++)
                     {
                        delete [] Fij[BigO][em][i][x][y][a];
                     }
                     delete [] Fij[BigO][em][i][x][y];
                  }
                  delete [] Fij[BigO][em][i][x];
               }
               delete [] Fij[BigO][em][i];
            }
            delete [] Fij[BigO][em];
         }
         delete [] Fij[BigO];
      }
      delete [] Fij;
   }

   //Definition for complex<double> AT/AS[ARRAY][ARRAY][3][3]
   //Dumb dummy variables, dummy
   void DeclDum()
   {
      AT = new complex<double>***[ARRAY];
      AS = new complex<double>***[ARRAY];
      for(int x=0; x<ARRAY; x++)
      {
         AT[x] = new complex<double>**[ARRAY];
         AS[x] = new complex<double>**[ARRAY];
         for(int y=0; y<ARRAY; y++)
         {
            AT[x][y] = new complex<double>*[3];
            AS[x][y] = new complex<double>*[3];
            for(int a=0; a<3; a++)
            {
               AT[x][y][a] = new complex<double>[3];
               AS[x][y][a] = new complex<double>[3];
               for(int b=0; b<3; b++)
               {
                  AT[x][y][a][b] = 0;
                  AS[x][y][a][b] = 0;
               }
            }
         }
      }
   }

   void RelDum()
   {
      for(int x=0; x<ARRAY; x++)
      {
         for(int y=0; y<ARRAY; y++)
         {
            for(int a=0; a<3; a++)
            {
               delete [] AT[x][y][a];
               delete [] AS[x][y][a];
            }
            delete [] AT[x][y];
            delete [] AS[x][y];
         }
         delete [] AT[x];
         delete [] AS[x];
      }
      delete [] AT;
      delete [] AS;
   }

   //Definition for complex<double> AvU[2][ARRAY][ARRAY][3][3]
   //Dumb dummy variables, dummy
   void DeclAvU()
   {
      AvU = new complex<double>****[2];
      for(int nuc=0; nuc<2; nuc++)
      {
         AvU[nuc] = new complex<double>***[ARRAY];
         for(int x=0; x<ARRAY; x++)
         {
            AvU[nuc][x] = new complex<double>**[ARRAY];
            for(int y=0; y<ARRAY; y++)
            {
               AvU[nuc][x][y] = new complex<double>*[3];
               for(int a=0; a<3; a++)
               {
                  AvU[nuc][x][y][a] = new complex<double>[3];
                  for(int b=0; b<3; b++)
                  {
                     AvU[nuc][x][y][a][b] = 0;
                  }
               }
            }
         }
      }
   }

   void RelAvU()
   {
      for(int nuc=0; nuc<2; nuc++)
      {
         for(int x=0; x<ARRAY; x++)
         {
            for(int y=0; y<ARRAY; y++)
            {
               for(int a=0; a<3; a++)
               {
                  delete [] AvU[nuc][x][y][a];
               }
               delete [] AvU[nuc][x][y];
            }
            delete [] AvU[nuc][x];
          }
         delete [] AvU[nuc];
      }
      delete [] AvU;
   }

   //Definition for static double T##[Order/2][ARRAY][ARRAY]
   //Stress-Energy Tensor
   //Need a better storage solution
   //Odds and evens are all sorts of f***ed up
   void DeclT()
   {
      T00 = new double**[(Order/2)+1];
      T11 = new double**[(Order/2)+1];
      T22 = new double**[(Order/2)+1];
      T33 = new double**[(Order/2)+1];
      T01 = new double**[(Order/2)+1];
      T02 = new double**[(Order/2)+1];
      T03 = new double**[(Order/2)+1];
      T12 = new double**[(Order/2)+1];
      T13 = new double**[(Order/2)+1];
      T23 = new double**[(Order/2)+1];
      
      for(int BigO=0; BigO<=(Order/2); BigO++)
      {
         T00[BigO] = new double*[ARRAY];
         T11[BigO] = new double*[ARRAY];
         T22[BigO] = new double*[ARRAY];
         T33[BigO] = new double*[ARRAY];
         T01[BigO] = new double*[ARRAY];
         T02[BigO] = new double*[ARRAY];
         T03[BigO] = new double*[ARRAY];
         T12[BigO] = new double*[ARRAY];
         T13[BigO] = new double*[ARRAY];
         T23[BigO] = new double*[ARRAY];
         
         for(int x=0; x<ARRAY; x++)
         {
            T00[BigO][x] = new double[ARRAY];
            T11[BigO][x] = new double[ARRAY];
            T22[BigO][x] = new double[ARRAY];
            T33[BigO][x] = new double[ARRAY];
            T01[BigO][x] = new double[ARRAY];
            T02[BigO][x] = new double[ARRAY];
            T03[BigO][x] = new double[ARRAY];
            T12[BigO][x] = new double[ARRAY];
            T13[BigO][x] = new double[ARRAY];
            T23[BigO][x] = new double[ARRAY];
            for(int y=0; y<ARRAY; y++)
            {
               T00[BigO][x][y] = 0;
               T11[BigO][x][y] = 0;
               T22[BigO][x][y] = 0;
               T33[BigO][x][y] = 0;
               T01[BigO][x][y] = 0;
               T02[BigO][x][y] = 0;
               T03[BigO][x][y] = 0;
               T12[BigO][x][y] = 0;
               T13[BigO][x][y] = 0;
               T23[BigO][x][y] = 0;
            }
         }
      }

      T03FT = new double**[2];

      for(int BigO=0; BigO<=1; BigO++)
      {
         T03FT[BigO] = new double*[(2*size+1)];

         for(int x=0; x<(2*size+1); x++)
         {
            T03FT[BigO][x] = new double[(2*size+1)];

            for(int y=0; y<(2*size+1); y++)
            {
               T03FT[BigO][x][y] = 0;
            }
         }
      }

      Recall = new double*** [(Order/2)+1];
      for(int BigO=0; BigO<=(Order/2); BigO++)
      {
         Recall[BigO] = new double**[10];
         for(int a=0; a<10; a++)
         {
            Recall[BigO][a] = new double*[ARRAY];
            for(int x=0; x<ARRAY; x++)
            {
               Recall[BigO][a][x] = new double[ARRAY];
               for(int y=0; y<ARRAY; y++)
               {
                  Recall[BigO][a][x][y] = 0;
               }               
            }
         }
      }

      EBGrid = new double** [(Order/2)+1];
      for(int BigO=0; BigO<=(Order/2); BigO++)
      {
         EBGrid[BigO] = new double*[ARRAY];
         for(int x=0; x<ARRAY; x++)
         {
            EBGrid[BigO][x] = new double[ARRAY];
            for(int y=0; y<ARRAY; y++)
            {
               EBGrid[BigO][x][y] = 0;
            }
         }
      }
   }

   void RelT()
   {
      for(int BigO=0; BigO<=(Order/2); BigO++)
      {
         for(int x=0; x<ARRAY; x++)
         {
            delete [] T00[BigO][x];
            delete [] T11[BigO][x];
            delete [] T22[BigO][x];
            delete [] T33[BigO][x];
            delete [] T01[BigO][x];
            delete [] T02[BigO][x];
            delete [] T03[BigO][x];
            delete [] T12[BigO][x];
            delete [] T13[BigO][x];
            delete [] T23[BigO][x];
         }
         delete [] T00[BigO];
         delete [] T11[BigO];
         delete [] T22[BigO];
         delete [] T33[BigO];
         delete [] T01[BigO];
         delete [] T02[BigO];
         delete [] T03[BigO];
         delete [] T12[BigO];
         delete [] T13[BigO];
         delete [] T23[BigO];
      }
      delete [] T00;
      delete [] T11;
      delete [] T22;
      delete [] T33;
      delete [] T01;
      delete [] T02;
      delete [] T03;
      delete [] T12;
      delete [] T13;
      delete [] T23;
   }

   //Definition for static complex<double> AxCor[2][ARRAY][ARRAY]
   //<AA> correlation
   
   void DeclAACor()
   {
      AxCor = new complex<double>**[2];
      AyCor = new complex<double>**[2];
      for(int nuc=0; nuc<2; nuc++)
      {
         AxCor[nuc] = new complex<double>*[ARRAY];
         AyCor[nuc] = new complex<double>*[ARRAY];
         
         for(int x=0; x<ARRAY; x++)
         {
            AxCor[nuc][x] = new complex<double>[ARRAY];
            AyCor[nuc][x] = new complex<double>[ARRAY];
            
            for(int y=0; y<ARRAY; y++)
            {
               AxCor[nuc][x][y] = 0;
               AyCor[nuc][x][y] = 0;
            }
         }
      }
   }

   void RelAACor()
   {
      for(int nuc=0; nuc=2; nuc++)
      {
         for(int x=0; x<ARRAY; x++)
         {
            delete [] AxCor[nuc][x];
            delete [] AyCor[nuc][x];
         }
         delete [] AxCor[nuc];
         delete [] AyCor[nuc];
      }
      delete [] AxCor;
      delete [] AyCor;
   }


   //Definition for static double AlphaCor[2][ARRAY][ARRAY]
   //<alpha-alpha> correlation
   
   void DeclAlphaCor()
   {
      AlphaCor = new double**[2];
      for(int nuc=0; nuc<2; nuc++)
      {
         AlphaCor[nuc] = new double*[ARRAY];
         
         for(int x=0; x<ARRAY; x++)
         {
            AlphaCor[nuc][x] = new double[ARRAY];

            for(int y=0; y<ARRAY; y++)
            {
               AlphaCor[nuc][x][y] = 0;
            }
         }
      }
   }

   void RelAlphaCor()
   {
      for(int nuc=0; nuc=2; nuc++)
      {
         for(int x=0; x<ARRAY; x++)
         {
            delete [] AlphaCor[nuc][x];
         }
         delete [] AlphaCor[nuc];
      }
      delete [] AlphaCor;
   }

   //Definition for static double T00Av[ARRAY][ARRAY]
   //T00(0) average
   
   void DeclT00Av()
   {
      T00Av = new double*[ARRAY];
      for(int x=0; x<ARRAY; x++)
      {
         T00Av[x] = new double[ARRAY];

         for(int y=0; y<ARRAY; y++)
         {
            T00Av[x][y] = 0;
         }
      }
   }

   void RelT00Av()
   {
      for(int x=0; x<ARRAY; x++)
      {
         delete [] T00Av[x];
      }
      delete [] T00Av;
   }

   //Definition for static double RhoCG[2][NG][etamax][ARRAY][ARRAY]
   //Hold onto your rho's!
   
//    static complex<double> AxCor[2][ARRAY][ARRAY];
//    static complex<double> AyCor[2][ARRAY][ARRAY];

   void DeclRho()
   {
      RhoCG = new double****[2];
      for(int nuc=0; nuc<2; nuc++)
      {
         RhoCG[nuc] = new double***[NG];

         for(int col=0; col<NG; col++)
         {
            RhoCG[nuc][col] = new double**[etamax];
            
            for(int eta=0; eta<etamax; eta++)
            {
               RhoCG[nuc][col][eta] = new double*[ARRAY];

               for(int x=0; x<ARRAY; x++)
               {
                  RhoCG[nuc][col][eta][x] = new double[ARRAY]; 
                  for(int y=0; y<ARRAY; y++)
                  {
                     RhoCG[nuc][col][eta][x][y] = 0.;
                  }
               }
            }
         }

      }
   }



   void RelRho()
   {
      for(int nuc=0; nuc<2; nuc++)
      {
         for(int col=0; col<NG; col++)
         {
            for(int eta=0; eta<etamax; eta++)
            {
               for(int x=0; x<ARRAY; x++)
               {
                  delete [] RhoCG[nuc][col][eta][x];
               }
               delete [] RhoCG[nuc][col][eta];
            }
            delete [] RhoCG[nuc][col];
         }
         delete [] RhoCG[nuc];
      }
      delete [] RhoCG;
   }
    /**********/
    /*End Decl*/
    /**********/
    void WSIntegration ()
    {
       for(int x = 0; x<ARRAY; x++)
       {
           for(int y = 0; y<ARRAY; y++)
           {
              double WSInt = 0;
              for(int z=0; z<ARRAY*10; z++)
              {
                 WSInt = WSInt + NuclDens*(h/10)/(1+exp((sqrt((x-HALF)*(x-HALF)+(y-HALF)*(y-HALF)+z*z/100)*h-WSr)/WSa));
              }
              MuShape[0][x][y]=WSInt*2.;
              MuShape[1][x][y]=WSInt*2.;
           }
       }
    }

/*    
    complex<double>** Comm(complex<double>** X, complex<double>** Y) //only for use with 3x3 complex matrix
    {
       complex<double> Z[3][3];
    }
    */
    
    void NucSampler ()
    {
       Nucl1 = new double*[3];
       for(int i=0;i<3;i++)
       {
          Nucl1[i] = new double[N1];
          for(int n=0; n<N1; n++)
          {
             Nucl1[i][n] = 0;
          }
       }
       NucGen(N1, Nucl1);//Nucleus 1

       Nucl2 = new double*[3];
       for(int i=0;i<3;i++)
       {
          Nucl2[i] = new double[N2];
          for(int n=0; n<N2; n++)
          {
             Nucl2[i][n] = 0;
          }
       }
       NucGen(N2,Nucl2);//Nucleus 2

       for(int x = 0; x<ARRAY; x++)
       {
           for(int y = 0; y<ARRAY; y++)
           {
              double Shape0 = 0;
              double Shape1 =0;
              int OSInt = (ARRAY-1)/2; //Because Nucl# has negative positions, but we don't.
              double offset = OSInt*h; //Actual, factual offset distance
              
              for(int z=0; z<N1; z++)
              {
                 Shape0 = Shape0 + (1/(2*PIE*gauss_n))*exp( -((x*h - (Nucl1[0][z] + offset) )*(x*h - (Nucl1[0][z] + offset) ) +
                                                             (y*h - (Nucl1[1][z] + offset) )*(y*h - (Nucl1[1][z] + offset) )
                                                            )/(2*gauss_n) );
              }
              
              for(int z=0; z<N2; z++)
              {
                 Shape1 = Shape1 + (1/(2*PIE*gauss_n))*exp( -((x*h - (Nucl2[0][z] + offset) )*(x*h - (Nucl2[0][z] + offset) ) +
                                                             (y*h - (Nucl2[1][z] + offset) )*(y*h - (Nucl2[1][z] + offset) )
                                                            )/(2*gauss_n) );
              }
              
              MuShape[0][x][y]=Shape0;
              MuShape[1][x][y]=Shape1;
              
              MuShapeCenter[0][0][x][y] = MuShapeCenter[0][0][x][y] + MuShape[0][x][y];
              MuShapeCenter[1][0][x][y] = MuShapeCenter[1][0][x][y] + MuShape[1][x][y];
           }
       }

       NCollGeo=0;
       int Part[2][N1] = {0};
       for(int z=0; z<N1; z++)
       {
          for(int z2=0; z2 <N2; z2++)
          {
             if( (Nucl1[0][z] - Nucl2[0][z2] + 2*imp)*(Nucl1[0][z] - Nucl2[0][z2] +2*imp) + (Nucl1[1][z] - Nucl2[1][z2])*(Nucl1[1][z] - Nucl2[1][z2]) <= 1.5*1.5 )
             {
                NCollGeo++;
                
                Part[0][z] = 1; Part[1][z2]=1;
             }
          }
       }
       
       NPartGeo=0;
       for(int z=0; z<N1; z++)
       {
          NPartGeo = NPartGeo + Part[0][z] + Part[1][z];
       }
       
       GlauberGeo << NPartGeo << "  " << NCollGeo << endl;
              
       for(int x=2+imp+5; x<ARRAY-2-imp-5; x++)
       {
          for(int y=2; y<ARRAY-2; y++)
          {
             MuShapeCenter[0][1][x][y] = MuShapeCenter[0][1][x][y] + MuShape[0][x-imp][y]*MuShape[1][x+imp][y];
             MuShapeCenter[1][1][x][y] = MuShapeCenter[1][1][x][y] +
                                         MuShape[1][x+imp][y]*
                                         (
                                            MuShape[0][x-imp-2][y]
                                          - MuShape[0][x-imp-1][y]*8
                                          + MuShape[0][x-imp+1][y]*8
                                          - MuShape[0][x-imp+2][y]
                                         )/(12*h);
                                         
          }
       }
    }

    void CDFTable ()
    {
        double denom = sqrt(2);
        for(int i=0; i<NUMSTEP; i++)
        {
            CDFTab[i][0] = i*RhoStep - RhoMax;
            CDFTab[i][1] = (1+erf(CDFTab[i][0]/denom))/2.;
        }
    }

    
    bool SplitSort (double goal, int floor, int ceiling)
    {
        int TargetPoint = ( (floor + ceiling)/2 );
        double TargetVal = CDFTab[TargetPoint][1];
        
        if(goal > TargetVal)
        {
            return (true);
        }
        else
        {
            return (false);
        }
    }
    
    void ReInit()
    {
         /*
       for(int BigO=0; BigO<=(Order/2);BigO++) //Showtime
       {
         for(int a=0; a<4; a++)
         {
            for(int b=0; b<4; b++)
            {
               for(int x=0; x<ARRAY; x++)
               {
                  for(int y=0; y<ARRAY; y++)
                  {
                     T[BigO][a][b][x][y] = 0;
				  }
			   }
            }
         }
      }
     */
       for(int BigO=0; BigO<=(Order/2);BigO++) //Showtime
       {
         for(int em=0; em<2; em++)
         {
            for(int dir=0; dir<3; dir++)
            {
               for(int x=0; x<ARRAY; x++)
               {
                  for(int y=0; y<ARRAY; y++)
                  {
                     for(int a=0; a<3; a++)
                     {
                        for(int b=0; b<3; b++)
                        {
                           Fij[BigO][em][dir][x][y][a][b] = 0;
                        }
                     }
				  }
			   }
            }
         }
      }

      for(int BigO=0; BigO<=(Order/2);BigO++) //Showtime
      {
         for(int x=0; x<ARRAY; x++)
         {
            for(int y=0; y<ARRAY; y++)
            {
               for(int a=0; a<3; a++)
               {
                  for(int b=0; b<3; b++)
                  {
                     A2[BigO][x][y][a][b] = 0;
                     for(int i=0; i<2; i++)
                     {
                        Ai[BigO][i][x][y][a][b] = 0;
                     }
                  }
               }
            }
         }
      }
    }
    
    
    int main()
    {
       DeclDelMe();

       static complex<double> ATest[6][2][3][3];
//       ofstream DerVsCom;
//       DerVsCom.open("DerVsCom.txt");

       ofstream Params;
       Params.open("Parameters.txt");
       Params << "Order Calcd: " << Order << endl;
       Params << "Events Calcd: " << N << endl;
       Params << "Nucl Model : " << NucMethod  << endl;
       if(NucMethod==1)
       {
          Params << "====================" << endl;
          Params << "Nucleonic Information" << endl;
          Params << "Collision Type: ";
          if(Atom==1)
          {
              Params << Atom << "-> Pb-Pb" << endl;
          }
          else
          {
              Params << Atom << "-> Au-Au (NOTE: This is the default collision type!)" << endl;
          }
          Params << "Energy Scale Factor: "  << NucScale << endl;
          Params << "Radius of Nucleons : " << r_core << endl;
          Params << "Gaus Profile Radius: " << gauss_n << endl;
          Params << "====================" << endl;
       }
       if(NucMethod==0)
       {
          Params << "====================" << endl;
          Params << "Woods-Saxon Info" << endl;
          Params << "Woods-Saxon Decay a: " << WSa << endl;
          Params << "Woods-Saxon Radius : " << WSr << endl;
          Params << "Nucleon Density    : " << NuclDens << endl;
          Params << "====================" << endl;
       }
       if(NucMethod==2)
       {
          Params << "====================" << endl;
          Params << "Flat mu Info" << endl;
          Params << "Value of mu: " << mu << endl;
          Params << "====================" << endl;
       }

       Params << "Step Size  : " << h << endl;
       Params << "Array Size : " << ARRAY << endl;
       Params << "Array Half : " << HALF << " (should be above-1/2)" << endl;
//       Params << "CArray Size: " << CGRho << " (should = Array Size)" << endl;
//       Params << "CArray Half: " << CGHalf << " (should = Array Half)" << endl;
       Params << "ImpactPrm b: " << impact << endl;
       Params << "Value of g : " << g << endl;
       Params << "IR Cutoff m: " << mIR << " fm^-1" << endl;
       Params << "           : " << mIR*GEVFM << " GeV" << endl;
       Params << "Grain Molds: " << FFTWGrain << " (# of modes kept on either end)" << endl;       
       Params << "UV Cutoff Q: " << 2*PIE*FFTWGrain/(h*ARRAY) << " fm^-1 (from grain, array size and step size)" << endl;
       Params << "           : " << 2*PIE*FFTWGrain*GEVFM/(h*ARRAY) << " GeV" << endl;
       
       Params << "Sat.ScaleQs: " << mIR*sqrt(exp(4*PIE*12*PIE/(g*g*(11*3-2*4)))) << " fm^-1 (from grain, array size and step size)" << endl;
       Params << "           : " << GEVFM*mIR*sqrt(exp(4*PIE*12*PIE/(g*g*(11*3-2*4)))) << " GeV" << endl;
       
       Params << "Eta Slices : " << etamax << endl;
       Params << "Gluon Types: " << NG << endl;
       Params << "Color Types: " << NC << endl;
       Params.close();


       ofstream WriteTime;
       WriteTime.open("TimeDiag.txt");
       
       runtime= time(0);
       timer=time(0);
       ofstream RhoPrint;
//       GridRat = (int)round(hCG/h);

        cout << "start!";
       srand(time(0));


       RhoStep = 2*RhoMax/(NUMSTEP-1);
       CDFTable();
//	DeclMuShape();

       cout << CDFTab[780][0] << endl;
       cout << CDFTab[780][1] << endl;

//Read in Mold
/*
ifstream MoldRead;
MoldRead.open("MoldHD2Chen.txt");

   DeclMold();
   for(int x = 0; x<MOLDARRAY; x++)
   {
      for(int y = 0; y<MOLDARRAY; y++)
      {
         MoldRead >> Mold[x][y];
      }
   }

MoldRead.close();
*/

DeclAACor();
DeclAlphaCor();
DeclT();
//DeclT00Av();

/***************************************************/
/***************************************************/
/******************STATISTICAL**********************/
/********************READ IN************************/
/*******************PROCEDURE***********************/
/******************INITIALIZING*********************/
/***************************************************/
/***************************************************/
/***************************************************/
/*
       ifstream ReadRh[4];
       ReadRh[0].open("RhoCentCor.txt");
       ReadRh[1].open("RhoOffCor.txt");
       ReadRh[2].open("RhoAv.txt");
       ReadRh[3].open("RhoSqAv.txt");
       
       for(int x = 0; x<CGRho; x++)
       {
           for(int y = 0; y<CGRho; y++)
           {
              ReadRh[0] >> RhoCGCor[0][x][y];
              ReadRh[1] >> RhoCGCor[1][x][y];
              ReadRh[2] >> RhoCGAv[0][x][y];
              ReadRh[3] >> RhoCGAv[1][x][y];
           }
       }
       ReadRh[0].close();
       ReadRh[3].close();


       ifstream ReadE[4];
       ReadE[0].open("E0CentCor.txt");
       ReadE[1].open("E0Cor.txt");
       ReadE[2].open("E0Av.txt");
       ReadE[3].open("E0SqAv.txt");

       ifstream ReadB[4];
       ReadB[0].open("B0CentCor.txt");
       ReadB[1].open("B0OffCor.txt");
       ReadB[2].open("B0Av.txt");
       ReadB[3].open("B0SqAv.txt");
*/
//       ifstream ReadAx[4];
//       ReadAx[0].open("AxCentCor.txt");
//       ReadAx[1].open("AxOffCor.txt");
//       ReadAx[2].open("AxAv.txt");
//       ReadAx[3].open("AxSqAv.txt");
       
//       ifstream ReadAy[4];
//       ReadAy[0].open("AyCentCor.txt");
//       ReadAy[1].open("AyOffCor.txt");
       
//       ReadAy[2].open("AyAv.txt");
//       ReadAy[3].open("AySqAv.txt");
/*
       ifstream ReadAlp[4];
       ReadAlp[0].open("AlpCentCor.txt");
       ReadAlp[1].open("AlpOffCor.txt");
       ReadAlp[2].open("AlpAv.txt");
       ReadAlp[3].open("AlpSqAv.txt");
       
       ifstream ReadED[4];
       ReadED[0].open("EDCentCor.txt");
       ReadED[1].open("EDOffCor.txt");
       ReadED[2].open("EDAv.txt");
       ReadED[3].open("EDSqAv.txt");
*/
//       for(int x = 0; x<ARRAY; x++)
//       {
//          for(int y = 0; y<ARRAY; y++)
//          {
/*             ReadE[0] >> E0Cor[0][x][y];
             ReadE[1] >> E0Cor[1][x][y];
             ReadE[2] >> E0Av[0][x][y];
             ReadE[3] >> E0Av[1][x][y];
             
             ReadB[0] >> B0Cor[0][x][y];
             ReadB[1] >> B0Cor[1][x][y];
             ReadB[2] >> B0Av[0][x][y];
             ReadB[3] >> B0Av[1][x][y];
*/
//             ReadAx[0] >> AxCor[0][x][y];
//             ReadAx[1] >> AxCor[1][x][y];
//             ReadAx[2] >> AxAv[0][x][y];  
//             ReadAx[3] >> AxAv[1][x][y];
             
//             ReadAy[0] >> AyCor[0][x][y];
//             ReadAy[1] >> AyCor[1][x][y];
//             ReadAy[2] >> AyAv[0][x][y];
//             ReadAy[3] >> AyAv[1][x][y];
/*
             ReadAlp[0] >> AlphaCor[0][x][y];
             ReadAlp[1] >> AlphaCor[1][x][y];
             ReadAlp[2] >> AlphaAv[0][x][y];
             ReadAlp[3] >> AlphaAv[1][x][y];
             
             ReadED[0] >> EDCor[0][x][y];
             ReadED[1] >> EDCor[1][x][y];
             ReadED[2] >> EDAv[0][x][y];
             ReadED[3] >> EDAv[1][x][y];
             */
//          }
//       }
       /*
          ReadE[0].close();
          ReadE[1].close();
          ReadE[2].close();
          ReadE[3].close();
             
          ReadB[0].close();
          ReadB[1].close();
          ReadB[2].close();
          ReadB[3].close();
*/
//          ReadAx[0].close();
//          ReadAx[1].close();
//          ReadAx[2].close();
//          ReadAx[3].close();
          
//          ReadAy[0].close();
//          ReadAy[1].close();
//          ReadAy[2].close();
//          ReadAy[3].close();
/*
          ReadAlp[0].close();
          ReadAlp[1].close();
          ReadAlp[2].close();
          ReadAlp[3].close();
          
          ReadED[0].close();
          ReadED[1].close();
          ReadED[2].close();
          ReadED[3].close();

*/
/*
       ifstream ReadNumb;
       ReadNumb.open("Numb.txt");
       ReadNumb >> Numb;
       ReadNumb.close();
       
       Numb = Numb + N;
       
       ofstream WriteNumb;
       WriteNumb.open("Numb.txt");
       WriteNumb << Numb;
       WriteNumb.close();
*/

/***************************************************/
/***************************************************/
/******************STATISTICAL**********************/
/********************READ IN************************/
/*******************PROCEDURE***********************/
/*******************FINALIZED***********************/
/***************************************************/
/***************************************************/
/***************************************************/










//ifstream RhoRead;
//RhoRead.open("RhoSameReal.txt");



cout << "Fix 5" << endl;



       
       
       

ofstream Histogram[16];
ofstream EBHistogram[4];
ofstream EccHist[2];

ofstream PeakHist[2];

ofstream GlaubOut;
/*
Histogram[0].open("d0T00Hist.txt");
Histogram[1].open("d1T01Hist.txt");
Histogram[2].open("d2T02Hist.txt");
Histogram[3].open("d3T03Hist.txt");

Histogram[4].open("d0T10Hist.txt");
Histogram[5].open("d1T11Hist.txt");
Histogram[6].open("d2T12Hist.txt");
Histogram[7].open("d3T13Hist.txt");

Histogram[8].open("d0T20Hist.txt");
Histogram[9].open("d1T21Hist.txt");
Histogram[10].open("d2T22Hist.txt");
Histogram[11].open("d3T23Hist.txt");

Histogram[12].open("d0T30Hist.txt");
Histogram[13].open("d1T31Hist.txt");
Histogram[14].open("d2T32Hist.txt");
Histogram[15].open("d3T33Hist.txt");
*/
Histogram[0].open("T00Hist.txt");
Histogram[1].open("T01Hist.txt");
Histogram[2].open("T02Hist.txt");
Histogram[3].open("T03Hist.txt");
Histogram[4].open("T11Hist.txt");
Histogram[5].open("T12Hist.txt");
Histogram[6].open("T13Hist.txt");
Histogram[7].open("T22Hist.txt");
Histogram[8].open("T23Hist.txt");
Histogram[9].open("T33Hist.txt");
Histogram[10].open("T11DerHist.txt");
Histogram[11].open("T22DerHist.txt");
Histogram[12].open("L2Hist.txt");
Histogram[13].open("CSCharge.txt");
Histogram[14].open("CSChargeSquare.txt");
Histogram[15].open("CSChargeAbs.txt");

EBHistogram[0].open("ELHist.txt");
EBHistogram[1].open("BLHist.txt");
EBHistogram[2].open("ETHist.txt");
EBHistogram[3].open("BTHist.txt");

EccHist[0].open("Ecc2.txt");
EccHist[1].open("Ecc3.txt");

PeakHist[0].open("MonopointEnergy.txt");
PeakHist[1].open("MultipointEnergy.txt");

GlaubOut.open("GlauberInfo.txt");
GlauberGeo.open("GlauberInfoGeo.txt");

       for(int n=0; n<N; n++)
       {
//          ReInit(); //Because of terrible memory management techniques

       if(NucMethod == 0)
       {
          WSIntegration();
       }
       
       if(NucMethod == 1)
       {
          NucSampler();
          
          /*
          ofstream MuShWrite[2];
          MuShWrite[0].open("ConMuSh0.txt");
          MuShWrite[1].open("ConMuSh1.txt");
          
          for(int nuc=0;nuc<2;nuc++)
          {
          for(int x=0;x<CGRho;x++)
          {
          for(int y=0;y<CGRho;y++)
          {
             MuShWrite[nuc] << MuShape[nuc][x][y] << "  ";
          }
             MuShWrite[nuc] << endl;
          }}
          
          MuShWrite[0].close();
          MuShWrite[1].close();
          */
       }

       if(NucMethod == 2)
       {
          for(int nuc=0;nuc<2;nuc++)
          {
          for(int x=0;x<ARRAY;x++)
          {
          for(int y=0;y<ARRAY;y++)
          {
          MuShape[nuc][x][y]=mu;
          }
          }
          }
       }
       

       NColl = 0;
       NPart = 0;
//       for(int x=0;x<ARRAY;x++)
       for(int x=imp; x<ARRAY-imp; x++)
       {
       for(int y=0;y<ARRAY;y++) // Calc at x,y
       {
       CollDens = XSect*MuShape[0][x-imp][y]*MuShape[1][x+imp][y];

       PartDens = MuShape[0][x-imp][y]*(1-pow(1-MuShape[1][x+imp][y]*XSect/N2,N2))
                + MuShape[1][x+imp][y]*(1-pow(1-MuShape[0][x-imp][y]*XSect/N1,N1));
            
       NColl = NColl + CollDens;
       NPart = NPart + PartDens;
       }
       }
       
       GlaubOut << NPart*h*h << "  " << NColl*h*h << endl;

       DeclCovPot();

       timer=time(0) - timer;
       
       WriteTime << "Starting run for N= " << n << endl;
       WriteTime << "Initial overhead: " << timer << endl;

       timer=time(0);
       
       for(int nuc=0; nuc<2; nuc++)
       {
       for(int eta=0; eta<etamax; eta++)
       {

       cout << "Ensemble: " << n << endl;

       DeclRho();

       for(int x = 0; x<ARRAY; x++)
//       for(int x = 0; x<(ARRAY-1)/2; x++)
       {
           for(int y = 0; y<ARRAY; y++)
           {
               for(int col = 0; col<8; col++)
               {
                   double RhoRoll = ran();

                   
                   if(RhoRoll < CDFTab[0][1] )
                   {
                      RhoCG[nuc][col][eta][x][y] = CDFTab[0][0];
                   }
                   
                   else if(RhoRoll > CDFTab[NUMSTEP-1][1] )
                   {
                      RhoCG[nuc][col][eta][x][y] = CDFTab[NUMSTEP-1][0];
                   }


                   else
                   {
                      flr = 0;
                      clg = NUMSTEP - 1;

                      for(int i = 0; i<10;i++)
                      {
                         if(SplitSort(RhoRoll, flr, clg)== true)
                         {
                            flr = ( (flr + clg)/2 );
                         }
            
                         else
                         {
                            clg = ( (flr + clg)/2 );
                         }
                      }

                      RhoCG[nuc][col][eta][x][y] = RhoStep*(RhoRoll - CDFTab[flr][1])/(CDFTab[clg][1] - CDFTab[flr][1]) + CDFTab[flr][0];


//Fabricating results of roll for testing
/*
                      if(x==CGHalf && y==CGHalf)
                      {
                          RhoCG[nuc][col][eta][x][y] = 1.0;
                      }
                      else
                      {
                          RhoCG[nuc][col][eta][x][y] = 0.0;
                      }
*/
//Fabrication Complete
                      
//                      RhoCGAv[col][x][y] = RhoCGAv[col][x][y] + RhoCG[col][x][y]/N;                           //Average for testing
//                      RhoCGVar[col][x][y] = RhoCGVar[col][x][y] + RhoCG[col][x][y]*RhoCG[col][x][y]/N;        //Variance for testing

//if(x<6 || x>12 || y<6 || y>12)
//{
//   RhoCG[col][x][y] = 0;
//}

                      //RhoCG[7][51][51]= 1.0;
                      /*DEBUG LINE*/
//                      RhoRead >> RhoCG[col][x][y];
                   }
                   //Scale by Sqrt(mu) to scale strength
                   //Scale by hCG to account for coarse graining (1/hCG) and include summation element for sum in alpha (hCG*hCG)
                   //
                   RhoCG[nuc][col][eta][x][y] = RhoCG[nuc][col][eta][x][y]*sqrt(NucScale*MuShape[nuc][x][y])*g*h;
//                   RhoCG[nuc][col][eta][ARRAY-1-x][y] = -RhoCG[nuc][col][eta][x][y];
               }
           }//close y loop
       }//close x loop
/*
ifstream RhoRead[2];
RhoRead[0].open("col0.txt");
RhoRead[1].open("col1.txt");

   for(int col = 0; col<8; col++)
   {
         for(int x = 0; x<CGRho; x++)
         {
            for(int y = 0; y<CGRho; y++)
            {
               if(col ==0)
               {
                  if(nuc==0)
                  {
                     RhoRead[0]>> RhoCG[nuc][col][eta][x][y];
                  }
                  else
                  {
                     RhoRead[1]>> RhoCG[nuc][col][eta][x][y];
                  }
               }
               if(col == 1)
               {
                  if(nuc==0)
                  {
                     RhoRead[1]>> RhoCG[nuc][col][eta][x][y];
                  }
                  else
                  {
                     RhoRead[0]>> RhoCG[nuc][col][eta][x][y];
                  }
               }
               if(col>1)
               {
                 if(col%2 ==0)
                  {
                    RhoCG[nuc][col][eta][x][y] = RhoCG[nuc][0][eta][x][y];
                 }
                 else
                 {
                    RhoCG[nuc][col][eta][x][y] = RhoCG[nuc][1][eta][x][y];
                 }
               }
//               else
//               {
//                  RhoCG[nuc][col][eta][x][y]=RhoCG[0][col][0][x][y];
//               }
            }
         }
   }

RhoRead[0].close();
RhoRead[1].close();
*/

/*
       if(nuc==0 && false)
       {
          ifstream Nuc0Write;
          if(eta==0)
          {
             Nuc0Write.open("ConNuc00.txt");
          }
          if(eta==1)
          {
             Nuc0Write.open("ConNuc01.txt");
          }
          if(eta==2)
          {
             Nuc0Write.open("ConNuc02.txt");
          }
          if(eta==3)
          {
             Nuc0Write.open("ConNuc03.txt");
          }
          if(eta==4)
          {
             Nuc0Write.open("ConNuc04.txt");
          }

          for(int x = 0; x<CGRho; x++)
          {
             for(int y = 0; y<CGRho; y++)
             {
                for(int col = 0; col<8; col++)
                {
                   Nuc0Write >> RhoCG[nuc][col][eta][x][y];// << "  ";
                }
//                Nuc0Write << endl;
             }
          }
          Nuc0Write.close();
       }

       if(nuc==1 && false)
       {
          ifstream Nuc1Write;
          if(eta==0)
          {
             Nuc1Write.open("ConNuc10.txt");
          }
          if(eta==1)
          {
             Nuc1Write.open("ConNuc11.txt");
          }
          if(eta==2)
          {
             Nuc1Write.open("ConNuc12.txt");
          }
          if(eta==3)
          {
             Nuc1Write.open("ConNuc13.txt");
          }
          if(eta==4)
          {
             Nuc1Write.open("ConNuc14.txt");
          }
          
          for(int x = 0; x<CGRho; x++)
          {
             for(int y = 0; y<CGRho; y++)
             {
                for(int col = 0; col<8; col++)
                {
                   Nuc1Write >> RhoCG[nuc][col][eta][x][y];// << "  ";
                }
//                Nuc1Write << endl;
             }
          }
          Nuc1Write.close();
       }
*/
       
/*
       for(int x = 0; x<CGRho; x++)
       {
           for(int y = 0; y<CGRho; y++)
           {
               for(int col = 0; col<8; col++)
               {
//                  RhoCGCor[0][x][y] = RhoCGCor[0][x][y] + RhoCG[nuc][col][eta][x][y]*RhoCG[nuc][col][eta][CGHalf][CGHalf];
//                  RhoCGCor[1][x][y] = RhoCGCor[1][x][y] + RhoCG[nuc][col][eta][x][y]*RhoCG[nuc][col][eta][(CGHalf/2)][(CGHalf/2)];
                  RhoCGAv[0][x][y] = RhoCGAv[0][x][y] + RhoCG[nuc][col][eta][x][y];
                  RhoCGAv[1][x][y] = RhoCGAv[1][x][y] + RhoCG[nuc][col][eta][x][y]*RhoCG[nuc][col][eta][x][y];
               }
           }
       }
*/
/*
ofstream RhoWrite;
RhoWrite.open("RhoSameReal.txt");
       /**Write Rho debugging**/
/*
       for(int x = 0; x<CGRho; x++)
       {
               RhoWrite << endl;
           for(int y = 0; y<CGRho; y++)
           {
               for(int col = 7; col<8; col++)
               {
                  RhoWrite << RhoCG[col][x][y] << "  ";
               }
           }
       }
*/
//RhoCG[7][CGHalf][CGHalf]= 1.0;
//RhoCG[8][CGHalf][CGHalf] = 1.0;
//RhoCG[9][CGHalf][CGHalf]= 1.0;
//Resetting for iteration loop
/*
       for(int x = 1; x<ARRAY-1; x++)
       {
           for(int y = 1; y<ARRAY-1; y++)
           {
               for(int col = 7; col<8; col++)
               {
                  Rho[col][x][y] = 0;
               }
           }
       }


// Charge Smearing
/*
       for(int x = 1; x<ARRAY-1; x++)
       {
           for(int y = 1; y<ARRAY-1; y++)
           {
               for(int col = 7; col<8; col++)
               {
                  for(int sx=0; sx<CGRho; sx++)
                  {
                     deltx = ((8+0-1)*h - sx*hCG);
                     for(int sy=0; sy<CGRho; sy++)
                     {
                        delty = ((y+0-1)*h - sy*hCG);
                        Rho[col][x][y] = Rho[col][x][y] + RhoCG[col][sx][sy]*exp( -(deltx*deltx + delty*delty)/(Lam*Lam) );
                     }
                  }
               }
           }
       }
*/       

          //Corelation Calculation
  /*        
          for(int x = 1; x<ARRAY-1; x++)
          {
             for(int y = 1; y<ARRAY-1; y++)
             {
                RhoCor[0][x][y]=RhoCor[0][x][y] + Rho[7][HALF][HALF]*Rho[7][x][y]/N;
                RhoCor[1][x][y]=RhoCor[1][x][y] + Rho[7][OFF][HALF]*Rho[7][x][y]/N;
             }
          }
       ofstream RhoPrint[9];
*/
//       RhoPrint[0].open("Rho0.txt");
//       RhoPrint[1].open("Rho1.txt");
//       RhoPrint[2].open("Rho2.txt");
//       RhoPrint[3].open("Rho3.txt");
//       RhoPrint[4].open("Rho4.txt");
//       RhoPrint[5].open("Rho5.txt");
//       RhoPrint[6].open("Rho6.txt");
/*       RhoPrint[7].open("Rho7-400.txt");
//       RhoPrint[8].open("Rho8Test.txt");

       for(int col = 7; col<8; col++)
       {
           for(int x = 1; x<ARRAY-1; x++)
           {
               RhoPrint[col] << endl;
               for(int y = 1; y<ARRAY-1; y++)
               {
                  RhoPrint[col] << Rho[col][x][y] << "   ";
               }
           }
       }
*/
//       RhoPrint[0].close();
//       RhoPrint[1].close();
//       RhoPrint[2].close();
//       RhoPrint[3].close();
//       RhoPrint[4].close();
//       RhoPrint[5].close();
//       RhoPrint[6].close();


//       Rho[8][41][51]= -1.0;
//       Rho[8][61][51]= 1.0;
//       Rho[8][HALF][HALF]= 1.0;

//Iterative Method (Lame)

if(Target == 1 )
{
/*   OLDE METHODE
       static int counter;
       static int count;
       for(int iter = 1; iter<1/*75000*//*; iter++)
       {
          counter++;
          if( counter == 100 )
          {
             count++;
             cout << "Iteration: " << count*100 << endl;
             counter=0;
          }
          for(int x = 1; x<ARRAY-1; x++)
          {
              for(int y = 1; y<ARRAY-1; y++)
              {
                  for(int col = 0; col<8; col++)
                  {
                     CovPot[0][col][x][y]=(1./4.)*( CovPot[1][col][x+1][y] + CovPot[1][col][x-1][y] +
                                                    CovPot[1][col][x][y+1] + CovPot[1][col][x][y-1] +
                                                    h*h*g*Rho[col][x][y]
                                                   );
                  }
              }
          }
          


          for(int x = 1; x<ARRAY-1; x++)
          {
              for(int y = 1; y<ARRAY-1; y++)
              {
                  for(int col = 0; col<9; col++)
                  {
                     CovPot[1][col][x][y]=CovPot[0][col][x][y];
                  }
              }
          }
       }
       
       */

// Point Charge superposition approximation
/*
          for(int x = 1; x<ARRAY-1; x++)
          {
              for(int y = 1; y<ARRAY-1; y++)
              {
                 for(int sx = 1; sx<ARRAY-1; sx++)
                 {
                    int dx = x - sx;
                    if(dx+HALF > 0 && dx+HALF<ARRAY)
                    {
                       for(int sy = 1; sy<ARRAY-1; sy++)
                       {
                          int dy = y - sy;
                          if(dy+HALF >=0 && dy+HALF<ARRAY-1)
                          CovPot[0][9][x][y]=CovPot[0][9][x][y] + Rho[7][sx][sy]*CovPot[0][8][dx+HALF][dy+HALF];
                       }
                    }
                 }
              }
          }
*/

// Reset Potential
/*Needed for olde methodology
       for(int x = 1; x<ARRAY-1; x++)
       {
           for(int y = 1; y<ARRAY-1; y++)
           {
              for(int col = 0; col<8; col++)
              {
                 CovPot[0][col][eta][x][y] = 0;
              }
           }
       }
*/

// INITIALIZE.

/*cout << "Starting CovPot initialization..." << endl;
          for(int x = 0; x<ARRAY; x++)
          {
              for(int y = 0; y<ARRAY; y++)
              {
                  for(int col = 0; col<8; col++)
                  {
                     CovPot[nuc][col][eta][x][y]=0;
                  }
              }
          }
*/


// Smeared Charge superposition approximation

/**************************************/

/**********AS IT WAS SO IT IS**********/

/**************************************/
/*
          cout << "Potential for nuc = " << nuc << " eta = " << eta << endl;
          for(int sx = 0; sx<CGRho; sx++)
          {
              int MaxX = GridRat*sx;
              int MinX = MaxX - MOLDRANGE; // Range limit for Green's f'n
              MaxX = MaxX + MOLDRANGE + 1; // Range limit for Green's f'n
              
              //check limits vs actual array and adjust
              if(MinX < 0)
              {
                 MinX = 0;
              }
              
              if(MaxX > ARRAY)
              {
                 MaxX = ARRAY;
              }
              
              for(int sy = 0; sy<CGRho; sy++)
              {
                 int MaxY = GridRat*sy;
                 int MinY = MaxY - MOLDRANGE; // Range limit for Green's f'n
                 MaxY = MaxY + MOLDRANGE + 1; // Range limit for Green's f'n
              
                 //check limits vs actual array and adjust
                 if(MinY < 0)
                 {
                    MinY = 0;
                 }
              
                 if(MaxY > ARRAY)
                 {
                    MaxY = ARRAY;
                 }
                 
                 if( (nuc==0 && MuShape[nuc][sx][sy]>(0.01*pow(N1, (1/3.))/(PIE*gauss_n)))
                   ||(nuc==1 && MuShape[nuc][sx][sy]>(0.01*pow(N2, (1/3.))/(PIE*gauss_n)))
                   )
                 {
                 for(int x = MinX; x<MaxX; x++)
                 {
                    int dx = ((x) - GridRat*sx);

                       for(int y = MinY; y<MaxY; y++)
                       {
                          int dy = ((y) - GridRat*sy);

                             for(int col = 0; col<8; col++)
                             {
                                CovPot[nuc][col][eta][x][y]=CovPot[nuc][col][eta][x][y] - RhoCG[nuc][col][eta][sx][sy]*Mold[dx+MOLDHALF][dy+MOLDHALF];
/*Calculating rho is costly!                          if(eta==0)
                          {
                             Rho[nuc][col][x][y] = Rho[nuc][col][x][y] + RhoCG[nuc][col][eta][sx][sy]*exp( -(deltx*deltx + delty*delty)/(Lam*Lam) );
                          }
*/



/* MIDDLE PIECE



                             }
//                       CovPotAv[x][y] = CovPotAv[x][y] + CovPot[0][9][x][y]/N;                                  // Average for testing (monochrome)
//                       CovPotVar[x][y] = CovPotVar[x][y] + CovPot[0][9][x][y]*CovPot[0][9][x][y]/N;             //Variance for testing (monochrome)
                       }
                 }
                 }//End if enough charge
              }
          }
DON'T FORGET MIDDLE PIECE*/

/*
          for(int sx = CGRho; sx<CGRho; sx++)
          {
              int MaxX = GridRat*sx;
              int MinX = 0;//MaxX - MOLDRANGE; // Range limit for Green's f'n
              MaxX = ARRAY;//MaxX + MOLDRANGE + 1; // Range limit for Green's f'n
              
              //check limits vs actual array and adjust
              if(MinX < 0)
              {
                 MinX = 0;
              }
              
              if(MaxX > ARRAY)
              {
                 MaxX = ARRAY;
              }
              
              for(int sy = 0; sy<CGRho; sy++)
              {
                 int MaxY = GridRat*sy;
                 int MinY = 0;//MaxY - MOLDRANGE; // Range limit for Green's f'n
                 MaxY = ARRAY;//MaxY + MOLDRANGE + 1; // Range limit for Green's f'n
              
                 //check limits vs actual array and adjust
                 if(MinY < 0)
                 {
                    MinY = 0;
                 }
              
                 if(MaxY > ARRAY)
                 {
                    MaxY = ARRAY;
                 }
                 
                 if( true//nuc==0 && MuShape[nuc][sx][sy]>(0.01*pow(N1, (1/3.))/(PIE*gauss_n)))
                   //||(nuc==1 && MuShape[nuc][sx][sy]>(0.01*pow(N2, (1/3.))/(PIE*gauss_n)))
                   )
                 {
                 for(int x = MinX; x<MaxX; x++)
                 {
                    int dx = ((x) - GridRat*sx);

                       for(int y = MinY; y<MaxY; y++)
                       {
                          int dy = ((y) - GridRat*sy);

                                CovPot[nuc][col][eta][x][y]=CovPot[nuc][col][eta][x][y] + RhoCG[nuc][col][eta][sx][sy]*Mold[dx+MOLDHALF][dy+MOLDHALF];
                             }
                       }
                 }
                 }//End if enough charge
              }
*/



       for(int col = 0; col<8; col++)
       {
          cout << "Potential for nuc = " << nuc << endl;

          if(nuc==0 && col==0 && eta==0)
          {
          /*
              ofstream WriteAlpha;
              WriteAlpha.open("XGreen.txt");

              for(int x = 0; x<ARRAY; x++)
              {
                 for(int y = 0; y<ARRAY; y++)
                 {
                    WriteAlpha << CovPot[nuc][col][eta][x][y] << "  ";
                 }
                 WriteAlpha << endl;
              }
           
              WriteAlpha.close();
           */
              for(int x = 0; x<ARRAY; x++)
              {
                 for(int y = 0; y<ARRAY; y++)
                 {
                    CovPot[nuc][col][eta][x][y] = 0.;
                 }
              }
          }
          
          
          
          fftw_complex *rhoalpha; //Will contain rho(x), rho(k), alpha(k), alpha(x)
//          fftw_complex *rhok;
//          double rho[CGRho*CGRho];
//          double alpha[CGRho*CGRho];
          fftw_plan p, q;


          rhoalpha = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ARRAY*ARRAY); // what the fuck is the 1d bullshit...


          for(int i=0; i<ARRAY; i++)
          {
             for(int j=0; j<ARRAY; j++)
             {
                rhoalpha[i*ARRAY + j][0] = RhoCG[nuc][col][eta][i][j]; //rho(x)
                rhoalpha[i*ARRAY + j][1] = 0;
             }
          }


          p = fftw_plan_dft_2d(ARRAY, ARRAY, rhoalpha, rhoalpha, FFTW_FORWARD, FFTW_ESTIMATE);
          q = fftw_plan_dft_2d(ARRAY, ARRAY, rhoalpha, rhoalpha, FFTW_BACKWARD, FFTW_ESTIMATE);
//          p = fftw_plan_dft_r2c_2d((CGRho+1)/2, (CGRho+1)/2, rho, rhok, FFTW_ESTIMATE);
//          q = fftw_plan_dft_c2r_2d((CGRho+1)/2, (CGRho+1)/2, rhok, alpha, FFTW_ESTIMATE);


          fftw_execute(p); //rhoalpha -> rho(k)

          double mterm = mIR*mIR*h*h;
                    
          for(int i=0; i<ARRAY; i++) //rhoalpha -> alpha(k)
          {
             double ki = 2*PIE*(min(i,ARRAY-i))/(ARRAY);
              for(int j=0; j<ARRAY; j++)
             {
/*
                if( i<=80 || CGRho-i <=80)
                {
                   if( j<=80 || CGRho-j <=80)
                   {
*/

                if( i<=FFTWGrain || ARRAY-i <=FFTWGrain )
                {
                   if( j<=FFTWGrain || ARRAY-j <=FFTWGrain )
                   {
                      double kj = 2*PIE*(min(j,ARRAY-j))/(ARRAY);
                      rhoalpha[i*ARRAY + j][0] = rhoalpha[i*ARRAY + j][0]/(ki*ki + kj*kj + mterm);
                      rhoalpha[i*ARRAY + j][1] = rhoalpha[i*ARRAY + j][1]/(ki*ki + kj*kj + mterm);
		        
                      if(i==j && j==0 && true)
                      {
		                 rhoalpha[i*ARRAY + j][0] = 0;
		                 rhoalpha[i*ARRAY + j][1] = 0;
		              }
                   }
                   else
                   {
                      rhoalpha[i*ARRAY + j][0] = 0;
                      rhoalpha[i*ARRAY + j][1] = 0;
                   }
                }
                else
                {
                   rhoalpha[i*ARRAY + j][0] = 0;
                   rhoalpha[i*ARRAY + j][1] = 0;
                }
             }
          }

          fftw_execute(q); //rhoalpha -> alpha(x)

/*
if(eta==0 && col==1)
{
          ofstream Reality;
	  ofstream Madness;
          Reality.open("FFTWRe.txt");
          Madness.open("FFTWIm.txt");
          for(int i=0; i<CGRho; i++)
          {
             for(int j=0; j<CGRho; j++)
             {
                Reality << rhoalpha[i*CGRho + j][0]/(CGRho*CGRho) << "  ";
                Madness << rhoalpha[i*CGRho + j][1]/(CGRho*CGRho) << "  ";
	     }
             Reality << endl;
             Madness << endl;
	  }
}
*/

          for(int i=0; i<ARRAY; i++)
          {
             for(int j=0; j<ARRAY; j++)
             {
                CovPot[nuc][col][eta][i][j] = rhoalpha[i*ARRAY + j][0]/(ARRAY*ARRAY);
             }
          }

/*
          if(nuc==0 && col==0 && eta==0)
          {
              ofstream WriteAlpha;
              WriteAlpha.open("PGreen.txt");

              for(int x = 0; x<ARRAY; x++)
              {
                 for(int y = 0; y<ARRAY; y++)
                 {
                    WriteAlpha << CovPot[nuc][col][eta][x][y] << "  ";
                 }
                 WriteAlpha << endl;
              }
              
              WriteAlpha.close();
              for(int x = 0; x<ARRAY; x++)
              {
                 for(int y = 0; y<ARRAY; y++)
                 {
                    CovPot[nuc][col][eta][x][y] = 0.;
                 }
              }
          }
*/
          fftw_free(rhoalpha);
          
          fftw_destroy_plan(p);
          fftw_destroy_plan(q);
}//end col loop


















          for(int x = 0; x<ARRAY; x++)
          {
             for(int y = 0; y<ARRAY; y++)
             {
                for(int col = 0; col<1; col++)
                {
                   AlphaCor[0][x][y] = AlphaCor[0][x][y] + CovPot[nuc][col][eta][x][y]*CovPot[nuc][col][eta][HALF][HALF];
//                   AlphaCor[1][x][y] = AlphaCor[1][x][y] + CovPot[nuc][col][eta][x][y]*CovPot[nuc][col][eta][(HALF/2)][(HALF/2)];
                   AlphaCor[1][x][y] = AlphaCor[1][x][y] + CovPot[nuc][col][eta][x][y];
//                   AlphaAv[0][x][y] = AlphaAv[0][x][y] + CovPot[nuc][col][eta][x][y];
//                   AlphaAv[1][x][y] = AlphaAv[1][x][y] + CovPot[nuc][col][eta][x][y]*CovPot[nuc][col][eta][x][y];
                }
             }
          }


} // If Targ 1

RelRho();

/*       ofstream WriteLap;
       WriteLap.open("LaplaceSolved.txt");
       ofstream WritePnt;
       WritePnt.open("PointSolved.txt");

        for(int x = 1; x<102; x++)
        {
            WritePnt << endl;
            WriteLap << endl;

            for(int y = 1; y<102; y++)
            {
               WriteLap << CovPot[0][7][x][y] << "   ";
               WritePnt << CovPot[0][9][x][y] << "   ";
            }
        }
*/
/*       ofstream WriteEx;
       WriteEx.open("ExactSolved.txt");

        for(int x = 0; x<102; x++)
        {
            WriteEx << endl;

            for(int y = 0; y<102; y++)
            {
               WriteEx << CovPotChk[x][y] << "   ";
            }
        }
*/

//Potential Corelation

/*//Only for testing
int countersame =0;
int counterdiff=0;

          for(int x = 0; x<ARRAY-1; x++)
          {
             for(int y = 0; y<ARRAY-1; y++)
             {
                PotCor[0][2][x][y]=PotCor[0][2][x][y] + CovPotChk[HALF][HALF]*CovPotChk[x][y]/N;
                PotCor[1][2][x][y]=PotCor[1][2][x][y] + CovPotChk[OFF][HALF]*CovPotChk[x][y]/N;
                if( x!=0 && y!=0)
                {
                    PotCor[0][0][x][y]=PotCor[0][0][x][y] + CovPot[0][7][HALF][HALF]*CovPot[0][7][x][y]/N;
                    PotCor[0][1][x][y]=PotCor[0][1][x][y] + CovPot[0][9][HALF][HALF]*CovPot[0][9][x][y]/N;
                    /*
                    if( ( (CovPot[0][9][HALF][HALF]>0) && (CovPot[0][9][x][y] <0) ) || 
                        ( (CovPot[0][9][HALF][HALF]<0) && (CovPot[0][9][x][y] >0) )
                      )
                    {
                      cout << "Yyyyep." << endl;
                       cout << CovPot[0][9][HALF][HALF] << endl;
                       cout << CovPot[0][9][x][y] << endl;
                       cout << CovPot[0][9][HALF][HALF]*CovPot[0][9][x][y]/N << endl;
                       cout << PotCor[0][1][x][y] << endl;
                    }
                    */
/*
                    PotCor[1][1][x][y]=PotCor[1][1][x][y] + CovPot[0][9][OFF][HALF]*CovPot[0][9][x][y]/N;
                    if(x == 51 && y == 51 ){
                    cout << "Center: " << CovPot[0][7][51][51] << endl;
                    cout << "Reference: " << CovPot[0][7][x][y] << endl;
                    cout << "Correlation: " << PotCor[0][1][x][y] << endl << endl;

                    }
                }

             }
          }
*/
          
          
cout << "N: " << n << " 's eta: " << eta  << "complete" << endl;

          } //end eta loop

/*
          ofstream CovPrint[4];
          CovPrint[0].open("CovPotSingle0.txt");
          CovPrint[1].open("CovPotSingle1.txt");
          CovPrint[2].open("CovPotFull0.txt");
          CovPrint[3].open("CovPotFull1.txt");
          for(int x = 0; x<ARRAY; x++)
          {
              for(int y = 0; y<ARRAY; y++)
              {
                 double Sing0 = 0.;
                 double Sing1 = 0.;
                 double Full0 = 0.;
                 double Full1 = 0.;
                 
                 for(int col = 0; col<8; col++)
                 {
                    Sing0 += CovPot[0][col][0][x][y];
                    Sing1 += CovPot[1][col][0][x][y];
                    
                    for(int eta = 0; eta<etamax; eta++)
                    {
                       Full0 += CovPot[0][col][eta][x][y];
                       Full1 += CovPot[1][col][eta][x][y];
                    }
                 }
                 CovPrint[0] << Sing0 << "  ";
                 CovPrint[1] << Sing1 << "  ";
                 CovPrint[2] << Full0 << "  ";
                 CovPrint[3] << Full1 << "  ";
              }
              CovPrint[0] << endl;
              CovPrint[1] << endl;
              CovPrint[2] << endl;
              CovPrint[3] << endl;
          }
          CovPrint[0].close();
          CovPrint[1].close();
          CovPrint[2].close();
          CovPrint[3].close();
*/
          
          } //end nuc loop
          

       timer=time(0) - timer;
       
       WriteTime << "Solving Alpha: " << timer << endl;

       timer=time(0);

       

//          RelMold();
/*
          ofstream PrintRho;
          ofstream PrintAlpha;
          PrintRho.open("RhoSample.txt");
          PrintAlpha.open("AlphaSample.txt");
          
          cout << "Printing Rho/Alpha..." << endl;
          for(int x = 0; x<ARRAY; x++)
          {
              for(int y = 0; y<ARRAY; y++)
              {
              //   PrintRho << Rho[0][0][x][y] << "  ";
              //   cout << "rho: " << Rho[0][0][x][y] << endl;
                 PrintAlpha << CovPot[0][0][0][x][y] << "  ";
              }
              //PrintRho << endl;
              PrintAlpha << endl;
          }
          PrintRho.close();
          PrintAlpha.close();
          cout << "Rho/Alpha Printed!" << endl;
          
          ofstream PrintRhoCG;
          PrintRhoCG.open("RhoCGSample.txt");
          for(int x = 0; x<CGRho; x++)
          {
              for(int y = 0; y<CGRho; y++)
              {
                 PrintRhoCG << RhoCG[0][0][0][x][y] << "  ";
                 cout << "rhoCG: " << RhoCG[0][0][0][x][y] << endl;
              }
              PrintRhoCG << endl;
          }
          PrintRhoCG.close();
*/
          //U=exp(-ig(n.alpha))
          //exp(i*theta*H)= (H^2-i*H*d/dtheta-I(1+d/dtheta^2))
          //DON'T DO ABOVE. USE EQ 4.
          //Tr[H^2]=2
          //Gellman Matricies
          
          

          complex<double> I(0.,1.);
          complex<double> t[8][3][3]; //8 3x3s
          //Tr[tatb]==2d_a,b
          int cl=0;
          t[cl][0][0] = 0.; t[cl][0][1] = 1.; t[cl][0][2] = 0.;
          t[cl][1][0] = 1.; t[cl][1][1] = 0.; t[cl][1][2] = 0.;
          t[cl][2][0] = 0.; t[cl][2][1] = 0.; t[cl][2][2] = 0.;

          cl=1;
          t[cl][0][0] = 0.; t[cl][0][1] = -I; t[cl][0][2] = 0.;
          t[cl][1][0] = I; t[cl][1][1] = 0.; t[cl][1][2] = 0.;
          t[cl][2][0] = 0.; t[cl][2][1] = 0.; t[cl][2][2] = 0.;

          cl=2;
          t[cl][0][0] = 1.; t[cl][0][1] = 0.; t[cl][0][2] = 0.;
          t[cl][1][0] = 0.; t[cl][1][1] = -1.; t[cl][1][2] = 0.;
          t[cl][2][0] = 0.; t[cl][2][1] = 0.; t[cl][2][2] = 0.;

          cl=3;
          t[cl][0][0] = 0.; t[cl][0][1] = 0.; t[cl][0][2] = 1.;
          t[cl][1][0] = 0.; t[cl][1][1] = 0.; t[cl][1][2] = 0.;
          t[cl][2][0] = 1.; t[cl][2][1] = 0.; t[cl][2][2] = 0.;

          cl=4;
          t[cl][0][0] = 0.; t[cl][0][1] = 0.; t[cl][0][2] = -I;
          t[cl][1][0] = 0.; t[cl][1][1] = 0.; t[cl][1][2] = 0.;
          t[cl][2][0] = I; t[cl][2][1] = 0.; t[cl][2][2] = 0.;

          cl=5;
          t[cl][0][0] = 0.; t[cl][0][1] = 0.; t[cl][0][2] = 0.;
          t[cl][1][0] = 0.; t[cl][1][1] = 0.; t[cl][1][2] = 1.;
          t[cl][2][0] = 0.; t[cl][2][1] = 1.; t[cl][2][2] = 0.;

          cl=6;
          t[cl][0][0] = 0.; t[cl][0][1] = 0.; t[cl][0][2] = 0.;
          t[cl][1][0] = 0.; t[cl][1][1] = 0.; t[cl][1][2] = -I;
          t[cl][2][0] = 0.; t[cl][2][1] = I; t[cl][2][2] = 0.;

          cl=7;
          t[cl][0][0] = 1.; t[cl][0][1] = 0.; t[cl][0][2] = 0.;
          t[cl][1][0] = 0.; t[cl][1][1] = 1.; t[cl][1][2] = 0.;
          t[cl][2][0] = 0.; t[cl][2][1] = 0.; t[cl][2][2] = -2.;
		  t[cl][0][0]=t[cl][0][0]/sqrt(3.);
		  t[cl][1][1]=t[cl][1][1]/sqrt(3.);
		  t[cl][2][2]=t[cl][2][2]/sqrt(3.);



          DeclU();
          
          for(int x = 0; x<ARRAY; x++)
          {
             for(int y = 0; y<ARRAY; y++)
             {
                U[0][x][y][0][0] = 1.; U[0][x][y][0][1] = 0.; U[0][x][y][0][2] = 0.;
                U[0][x][y][1][0] = 0.; U[0][x][y][1][1] = 1.; U[0][x][y][1][2] = 0.;
                U[0][x][y][2][0] = 0.; U[0][x][y][2][1] = 0.; U[0][x][y][2][2] = 1.;
                
                U[1][x][y][0][0] = 1.; U[1][x][y][0][1] = 0.; U[1][x][y][0][2] = 0.;
                U[1][x][y][1][0] = 0.; U[1][x][y][1][1] = 1.; U[1][x][y][1][2] = 0.;
                U[1][x][y][2][0] = 0.; U[1][x][y][2][1] = 0.; U[1][x][y][2][2] = 1.;
             }
          }

          cout << "Forming U for N = " << n << endl;
          for(int nuc=0; nuc<2; nuc++)
          {
          for(int eta=0; eta<etamax; eta++)
          {
             for(int x = 0; x<ARRAY; x++)
             {
                for(int y = 0; y<ARRAY; y++)
                { //at each point...

                   //Initializations- H as in Curtright paper U is Unitary Trasnformation
                   complex<double> H[3][3];
                   H[0][0] = 0.; H[0][1] = 0.; H[0][2] = 0.;
                   H[1][0] = 0.; H[1][1] = 0.; H[1][2] = 0.;
                   H[2][0] = 0.; H[2][1] = 0.; H[2][2] = 0.;

                   complex<double> pU[3][3]; //U for this x,y,eta
                   pU[0][0] = 0.; pU[0][1] = 0.; pU[0][2] = 0.;
                   pU[1][0] = 0.; pU[1][1] = 0.; pU[1][2] = 0.;
                   pU[2][0] = 0.; pU[2][1] = 0.; pU[2][2] = 0.;

                   double alpha = 0;


/*********************************************/
/*ARNING*******WARNING*******WARNING*******WA*/
/*********************************************/
/************COVPOT IS BEING SHIFTED**********/
/*********************************************/
/***********MAKE SURE THIS MAKES SENSE!*******/
/*********************************************/
/*********************************************/


                   for(int col = 0; col<8; col++) // alpha^2
                   {
                      alpha = alpha + CovPot[nuc][col][eta][x][y]*CovPot[nuc][col][eta][x][y];
                   }
                   
                   alpha = sqrt(alpha); //Magnitude of alpha

                   if(alpha != 0)
                   {                  
                   for(int a=0;a<3;a++)
                   {
                      for(int b=0;b<3;b++)
                      {
                         for(int col=0; col<8; col++)
                         { //Add all colors to get H
                            H[a][b] = H[a][b] + /*4.*/CovPot[nuc][col][eta][x][y]*t[col][a][b]/alpha;
                         }
                      }
                   }



                   //We have H for this position/eta next need theta/phi
                   double theta = (g/*/4.*/)*alpha/sqrt(etamax); //easy.

                   //phi
                   double DetH; //Determinant of H for phi; presuming DetH is real
                   DetH = real(H[0][0]*H[1][1]*H[2][2] +
                          H[0][1]*H[1][2]*H[2][0] +
                          H[0][2]*H[1][0]*H[2][1] -
                          H[0][2]*H[1][1]*H[2][0] -
                          H[0][0]*H[1][2]*H[2][1] -
                          H[0][1]*H[1][0]*H[2][2]); //Compiler needs to be told use only real part
                          //Likely DetH is currently imaginary :(
                          complex<double> DetH2; //Determinant of H for phi; presuming DetH is real
/*                   DetH2 = H[0][0]*H[1][1]*H[2][2] +
                          H[0][1]*H[1][2]*H[2][0] +
                          H[0][2]*H[1][0]*H[2][1] -
                          H[0][2]*H[1][1]*H[2][0] -
                          H[0][0]*H[1][2]*H[2][1] -
                          H[0][1]*H[1][0]*H[2][2];
*/
/*
                         if( (x==144 && y==215) || (x==54 && y==514) || (x==341 && y==300) || (x==98 && y==277) || (x==144 && y==215) || (x==300 && y==444) )
                         {
                          cout << "Here's  H, DetH= " << DetH << endl;
                          cout << H[0][0] << "   " << H[0][1] << "   " << H[0][2] << endl;
                          cout << H[1][0] << "   " << H[1][1] << "   " << H[1][2] << endl;
                          cout << H[2][0] << "   " << H[2][1] << "   " << H[2][2] << endl << endl;
                          }
/*                          double stop;
                          cin >> stop;
*/
                          
                   double acosarg=3.*DetH*sqrt(3.)/2.; //argument for acos
                   //Fit into proper domain (-1 to 1)              
/*     
                   if(acosarg>1.)
                   {
                      while(acosarg>1.)
                      {
                         acosarg = acosarg -1.;
                      }
                   }
                   if(acosarg<-1.)
                   {
                      while(acosarg<-1.)
                      {
                         acosarg = acosarg + 1.;
                      }
                   }
*/

                   if(abs(acosarg) > 1.0)
                   {
                      if(abs(acosarg) < 1.0000001)
                      {
                         if(acosarg < -1)
                         {
                            acosarg = -1.;
                         }
                         else
                         {
                            acosarg = 1.;
                         }
                      }
                      else
                      {
                         cout << "Warning! Acosarg too large!  " << acosarg << " x: " << x << " y: " << y << endl;
                      }
                      //double errval;
                      //cin >> errval;
                   } 
                   double phi = (acos(acosarg)-(PIE/2.))/3.;
                   
                   //Now Actually Exponentiate
                   //Sum over K
                   
                   for(int k=0; k<3; k++) // 0,1,2
                   {    
                      double cangle = phi+2.*PIE*k/3.; //angle for common factor
                      complex<double> cfactor; //factor for all k-th terms in sum
                      cfactor = exp((2./sqrt(3.))*I*theta*sin(cangle) );
                      cfactor = cfactor/(1 - 2*cos(2*cangle));

                      complex<double> Uk[3][3];
                      Uk[0][0] = 0.; Uk[0][1] = 0.; Uk[0][2] = 0.;
                      Uk[1][0] = 0.; Uk[1][1] = 0.; Uk[1][2] = 0.;
                      Uk[2][0] = 0.; Uk[2][1] = 0.; Uk[2][2] = 0.;
                      // Uk[a][b]
                      for(int a=0; a<3; a++)
                      {
                         for(int b=0; b<3; b++)
                         {
                            // H^2 term
                            for(int m=0; m<3; m++)
                            {
                               Uk[a][b]=Uk[a][b]+H[a][m]*H[m][b];
                            }

                            //H^1 Term
                            Uk[a][b]=Uk[a][b]+H[a][b]*2.*sin(cangle)/sqrt(3.);

                            
                            //H^0 Term
                            if(a==b)
                            {
                               Uk[a][b]=Uk[a][b]-(1.+2.*cos(2.*cangle))/3.;
                            }
                         }
                      }
                      //multiply by common factor for kth terms
                      for(int a=0; a<3; a++)
                      {
                         for(int b=0; b<3; b++)
                         {
                            pU[a][b]=pU[a][b]+Uk[a][b]*cfactor;
                         }
                      }
                   } // end kth sum

                   //Should be Unitary- give Det = 1 +0*i
                   complex<double>DetU = pU[0][0]*pU[1][1]*pU[2][2] +
                                         pU[0][1]*pU[1][2]*pU[2][0] +
                                         pU[0][2]*pU[1][0]*pU[2][1] -
                                         pU[0][2]*pU[1][1]*pU[2][0] -
                                         pU[0][0]*pU[1][2]*pU[2][1] -
                                         pU[0][1]*pU[1][0]*pU[2][2];
/*
                   cout << "x: " << x << " y: " << y <<endl;
                   cout << "Det U: " << DetU << endl << endl;

                   //Unitary Test
                   cout << "Uk unitarity test:" << endl;
*/

                   complex<double> UdU[3][3];
                   UdU[0][0] = 0.; UdU[0][1] = 0.; UdU[0][2] = 0.;
                   UdU[1][0] = 0.; UdU[1][1] = 0.; UdU[1][2] = 0.;
                   UdU[2][0] = 0.; UdU[2][1] = 0.; UdU[2][2] = 0.;
                   complex<double> UUd[3][3];
                   UUd[0][0] = 0.; UUd[0][1] = 0.; UUd[0][2] = 0.;
                   UUd[1][0] = 0.; UUd[1][1] = 0.; UUd[1][2] = 0.;
                   UUd[2][0] = 0.; UUd[2][1] = 0.; UUd[2][2] = 0.;
                   
                   for(int a=0; a<3; a++)
                      {
                         for(int b=0; b<3; b++)
                         {
                            for(int c=0; c<3; c++)
                            {
                               UdU[a][c]=UdU[a][c] + conj(pU[b][a])* pU[b][c];
                               UUd[a][c]=UUd[a][c] + pU[a][b]*conj(pU[c][b]);
                            }
                         }
                      }
                   complex<double>DetUdU = UdU[0][0]*UdU[1][1]*UdU[2][2] +
                                           UdU[0][1]*UdU[1][2]*UdU[2][0] +
                                           UdU[0][2]*UdU[1][0]*UdU[2][1] -
                                           UdU[0][2]*UdU[1][1]*UdU[2][0] -
                                           UdU[0][0]*UdU[1][2]*UdU[2][1] -
                                           UdU[0][1]*UdU[1][0]*UdU[2][2];
                                           
                   complex<double>DetUUd = UUd[0][0]*UUd[1][1]*UUd[2][2] +
                                           UUd[0][1]*UUd[1][2]*UUd[2][0] +
                                           UUd[0][2]*UUd[1][0]*UUd[2][1] -
                                           UUd[0][2]*UUd[1][1]*UUd[2][0] -
                                           UUd[0][0]*UUd[1][2]*UUd[2][1] -
                                           UUd[0][1]*UUd[1][0]*UUd[2][2];
                      
                      
                   //SEE! Is Unitary!

//                         if( (x==144 && y==215) || (x==54 && y==514) || (x==341 && y==300) || (x==98 && y==277) || (x==144 && y==215) || (x==300 && y==444) )
/*                         if((x==300 && y==444 && nuc==1 && col <7))
                         {
                   cout << "This test..." << endl;
                   cout << "At eta... " << eta << endl;
                   cout << "Not even UdU's final forme! Det is:  "  << DetUdU << endl;
                   cout << UdU[0][0] << "   " << UdU[0][1] << "   " << UdU[0][2] << endl;
                   cout << UdU[1][0] << "   " << UdU[1][1] << "   " << UdU[1][2] << endl;
                   cout << UdU[2][0] << "   " << UdU[2][1] << "   " << UdU[2][2] << endl << endl;
                   cout << "Not even UUd's final forme! Det is:  "  << DetUUd << endl;
                   cout << UUd[0][0] << "   " << UUd[0][1] << "   " << UUd[0][2] << endl;
                   cout << UUd[1][0] << "   " << UUd[1][1] << "   " << UUd[1][2] << endl;
                   cout << UUd[2][0] << "   " << UUd[2][1] << "   " << UUd[2][2] << endl << endl;

                   }
*/
                   complex<double> Uk[3][3];
                   Uk[0][0] = 0.; Uk[0][1] = 0.; Uk[0][2] = 0.;
                   Uk[1][0] = 0.; Uk[1][1] = 0.; Uk[1][2] = 0.;
                   Uk[2][0] = 0.; Uk[2][1] = 0.; Uk[2][2] = 0.;

                   //Product over Etas
                   for(int a=0; a<3; a++)
                      {
                         for(int b=0; b<3; b++)
                         {
                            for(int c=0; c<3; c++)
                            {//New U needs a temporary name...
                               //U after next Mult = sum of Old U *new U
                               Uk[a][c]=Uk[a][c] + U[nuc][x][y][a][b]*pU[b][c];
                            }
                         }
                      }
                   //Uk is new Uproduct
                   for(int a=0; a<3; a++)
                      {
                         for(int b=0; b<3; b++)
                         {
                            U[nuc][x][y][a][b] = Uk[a][b];
                         }
                      }

                   if(eta == etamax-1)
                   {
/*                      Numb++;
                      for(int a=0; a<3; a++)
                      {
                         for(int b=0; b<3; b++)
                         {
                            AvU[nuc][x][y][a][b] = AvU[nuc][x][y][a][b] + U[nuc][x][y][a][b];
                         }
                      }
*/                      

                      UdU[0][0] = 0.; UdU[0][1] = 0.; UdU[0][2] = 0.;
                      UdU[1][0] = 0.; UdU[1][1] = 0.; UdU[1][2] = 0.;
                      UdU[2][0] = 0.; UdU[2][1] = 0.; UdU[2][2] = 0.;
                      
                      UUd[0][0] = 0.; UUd[0][1] = 0.; UUd[0][2] = 0.;
                      UUd[1][0] = 0.; UUd[1][1] = 0.; UUd[1][2] = 0.;
                      UUd[2][0] = 0.; UUd[2][1] = 0.; UUd[2][2] = 0.;
                      
                      for(int a=0; a<3; a++)
                      {
                         for(int b=0; b<3; b++)
                         {
                            for(int c=0; c<3; c++)
                            {
                               UdU[a][c]=UdU[a][c] + conj(U[nuc][x][y][b][a])* U[nuc][x][y][b][c];
                               UUd[a][c]=UUd[a][c] + U[nuc][x][y][a][b]*conj(U[nuc][x][y][c][b]);
                            }
                         }
                      }
                      //SEE! Is Unitary!
                      
                   complex<double>DetUdU = UdU[0][0]*UdU[1][1]*UdU[2][2] +
                                           UdU[0][1]*UdU[1][2]*UdU[2][0] +
                                           UdU[0][2]*UdU[1][0]*UdU[2][1] -
                                           UdU[0][2]*UdU[1][1]*UdU[2][0] -
                                           UdU[0][0]*UdU[1][2]*UdU[2][1] -
                                           UdU[0][1]*UdU[1][0]*UdU[2][2];
                                           
                   complex<double>DetUUd = UUd[0][0]*UUd[1][1]*UUd[2][2] +
                                           UUd[0][1]*UUd[1][2]*UUd[2][0] +
                                           UUd[0][2]*UUd[1][0]*UUd[2][1] -
                                           UUd[0][2]*UUd[1][1]*UUd[2][0] -
                                           UUd[0][0]*UUd[1][2]*UUd[2][1] -
                                           UUd[0][1]*UUd[1][0]*UUd[2][2];

/*
                         if( (x==144 && y==215) || (x==341 && y==300)|| (x==144 && y==215) || (x==300 && y==444) )
                         {
                         cout << "FOR NUC = " << nuc << " at x " << x << " and y " << y << endl;
                      cout << "UdU's true final forme of final U! Det is: " << DetUdU << endl;
                      cout << UdU[0][0] << "   " << UdU[0][1] << "   " << UdU[0][2] << endl;
                      cout << UdU[1][0] << "   " << UdU[1][1] << "   " << UdU[1][2] << endl;
                      cout << UdU[2][0] << "   " << UdU[2][1] << "   " << UdU[2][2] << endl << endl;
                      cout << "UUd's true final forme of final U! Det is: " << DetUUd << endl;
                      cout << UUd[0][0] << "   " << UUd[0][1] << "   " << UUd[0][2] << endl;
                      cout << UUd[1][0] << "   " << UUd[1][1] << "   " << UUd[1][2] << endl;
                      cout << UUd[2][0] << "   " << UUd[2][1] << "   " << UUd[2][2] << endl << endl;
                      }
*/

/*                    double trashval;
                      cin >> trashval;                      
*/
                   }

                   } 
                }
             }//ends transverse loops
          }//ends eta loop
          }//end nuc loop



          RelCovPot();


          
          //FS Potential
/*
          ofstream ForwardOrder;
          ofstream ForwardOrderR;
          ofstream ForwardOrderI;
          ForwardOrder.open("FourthOrder.txt");
          ForwardOrderR.open("FourthOrderRe.txt");
          ForwardOrderI.open("FourthOrderIm.txt");
          ofstream ForwardOrderY;
          ofstream ForwardOrderRY;
          ofstream ForwardOrderIY;
          ForwardOrder.open("FourthOrderY.txt");
          ForwardOrderR.open("FourthOrderReY.txt");
          ForwardOrderI.open("FourthOrderImY.txt");
          
          ofstream BackwardOrder;
          ofstream SecOrder;
          ofstream FrtOrder;
          ofstream BackwardOrderR;
          ofstream SecOrderR;
          ofstream FrtOrderR;
          ofstream BackwardOrderI;
          ofstream SecOrderI;
          ofstream FrtOrderI;
          BackwardOrder.open("FirstOrderB.txt");
          SecOrder.open("SecondOrder.txt");
          FrtOrder.open("FourthOrder.txt");
          BackwardOrderR.open("FirstOrderBRe.txt");
          SecOrderR.open("SecondOrderRe.txt");
          FrtOrderR.open("FourthOrderRe.txt");
          BackwardOrderI.open("FirstOrderBIm.txt");
          SecOrderI.open("SecondOrderIm.txt");
          FrtOrderI.open("FourthOrderIm.txt");
          */

          DeclA();


          cout << "Calculation A for N = " << n << endl;
          for(int nuc = 0; nuc<2; nuc++)
          {
          
          for(int x = 2; x<ARRAY-2; x++)
             {
                for(int y = 2; y<ARRAY-2; y++)
                {
                   for(int a=0; a<3; a++)
                   {
                      for(int b=0; b<3; b++)
                      {
                         complex<double> Id=0.;//Identity Element
                          if(a==b)
                          {
                             Id=1.;
                          }
/*
                         //Forward Derivative
                         complex<double> UUdx; //Element of UUd that are relevant- "a, bth element"
                         UUdx = 0.;
                         
                         complex<double> UUdy;
                         UUdy = 0.;

                         //Backward Derivative
                         complex<double> UUdxB; //Element of UUd that are relevant- "a, bth element"
                         UUdxB = 0.;
                         
                         complex<double> UUdyB;
                         UUdyB = 0.;
*/
                         //Centered Derivative (2nd order in h)
                         complex<double> UUdx2; //Element of UUd that are relevant- "a, bth element"
                         UUdx2 = 0.;
                         
                         complex<double> UUdy2;
                         UUdy2 = 0.;

                         //Centered Derivative (2nd order in h)

                         complex<double> UUdx4[4]; //Element of UUd that are relevant- "a, bth element"
                         UUdx4[0] = 0.;
                         UUdx4[1] = 0.;
                         UUdx4[2] = 0.;
                         UUdx4[3] = 0.;
                         
                         complex<double> UUdy4[4];
                         UUdy4[0] = 0.;
                         UUdy4[1] = 0.;
                         UUdy4[2] = 0.;
                         UUdy4[3] = 0.;
                         
                         complex<double> UUdx4dagger; //Element of UUd that are relevant- "a, bth element"
                         UUdx4dagger = 0.;
                         
                         complex<double> UUdy4dagger;
                         UUdy4dagger = 0.;

                         for(int c=0; c<3; c++)
                         {
                         /*
                            //Forward Derivative
                            UUdx=UUdx + U[nuc][x][y][a][c]*conj(U[nuc][x+1][y][b][c]);
                            UUdy=UUdy + U[nuc][x][y][a][c]*conj(U[nuc][x][y+1][b][c]);
                         
                            //Backward
                            UUdxB=UUdxB + U[x][y][a][c]*conj(U[x-1][y][b][c]);
                            UUdyB=UUdyB + U[x][y][a][c]*conj(U[x][y-1][b][c]);
                         */   
                            //2nd order
                            UUdx2=UUdx2 - U[nuc][x][y][a][c]*conj(U[nuc][x-1][y][b][c])
                                        + U[nuc][x][y][a][c]*conj(U[nuc][x+1][y][b][c]);
                            UUdy2=UUdy2 - U[nuc][x][y][a][c]*conj(U[nuc][x][y-1][b][c])
                                        + U[nuc][x][y][a][c]*conj(U[nuc][x][y+1][b][c]);
                         
                            //4th order

                            UUdx4[0]=UUdx4[0] + U[nuc][x][y][a][c]*conj(U[nuc][x-2][y][b][c]);
                            UUdx4[1]=UUdx4[1] - U[nuc][x][y][a][c]*conj(U[nuc][x-1][y][b][c])*8.;
                            UUdx4[2]=UUdx4[2] + U[nuc][x][y][a][c]*conj(U[nuc][x+1][y][b][c])*8.;
                            UUdx4[3]=UUdx4[3] - U[nuc][x][y][a][c]*conj(U[nuc][x+2][y][b][c]);
                            UUdy4[0]=UUdy4[0] + U[nuc][x][y][a][c]*conj(U[nuc][x][y-2][b][c]);
                            UUdy4[1]=UUdy4[1] - U[nuc][x][y][a][c]*conj(U[nuc][x][y-1][b][c])*8.;
                            UUdy4[2]=UUdy4[2] + U[nuc][x][y][a][c]*conj(U[nuc][x][y+1][b][c])*8.;
                            UUdy4[3]=UUdy4[3] - U[nuc][x][y][a][c]*conj(U[nuc][x][y+2][b][c]);
                                        
                                        
                            
                            
                            UUdx4dagger=UUdx4dagger + U[nuc][x-2][y][a][c]*conj(U[nuc][x][y][b][c])
                                                    - U[nuc][x-1][y][a][c]*conj(U[nuc][x][y][b][c])*8.
                                                    + U[nuc][x+1][y][a][c]*conj(U[nuc][x][y][b][c])*8.
                                                    - U[nuc][x+2][y][a][c]*conj(U[nuc][x][y][b][c]);
                            UUdy4dagger=UUdy4dagger + U[nuc][x][y-2][a][c]*conj(U[nuc][x][y][b][c])
                                                    - U[nuc][x][y-1][a][c]*conj(U[nuc][x][y][b][c])*8.
                                                    + U[nuc][x][y+1][a][c]*conj(U[nuc][x][y][b][c])*8.
                                                    - U[nuc][x][y+2][a][c]*conj(U[nuc][x][y][b][c]);
                         
                         }

                          A[nuc][0][x][y][a][b] = -I*(UUdx4[0]+UUdx4[1]+UUdx4[2]+UUdx4[3])/(g*12*h);
                          A[nuc][1][x][y][a][b] = -I*(UUdy4[0]+UUdy4[1]+UUdy4[2]+UUdy4[3])/(g*12*h);
                          
//                          A[nuc][0][x][y][a][b] = -I*(UUdx2)/(g*2*h);
//                          A[nuc][1][x][y][a][b] = -I*(UUdy2)/(g*2*h);
                          
                         /*
                          A[nuc][0][x][y][a][b] = -I*(Id-UUdx)/(g*h);
                          A[nuc][1][x][y][a][b] = -I*(Id-UUdy)/(g*h);
                          
                          AB[0][x][y][a][b] = I*(Id-UUdxB)/(g*h);
                          AB[1][x][y][a][b] = I*(Id-UUdyB)/(g*h);
                          
                          A2nd[0][x][y][a][b] = I*(UUdx2)/(g*2*h);
                          A2nd[1][x][y][a][b] = I*(UUdy2)/(g*2*h);
                          
                          A4th[0][x][y][a][b] = I*(UUdx4)/(g*12*h);
                          A4th[1][x][y][a][b] = I*(UUdy4)/(g*12*h);
                          */

//endl at y-loop

                         if(a==b && a==2)
                         {

/*                         if( (x==144 && y==215) || (x==54 && y==514) || (x==341 && y==300) || (x==98 && y==277) || (x==144 && y==215) || (x==300 && y==444) )
                         {
                          cout << "U,x,y:" << endl;
                          cout << U[nuc][x][y][0][0] << "   " << U[nuc][x][y][0][1] << "   " << U[nuc][x][y][0][2] << endl;
                          cout << U[nuc][x][y][1][0] << "   " << U[nuc][x][y][1][1] << "   " << U[nuc][x][y][1][2] << endl;
                          cout << U[nuc][x][y][2][0] << "   " << U[nuc][x][y][2][1] << "   " << U[nuc][x][y][2][2] << endl << endl;
                          
                          cout << "U,x+1,y:" << endl;
                          cout << U[nuc][x+1][y][0][0] << "   " << U[nuc][x+1][y][0][1] << "   " << U[nuc][x+1][y][0][2] << endl;
                          cout << U[nuc][x+1][y][1][0] << "   " << U[nuc][x+1][y][1][1] << "   " << U[nuc][x+1][y][1][2] << endl;
                          cout << U[nuc][x+1][y][2][0] << "   " << U[nuc][x+1][y][2][1] << "   " << U[nuc][x+1][y][2][2] << endl << endl;
                          
                          cout << "A^1:" << endl;
                          cout << A[nuc][0][x][y][0][0] << "   " << A[nuc][0][x][y][0][1] << "   " << A[nuc][0][x][y][0][2] << endl;
                          cout << A[nuc][0][x][y][1][0] << "   " << A[nuc][0][x][y][1][1] << "   " << A[nuc][0][x][y][1][2] << endl;
                          cout << A[nuc][0][x][y][2][0] << "   " << A[nuc][0][x][y][2][1] << "   " << A[nuc][0][x][y][2][2] << endl << endl;
                          }                          

               //           cout << "Color Elements:" << endl;
*/


//XYSYM
/*
if(nuc==0 && x==144 && y==215)
{
          cout << "A print " << endl;
          cout << A[0][0][x][y][0][0] << "   " << A[0][0][x][y][0][1] << "   " << A[0][0][x][y][0][2] << endl;
          cout << A[0][0][x][y][1][0] << "   " << A[0][0][x][y][1][1] << "   " << A[0][0][x][y][1][2] << endl;
          cout << A[0][0][x][y][2][0] << "   " << A[0][0][x][y][2][1] << "   " << A[0][0][x][y][2][2] << endl << endl;

          for(int q=0;q<3;q++)
          {
          for(int r=0;r<3;r++)
          {
             ATest[0][0][q][r] = A[0][0][x][y][q][r];
          }
          }
}
if(nuc==0 && y==144 && x==215)
{
          for(int q=0;q<3;q++)
          {
          for(int r=0;r<3;r++)
          {
             ATest[0][1][q][r] = A[0][1][x][y][q][r];
          }
          }
}



if(nuc==0 && x==54 && y==514)
{
          cout << "A print " << endl;
          cout << A[0][0][x][y][0][0] << "   " << A[0][0][x][y][0][1] << "   " << A[0][0][x][y][0][2] << endl;
          cout << A[0][0][x][y][1][0] << "   " << A[0][0][x][y][1][1] << "   " << A[0][0][x][y][1][2] << endl;
          cout << A[0][0][x][y][2][0] << "   " << A[0][0][x][y][2][1] << "   " << A[0][0][x][y][2][2] << endl << endl;

          for(int q=0;q<3;q++)
          {
          for(int r=0;r<3;r++)
          {
             ATest[1][0][q][r] = A[0][0][x][y][q][r];
          }
          }
}
if(nuc==0 && y==54 && x==514)
{
          for(int q=0;q<3;q++)
          {
          for(int r=0;r<3;r++)
          {
             ATest[1][1][q][r] = A[0][1][x][y][q][r];
          }
          }
}

if(nuc==0 && x==341 && y==300)
{
          cout << "A print " << endl;
         cout << A[0][0][x][y][0][0] << "   " << A[0][0][x][y][0][1] << "   " << A[0][0][x][y][0][2] << endl;
          cout << A[0][0][x][y][1][0] << "   " << A[0][0][x][y][1][1] << "   " << A[0][0][x][y][1][2] << endl;
          cout << A[0][0][x][y][2][0] << "   " << A[0][0][x][y][2][1] << "   " << A[0][0][x][y][2][2] << endl << endl;

          for(int q=0;q<3;q++)
          {
          for(int r=0;r<3;r++)
          {
             ATest[2][0][q][r] = A[0][0][x][y][q][r];
          }
          }
}
if(nuc==0 && y==341 && x==300)
{
          for(int q=0;q<3;q++)
          {
          for(int r=0;r<3;r++)
          {
             ATest[2][1][q][r] = A[0][1][x][y][q][r];
          }
          }
}

if(nuc==0 && x==98 && y==277)
{
          cout << "A print " << endl;
         cout << A[0][0][x][y][0][0] << "   " << A[0][0][x][y][0][1] << "   " << A[0][0][x][y][0][2] << endl;
          cout << A[0][0][x][y][1][0] << "   " << A[0][0][x][y][1][1] << "   " << A[0][0][x][y][1][2] << endl;
          cout << A[0][0][x][y][2][0] << "   " << A[0][0][x][y][2][1] << "   " << A[0][0][x][y][2][2] << endl << endl;

          for(int q=0;q<3;q++)
          {
          for(int r=0;r<3;r++)
          {
             ATest[3][0][q][r] = A[0][0][x][y][q][r];
          }
          }
}
if(nuc==0 && y==98 && x==277)
{
          for(int q=0;q<3;q++)
          {
          for(int r=0;r<3;r++)
          {
             ATest[3][1][q][r] = A[0][1][x][y][q][r];
          }
          }
}

if(nuc==0 && x==144 && y==215)
{
          cout << "A print " << endl;
         cout << A[0][0][x][y][0][0] << "   " << A[0][0][x][y][0][1] << "   " << A[0][0][x][y][0][2] << endl;
          cout << A[0][0][x][y][1][0] << "   " << A[0][0][x][y][1][1] << "   " << A[0][0][x][y][1][2] << endl;
          cout << A[0][0][x][y][2][0] << "   " << A[0][0][x][y][2][1] << "   " << A[0][0][x][y][2][2] << endl << endl;

          for(int q=0;q<3;q++)
          {
          for(int r=0;r<3;r++)
          {
             ATest[4][0][q][r] = A[0][0][x][y][q][r];
          }
          }
}
if(nuc==0 && y==144 && x==215)
{
          for(int q=0;q<3;q++)
          {
          for(int r=0;r<3;r++)
          {
             ATest[4][1][q][r] = A[0][1][x][y][q][r];
          }
          }
}

if(nuc==0 && x==300 && y==454)
{
          cout << "A print " << endl;
          cout << A[0][0][x][y][0][0] << "   " << A[0][0][x][y][0][1] << "   " << A[0][0][x][y][0][2] << endl;
          cout << A[0][0][x][y][1][0] << "   " << A[0][0][x][y][1][1] << "   " << A[0][0][x][y][1][2] << endl;
          cout << A[0][0][x][y][2][0] << "   " << A[0][0][x][y][2][1] << "   " << A[0][0][x][y][2][2] << endl << endl;

          for(int q=0;q<3;q++)
          {
          for(int r=0;r<3;r++)
          {
             ATest[5][0][q][r] = A[0][0][x][y][q][r];
          }
          }
}
if(nuc==0 && y==300 && x==454)
{
          for(int q=0;q<3;q++)
          {
          for(int r=0;r<3;r++)
          {
             ATest[5][1][q][r] = A[0][1][x][y][q][r];
          }
          }
}
*/

/*
                         if(x==256 && y==332)
                         {
                            if(a==2&&b==2 && nuc==0)
                            {
                            cout << "Pre-AFix:" << endl;
                            cout << "A^x_1:" << endl;
                            for(int aa=0; aa<3; aa++)
                            {
                               for(int bb=0; bb<3; bb++)
                               {
                                  cout << A[0][0][x][y][aa][bb] << "  ";
                               }
                               cout << endl;
                            }

                            cout << "A^y_1:" << endl;
                            for(int aa=0; aa<3; aa++)
                            {
                               for(int bb=0; bb<3; bb++)
                               {
                                  cout << A[0][1][x][y][aa][bb] << "  ";
                               }
                               cout << endl;
                            }
                            }
                            

                            if(a==2&&b==2 && nuc==1)
                            {
                            cout << "Pre-AFix:" << endl;
                            
                            cout << "A^x_2:" << endl;
                            for(int aa=0; aa<3; aa++)
                            {
                               for(int bb=0; bb<3; bb++)
                               {
                                  cout << A[1][0][x][y][aa][bb] << "  ";
                               }
                               cout << endl;
                            }
                            cout << "A^y_2:" << endl;
                            for(int aa=0; aa<3; aa++)
                            {
                               for(int bb=0; bb<3; bb++)
                               {
                                  cout << A[1][1][x][y][aa][bb] << "  ";
                               }
                               cout << endl;
                            }
                            }



                         }
*/




                          complex<double> color[2][8];
                          for(int col=0; col<8; col++)
                          {
                             for(int dir=0; dir<2; dir++)
                             {
                                color[dir][col]=0.;
                             }
                          }
                          for(int col=0; col<8; col++)
                          {
                          for(int dir=0; dir<2; dir++)
                          {
                             /*
                             complex<double> color=0.;
                             complex<double> colorb=0.;
                             complex<double> color2=0.;
                             complex<double> color4=0.;
                             */


                             for(int c=0; c<3; c++)
                             {
                                color[dir][col]= color[dir][col] + (t[col][0][c]*A[nuc][dir][x][y][c][0] + t[col][1][c]*A[nuc][dir][x][y][c][1] + t[col][2][c]*A[nuc][dir][x][y][c][2])/2.;
                                
                                
                                /*
                                colorb= colorb + (t[col][0][c]*AB[0][x][y][c][0] + t[col][1][c]*AB[0][x][y][c][1] + t[col][2][c]*AB[0][x][y][c][2])/2.;
                                
                                color2= color2 + (t[col][0][c]*A2nd[0][x][y][c][0] + t[col][1][c]*A2nd[0][x][y][c][1] + t[col][2][c]*A2nd[0][x][y][c][2])/2.;
                                
                                color4= color4 + (t[col][0][c]*A4th[0][x][y][c][0] + t[col][1][c]*A4th[0][x][y][c][1] + t[col][2][c]*A4th[0][x][y][c][2])/2.;
                                */

                             }
                             
                             
                             
                             /**Faking color**/
                             
                             bool FakeCol = false;
                             
                             if(FakeCol)
                             {
                             if(col == nuc && dir == 0)
                             {
                   //             cout << "Used color: " << col << endl;
                                color[dir][col] = 10.;
                             }
                             else
                             {
                                if(nuc == 1 && dir == 1 && col == 1)
                                {
                      //             cout << "Used color: " << col << " in second chance for nucleus: " << nuc << " in direction: " << dir << endl;
                                   color[dir][col] = 10.;
                                }
                                else
                                {
                                   color[dir][col] = 0.;
                                }
                             }
                             }
                             
                             
                             //Alternate Method
                             FakeCol = false;
                             
                             if(FakeCol)
                             {
                                if(nuc == 0 && dir == 0)
                                {
                                   color[dir][col] = (col+1)/2.;
                                }

                                if(nuc == 1 && dir == 0)
                                {
                                   color[dir][col] = (3+col%2)/sqrt(3.);
                                }

                                if(nuc == 0 && dir == 1)
                                {
                                   if(col%3 == 0)
                                   {
                                      color[dir][col] = 5/sqrt(5.);
                                   }
                                   if(col%3 == 1)
                                   {
                                      color[dir][col] = 9/sqrt(5.);
                                   }
                                   if(col%3 == 2)
                                   {
                                      color[dir][col] = 2/sqrt(5.);
                                   }
                                }

                                if(nuc == 1 && dir == 1)
                                {
                                   color[dir][col] = (8-col)/5.;
                                }
                                
                                if(nuc == 1)
                                {color[dir][col] = 0.;}

                             }

                             
                             
                             
                             
                             
                             /*
                             ForwardOrder << real(color[0][col]) << "  " << imag(color[0][col]) << "  ";
                             ForwardOrderR << real(color[0][col]) << "  ";
                             ForwardOrderI << imag(color[0][col]) << "  ";

                             ForwardOrderY << real(color[1][col]) << "  " << imag(color[1][col]) << "  ";
                             ForwardOrderRY << real(color[1][col]) << "  ";
                             ForwardOrderIY << imag(color[1][col]) << "  ";
                             
                             cout << col << " 4th order  forward: " << color[dir][col] << endl;
                             cout << col << " 1st order backward: " << colorb << endl;
                             cout << col << " 2nd order centered: " << color2 << endl;
                             cout << col << " 4th order centered: " << color4 << endl;
                             
                             
                             BackwardOrder << real(colorb) << "  " << imag(colorb) << "  ";
                             SecOrder << real(color2) << "  " << imag(color2) << "  ";
                             FrtOrder << real(color4) << "  " << imag(color4) << "  ";
                             
                             
                             
                             BackwardOrderR << real(colorb) << "  ";
                             SecOrderR << real(color2) << "  ";
                             FrtOrderR << real(color4) << "  ";
                             
                             
                             BackwardOrderI << imag(colorb) << "  ";
                             SecOrderI << imag(color2) << "  ";
                             FrtOrderI << imag(color4) << "  ";
                      */

                          }
                          }
                          
                          
/*                          
                          cout << "Color Check, get real..." << endl;
                          cout << "For Ax as: " << endl;
                          cout << A[nuc][0][x][y][0][0] << "  " << A[nuc][0][x][y][0][1] << "  " << A[nuc][0][x][y][0][2] << endl;
                          cout << A[nuc][0][x][y][1][0] << "  " << A[nuc][0][x][y][1][1] << "  " << A[nuc][0][x][y][1][2] << endl;
                          cout << A[nuc][0][x][y][2][0] << "  " << A[nuc][0][x][y][2][1] << "  " << A[nuc][0][x][y][2][2] << endl;
                          cout << "Calculated Conjugate: " << endl;
                          cout << Ad[nuc][0][x][y][0][0] << "  " << Ad[nuc][0][x][y][0][1] << "  " << Ad[nuc][0][x][y][0][2] << endl;
                          cout << Ad[nuc][0][x][y][1][0] << "  " << Ad[nuc][0][x][y][1][1] << "  " << Ad[nuc][0][x][y][1][2] << endl;
                          cout << Ad[nuc][0][x][y][2][0] << "  " << Ad[nuc][0][x][y][2][1] << "  " << Ad[nuc][0][x][y][2][2] << endl;
                          cout << "Individual Conj of Calculated Conjugate (test of conj() function)" << endl;
                          cout << conj(Ad[nuc][0][x][y][0][0]) << "  " << conj(Ad[nuc][0][x][y][0][1]) << "  " << conj(Ad[nuc][0][x][y][0][2]) << endl;
                          cout << conj(Ad[nuc][0][x][y][1][0]) << "  " << conj(Ad[nuc][0][x][y][1][1]) << "  " << conj(Ad[nuc][0][x][y][1][2]) << endl;
                          cout << conj(Ad[nuc][0][x][y][2][0]) << "  " << conj(Ad[nuc][0][x][y][2][1]) << "  " << conj(Ad[nuc][0][x][y][2][2]) << endl;
                          cout << "gets colors" << endl;
                          cout << color[0][0] << "  " << color[0][1] << "  " << color[0][2] << "  " << color[0][3] << endl;
                          cout << color[0][4] << "  " << color[0][5] << "  " << color[0][6] << "  " << color[0][7] << endl;
                          cout << "For Ay as: " << endl;
                          cout << A[nuc][1][x][y][0][0] << "  " << A[nuc][1][x][y][0][1] << "  " << A[nuc][1][x][y][0][2] << endl;
                          cout << A[nuc][1][x][y][1][0] << "  " << A[nuc][1][x][y][1][1] << "  " << A[nuc][1][x][y][1][2] << endl;
                          cout << A[nuc][1][x][y][2][0] << "  " << A[nuc][1][x][y][2][1] << "  " << A[nuc][1][x][y][2][2] << endl;
                          cout << "Calculated Conjugate: " << endl;
                          cout << Ad[nuc][1][x][y][0][0] << "  " << Ad[nuc][1][x][y][0][1] << "  " << Ad[nuc][1][x][y][0][2] << endl;
                          cout << Ad[nuc][1][x][y][1][0] << "  " << Ad[nuc][1][x][y][1][1] << "  " << Ad[nuc][1][x][y][1][2] << endl;
                          cout << Ad[nuc][1][x][y][2][0] << "  " << Ad[nuc][1][x][y][2][1] << "  " << Ad[nuc][1][x][y][2][2] << endl;
                          cout << "gets colors" << endl;
                          cout << color[1][0] << "  " << color[1][1] << "  " << color[1][2] << "  " << color[1][3] << endl;
                          cout << color[1][4] << "  " << color[1][5] << "  " << color[1][6] << "  " << color[1][7] << endl << endl;
*/
                          
                          //Forget old Transformation...

                          for(int q=0; q<3; q++)
                          {
                             for(int r=0; r<3; r++)
                             {
                                A[nuc][0][x][y][q][r] = 0;
                                A[nuc][1][x][y][q][r] = 0;
                             }
                          }
                          //Remake with real part of color calculated
                          for(int q=0; q<3; q++)
                          {
                             for(int r=0; r<3; r++)
                             {
                                for(int col=0; col<8; col++)
                                {
                                    A[nuc][0][x][y][q][r] = A[nuc][0][x][y][q][r] +
                                                            real(color[0][col])*t[col][q][r];
                                    A[nuc][1][x][y][q][r] = A[nuc][1][x][y][q][r] +
                                                            real(color[1][col])*t[col][q][r];
/*
if(nuc==0 && x==144 && y==215 && nuc==0 && q==0 && r==2)
{
          cout << "A test PRIME for " << col << endl;
          cout << "x: " << x << " y: " << y << " a= " << q << " b= " << r << endl;
          cout <<  real(color[0][col]) << " vs " << imag(color[0][col]) << endl;
}

if(nuc==0 && x==54 && y==514 && nuc==0 && q==1 && r==2)
{
          cout << "A test PRIME for " << col << endl;
          cout << "x: " << x << " y: " << y << " a= " << q << " b= " << r << endl;
          cout << real(color[0][col]) << " vs " << imag(color[0][col]) << endl;
}
if(nuc==0 && x==341 && y==300 && nuc==0 && q==2 && r==2)
{
          cout << "A test PRIME for " << col << endl;
          cout << "x: " << x << " y: " << y << " a= " << q << " b= " << r << endl;
          cout << real(color[0][col]) << " vs " << imag(color[0][col]) << endl;
}
if(nuc==0 && x==98 && y==277 && nuc==0 && q==0 && r==2)
{
          cout << "A test PRIME for " << col << endl;
          cout << "x: " << x << " y: " << y << " a= " << q << " b= " << r << endl;
          cout << real(color[0][col]) << " vs " << imag(color[0][col]) << endl;
}
if(nuc==0 && x==144 && y==215 && nuc==0 && q==1 && r==2)
{
          cout << "A test PRIME for " << col << endl;
          cout << "x: " << x << " y: " << y << " a= " << q << " b= " << r << endl;
          cout << real(color[0][col]) << " vs " << imag(color[0][col]) << endl;
}
if(nuc==0 && x==300 && y==444 && nuc==0 && q==2 && r==2)
{
          cout << "A test PRIME for " << col << endl;
          cout << "x: " << x << " y: " << y << " a= " << q << " b= " << r << endl;
          cout << real(color[0][col]) << " vs " << imag(color[0][col]) << endl << endl;
}
*/

                                }
                             }
                          }

                          /*
                          ForwardOrder << endl;
                          BackwardOrder << endl;
                          SecOrder << endl;
                          FrtOrder << endl;
                          
                          ForwardOrderR << endl;
                          BackwardOrderR << endl;
                          SecOrderR << endl;
                          FrtOrderR << endl;
                          
                          ForwardOrderI << endl;
                          BackwardOrderI << endl;
                          SecOrderI << endl;
                          FrtOrderI << endl;
                          */

                          //cout << endl;
                          }
                      }
                   }
                }// end y loop


             }
             
/*********/
/*CLEAN B*/
/*********/
if(false)
{
             for(int x = 2; x<ARRAY-2; x++)
             {
                for(int y = 2; y<ARRAY-2; y++)
                {
                   static complex<double> di[3][3][2];
                   for(int a=0; a<3; a++)
                   {
                      for(int b=0; b<3; b++)
                      {
                         for(int c=0; c<2; c++)
                         {
                            di[a][b][c] = 0.;
                         }
                      }
                   }
                   
                   for(int a=0; a<3; a++)
                   {
                      for(int b=0; b<3; b++)
                      {
                         di[a][b][0] = 
                           A[nuc][0][x][y-2][a][b]
                         - A[nuc][0][x][y-1][a][b]*8
                         + A[nuc][0][x][y+1][a][b]*8.
                         - A[nuc][0][x][y+2][a][b]

                         - A[nuc][1][x-2][y][a][b]
                         + A[nuc][1][x-1][y][a][b]*8.
                         - A[nuc][1][x+1][y][a][b]*8.
                         + A[nuc][1][x+2][y][a][b];
                         
                         di[a][b][0] = -di[a][b][0]/(12.*h);

                         for(int c=0; c<3; c++)
                         {
                            di[a][b][1] = di[a][b][1]  -I*g*(A[nuc][1][x][y][a][c]*A[nuc][0][x][y][c][b]
                                                       - A[nuc][0][x][y][a][c]*A[nuc][1][x][y][c][b]);
                         }
                      }
                   }

                   complex<double> color[2][8];

                   for(int col=0; col<8; col++)
                   {
                      for(int dir=0; dir<2; dir++)
                      {
                         color[dir][col]=0.;
                      }
                   }

                   for(int col=0; col<8; col++)
                   {
                      for(int c=0; c<3; c++)
                      {
                         for(int dir=0; dir<2; dir++)
                         {
                            color[dir][col]= color[dir][col] + (t[col][0][c]*di[c][0][dir] + t[col][1][c]*di[c][1][dir] + t[col][2][c]*di[c][2][dir])/2.;
                         }
                      }
                      // AScale[col][x][y] = -color[0][col]/color[1][col];
                   }
                }
             }
 /***AScale Now Generated for all (x,y) of this nuc, remake A yet again:***/

             for(int x = 2; x<ARRAY-2; x++)
             {
                for(int y = 2; y<ARRAY-2; y++)
                {
                   complex<double> color[2][8];
                   for(int col=0; col<8; col++)
                   {
                      for(int dir=0; dir<2; dir++)
                      {
                         color[dir][col]=0.;
                      }
                   }
                   for(int col=0; col<8; col++)
                   {
                      for(int dir=0; dir<2; dir++)
                      {
                         for(int c=0; c<3; c++)
                         {
                            color[dir][col]= color[dir][col] + (t[col][0][c]*A[nuc][dir][x][y][c][0] + t[col][1][c]*A[nuc][dir][x][y][c][1] + t[col][2][c]*A[nuc][dir][x][y][c][2])/2.;
                         }
                      }
                   }

                   //Forget old Transformation...

                   for(int q=0; q<3; q++)
                   {
                      for(int r=0; r<3; r++)
                      {
                         A[nuc][0][x][y][q][r] = 0;
                         A[nuc][1][x][y][q][r] = 0;
                      }
                   }

                   //Remake with real part of color calculated
                   for(int q=0; q<3; q++)
                   {
                      for(int r=0; r<3; r++)
                      {
                         for(int col=0; col<8; col++)
                         {
                            A[nuc][0][x][y][q][r] = A[nuc][0][x][y][q][r] +
                                                    real(color[0][col])*t[col][q][r];
                            A[nuc][1][x][y][q][r] = A[nuc][1][x][y][q][r] +
                                                    real(color[1][col])*t[col][q][r];
                         }
                      }
                   }
                   /*
                   if(x==256 && y==332)
                   {
                      cout << "Color Fix Print: " << endl;
                         for(int cl=0; cl<8;cl++)
                         {
                            cout << color[0][cl] << "  ";
                         }
                      cout << endl << endl;
                      
                      cout << "Color Scale Print: " << endl;
                         for(int cl=0; cl<8;cl++)
                         {
                            cout << AScale[cl][x][y] << "  ";
                         }
                      cout << endl << endl;
                   }
                   */
                }
             }
}
/***********/
/*CLEANED B*/
/***********/


/*
             int xx=256; int yy=332;
             
             complex<double> di[4][3][3];

             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   //d2 U
                   di[0][a][b] = U[nuc][xx][yy-2][a][b]
                               - U[nuc][xx][yy-1][a][b]*8
                               + U[nuc][xx][yy+1][a][b]*8.
                               - U[nuc][xx][yy+2][a][b];
                   //d1 U
                   di[1][a][b] = U[nuc][xx-2][yy][a][b]
                               - U[nuc][xx-1][yy][a][b]*8
                               + U[nuc][xx+1][yy][a][b]*8.
                               - U[nuc][xx+2][yy][a][b];
                   //d2 Udag
                   di[2][a][b] = conj(U[nuc][xx][yy-2][b][a])
                               - conj(U[nuc][xx][yy-1][b][a])*8
                               + conj(U[nuc][xx][yy+1][b][a])*8.
                               - conj(U[nuc][xx][yy+2][b][a]);
                   //d1 Udag
                   di[3][a][b] = conj(U[nuc][xx-2][yy][b][a])
                               - conj(U[nuc][xx-1][yy][b][a])*8
                               + conj(U[nuc][xx+1][yy][b][a])*8.
                               - conj(U[nuc][xx+2][yy][b][a]);
                }
             }

/*
             complex<double> UTest[3][3];
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   UTest[a][b] = 0.;
                }
             }

             cout << "UTest" << endl << endl;
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   for(int c = 0; c<3; c++)
                   {
                      //                                 d2 U d1 Udag           d1 U d2 Udag
                      UTest[a][b] = UTest[a][b] - I*(di[0][a][c]*di[3][c][b] - di[1][a][c]*di[2][c][b])/(144.*g*h*h);
                   }
                   cout << UTest[a][b] << "  ";
                }
                cout << endl;
             }
             cout << endl;

/**********************************************/
/*
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   //d2 U
                   di[0][a][b] = - U[nuc][xx][yy-1][a][b]
                                 + U[nuc][xx][yy+1][a][b];
                   //d1 U
                   di[1][a][b] = - U[nuc][xx-1][yy][a][b]
                                 + U[nuc][xx+1][yy][a][b];
                   //d2 Udag
                   di[2][a][b] = - conj(U[nuc][xx][yy-1][b][a])
                                 + conj(U[nuc][xx][yy+1][b][a]);
                   //d1 Udag
                   di[3][a][b] = - conj(U[nuc][xx-1][yy][b][a])
                                 + conj(U[nuc][xx+1][yy][b][a]);
                }
             }

             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   UTest[a][b] = 0.;
                }
             }
                                                                                                                                        
             cout << "UTest 1st order acc7" << endl << endl;
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   for(int c = 0; c<3; c++)
                   {
                      //                                 d2 U d1 Udag           d1 U d2 Udag
                      UTest[a][b] = UTest[a][b] + I*(di[0][a][c]*di[3][c][b] - di[1][a][c]*di[2][c][b])/(4.*h*h*g);
                   }
                   cout << UTest[a][b] << "  ";
                }
                cout << endl;
             }
             cout << endl;
             
/**********************************************/
/*
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   di[0][a][b] = 0.;
                   di[1][a][b] = 0.;
                }
             }

             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   for(int c=0; c<3; c++)
                   {
                      //d2 UUdag
                      di[0][a][b] = di[0][a][b] + U[nuc][xx][yy-2][a][c]*conj(U[nuc][xx][yy-2][b][c])
                                  - U[nuc][xx][yy-1][a][c]*conj(U[nuc][xx][yy-1][b][c])*8
                                  + U[nuc][xx][yy+1][a][c]*conj(U[nuc][xx][yy+1][b][c])*8.
                                  - U[nuc][xx][yy+2][a][c]*conj(U[nuc][xx][yy+2][b][c]);
                      //d1 UUdag
                      di[1][a][b] = di[1][a][b] + U[nuc][xx-2][yy][a][c]*conj(U[nuc][xx-2][yy][b][c])
                                  - U[nuc][xx-1][yy][a][c]*conj(U[nuc][xx-1][yy][b][c])*8
                                  + U[nuc][xx+1][yy][a][c]*conj(U[nuc][xx+1][yy][b][c])*8.
                                  - U[nuc][xx+2][yy][a][c]*conj(U[nuc][xx+2][yy][b][c]);
                   }
                }
             }

             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   UTest[a][b] = 0.;
                }
             }

             cout << "d2UUdag Test" << endl << endl;
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   cout << di[0][a][b] << "  ";
                }
                cout << endl;
             }
             cout << endl;

             cout << "d1UUdag Test" << endl << endl;
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   cout << di[1][a][b] << "  ";
                }
                cout << endl;
             }
             cout << endl;

/**********************************************/
/*
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   //d1d2 U
                   di[0][a][b] = conj(U[nuc][xx-2][yy-2][b][a] - U[nuc][xx-2][yy-1][b][a]*8 + U[nuc][xx-2][yy+1][b][a]*8. - U[nuc][xx-2][yy+2][b][a])
                               - conj(U[nuc][xx-1][yy-2][b][a] - U[nuc][xx-1][yy-1][b][a]*8 + U[nuc][xx-1][yy+1][b][a]*8. - U[nuc][xx-1][yy+2][b][a])*8.
                               + conj(U[nuc][xx+1][yy-2][b][a] - U[nuc][xx+1][yy-1][b][a]*8 + U[nuc][xx+1][yy+1][b][a]*8. - U[nuc][xx+1][yy+2][b][a])*8.
                               - conj(U[nuc][xx+2][yy-2][b][a] - U[nuc][xx+2][yy-1][b][a]*8 + U[nuc][xx+2][yy+1][b][a]*8. - U[nuc][xx+2][yy+2][b][a]);
                   //d2d1 U
                   di[1][a][b] = conj(U[nuc][xx-2][yy-2][b][a] - U[nuc][xx-1][yy-2][b][a]*8 + U[nuc][xx+1][yy-2][b][a]*8. - U[nuc][xx+2][yy-2][b][a])
                               - conj(U[nuc][xx-2][yy-1][b][a] - U[nuc][xx-1][yy-1][b][a]*8 + U[nuc][xx+1][yy-1][b][a]*8. - U[nuc][xx+2][yy-1][b][a])*8.
                               + conj(U[nuc][xx-2][yy+1][b][a] - U[nuc][xx-1][yy+1][b][a]*8 + U[nuc][xx+1][yy+1][b][a]*8. - U[nuc][xx+2][yy+1][b][a])*8.
                               - conj(U[nuc][xx-2][yy+2][b][a] - U[nuc][xx-1][yy+2][b][a]*8 + U[nuc][xx+1][yy+2][b][a]*8. - U[nuc][xx+2][yy+2][b][a]);
                }
             }

             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   UTest[a][b] = 0.;
                }
             }

             cout << "Ud1d2Udag" << endl << endl;
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   for(int c = 0; c<3; c++)
                   {
                      //                                 U d1d2 Udag
                      UTest[a][b] = UTest[a][b] - (U[nuc][xx][yy][a][c]*di[0][c][b])/(144.*g*h*h);
                   }
                   cout << UTest[a][b] << "  ";
                }
                cout << endl;
             }
             cout << endl;

             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   UTest[a][b] = 0.;
                }
             }

             cout << "Ud2d1Udag" << endl << endl;
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   for(int c = 0; c<3; c++)
                   {
                      //                                 U d1d2 Udag
                      UTest[a][b] = UTest[a][b] - (U[nuc][xx][yy][a][c]*di[1][c][b])/(144.*g*h*h);
                   }
                   cout << UTest[a][b] << "  ";
                }
                cout << endl;
             }
             cout << endl;
/**********************************************/
/*
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   //d2 U
                   di[0][a][b] = U[nuc][xx][yy-2][a][b]
                               - U[nuc][xx][yy-1][a][b]*8
                               + U[nuc][xx][yy+1][a][b]*8.
                               - U[nuc][xx][yy+2][a][b];
                   //d1 U
                   di[1][a][b] = U[nuc][xx-2][yy][a][b]
                               - U[nuc][xx-1][yy][a][b]*8
                               + U[nuc][xx+1][yy][a][b]*8.
                               - U[nuc][xx+2][yy][a][b];
                   //d2 Udag
                   di[2][a][b] = conj(U[nuc][xx][yy-2][b][a])
                               - conj(U[nuc][xx][yy-1][b][a])*8
                               + conj(U[nuc][xx][yy+1][b][a])*8.
                               - conj(U[nuc][xx][yy+2][b][a]);
                   //d1 Udag
                   di[3][a][b] = conj(U[nuc][xx-2][yy][b][a])
                               - conj(U[nuc][xx-1][yy][b][a])*8
                               + conj(U[nuc][xx+1][yy][b][a])*8.
                               - conj(U[nuc][xx+2][yy][b][a]);
                }
             }

             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   UTest[a][b] = 0.;
                }
             }

             cout << "(d1U)Udag" << endl << endl;
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   for(int c = 0; c<3; c++)
                   {
                      //                              d1 U         Udag
                      UTest[a][b] = UTest[a][b] - I*(di[1][a][c]*conj(U[nuc][xx][yy][b][c]))/(12*g*h);
                   }
                   cout << UTest[a][b] << "  ";
                }
                cout << endl;
             }
             cout << endl;

             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   UTest[a][b] = 0.;
                }
             }

             cout << "U(d1Udag)" << endl << endl;
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   for(int c = 0; c<3; c++)
                   {
                      //                                         U         d1Udag
                      UTest[a][b] = UTest[a][b] - I*(U[nuc][xx][yy][a][c]*di[3][c][b])/(12*g*h);
                   }
                   cout << UTest[a][b] << "  ";
                }
                cout << endl;
             }
             cout << endl;

             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   UTest[a][b] = 0.;
                }
             }

             cout << "(d2U)Udag" << endl << endl;
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   for(int c = 0; c<3; c++)
                   {
                      //                              d2 U         Udag
                      UTest[a][b] = UTest[a][b] - I*(di[0][a][c]*conj(U[nuc][xx][yy][b][c]))/(12*g*h);
                   }
                   cout << UTest[a][b] << "  ";
                }
                cout << endl;
             }
             cout << endl;

             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   UTest[a][b] = 0.;
                }
             }

             cout << "U(d2Udag)" << endl << endl;
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   for(int c = 0; c<3; c++)
                   {
                      //                                         U         d2Udag
                      UTest[a][b] = UTest[a][b] - I*(U[nuc][xx][yy][a][c]*di[2][c][b])/(12*g*h);
                   }
                   cout << UTest[a][b] << "  ";
                }
                cout << endl;
             }
             cout << endl;

/**********************************************/

             for(int x = 0; x<ARRAY; x++)
             {
                for(int y = 0; y<ARRAY; y++)
                {
                      for(int c = 0; c<3; c++)
                      {
                         AxCor[0][x][y] = AxCor[0][x][y] + A[nuc][0][x][y][0][c]*A[nuc][0][HALF][HALF][c][0]
                                                         + A[nuc][0][x][y][1][c]*A[nuc][0][HALF][HALF][c][1]
                                                         + A[nuc][0][x][y][2][c]*A[nuc][0][HALF][HALF][c][2];
                         AxCor[1][x][y] = AxCor[1][x][y] + A[nuc][0][x][y][0][c]*A[nuc][0][(HALF/2)][(HALF/2)][c][0]
                                                         + A[nuc][0][x][y][1][c]*A[nuc][0][(HALF/2)][(HALF/2)][c][1]
                                                         + A[nuc][0][x][y][2][c]*A[nuc][0][(HALF/2)][(HALF/2)][c][2];
//                         AxAv[0][x][y] = AxAv[0][x][y] + A[nuc][0][x][y][c][c];
//                         AxAv[1][x][y] = AxAv[1][x][y] + A[nuc][0][x][y][0][c]*A[nuc][0][x][y][c][0]
//                                                       + A[nuc][0][x][y][1][c]*A[nuc][0][x][y][c][1]
//                                                       + A[nuc][0][x][y][2][c]*A[nuc][0][x][y][c][2];

                         AyCor[0][x][y] = AyCor[0][x][y] + A[nuc][1][x][y][0][c]*A[nuc][1][HALF][HALF][c][0]
                                                         + A[nuc][1][x][y][1][c]*A[nuc][1][HALF][HALF][c][1]
                                                         + A[nuc][1][x][y][2][c]*A[nuc][1][HALF][HALF][c][2];
                         AyCor[1][x][y] = AyCor[1][x][y] + A[nuc][1][x][y][0][c]*A[nuc][1][(HALF/2)][(HALF/2)][c][0]
                                                         + A[nuc][1][x][y][1][c]*A[nuc][1][(HALF/2)][(HALF/2)][c][1]
                                                         + A[nuc][1][x][y][2][c]*A[nuc][1][(HALF/2)][(HALF/2)][c][2];
//                         AyAv[0][x][y] = AyAv[0][x][y] + A[nuc][1][x][y][c][c];
//                         AyAv[1][x][y] = AyAv[1][x][y] + A[nuc][1][x][y][0][c]*A[nuc][1][x][y][c][0]
//                                                       + A[nuc][1][x][y][1][c]*A[nuc][1][x][y][c][1]
//                                                       + A[nuc][1][x][y][2][c]*A[nuc][1][x][y][c][2];
                      }
                }
             }
             }//end nuc loop
          // XYSYM
/*          cout << "A test" << endl;
          int xtemp = 144; int ytemp = 215; int dirtemp = 0; int nuctemp=0; int atemp=0; int btemp=2;
          cout <<  A[nuctemp][0][xtemp][ytemp][atemp][btemp] << " vs " << A[nuctemp][1][ytemp][xtemp][atemp][btemp] << endl;
          xtemp = 54; ytemp = 514; dirtemp = 0; nuctemp=0; atemp=1; btemp=2;
          cout << A[nuctemp][0][xtemp][ytemp][atemp][btemp] << " vs " << A[nuctemp][1][ytemp][xtemp][atemp][btemp] << endl;
          xtemp = 341; ytemp = 300; dirtemp = 0; nuctemp=0; atemp=2; btemp=2;
          cout << A[nuctemp][0][xtemp][ytemp][atemp][btemp] << " vs " << A[nuctemp][1][ytemp][xtemp][atemp][btemp] << endl;
          xtemp = 98; ytemp = 277; dirtemp = 1; nuctemp=0; atemp=0; btemp=2;
          cout << A[nuctemp][0][xtemp][ytemp][atemp][btemp] << " vs " << A[nuctemp][1][ytemp][xtemp][atemp][btemp] << endl;
          xtemp = 144; ytemp = 215; dirtemp = 1; nuctemp=0; atemp=1; btemp=2;
          cout << A[nuctemp][0][xtemp][ytemp][atemp][btemp] << " vs " << A[nuctemp][1][ytemp][xtemp][atemp][btemp] << endl;
          xtemp = 300; ytemp = 444; dirtemp = 1; nuctemp=0; atemp=2; btemp=2;
          cout << A[nuctemp][0][xtemp][ytemp][atemp][btemp] << " vs " << A[nuctemp][1][ytemp][xtemp][atemp][btemp] << endl << endl;
*/


             RelU();
             
/*
             ofstream APrint[4];
             APrint[0].open("Ax0.txt");
             APrint[1].open("Ay0.txt");
             APrint[2].open("Ax1.txt");
             APrint[3].open("Ay1.txt");
             
             for(int x = 0; x<ARRAY; x++)
             {
                for(int y = 0; y<ARRAY; y++)
                {
                   complex<double> Ax0=0.;
                   complex<double> Ay0=0.;
                   complex<double> Ax1=0.;
                   complex<double> Ay1=0.;
                   
                   for(int col = 0; col<8; col++)
                   {
                      for(int c = 0; c<3; c++)
                      {
                         Ax0 += A[0][0][x][y][c][c];
                         Ay0 += A[0][1][x][y][c][c];
                         Ax1 += A[1][0][x][y][c][c];
                         Ay1 += A[1][1][x][y][c][c];
                      }
                   }
                   APrint[0] << real(Ax0) << "  ";
                   APrint[1] << real(Ay0) << "  ";
                   APrint[2] << real(Ax1) << "  ";
                   APrint[3] << real(Ay1) << "  ";
                }
                APrint[0] << endl;
                APrint[1] << endl;
                APrint[2] << endl;
                APrint[3] << endl;
             }
             
             APrint[0].close();
             APrint[1].close();
             APrint[2].close();
             APrint[3].close();
*/

/*
int xm[6]; int ym[6];
xm[0] = 144;  ym[0] = 215;
xm[1] = 54;   ym[1] = 514;
xm[2] = 341;  ym[2] = 300;
xm[3] = 98;   ym[3] = 277;
xm[4] = 144;  ym[4] = 215;
xm[5] = 300;  ym[5] = 454;

          for(int samp=0;samp<6;samp++)
          {
             cout << "AXYSYM " << samp << endl;
             cout << "x: " << xm[samp] << " y: " << ym[samp] << endl;
             cout << ATest[samp][0][0][0] << "   " << ATest[samp][0][0][1] << "   " << ATest[samp][0][0][2] << endl;
             cout << ATest[samp][0][1][0] << "   " << ATest[samp][0][1][1] << "   " << ATest[samp][0][1][2] << endl;
             cout << ATest[samp][0][2][0] << "   " << ATest[samp][0][2][1] << "   " << ATest[samp][0][2][2] << endl;
             cout << "VERSUS" << endl;
             cout << ATest[samp][1][0][0] << "   " << ATest[samp][1][0][1] << "   " << ATest[samp][1][0][2] << endl;
             cout << ATest[samp][1][1][0] << "   " << ATest[samp][1][1][1] << "   " << ATest[samp][1][1][2] << endl;
             cout << ATest[samp][1][2][0] << "   " << ATest[samp][1][2][1] << "   " << ATest[samp][1][2][2] << endl << endl;
          }
*/


             DeclFij();
             DeclA2();
             DeclAi();

             //Initialize Combined Field
/*             
             ofstream PrintAi[4];
             PrintAi[0].open("AxRe.txt");
             PrintAi[1].open("AxIm.txt");
             PrintAi[2].open("AyRe.txt");
             PrintAi[3].open("AyIm.txt");

             ofstream PrintA2[2];
             PrintA2[0].open("AzRe.txt");
             PrintA2[1].open("AzIm.txt");

             ofstream PrintFij[2];
             PrintFij[0].open("EzRe.txt");
             PrintFij[1].open("EzIm.txt");
*/

             for(int x = 0; x<ARRAY; x++)
             {
                for(int y = 0; y<ARRAY; y++)
                {
                   for(int a=0; a<3; a++)
                   {
                      for(int b=0; b<3; b++)
                      {
                         // A is of form A[nuc][dir]
                         if(x-imp < 0)
                         {
                            Ai[0][0][x][y][a][b] = (0. + C*A[1][0][x+imp][y][a][b] + 2.*pedestal );
                            Ai[0][1][x][y][a][b] = (0. + C*A[1][1][x+imp][y][a][b] + 2.*pedestal);
                         }
                         else if(x+imp > ARRAY-1)
                         {
                            Ai[0][0][x][y][a][b] = (C*A[0][0][x-imp][y][a][b] + 2*pedestal);
                            Ai[0][1][x][y][a][b] = (C*A[0][1][x-imp][y][a][b] + 2*pedestal);
                         }
                         else{
                            Ai[0][0][x][y][a][b] = (C*A[0][0][x-imp][y][a][b] + C*A[1][0][x+imp][y][a][b]  + 2.*pedestal );
                            Ai[0][1][x][y][a][b] = (C*A[0][1][x-imp][y][a][b] + C*A[1][1][x+imp][y][a][b]  + 2.*pedestal );
                         }
                         /*
                         if(a==b && a==2)
                         {
                            PrintAi[0] << real(Ai[0][0][x][y][a][b]) << "  ";
                            PrintAi[1] << imag(Ai[0][0][x][y][a][b]) << "  ";
                            PrintAi[2] << real(Ai[0][1][x][y][a][b]) << "  ";
                            PrintAi[3] << imag(Ai[0][1][x][y][a][b]) << "  ";
                         }
                         */
                         for(int c=0; c<3; c++)
                         {
                            if(x-imp < 0)
                            {
                               A2[0][x][y][a][b] = A2[0][x][y][a][b] +
                                                   pedestal*(C*A[1][0][x+imp][y][c][b] + pedestal) -
                                                   (C*A[1][0][x+imp][y][a][c] + pedestal)*pedestal +
                                                   pedestal*(C*A[1][1][x+imp][y][c][b] + pedestal) -
                                                   (C*A[1][1][x+imp][y][a][c] + pedestal)*pedestal;
                            }
                            else if(x+imp > ARRAY-1)
                            {
                               A2[0][x][y][a][b] = A2[0][x][y][a][b] +
                                                   (C*A[0][0][x-imp][y][a][c] + pedestal)*pedestal -
                                                   pedestal*(C*A[0][0][x-imp][y][c][b] + pedestal) +
                                                   (C*A[0][1][x-imp][y][a][c] + pedestal)*pedestal -
                                                   pedestal*(C*A[0][1][x-imp][y][c][b] + pedestal);
                               }
                            else
                            {
                               A2[0][x][y][a][b] = A2[0][x][y][a][b] +
                                                   (C*A[0][0][x-imp][y][a][c] + pedestal)*(C*A[1][0][x+imp][y][c][b] + pedestal) -
                                                   (C*A[1][0][x+imp][y][a][c] + pedestal)*(C*A[0][0][x-imp][y][c][b] + pedestal) +
                                                   (C*A[0][1][x-imp][y][a][c] + pedestal)*(C*A[1][1][x+imp][y][c][b] + pedestal) -
                                                   (C*A[1][1][x+imp][y][a][c] + pedestal)*(C*A[0][1][x-imp][y][c][b] + pedestal);
                            }
                         }
                         A2[0][x][y][a][b] = (-I*g/2.)*A2[0][x][y][a][b];
                         /*
                         if(a==1 && b==2)
                         {
                            PrintA2[0] << real(A2[0][x][y][a][b]) << "  ";
                            PrintA2[1] << imag(A2[0][x][y][a][b]) << "  ";
                         }
                         */
                      }
                   } 
                } // end y
                /*
                PrintAi[0] << endl;
                PrintAi[1] << endl;
                PrintAi[2] << endl;
                PrintAi[3] << endl;
                
                PrintA2[0] << endl;
                PrintA2[1] << endl;
                
                PrintFij[0] << endl;
                PrintFij[1] << endl;
                */
             } // end x



/* Here */
{
             for(int x = 2; x<ARRAY-2; x++)
             {
                for(int y = 2; y<ARRAY-2; y++)
                {
                   for(int a=0; a<3; a++)
                   {
                      for(int b=0; b<3; b++)
                      {
                         for(int c=0; c<3; c++)
                         {
                            if(x-imp < 0)
                            {
                               //E0 = [Ai1,Ai2]
                               Fij[0][0][0][x][y][a][b] = Fij[0][0][0][x][y][a][b] 
                                                       + I*g*(pedestal*(C*A[1][0][x+imp][y][c][b] + pedestal)
                                                            - (C*A[1][0][x+imp][y][a][c] + pedestal)*pedestal
                                                            + pedestal*(C*A[1][1][x+imp][y][c][b] + pedestal)
                                                            - (C*A[1][1][x+imp][y][a][c] + pedestal)*pedestal
                                                             );
                               //B0 = eij Ig[Ai1,Aj2]
                               Fij[0][1][0][x][y][a][b] = Fij[0][1][0][x][y][a][b] 
                                                        + I*g*(pedestal*(C*A[1][1][x+imp][y][c][b] + pedestal)
                                                             - (C*A[1][1][x+imp][y][a][c] + pedestal)*pedestal
                                                             - pedestal*(C*A[1][0][x+imp][y][c][b] + pedestal)
                                                             + (C*A[1][0][x+imp][y][a][c] + pedestal)*pedestal
                                                              );


                          }
                            else if(x+imp > ARRAY-1)
                            {
                               //E0 = [Ai1,Ai2]
                               Fij[0][0][0][x][y][a][b] = Fij[0][0][0][x][y][a][b] 
                                                       + I*g*((C*A[0][0][x-imp][y][a][c] + pedestal)*pedestal
                                                            - pedestal*(C*A[0][0][x-imp][y][c][b] + pedestal)
                                                            + (C*A[0][1][x-imp][y][a][c] + pedestal)*pedestal
                                                            - pedestal*(C*A[0][1][x-imp][y][c][b] + pedestal)
                                                             );
                               //B0 = eij Ig[Ai1,Aj2]
                               Fij[0][1][0][x][y][a][b] = Fij[0][1][0][x][y][a][b] 
                                                        + I*g*((C*A[0][0][x-imp][y][a][c] + pedestal)*pedestal
                                                             - pedestal*(C*A[0][0][x-imp][y][c][b] + pedestal)
                                                             - (C*A[0][1][x-imp][y][a][c] + pedestal)*pedestal
                                                             + pedestal*(C*A[0][1][x-imp][y][c][b] + pedestal)
                                                              );

                            }
                            else
                            {
                               //E0 = [Ai1,Ai2]
                               Fij[0][0][0][x][y][a][b] = Fij[0][0][0][x][y][a][b] 
                                                       + I*g*((C*A[0][0][x-imp][y][a][c] + pedestal)*(C*A[1][0][x+imp][y][c][b] + pedestal)
                                                            - (C*A[1][0][x+imp][y][a][c] + pedestal)*(C*A[0][0][x-imp][y][c][b] + pedestal)
                                                            + (C*A[0][1][x-imp][y][a][c] + pedestal)*(C*A[1][1][x+imp][y][c][b] + pedestal)
                                                            - (C*A[1][1][x+imp][y][a][c] + pedestal)*(C*A[0][1][x-imp][y][c][b] + pedestal)
                                                             );
                               //B0 = eij Ig[Ai1,Aj2]
                               Fij[0][1][0][x][y][a][b] = Fij[0][1][0][x][y][a][b] 
                                                        + I*g*((C*A[0][0][x-imp][y][a][c] + pedestal)*(C*A[1][1][x+imp][y][c][b] + pedestal)
                                                             - (C*A[1][1][x+imp][y][a][c] + pedestal)*(C*A[0][0][x-imp][y][c][b] + pedestal)
                                                             - (C*A[0][1][x-imp][y][a][c] + pedestal)*(C*A[1][0][x+imp][y][c][b] + pedestal)
                                                             + (C*A[1][0][x+imp][y][a][c] + pedestal)*(C*A[0][1][x-imp][y][c][b] + pedestal)
                                                              );


                            }
                         }// end c
                         /*
                         static complex<double> di[3][3];
                            
                         di[a][b] = Ai[0][0][x][y-2][a][b]
                         - Ai[0][0][x][y-1][a][b]*8
                         + Ai[0][0][x][y+1][a][b]*8.
                         - Ai[0][0][x][y+2][a][b]

                         - Ai[0][1][x-2][y][a][b]
                         + Ai[0][1][x-1][y][a][b]*8.
                         - Ai[0][1][x+1][y][a][b]*8.
                         + Ai[0][1][x+2][y][a][b];
                         
                         di[a][b] = -di[a][b]/(12.*h);

                         for(int c=0; c<3; c++)
                         {
                            di[a][b] = di[a][b]  -I*g*(Ai[0][1][x][y][a][c]*Ai[0][0][x][y][c][b]
                                                      - Ai[0][0][x][y][a][c]*Ai[0][1][x][y][c][b]);
                         }
                         Fij[0][1][0][x][y][a][b] = di[a][b];
                         */
                      }
                   }

/*V Other B V*/
/*
                   for(int a=0; a<3; a++)
                   {
                      for(int b=0; b<3; b++)
                      {
                              
                         if(x==256 && y==332)
                         {
                         static complex<double> di[3][3];
                         static complex<double> di2[3][3][4];
                         static complex<double> di3[3][3][4];

                            di[a][b] = Ai[0][0][x][y-2][a][b]
                            - Ai[0][0][x][y-1][a][b]*8
                            + Ai[0][0][x][y+1][a][b]*8.
                            - Ai[0][0][x][y+2][a][b]

                            - Ai[0][1][x-2][y][a][b]
                            + Ai[0][1][x-1][y][a][b]*8.
                            - Ai[0][1][x+1][y][a][b]*8.
                            + Ai[0][1][x+2][y][a][b];
                            
                            di2[a][b][1] = C*A[0][0][x-imp][y-2][a][b]
                            - C*A[0][0][x-imp][y-1][a][b]*8.
                            + C*A[0][0][x-imp][y+1][a][b]*8.
                            - C*A[0][0][x-imp][y+2][a][b];

                            di2[a][b][2] = - C*A[0][1][x-2-imp][y][a][b]
                            + C*A[0][1][x-1-imp][y][a][b]*8.
                            - C*A[0][1][x+1-imp][y][a][b]*8.
                            + C*A[0][1][x+2-imp][y][a][b];
                            
                            di3[a][b][1] = C*A[1][0][x+imp][y-2][a][b]
                            - C*A[1][0][x+imp][y-1][a][b]*8.
                            + C*A[1][0][x+imp][y+1][a][b]*8.
                            - C*A[1][0][x+imp][y+2][a][b];

                            di3[a][b][2] = - C*A[1][1][x-2+imp][y][a][b]
                            + C*A[1][1][x-1+imp][y][a][b]*8.
                            - C*A[1][1][x+1+imp][y][a][b]*8.
                            + C*A[1][1][x+2+imp][y][a][b];
                            
                            di[a][b] = -di[a][b]/(12.*h);
                            di2[a][b][0] = -(di2[a][b][1]+di2[a][b][2])/(12.*h);
                            di3[a][b][0] = -(di3[a][b][1]+di3[a][b][2])/(12.*h);
                            
                         for(int c=0; c<3; c++)
                         {
                            di[a][b] = di[a][b]  -I*g*(Ai[0][1][x][y][a][c]*Ai[0][0][x][y][c][b]
                                                      - Ai[0][0][x][y][a][c]*Ai[0][1][x][y][c][b]);
                            di2[a][b][3] = di2[a][b][3]  -I*g*(C*A[0][1][x-imp][y][a][c]*C*A[0][0][x-imp][y][c][b]
                                                      - C*A[0][0][x-imp][y][a][c]*C*A[0][1][x-imp][y][c][b]);
                            di3[a][b][3] = di3[a][b][3]  -I*g*(C*A[1][1][x+imp][y][a][c]*C*A[1][0][x+imp][y][c][b]
                                                      - C*A[1][0][x+imp][y][a][c]*C*A[1][1][x+imp][y][c][b]);
                         }
                         di2[a][b][0] = di2[a][b][0]+di2[a][b][3];
                         di3[a][b][0] = di3[a][b][0]+di3[a][b][3];

                            if(a==2&&b==2)
                            {
                            cout << "ACheck:" << endl;
                            cout << "A^x_1:" << endl;
                            for(int aa=0; aa<3; aa++)
                            {
                               for(int bb=0; bb<3; bb++)
                               {
                                  cout << C*A[0][0][x][y][aa][bb] << "  ";
                               }
                               cout << endl;
                            }
                            cout << "A^x_2:" << endl;
                            for(int aa=0; aa<3; aa++)
                            {
                               for(int bb=0; bb<3; bb++)
                               {
                                  cout << C*A[1][0][x][y][aa][bb] << "  ";
                               }
                               cout << endl;
                            }
                            cout << "A^y_1:" << endl;
                            for(int aa=0; aa<3; aa++)
                            {
                               for(int bb=0; bb<3; bb++)
                               {
                                  cout << C*A[0][1][x][y][aa][bb] << "  ";
                               }
                               cout << endl;
                            }
                            cout << "A^y_2:" << endl;
                            for(int aa=0; aa<3; aa++)
                            {
                               for(int bb=0; bb<3; bb++)
                               {
                                  cout << C*A[1][1][x][y][aa][bb] << "  ";
                               }
                               cout << endl;
                            }
                            
                            cout << "Full:" << endl;
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  cout << di[a][b] << "  ";
                               }
                               cout << endl;
                            }
                            cout << "Simple:" << endl;
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  cout << Fij[0][1][0][x][y][a][b] << "  ";
                               }
                               cout << endl;
                            }
                            cout << "Single:" << endl;
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  cout << di2[a][b][0] << "  ";
                               }
                               cout << endl;
                            }
                            cout << "Single2:" << endl;
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  cout << di3[a][b][0] << "  ";
                               }
                               cout << endl;
                            }
                            
                            cout << "References: " << endl;
                            cout << "y-der" << endl;
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  cout << di2[a][b][1] << "  ";
                               }
                               cout << endl;
                            }
                            cout << endl;
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  cout << di3[a][b][1] << "  ";
                               }
                               cout << endl;
                            }
                            cout << endl;
                            cout << "x-der" << endl;
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  cout << di2[a][b][2] << "  ";
                               }
                               cout << endl;
                            }
                            cout << endl;
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  cout << di3[a][b][2] << "  ";
                               }
                               cout << endl;
                            }
                            cout << "diAj" << endl;
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  cout << di2[a][b][1]+di2[a][b][2] << "  ";
                               }
                               cout << endl;
                            }
                            cout << endl;
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  cout << di3[a][b][1]+di3[a][b][2] << "  ";
                               }
                               cout << endl;
                            }
                            cout << endl;
                            cout << "comm" << endl;
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  cout << di2[a][b][3] << "  ";
                               }
                               cout << endl;
                            }
                            cout << endl;
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  cout << di3[a][b][3] << "  ";
                               }
                               cout << endl;
                            }
                            cout << endl;
                            }

                         }
                      }
                   }
*/
                }
             }
             /*  eps ijDiDj B0 =?= 0
             di[a][b] =
              (Fij[0][1][0][x-2][y-2][a][b] - Fij[0][1][0][x-1][y-2][a][b]*8 + Fij[0][1][0][x+1][y-2][a][b]*8. - Fij[0][1][0][x+2][y-2][a][b])
            - (Fij[0][1][0][x-2][y-1][a][b] - Fij[0][1][0][x-1][y-1][a][b]*8 + Fij[0][1][0][x+1][y-1][a][b]*8. - Fij[0][1][0][x+2][y-1][a][b])*8.
            + (Fij[0][1][0][x-2][y+1][a][b] - Fij[0][1][0][x-1][y+1][a][b]*8 + Fij[0][1][0][x+1][y+1][a][b]*8. - Fij[0][1][0][x+2][y+1][a][b])*8.
            - (Fij[0][1][0][x-2][y+2][a][b] - Fij[0][1][0][x-1][y+2][a][b]*8 + Fij[0][1][0][x+1][y+2][a][b]*8. - Fij[0][1][0][x+2][y+2][a][b]);
             
             for(int c=0; c<3; c++)
             {
             di[a][b] = di[a][b]  -I*g*(Ai[0][0][x][y][a][c]*Ai[0][0][x][y][c][b]
                                                      - Ai[0][0][x][y][a][c]*Ai[0][1][x][y][c][b]);
             }
             */
}


             RelA();

	     // Matt: Just did for loops because I'm lazy and it's not part of main sequence
             for(int x = 2; x<ARRAY-2; x++)
             {
                for(int y = 2; y<ARRAY-2; y++)
                {
                   /* This appears to be testing garbage -SR
                   for(int a=0; a<3; a++)
                   {
                      for(int b=0; b<3; b++)
                      {


                         for(int c=0; c<3; c++)
                         {













                            Fij[1][1][0][x][y][a][b] = Fij[1][1][0][x][y][a][b]
                                                     -I*g*(Ai[0][1][x][y][a][c]*Ai[0][0][x][y][c][b]
                                                         - Ai[0][0][x][y][a][c]*Ai[0][1][x][y][c][b]);
                            Fij[2][1][0][x][y][a][b] = Fij[2][1][0][x][y][a][b]
                                                     -I*g*(A[0][1][x][y][a][c]*A[0][0][x][y][c][b]
                                                         - A[0][0][x][y][a][c]*A[0][1][x][y][c][b]);
                            Fij[3][1][0][x][y][a][b] = Fij[3][1][0][x][y][a][b]
                                                     -I*g*(A[1][1][x][y][a][c]*A[1][0][x][y][c][b]
                                                         - A[1][0][x][y][a][c]*A[1][1][x][y][c][b]);
                         }
                         
                         static complex<double> di; //elemnt of derivative from l
                         static complex<double> di2; //elemnt of derivative from l
                         static complex<double> di3; //elemnt of derivative from l
                         
                         di = Ai[0][0][x][y-2][a][b]
                            - Ai[0][0][x][y-1][a][b]*8.
                            + Ai[0][0][x][y+1][a][b]*8.
                            - Ai[0][0][x][y+2][a][b]

                            - Ai[0][1][x-2][y][a][b]
                            + Ai[0][1][x-1][y][a][b]*8.
                            - Ai[0][1][x+1][y][a][b]*8.
                            + Ai[0][1][x+2][y][a][b];
                            
                         di2 = A[0][0][x][y-2][a][b]
                            - A[0][0][x][y-1][a][b]*8.
                            + A[0][0][x][y+1][a][b]*8.
                            - A[0][0][x][y+2][a][b]

                            - A[0][1][x-2][y][a][b]
                            + A[0][1][x-1][y][a][b]*8.
                            - A[0][1][x+1][y][a][b]*8.
                            + A[0][1][x+2][y][a][b];
                            
                         di3 = A[1][0][x][y-2][a][b]
                            - A[1][0][x][y-1][a][b]*8.
                            + A[1][0][x][y+1][a][b]*8.
                            - A[1][0][x][y+2][a][b]

                            - A[1][1][x-2][y][a][b]
                            + A[1][1][x-1][y][a][b]*8.
                            - A[1][1][x+1][y][a][b]*8.
                            + A[1][1][x+2][y][a][b];


                         di = -di/(12.*h);
                         di2 = -di2/(12.*h); 
                         di3 = -di3/(12.*h); 
                         Fij[1][1][0][x][y][a][b] = Fij[1][1][0][x][y][a][b] + di;
                         Fij[2][1][0][x][y][a][b] = Fij[2][1][0][x][y][a][b] + di2;
                         Fij[3][1][0][x][y][a][b] = Fij[3][1][0][x][y][a][b] + di3;
                         
                      }
                   }
                   */

                   for(int a=0; a<3; a++)
                   {
                      for(int c=0; c<3; c++)
                      {
                         T00[0][x][y] = T00[0][x][y]      // EzEz+BzBz  No transverse fields
                                      + real(Fij[0][0][0][x][y][a][c]*Fij[0][0][0][x][y][c][a]
                                      +  Fij[0][1][0][x][y][a][c]*Fij[0][1][0][x][y][c][a])/4.;
                      }
                   }
                   
//                   T00Av[x][y] = T00Av[x][y] + real(T[0][0][0][x][y]);
                   
                   
                   
                   T11[0][x][y] = T00[0][x][y];
                   T22[0][x][y] = T00[0][x][y];
                   T33[0][x][y] =-T00[0][x][y];
                }
             }
             //Matt Der. Procedure
             
             //But why? No derivative for T00... -SR
/*
	     for(int x = 0; x<2; x++)
             {
                for(int y = 0; y<ARRAY; y++)
                {
                   for(int a=0; a<4; a++)
                   {
					   T[0][a][a][x][y] = T[0][a][a][2][y];
				   }
				}
			 }
			 for(int x = ARRAY - 2; x<ARRAY; x++)
             {
                for(int y = 0; y<ARRAY; y++)
                {
                   for(int a=0; a<4; a++)
                   {
					   T[0][a][a][x][y] = T[0][a][a][ARRAY - 3][y];
				   }
				}
			 }
			 for(int y = 0; y<2; y++)
             {
                for(int x = 0; x<ARRAY; x++)
                {
                   for(int a=0; a<4; a++)
                   {
					   T[0][a][a][x][y] = T[0][a][a][x][2];
				   }
				}
			 }
			 for(int y = ARRAY - 2; y<ARRAY; y++)
             {
                for(int x = 0; x<ARRAY; x++)
                {
                   for(int a=0; a<4; a++)
                   {
					   T[0][a][a][x][y] = T[0][a][a][x][ARRAY - 3];
				   }
				}
			 }
*/

/* This is a removal of the testing garbage. These Fij are never used in anything, looks like a test remnant
             for(int x = 4; x<ARRAY-4; x++)
             {
                for(int y = 4; y<ARRAY-4; y++)
                {
                   for(int a=0; a<3; a++)
                   {
                      for(int b=0; b<3; b++)
                      {
                         Fij[1][1][0][x][y][a][b] = 0.;
                         Fij[2][1][0][x][y][a][b] = 0.;
                         Fij[3][1][0][x][y][a][b] = 0.;
                      }
                   }
                }
             }
*/


//TempDecl             RelA();


       timer=time(0) - timer;
       
       WriteTime << "Gauge to A + Initial conditions: " << timer << endl;

       timer=time(0);

//cout << "Starting loop procedure..." << endl;
             for(int BigO=1; BigO<=(Order/2);BigO++) //Showtime
             {
                int klim;
                if(!LKApprox)
                {
                   klim=BigO-1;
                }
                else
                {
                   klim=0;
                }

                for(int k=0; k<=klim; k++)
                {
                   int llim;
                   if(!LKApprox)
                   {
                      llim=BigO-k-1;
                   }
                   else
                   {
                      llim=0;
                   }
                   for(int l=0; l<=llim; l++)
                   {
                      int m = BigO-k-l-1;
                      for(int i =0; i<2; i++)
                      {
                         /**Inner Comm**/
                         DeclDum(); //Temp arrays for calclation- need to hold x-y array info
						 {
						 int x = BigO*0 + 2;
                         for(int matt_index_x = 0; matt_index_x <ARRAY; matt_index_x++) // Starts at 2 + last order lowest (with x=2 for BigO=0)
                         {                                                      // Ends ar -2 + last order highest (With x=Array-2 for BigO=0)
                            x = x%ARRAY;
						if(x >= BigO*0 +2 && x < ARRAY - (2 + BigO*0))
							{
							int y = BigO*0 + 2;
							for(int matt_index_y = 0; matt_index_y <ARRAY; matt_index_y++)
                            {
							y = y%ARRAY;
                           /**********************/
                           /**********************/
                           /******Long. Recur.****/
                           /**********************/
                           /**********************/
						   if(y>= BigO*0 + 2 && y< ARRAY - (BigO*0 + 2))
						   {
                               for(int a=0; a<3; a++)
                               {
                                  for(int b=0; b<3; b++)
                                  {
                                     // [Ai,A]
                                     for(int c=0; c<3; c++)
                                     {
                                        AT[x][y][a][b] = AT[x][y][a][b]
                                                       + Ai[l][i][x][y][a][c]*A2[m][x][y][c][b]
                                                       - A2[m][x][y][a][c]*Ai[l][i][x][y][c][b];
                                     }
                                     
                                     AT[x][y][a][b] = (-I*g)*AT[x][y][a][b];
                                     
                                     // for 0th Ai term, also get derivative because really is [Di,A]
                                     if(l==0)
                                     {
                                        complex<double> di; //elemnt of derivative from l
                                        if(i==0)
                                        {
                                           di = A2[m][x-2][y][a][b]
                                              - A2[m][x-1][y][a][b]*8.
                                              + A2[m][x+1][y][a][b]*8.
                                              - A2[m][x+2][y][a][b];
                                        }
                                        else
                                        {
                                           di = A2[m][x][y-2][a][b]
                                              - A2[m][x][y-1][a][b]*8.
                                              + A2[m][x][y+1][a][b]*8.
                                              - A2[m][x][y+2][a][b];
                                        }
                                       di = -di/(12.*h); 
                                       AT[x][y][a][b] = AT[x][y][a][b] + di;
                                     }
                                     
                                     //Di = di -igAperpn
                                     
                                  } //end b loop
                               } //end a loop
							   
						   }
						   if(y< BigO*0 +2 )
						   {
							   for(int a=0; a<3; a++)
                               {
                                  for(int b=0; b<3; b++)
                                  {
									AT[x][y][a][b]  = AT[x][BigO*0 +2][a][b];
									  
								  }
							   
							   }
							   
						   }
						   
						    if(y>= ARRAY - BigO*0 -2 )
						   {
							   for(int a=0; a<3; a++)
                               {
                                  for(int b=0; b<3; b++)
                                  { 
									AT[x][y][a][b]  = AT[x][ARRAY - BigO*0 -3][a][b];									  
								  }
							   
							   }
							   
						   }
						   
						   
						   
							y++;
                            } //end y loop
						 }
						if(x < BigO*0 +2 )
						{
							for( int y = 0; y < ARRAY; y++)
							{
								for (int a = 0; a < 3; a++)
								{
									for (int b =0; b < 3; b++)
									{						
									AT[x][y][a][b]  = AT[BigO*0 +2][y][a][b];										
									}
								}
							}
							
								
						}
							
						if(x >= ARRAY - BigO*0 -2)
						{
							for( int y = 0; y < ARRAY; y++)
							{
								for (int a = 0; a < 3; a++)
								{
									for (int b =0; b < 3; b++)
									{						
									AT[x][y][a][b]  = AT[ARRAY - BigO*0 - 3][y][a][b];										
									}
								}
							}
								
						}
														
							x++;
                         } //end x loop
					  }
                        // cout << "Inner commutator 1 done" << endl;

                   
                         /**Outter Comm**/
						 {
						 int x = 2 + BigO*0;
						for(int matt_index_x = 0; matt_index_x <ARRAY; matt_index_x++) // Starts at 4 + last order lowest (with x=2 for BigO=0)
                         {                                              // Ends at -4 + last order highest (with x=Array-2 for BigO=0)
                            x = x% ARRAY;
							if(x >= 2 + BigO*0 && x < ARRAY - (2 + BigO*0))
							{
							int y = 2 + BigO*0;
							for(int matt_index_y = 0; matt_index_y <ARRAY; matt_index_y++)
                            {
								y = y%ARRAY;
								if (y >= 2 + BigO*0 && y < ARRAY - (2 + BigO*0))
								{
                               for(int a=0; a<3; a++)
                               {
                                  for(int b=0; b<3; b++)
                                  {
                                     for(int c=0; c<3; c++)
                                     {
                                        AS[x][y][a][b] = AS[x][y][a][b]
                                                       + Ai[k][i][x][y][a][c]*AT[x][y][c][b]
                                                       - AT[x][y][a][c]*Ai[k][i][x][y][c][b];
                                     }
                                     
                                     AS[x][y][a][b] = (-I*g)*AS[x][y][a][b];
                                     
                                     if(k==0)
                                     {
                                        complex<double> di; //elemnt of derivative from l
                                        if(i==0)
                                        {
                                           di  = AT[x-2][y][a][b]
                                               - AT[x-1][y][a][b]*8.
                                               + AT[x+1][y][a][b]*8.
                                               - AT[x+2][y][a][b];
                                        }
                                        else
                                        {
                                           di = AT[x][y-2][a][b]
                                              - AT[x][y-1][a][b]*8.
                                              + AT[x][y+1][a][b]*8.
                                              - AT[x][y+2][a][b];
                                        }
                                       di = -di/(12.*h); 
                                       AS[x][y][a][b] = AS[x][y][a][b] + di;
                                     }
                                     
                                     //Di = di -igAperpn
                                     
                                  }
                               }// end outter comm comps 
                                  for(int a=0; a<3; a++)
                                  {
                                     for(int b=0; b<3; b++)
                                     {
                                        A2[BigO][x][y][a][b] = A2[BigO][x][y][a][b] + AS[x][y][a][b]/(BigO*2.*(BigO*2.+2));
                                     }
                                  }
								}
								if(y < 2+ BigO*0)
								{
									 for(int a=0; a<3; a++)
									{
										for(int b=0; b<3; b++)
										{
										 A2[BigO][x][y][a][b] =  A2[BigO][x][2 + BigO*0][a][b];
									
										}
										
									}
									
								}
								if(y >= ARRAY - (2+ BigO*0))
								{
									 for(int a=0; a<3; a++)
									{
										for(int b=0; b<3; b++)
										{
										 A2[BigO][x][y][a][b] =  A2[BigO][x][ARRAY - (3+ BigO*0)][a][b];
									
										}
										
									}
									
								}
								
								y++;
                            } // end y loop
							}
						if(x < 2+BigO*0)
						{
							for( int y = 0; y < ARRAY; y++)
							{
								for (int a = 0; a < 3; a++)
								{
									for (int b =0; b < 3; b++)
									{						
									A2[BigO][x][y][a][b]  = A2[BigO][BigO*0+2][y][a][b];										
									}
								}
							}
								
						}
							
						if(x >= ARRAY - BigO*0 -2)
						{
							for( int y = 0; y < ARRAY; y++)
							{
								for (int a = 0; a < 3; a++)
								{
									for (int b =0; b < 3; b++)
									{						
									A2[BigO][x][y][a][b]  = A2[BigO][ARRAY - BigO*0 - 3][y][a][b];										
									}
								}
							}
								
						}							
						x++;	
							
							
							
							
							
							
                         } // end x loop
					  }
                        // cout << "Outter commutator 1 done" << endl;
                      RelDum();
                      } //end i loop
                      
                   }//end l loop
                }//end k loop

                      /**********************/
                      /**********************/
                      /******Perp. Recur.****/
                      /**********************/
	{             	  /**********************/                  
		int x = 2 + BigO*0;
		for(int matt_index_x = 0; matt_index_x <ARRAY; matt_index_x++) // Starts at 4 + last order lowest (with x=2 for BigO=0)
		 {                                              // Ends at -4 + last order highest (with x=Array-2 for BigO=0)
			x = x% ARRAY;
		 if(x < 2 + BigO*0)
		 {
				for(int y = 0; y < ARRAY; y++)
				{
					for (int i = 0; i <2; i++)
					{
					 for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
						 Ai[BigO][i][x][y][a][b] =  Ai[BigO][i][2+ BigO*0][y][a][b];
					
						}
						
					}
					}
				}			 
			 
			 
			 
			 
			 
			 
		 }
		 if(x>= ARRAY - (2+ BigO*0))
		 {		 
				for(int y = 0; y < ARRAY; y++)
				{
					for (int i = 0; i <2; i++)
					{
					 for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
						 Ai[BigO][i][x][y][a][b] =  Ai[BigO][i][ARRAY - (3+ BigO*0)][y][a][b];
						
						}
						
					}
					}
				}			 
			 
			 
			 
			 
			 
			 
		 }				
			else if(x >= 2 + BigO*0 && x < ARRAY - (2 + BigO*0))
			{
				//std::cout<< "matt x = " << x << std::endl;
				if(x==50)
				{
				//std::cout<< "matt BigO = " << BigO << std::endl;	
				}
			int y = 2 + BigO*0;
			for(int matt_index_y = 0; matt_index_y <ARRAY; matt_index_y++)
			{
				y = y%ARRAY;
				if (y >= 2 + BigO*0 && y < ARRAY - (2 + BigO*0))
				{
				if(x==50)
				{
				//std::cout<< "matt BigO = " << BigO <<"   matt y =" << y << std::endl;	
				}
                   for(int i=0; i<2; i++)
                   {
                      int j = (i+1)%2;
                      double eps = (j*1.-0.5)*2.; // sign of non-zero Fji
                         
                      
                      if(!LKApprox)
                      {
                         klim=BigO-1;
                      }
                      else
                      {
                         klim=0;
                      }

                      for(int k=0; k<=klim; k++)
                      {
                         int l = BigO-k-1;
                            
                         for(int a=0; a<3; a++)
                         {
                            for(int b=0; b<3; b++)
                            {
                               for(int c=0; c<3; c++)
                               {
                                  //F21 = Bz = Fij[Order][1][0][...]
                                  Ai[BigO][i][x][y][a][b] = Ai[BigO][i][x][y][a][b]
                                                          + (Ai[k][j][x][y][a][c]*Fij[l][1][0][x][y][c][b]
                                                          - Fij[l][1][0][x][y][a][c]*Ai[k][j][x][y][c][b])*(-I*g*eps);
                               }
                               
                               if(k==0)
                               {
                                  complex<double> di; //elemnt of derivative from l
                                  if(j==0)
                                  {
                                     di = Fij[l][1][0][x-2][y][a][b]
                                        - Fij[l][1][0][x-1][y][a][b]*8.
                                        + Fij[l][1][0][x+1][y][a][b]*8.
                                        - Fij[l][1][0][x+2][y][a][b];
                                  }
                                  else
                                  {
                                     di = Fij[l][1][0][x][y-2][a][b]
                                        - Fij[l][1][0][x][y-1][a][b]*8.
                                        + Fij[l][1][0][x][y+1][a][b]*8.
                                        - Fij[l][1][0][x][y+2][a][b];
                                  }
                                  di = -di/(12.*h); 
                                     
                                  Ai[BigO][i][x][y][a][b] = Ai[BigO][i][x][y][a][b] + di*eps; //eps gives overall sign of Fji
                               }
                            } //end b loop
                         } //end a loop
                      } //end k loop for first part of Aperp
                      
                      static complex<double> AT[3][3];
                      
                      int klim;
                      if(!LKApprox)
                      {
                         klim=BigO-2;
                      }
                      else
                      {
                         klim=-2;
                      }
                      for(int k=0; k<=klim; k++) //Should NOT run for BigO=1!
                      {
                         for(int l=0; l<=BigO-k-2; l++)
                         {
                            int m = BigO-k-l-2;
                            
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  AT[a][b] = 0.;
                               }
                            }
                            
                            //static complex<double> AT[3][3];
                            //static complex<double> AS[3][3];
                            
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  for(int c=0; c<3; c++)
                                  {
                                     AT[a][b] = AT[a][b] +
                                                Ai[l][i][x][y][a][c]*A2[m][x][y][c][b] -
                                                A2[m][x][y][a][c]*Ai[l][i][x][y][c][b];
                                  }
                                     
                                  AT[a][b] = (-I*g)*AT[a][b];
                                  
                                  //Di = di -igAperpn
                                  if(l==0)
                                  {
                                     complex<double> di; //elemnt of derivative from l
                                     if(i==0)
                                     {
                                        di = A2[m][x-2][y][a][b]
                                           - A2[m][x-1][y][a][b]*8.
                                           + A2[m][x+1][y][a][b]*8.
                                           - A2[m][x+2][y][a][b];
                                     }
                                     else
                                     {
                                        di = A2[m][x][y-2][a][b]
                                           - A2[m][x][y-1][a][b]*8.
                                           + A2[m][x][y+1][a][b]*8.
                                           - A2[m][x][y+2][a][b];
                                     }
                                     di = -di/(12.*h); 
                                     AT[a][b] = AT[a][b] + di;
                                  }
                               }
                            }
                               
                            //Outter Comm
                            for(int a=0; a<3; a++)
                            {
                               for(int b=0; b<3; b++)
                               {
                                  for(int c=0; c<3; c++)
                                  {
                                     Ai[BigO][i][x][y][a][b] = Ai[BigO][i][x][y][a][b] +
                                                        (I*g)*(A2[k][x][y][a][c]*AT[c][b] -
                                                               AT[a][c]*A2[k][x][y][c][b]);
                                  }
                               }
                            }
                         } // end l loop
                      } // end k loop
                         
                      for(int a=0; a<3; a++)
                      {
                         for(int b=0; b<3; b++)
                         {
                            Ai[BigO][i][x][y][a][b] = Ai[BigO][i][x][y][a][b]/(BigO*BigO*4.);
                         }
                      }
                   } //end i loop for perp   



				   
				}
				if(y < 2+ BigO*0)
				{
					for (int i = 0; i <2; i++)
					{
					 for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
						 Ai[BigO][i][x][y][a][b] =  Ai[BigO][i][x][2+ BigO*0][a][b];
					
						}
						
					}
					}
				}
				if(y >= ARRAY - (2+ BigO*0))
				{
					for (int i = 0; i <2; i++)
					{
					 for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
						 Ai[BigO][i][x][y][a][b] =  Ai[BigO][i][x][ARRAY - (3+ BigO*0)][a][b];
					
						}
						
					}
					}
					
				}								
								
								
								
								
				y++;			
                } // end y for A
		 }				
				
				x++;
                }// end x for A
				
				
			 }
              //  cout << "Perp fields done" << endl;

                      /**********************/
                      /**********************/
                      /******Chrm. Fields****/
                      /**********************/
                      /**********************/
	{
		int x = 2 + BigO*0;
		for(int matt_index_x = 0; matt_index_x <ARRAY; matt_index_x++) // Starts at 4 + last order lowest (with x=2 for BigO=0)
		 {                                              // Ends at -4 + last order highest (with x=Array-2 for BigO=0)
			x = x% ARRAY;
			if(x >= 2 + BigO*0 && x < ARRAY - (2 + BigO*0))
			{
			int y = 2 + BigO*0;
			for(int matt_index_y = 0; matt_index_y <ARRAY; matt_index_y++)
			{
				y = y%ARRAY;
				if (y >= 2 + BigO*0 && y < ARRAY - (2 + BigO*0))
				{
                      for(int a=0; a<3; a++)
                      {
                         for(int b=0; b<3; b++)
                         {
                            //E0 (F+- = E0)
                            //Recent Change: if above is true, then divide by -2 to get E0
                            Fij[BigO][0][0][x][y][a][b] = -(2.*(BigO+1))*A2[BigO][x][y][a][b];
                            
                            //Ei ( (Fi+ + Fi-)/sqrt(2) ) = Ei
                            for(int i=0; i<2;i++)
                            {
                               //BigO*2th order of Ai gives (BigO-1)*2th order of Ei
                               //eg. 2nd in Ai (BigO -> BigO*2th order for even terms) gives 1st in Ei (BigO -> BigO*2+1th order for odd)
                               Fij[BigO][0][i+1][x][y][a][b] = -(2.*BigO*Ai[BigO][i][x][y][a][b]);
                            }
                            
                            //B0 from Fij, with ij=21 (=10 indexes)
                            for(int k=0; k<=BigO; k++)
                            {
                               int l = BigO-k;
                               for(int c=0; c<3; c++)
                               {
                               //F21 (10)
                               Fij[BigO][1][0][x][y][a][b] = Fij[BigO][1][0][x][y][a][b]
                                                           + (Ai[k][1][x][y][a][c]*Ai[l][0][x][y][c][b]
                                                           -  Ai[l][0][x][y][a][c]*Ai[k][1][x][y][c][b])*(-I*g);
                               }
//                               DerVsCom << "["<< (Ai[k][1][x][y][a][0]*Ai[l][0][x][y][0][b]
//                                          -  Ai[l][0][x][y][a][0]*Ai[k][1][x][y][0][b]
//                                          +  Ai[k][1][x][y][a][1]*Ai[l][0][x][y][1][b]
//                                          -  Ai[l][0][x][y][a][1]*Ai[k][1][x][y][1][b]
//                                          +  Ai[k][1][x][y][a][2]*Ai[l][0][x][y][2][b]
//                                          -  Ai[l][0][x][y][a][2]*Ai[k][1][x][y][2][b]
//                                            )*(-I*g) << "  , ";
                            }
                            
                            //only terms of same order for derivative
                            complex<double> di; //elemnt of derivative from l
                            //diAj
                            di = Ai[BigO][0][x][y-2][a][b]
                               - Ai[BigO][0][x][y-1][a][b]*8.
                               + Ai[BigO][0][x][y+1][a][b]*8.
                               - Ai[BigO][0][x][y+2][a][b]
                            //-djAi
                               - Ai[BigO][1][x-2][y][a][b]
                               + Ai[BigO][1][x-1][y][a][b]*8.
                               - Ai[BigO][1][x+1][y][a][b]*8.
                               + Ai[BigO][1][x+2][y][a][b];
                            di = -di/(12.*h);
                            
//                            DerVsCom << di << "] ";
                            
                            Fij[BigO][1][0][x][y][a][b] = Fij[BigO][1][0][x][y][a][b] + di;
                            
                            //Bi 
                            //(Fi+ - Fi- = -*sqrt(2)epsji3*Bj)
                            //Fi3 = tau*[D^i,A]
                            //Bj = eps(ij)*tau*[D^i,A]    NO FACTOR 2, BECAUSE WE CAN'T SUM OVER BOTH ENTRIES
                            //Runs to one smaller order because tau included
                            for(int i=0; i<2;i++)
                            {
                               double eps;// Index of D^i to calculata Bj
                               if( i == 0 ) // F13 = By
                               {
                                  eps = 1.;
                               }
                               
                               else // F23 = -Bx
                               {
                                  eps = -1.;
                               }
                               for(int k=0; k<=BigO-1; k++)
                               {
                                  
                                  int m = BigO-1-k;
                                  
                                  for(int c=0; c<3; c++) //[2-i gets correct x/y based on derivative D^i]
                                  {
                                      Fij[BigO][1][2-i][x][y][a][b] = Fij[BigO][1][2-i][x][y][a][b]+
                                                                     (Ai[k][i][x][y][a][c]*A2[m][x][y][c][b]
                                                                    - A2[m][x][y][a][c]*Ai[k][i][x][y][c][b])*(-I*g*eps);
                                      if(i==1 && BigO==2 && x==50 && y==50)
                                      {
                                      Bx[k][a][b] = Bx[k][a][b]+
                                                (Ai[k][i][x][y][a][c]*A2[m][x][y][c][b]
                                               - A2[m][x][y][a][c]*Ai[k][i][x][y][c][b])*(-I*g*eps);
                                      }
                                  }
                                  
                                  if(k==0)
                                  {
                                     complex<double> di; //elemnt of derivative from l
                                     if(i==0)
                                     {
                                        di = A2[m][x-2][y][a][b]
                                           - A2[m][x-1][y][a][b]*8.
                                           + A2[m][x+1][y][a][b]*8.
                                           - A2[m][x+2][y][a][b];
                                     }
                                     else
                                     {
                                        di = A2[m][x][y-2][a][b]
                                           - A2[m][x][y-1][a][b]*8.
                                           + A2[m][x][y+1][a][b]*8.
                                           - A2[m][x][y+2][a][b];
                                     }
                                     di = -di/(12.*h);
                                     
                                     if(i==1 && x==50 && y==50 && BigO==2)
                                     {
                                     Bx[2][a][b] = Bx[2][a][b]+
                                                   di*eps;
                                     }
                                      
                                     Fij[BigO][1][2-i][x][y][a][b] = Fij[BigO][1][2-i][x][y][a][b] + di*eps;
                                  }
                               }//end k loop
                            }// end i loop
                         }
                      }
					  
				}
				if(y < 2+ BigO*0)
				{
					for (int i = 0; i <2; i++)
					{
					 for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
						Fij[BigO][0][0][x][y][a][b] = Fij[BigO][0][0][x][2 + BigO*0][a][b];
						Fij[BigO][0][1][x][y][a][b] = Fij[BigO][0][1][x][2 + BigO*0][a][b];
						Fij[BigO][0][2][x][y][a][b] = Fij[BigO][0][2][x][2 + BigO*0][a][b];
						Fij[BigO][1][0][x][y][a][b] = Fij[BigO][1][0][x][2 + BigO*0][a][b];
						Fij[BigO][1][1][x][y][a][b] = Fij[BigO][1][1][x][2 + BigO*0][a][b];
						Fij[BigO][1][2][x][y][a][b] = Fij[BigO][1][2][x][2 + BigO*0][a][b];
					
						}
						
					}
					}
				}
				if(y >= ARRAY - (2+ BigO*0))
				{
					for (int i = 0; i <2; i++)
					{
					 for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
						Fij[BigO][0][0][x][y][a][b] = Fij[BigO][0][0][x][ARRAY - (3+ BigO*0)][a][b];
						Fij[BigO][0][1][x][y][a][b] = Fij[BigO][0][1][x][ARRAY - (3+ BigO*0)][a][b];
						Fij[BigO][0][2][x][y][a][b] = Fij[BigO][0][2][x][ARRAY - (3+ BigO*0)][a][b];
						Fij[BigO][1][0][x][y][a][b] = Fij[BigO][1][0][x][ARRAY - (3+ BigO*0)][a][b];
						Fij[BigO][1][1][x][y][a][b] = Fij[BigO][1][1][x][ARRAY - (3+ BigO*0)][a][b];
						Fij[BigO][1][2][x][y][a][b] = Fij[BigO][1][2][x][ARRAY - (3+ BigO*0)][a][b];
					
						}
						
					}
					}
					
				}
				
				y++;
				} // end y loop
			}
		 if(x < 2 + BigO*0)
		 {
				for(int y = 0; y < ARRAY; y++)
				{
					for (int i = 0; i <2; i++)
					{
					 for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
						Fij[BigO][0][0][x][y][a][b] = Fij[BigO][0][0][2 + BigO*0][y][a][b];
						Fij[BigO][0][1][x][y][a][b] = Fij[BigO][0][1][2 + BigO*0][y][a][b];
						Fij[BigO][0][2][x][y][a][b] = Fij[BigO][0][2][2 + BigO*0][y][a][b];
						Fij[BigO][1][0][x][y][a][b] = Fij[BigO][1][0][2 + BigO*0][y][a][b];
						Fij[BigO][1][1][x][y][a][b] = Fij[BigO][1][1][2 + BigO*0][y][a][b];
						Fij[BigO][1][2][x][y][a][b] = Fij[BigO][1][2][2 + BigO*0][y][a][b];
					
						}
						
					}
					}
				}			 
			 
			 
			 
			 
			 
			 
		 }
		 if(x>= ARRAY - (2+ BigO*0))
		 {
				for(int y = 0; y < ARRAY; y++)
				{
					for (int i = 0; i <2; i++)
					{
					 for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
						Fij[BigO][0][0][x][y][a][b] = Fij[BigO][0][0][ARRAY - (3+ BigO*0)][y][a][b];
						Fij[BigO][0][1][x][y][a][b] = Fij[BigO][0][1][ARRAY - (3+ BigO*0)][y][a][b];
						Fij[BigO][0][2][x][y][a][b] = Fij[BigO][0][2][ARRAY - (3+ BigO*0)][y][a][b];
						Fij[BigO][1][0][x][y][a][b] = Fij[BigO][1][0][ARRAY - (3+ BigO*0)][y][a][b];
						Fij[BigO][1][1][x][y][a][b] = Fij[BigO][1][1][ARRAY - (3+ BigO*0)][y][a][b];
						Fij[BigO][1][2][x][y][a][b] = Fij[BigO][1][2][ARRAY - (3+ BigO*0)][y][a][b];
					
						}
						
					}
					}
				}			 
			 
			 
			 
			 
			 
			 
		 }
				
				
			x++;
//            DerVsCom<< endl;
                } // end x loop

				/* CHECKS BY RJF */
/*				
				cout << "Yang-Mills Tests" << endl << endl;
				if(BigO==1)
				{
                    complex<double> Summed[3][3];
                    for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
                           Summed[a][b] = 0;
                        }
                    }
                    cout << "nu = 0" << endl << endl;
                    
                    int x=298; int y=298;
					
					cout << "ig[Ai,Ei]" << endl;
                    for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p2;
							p2=0;
							
							for(int c=0; c<3; c++) 
							{
							    p2=p2+(Ai[0][0][x][y][a][c]*Fij[1][0][1][x][y][c][b]
                                - Fij[1][0][1][x][y][a][c]*Ai[0][0][x][y][c][b]
								+Ai[0][1][x][y][a][c]*Fij[1][0][2][x][y][c][b]
                                - Fij[1][0][2][x][y][a][c]*Ai[0][1][x][y][c][b])*(I*g);
							}
							cout << p2 << ", " ;
							Summed[a][b] = Summed[a][b] + p2;
						}
						cout << endl;
					}
					cout << endl;
					
					cout << "ig[t*A,Ez]" << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p3;
							p3=0;
							
							for(int c=0; c<3; c++) 
							{
								p3=p3+(A2[0][x][y][a][c]*Fij[0][0][0][x][y][c][b]
                                - Fij[0][0][0][x][y][a][c]*A2[0][x][y][c][b])*(I*g);
							}
							cout << p3 << ", " ;
							Summed[a][b] = Summed[a][b] + p3;
						}
						cout << endl;
					}
					cout << endl;
					
					cout << "n-1 summed:" << endl;
					
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							cout << Summed[a][b] << ", " ;
						}
						cout << endl;
					}
					cout << endl;

					cout << "diEi" << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p1;
							p1=0;
							
							p1 = Fij[1][0][1][x-2][y][a][b]
                               - Fij[1][0][1][x-1][y][a][b]*8.
                               + Fij[1][0][1][x+1][y][a][b]*8.
                               - Fij[1][0][1][x+2][y][a][b]

                               + Fij[1][0][2][x][y-2][a][b]
                               - Fij[1][0][2][x][y-1][a][b]*8.
                               + Fij[1][0][2][x][y+1][a][b]*8.
                               - Fij[1][0][2][x][y+2][a][b];
							p1=p1/(12*h);
							
							cout << p1 << ", " ;
							Summed[a][b] = Summed[a][b] + p1;
						}
						cout << endl;
					}
					cout << endl;
					
					/***/
/*
					cout << "dA check " << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p1;
							p1=0;
							
							p1 = - Ai[1][0][x-2][y][a][b]
                                 + Ai[1][0][x-1][y][a][b]*8.
                                 - Ai[1][0][x+1][y][a][b]*8.
                                 + Ai[1][0][x+2][y][a][b]

                                 - Ai[1][1][x][y-2][a][b]
                                 + Ai[1][1][x][y-1][a][b]*8.
                                 - Ai[1][1][x][y+1][a][b]*8.
                                 + Ai[1][1][x][y+2][a][b];

							p1=p1/(12*h);
							
							cout << 2*p1 << ", " ;
						}
						cout << endl;
					}
					cout << endl;
                    /***/
/*
					cout << "n summed:" << endl;
					
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							cout << Summed[a][b] << ", " ;
						}
						cout << endl;
					}
					cout << endl;
					
					/*
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{			
							complex<double> p1;
							p1=0;
							
							p1 = Fij[0][1][0][x][y-2][a][b]
                                           - Fij[0][1][0][x][y-1][a][b]*8.
                                           + Fij[0][1][0][x][y+1][a][b]*8.
                                           - Fij[0][1][0][x][y+2][a][b];
							p1=0.5*p1/(12*h);

					
							cout << p1 << " ..... " << Fij[1][0][1][x][y][a][b] << ", " ;
						}
						cout << endl;
					}
					*/
/*
					cout << endl << "nu = 1" << endl << endl;

                    for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
                           Summed[a][b] = 0;
                        }
                    }
                    
					cout << "d0E1" << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p3;
							p3=0;

								p3=-Fij[1][0][1][x][y][a][b];

							cout << p3 << ", " ;
							Summed[a][b] = Summed[a][b] + p3;
						}
						cout << endl;
					}
					cout << endl;
					
					cout << "-d2B0 eps21" << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p1;
							p1=0;
							
							p1 =   Fij[0][1][0][x][y-2][a][b]
							     - Fij[0][1][0][x][y-1][a][b]*8
							     + Fij[0][1][0][x][y+1][a][b]*8
							     - Fij[0][1][0][x][y+2][a][b];

							p1=p1/(12*h);
							
							cout << p1 << ", " ;
							Summed[a][b] = Summed[a][b] + p1;
						}
						cout << endl;
					}
					cout << endl;
					
					cout << "ig[A^2,-B0] eps21" << endl;
                    for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p2;
							p2=0;
							
							for(int c=0; c<3; c++) 
							{
                                p2=p2+(
                                  Ai[0][1][x][y][a][c]*Fij[0][1][0][x][y][c][b]
                                - Fij[0][1][0][x][y][a][c]*Ai[0][1][x][y][c][b] 
                                      )*I*g;
							}
							cout << p2 << ", " ;
							Summed[a][b] = Summed[a][b] + p2;
						}
						cout << endl;
					}
					cout << endl;

					cout << "n-1 summed:" << endl;
					
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							cout << Summed[a][b] << ", " ;
						}
						cout << endl;
					}
					cout << endl;

					cout << "(1/t) d3F31=F01" << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p3;
							p3=0;
							
							p3 = -Fij[1][0][1][x][y][a][b];
							cout << p3 << ", " ;
							Summed[a][b] = Summed[a][b] + p3;
						}
						cout << endl;
					}
					cout << endl;
					
					cout << "n summed:" << endl;
					
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							cout << Summed[a][b] << ", " ;
						}
						cout << endl;
					}
					cout << endl;

/*************/
/*
					cout << endl << "nu = 2" << endl << endl;
					
                    for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
                           Summed[a][b] = 0;
                        }
                    }
                    
					cout << "d0E2" << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p3;
							p3=0;

								p3=-Fij[1][0][2][x][y][a][b];

							cout << p3 << ", " ;
							Summed[a][b] = Summed[a][b] + p3;
						}
						cout << endl;
					}
					cout << endl;

					cout << "-d1B0 eps12" << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p1;
							p1=0;
							
							p1 = - Fij[0][1][0][x-2][y][a][b]
							     + Fij[0][1][0][x-1][y][a][b]*8
							     - Fij[0][1][0][x+1][y][a][b]*8
							     + Fij[0][1][0][x+2][y][a][b];

							p1=p1/(12*h);
							
							cout << p1 << ", " ;
							Summed[a][b] = Summed[a][b] + p1;
						}
						cout << endl;
					}
					cout << endl;
					
					cout << "ig[A^1,-B0] eps12" << endl;
                    for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p2;
							p2=0;
							
							for(int c=0; c<3; c++) 
							{
                                p2=p2+(
                                - Ai[0][0][x][y][a][c]*Fij[0][1][0][x][y][c][b]
                                + Fij[0][1][0][x][y][a][c]*Ai[0][0][x][y][c][b] 
                                      )*I*g;
							}
							cout << p2 << ", " ;
							Summed[a][b] = Summed[a][b] + p2;
						}
						cout << endl;
					}
					cout << endl;

					cout << "n-1 summed:" << endl;
					
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							cout << Summed[a][b] << ", " ;
						}
						cout << endl;
					}
					cout << endl;

					cout << "(1/t) d3F32=F02" << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p3;
							p3=0;
							
							p3 = -Fij[1][0][2][x][y][a][b];
							cout << p3 << ", " ;
							Summed[a][b] = Summed[a][b] + p3;
						}
						cout << endl;
					}
					cout << endl;

					cout << "n summed:" << endl;
					
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							cout << Summed[a][b] << ", " ;
						}
						cout << endl;
					}
					cout << endl;

/*************/
/*
					cout << endl << "nu = 3" << endl << endl;
					
                    for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
                           Summed[a][b] = 0;
                        }
                    }
                    
					cout << "d0E3" << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p3;
							p3=0;

								p3=2*Fij[1][0][0][x][y][a][b];

							cout << p3 << ", " ;
							Summed[a][b] = Summed[a][b] + p3;
						}
						cout << endl;
					}
					cout << endl;

					cout << "-diBj epsij " << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p1;
							p1=0;
							
							p1 = - Fij[1][1][2][x-2][y][a][b]
							     + Fij[1][1][2][x-1][y][a][b]*8
							     - Fij[1][1][2][x+1][y][a][b]*8
							     + Fij[1][1][2][x+2][y][a][b]
							     
                                 + Fij[1][1][1][x][y-2][a][b]
							     - Fij[1][1][1][x][y-1][a][b]*8
							     + Fij[1][1][1][x][y+1][a][b]*8
							     - Fij[1][1][1][x][y+2][a][b];

							p1=p1/(12*h);
							
							cout << p1 << ", " ;
							Summed[a][b] = Summed[a][b] + p1;
						}
						cout << endl;
					}
					cout << endl;

					cout << "n-1 summed:" << endl;
					
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							cout << Summed[a][b] << ", " ;
						}
						cout << endl;
					}
					cout << endl;

					cout << "-ig[Ai,Bj] epsij" << endl;
                    for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p2;
							p2=0;
							
							for(int c=0; c<3; c++) 
							{
							    p2=p2+(
                                - Ai[0][0][x][y][a][c]*Fij[1][1][2][x][y][c][b]
                                + Fij[1][1][2][x][y][a][c]*Ai[0][0][x][y][c][b]
								+ Ai[0][1][x][y][a][c]*Fij[1][1][1][x][y][c][b]
                                - Fij[1][1][1][x][y][a][c]*Ai[0][1][x][y][c][b])*(I*g);
							}
							cout << p2 << ", " ;
							Summed[a][b] = Summed[a][b] + p2;
						}
						cout << endl;
					}
					cout << endl;

					cout << "n summed:" << endl;
					
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							cout << Summed[a][b] << ", " ;
						}
						cout << endl;
					}
					cout << endl;




					cout << "dA check " << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
							complex<double> p1;
							p1=0;
							
							p1 = - Ai[1][0][x-2][y][a][b]
                                 + Ai[1][0][x-1][y][a][b]*8.
                                 - Ai[1][0][x+1][y][a][b]*8.
                                 + Ai[1][0][x+2][y][a][b]

                                 - Ai[1][1][x][y-2][a][b]
                                 + Ai[1][1][x][y-1][a][b]*8.
                                 - Ai[1][1][x][y+1][a][b]*8.
                                 + Ai[1][1][x][y+2][a][b];

							p1=p1/(12*h);
							
							cout << 2*p1 << ", " ;
						}
						cout << endl;
					}
					cout << endl;

					cout << "[D,A] check " << endl;
					for(int a=0; a<3; a++)
					{
						for(int b=0; b<3; b++)
						{
                            complex<double> p2;
                            p2=0;
							for(int c=0; c<3; c++) 
							{
							    p2=p2+(
                                - Ai[0][0][x][y][a][c]*Ai[1][0][x][y][c][b]
                                + Ai[1][0][x][y][a][c]*Ai[0][0][x][y][c][b]
								- Ai[0][1][x][y][a][c]*Ai[1][1][x][y][c][b]
                                + Ai[1][1][x][y][a][c]*Ai[0][1][x][y][c][b])*(I*g);
							}
							cout << 2*p2 << ", " ;
						}
						cout << endl;
					}
					cout << endl;


				}
*/

				
	}
                      /**********************/
                      /**********************/
                      /*******TUV Fields*****/
                      /**********************/
                      /**********************/

	{
		int x = 2 + BigO*0;
		for(int matt_index_x = 0; matt_index_x <ARRAY; matt_index_x++) // Starts at 4 + last order lowest (with x=2 for BigO=0)
		 {                                              // Ends at -4 + last order highest (with x=Array-2 for BigO=0)
			x = x% ARRAY;
			if(x >= 2 + BigO*0 && x < ARRAY - (2 + BigO*0))
			{
			int y = 2 + BigO*0;
			for(int matt_index_y = 0; matt_index_y <ARRAY; matt_index_y++)
			{
				y = y%ARRAY;
				if (y >= 2 + BigO*0 && y < ARRAY - (2 + BigO*0))
				{
                      for(int a=0; a<3; a++)
                      {
                         //Off-diagonal odds
                         //Old
/*                       T[BigO][0][1][x][y][a][b] = -Ey*Bz + -Ez*-By
/*                       T[BigO][0][2][x][y][a][b] = -Ex*-Bz + -Ez*Bx
/*                       T[BigO][1][3][x][y][a][b] = Ex*Ez + -Bz*-Bx
/*                       T[BigO][2][3][x][y][a][b] = Ey*Ez + Bz*By
                         New
/*                       T[BigO][0][1][x][y][a][b] = Sx
/*                       T[BigO][0][2][x][y][a][b] = Sy
/*                       T[BigO][1][3][x][y][a][b] = -Sig_xz
/*                       T[BigO][2][3][x][y][a][b] = -Sig_yz
*/
                            
                         for(int k=1; k<=BigO; k++) // Use for ODD counter (k -> order 2k-1)
                         {
                             int l = BigO-k;// Use for EVEN counter (l -> order 2l)
                                
                             for(int c=0; c<3; c++)
                             {
                                T01[BigO][x][y] = T01[BigO][x][y]
                                                + real(Fij[k][0][2][x][y][a][c]*Fij[l][1][0][x][y][c][a]
                                                -  Fij[l][0][0][x][y][a][c]*Fij[k][1][2][x][y][c][a])/2.;
                                T02[BigO][x][y] = T02[BigO][x][y]
                                                + real(Fij[l][0][0][x][y][a][c]*Fij[k][1][1][x][y][c][a]
                                                -  Fij[k][0][1][x][y][a][c]*Fij[l][1][0][x][y][c][a])/2.;
                                T13[BigO][x][y] =  T13[BigO][x][y]
                                                - real(Fij[k][0][1][x][y][a][c]*Fij[l][0][0][x][y][c][a]
                                                +  Fij[l][1][0][x][y][a][c]*Fij[k][1][1][x][y][c][a])/2.;
                                T23[BigO][x][y] =  T23[BigO][x][y]
                                                - real(Fij[k][0][2][x][y][a][c]*Fij[l][0][0][x][y][c][a]
                                                +  Fij[l][1][0][x][y][a][c]*Fij[k][1][2][x][y][c][a])/2.;
                             }
                         }
                            
                         //Off Diagonol evens
                         //Old
/*                       T[BigO][0][3][x][y][a][b] = -Ex*By + -Ey*-Bx
/*                       T[BigO][1][2][x][y][a][b] = Ex*Ey + By*Bx
                         New
/*                       T[BigO][0][3][x][y][a][b] = Sz
/*                       T[BigO][1][2][x][y][a][b] = -Sig_xy
*/
                         for(int k=1; k<=BigO;k++) // Both odd
                         {
                             int l = BigO-k+1;
                                
                             for(int c=0; c<3; c++)
                             {
                                T03[BigO][x][y] = T03[BigO][x][y]
                                                + real(Fij[k][0][1][x][y][a][c]*Fij[l][1][2][x][y][c][a]
                                                -  Fij[l][0][2][x][y][a][c]*Fij[k][1][1][x][y][c][a])/2.;
                                T12[BigO][x][y] = T12[BigO][x][y]
                                                - real(Fij[k][0][1][x][y][a][c]*Fij[l][0][2][x][y][c][a]
                                                +  Fij[l][1][2][x][y][a][c]*Fij[k][1][1][x][y][c][a])/2.;
                             }
                         }
                            
                       //Diagonals (even)
                       //Old
                       //T[BigO][0][0][x][y][a][b] = -Ex*Ex + -Ey*Ey -Ez*Ez + (+1)*(-2)*(0.25)*(Ex*Ex + Ey*Ey + Ez*Ez + Bx*Bx + By*By + Bz*Bz)
                       //T[BigO][1][1][x][y][a][b] =  Ex*Ex + -By*By -Bz*Bz + (-1)*(-2)*(0.25)*(Ex*Ex + Ey*Ey + Ez*Ez + Bx*Bx + By*By + Bz*Bz)
                       //T[BigO][2][2][x][y][a][b] = -Bx*Bx +  Ey*Ey -Bz*Bz + (-1)*(-2)*(0.25)*(Ex*Ex + Ey*Ey + Ez*Ez + Bx*Bx + By*By + Bz*Bz)
                       //T[BigO][3][3][x][y][a][b] = -Bx*Bx + -By*By -Ez*Ez + (-1)*(-2)*(0.25)*(Ex*Ex + Ey*Ey + Ez*Ez + Bx*Bx + By*By + Bz*Bz)
                       //New
                       //T[BigO][0][0][x][y][a][b] = (E^2 + B^2) /2
                       //T[BigO][1][1][x][y][a][b] =  -Sig_xx
                       //T[BigO][2][2][x][y][a][b] = -Sig_yy
                       //T[BigO][3][3][x][y][a][b] = -Sig_zz
                         //Even Diagonal Pieces
                         for(int k=0; k<=BigO; k++)
                         {
                            int l = BigO-k;
                            
                            for(int c=0; c<3; c++)
                            {
                               T00[BigO][x][y] = T00[BigO][x][y]      // EzEz+BzBz
                                               + real(Fij[k][0][0][x][y][a][c]*Fij[l][0][0][x][y][c][a]
                                               +  Fij[k][1][0][x][y][a][c]*Fij[l][1][0][x][y][c][a])/4.;
                               T33[BigO][x][y] = T33[BigO][x][y]      //-EzEz-BzBz
                                               - real(Fij[k][0][0][x][y][a][c]*Fij[l][0][0][x][y][c][a]
                                               +  Fij[k][1][0][x][y][a][c]*Fij[l][1][0][x][y][c][a])/2.;
                                               
                               T11[BigO][x][y] = T11[BigO][x][y]      // EzEz+BzBz (for dij term)
                                               + real(Fij[k][0][0][x][y][a][c]*Fij[l][0][0][x][y][c][a]
                                               +  Fij[k][1][0][x][y][a][c]*Fij[l][1][0][x][y][c][a])/4.;
                               T22[BigO][x][y] = T22[BigO][x][y]      // EzEz+BzBz (for dij term)
                                               + real(Fij[k][0][0][x][y][a][c]*Fij[l][0][0][x][y][c][a]
                                               +  Fij[k][1][0][x][y][a][c]*Fij[l][1][0][x][y][c][a])/4.;
                               T33[BigO][x][y] = T33[BigO][x][y]      // EzEz+BzBz (for dij term)
                                               + real(Fij[k][0][0][x][y][a][c]*Fij[l][0][0][x][y][c][a]
                                               +  Fij[k][1][0][x][y][a][c]*Fij[l][1][0][x][y][c][a])/4.;
                            }
                         }// Even Diagonals 
                            
                         //Odd Diagonal Pieces
                         for(int k=1; k<=BigO; k++)
                         {
                            int l = BigO-k+1;
                               
                            for(int c=0; c<3; c++)
                            {
                               T00[BigO][x][y] = T00[BigO][x][y]      // ExEx + EyEy + BxBx + ByBy
                                               + real(Fij[k][0][1][x][y][a][c]*Fij[l][0][1][x][y][c][a]
                                               +  Fij[k][0][2][x][y][a][c]*Fij[l][0][2][x][y][c][a]
                                               +  Fij[k][1][1][x][y][a][c]*Fij[l][1][1][x][y][c][a]
                                               +  Fij[k][1][2][x][y][a][c]*Fij[l][1][2][x][y][c][a])/4.;
                               T11[BigO][x][y] = T11[BigO][x][y]      //-ExEx - BxBx
                                               - real(Fij[k][0][1][x][y][a][c]*Fij[l][0][1][x][y][c][a]
                                               +  Fij[k][1][1][x][y][a][c]*Fij[l][1][1][x][y][c][a])/2.;
                               T22[BigO][x][y] = T22[BigO][x][y]      //-ByBy - EyEy
                                               - real(Fij[k][1][2][x][y][a][c]*Fij[l][1][2][x][y][c][a]
                                               +  Fij[k][0][2][x][y][a][c]*Fij[l][0][2][x][y][c][a])/2.;
                               
                               T11[BigO][x][y] = T11[BigO][x][y]      // ExEx + EyEy + BxBx + ByBy (for dij term)
                                               + real(Fij[k][0][1][x][y][a][c]*Fij[l][0][1][x][y][c][a]
                                               +  Fij[k][0][2][x][y][a][c]*Fij[l][0][2][x][y][c][a]
                                               +  Fij[k][1][1][x][y][a][c]*Fij[l][1][1][x][y][c][a]
                                               +  Fij[k][1][2][x][y][a][c]*Fij[l][1][2][x][y][c][a])/4.;
                               T22[BigO][x][y] = T22[BigO][x][y]      // ExEx + EyEy + BxBx + ByBy (for dij term)
                                               + real(Fij[k][0][1][x][y][a][c]*Fij[l][0][1][x][y][c][a]
                                               +  Fij[k][0][2][x][y][a][c]*Fij[l][0][2][x][y][c][a]
                                               +  Fij[k][1][1][x][y][a][c]*Fij[l][1][1][x][y][c][a]
                                               +  Fij[k][1][2][x][y][a][c]*Fij[l][1][2][x][y][c][a])/4.;
                               T33[BigO][x][y] = T33[BigO][x][y]      // ExEx + EyEy + BxBx + ByBy (for dij term)
                                               + real(Fij[k][0][1][x][y][a][c]*Fij[l][0][1][x][y][c][a]
                                               +  Fij[k][0][2][x][y][a][c]*Fij[l][0][2][x][y][c][a]
                                               +  Fij[k][1][1][x][y][a][c]*Fij[l][1][1][x][y][c][a]
                                               +  Fij[k][1][2][x][y][a][c]*Fij[l][1][2][x][y][c][a])/4.;
                            }
                         }
                      }// end a loop
				}
				
				/*WEIRD MATT STUFF
				
				if(y < 2+ BigO*0)
				{
					for (int i = 0; i <2; i++)
					{
					 for(int a=0; a<4; a++)
					{
						for(int b=a; b<4; b++)
						{
						T[BigO][a][b][x][y] = T[BigO][a][b][x][2+ BigO*0];
					
						}
						
					}
					}
				}
				if(y >= ARRAY - (2+ BigO*0))
				{
					for (int i = 0; i <2; i++)
					{
					 for(int a=0; a<4; a++)
					{
						for(int b=a; b<4; b++)
						{
						T[BigO][a][b][x][y] = T[BigO][a][b][x][ARRAY - (3+ BigO*0)];
					
						}
						
					}
					}
					
				}
				
				
				*/
				
                   y++;   
                   
                   }// end y loop
			}
			/*WEIRD MATT STUFF
		 if(x < 2 + BigO*0)
		 {
				for(int y = 0; y < ARRAY; y++)
				{
					for (int i = 0; i <2; i++)
					{
					 for(int a=0; a<4; a++)
					{
						for(int b=0; b<4; b++)
						{
						T[BigO][a][b][x][y] = T[BigO][a][b][2+ BigO*0][y];
					
						}
						
					}
					}
				}			 
			 

			
		 }
		 if(x>= ARRAY - (2+ BigO*0)) 
		 {
				for(int y = 0; y < ARRAY; y++)
				{
					for (int i = 0; i <2; i++)
					{
					 for(int a=0; a<4; a++)
					{
						for(int b=0; b<4; b++)
						{
						T[BigO][a][b][x][y] = T[BigO][a][b][ARRAY - (3+ BigO*0)][y];
					
						}
						
					}
					}
				}			 
			 
			 
			 
			 
			 
			 
		 }				   
				   
*/
				   
				   
				   
				   
				   
				   
				   x++;
                }// end x loop
		 }
		 
       timer=time(0) - timer;
       
       WriteTime << "BigO = " << BigO << " calculation: " << timer << endl << endl;

       timer=time(0);

             } // end BigO loop

//XYSYM
/*
             ofstream PrintAt[2][4];
             PrintAt[0][0].open("AxRe2.txt");
             PrintAt[0][1].open("AxIm2.txt");
             PrintAt[0][2].open("AyRe2.txt");
             PrintAt[0][3].open("AyIm2.txt");
             PrintAt[1][0].open("AxRe4.txt");
             PrintAt[1][1].open("AxIm4.txt");
             PrintAt[1][2].open("AyRe4.txt");
             PrintAt[1][3].open("AyIm4.txt");

             ofstream PrintAl[2][2];
             PrintAl[0][0].open("AzRe2.txt");
             PrintAl[0][1].open("AzIm2.txt");
             PrintAl[1][0].open("AzRe4.txt");
             PrintAl[1][1].open("AzIm4.txt");

             ofstream PrintB[2][2];
             PrintB[0][0].open("BzRe2.txt");
             PrintB[0][1].open("BzIm2.txt");
             PrintB[1][0].open("BzRe4.txt");
             PrintB[1][1].open("BzIm4.txt");

             ofstream PrintE[2][2];
             PrintE[0][0].open("EzRe2.txt");
             PrintE[0][1].open("EzIm2.txt");
             PrintE[1][0].open("EzRe4.txt");
             PrintE[1][1].open("EzIm4.txt");
             
             ofstream PrintBt[2][4];
             PrintBt[0][0].open("BxRe2.txt");
             PrintBt[0][1].open("BxIm2.txt");
             PrintBt[0][2].open("ByRe2.txt");
             PrintBt[0][3].open("ByIm2.txt");
             PrintBt[1][0].open("BxRe4.txt");
             PrintBt[1][1].open("BxIm4.txt");
             PrintBt[1][2].open("ByRe4.txt");
             PrintBt[1][3].open("ByIm4.txt");             

             ofstream PrintEt[2][4];
             PrintEt[0][0].open("ExRe2.txt");
             PrintEt[0][1].open("ExIm2.txt");
             PrintEt[0][2].open("EyRe2.txt");
             PrintEt[0][3].open("EyIm2.txt");
             PrintEt[1][0].open("ExRe4.txt");
             PrintEt[1][1].open("ExIm4.txt");
             PrintEt[1][2].open("EyRe4.txt");
             PrintEt[1][3].open("EyIm4.txt");             

             
             for(int x=0; x<ARRAY; x++)
             {
                for(int y=0; y<ARRAY; y++)
                {
                   PrintAt[0][0] << real(Ai[1][0][x][y][2][2]) << " ";
                   PrintAt[0][1] << imag(Ai[1][0][x][y][2][2]) << " ";
                   PrintAt[0][2] << real(Ai[1][1][x][y][2][2]) << " ";
                   PrintAt[0][3] << imag(Ai[1][1][x][y][2][2]) << " ";
                   PrintAt[1][0] << real(Ai[2][0][x][y][2][2]) << " ";
                   PrintAt[1][1] << imag(Ai[2][0][x][y][2][2]) << " ";
                   PrintAt[1][2] << real(Ai[2][1][x][y][2][2]) << " ";
                   PrintAt[1][3] << imag(Ai[2][1][x][y][2][2]) << " ";

                   PrintAl[0][0] << real(A2[1][x][y][2][2]) << " ";
                   PrintAl[0][1] << imag(A2[1][x][y][2][2]) << " ";
                   PrintAl[1][0] << real(A2[2][x][y][2][2]) << " ";
                   PrintAl[1][1] << imag(A2[2][x][y][2][2]) << " ";

                   PrintB[0][0] << real(Fij[1][1][0][x][y][2][2]) << " ";
                   PrintB[0][1] << imag(Fij[1][1][0][x][y][2][2]) << " ";
                   PrintB[1][0] << real(Fij[2][1][0][x][y][2][2]) << " ";
                   PrintB[1][1] << imag(Fij[2][1][0][x][y][2][2]) << " ";

                   PrintE[0][0] << real(Fij[1][0][0][x][y][2][2]) << " ";
                   PrintE[0][1] << imag(Fij[1][0][0][x][y][2][2]) << " ";
                   PrintE[1][0] << real(Fij[2][0][0][x][y][2][2]) << " ";
                   PrintE[1][1] << imag(Fij[2][0][0][x][y][2][2]) << " ";

                   PrintBt[0][0] << real(Fij[1][1][1][x][y][2][2]) << " ";
                   PrintBt[0][1] << imag(Fij[1][1][1][x][y][2][2]) << " ";
                   PrintBt[0][2] << real(Fij[1][1][2][x][y][2][2]) << " ";
                   PrintBt[0][3] << imag(Fij[1][1][2][x][y][2][2]) << " ";
                   PrintBt[1][0] << real(Fij[2][1][1][x][y][2][2]) << " ";
                   PrintBt[1][1] << imag(Fij[2][1][1][x][y][2][2]) << " ";
                   PrintBt[1][2] << real(Fij[2][1][2][x][y][2][2]) << " ";
                   PrintBt[1][3] << imag(Fij[2][1][2][x][y][2][2]) << " ";

                   PrintEt[0][0] << real(Fij[1][0][1][x][y][2][2]) << " ";
                   PrintEt[0][1] << imag(Fij[1][0][1][x][y][2][2]) << " ";
                   PrintEt[0][2] << real(Fij[1][0][2][x][y][2][2]) << " ";
                   PrintEt[0][3] << imag(Fij[1][0][2][x][y][2][2]) << " ";
                   PrintEt[1][0] << real(Fij[2][0][1][x][y][2][2]) << " ";
                   PrintEt[1][1] << imag(Fij[2][0][1][x][y][2][2]) << " ";
                   PrintEt[1][2] << real(Fij[2][0][2][x][y][2][2]) << " ";
                   PrintEt[1][3] << imag(Fij[2][0][2][x][y][2][2]) << " ";
                }
                PrintAt[0][0] << endl;
                PrintAt[0][1] << endl;
                PrintAt[0][2] << endl;
                PrintAt[0][3] << endl;
                PrintAt[1][0] << endl;
                PrintAt[1][1] << endl;
                PrintAt[1][2] << endl;
                PrintAt[1][3] << endl;

                PrintAl[0][0] << endl;
                PrintAl[0][1] << endl;
                PrintAl[1][0] << endl;
                PrintAl[1][1] << endl;

                PrintB[0][0] << endl;
                PrintB[0][1] << endl;
                PrintB[1][0] << endl;
                PrintB[1][1] << endl;

                PrintE[0][0] << endl;
                PrintE[0][1] << endl;
                PrintE[1][0] << endl;
                PrintE[1][1] << endl;

                PrintBt[0][0] << endl;
                PrintBt[0][1] << endl;
                PrintBt[0][2] << endl;
                PrintBt[0][3] << endl;
                PrintBt[1][0] << endl;
                PrintBt[1][1] << endl;
                PrintBt[1][2] << endl;
                PrintBt[1][3] << endl;

                PrintEt[0][0] << endl;
                PrintEt[0][1] << endl;
                PrintEt[0][2] << endl;
                PrintEt[0][3] << endl;
                PrintEt[1][0] << endl;
                PrintEt[1][1] << endl;
                PrintEt[1][2] << endl;
                PrintEt[1][3] << endl;
             }

             PrintAt[0][0].close();
             PrintAt[0][1].close();
             PrintAt[0][2].close();
             PrintAt[0][3].close();
             PrintAt[1][0].close();
             PrintAt[1][1].close();
             PrintAt[1][2].close();
             PrintAt[1][3].close();

             PrintAl[0][0].close();
             PrintAl[0][1].close();
             PrintAl[1][0].close();
             PrintAl[1][1].close();

             PrintB[0][0].close();
             PrintB[0][1].close();
             PrintB[1][0].close();
             PrintB[1][1].close();

             PrintE[0][0].close();
             PrintE[0][1].close();
             PrintE[1][0].close();
             PrintE[1][1].close();

             PrintBt[0][0].close();
             PrintBt[0][1].close();
             PrintBt[0][2].close();
             PrintBt[0][3].close();
             PrintBt[1][0].close();
             PrintBt[1][1].close();
             PrintBt[1][2].close();
             PrintBt[1][3].close();

             PrintEt[0][0].close();
             PrintEt[0][1].close();
             PrintEt[0][2].close();
             PrintEt[0][3].close();
             PrintEt[1][0].close();
             PrintEt[1][1].close();
             PrintEt[1][2].close();
             PrintEt[1][3].close();
*/

            for(int BigO=0; BigO<=Order/2; BigO++)
            {
               CSCharge[BigO][0] = 0;
               CSCharge[BigO][1] = 0;
               CSCharge[BigO][2] = 0;
               
               EBsq[BigO][0] = 0;
               EBsq[BigO][1] = 0;
               EBsq[BigO][2] = 0;
               EBsq[BigO][3] = 0;
               
               double tempcharge;

               for(int x=2; x<ARRAY-2; x++)
               {
                  for(int y=2; y<ARRAY-2; y++)
                  {
                     for(int a=0; a<3; a++)
                     {
                        for(int c=0; c<3; c++)
                        {
                        
                           //Longitudinal Contribution to Chern-Simmons Charge (Ez*Bz)
                           for(int k=0; k<=BigO; k++)
                           {
                              int l = BigO-k;
                               
                              CSCharge[BigO][0] = CSCharge[BigO][0] + real(Fij[k][0][0][x][y][a][c]*Fij[l][1][0][x][y][c][a]
                                                                    )/2.;
                              EBsq[BigO][0] = EBsq[BigO][0] + real(Fij[k][0][0][x][y][a][c]*Fij[l][0][0][x][y][c][a]
                                                            )/2.;
                              EBsq[BigO][1] = EBsq[BigO][1] + real(Fij[k][1][0][x][y][a][c]*Fij[l][1][0][x][y][c][a]
                                                            )/2.;

                              tempcharge = tempcharge + real(Fij[k][0][0][x][y][a][c]*Fij[l][1][0][x][y][c][a]
                                                                    )/2.;
                              EBGrid[BigO][x][y] = EBGrid[BigO][x][y] + real(Fij[k][0][0][x][y][a][c]*Fij[l][1][0][x][y][c][a]
                                                                            )/2.;
                           }
                           //Transverse Contribution to Chern-Simmons Charge (Ei*Bi)
                           for(int k=1; k<=BigO; k++)
                           {
                              int l = BigO-k+1;

                              CSCharge[BigO][0] = CSCharge[BigO][0] + real(Fij[k][0][1][x][y][a][c]*Fij[l][1][1][x][y][c][a] +
                                                                     Fij[k][0][2][x][y][a][c]*Fij[l][1][2][x][y][c][a]
                                                                     )/2.;
                              EBsq[BigO][2] = EBsq[BigO][2] + real(Fij[k][0][1][x][y][a][c]*Fij[l][0][1][x][y][c][a] +
                                                              Fij[k][0][2][x][y][a][c]*Fij[l][0][2][x][y][c][a]
                                                            )/2.;
                              EBsq[BigO][3] = EBsq[BigO][3] + real(Fij[k][1][1][x][y][a][c]*Fij[l][1][1][x][y][c][a] +
                                                              Fij[k][1][2][x][y][a][c]*Fij[l][1][2][x][y][c][a]
                                                            )/2.;
                              tempcharge = tempcharge + real(Fij[k][0][1][x][y][a][c]*Fij[l][1][1][x][y][c][a] +
                                                                     Fij[k][0][2][x][y][a][c]*Fij[l][1][2][x][y][c][a]
                                                                     )/2.;
                              EBGrid[BigO][x][y] = EBGrid[BigO][x][y] + real(Fij[k][0][1][x][y][a][c]*Fij[l][1][1][x][y][c][a] +
                                                                     Fij[k][0][2][x][y][a][c]*Fij[l][1][2][x][y][c][a]
                                                                     )/2.;
                           }
                        }
                     }
                     CSCharge[BigO][1] = CSCharge[BigO][1] + tempcharge*tempcharge;
                     CSCharge[BigO][2] = CSCharge[BigO][2] + abs(tempcharge);
                  }
               }
            }
             RelFij();
             RelA2();
             RelAi();
             
// INITIALIZE E0, B0, EDENS
       //E0[x][y][a][b]
       //B0[x][y][a][b]
       //edens[x][y][a][b]
//TempStop
/*
             for(int x = 0; x<ARRAY; x++)
             {
                for(int y = 0; y<ARRAY; y++)
                {
                   edens[x][y] = 0.;
                   
                   for(int a = 0; a<3; a++)
                   {
                      for(int b = 0; b<3; b++)
                      {
                         E0[x][y][a][b] = 0;
                         B0[x][y][a][b] = 0;
                      }
                   }
                }
             }

             for(int x = 3; x<ARRAY-3; x++)
             {
                cout << "Field Calculation x: " << x << endl;
                for(int y = 3; y<ARRAY-3; y++)
                {
                   for(int a = 0; a<3; a++)
                   {
                      for(int b = 0; b<3; b++)
                      {
                         for(int c = 0; c<3; c++)
                         {
                                 //cout << "Writing E's and B's" << endl;
                            E0[x][y][a][b] = E0[x][y][a][b]+ I*g*(
                                       A[0][0][x][y][a][c]*A[1][0][x][y][c][b] -
                                       A[1][0][x][y][a][c]*A[0][0][x][y][c][b] +
                                       A[0][1][x][y][a][c]*A[1][1][x][y][c][b] -
                                       A[1][1][x][y][a][c]*A[0][1][x][y][c][b]
                                           );
                            B0[x][y][a][b] = B0[x][y][a][b] + I*g*(
                                       A[0][0][x][y][a][c]*A[1][1][x][y][c][b] -
                                       A[1][1][x][y][a][c]*A[0][0][x][y][c][b] +
                                       A[1][0][x][y][a][c]*A[0][1][x][y][c][b] -
                                       A[0][1][x][y][a][c]*A[1][0][x][y][c][b]
                                           );
                         }
                      }
                   }
                }
             }
             
             
             for(int x = 3; x<ARRAY-3; x++)
             {
                cout << "Field Calculation x: " << x << endl;
                for(int y = 3; y<ARRAY-3; y++)
                {
                   for(int a = 0; a<3; a++)
                   {
                      for(int c = 0; c<3; c++)
                      {
                         //cout << "Writing E's and B's" << endl;
                         edens[x][y] = edens[x][y] + 
                                       (E0[x][y][a][c]*E0[x][y][c][a] +
                                       B0[x][y][a][c]*B0[x][y][c][a])/2.;
                      }
                   }
                }
             }
             
             for(int x = 0; x<ARRAY; x++)
             {
                for(int y = 0; y<ARRAY; y++)
                {
                   EDCor[0][x][y] = EDCor[0][x][y] + edens[x][y]*edens[HALF][HALF];
                   EDAv[0][x][y] = EDAv[0][x][y] + edens[x][y];
                   EDAv[1][x][y] = EDAv[1][x][y] + edens[x][y]*edens[x][y];

                   for(int c = 0; c<3; c++)
                   {
                      E0Cor[0][x][y] = E0Cor[0][x][y] + E0[x][y][0][c]*E0[HALF][HALF][c][0]
                                                      + E0[x][y][0][c]*E0[HALF][HALF][c][0]
                                                      + E0[x][y][0][c]*E0[HALF][HALF][c][0];
                      E0Cor[1][x][y] = E0Cor[1][x][y] + E0[x][y][0][c]*E0[(HALF/2)][(HALF/2)][c][0]
                                                      + E0[x][y][0][c]*E0[(HALF/2)][(HALF/2)][c][1]
                                                      + E0[x][y][0][c]*E0[(HALF/2)][(HALF/2)][c][2];
                      E0Av[0][x][y] = E0Av[0][x][y] + E0[x][y][c][c];
                      E0Av[1][x][y] = E0Av[1][x][y] + E0[x][y][0][c]*E0[x][y][c][0]
                                                    + E0[x][y][1][c]*E0[x][y][c][1]
                                                    + E0[x][y][2][c]*E0[x][y][c][2];

                      B0Cor[0][x][y] = B0Cor[0][x][y] + B0[x][y][0][c]*B0[HALF][HALF][c][0]
                                                      + B0[x][y][0][c]*B0[HALF][HALF][c][0]
                                                      + B0[x][y][0][c]*B0[HALF][HALF][c][0];
                      B0Cor[1][x][y] = B0Cor[1][x][y] + B0[x][y][0][c]*B0[(HALF/2)][(HALF/2)][c][0]
                                                      + B0[x][y][0][c]*B0[(HALF/2)][(HALF/2)][c][1]
                                                      + B0[x][y][0][c]*B0[(HALF/2)][(HALF/2)][c][2];
                      B0Av[0][x][y] = B0Av[0][x][y] + B0[x][y][c][c];
                      B0Av[1][x][y] = B0Av[1][x][y] + B0[x][y][0][c]*B0[x][y][c][0]
                                                    + B0[x][y][1][c]*B0[x][y][c][1]
                                                    + B0[x][y][2][c]*B0[x][y][c][2];
                   }
                }
             }

             for(int x = 0; x<ARRAY; x++)
             {
                for(int y = 0; y<ARRAY; y++)
                {
                   Flow[0][0][x][y] = 0.;
                   Flow[0][1][x][y] = 0.;
                   Flow[1][0][x][y] = 0.;
                   Flow[1][1][x][y] = 0.;
                }
             }
             
             for(int x = 3; x<ARRAY-3; x++)
             {
                cout << "Field Calculation x: " << x << endl;
                for(int y = 3; y<ARRAY-3; y++)
                {
                   Flow[0][0][x][y]=-(
                                       edens[x-2][y] -
                                       edens[x-1][y]*8. +
                                       edens[x+1][y]*8. -
                                       edens[x+2][y]
                                    )/(12*h);
                   
                   Flow[0][1][x][y]=-(
                                       edens[x][y-2] -
                                       edens[x][y-1]*8. +
                                       edens[x][y+1]*8. -
                                       edens[x][y+2]
                                    )/(12*h);

                   //Grads if E and B @ this x,y
                   complex<double> GE[2][3][3];complex<double> GB[2][3][3];

                   //Ai @ this x,y
                   complex<double> Ax[3][3]; complex<double> Ay[3][3];
                   
                   for(int a=0;a<3;a++)
                   {
                      for(int b=0;b<3;b++)  //Matrix form or nonsense
                      {
                         GE[0][a][b] =-( E0[x-2][y][a][b] -
                                        E0[x-1][y][a][b]*8. +
                                        E0[x+1][y][a][b]*8. -
                                        E0[x+2][y][a][b]
                                      )/(12*h);
                          GE[1][a][b] =-( E0[x][y-2][a][b] -
                                         E0[x][y-1][a][b]*8. +
                                         E0[x][y+1][a][b]*8. -
                                         E0[x][y+2][a][b]
                                       )/(12*h);

                          GB[0][a][b] =-( B0[x-2][y][a][b] -
                                         B0[x-1][y][a][b]*8. +
                                         B0[x+1][y][a][b]*8. -
                                         B0[x+2][y][a][b]
                                       )/(12*h);
                           GB[1][a][b] =-( B0[x][y-2][a][b] -
                                          B0[x][y-1][a][b]*8. +
                                          B0[x][y+1][a][b]*8. -
                                          B0[x][y+2][a][b]
                                        )/(12*h);

                         Ax[a][b]= A[0][0][x][y][a][b]+A[1][0][x][y][a][b];

                         Ay[a][b]= A[0][1][x][y][a][b]+A[1][1][x][y][a][b];
                      }
                   }
                   //Di = di - igAi(0)
                   //Bi = epsij( [Dj,B0]E0-[Dj,E0]B0 )
                   for(int a=0;a<3;a++)
                   {
                      for(int c=0;c<3;c++)
                      {
                         Flow[1][0][x][y]= Flow[1][0][x][y] +
                                           GB[1][a][c]*E0[x][y][c][a] -
                                           
                                           GE[1][a][c]*B0[x][y][c][a];
                                           
                         Flow[1][1][x][y]= Flow[1][1][x][y] -
                                           GB[0][a][c]*E0[x][y][c][a] +
                                           GE[0][a][c]*B0[x][y][c][a];
                                           
                         for(int b=0; b<3;b++)
                         {
                            Flow[1][0][x][y]= Flow[1][0][x][y] -
                                              (I*g*Ay[a][c]*B0[x][y][c][b]*E0[x][y][b][a] +
                                              B0[x][y][a][c]*I*g*Ay[c][b]*E0[x][y][b][a] +
                                     
                                              I*g*Ay[a][c]*E0[x][y][c][b]*B0[x][y][b][a] -
                                              E0[x][y][a][c]*I*g*Ay[c][b]*B0[x][y][b][a])/2.;
                            
                            Flow[1][1][x][y]= Flow[1][1][x][y] +
                                              (I*g*Ax[a][c]*B0[x][y][c][b]*E0[x][y][b][a] -
                                              B0[x][y][a][c]*I*g*Ax[c][b]*E0[x][y][b][a] -
                                     
                                              I*g*Ax[a][c]*E0[x][y][c][b]*B0[x][y][b][a] +
                                              E0[x][y][a][c]*I*g*Ax[c][b]*B0[x][y][b][a])/2.;
                         }
                      }//c
                   }//a
                }//y
             }//x

             for(int x = 0; x<ARRAY; x++)
             {
                for(int y = 0; y<ARRAY; y++)
                { 
                      FAlphAv[0][x][y] = FAlphAv[0][x][y] + Flow[0][0][x][y];
                      FAlphAv[1][x][y] = FAlphAv[1][x][y] + Flow[0][1][x][y];

                      FDirAv[0][x][y] = FDirAv[0][x][y] + Flow[1][0][x][y];
                      FDirAv[1][x][y] = FDirAv[1][x][y] + Flow[1][1][x][y];
                }
             }

             for(int x = 3; x<ARRAY-3; x++)
             {
                for(int y = 3; y<ARRAY-3; y++)
                {
                   //Grads if E and B @ this x,y
                   complex<double> GE[2][3][3];complex<double> GB[2][3][3];
                   
                   //[D,{E/B}]
                   complex<double> DE[2][3][3];complex<double> DB[2][3][3];

                   //Ai @ this x,y
                   complex<double> Ax[3][3]; complex<double> Ay[3][3];
                   
                   for(int a=0;a<3;a++)
                   {
                      for(int b=0;b<3;b++)  //Matrix form or nonsense
                      {
                         GE[0][a][b] =-( E0[x-2][y][a][b] -
                                        E0[x-1][y][a][b]*8. +
                                        E0[x+1][y][a][b]*8. -
                                        E0[x+2][y][a][b]
                                      )/(12*h);
                          GE[1][a][b] =-( E0[x][y-2][a][b] -
                                         E0[x][y-1][a][b]*8. +
                                         E0[x][y+1][a][b]*8. -
                                         E0[x][y+2][a][b]
                                       )/(12*h);

                          GB[0][a][b] =-( B0[x-2][y][a][b] -
                                         B0[x-1][y][a][b]*8. +
                                         B0[x+1][y][a][b]*8. -
                                         B0[x+2][y][a][b]
                                       )/(12*h);
                           GB[1][a][b] =-( B0[x][y-2][a][b] -
                                          B0[x][y-1][a][b]*8. +
                                          B0[x][y+1][a][b]*8. -
                                          B0[x][y+2][a][b]
                                        )/(12*h);

                         Ax[a][b]= A[0][0][x][y][a][b]+A[1][0][x][y][a][b];

                         Ay[a][b]= A[0][1][x][y][a][b]+A[1][1][x][y][a][b];
                      }//b
                   }//a

                   //[D,E] like terms
                   //Non multiplied term
                   for(int a=0;a<3;a++)
                   {
                      for(int b=0;b<3;b++)
                      {
                         for(int k=0;k<2;k++)
                         {
                            DE[k][a][b] = GE[k][a][b];
                            DB[k][a][b] = GB[k][a][b];
                         }
                      }
                   }
                   //plus multiplying terms
                   for(int a=0;a<3;a++)
                   {
                      for(int b=0;b<3;b++)
                      {
                         for(int c=0;c<3;c++)
                         {
                            DE[0][a][b] = DE[0][a][b] -I*g*(Ax[a][c]*E0[x][y][c][b] - E0[x][y][a][c]*Ax[c][b]);
                            DB[0][a][b] = DE[0][a][b] -I*g*(Ax[a][c]*B0[x][y][c][b] - B0[x][y][a][c]*Ax[c][b]);
                            DE[1][a][b] = DE[1][a][b] -I*g*(Ay[a][c]*E0[x][y][c][b] - E0[x][y][a][c]*Ay[c][b]);
                            DB[1][a][b] = DE[1][a][b] -I*g*(Ay[a][c]*B0[x][y][c][b] - B0[x][y][a][c]*Ay[c][b]);
                         }//c
                         delt[x][y] = 0;
                      }//b
                   }//a
                   
                   for(int a=0;a<3;a++)
                   {
                      for(int c=0;c<3;c++)
                      {
                         delt[x][y] = delt[x][y]              +
                                      (DE[0][a][c]*DE[0][c][a]+
                                      DE[1][a][c]*DE[1][c][a] +
                                      DB[0][a][c]*DB[0][c][a] +
                                      DB[1][a][c]*DB[1][c][a])/2.;
                      }//c
                   }//a
                }//y
             }//x

             for(int x = 0; x<ARRAY; x++)
             {
                for(int y = 0; y<ARRAY; y++)
                { 
                      deltav[x][y] = deltav[x][y] + delt[x][y];
                }
             }

             for(int x = 3; x<ARRAY-3; x++)
             {
                for(int y = 3; y<ARRAY-3; y++)
                {
                   //Div if Alpha (A) and Beta (B) flows @ this x,y
                   complex<double> DivA;complex<double> DivB;
                   DivA = ( Flow[0][0][x-2][y] -
                              Flow[0][0][x-1][y]*8. +
                              Flow[0][0][x+1][y]*8. -
                              Flow[0][0][x+2][y] +
                              
                              Flow[0][1][x][y-2] -
                              Flow[0][1][x][y-1]*8. +
                              Flow[0][1][x][y+1]*8. -
                              Flow[0][1][x][y+2]
                            )/(12*h);
                                 
                   DivB = ( Flow[1][0][x-2][y] -
                              Flow[1][0][x-1][y]*8. +
                              Flow[1][0][x+1][y]*8. -
                              Flow[1][0][x+2][y] +
                              
                              Flow[1][1][x][y-2] -
                              Flow[1][1][x][y-1]*8. +
                              Flow[1][1][x][y+1]*8. -
                              Flow[1][1][x][y+2]
                            )/(12*h);

                   // Averaging FlowDivergens
                   DivAAv[x][y] = DivAAv[x][y] + DivA;
                   DivBAv[x][y] = DivBAv[x][y] + DivB;
                   //
                   
                      edens2[0][x][y] = -(1./4.)*(DivA + delt[x][y])
                                       -(1./8.)*( DivB*sinh(2.*0.) ) //eta = 0
                                       +(1./8.)*(delt[x][y]*cosh(2.*0.));
                      edens2[1][x][y] = -(1./4.)*(DivA + delt[x][y])
                                       -(1./8.)*( DivB*sinh(2.*1.) ) //eta = 1
                                       +(1./8.)*(delt[x][y]*cosh(2.*1.));

                      PL2[0][x][y] = (1./4.)*(DivA + delt[x][y])
                                    -(1./8.)*( DivB*sinh(2.*0.) ) //eta = 0
                                    +(1./8.)*(delt[x][y]*cosh(2.*0.));
                      PL2[1][x][y] = (1./4.)*(DivA + delt[x][y])
                                    -(1./8.)*( DivB*sinh(2.*1.) ) //eta = 1
                                    +(1./8.)*(delt[x][y]*cosh(2.*1.));
                }//y
             }//x

             for(int x = 0; x<ARRAY; x++)
             {
                for(int y = 0; y<ARRAY; y++)
                { 
                      ED2Av[0][x][y] = ED2Av[0][x][y] + edens2[0][x][y];
                      ED2Av[1][x][y] = ED2Av[1][x][y] + edens2[1][x][y];

                      PL2Av[0][x][y] = PL2Av[0][x][y] + PL2[0][x][y];
                      PL2Av[1][x][y] = PL2Av[1][x][y] + PL2[1][x][y];
                }
             }
             
*/
//TempStop
             
             

/*
             cout << "Print E's and B's..." << endl;
             ofstream E0Print;
             ofstream B0Print;
             E0Print.open("E0.txt");
             B0Print.open("B0.txt");
             for(int x = 3; x<ARRAY-3; x++)
             {
                for(int y = 3; y<ARRAY-3; y++)
                {
                   double E0Val=0.;
                   double B0Val=0.;
                   for(int a = 0; a<3; a++)
                   {
                      E0Val = E0Val + real(E0[x][y][0][a]*E0[x][y][a][0] +
                                           E0[x][y][1][a]*E0[x][y][a][1] +
                                           E0[x][y][2][a]*E0[x][y][a][2] )/2.;
                                           
                      B0Val = B0Val + real(B0[x][y][0][a]*B0[x][y][a][0] +
                                           B0[x][y][1][a]*B0[x][y][a][1] +
                                           B0[x][y][2][a]*B0[x][y][a][2] )/2.;
                   }
                   E0Print << E0Val << "  ";
                       
                   B0Print << B0Val << "  ";
                }
                
                E0Print << endl;
                
                B0Print << endl;
             }
             E0Print.close();
             B0Print.close();
             cout << "E's and B's Printed!" << endl;
*/



            double CenterT00[2][Order/2 + 1];
            
            for(int BigO=0; BigO<=Order/2; BigO++)//Usual
            {
               CenterT00[0][BigO] = 0;
               CenterT00[1][BigO] = 0;
            }
            
            double CMX = 0;
            double CMY = 0;
            double TE = 0;

            //Need CoM to do Eccentricity
            // COM IS CALCULATED IN THE GRID SPACE CONTINUUM! :D
            for(int x=2+imp+Order/2+5; x<ARRAY-2-imp-Order/2-5; x++)
            {
               for(int y=2+Order/2+5; y<ARRAY-2-Order/2-5; y++)
               {
                  CMX = CMX + (x-HALF)*(T00[0][x][y] - Recall[0][0][x][y]);
                  CMY = CMY + (y-HALF)*(T00[0][x][y] - Recall[0][0][x][y]);
                  TE  = TE  + (T00[0][x][y] - Recall[0][0][x][y]);
               }
            }
            
            CMX = CMX/TE;
            CMY = CMY/TE;
            double Ecc2[3];
               double Ecc3[3];

//            for(int BigO=0; BigO<Order/2; BigO++)//ALTERNATE
            for(int BigO=0; BigO<=Order/2; BigO++)//Usual
            {
               double L0=0;
               double L1=0;
               double L2=0;
               double L3=0;
               double L4=0;
               double L5=0;
               double L6=0;
               double L7=0;
               double L8=0;
               double L9=0;
               double L10=0;
               double L11=0;
               double L12=0;
               double L13=0;
               double L14=0;
               double L15=0;

               double L01=0;
               double L02=0;
               ////
               
               double EBEL = 0;
               double EBBL = 0;
               double EBET = 0;
               double EBBT = 0;
               
               for(int x=0; x<3;x++)
               {
                  Ecc2[x] = 0.;
                  Ecc3[x] = 0.;
               }
               for(int x=2+imp+Order/2+5; x<ARRAY-2-imp-Order/2-5; x++)
               {
                  for(int y=2+Order/2+5; y<ARRAY-2-Order/2-5; y++)
                  {

                     if( abs(x-HALF) < 11 && abs(y-HALF) < 11 )
                     {
                        L02 = L02 + (T00[BigO][x][y] - Recall[BigO][0][x][y])/441;
                     }
                     if( x==HALF && y==HALF )
                     {
                        L01 = (T00[BigO][HALF][HALF] - Recall[BigO][0][HALF][HALF]);
                     }
                     L0 = L0 + (T00[BigO][x][y] - Recall[BigO][0][x][y]);
                     L1 = L1 + (T01[BigO][x][y] - Recall[BigO][1][x][y]);
                     L2 = L2 + (T02[BigO][x][y] - Recall[BigO][2][x][y]);
                     L3 = L3 + (T03[BigO][x][y] - Recall[BigO][3][x][y]);
                     L4 = L4 + (T11[BigO][x][y] - Recall[BigO][4][x][y]);
                     L5 = L5 + (T12[BigO][x][y] - Recall[BigO][5][x][y]);
                     L6 = L6 + (T13[BigO][x][y] - Recall[BigO][6][x][y]);
                     L7 = L7 + (T22[BigO][x][y] - Recall[BigO][7][x][y]);
                     L8 = L8 + (T23[BigO][x][y] - Recall[BigO][8][x][y]);
                     L9 = L9 + (T33[BigO][x][y] - Recall[BigO][9][x][y]);


                     L10 = L10 +                ((T11[BigO][x-2][y] - Recall[BigO][4][x-2][y] )
                                               - (T11[BigO][x-1][y] - Recall[BigO][4][x-1][y] )*8.
                                               + (T11[BigO][x+1][y] - Recall[BigO][4][x+1][y] )*8.
                                               - (T11[BigO][x+2][y] - Recall[BigO][4][x+2][y] )
                                                );

                     L11 = L11 +                ((T22[BigO][x][y-2] - Recall[BigO][7][x][y-2] )
                                               - (T22[BigO][x][y-1] - Recall[BigO][7][x][y-1] )*8.
                                               + (T22[BigO][x][y+1] - Recall[BigO][7][x][y+1] )*8.
                                               - (T22[BigO][x][y+2] - Recall[BigO][7][x][y+2] )
                                                );
                     L12 = L12 - (x-HALF)*h*(T03[BigO][x][y] - Recall[BigO][3][x][y]);
                     
                     L13 = CSCharge[BigO][0];
                     
                     EBEL = EBsq[BigO][0];
                     
                     EBBL = EBsq[BigO][1];
                     
                     EBET = EBsq[BigO][2];
                     
                     EBBT = EBsq[BigO][3];
                     
                     double rd = sqrt((x-HALF-CMX)*(x-HALF-CMX)*1. + (y-HALF-CMY)*(y-HALF-CMY)*1.)*h;
                     
                     if(rd != 0)
                     {
                        double tht = atan2((y-HALF-CMY)*h,(x-HALF-CMX)*h);
/*                        
                        Ecc2 = Ecc2 + (T00[BigO][x][y] - Recall[BigO][0][x][y])*sqrt( pow(rd*rd*cos(2*tht),2.0) + pow(rd*rd*sin(2*tht),2.0) )/pow(rd,2.0);
                        Ecc3 = Ecc3 + (T00[BigO][x][y] - Recall[BigO][0][x][y])*sqrt( pow(rd*rd*rd*cos(3*tht),2.0) + pow(rd*rd*rd*sin(3*tht),2.0) )/pow(rd,3.0);
*/
                        /**Idiot Fix**/
                        Ecc2[0] = Ecc2[0] + (T00[BigO][x][y] - Recall[BigO][0][x][y])*rd*rd;
                        Ecc2[1] = Ecc2[1] + (T00[BigO][x][y] - Recall[BigO][0][x][y])*rd*rd*cos(2.*tht);
                        Ecc2[2] = Ecc2[2] + (T00[BigO][x][y] - Recall[BigO][0][x][y])*rd*rd*sin(2.*tht);
                        
                        Ecc3[0] = Ecc3[0] + (T00[BigO][x][y] - Recall[BigO][0][x][y])*rd*rd*rd;
                        Ecc3[1] = Ecc3[1] + (T00[BigO][x][y] - Recall[BigO][0][x][y])*rd*rd*rd*cos(3.*tht);
                        Ecc3[2] = Ecc3[2] + (T00[BigO][x][y] - Recall[BigO][0][x][y])*rd*rd*rd*sin(3.*tht);
                     }
                  }
               }
               if(BigO==1)
               {
                  fftw_complex *T03Trans; //Will contain rho(x), rho(k), alpha(k), alpha(x)
                  fftw_plan p2;

                  T03Trans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*size+1)*(2*size+1)); // what the fuck is the 1d bullshit...


                  for(int i=0; i<2*size+1; i++)
                  {
                     for(int j=0; j<2*size+1; j++)
                     {
                        T03Trans[i*(2*size+1) + j][0] = T03[0][HALF-size+i][HALF-size+j]; //rho(x)
                        T03Trans[i*(2*size+1) + j][1] = 0;
                     }
                  }

                  p2 = fftw_plan_dft_2d((2*size+1), (2*size+1), T03Trans, T03Trans, FFTW_FORWARD, FFTW_ESTIMATE);

                  fftw_execute(p2); //rhoalpha -> rho(k)

                  for(int i=0; i<(2*size+1); i++)
                  {
                     for(int j=0; j<(2*size+1); j++)
                     {
                        T03FT[0][i][j] = T03FT[0][i][j] + T03Trans[i*(2*size+1) + j][0]*T03Trans[i*(2*size+1) + j][0];
                        T03FT[1][i][j] = T03FT[1][i][j] + T03Trans[i*(2*size+1) + j][1]*T03Trans[i*(2*size+1) + j][1];
                     }
                  }

                  fftw_free(T03Trans);
          
                  fftw_destroy_plan(p2);
               }

               for(int x=0; x<ARRAY; x++)
               {
                  for(int y=0; y<ARRAY; y++)
                  {
                     Recall[BigO][0][x][y] = T00[BigO][x][y];
                     Recall[BigO][1][x][y] = T01[BigO][x][y];
                     Recall[BigO][2][x][y] = T02[BigO][x][y];
                     Recall[BigO][3][x][y] = T03[BigO][x][y];
                     Recall[BigO][4][x][y] = T11[BigO][x][y];
                     Recall[BigO][5][x][y] = T12[BigO][x][y];
                     Recall[BigO][6][x][y] = T13[BigO][x][y];
                     Recall[BigO][7][x][y] = T22[BigO][x][y];
                     Recall[BigO][8][x][y] = T23[BigO][x][y];
                     Recall[BigO][9][x][y] = T33[BigO][x][y];
                  }
               }

               Hist[0][BigO] = L0;
               Hist[1][BigO] = L1;
               Hist[2][BigO] = L2;
               Hist[3][BigO] = L3;
               Hist[4][BigO] = L4;
               Hist[5][BigO] = L5;
               Hist[6][BigO] = L6;
               Hist[7][BigO] = L7;
               Hist[8][BigO] = L8;
               Hist[9][BigO] = L9;
               Hist[10][BigO] = L10/(12.*h);
               Hist[11][BigO] = L11/(12.*h);
               Hist[12][BigO] = L12;
               Hist[13][BigO] = L13;

               CenterT00[0][BigO] = L01;
               CenterT00[1][BigO] = L02;
               
               EBHist[0][BigO] = EBEL;
               EBHist[1][BigO] = EBBL;
               EBHist[2][BigO] = EBET;
               EBHist[3][BigO] = EBBT;
                              
               Ecc[0][BigO] = sqrt(Ecc2[1]*Ecc2[1] + Ecc2[2]*Ecc2[2])/Ecc2[0];
               Ecc[1][BigO] = sqrt(Ecc3[1]*Ecc3[1] + Ecc3[2]*Ecc3[2])/Ecc3[0];
               
/*               Hist[14][BigO] = L14;
               Hist[15][BigO] = L15;*/
            }




/*
CenterT00[1][BigO] = L02
*/

for(int BigO=0; BigO<=Order/2; BigO++)
{
Histogram[0] << (Hist[0][BigO]*h*h) << " ";
Histogram[1] << (Hist[1][BigO]*h*h) << " ";
Histogram[2] << (Hist[2][BigO]*h*h) << " ";
Histogram[3] << (Hist[3][BigO]*h*h) << " ";
Histogram[4] << (Hist[4][BigO]*h*h) << " ";
Histogram[5] << (Hist[5][BigO]*h*h) << " ";
Histogram[6] << (Hist[6][BigO]*h*h) << " ";
Histogram[7] << (Hist[7][BigO]*h*h) << " ";
Histogram[8] << (Hist[8][BigO]*h*h) << " ";
Histogram[9] << (Hist[9][BigO]*h*h) << " ";
Histogram[10] << (Hist[10][BigO]*h*h) << " ";
Histogram[11] << (Hist[11][BigO]*h*h) << " ";
Histogram[12] << (Hist[12][BigO]*h*h) << " ";
/*
//Histogram[13] << (Hist[13][BigO]*h*h) << " ";
Histogram[13] << (CSCharge[BigO][0]*h*h) << " ";
Histogram[14] << (CSCharge[BigO][1]*h*h) << " ";
Histogram[15] << (CSCharge[BigO][2]*h*h) << " ";
*/
Histogram[13] << (Hist[13][BigO]*h*h) << " ";
Histogram[14] << (Hist[14][BigO]*h*h) << " ";
Histogram[15] << (Hist[15][BigO]*h*h) << " ";

EBHistogram[0] << (EBHist[0][BigO]*h*h) << "  ";
EBHistogram[1] << (EBHist[1][BigO]*h*h) << "  ";
EBHistogram[2] << (EBHist[2][BigO]*h*h) << "  ";
EBHistogram[3] << (EBHist[3][BigO]*h*h) << "  ";


EccHist[0] << (Ecc[0][BigO]) << "  ";
EccHist[1] << (Ecc[1][BigO]) << "  ";

PeakHist[0] << CenterT00[0][BigO] << "  ";
PeakHist[1] << CenterT00[1][BigO] << "  ";
}
Histogram[0] << endl;
Histogram[1] << endl;
Histogram[2] << endl;
Histogram[3] << endl;
Histogram[4] << endl;
Histogram[5] << endl;
Histogram[6] << endl;
Histogram[7] << endl;
Histogram[8] << endl;
Histogram[9] << endl;
Histogram[10] << endl;
Histogram[11] << endl;
Histogram[12] << endl;
Histogram[13] << endl;
Histogram[14] << endl;
Histogram[15] << endl;

EBHistogram[0] << endl;
EBHistogram[1] << endl;
EBHistogram[2] << endl;
EBHistogram[3] << endl;

EccHist[0] << " CM: " << CMX << "  " << CMY;
EccHist[0] << endl;
EccHist[1] << endl;

PeakHist[0] << endl;
PeakHist[1] << endl;


if(NucMethod == 1)
{
   for(int i=0; i<3; i++)
   {
      delete [] Nucl1[i];
      delete [] Nucl2[i];
   }
   delete [] Nucl1;
   delete [] Nucl2;
}

            } //CLOSES N ITERATIONS
            
Histogram[0].close();
Histogram[1].close();
Histogram[2].close();
Histogram[3].close();
Histogram[4].close();
Histogram[5].close();
Histogram[6].close();
Histogram[7].close();
Histogram[8].close();
Histogram[9].close();
Histogram[10].close();
Histogram[11].close();
Histogram[12].close();
Histogram[13].close();
Histogram[14].close();
Histogram[15].close();

EBHistogram[0].close();
EBHistogram[1].close();
EBHistogram[2].close();
EBHistogram[3].close();

EccHist[0].close();
EccHist[1].close();

/*
       ofstream PrintRh[4];
       PrintRh[0].open("RhoCentCor.txt");
       PrintRh[1].open("RhoOffCor.txt");
       PrintRh[2].open("RhoAv.txt");
       PrintRh[3].open("RhoSqAv.txt");
       
       for(int x = 0; x<CGRho; x++)
       {
           for(int y = 0; y<CGRho; y++)
           {
//              PrintRh[0] << RhoCGCor[0][x][y] << "  ";
//              PrintRh[1] << RhoCGCor[1][x][y] << "  ";
              PrintRh[2] << RhoCGAv[0][x][y] << "  ";
              PrintRh[3] << RhoCGAv[1][x][y] << "  ";
           }
           
           PrintRh[0] << endl;
           PrintRh[1] << endl;
           PrintRh[2] << endl;
           PrintRh[3] << endl;
       }
       PrintRh[0].close();
       PrintRh[1].close();
       PrintRh[2].close();
       PrintRh[3].close();


       ofstream PrintE[4];
       PrintE[0].open("E0CentCor.txt");
       PrintE[1].open("E0OffCor.txt");
       PrintE[2].open("E0Av.txt");
       PrintE[3].open("E0SqAv.txt");
       
       ofstream PrintB[4];
       PrintB[0].open("B0CentCor.txt");
       PrintB[1].open("B0OffCor.txt");
       PrintB[2].open("B0Av.txt");
       PrintB[3].open("B0SqAv.txt");
       
       ofstream PrintED[4];
       PrintED[0].open("EDCentCor.txt");
       PrintED[1].open("EDOffCor.txt");
       PrintED[2].open("EDAv.txt");
       PrintED[3].open("EDSqAv.txt");
       
       ofstream PrintED2[2];
       PrintED2[0].open("ED2Av0.txt");
       PrintED2[1].open("ED2Av1.txt");

       ofstream PrintPL2[2];
       PrintPL2[0].open("PL2Av0.txt");
       PrintPL2[1].open("PL2Av1.txt");

       ofstream PrintAx[4];
       PrintAx[0].open("AxCentCor.txt");
       PrintAx[1].open("AxOffCor.txt");
       PrintAx[2].open("AxAv.txt");
       PrintAx[3].open("AxSqAv.txt");
       
       ofstream PrintAy[4];
       PrintAy[0].open("AyCentCor.txt");
       PrintAy[1].open("AyOffCor.txt");
       PrintAy[2].open("AyAv.txt");
       PrintAy[3].open("AySqAv.txt");
       
       ofstream PrintAlp[4];
       PrintAlp[0].open("AlpCentCor.txt");
       PrintAlp[1].open("AlpOffCor.txt");
       PrintAlp[2].open("AlpAv.txt");
       PrintAlp[3].open("AlpSqAv.txt");

       ofstream PrintFlow[4];
       PrintFlow[0].open("RadFlowX.txt");
       PrintFlow[1].open("RadFlowY.txt");
       PrintFlow[2].open("DirFlowX.txt");
       PrintFlow[3].open("DirFlowY.txt");
       
       ofstream PrintDiv[2];
       PrintDiv[0].open("DivA.txt");
       PrintDiv[1].open("DivB.txt");
       
       ofstream PrintDelt;
       PrintDelt.open("deltav.txt");

       
       
       for(int x = 0; x<ARRAY; x++)
       {
          for(int y = 0; y<ARRAY; y++)
          {
//             PrintE[0] << real(E0Cor[0][x][y]) << "  ";
//             PrintE[1] << real(E0Cor[1][x][y]) << "  ";
             PrintE[2] << real(E0Av[0][x][y]) << "  ";
             PrintE[3] << real(E0Av[1][x][y]) << "  ";
             
//             PrintB[0] << real(B0Cor[0][x][y]) << "  ";
//             PrintB[1] << real(B0Cor[1][x][y]) << "  ";
             PrintB[2] << real(B0Av[0][x][y]) << "  ";
             PrintB[3] << real(B0Av[1][x][y]) << "  ";
             
//             PrintED[0] << real(EDCor[0][x][y]) << "  ";
//             PrintED[1] << real(EDCor[1][x][y]) << "  ";
             PrintED[2] << real(EDAv[0][x][y]) << "  ";
             PrintED[3] << real(EDAv[1][x][y]) << "  ";
             
             PrintED2[0] << real(ED2Av[0][x][y]) << "  ";
             PrintED2[1] << real(ED2Av[1][x][y]) << "  ";

             PrintPL2[0] << real(PL2Av[0][x][y]) << "  ";
             PrintPL2[1] << real(PL2Av[1][x][y]) << "  ";

//             PrintAx[0] << real(AxCor[0][x][y]) << "  ";
//             PrintAx[1] << real(AxCor[1][x][y]) << "  ";
             PrintAx[2] << real(AxAv[0][x][y]) << "  ";
             PrintAx[3] << real(AxAv[1][x][y]) << "  ";
//             PrintAy[0] << real(AyCor[0][x][y]) << "  ";
//             PrintAy[1] << real(AyCor[1][x][y]) << "  ";
             PrintAy[2] << real(AyAv[0][x][y]) << "  ";
             PrintAy[3] << real(AyAv[1][x][y]) << "  ";

//             PrintAlp[0] << AlphaCor[0][x][y] << "  ";
//             PrintAlp[1] << AlphaCor[1][x][y] << "  ";
             PrintAlp[2] << AlphaAv[0][x][y] << "  ";
             PrintAlp[3] << AlphaAv[1][x][y] << "  ";

             PrintFlow[0] << real(FAlphAv[0][x][y]) << "  ";
             PrintFlow[1] << real(FAlphAv[1][x][y]) << "  ";
             PrintFlow[2] << real(FDirAv[0][x][y]) << "  ";
             PrintFlow[3] << real(FDirAv[1][x][y]) << "  ";
             
             PrintDiv[0] << real(DivAAv[x][y]) << "  ";
             PrintDiv[1] << real(DivBAv[x][y]) << "  ";
             
             PrintDelt << real(deltav[x][y]) << "  ";
          }
          PrintE[0] << endl;
          PrintE[1] << endl;
          PrintE[2] << endl;
          PrintE[3] << endl;
             
          PrintB[0] << endl;
          PrintB[1] << endl;
          PrintB[2] << endl;
          PrintB[3] << endl;

          PrintED[0] << endl;
          PrintED[1] << endl;
          PrintED[2] << endl;
          PrintED[3] << endl;
          
          PrintED2[0] << endl;
          PrintED2[0] << endl;

          PrintPL2[0] << endl;
          PrintPL2[0] << endl;

          PrintAx[0] << endl;
          PrintAx[1] << endl;
          PrintAx[2] << endl;
          PrintAx[3] << endl;
          PrintAy[0] << endl;
          PrintAy[1] << endl;
          PrintAy[2] << endl;
          PrintAy[3] << endl;

          PrintAlp[0] << endl;
          PrintAlp[1] << endl;
          PrintAlp[2] << endl;
          PrintAlp[3] << endl;

          PrintFlow[0] << endl;
          PrintFlow[1] << endl;
          PrintFlow[2] << endl;
          PrintFlow[3] << endl;

          PrintDiv[0] << endl;
          PrintDiv[1] << endl;
             
          PrintDelt << endl;
       }
          PrintE[0].close();
          PrintE[1].close();
          PrintE[2].close();
          PrintE[3].close();
             
          PrintB[0].close();
          PrintB[1].close();
          PrintB[2].close();
          PrintB[3].close();

          PrintED[0].close();
          PrintED[1].close();
          PrintED[2].close();
          PrintED[3].close();

          PrintED2[0].close();
          PrintED2[1].close();

          PrintPL2[0].close();
          PrintPL2[1].close();

          PrintAx[0].close();
          PrintAx[1].close();
          PrintAx[2].close();
          PrintAx[3].close();
          PrintAy[0].close();
          PrintAy[1].close();
          PrintAy[2].close();
          PrintAy[3].close();

          PrintAlp[0].close();
          PrintAlp[1].close();
          PrintAlp[2].close();
          PrintAlp[3].close();

          PrintFlow[0].close();
          PrintFlow[1].close();
          PrintFlow[2].close();
          PrintFlow[3].close();

          PrintDiv[0].close();
          PrintDiv[1].close();
             
          PrintDelt.close();
*/
/*
       ofstream PrintNumb;
       PrintNumb.open("Number.txt");
       PrintNumb << Numb;
       PrintNumb.close();
       
       ofstream PrintU[4];
       PrintU[0].open("UAv0Re.txt");
       PrintU[1].open("UAv1Re.txt");
       PrintU[2].open("UAv0Im.txt");
       PrintU[3].open("UAv1Im.txt");
       
       for(int x = 0; x<ARRAY; x++)
       {
          for(int y = 0; y<ARRAY; y++)
          {
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   PrintU[0] << real(AvU[0][x][y][a][b]) << "  ";
                   PrintU[1] << real(AvU[1][x][y][a][b]) << "  ";
                   PrintU[2] << imag(AvU[0][x][y][a][b]) << "  ";
                   PrintU[3] << imag(AvU[1][x][y][a][b]) << "  ";
                }
             }
          }

          PrintU[0] << endl;
          PrintU[1] << endl;
          PrintU[2] << endl;
          PrintU[3] << endl;
       }
*/
//       RelAi();

cout << "Out of N" << endl;

      ofstream PrintMu[2][2];
      PrintMu[0][0].open("Mu0.txt");
      PrintMu[0][1].open("Mu1.txt");
      PrintMu[1][0].open("MuMu.txt");
      PrintMu[1][1].open("MuDMu.txt");
      for(int x=0; x<ARRAY; x++)
      {
         for(int y=0; y<ARRAY; y++)
         {
            PrintMu[0][0] << MuShapeCenter[0][0][x][y] << "  ";
            PrintMu[0][1] << MuShapeCenter[1][0][x][y] << "  ";
            PrintMu[1][0] << MuShapeCenter[0][1][x][y] << "  ";
            PrintMu[1][1] << MuShapeCenter[1][1][x][y] << "  ";
         }
         PrintMu[0][0] << endl;
         PrintMu[0][1] << endl;
         PrintMu[1][0] << endl;
         PrintMu[1][1] << endl;
      }

      ofstream PrintAx[2];
      PrintAx[0].open("AxCentCor.txt");
      PrintAx[1].open("AxOffCor.txt");

      ofstream PrintAy[2];
      PrintAy[0].open("AyCentCor.txt");
      PrintAy[1].open("AyOffCor.txt");
      
      ofstream PrintAlphaCor[2];
      PrintAlphaCor[0].open("AlphaCorCent.txt");
      PrintAlphaCor[1].open("AlphaCorOff.txt");
      
//      ofstream PrintFij[2];
//      PrintFij[0].open("Fij0.txt");
//      PrintFij[1].open("Fij1.txt");

//      ofstream PrintT00Av;
//      PrintT00Av.open("T00Av.txt");

       for(int x = 0; x<ARRAY; x++)
       {
          for(int y = 0; y<ARRAY; y++)
          {
             PrintAx[0] << real(AxCor[0][x][y]) << "  ";
             PrintAx[1] << real(AxCor[1][x][y]) << "  ";

             PrintAy[0] << real(AyCor[0][x][y]) << "  ";
             PrintAy[1] << real(AyCor[1][x][y]) << "  ";
             
             PrintAlphaCor[0] << AlphaCor[0][x][y] << "  ";
             PrintAlphaCor[1] << AlphaCor[1][x][y] << "  ";
             
//             PrintFij[0] << sqrt( real(Fij[0][0][0][x][y][1][2])*real(Fij[0][0][0][x][y][1][2])+imag(Fij[0][0][0][x][y][1][2])*imag(Fij[0][0][0][x][y][1][2]) ) << "  ";
//             PrintFij[1] << sqrt( real(Fij[0][1][0][x][y][0][0])*real(Fij[0][1][0][x][y][0][0])+imag(Fij[0][1][0][x][y][0][0])*imag(Fij[0][1][0][x][y][0][0]) ) << "  ";

//             PrintT00Av << T00Av[x][y] << "  ";
          }
          PrintAx[0] << endl;
          PrintAx[1] << endl;
          PrintAy[0] << endl;
          PrintAy[1] << endl;
          PrintAlphaCor[0] << endl;
          PrintAlphaCor[1] << endl;
//          PrintFij[0] << endl;
//          PrintFij[1] << endl;
//          PrintT00Av << endl;
       }


       ofstream PrintT[10];
       ofstream PrintT03FT[2];
       ofstream PrintCS;
       for(int BigO=0; BigO<=(Order/2); BigO++)
       {
          ostringstream TheDamnName[10];
          ostringstream TheFTName[2];
          ostringstream CSName;
                   
          if(BigO < 10)
          {
             TheDamnName[0] << "T" << 0 << 0 << "_" << "0" << BigO << ".txt";
             TheDamnName[1] << "T" << 0 << 3 << "_" << "0" << BigO << ".txt";
             TheDamnName[2] << "T" << 1 << 3 << "_" << "0" << BigO << ".txt";

             TheDamnName[3] << "T" << 0 << 1 << "_" << "0" << BigO << ".txt";
             TheDamnName[4] << "T" << 0 << 2 << "_" << "0" << BigO << ".txt";
             TheDamnName[5] << "T" << 1 << 1 << "_" << "0" << BigO << ".txt";
             TheDamnName[6] << "T" << 1 << 2 << "_" << "0" << BigO << ".txt";
             TheDamnName[7] << "T" << 2 << 2 << "_" << "0" << BigO << ".txt";
             TheDamnName[8] << "T" << 2 << 3 << "_" << "0" << BigO << ".txt";
             TheDamnName[9] << "T" << 3 << 3 << "_" << "0" << BigO << ".txt";
       
             CSName << "EB_" << "0" << BigO << ".txt";

          }
          else
          {
             TheDamnName[0] << "T" << 0 << 0 << "_" << BigO << ".txt";
             TheDamnName[1] << "T" << 0 << 3 << "_" << BigO << ".txt";
             TheDamnName[2] << "T" << 1 << 3 << "_" << BigO << ".txt";

             TheDamnName[3] << "T" << 0 << 1 << "_" << BigO << ".txt";
             TheDamnName[4] << "T" << 0 << 2 << "_" << BigO << ".txt";
             TheDamnName[5] << "T" << 1 << 1 << "_" << BigO << ".txt";
             TheDamnName[6] << "T" << 1 << 2 << "_" << BigO << ".txt";
             TheDamnName[7] << "T" << 2 << 2 << "_" << BigO << ".txt";
             TheDamnName[8] << "T" << 2 << 3 << "_" << BigO << ".txt";
             TheDamnName[9] << "T" << 3 << 3 << "_" << BigO << ".txt";
             
             CSName << "EB_" << BigO << ".txt";
          }

          TheFTName[0] << "T" << 0 << 3 << "_" << "FT_ReSq" << ".txt";
          TheFTName[1] << "T" << 0 << 3 << "_" << "FT_ImSq" << ".txt";

          PrintT[0].open(TheDamnName[0].str());
          PrintT[1].open(TheDamnName[1].str());
          PrintT[2].open(TheDamnName[2].str());

          PrintT[3].open(TheDamnName[3].str());
          PrintT[4].open(TheDamnName[4].str());
          PrintT[5].open(TheDamnName[5].str());
          PrintT[6].open(TheDamnName[6].str());
          PrintT[7].open(TheDamnName[7].str());
          PrintT[8].open(TheDamnName[8].str());
          PrintT[9].open(TheDamnName[9].str());
          
          PrintT03FT[0].open(TheFTName[0].str());
          PrintT03FT[1].open(TheFTName[1].str());
          
          PrintCS.open(CSName.str());
                
          for(int x=0; x<ARRAY; x++)
          {
             for(int y=0; y<ARRAY; y++)
             {
                PrintT[0] << T00[BigO][x][y] << " ";
                PrintT[1] << T03[BigO][x][y] << " ";
                PrintT[2] << T13[BigO][x][y] << " ";

                PrintT[3] << T01[BigO][x][y] << " ";
                PrintT[4] << T02[BigO][x][y] << " ";
                PrintT[5] << T11[BigO][x][y] << " ";
                PrintT[6] << T12[BigO][x][y] << " ";
                PrintT[7] << T22[BigO][x][y] << " ";
                PrintT[8] << T23[BigO][x][y] << " ";
                PrintT[9] << T33[BigO][x][y] << " ";
                
                PrintCS   << EBGrid[BigO][x][y] << " ";
             }
             PrintT[0] << endl;
             PrintT[1] << endl;
             PrintT[2] << endl;

             PrintT[3] << endl;
             PrintT[4] << endl;
             PrintT[5] << endl;
             PrintT[6] << endl;
             PrintT[7] << endl;
             PrintT[8] << endl;
             PrintT[9] << endl;
             
             PrintCS << endl;
          }

          for(int x=0; x<(2*size+1); x++)
          {
             for(int y=0; y<(2*size+1); y++)
             {
                PrintT03FT[0] << T03FT[0][x][y] << " ";
                PrintT03FT[1] << T03FT[1][x][y] << " ";
             }
             PrintT03FT[0] << endl;
             PrintT03FT[1] << endl;
          }
                   
          PrintT[0].close();
          TheDamnName[0].str(std::string());
          PrintT[1].close();
          TheDamnName[1].str(std::string());
          PrintT[2].close();
          TheDamnName[2].str(std::string());

          PrintT[3].close();
          TheDamnName[3].str(std::string());
          PrintT[4].close();
          TheDamnName[4].str(std::string());
          PrintT[5].close();
          TheDamnName[5].str(std::string());
          PrintT[6].close();
          TheDamnName[6].str(std::string());
          PrintT[7].close();
          TheDamnName[7].str(std::string());
          PrintT[8].close();
          TheDamnName[8].str(std::string());
          PrintT[9].close();
          TheDamnName[9].str(std::string());
          PrintCS.close();
          CSName.str(std::string());
          
          PrintT03FT[0].close();
          TheFTName[0].str(std::string());
          PrintT03FT[1].close();
          TheFTName[1].str(std::string());
       }
       
cout << "All writted!" << endl;

       if(false)
       {
       RhoPrint.precision(12);
       for(int nuc=0; nuc<2; nuc++)
       {
          for(int col=0;col<8;col++)
          {
             for(int eta=0; eta<etamax; eta++)
             {
                for(int x = 0; x<ARRAY; x++)
                {
                   for(int y = 0; y<ARRAY; y++)
                   {
                      RhoPrint << RhoCG[nuc][col][eta][x][y] << "   ";
                   if(x==y && y==(ARRAY-1)/2)
                   {
                      cout << "Rho for nuc " << nuc << " with col " << col << ": " << RhoCG[nuc][col][eta][x][y] << endl;
                   }
                   }
                   RhoPrint << endl;
                }
             }
          }
       }
       }
/*       
       ofstream EMField[3][2];
       EMField[0][0].open("E02ndRe.txt");
       EMField[1][0].open("E04thRe.txt");
       EMField[2][0].open("Ex3rdRe.txt");
       
       EMField[0][1].open("E02ndIm.txt");
       EMField[1][1].open("E04thIm.txt");
       EMField[2][1].open("Ex3rdIm.txt");
       
       ofstream Zeroth[2][2];
       Zeroth[0][0].open("E0Re.txt");
       Zeroth[0][1].open("E0Im.txt");
       Zeroth[1][0].open("B0Re.txt");
       Zeroth[1][1].open("B0Im.txt");
       
       ofstream AiPrint[2][2];
       AiPrint[0][0].open("AxRe.txt");
       AiPrint[0][1].open("AxIm.txt");
       AiPrint[1][0].open("AyRe.txt");
       AiPrint[1][1].open("AyIm.txt");
       
       for(int x=0; x<ARRAY; x++)
       {
	  cout << "EMF x: " << x << endl;
          for(int y=0; y<ARRAY; y++)
          {
             for(int a=0; a<3; a++)
             {
                for(int b=0; b<3; b++)
                {
                   EMField[0][0] << real(Fij[1][0][0][x][y][a][b]) << "   ";
                   EMField[1][0] << real(Fij[2][0][0][x][y][a][b]) << "   ";
                   EMField[2][0] << real(Fij[2][0][1][x][y][a][b]) << "   ";
                   
                   EMField[0][1] << imag(Fij[1][0][0][x][y][a][b]) << "   ";
                   EMField[1][1] << imag(Fij[2][0][0][x][y][a][b]) << "   ";
                   EMField[2][1] << imag(Fij[2][0][1][x][y][a][b]) << "   ";
                   
                   Zeroth[0][0] << real(Fij[0][0][0][x][y][a][b]) << "   ";
                   Zeroth[0][1] << imag(Fij[0][0][0][x][y][a][b]) << "   ";
                   
                   Zeroth[1][0] << real(Fij[0][1][0][x][y][a][b]) << "   ";
                   Zeroth[1][1] << imag(Fij[0][1][0][x][y][a][b]) << "   ";
                   
                   AiPrint[0][0] << real(Ai[0][0][x][y][a][b]) << "   ";
                   AiPrint[0][1] << imag(Ai[0][0][x][y][a][b]) << "   ";
                   
                   AiPrint[1][0] << real(Ai[0][1][x][y][a][b]) << "   ";
                   AiPrint[1][1] << imag(Ai[0][1][x][y][a][b]) << "   ";
                }
                EMField[0][0] << endl;
                EMField[1][0] << endl;
                EMField[2][0] << endl;
                
                EMField[0][1] << endl;
                EMField[1][1] << endl;
                EMField[2][1] << endl;
                
                Zeroth[0][0] << endl;
                Zeroth[0][1] << endl;
                
                Zeroth[1][0] << endl;
                Zeroth[1][1] << endl;

                AiPrint[0][0] << endl;
                AiPrint[0][1] << endl;
                
                AiPrint[1][0] << endl;
                AiPrint[1][1] << endl;
             }
          }
       }
*/       


//RelAvU();

cout << endl;


runtime= time(0) - runtime;

cout << "Total runtime: " << runtime << endl << endl;

       timer=time(0) - timer;
       
       WriteTime << "Closing up: " << timer << endl;

cout << "~~fin~~";

       WriteTime << endl << "~~fin~~";
    }
