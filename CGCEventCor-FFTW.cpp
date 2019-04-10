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
    double Hist[16][200];
    double Ecc[2][200];
    double EBHist[4][200];
    double CSCharge[200][3];
    double EBsq[100][4];
    double RhoStep;
    double***** CovPot;
    static double*** AlphaCor;


      static double**** MuShapeCenter;

       static int Numb;

//    static complex<double>*** AScale;

    complex<double>**** AT;
    complex<double>**** AS;

    complex<double>***** AvU;
    
    static complex<double> Bx[3][3][3];

std::mt19937_64 eng(std::random_device{}());
//std::mt19937_64 eng(145709713); // run with set seed: 145709713
double ran() {std::uniform_real_distribution<double> uniran(0.0,1.0); return uniran(eng);}


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
              
              MuShape[0][x][y]=NucScale*Shape0;
              MuShape[1][x][y]=NucScale*Shape1;
              
              MuShapeCenter[0][0][x][y] = MuShapeCenter[0][0][x][y] + MuShape[0][x][y];
              MuShapeCenter[1][0][x][y] = MuShapeCenter[1][0][x][y] + MuShape[1][x][y];
           }
       }
       
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
    
    int main()
    {
       DeclDelMe();
       
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
       
       timer=time(0);

       DeclAlphaCor();

       ofstream RhoPrint;

       srand(time(0));

       RhoStep = 2*RhoMax/(NUMSTEP-1);
       CDFTable();

       cout << CDFTab[780][0] << endl;
       cout << CDFTab[780][1] << endl;

       for(int n=0; n<N; n++)
       {
       if(NucMethod == 0)
       {
          WSIntegration();
       }
       
       if(NucMethod == 1)
       {
          NucSampler();
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
       
       DeclCovPot();

          cout << "Sample'd: " << endl;
          cout << "Runtime : " << time(0) - timer << endl;
          cout << "Event No: " << n << endl;
       for(int nuc=0; nuc<2; nuc++)
       {
       for(int eta=0; eta<etamax; eta++)
       {

/*** The Rho Sampler ***/

       DeclRho();

       for(int x = 0; x<ARRAY; x++)
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
                   }
                   //Scale by Sqrt(mu) to scale strength
                   //Scale by hCG to account for coarse graining (1/hCG) and include summation element for sum in alpha (hCG*hCG)
                   RhoCG[nuc][col][eta][x][y] = RhoCG[nuc][col][eta][x][y]*sqrt(MuShape[nuc][x][y])*g*h;
               }
           }//close y loop
       }//close x loop

   /*** Solving Poisson ***/


       for(int col = 0; col<8; col++)
       {
          cout << "Potential for nuc = " << nuc << endl;

          if(nuc==0 && col==0 && eta==0)
          {
              for(int x = 0; x<ARRAY; x++)
              {
                 for(int y = 0; y<ARRAY; y++)
                 {
                    CovPot[nuc][col][eta][x][y] = 0.; //to be safe
                 }
              }
          }
          
          
          
          fftw_complex *rhoalpha; //Will contain rho(x), rho(k), alpha(k), alpha(x)
          fftw_plan p, q; //Is how FFTW knows what to do


          rhoalpha = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ARRAY*ARRAY); // This is a linear list of ARRAY*ARRAY size- it's what FFTW expects


          for(int i=0; i<ARRAY; i++) //real part
          {
             for(int j=0; j<ARRAY; j++) // imaginary part
             {
                rhoalpha[i*ARRAY + j][0] = RhoCG[nuc][col][eta][i][j]; //rho(x)
                rhoalpha[i*ARRAY + j][1] = 0;
             }
          }


          p = fftw_plan_dft_2d(ARRAY, ARRAY, rhoalpha, rhoalpha, FFTW_FORWARD, FFTW_ESTIMATE); //Forward transformation, overwriting rhoalpha
          q = fftw_plan_dft_2d(ARRAY, ARRAY, rhoalpha, rhoalpha, FFTW_BACKWARD, FFTW_ESTIMATE); //Backward transformation, overwriting rhoalpha



          fftw_execute(p); //rhoalpha -> rho(k) - It's just that easy!

          double mterm = mIR*mIR*h*h; //Mass term
                    
          for(int i=0; i<ARRAY; i++) //rhoalpha -> alpha(k)
          {
             double ki = 2*PIE*(min(i,ARRAY-i))/(ARRAY);
              for(int j=0; j<ARRAY; j++)
             {

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


          for(int i=0; i<ARRAY; i++)
          {
             for(int j=0; j<ARRAY; j++)
             {
                CovPot[nuc][col][eta][i][j] = rhoalpha[i*ARRAY + j][0]/(ARRAY*ARRAY); //Forward then backward transform is renormalized
             }
          }

          fftw_free(rhoalpha); // Discard variable
          
          fftw_destroy_plan(p); // Discard plans
          fftw_destroy_plan(q);
}//end col loop

/*** This is what my correlation calculations look like ***/

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




RelRho();
          } //end eta loop
          } //end nuc loop

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
      
      ofstream PrintAlphaCor[2];
      PrintAlphaCor[0].open("AlphaCorCent.txt");
      PrintAlphaCor[1].open("AlphaCorOff.txt");

       for(int x = 0; x<ARRAY; x++)
       {
          for(int y = 0; y<ARRAY; y++)
          {             
             PrintAlphaCor[0] << AlphaCor[0][x][y] << "  ";
             PrintAlphaCor[1] << AlphaCor[1][x][y] << "  ";
          }
          PrintAlphaCor[0] << endl;
          PrintAlphaCor[1] << endl;
       }
       
       
cout << "All writted!" << endl;

cout << endl;

cout << "Final time: " << time(0) - timer << endl;

cout << "~~fin~~";
    }
