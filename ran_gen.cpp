/*ingw*/
#include<random>
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

using namespace std;

double temp[8192];
double r_sim;
double val1,val2,val3;
int i_sort,i_sign,i_fail;
//int i,j,k;
int iseed;
double cdf_au[2][8192];

std::mt19937_64 engg(std::random_device{}());
double rann() {std::uniform_real_distribution<double> uniran(0.0,1.0); return uniran(engg);}

      void hpsort(unsigned long n, double ra[])
      {
	     unsigned long i,ir,j,l;
	     float rra;

	     if (n < 2) return;
	     l=(n >> 1)+1;
	     ir=n;
	     for (;;) 
         {
		    if (l > 1) 
           {
		      rra=ra[--l];
           }
           else 
           {
		      rra=ra[ir];
			  ra[ir]=ra[1];
			  if (--ir == 1)
              {
                 ra[1]=rra;
				 break;
			  }
		   }
		   i=l;
		   j=l+l;
		   while (j <= ir)
           {
              if (j < ir && ra[j] < ra[j+1]) j++;
			  if (rra < ra[j]) 
              {
                 ra[i]=ra[j];
				 i=j;
				 j <<= 1;
   	          }
              else j=ir+1;
		   }
	       ra[i]=rra;
	    }
      }

      double ran2 () 
      {
         int idum = -1*rann();
         const int IM1 = 2147483563, IM2 = 2147483399;
         const double AM=(1.0/IM1);
         const int IMM1 = IM1-1;
         const int IA1 = 40014, IA2 = 40692, IQ1 = 53668, IQ2 = 52774;
         const int IR1 = 12211, IR2 = 3791, NTAB = 32;
         const int NDIV = 1+IMM1/NTAB;
         const double EPS = 3.0e-16, RNMX = 1.0-EPS;

         int j, k;
         static int idum2=123456789, iy = 0;
         static int iv[NTAB];
         double temp;

         if (idum <= 0)
         {
            idum = (idum == 0 ? 1 : -idum);
            idum2=idum;
            for (j=NTAB+7;j>=0;j--)
            {
               k=idum/IQ1;
               idum=IA1*(idum-k*IQ1)-k*IR1;
               if (idum < 0) idum += IM1;
               if (j < NTAB) iv[j] = idum;
            }
            iy=iv[0];
         }
         k=idum/IQ1;
         idum=IA1*(idum-k*IQ1)-k*IR1;
         if (idum < 0) idum += IM1;
         k=idum2/IQ2;
         idum2=IA2*(idum2-k*IQ2)-k*IR2;
         if (idum2 < 0) idum2 += IM2;
         j=iy/NDIV;
         iy=iv[j]-idum2;
         iv[j] = idum;
         if (iy < 1) iy += IMM1;
         if ((temp=AM*iy) > RNMX)return RNMX;
         else return temp;
      }

      void NucGen(int NN, double** cartisan)
      {
          double ran_sample_r[NN],r_sample[NN],theta_sample[NN],phi_sample[NN];
          
/*          double** cartisan;
          cartisan = new double*[3];
          for(int i=0;i<3;i++)
          {
             cartisan[i] = new double[NN];
             for(int n=0; n<NN; n++)
             {
             cartisan[i][NN] = 0;
             }
          }
*/

          //srand(time(0));          
          ifstream CDFRead;
          if(Atom==1)
          {
              CDFRead.open("CDF_Pb.dat");
          }
          else
          {
              CDFRead.open("CDF_Au.dat");
          }


          for(int j = 0; j<8192; j++)
          {
             for(int i = 0; i<2; i++)
             {
                CDFRead >> cdf_au[i][j];
             }
          }
          
          
          CDFRead.close();
          while(1==0)
          {
             cout << "i?" << endl;
             int k;
             cin >> k;
             cout << "j?" << endl;
             int p;
             cin >> p;
             cout << cdf_au[k][p]<< endl << endl;
          }

          //here then asks for iseed.
          iseed= -1*rann();
          val1=rann();

          bool Succ = false; //Are all Nucleons placed successfully so far?
          while(Succ==false)
          {
             Succ=true; //Well it is true so far! By default we will escape the loop.
             
             for(int i=0; i<NN;i++) // Sample random # from 0 ->1
             {
                val1=rann(); //Rand #
                ran_sample_r[i]=val1; //Rand # Array
             }

             cout << "****************************" << endl;
          
             for(int n=0; n<NN;n++) //Turn sampled # into sampled radius
             {
                for(int i=0; i<8192;i++) //Abs Diff between each CDF point and sampled val- will want to find where 0 is
                {
                   temp[i]=abs(ran_sample_r[n]-cdf_au[1][i]);
                }

             
                i_sort=0;
                double min=temp[0];
                for(int j=1; j<8192; j++)
                {
                   if(min>temp[j]) //Find smallest deviation from tabulated CDF
                   {
                      i_sort = j;
                      min=temp[j];
                   }
                }

                if( (ran_sample_r[n]-cdf_au[1][i_sort]) >0 ) //Are we above or below at min diff?
                {
                   i_sign=1.0;
                }

                else
                {
                   i_sign=-1.0;
                }

                if(i_sign<0) //Interpolate radius for placement
                {
                   r_sample[n]=2.E-3*i_sort+2.E-3*(ran_sample_r[n]-cdf_au[1][i_sort-1])/(cdf_au[1][i_sort]-cdf_au[1][i_sort-1]);
                }
                else
                {
                   r_sample[n]=2.E-3*i_sort+2.E-3*temp[i_sort]/(cdf_au[1][i_sort+1]-cdf_au[1][i_sort]);
                }
             } // All nucleons are given a radius

             hpsort(NN,r_sample);

             ofstream RSampleWrite;
             if(Atom==1)
             {
                 RSampleWrite.open("rsample_Pb.dat");
             }
             else
             {
                 RSampleWrite.open("rsample_Au.dat");
             }
             

             for(int i = 0; i<NN; i++)
             {
                RSampleWrite << r_sample[i] << endl;
             }
          
             RSampleWrite.close();

             val1=rann();
             val1=2.0*(val1-0.5);
             theta_sample[0]=acos(val1);
          
             val1=rann();
             phi_sample[0]=2.0*3.1415926*val1;
          
             cartisan[0][0]=r_sample[0]*sin(theta_sample[0])*cos(phi_sample[0]);
             cartisan[1][0]=r_sample[0]*sin(theta_sample[0])*sin(phi_sample[0]);
             cartisan[2][0]=r_sample[0]*cos(theta_sample[0]);

             i_fail=0;
             for(int i=1;i<NN;i++)
             {
                cout << "i_fail: " << i_fail << endl;
                
                i_fail=0;


                bool Nuc=false; //Has this nucleon beeen placed?
                while(Nuc==false)
                {
                   val1=rann();
                   val1=2.0*(val1-0.5);
                   theta_sample[i]=acos(val1);
//test                   cout << "ran val for theta input: " << val1 << endl;
             
                   val1=rann();
                   phi_sample[i]=2.0*3.1415926*val1;

                   cartisan[0][i]=r_sample[i]*sin(theta_sample[i])*cos(phi_sample[i]);
                   cartisan[1][i]=r_sample[i]*sin(theta_sample[i])*sin(phi_sample[i]);
                   cartisan[2][i]=r_sample[i]*cos(theta_sample[i]);
             
                   for(int j=0;j<i;j++)
                   {
                      if(abs(r_sample[i]-r_sample[j])<r_core)
                      {
                         r_sim=sqrt((cartisan[0][j]-cartisan[0][i])*(cartisan[0][j]-cartisan[0][i])+
                                    (cartisan[1][j]-cartisan[1][i])*(cartisan[1][j]-cartisan[1][i])+
                                    (cartisan[2][j]-cartisan[2][i])*(cartisan[2][j]-cartisan[2][i])
                                   );
                   
                         if(r_sim < r_core)
                         {
                            i_fail=i_fail+1; // Count Fails
//                            cout << "nucleon index = " << i << endl;
//                            cout << "i_fail = " << i_fail;
                         }
                         else{Nuc=true;}
                      }
                      else{Nuc=true;}
                   }
                   if(i_fail==20)
                   {
                      Nuc=true; //Will escape Nuc Loop in shame
                      i=2*NN; // Escape the i loop, don't bother placing more nuclei
                      Succ=false; //Do NOT escape the Succ loop! Restart the process with radii sampling
                   }
                }//While Nuc
             }//i loop
          }//while succ

          ofstream SampleWrite;
          if(Atom==1)
          {
              SampleWrite.open("sample_Pb.dat");
          }
          else
          {
              SampleWrite.open("sample_Au.dat");
          }

          for(int i = 0; i<NN; i++)
          {
             SampleWrite << cartisan[0][i] << "  " << cartisan[1][i] << "  " << cartisan[2][i] << endl;
          }
          
          SampleWrite.close();
          
          ofstream BubbleWrite;
          if(Atom==1)
          {
              BubbleWrite.open("sample_Pb_bubble.dat");
          }
          else
          {
              BubbleWrite.open("sample_Au_bubble.dat");
          }

          for(int i = 0; i<NN; i++)
          {
             BubbleWrite << cartisan[0][i] << "  " << cartisan[1][i] << "  " << cartisan[2][i] << "  " << 0.7 << endl;
          }
          
          BubbleWrite.close();
      } 
