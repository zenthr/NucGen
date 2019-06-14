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
#include <complex>
#include <cmath>
#include <fftw3.h>

using namespace std;

int main()
{

int N = 125;
//   static double Input[125];
//   static double Output[125];



fftw_complex *in, *out;
fftw_plan p;
//...
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

for(int i=0; i<N; i++)
{

      in[i][0] = 0.;
      in[i][1] = 0;
}

in[1][0] = 10.;
//Input[0] = 0.;
fftw_execute(p); /*repeat as needed*/
fftw_destroy_plan(p);

ofstream Sin;
Sin.open("Pride.txt");
for(int i=0; i<N; i++)
{
   Sin << out[i][0]/(125*125) << endl;
//   Sin << Output[i] << endl;
}

Sin.close();

Sin.open("PrideIm.txt");
for(int i=0; i<N; i++)
{
   Sin << out[i][1] << endl;
}

Sin.close();


in[0][0] = 0;
//Input[0] = 0;

fftw_free(in); fftw_free(out);

cout << "Pride" << endl;


/********************************************************************/


in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

for(int i=0; i<N; i++)
{
   for(int j=0; j<2; j++)
   {
      in[i][j] = 0;
//      Input[i] = 0;
   }
}

in[2][0] = 10.;
//Input[2] = 10.;
fftw_execute(p); /*repeat as needed*/
fftw_destroy_plan(p);

Sin.open("Greed.txt");
for(int i=0; i<N; i++)
{
   Sin << out[i][0] << endl;
//   Sin << Output[i] << endl;
}

Sin.close();

Sin.open("GreedIm.txt");
for(int i=0; i<N; i++)
{
   Sin << out[i][1] << endl;
}

Sin.close();

in[3][0] = 0;
//Input[3] = 0;

fftw_free(in); fftw_free(out);

cout << "Greed" << endl;


/********************************************************************/


in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

for(int i=0; i<N; i++)
{
   for(int j=0; j<2; j++)
   {
      in[i][j] = 0;
//      Input[i] = 0;
   }
}

in[62][0] = 10.;
//Input[62] = 10.;
fftw_execute(p); /*repeat as needed*/
fftw_destroy_plan(p);

Sin.open("Sloth.txt");
for(int i=0; i<N; i++)
{
   Sin << out[i][0] << endl;
//   Sin << Output[i] << endl;
}

Sin.close();

Sin.open("SlothIm.txt");
for(int i=0; i<N; i++)
{
   Sin << out[i][1] << endl;
}

Sin.close();

in[62][0] = 0;
//Input[62] = 0;

fftw_free(in); fftw_free(out);

cout << "Sloth" << endl;


/********************************************************************/


in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

for(int i=0; i<N; i++)
{
   for(int j=0; j<2; j++)
   {
      in[i][j] = 0;
//      Input[i] = 0;
   }
}

in[123][0] = 10.;
//Input[123] = 10.;
fftw_execute(p); /*repeat as needed*/
fftw_destroy_plan(p);

Sin.open("Envy.txt");
for(int i=0; i<N; i++)
{
   Sin << out[i][0] << endl;
//   Sin << Output[i] << endl;
}

Sin.close();

Sin.open("EnvyIm.txt");
for(int i=0; i<N; i++)
{
   Sin << out[i][1] << endl;
}

Sin.close();

in[124][0] = 0;
//Input[124] = 0;

fftw_free(in); fftw_free(out);

cout << "Envy" << endl;


/********************************************************************/

N=63;

in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

for(int i=0; i<N; i++)
{
   for(int j=0; j<2; j++)
   {
      in[i][j] = 0;
//      Input[i] = 0;
   }
}


in[62][0] = 10.;
//Input[62] = 10.;
fftw_execute(p); /*repeat as needed*/
fftw_destroy_plan(p);

Sin.open("Wrath.txt");
for(int i=0; i<N; i++)
{
   Sin << out[i][0] << " " << endl;
//   Sin << Output[i] << " " << endl;
}

Sin.close();

Sin.open("WrathIm.txt");
for(int i=0; i<N; i++)
{
   Sin << out[i][1] << " " << endl;
}

Sin.close();

in[62][0] = 0;
//Input[62] = 0;

fftw_free(in); fftw_free(out);

cout << "Wrath" << endl;


/********************************************************************/

N=125;

in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

for(int i=0; i<N; i++)
{
   for(int j=0; j<2; j++)
   {
      in[i][j] = 0;
//      Input[i] = 0;
   }
}

in[100][0] = 5.;
//Input[100] = 10.;
fftw_execute(p); /*repeat as needed*/
fftw_destroy_plan(p);

p = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

fftw_execute(p); /*repeat as needed*/
fftw_destroy_plan(p);

Sin.open("Lust.txt");
for(int i=0; i<N; i++)
{
   Sin << in[i][0] << endl;
//   Sin << Input[i] << endl;
}

Sin.close();

Sin.open("LustIm.txt");
for(int i=0; i<N; i++)
{
   Sin << in[i][1] << endl;
}

Sin.close();

in[100][0] = 0;
//Input[100][0] = 0;

fftw_free(in); fftw_free(out);

cout << "Lust" << endl;


/********************************************************************/

fftw_complex *in2, *out2;

in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N); // what the fuck is the 1d bullshit...
out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N);

p = fftw_plan_dft_2d(N, N, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);

for(int i=0; i<N*N; i++)
{
   for(int j=0; j<2; j++)
   {
      in2[i][j] = 0;
//      Input[i] = 0;
   }
}

//for(int j=0; j<N; j++)
//{
//   int i=5;
//   in2[i*N + j][0] = 10.;
//}

in2[0*N + 5][0] = 10.;
//Input[123] = 10.;
fftw_execute(p); /*repeat as needed*/

fftw_destroy_plan(p);

Sin.open("Gluttony.txt");
for(int i=0; i<N; i++)
{
   for(int j=0; j<N; j++)
   {
      Sin << out2[i*N + j][0] << "  ";
   }
   Sin << endl;
//   Sin << Output[i] << endl;
}

Sin.close();

Sin.open("GluttonyIm.txt");
for(int i=0; i<N; i++)
{
   for(int j=0; j<N; j++)
   {
      Sin << out2[i*N + j][1] << "  ";
   }
   Sin << endl;
}

Sin.close();

Sin.open("GluttonyIn.txt");
for(int i=0; i<N; i++)
{
   for(int j=0; j<N; j++)
   {
      Sin << in2[i*N + j][0] << "  ";
   }
   Sin << endl;
//   Sin << Output[i] << endl;
}

Sin.close();

fftw_free(in2); fftw_free(out2);

cout << "Gluttony" << endl;
 
   printf("Hello cruel world\n");
   
//   cout << Input[62] << endl;
   
   cout << "Goodbye, cruel world..." << endl;
   printf("... b a n g ....");
}
