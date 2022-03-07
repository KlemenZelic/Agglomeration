#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**************************************************************************************************/
/**Parameters**/

#define a       200.0       //Primary Particle principal semi - axis in x direction (nm)
#define b       150.0       //Primary Particle principal semi - axis in y direction (nm)
#define c       100.0       //Primary Particle principal semi - axis in z direction (nm)
#define NP      1000         //Number of particle in agglomerate
#define R2      0.5           //ratio between agglomerate semi-axis in direction y and x (r2/r1)
#define R3      0.5           //ratio between agglomerate semi-axis in direction z and x (r3/r1)
#define S       0.4         //particle sizes dispersion (-)
#define TH      1.5         //Primary particles are randomly rotated around random direction for angles between -TH and TH (-)
#define Dxx     1e-17       //diffusion tensor component 11 in particle crystal lattice coordinate system
#define Dyy     1e-15       //diffusion tensor component 22 in particle crystal lattice coordinate system
#define Dzz     1e-17       //diffusion tensor component 33 in particle crystal lattice coordinate system

                            //All other components of diffusion coefficient tensor equal 0 in self coordinate system

/**************************************************************************************************/
/**Functions**/

/*funcrion that compare two double values by magnitude. Serving as a tool for sorting the lists from higher to lower values*/
int cmpfunc (const void * aa, const void * bb)
{
  if (*(double*)aa < *(double*)bb)
    return 1;
  else if (*(double*)aa > *(double*)bb)
    return -1;
  else
    return 0;
}

/*Same function as previous one but with the inverted inequality sign. Serving as a tool for sorting the lists from lower to higher values */
int cmpfunc2 (const void * aa, const void * bb)
{
  if (*(double*)aa > *(double*)bb)
    return 1;
  else if (*(double*)aa < *(double*)bb)
    return -1;
  else
    return 0;
}

/*Function that compares tow lists and tels if they equal*/
int compareArrays(double array1[], double array2[], int N) {
  int ii;
  for(ii = 1; ii <= N; ii++) {
    if (array1[ii] != array2[ii]) return 0;
  }
  return 1;
}

/*Function that tests of value A is an element of interval (min,max)*/
int Outside_interval(double A, double min, double max){
if(A<min) return 1;
else if (A>max) return 1;
else return -1;
}

/*Nine functions that calculate components of rotational matrix from normalized axial vector and rotational angle*/

double R11(double PhiRot, double ux, double uy, double uz){
return cos(PhiRot) + ux*ux*(1 - cos(PhiRot));
}

double R12(double PhiRot, double ux, double uy, double uz){
return ux*uy*(1 - cos(PhiRot)) - uz*sin(PhiRot);
}

double R13(double PhiRot, double ux, double uy, double uz){
return ux*uz*(1 - cos(PhiRot)) + uy*sin(PhiRot);
}

double R21(double PhiRot, double ux, double uy, double uz){
return uy*ux*(1 - cos(PhiRot)) + uz*sin(PhiRot);
}

double R22(double PhiRot, double ux, double uy, double uz){
return cos(PhiRot) + uy*uy*(1 - cos(PhiRot));
}

double R23(double PhiRot, double ux, double uy, double uz){
return uy*uz*(1 - cos(PhiRot)) - ux*sin(PhiRot);
}

double R31(double PhiRot, double ux, double uy, double uz){
return uz*ux*(1 - cos(PhiRot)) - uy*sin(PhiRot);
}

double R32(double PhiRot, double ux, double uy, double uz){
return uz*uy*(1 - cos(PhiRot)) + ux*sin(PhiRot);
}

double R33(double PhiRot, double ux, double uy, double uz){
return cos(PhiRot) + uz*uz*(1 - cos(PhiRot));
}

/*****************************************************************************************/
/**main**/

int main()
{

/*files for exporting data*/
FILE *data;
FILE *parameters;

/*Datoteke za uvažanje*/
FILE *SurfaceAreas12;
FILE *NearesNeighbours12;
FILE *SurfaceAreas3;
FILE *NearesNeighbours3;

int M = 20; //maximal possible size of agglomerate in a1, a2 and a3 direction
//ARBITRARY SET VALUE FOR THE MAXIMAL NUMBER OF PARTICLES IN DIRECTIONS a1,a2 AND a3.

srand((unsigned int)time(NULL));

//Definitions of base vector of body centered orthorombic Bravais lattice.
double a1[] = {a, b, -c};
double a2[] = {-a, b, c};
double a3[] = {a, -b, c};
printf("%f, %f, %f\n", a1[1], a2[1], a3[1]);
//Initialization of different variables
double RAn;         //particle Rotation angle (between -TH and TH)
double RAx[3];      //particle Rotation axis (normal vector along axis of rotation)
double CosPhi[3];   //particle rotation axis vector non-normalized components
double NormRAx;     //Normalization of rotation axis vector
double BufV[3];     //Particle location (vector to the center of particle)
double r1;          //agglomerate semi-axis in x direction
double r2;          //agglomerate semi-axis in y direction
double r3;          //agglomerate semi-axis in z direction
double Vp;          //average sized particle volume calculated as (a1 cross a2).a3

Vp = 2.0*fabs((a1[1]*a2[2]+a1[2]*a2[1])*a3[0] + (a1[2]*a2[0]+a1[0]*a2[2])*a3[1] + (a1[0]*a2[1]+a1[1]*a2[0])*a3[2]); //Particle volume defined as mixed product of
r1 = ceil(cbrt(3.0*NP*Vp/(4.0*3.14159*R2*R3))/a)*a; //Agglomerate semi-axis.



r2 = R2*r1;     //r2 in r3 calculated from R2 and R3 ratios
r3 = R3*r1;

double AggData[2*(M*M*M) + 3][9]; //Initialization of array to store results - length of an array is maximal possible number of particle

//iterators that will run in the directions of base vectors a1, a2 and a3 and N of particles
int N = 0;

/*****Agglomerate construction******/

for(int i=-M; i<M; i++){
    for(int j=-M; j<M; j++){
        for(int k=-M; k<M; k++){
            BufV[0] = i*a1[0]+j*a2[0]+k*a3[0];
            BufV[1] = i*a1[1]+j*a2[1]+k*a3[1];
            BufV[2] = i*a1[2]+j*a2[2]+k*a3[2];
            if((BufV[0]*BufV[0])/(r1*r1) + (BufV[1]*BufV[1])/(r2*r2) + (BufV[2]*BufV[2])/(r3*r3) <= 1.0){
                RAn = ((double)rand()/(double)(RAND_MAX)) * 2*TH - TH;
                CosPhi[0] = cos(((double)rand()/(double)(RAND_MAX)) * 2*3.14159 - 3.14159);
                CosPhi[1] = cos(((double)rand()/(double)(RAND_MAX)) * 2*3.14159 - 3.14159);
                CosPhi[2] = cos(((double)rand()/(double)(RAND_MAX)) * 2*3.14159 - 3.14159);
                NormRAx = sqrt(CosPhi[0]*CosPhi[0] + CosPhi[1]*CosPhi[1] + CosPhi[2]*CosPhi[2]);
                RAx[0] = CosPhi[0]/NormRAx;
                RAx[2] = CosPhi[2]/NormRAx;

                AggData[N][0] = sqrt(BufV[0]*BufV[0]/(a*a) + BufV[1]*BufV[1]/(b*b) + BufV[2]*BufV[2]/(c*c));
                AggData[N][1] = BufV[0]/a;
                AggData[N][2] = BufV[1]/b;
                AggData[N][3] = BufV[2]/c;
                AggData[N][4] = 1.0 + ((double)rand()/(double)(RAND_MAX)) * 2*S - S;
                AggData[N][5] = RAn;
                AggData[N][6] = RAx[0];
                AggData[N][7] = RAx[1];
                AggData[N][8] = RAx[2];
                N = N + 1;
            }
        }
    }
}

printf("Number of particles in aglomerate 1: %i\n", N);

if(N<NP){
    r1=r1+a;
    r2=r1*R2;
    r3=r1*R3;

    N = 0;

    for(int i=-M; i<M; i++){
        for(int j=-M; j<M; j++){
            for(int k=-M; k<M; k++){
                BufV[0] = i*a1[0]+j*a2[0]+k*a3[0];
                BufV[1] = i*a1[1]+j*a2[1]+k*a3[1];
                BufV[2] = i*a1[2]+j*a2[2]+k*a3[2];
                if((BufV[0]*BufV[0])/(r1*r1) + (BufV[1]*BufV[1])/(r2*r2) + (BufV[2]*BufV[2])/(r3*r3) <= 1.0){
                    RAn = ((double)rand()/(double)(RAND_MAX)) * 2*TH - TH;
                    CosPhi[0] = cos(((double)rand()/(double)(RAND_MAX)) * 2*3.14159 - 3.14159);
                    CosPhi[1] = cos(((double)rand()/(double)(RAND_MAX)) * 2*3.14159 - 3.14159);
                    CosPhi[2] = cos(((double)rand()/(double)(RAND_MAX)) * 2*3.14159 - 3.14159);
                    NormRAx = sqrt(CosPhi[0]*CosPhi[0] + CosPhi[1]*CosPhi[1] + CosPhi[2]*CosPhi[2]);
                    RAx[0] = CosPhi[0]/NormRAx;
                    RAx[1] = CosPhi[1]/NormRAx;
                    RAx[2] = CosPhi[2]/NormRAx;

                    AggData[N][0] = sqrt(BufV[0]*BufV[0]/(a*a) + BufV[1]*BufV[1]/(b*b) + BufV[2]*BufV[2]/(c*c));
                    AggData[N][1] = BufV[0]/a;
                    AggData[N][2] = BufV[1]/b;
                    AggData[N][3] = BufV[2]/c;
                    AggData[N][4] = 1.0 + ((double)rand()/(double)(RAND_MAX)) * 2*S - S;
                    AggData[N][5] = RAn;
                    AggData[N][6] = RAx[0];
                    AggData[N][7] = RAx[1];
                    AggData[N][8] = RAx[2];
                    N = N + 1;
                }
            }
        }
    }
}

qsort(AggData, N, sizeof(*AggData), cmpfunc2);

N=NP;

int Position[N][3];
double Orientation[N][4];
double Size[N];

//Printing particle parameters to file
parameters = fopen("parameters.txt","w");
fprintf(parameters, "%f %f %f", a, b, c);
fclose(parameters);

//Printing results to file
data = fopen("data.txt","w");

for(int i=0; i<N; i++){
    fprintf(data, "%f %f %f %f %f %f %f %f\n",AggData[i][1], AggData[i][2], AggData[i][3], AggData[i][4], AggData[i][5], AggData[i][6], AggData[i][7], AggData[i][8]);

    Position[i][0] = (int)AggData[i][1];
    Position[i][1] = (int)AggData[i][2];
    Position[i][2] = (int)AggData[i][3];
    Size[i] = AggData[i][4];
    Orientation[i][0] = AggData[i][5];
    Orientation[i][1] = AggData[i][6];
    Orientation[i][2] = AggData[i][7];
    Orientation[i][3] = AggData[i][8];
}
fclose(data);

printf("sizeof position: %i\n", sizeof(Position)/sizeof(Position[0]));
/********Look up tables for areas and neighbours*********************/
//Routine here needs look up tables for the sizes of interfaces, that were previously calculated in order to boost computational time. Look up tables are atached
//SurfaceAreas1.txt, NearestNeighbours1.txt, SurfaceAreas31.txt and NearestNeighbours31.txt


double ParamBuf[]={a,b,c};
qsort(ParamBuf, 3, sizeof(double), cmpfunc);
double NormParam[]={1,ParamBuf[1]/ParamBuf[0],ParamBuf[2]/ParamBuf[0]};

double Permutations[6][3]={{a,b,c},{a,c,b},{b,a,c},{b,c,a},{c,a,b},{c,b,a}};

int Cyc[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
int CC=0;

for(CC=0; CC<6; CC++){
    if(compareArrays(ParamBuf, Permutations[CC], 2)){
        break;
    }
}

SurfaceAreas12 = fopen("SurfaceAreas1.txt", "r");
NearesNeighbours12 = fopen("NearestNeighbours1.txt", "r");
SurfaceAreas3 = fopen("SurfaceAreas31.txt", "r");
NearesNeighbours3 = fopen("NearestNeighbours31.txt", "r");

double SizeRatios[2];
double SABuf12[14];
double NNBuf12[42];
double SABuf3[12];
double NNBuf3[36];
char line[500];
char line1[80000];
int LineIndex = 0;

while(fgets(line, sizeof(line), SurfaceAreas12) != NULL){      //searching for the last line of file

    for(int p=0; p<16; p++){
        if(p<2)
        fscanf(SurfaceAreas12, "%lf", &SizeRatios[p]);
        else
        fscanf(SurfaceAreas12, "%lf", &SABuf12[p-2]);
        }
    if(SizeRatios[0]>=NormParam[1] && SizeRatios[1]>=NormParam[2])
        break;
    LineIndex++;
}

int LineIndex1=0;
while(fgets(line1, sizeof(line1), NearesNeighbours12) != NULL){      //searching for the last line of file

    for(int r=0; r<44; r++){
            fscanf(NearesNeighbours12, "%lf", &NNBuf12[r-2]);
        }
    if(LineIndex1==LineIndex)
        break;
    LineIndex1++;
}

int LineIndex2=0;
while(fgets(line, sizeof(line), SurfaceAreas3) != NULL){      //searching for the last line of file

    for(int p=0; p<14; p++){
        fscanf(SurfaceAreas3, "%lf", &SABuf3[p-2]);
        }
    if(LineIndex2 == LineIndex){
        break;
    }
    LineIndex2++;
}


int LineIndex3=0;
while(fgets(line1, sizeof(line1), NearesNeighbours3) != NULL){      //searching for the last line of file
    for(int r=0; r<38; r++){
        fscanf(NearesNeighbours3, "%lf", &NNBuf3[r-2]);                  //reading last line from Result.txt and saving it to vector x
        }
    if(LineIndex3==LineIndex){
        break;
    }
    LineIndex3++;
}

fclose(SurfaceAreas12);
fclose(NearesNeighbours12);
fclose(SurfaceAreas3);
fclose(NearesNeighbours3);

/*
for(int s=0; s<14; s++){
        printf("%f %f %f         %f\n", NNBuf12[3*s],  NNBuf12[3*s+1],  NNBuf12[3*s+2], SABuf12[s]);  //FOR DEBUGGING
}

for(int s=0; s<12; s++){
        printf("%f %f %f         %f\n", NNBuf3[3*s],  NNBuf3[3*s+1],  NNBuf3[3*s+2], SABuf3[s]);   //FOR DEBUGGING
}
*/
/**************Normalization of imported data to the specific case****************/

int NNN = 0;
double PSA=0;
double PSA3=0;
double AvgS=0;

for(int i=0; i<N; i++){
    AvgS=AvgS+Size[i];
}
AvgS=AvgS/N;
printf("%f\n", AvgS);

for(int i=0; i<14; i++){
    PSA = PSA + ParamBuf[0]*AvgS*ParamBuf[0]*AvgS*SABuf12[i];
    PSA3 = PSA3 + ParamBuf[0]*AvgS*ParamBuf[0]*AvgS*SABuf3[i];
    if(SABuf12[i]>0){
        NNN = NNN + 1;
    }
}

int NN12[NNN][3];
double SA12[NNN];
int NN3[12][3];
double SA3[12];

for(int i=0; i<NNN; i++){
        for(int j=0; j<3; j++){
            NN12[i][(Cyc[CC][j])]=(int)(NNBuf12[3*i+j]/NormParam[j]);
        }
        SA12[i] = ParamBuf[0]*ParamBuf[0]*SABuf12[i]/PSA;
}

for(int i=0; i<12; i++){
        for(int j=0; j<3; j++){
            NN3[i][(Cyc[CC][j])]=(int)(NNBuf3[3*i+j]/NormParam[j]);
        }
        SA3[i] = ParamBuf[0]*ParamBuf[0]*SABuf3[i]/PSA3;
}


/*
printf("\n Nearest Neighbors 3index                  Normalized surface areas\n"); //FOR DEBUGGING
for(int s=0; s<NNN; s++){
        printf("    %d %d %d                                     %f\n", NN12[s][0],  NN12[s][1],  NN12[s][2], SA12[s]); //FOR DEBUGGING
}


printf("\n\n Nearest Neighbors 3index                  Normalized surface areas\n");
for(int s=0; s<12; s++){
        printf("    %d %d %d                                     %f\n", NN3[s][0],  NN3[s][1],  NN3[s][2], SA3[s]); //FOR DEBUGGING
}
*/
/***********Sorting neighbours and surfacces in arays***************/


int Neighbours12[N][NNN];
int Neighbours3[N][12];

for(int i=0; i<N; i++){
    for(int j=0; j<NNN; j++){
        Neighbours12[i][j]=-1;
    }
    for(int j=0; j<12; j++){
        Neighbours3[i][j]=-1;
    }
}

for(int i=0; i<N;i++){
    for(int j=0; j<NNN; j++){
        for(int k=0; k<N; k++){
            if(Position[i][0] + NN12[j][0] == Position[k][0] && Position[i][1] + NN12[j][1] == Position[k][1] && Position[i][2] + NN12[j][2] == Position[k][2] && SA12[j]>0.0){
                Neighbours12[i][j]=k;
                break;
            }
        }
    }
    for(int j=0; j<12; j++){
        for(int k=0; k<N; k++){
            if(Position[i][0] + NN3[j][0] == Position[k][0] && Position[i][1] + NN3[j][1] == Position[k][1] && Position[i][2] + NN3[j][2] == Position[k][2]){
                Neighbours3[i][j]=k;
                break;
            }
        }
    }
}

double A12ij[N][NNN];
double A3ij[N][12];

for(int i=0; i<N; i++){
    for(int j=0; j<NNN; j++){
        if(Neighbours12[i][j]<0){
            A12ij[i][j]=SA12[j]*Size[i]*Size[i];
        }
        else{
            A12ij[i][j] = pow(fmin(Size[i],Size[Neighbours12[i][j]]),2)*SA12[j];
        }
    }
}

for(int i=0; i<N; i++){
    for(int j=0; j<12; j++){
        if(Neighbours3[i][j]<0){
            A3ij[i][j]=pow(Size[i],2)*SA3[j];
        }
        else{
            A3ij[i][j] = pow(fmin(Size[i],Size[Neighbours3[i][j]]),2)*SA3[j];
        }
    }
}

/********************Permeability*************/

double D_coeff[3] = {Dxx,Dyy,Dzz};
qsort(D_coeff, 3, sizeof(double), cmpfunc);

double Da = Dxx/D_coeff[0];
double Db = Dyy/D_coeff[0];
double Dc = Dzz/D_coeff[0];


double Dij12;
double Dji12;
double Dij3;
double Dji3;
double d;
double nx;
double ny;
double nz;
double P12[N][NNN];
double P3[N][12];


for(int i = 0; i<N; i++){
    for(int j = 0; j<NNN; j++){

        d = sqrt(NN12[j][0]*a*NN12[j][0]*a + NN12[j][1]*b*NN12[j][1]*b + NN12[j][2]*c*NN12[j][2]*c);

        nx = NN12[j][0]*a/d;
        ny = NN12[j][1]*b/d;
        nz = NN12[j][2]*c/d;

     Dij12 = nx*(nx*(Da*R11(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R11(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R21(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R21(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R31(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R31(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])) +

                 ny*(Da*R11(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R12(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R21(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R22(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R31(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R32(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])) +

                 nz*(Da*R11(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R13(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R21(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R23(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R31(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R33(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]))) +

             ny*(nx*(Da*R11(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R12(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R21(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R22(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R31(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R32(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])) +

                 ny*(Da*R12(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R12(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R22(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R22(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R32(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R32(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])) +

                 nz*(Da*R12(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R13(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R22(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R23(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R32(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R33(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]))) +

             nz*(nx*(Da*R11(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R13(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R21(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R23(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R31(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R33(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])) +

                 ny*(Da*R12(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R13(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R22(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R23(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R32(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R33(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])) +

                 nz*(Da*R13(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R13(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R23(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R23(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R33(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R33(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])));


        if(Neighbours12[i][j] >= 0){
      Dji12 = nx*(nx*(Da*R11(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R11(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Db*R21(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R21(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Dc*R31(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R31(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])) +

                 ny*(Da*R11(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R12(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Db*R21(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R22(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Dc*R31(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R32(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])) +

                 nz*(Da*R11(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R13(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Db*R21(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R23(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Dc*R31(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R33(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]))) +

             ny*(nx*(Da*R11(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R12(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Db*R21(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R22(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Dc*R31(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R32(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])) +

                 ny*(Da*R12(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R12(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Db*R22(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R22(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Dc*R32(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R32(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])) +

                 nz*(Da*R12(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R13(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Db*R22(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R23(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Dc*R32(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R33(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]))) +

             nz*(nx*(Da*R11(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R13(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Db*R21(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R23(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Dc*R31(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R33(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])) +

                 ny*(Da*R12(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R13(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Db*R22(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R23(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Dc*R32(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R33(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])) +

                 nz*(Da*R13(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R13(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Db*R23(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R23(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3]) +
                     Dc*R33(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])*R33(Orientation[Neighbours12[i][j]][0],Orientation[Neighbours12[i][j]][1],Orientation[Neighbours12[i][j]][2],Orientation[Neighbours12[i][j]][3])));


            P12[i][j] = D_coeff[0]/(0.5*d*((1/Dij12)+(1/Dji12))*1e-9);

        }

        else{
            P12[i][j] = D_coeff[0]/(0.5*d*(1/Dij12)*1e-9);
        }
    }
}

for(int i = 0; i<N; i++){
    for(int j = 0; j<12; j++){

        d = sqrt(NN3[j][0]*a*NN3[j][0]*a + NN3[j][1]*b*NN3[j][1]*b + NN3[j][2]*c*NN3[j][2]*c);
        nx = NN3[j][0]*a/d;
        ny = NN3[j][1]*b/d;
        nz = NN3[j][2]*c/d;


     Dij3 = nx*(nx*(Da*R11(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R11(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R21(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R21(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R31(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R31(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])) +

                 ny*(Da*R11(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R12(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R21(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R22(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R31(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R32(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])) +

                 nz*(Da*R11(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R13(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R21(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R23(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R31(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R33(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]))) +

             ny*(nx*(Da*R11(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R12(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R21(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R22(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R31(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R32(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])) +

                 ny*(Da*R12(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R12(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R22(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R22(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R32(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R32(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])) +

                 nz*(Da*R12(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R13(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R22(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R23(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R32(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R33(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]))) +

             nz*(nx*(Da*R11(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R13(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R21(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R23(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R31(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R33(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])) +

                 ny*(Da*R12(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R13(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R22(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R23(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R32(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R33(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])) +

                 nz*(Da*R13(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R13(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Db*R23(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R23(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3]) +
                     Dc*R33(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])*R33(Orientation[i][0],Orientation[i][1],Orientation[i][2],Orientation[i][3])));


        if(Neighbours3[i][j] >= 0){
      Dji3 = nx*(nx*(Da*R11(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R11(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Db*R21(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R21(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Dc*R31(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R31(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])) +

                 ny*(Da*R11(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R12(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Db*R21(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R22(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Dc*R31(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R32(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])) +

                 nz*(Da*R11(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R13(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Db*R21(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R23(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Dc*R31(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R33(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]))) +

             ny*(nx*(Da*R11(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R12(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Db*R21(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R22(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Dc*R31(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R32(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])) +

                 ny*(Da*R12(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R12(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Db*R22(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R22(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Dc*R32(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R32(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])) +

                 nz*(Da*R12(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R13(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Db*R22(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R23(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Dc*R32(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R33(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]))) +

             nz*(nx*(Da*R11(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R13(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Db*R21(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R23(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Dc*R31(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R33(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])) +

                 ny*(Da*R12(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R13(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Db*R22(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R23(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Dc*R32(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R33(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])) +

                 nz*(Da*R13(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R13(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Db*R23(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R23(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3]) +
                     Dc*R33(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])*R33(Orientation[Neighbours3[i][j]][0],Orientation[Neighbours3[i][j]][1],Orientation[Neighbours3[i][j]][2],Orientation[Neighbours3[i][j]][3])));

            P3[i][j] = D_coeff[0]/(0.5*d*((1/Dij3)+(1/Dji3))*1e-9);
        }

        else{
            P3[i][j] = D_coeff[0]/(0.5*d*(1/Dij3)*1e-9);
        }
    }
}

/*******Final determination of neibours and areas***************/

int Neighbours[N][NNN + 12];
double SurfaceAreas[N][NNN + 12];
double Permeability[N][NNN + 12];

double DC12;
double BV12;
double DC3;
double BV3;


double ST = 1.0 + S;
double SR = 0.1;
double SRold = 1;
double IntervalSize = 0.01;
int counter = 0;

while(Outside_interval(SR,1.0-IntervalSize,1.0+IntervalSize)>0){

    DC12 = 0;
    BV12 = 0;
    DC3 = 0;
    BV3 = 0;

    for(int i=0; i<N; i++){
        for(int j=0; j<NNN + 12; j++){
            Neighbours[i][j] = -2;
            SurfaceAreas[i][j] = -2.0;
            Permeability[i][j] = -2.0;
        }
    }


    for(int i=0; i<N; i++){
        for(int j=0; j<NNN; j++){
            Neighbours[i][j] = Neighbours12[i][j];
            SurfaceAreas[i][j] = A12ij[i][j];
            Permeability[i][j] = P12[i][j];
            if(Neighbours12[i][j]<0)
                BV12=BV12 + A12ij[i][j]/N;
            else
                DC12=DC12 + A12ij[i][j]/N;
        }
         if(Size[i]>ST){
            for(int k=0; k<12; k++){
                if(Neighbours3[i][k]<0){
                    BV3=BV3 + A3ij[i][k]/N;
                    Neighbours[i][NNN + k] = Neighbours3[i][k];
                    SurfaceAreas[i][NNN + k] = A3ij[i][k];
                    Permeability[i][NNN + k] = P3[i][k];
                }
                else{
                    if(Size[Neighbours3[i][k]]>ST){
                        DC3=DC3 + A3ij[i][k]/N;
                        Neighbours[i][NNN + k] = Neighbours3[i][k];
                        SurfaceAreas[i][NNN + k] = A3ij[i][k];
                        Permeability[i][NNN + k] = P3[i][k];
                    }
                }
            }
        }
    }

    printf("povrsine - DC12 = %f,  BV12 = %f, DC3 = %f,  BV3 = %f,  Total = %f ST = %f\n", DC12, BV12, DC3, BV3, DC12 + BV12 + DC3 +  BV3, ST);

    SR = DC12 + BV12 + DC3 + BV3;
    ST = 0.5*ST + 0.5*ST*SR;

    counter = counter + 1;
    if(counter == 1000){
        IntervalSize = IntervalSize + 0.01;
        counter = 0;
    }
}

printf("\n");
for(int i=0; i<N; i++){
    for(int j=0; j<NNN + 12; j++){
        printf("%g ", SurfaceAreas[i][j]);
    }
    printf("\n");
}

return 0;
}
