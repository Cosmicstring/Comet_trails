#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

/**

THIS IS A CODE WRITTEN FOR MAKING A BINNED IMAGE FROM THE XYZ DATA OF DUST PARTICLES
EJECTED FROM THE COMET SURFACE. THE IMAGE IS LATER ON FORWARDED TO MATLAB CODES WHICH
PROCESS THEM FURTHER. 

**/


#define MaxN 2001
#define DEBUG 0


long A[MaxN*MaxN];
long Test[MaxN];

FILE* SaSOHO;
FILE* Testing;
FILE* TestingOutput;
FILE* DebugText;

int main(){

    SaSOHO = fopen("SaSOHO.txt", "r");
    Testing = fopen("Test.txt", "r");
    TestingOutput = fopen("TestingOutput.txt", "w");
    DebugText = fopen("Debug.txt", "w");
    double Testx, Testy, x,y, MinX, MinY, MaxX, MaxY, binstep;
    long i,j;
    char input[256];

    MaxX = 0.2; MaxY = 0.2;
    MinX = -0.2; MinY = -0.2; /** Vrednost celog vidnog polja sa Suncem u Centru */
    memset(Test, sizeof(Test), 0);
    memset(A, sizeof(A), 0);
    binstep = (MaxX-MinX)/2000;

    if(DEBUG) printf("%lf\n", binstep);

    /*
    TESTING ROUTINE 

    for(i=0; i<10; i++){
        /*fgets(input, sizeof input, Testing);
        sscanf(input, "%lf %lf", &Testx, &Testy);
        scanf("%lf%lf", &Testx, &Testy);
        printf("%lf\t%ld\t%lf\t%ld\n", fabs(Testx)/binstep, lround(fabs(Testx)/binstep), fabs(Testy)/binstep, lround(fabs(Testy)/binstep));
        printf("1: %ld\n", lround(fabs(Testx)/binstep) + lround(fabs(Testy)/binstep)*10);
        Test[lround(fabs(Testx)/binstep) + lround(fabs(Testy)/binstep)*10]++;
    }

    for(i=0; i<10; i++)
        for(j=0; j<10; j++){
            if(j==9) printf("%ld\n", Test[i+j*10]);
            else printf("%ld ", Test[i+j*10]);
        }*/
    i=0;
    while(!feof(SaSOHO)){
        fgets(input, sizeof input, SaSOHO);
        sscanf(input, "%lf %lf", &x, &y);
        if(DEBUG) printf("Ok: %ld %ld %ld\n", A[lround(fabs(x)/binstep) + lround((fabs(y)/binstep))*2000], lround(fabs(x)/binstep) + lround((fabs(y)/binstep))*2000, i);
        A[lround(fabs(x)/binstep) + lround((fabs(y)/binstep))*2000]++;
        if(DEBUG) printf("Ok: %ld %ld %ld\n", A[lround(fabs(x)/binstep) + lround((fabs(y)/binstep))*2000], lround(fabs(x)/binstep) + lround((fabs(y)/binstep))*2000, i);
        i++;
    }

    for(i=0; i<2000*2000; i++){
        if((i%2000==0)) fprintf(TestingOutput, "%ld\n", A[i]);
        else fprintf(TestingOutput, "%ld\t", A[i]);
    }

    if(DEBUG) printf("OkIzaso %ld\n", A[1206304]);
/*
    for(i=0; i<2000*2000; i++){
        if(i%2000==0) fprintf(TestingOutput, "%ld\n", A[i]);
        else fprintf(TestingOutput, "%ld\t", A[i]);
    }*/

    if(DEBUG) printf("Nema SegFault\n");

    fclose(SaSOHO);
    fclose(Testing);
    fclose(DebugText);
    fclose(TestingOutput);
    return 0;
}
