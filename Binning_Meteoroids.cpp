#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MaxN 2001
#define resolution 1024
#define DEBUG 1


long A[MaxN*MaxN];
long P1[MaxN*MaxN], P2[MaxN*MaxN];
long Test[MaxN];

FILE* SaSOHO;
FILE* Presek1;
FILE* Presek2;
FILE* TestingOutput;
FILE* DebugText;
FILE* Presek1Ispis;
FILE* Presek2Ispis;

int main(){

    Presek1 = fopen("Trenutak1.txt", "r");
    Presek2 = fopen("Trenutak2.txt", "r");
    Presek1Ispis = fopen("Presek1Ispis.txt", "w");
    Presek2Ispis = fopen("Presek2Ispis.txt", "w");
    SaSOHO = fopen("SaSOHO.txt", "r");
    TestingOutput = fopen("TestingOutput_Beta06.txt", "w"); /**Od najvece dimenzije do najmanje 3-0 {0.05, 0.5, 0.6, 1}*/
    DebugText = fopen("Debug.txt", "w");
    double foo, xi, yi, Testx, Testy, x,y, MinX, MinY, MaxX, MaxY, binstep, binstepP1;
    long i,j;
    char input[256];

    MaxX = 17; MaxY = 17;
    MinX = -17; MinY = -17; /** [stepeni] Vrednost celog vidnog polja sa Suncem u Centru, Sekanin rad za velicinu piksela(Uzeo sam kameru C3)*/
    memset(Test, sizeof(Test), 0);
    memset(A, sizeof(A), 0);
    binstep = (MaxX-MinX)/resolution; /**Kamera C3 ima 1024*1024 piksela*/

    if(DEBUG) printf("%lf\n", binstep);

    /*for(i=0; i<10; i++){
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

    while(!feof(SaSOHO)){
        fgets(input, sizeof input, SaSOHO);
        sscanf(input, "%lf %lf", &x, &y);
       // if(DEBUG) printf("Ok: %ld %ld %ld\n", A[lround(fabs(x)/binstep) + lround((fabs(y)/binstep))*2000], lround(fabs(x)/binstep) + lround((fabs(y)/binstep))*2000, i);
        if(x<0) xi = MaxX - fabs(x);
        else
            if(x>0) xi = MaxX+x;
            else xi = x;
        if(y<0) yi = MaxY+fabs(y);
        else
            if(y>0) yi = MaxY - y;
            else yi = y;
        A[lround(xi/binstep) + lround(yi/binstep)*resolution]++;
       // if(DEBUG) printf("Ok: %ld %ld %ld\n", A[lround(fabs(x)/binstep) + lround((fabs(y)/binstep))*2000], lround(fabs(x)/binstep) + lround((fabs(y)/binstep))*2000, i);
    }


    for(i=0; i<resolution*resolution; i++){
        if(i%resolution==0) fprintf(TestingOutput, "%ld\n", A[i]);
        else fprintf(TestingOutput, "%ld\t", A[i]);
    }

    /** Presek 1 */

    memset(input, sizeof input, 0);
    memset(P1, sizeof P1, 0);
    double MaxXP1, MinXP1, MaxYP1, MinYP1;
    double binstepP1X, binstepP1Y;

    MaxXP1 = -100; MaxYP1 = -100;
    MinXP1 = 100; MinYP1 = 100;

    while(!feof(Presek1)){
        fgets(input, sizeof input, Presek1);
        sscanf(input, "%lf %lf %lf", &foo, &x, &y);
        if(MaxXP1<fabs(x)) MaxXP1 = fabs(x);
        if(MaxYP1<fabs(y)) MaxYP1 = fabs(y);
        if(MinXP1>fabs(x)) MinXP1 = fabs(x);
        if(MinYP1>fabs(y)) MinYP1 = fabs(y);
    }

    binstepP1X = (fabs(MaxXP1)-fabs(MinXP1))/resolution;
    binstepP1Y = (fabs(MaxYP1)-fabs(MinYP1))/resolution;
    if(DEBUG) printf("%.10lf %lf %lf %lf %lf %lf %lf %lf\n", binstepP1X, fabs(MaxXP1)-fabs(MinXP1), (fabs(MaxYP1) - fabs(MinYP1))/binstepP1X, (fabs(MaxXP1)-fabs(MinXP1))/binstepP1Y, MaxXP1, MaxYP1, MinXP1, MinYP1);
    rewind(Presek1);


    i=0;
    while(!feof(Presek1)){
        fgets(input, sizeof input, Presek1);
        sscanf(input, "%lf %lf %lf", &foo, &x, &y);
        P1[lround((fabs(x)-fabs(MinXP1))/binstepP1X) + lround((fabs(y)-fabs(MinYP1))/binstepP1Y)*resolution]++;
    }

    for(i=0; i<resolution*resolution; i++)
        if(i%resolution==0) fprintf(Presek1Ispis, "%ld\n", P1[i]);
        else fprintf(Presek1Ispis, "%ld\t", P1[i]);

    /** Presek 2*/

    memset(input, sizeof input, 0);
    memset(P2, sizeof P2, 0);
    double MaxXP2, MinXP2, MaxYP2, MinYP2;
    double binstepP2X, binstepP2Y;

    MaxXP2 = -100; MaxYP2 = -100;
    MinXP2 = 100; MinYP2 = 100;

    while(!feof(Presek2)){
        fgets(input, sizeof input, Presek2);
        sscanf(input, "%lf %lf %lf", &foo, &x, &y);
        if(MaxXP2<fabs(x)) MaxXP2 = fabs(x);
        if(MaxYP2<fabs(y)) MaxYP2 = fabs(y);
        if(MinXP2>fabs(x)) MinXP2 = fabs(x);
        if(MinYP2>fabs(y)) MinYP2 = fabs(y);
    }

    if(DEBUG) printf("%lf %lf %lf %lf\n", MaxXP2, MinXP2, MaxYP2, MinYP2);

    binstepP2X = (fabs(MaxXP2)-fabs(MinXP2))/resolution;
    binstepP2Y = (fabs(MaxYP2)-fabs(MinYP2))/resolution;
    if(DEBUG) printf("%.10lf %lf %lf %lf %lf %lf %lf %lf\n", binstepP2X, fabs(MaxXP2)-fabs(MinXP2), (fabs(MaxYP2) - fabs(MinYP2))/binstepP2X, (fabs(MaxXP2)-fabs(MinXP2))/binstepP2Y, MaxXP2, MaxYP2, MinXP2, MinYP2);

    rewind(Presek2);

    while(!feof(Presek2)){
        fgets(input, sizeof input, Presek2);
        sscanf(input, "%lf %lf %lf", &foo, &x, &y);
        P2[lround((fabs(x)-fabs(MinXP2))/binstepP2X) + lround((fabs(y)-fabs(MinYP2))/binstepP2Y)*resolution]++;
    }

    for(i=0; i<resolution*resolution; i++)
        if(i%resolution==0) fprintf(Presek2Ispis, "%ld\n", P2[i]);
        else fprintf(Presek2Ispis, "%ld\t", P2[i]);

    fclose(SaSOHO);
    fclose(Presek1Ispis);
    fclose(Presek1);
    fclose(Presek2);
    fclose(DebugText);
    fclose(TestingOutput);
    return 0;
}
