#include "Generate_SOHO_Image.h"


FILE* Earthcoordinates = fopen("EARTHMOO.aei","r");
FILE* SohoAndEarthcoordinates = fopen("SohoAndEarthcoordinates.txt", "w");
FILE* KoordinateSvaTriTela = fopen("Orijentacija_i_koordinate.txt", "w");
FILE* CheckIfNumbersWork = fopen("CheckIfNumbersWork.txt", "w");
FILE* CheckForMisinPredlog = fopen("CheckForMisinPredlog.txt", "w");
FILE* TMoment_slika = fopen("TMoment_slika.txt", "r");         

//Prefix for the name of the TMoment Files
char Prefix[20];
char input[256];
int main()
{


    //Pics and Vids
    FILE* VideoFile;
    FILE* TMoment;

    vectorr Sunce;
    Sunce.x =0; Sunce.y=0; Sunce.z=0;

    OrbitalElements Orbitalfoo,o,o1;
    PhaseState foo, p;
    double M,M0;
    double x,y,z;
    double t, tobserved; /** u zavisnosti, koji moment gledam dobija vrednost */
    double mi;
    int n1=0;

/** Podaci za Zemlju http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
    A pozicija zemlje 26. og novembra 2013 sa http://omniweb.gsfc.nasa.gov/coho/helios/planet.html
*/

    o.a = 1.00000011*AU;
    o.e = 0.01671022;
    o.i = ConvertToRadians(0.00005);
    o.longnode = ConvertToRadians(348.73936);
    o.longperi = ConvertToRadians(102.94719);
    o.argperi = 2*pi + (o.longperi - o.longnode);
    o.meanlongitude = ConvertToRadians(100.46435);

    cartesian(GMs, o, &foo);
    keplerian(GMs, foo, &Orbitalfoo);

    double MeanAnom_At_28_November_2013 = ConvertToRadians(323.194966667);// degrees
    //Mean anomaly for correction
    M0 = MeanAnom_At_28_November_2013;

    tobserved = T_ISON_perihel-30;
    //Mean anomaly for the seeked position
    M = M0 + n*(tobserved - T_ISON_perihel);

    if(DEBUGzaKonverzije) printf("MeanAnom: %lf\t%lf\t%lf\n", ConvertToDegrees(M0), ConvertToDegrees(M), n*(tobserved - T_Rosetta_perihel)*180/pi);//(ConvertToDegrees(M) - 360*((int)ConvertToDegrees(M)/360)));
    else
    {
      /*

        SORTIRATI PO VREMENU SVE CESTICE I ONDA UCITAVATI ZA SVAKI TRENUTAK KOORDINATE ZEMLJE I KOORDINATE CESTICA, PROJEKTOVATI 
        I FORMIRATI REP UZ POMOC MATLABA,
        VIDETI DA LI LEPO RADI MATLAB SA VISE ULAZNIH FAJLOVA (TREBALO BI DA SAM TO ODRADIO RANIJE)

      */
      
      // Vectors needed for projection
      vectorr particlenorm, soho /** Sunce-SOHO*/, particle /** SOHO - Cestica*/, CesticaHC/** u odnosu na Sunce*/, Xprojekcije, Yprojekcije, Z0/** Jedinicni izabrani vektor za pozitivni smer Z ose*/;

      //Number of Files
      int NFiles;
      int nfile = 0;
      int color;
      int whichone;

      printf("What do you want to project?\n For Video: 1\t For Pics: 2\n");
      scanf("%d", &whichone);

      switch(whichone)
      {
        case 1: 
          NFiles = 100;
          break;
        case 2:
          NFiles = 61;
          break;
      }

      printf("NFiles: %d\n", NFiles);

      double foo, cosTheta, fi, fis, theta;
      double r,rs;
      double X, Y, Particlexy, Sohoxy,SohoParticlexy, RParticleSoho;
      double cosAzimut, cosDeklinacija,Azimut, Deklinacija;

      while(nfile<NFiles)
      {
        
        int nline = 0;
        
        switch(whichone)
        {
          case 1:
          {
            memset(Prefix, 0, sizeof Prefix);
            sprintf(Prefix, "VideoFile_%d.txt", nfile);
            VideoFile = fopen(Prefix, "r");
            break;
          }
          case 2:
          {
            memset(Prefix, 0, sizeof Prefix);
            sprintf(Prefix, "TMoment_%d.txt", nfile);
            printf("%s\n", Prefix);
            TMoment = fopen(Prefix, "r");
            break;
          }
        }
        
        memset(Prefix, 0, sizeof Prefix);
        sprintf(Prefix, "SaSOHO_%d.txt", nfile);
        FILE* SaSOHO = fopen(Prefix, "w");
 
        while(1)
        {

            //fscanf(TMoment, "%lf %lf %lf %lf %lf %lf %lf %lf", &foo, &foo, &foo, &foo, &x, &y, &z, &foo);
            if(whichone==1) fgets(input, sizeof input, VideoFile);
            if(whichone==2) fgets(input, sizeof input, TMoment);
            nline ++;
          
            if(whichone==1) 
                if(feof(VideoFile)) break;
            if(whichone==2) 
                if(feof(TMoment)) break;
            
            if(whichone==1)
              sscanf(input, "%lf %lf %lf %lf %lf", &foo, &foo, &x, &y, &z);
            if(whichone==2)
              sscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf %d", &foo, &foo, &foo, &foo, &x, &y, &z, &foo, &color);

            CesticaHC.x = x; CesticaHC.y = y; CesticaHC.z = z;

            n1++;

            r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);

            rs = r - 1.5e6; /*km*/

            soho.x = p.x*rs/r;
            soho.y = p.y*rs/r;
            soho.z = p.z*rs/r;

            if(n1==1) fprintf(SohoAndEarthcoordinates, "%.10lf\t%.10lf\t%.10lf\n%.10lf\t%.10lf\t%.10lf\n", soho.x, soho.y, soho.z, p.x, p.y, p.z);

            particle.x = x-soho.x;
            particle.y = y-soho.y;
            particle.z = z-soho.z;

            Z0.x = 0; Z0.y = 0; Z0.z = 1; /** Vrednost za Z koordinatu rucno namestiti?*/

            Xprojekcije = vp(Z0, soho);
            Yprojekcije = vp(soho, Xprojekcije);

            X = intenzitet(Xprojekcije);
            Xprojekcije.x /= X; Xprojekcije.y /= X; Xprojekcije.z /= X;

            Y = intenzitet(Yprojekcije);
            Yprojekcije.x /= Y; Yprojekcije.y /= Y; Yprojekcije.z /= Y;

            double d;
            vectorr SohoTackaURavni;
           
            SohoTackaURavni = PresekRavniIPrave(Sunce, Norm(particle), soho, soho); /** Vektor Sunce - Tacka u Ravni polja*/
           
            double ProjekcijaX, ProjekcijaY;

            ProjekcijaX = sp(SohoTackaURavni, Xprojekcije);
            ProjekcijaY = sp(SohoTackaURavni, Yprojekcije);

            if(DEBUGzaMisinPredlog) fprintf(CheckForMisinPredlog, "%.10lf\t%.10lf\n", ProjekcijaX, ProjekcijaY);

            Azimut = atan(ProjekcijaX/intenzitet(soho));
            Deklinacija = atan(ProjekcijaY/intenzitet(soho));
           // if(DEBUGzaMisinPredlog) fprintf(SaSOHO, "%lf\t%lf\t%lf\t%lf\n", d, sp(vp(SohoTackaURavni, Xprojekcije), Yprojekcije), ProjekcijaX, ProjekcijaY);

            fprintf(SaSOHO, "%.10lf\t%.10lf\n", ConvertToDegrees(Azimut), ConvertToDegrees(Deklinacija));//, color);
        }
        switch(whichone)
        {
          case 1:
            fclose(VideoFile);
            break;
          case 2:
            fclose(TMoment);
            break;
        }
        fclose(SaSOHO);
        nfile++;
        printf("%d\n", nfile);
      }
  }

  fclose(CheckForMisinPredlog);
  fclose(SohoAndEarthcoordinates);
  fclose(KoordinateSvaTriTela);
  return 0;
}
