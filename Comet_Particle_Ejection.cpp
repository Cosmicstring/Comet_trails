#include "Comet_Particle_Ejection.h"

    double Z1,ZBB,TBB,PvBB,PvZ1,T1;
    double molmasavode = (18.02*1.661*1e-27);
    double BolzConst = 1.3806488e-23;
    // Just to eject one particle at a time which is according to Whipple formulae
    int ni=0;

    //To calculate the total production, i.e. total ejected mass
    double Zuu = 0;

int main()
{
    //Declaration of I/O files
    FILE* Provera = fopen("Provera.txt","w");
    FILE* CheckInput = fopen("CheckInput.txt","w");
    FILE* Vecheck = fopen("Vecheck.txt","w");
    FILE* BrzinaiPolozaj_MercuryCestice = fopen("beta0_5_small.in","w");
    FILE* Churyumov_Gerasimenko = fopen("67PComet.aei", "r");
    FILE* Comet_Cartesian = fopen("Cartesian_Comet.txt", "w");

    OrbitalElements p, o1,o;
    PhaseState p1,p2,pcomet;
    
    comet Comet_CG;
    
    vector r,v;
    
    int i,k,kt=1,j;
    int MeteorNum = 0;

    char input[200];

    double Te = 365.256;
    double ae = 149.6e9;
    double K= Te/sqrt(ae*ae*ae);
    double N= 5120;
    double u1,u,u0,eps,mi,m;
    double f0,f1,f2,f3,d1,d2,d3;
    double TA,T,TA0,ta,ta0;
    double t, t0;
    double rp,x3,y3,z3,x2,y2,z2,x1,y1,z1,x,y,z;
    double lambda,beta;
    double Q;
    double timeobs = 365.256;
    double h = 0.01;
    double dt1, dt2;
    double Ugao,Ri;

    double Time;

    srand(time(NULL));

    /*From paper: A homogeneous nucleus for comet 67P/Churyumov–Gerasimenko from its gravity field*/
    
    Comet_CG.R = 2.65; 
    Comet_CG.ro = 0.553e12; 
    Comet_CG.Mass = 9.982e9;

    /*From JPL*/
    p.e = 0.64093009801435;
    p.a = 3.462505635805988; /**[AU]*/
    p.i = ConvertToRadians(7.040168070284028);
    p.argperi = ConvertToRadians(12.80065127074631);
    p.longnode = ConvertToRadians(50.13500663540592);

    cartesian(2.959122082855911e-4, p, &pcomet);

    fprintf(Comet_Cartesian, "%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", pcomet.x, pcomet.y, pcomet.z, pcomet.xd, pcomet.yd, pcomet.zd);

    /*Skip first 4 lines*/
    for(int i=1; i<5; i++)   fgets(input, sizeof input, Churyumov_Gerasimenko);
    
    //Set up the file for mercury
    fprintf(BrzinaiPolozaj_MercuryCestice, ")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)\nstyle (Cartesian, Asteroidal, Cometary) = Ast\n)---------------------------------------------------------------------\n");
    
while(!feof(Churyumov_Gerasimenko)){ 
    
    fgets(input, sizeof input, Churyumov_Gerasimenko);
    if(DEBUG) printf("Ok1\n");

    sscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &Time, &r.x, &r.y, &r.z, &v.x, &v.y, &v.z, &Comet_CG.a, &Comet_CG.e, &Comet_CG.i, &Comet_CG.argperi, &Comet_CG.longnode, &Comet_CG.meananom);
    if(DEBUG) printf("Ok2\n");

    fprintf(CheckInput, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Time, r.x, r.y, r.z, v.x, v.y, v.z, Comet_CG.a, Comet_CG.e, Comet_CG.i, Comet_CG.argperi, Comet_CG.longnode, Comet_CG.meananom);

    /*Vracam u km i km/d posto se sa tim radi u GiveMeParticles*/
    double AUkm = AU * 1e-3;
    r.x = r.x * AUkm;
    r.y = r.y * AUkm;
    r.z = r.z * AUkm;

    v.x = v.x * AUkm;
    v.y = v.y * AUkm;
    v.z = v.z * AUkm;
    
    //Time in days so it can be added to Epoch of the comet
    GiveMeParticles(BrzinaiPolozaj_MercuryCestice, Vecheck, Provera, MeteorNum, Zuu, ni, Comet_CG.Mass, Comet_CG.R, Time*365.25, p, r, v); 

}
    
    if(DEBUG) printf("%.10lf\n", Comet_CG.Mass);

    fclose(Provera);
    fclose(Vecheck);
    fclose(CheckInput);
    fclose(Churyumov_Gerasimenko);
    fclose(BrzinaiPolozaj_MercuryCestice);
    fclose(Comet_Cartesian);

    return 0;
}