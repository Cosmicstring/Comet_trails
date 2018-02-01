#include "Particle_Orbit.h"


double mi,m2;
double Es;
double ro; // [kg/km^3]
double R; //km

OrbitalElements o,p,o0, O_preseka;
FILE* TailFiles[61];

    double Time;
    double Tp = 2456625.27616;
    double n=0;
    double MaxEffect=-1;

int main(){

	vectorr r, r0, sam, v2, v,v0;
	PhaseState p1,p2,p, P_preseka;
	meteoroid Meteor;

    int i;
    
    char Prefix[20];
    
    memset(Prefix, sizeof Prefix, 0);

    for(i=0; i<61; i++)
    {
        sprintf(Prefix, "TMoment_%d.txt", i);
        TailFiles[i] = fopen(Prefix, "w");
    }

    // Setup for the vid files 
    int NFrames = 100; 
    FILE* VideoFile[NFrames];
    
    memset(Prefix, sizeof Prefix, 0);
            
    for(int NFrame=0; NFrame<NFrames; NFrame++)
    {
        sprintf(Prefix, "VideoFile_%d.txt", NFrame);
        VideoFile[NFrame] = fopen(Prefix, "w");
    }
            

    //Input files
    FILE* BrzinaiPolozaj = fopen("BrzinaiPolozaj.txt", "r");

    //Ouput Files
    FILE* Trenutak1 = fopen("Trenutak1.txt", "w");
	FILE* Trenutak2 = fopen("Trenutak2.txt", "w");
	
    FILE* Checkinput = fopen("Checkinput.txt", "w");
    FILE* TrenutakZ = fopen("TrenutakZ.txt","w");
    FILE* Checkorbit = fopen("Checkorbit.txt","w");
    FILE* Odnos = fopen("Odnos.txt", "w");
    FILE* Checking = fopen("CheckIntersection.txt", "w");
    
    FILE* Nodes1 = fopen("Nodes.txt", "w");
    FILE* Nodesa = fopen("NodesAscending.txt", "w");
    FILE* Nodesd = fopen("NodesDescending.txt", "w");

    FILE* TMoment_slika_prim = fopen("TMoment_slika_prim.txt", "w");
    FILE* TMoment_slika_1= fopen("TMoment_slika_1.txt", "w");
    FILE* TMoment_slika_2 = fopen("TMoment_slika_2.txt", "w");

    char input[512];

	vectorr k1, k2, k3, k4;
	vectorr l1, l2, l3, l4;
	vectorr temp;
	vectorr Vs,vc,Mi;

	double h;
	double E, Ek, Ep;
	double rp, V, betac, ratio;
	double Te = 365.256;
	double ae = 149.6e6;
    double K= Te/sqrt(ae*ae*ae);
    double N=5120;
    double ra, rd0, rd;
    double S;
    double dS=0;
    double n0=0.01;
    long hi=0;
    long he=0;
    long hp=0;
    int ti=0;
    int kt=1;
    int u1, bla =0, foo = 0;

    o.a = AU; /**[m]*/
    o.e = 0.00001;
    o.i = 0;
    o.argperi = 0;
    o.longnode = 0;

while(1){

    int FrameCount = 0;

    int nf, BI=0;
    int EjectionColor;

    double timep = 0;
    double T0;
    double Tslike_prim = Tp - 4.2*sat;
    double Tslike_1 = Tp + 16.4*sat;
    double Tslike_2 = Tp + 25.8*sat;
    double TVideo0 = Tp - 1;
    double TVideo1 = Tp + 3;

    double IntVe, TA0;
    double T00 = Tp - 30;
    double dt;
    double eps = 1e-3;
    double rocgs,Rcgs;
    
    bool PresekaoEkliptiku = false;
    
    vectorr TestUbrzanje, ri,rmax,rmin;
    
    fgets(input, sizeof input, BrzinaiPolozaj);
    sscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", &T0, &Meteor.ro, &Meteor.R, &r.x, &r.y, &r.z, &v.x, &v.y, &v.z, &v2.x, &v2.y, &v2.z, &IntVe, &TA0, &EjectionColor);
  //  fprintf(Checkinput, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n", T0, Meteor.ro, Meteor.R, r.x, r.y, r.z, v.x, v.y, v.z, v2.x, v2.y, v2.z, IntVe, TA0, EjectionColor);


   	ro = Meteor.ro;
	R = Meteor.R;
	rp = intenzitet(r);
	Time = T0;
	dt = 0;
	nf = -1;

    while(T0>(T00+dt))
    {
        dt++;
    }

    n++;

    p.x = r.x;
    p.y = r.y;
    p.z = r.z;

    p.xd = v.x;
    p.yd = v.y;
    p.zd = v.z;

    keplerian(gama*msol, p, &o0);

   /* r.x = -9070435.45614153;
    r.y = 11538032.6357978;
    r.z = -6083952.81855579;

    v.x = 5238017.67614801;
    v.y = -9729772.53481191;
    v.z = 898543.194909737;*/

    Ep = -((gama * msol*m2)/rp);
    Ek = (m2*V*V)/2.0;
    E = Ek + Ep;

    betac = 3.0*L/(16.0*pi*ro*R);

	m2 = ro * (4.0/3.0) * pi * R*R*R;
	mi = gama * msol - betac;
    ratio = betac/(gama*msol);
 
    h = Adapth(h,K,N,r);
   
	while(Time<=Tp+30)
    {

        if(Video) h = 0.5*sat;
        else
        {
            if((Time>=Tp-1) && (Time<=Tp+1)) h = 2*minut;
            else h = Adapth(h,K,N,r);
        }

        r0 = r;

		l1.x = v.x;
		l1.y = v.y;
		l1.z = v.z;

        k1 = ubrzanja(r,v,R,m2);

		l2.x = l1.x + k1.x * h/2.0;
		l2.y = l1.y + k1.y * h/2.0;
		l2.z = l1.z + k1.z * h/2.0;

		temp.x = r.x + l1.x * h/2.0;
		temp.y = r.y + l1.y * h/2.0;
		temp.z = r.z + l1.z * h/2.0;

		k2 = ubrzanja(temp,l2,R,m2);

		l3.x = l1.x + k2.x * h/2.0;
		l3.y = l1.y + k2.y * h/2.0;
		l3.z = l1.z + k2.z * h/2.0;

		temp.x = r.x + l2.x * (h/2.0);
		temp.y = r.y + l2.y * (h/2.0);
		temp.z = r.z + l2.z * (h/2.0);

		k3 = ubrzanja(temp, l3,R,m2);

		l4.x = l1.x + k3.x * h;
		l4.y = l1.y + k3.y * h;
		l4.z = l1.z + k3.z * h;

		temp.x = r.x + l3.x * h;
		temp.y = r.y + l3.y * h;
		temp.z = r.z + l3.z * h;

        k4 = ubrzanja(temp, l4,R,m2);

		v.x = v.x + (h/6.0)*(k1.x + 2*k2.x + 2*k3.x + k4.x);
		v.y = v.y + (h/6.0)*(k1.y + 2*k2.y + 2*k3.y + k4.y);
		v.z = v.z + (h/6.0)*(k1.z + 2*k2.z + 2*k3.z + k4.z);

		r.x = r.x + (h/6.0)*(l1.x + 2*l2.x + 2*l3.x + l4.x);
		r.y = r.y + (h/6.0)*(l1.y + 2*l2.y + 2*l3.y + l4.y);
		r.z = r.z + (h/6.0)*(l1.z + 2*l2.z + 2*l3.z + l4.z);

        // Write out to a file for the date and time of interest
        if((Time<=Tslike_prim+10*minut) && (Time>=Tslike_prim-10*minut)) {
            fprintf(TMoment_slika_prim, "%.10lf\t%.10lf\t%lf\t%.10lf\t%.10lf\t%.10lf\n", T0 - T00, Time-T0, n, r.x, r.y, r.z);
        }

        if((Time<=Tslike_1+10*minut) && (Time>=Tslike_1-10*minut)) {
            fprintf(TMoment_slika_1, "%.10lf\t%.10lf\t%lf\t%.10lf\t%.10lf\t%.10lf\n", T0 - T00, Time-T0, n, r.x, r.y, r.z);
        }

        if((Time<=Tslike_2+10*minut) && (Time>=Tslike_2-10*minut)) {
            fprintf(TMoment_slika_2, "%.10lf\t%.10lf\t%lf\t%.10lf\t%.10lf\t%.10lf\n", T0 - T00, Time-T0, n, r.x, r.y, r.z);
        }

        // Write out for the frames for the video
        if(Video)
        {
            if((Time>=TVideo0) && (Time<=TVideo1))
            {
                fprintf(VideoFile[FrameCount], "%.2lf\t%.10lf\t%d\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%d\n", T0 - T00, Time-T0, nf, n, r.x, r.y, r.z, intenzitet(r)- 10.0*Rs, EjectionColor);
                FrameCount++;
            }
        }

      //  if(DEBUG) printf("%lf\t%lf\n", Time, intenzitet(r));
        if( ((Time>=(T00+dt)) && (T0<=(T00+dt))) && (intenzitet(r)> Rs))
        {
            nf = dt;
            fprintf(TailFiles[nf], "%.10lf\t%.10lf\t%d\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%d\n", T0 - T00, Time-T0, nf, n, r.x, r.y, r.z, intenzitet(r)- 10.0*Rs, EjectionColor);
            dt=dt+1;
        }
        Time = Time+h;

		p1.x = r.x;
        p1.y = r.y;
        p1.z = r.z;

        p1.xd = v.x;
        p1.yd = v.y;
        p1.zd = v.z;

		keplerian(gama*msol,p1,&o);
        
        ra = o.a*(1-o.e*o.e)/(1+o.e*cos(o.argperi));
        rd = o.a*(1-o.e*o.e)/(1-o.e*cos(o.argperi));
        
        if( ((r0.z<0) && (r.z>0)) || ((r0.z>0) && (r.z<0)) )
        {
            PresekaoEkliptiku = true;
            if(DEBUG){
                //printf("%.10lf\t%.10lf\n", r0.z, r.z);
                foo++;
            }
            if(fabs(r.z) <= fabs(r0.z))
            {
                rmax = r0;
                rmin = r;
            }
            else
            {
                rmax = r;
                rmin = r0;
            }
            ri.x = rmin.x*(fabs(rmax.z)/(fabs(r0.z) + fabs(r.z))) + rmax.x*(fabs(rmin.z)/(fabs(r0.z) + fabs(r.z)));
            ri.y = rmin.y*(fabs(rmax.z)/(fabs(r0.z) + fabs(r.z))) + rmax.y*(fabs(rmin.z)/(fabs(r0.z) + fabs(r.z)));
            ri.z = 0;

            if( (intenzitet(ri)<= AU + 0.0005*AU) && (intenzitet(ri)>= AU - 0.0005*AU) )
            {
                fprintf(TrenutakZ, "%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",n, Time, rp, ra, rd);
            }
            keplerian(gama*msol, p1, &O_preseka);

            if( (r0.z>0) && (r.z<0) ){
                rd = O_preseka.a*(1-O_preseka.e*O_preseka.e)/(1-O_preseka.e*cos(O_preseka.argperi));
                fprintf(Checking, "%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", ri.x, ri.y, ri.z, intenzitet(ri)/AU, ra/AU, rd/AU);
                fprintf(Nodesd, "%d\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", PresekaoEkliptiku, n, O_preseka.a, O_preseka.e, O_preseka.argperi, TA0, IntVe, rd/AU, rd/AU-0.7854238483);
                fprintf(Trenutak1, "%.10lf\t%.10lf\t%.10lf\t%d\n", Time - T0, ri.x/AU, ri.y/AU, EjectionColor);
            }
            else{
                if((r0.z<0) && (r.z>0)){
                    ra = O_preseka.a*(1-O_preseka.e*O_preseka.e)/(1+O_preseka.e*cos(O_preseka.argperi));
                    fprintf(Nodesa, "%d\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", PresekaoEkliptiku, n, O_preseka.a, O_preseka.e, O_preseka.argperi, TA0, IntVe, ra/AU, ra/AU-0.0126468222);
                    fprintf(Trenutak2, "%.10lf\t%.10lf\t%.10lf\t%d\n", Time - T0, ri.x/AU, ri.y/AU, EjectionColor);
                }
            }
            u1 = 0;
        }

		V = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);

		rp = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);

    }
fprintf(Nodes1, "%d\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", PresekaoEkliptiku, n, o.a, o.e, o.argperi, TA0, IntVe, ra/AU, rd/AU, ra/AU-0.01265, rd/AU-0.7854238483);

if(n/60010 > n0){
    printf("%.2lf\n", n/60010);
    printf("%d\n", foo);

    n0+=0.01;
}
    if(o0.e>1)
    {
        hi +=1;
    }
    if(o0.e<1)
    {
        he +=1;
    }
    if(o0.e==1)
    {
        hp +=1;
    }
    Time = 0;
    PresekaoEkliptiku = false;
        if(feof(BrzinaiPolozaj)){
            printf("Upao u FEOF\n");
            printf("%.10lf\n", MaxEffect);
            break;
        }
}
 	fclose(Nodes1);
	fclose(Checkorbit);
	fclose(Checkinput);
	fclose(Nodesa);
	fclose(Nodesd);
	fclose(BrzinaiPolozaj);
	fclose(Trenutak1);
	fclose(Trenutak2);
	fclose(Checking);
	fclose(TMoment_slika_2);
    fclose(TMoment_slika_1);
    fclose(TMoment_slika_prim);
	fclose(TrenutakZ);

 return 0;
}






