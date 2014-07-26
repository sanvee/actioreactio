#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_const_mksa.h>
#include<math.h>

int main (int argc, char **argv)
{
  if (argc !=4)
    {
      printf("Please Enter  [number of samples] [delta] [flag] as argument");
      printf("Flag = 0 => trajektorien berechnen");
      printf("Flag = 1 => energie berechnen");
      printf("Flag = 2 => energie iteration");
      return 0;
    }
  unsigned long  nbr_of_smpls= strtoull((argv[1]),NULL,10);
  double delta = atof(argv[2]); 
  int fl = atoi(argv[3]);

  double masse_a = 1.0;
  double R_a = 1.0;
  double masse_b = 1.0;
  double R_b = 1.0;
  double sparam;
  double delta_t = 0.0005;
  double k=0.25;
  double epsilon = (1+delta)/(1-delta);
  double M = (masse_a+epsilon*masse_b);
  double red_masse = (double)(epsilon*masse_a*masse_b)/(double)(masse_a+epsilon*masse_b);
  double therm_en1=0.01;
  double therm_en2=0.01;

  double sigma_a,sigma_b;
  double e_kin,e_pot,e_tot;
  int i;
  FILE *fpa ,*fpb,*fpc, *fpd, *fpe, *fpf, *fpg;
  char filename[16];
  int double_sz = sizeof(double);

  double va[2]={0,0}, vb[2]={0,0}, vab[2]={0,0};
  double phi=0,theta=0;
  double vb_norm=0, vr_norm=0, ar_norm2=0,r_norm=0;
  double vr[2]={0,0}, vr2[2]={0,0},ar[2]={0,0},ar2[2]={0,0},v1_l[2]={0,0},v2_l[2]={0,0},v_l[2]={0,0};  

  double omega=0;
  double V[2]={0,0};
  double R[2] = {0,0};
  double r[2]={0,0};
  double r_1[2]={0,0};
  double r_2[2]={0,0};
  double r1_l[2]={0,0};
  double r2_l[2]={0,0};
  double r_l[2]={0,0};

  double delta_E1=0,delta_E2=0,delta_Etot=0;
  double sum_delta_E1=0, sum_delta_E2=0, sum_delta_Etot=0;
  double mean_delta_E1=0, mean_delta_E2=0,mean_delta_Etot=0;

  double r_tot = R_a+R_b;
  gsl_rng *rand;  
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  
//reading from /dev/urandom to inntialise seed---------------------------------
  fpa = fopen("/dev/urandom","r");
  fread (&gsl_rng_default_seed,sizeof(gsl_rng_default_seed),1,fpa);
  fclose(fpa);
  T = gsl_rng_mt19937; // default generator
  rand = gsl_rng_alloc(T);

   fpb = fopen("data/stoß.bin","wb");
   fpc = fopen("data/labor.bin","wb");
   fpg = fopen("data/finaltem.bin","wb");

   sprintf(filename,"data/statenerg%i.bin",(int)(delta*101));
   fpe = fopen(filename,"wb");
   int counter =0;
   int max =30;
   for ( counter =0;counter <= max; counter++)
     {
      if (fl ==1 || counter ==0 )
	    {
	       fpd = fopen("data/energie.bin","wb");
	    }
      if ( counter == max)
	    {
	      fpf = fopen("data/energie_stat.bin","wb");
	    }
       printf("%i\n",counter);
      
       therm_en1+=mean_delta_E1;
       therm_en2+=mean_delta_E2;

       mean_delta_E1=0;
       mean_delta_E2=0;
       mean_delta_Etot=0;

       sigma_a = sqrt((fabs(therm_en1))/masse_a);
       sigma_b = sqrt((fabs(therm_en2))/masse_b);
       
       for (i=0; i<=nbr_of_smpls-1;i++)
	 {
	   va[0]= sigma_a*sqrt(-2*log(gsl_rng_uniform(rand)));
	   va[1]=0;

	   vb_norm =sigma_b*sqrt(-2*log(gsl_rng_uniform(rand)));
	   phi = 2*M_PI*gsl_rng_uniform(rand);//Winkel der Geschwindigkeit des 2ten Teilchen
	  
	   vb[0] = vb_norm*cos(phi);// Vb_x
	   vb[1] = vb_norm*sin(phi);// Vb_Y
	   
	   vab[0]=(vb[0]-va[0]);
	   vab[1]=(vb[1]-va[1]);
       
	   vr[0]=sqrt(vab[0]*vab[0]+vab[1]*vab[1]); // Geschindigkeit parralel zu Vb-Va im Schwerpunktsystem 
	   vr[1] = 0;

	   theta= atan2(vab[1],vab[0]);// Winkelzwichen Va und Vr
	   sparam =r_tot- 2*r_tot*gsl_rng_uniform_pos(rand); // <---- Stoßparameter wird hier gezogen.
	   omega = asin(sparam/r_tot);

	   r[0] = -r_tot*cos(omega);
	   r[1] = sparam;
	   r1_l[0]=0;
	   r1_l[1]=0;
	   r2_l[0]=r[0]*cos(-theta)+r[1]*sin(-theta);
	   r2_l[1]=-r[0]*sin(-theta)+r[1]*cos(-theta);
	   r_norm = r_tot;
	   ar[0] = 0;
	   ar[1] = 0;
	   V[0]=(masse_a*va[0]+epsilon*masse_b*vb[0])/M;
	   V[1]=(masse_a*va[1]+epsilon*masse_b*vb[1])/M;
	   R[0]=(epsilon*masse_b*r2_l[0])/M;
	   R[1]=(epsilon*masse_b*r2_l[1])/M;
       //velocity verlet:
	   do
	     { 
	       r[0]=r[0]+vr[0]*delta_t+0.5*ar[0]*delta_t*delta_t; //Berechnet neue Position
	       r[1]=r[1]+vr[1]*delta_t+0.5*ar[1]*delta_t*delta_t;
	      
	       r_norm = sqrt(r[0]*r[0]+r[1]*r[1]);

	       ar_norm2 =(epsilon*k*(1+delta)*(r_tot-r_norm))/red_masse; //Berechnet Beschleunigung	       
	       ar2[0]=(ar_norm2/r_norm)*r[0];
	       ar2[1]=(ar_norm2/r_norm)*r[1];

	       vr2[0]=vr[0]+0.5*(ar[0]+ar2[0])*delta_t; //Berechnet Geschwingkeit zu t+Dt
	       vr2[1]=vr[1]+0.5*(ar[1]+ar2[1])*delta_t;
	  
	       vr[0]=vr2[0];
	       vr[1]=vr2[1];
	       ar[0]=ar2[0];
	       ar[1]=ar2[1];
	  if(fl==0)
		 {
		   fwrite(&r[0], double_sz,1,fpb);
		   fwrite(&r[1], double_sz,1,fpb);
		   fwrite(&e_tot, double_sz,1,fpb);
		 }
	       // umrechnung ins laborsystem:
	   if (fl==0)
	     {
	       R[0]=R[0]+V[0]*delta_t;
	       R[1]=R[1]+V[1]*delta_t;

	       r_l[0]=r[0]*cos(-theta)+r[1]*sin(-theta);
	       r_l[1]=-r[0]*sin(-theta)+r[1]*cos(-theta);
	  
	       r1_l[0]=R[0]-(red_masse/masse_a)*r_l[0];
	       r1_l[1]=R[1]-(red_masse/masse_a)*r_l[1];
	  
	       r2_l[0]=R[0]+(red_masse/(epsilon*masse_b))*r_l[0];
	       r2_l[1]=R[1]+(red_masse/(epsilon*masse_b))*r_l[1];
	
	       fwrite(&r1_l[0], double_sz,1,fpc);
	       fwrite(&r1_l[1], double_sz,1,fpc);	
	       fwrite(&r2_l[0], double_sz,1,fpc);
	       fwrite(&r2_l[1], double_sz,1,fpc);
	       fwrite(&R[0], double_sz,1,fpc);
	       fwrite(&R[1], double_sz,1,fpc);
	     }
	  }while(r_tot-r_norm>=0);

      // berechnung der Geschwindigkeiten nach dem stoß:

	   v_l[0]=vr[0]*cos(-theta)+vr[1]*sin(-theta);
	   v_l[1]=-vr[0]*sin(-theta)+vr[1]*cos(-theta);
       
	   v1_l[0]=V[0]-(red_masse/masse_a)*v_l[0];
	   v1_l[1]=V[1]-(red_masse/masse_a)*v_l[1];
	  
	   v2_l[0]=V[0]+(red_masse/(epsilon*masse_b))*v_l[0];
	   v2_l[1]=V[1]+(red_masse/(epsilon*masse_b))*v_l[1];

	   delta_E1 = masse_a*0.5*(v1_l[0]*v1_l[0]+v1_l[1]*v1_l[1]-va[0]*va[0]-va[1]*va[1]);
	   delta_E2 = masse_b*0.5*(v2_l[0]*v2_l[0]+v2_l[1]*v2_l[1]-vb[0]*vb[0]-vb[1]*vb[1]);
	   delta_Etot=delta_E1+delta_E2;

	   sum_delta_E1 += delta_E1;
	   sum_delta_E2 += delta_E2;
	   sum_delta_Etot += delta_Etot;
	   
	   if (fl ==1 || counter ==0 )
	     {
	       fwrite(&delta_E1, double_sz,1,fpd);
	       fwrite(&delta_E2, double_sz,1,fpd);
	       fwrite(&delta_Etot, double_sz,1,fpd);
	     }
	   if (counter == max)
	     {
	       fwrite(&delta_E1, double_sz,1,fpf);
	       fwrite(&delta_E2, double_sz,1,fpf);
	       fwrite(&delta_Etot, double_sz,1,fpf);
	     }
       }

       mean_delta_E1 = sum_delta_E1 / (double) nbr_of_smpls;
       mean_delta_E2 = sum_delta_E2 / (double) nbr_of_smpls;
       mean_delta_Etot = sum_delta_Etot / (double) nbr_of_smpls;
   
       printf("mean_delta_E1=%F\n",mean_delta_E1);
       printf("mean_delta_E2=%F\n",mean_delta_E2);
       printf("mean_delta_Etot=%F\n",mean_delta_Etot);
       fwrite(&therm_en1,double_sz,1,fpe);
       fwrite(&therm_en2,double_sz,1,fpe);
       
       delta_E1=0;
       delta_E2=0;
       delta_Etot=0;
       sum_delta_E1 = 0;
       sum_delta_E2 =0;
       sum_delta_Etot =0;



       if (fl ==1 || counter ==0)
	 {
	       fclose(fpd);
	 }
       if ( counter == max)
	 {
	       fclose(fpf);
	 }
       if (fl != 2){break ;}
     }
   fclose(fpb);
   fclose(fpc);  
   fclose(fpe);
   fclose(fpg);

  return 0;
}
