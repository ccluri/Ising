// This code is available under GNU GPL 3 or later
// This was written by C.Hanuma Chaitanya and Angad S. Gill 

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <math.h>
#define MODULUS    2147483647 
#define MULTIPLIER 48271      
#define CHECK      399268537  
#define STREAMS    256        
#define A256       22925      
#define DEFAULT    123456789  
#if !defined( _RNGS_ )
#define _RNGS_
double Random(void);
void   PlantSeeds(long x);
void   GetSeed(long *x);
void   PutSeed(long x);
void   SelectStream(int index);
void   TestRandom(void);

#endif

main()
{
	int conf[52][52],confp,e[2501];unsigned int l;
	int r,i,j,s,n,p,q,J=1;
	int a,b,z;
	float spin=0,kbT,m1[2501],m2[2501],x[2501],h[2501],spin1=0,spin2=0,spin3=0;
	double r1,r2,energy,power,chi,ham,fef=0,ftef=0,ef=0,tef=0;

	FILE *obj;
	obj=fopen("3potts.dat","w");

	clrscr ();
	randomize();
	printf ("enter the size of lattice:    ");
	scanf ("%d",&n);

	for (j=1;j<=n;j++)		/*random matrix*/
	{
	 	for (i=1;i<=n;i++)/*printing the first conf and calculating the first spin*/
		{
		 	r=rand();
		 	if (r%3==0)
	  		{
				printf (" 1");
				conf[j][i]=1;
				spin1=spin1+1;
			}
			else if (r%3==1)
			{
				printf (" 2");
				conf[j][i]=2;
				spin2=spin2-0.5;/*REAL PART*/
				spin3=spin3+0.8660;/*IMAGINARY PART*/
			}
			else
			{
				printf (" 3");
				conf[j][i]=3;
				spin2=spin2-0.5;/*REAL PART*/
				spin3=spin3-0.8660;/*IMAGINARY PART*/
			}
		}
		printf ("\n");
	}
					/*SPIN CALCULATION*/
	s=n*n;
	spin=spin1+spin2;
	printf("net spin of start config is %f+i%f\n",spin/s,spin3/s);
	{
		for(i=1;i<n+1;i++)
		{
			for (j=1;j<n+1;j++)
			{	if(i==1||i==n||j==1||j==n)
				{	conf[j][0]=conf[j][n];/*boundary conditions*/
					conf[j][n+1]=conf[j][1];
					conf[0][i]=conf[n][i];
					conf[n+1][i]=conf[1][i];
				}
				ef=0;
				if(conf[j-1][i]-conf[j][i]==0)
					{ef=1;}
				if(conf[j][i+1]-conf[j][i]==0)
					{ef=ef+1;}
				if(conf[j+1][i]-conf[j][i]==0)
					{ef=ef+1;}
				if(conf[j][i-1]-conf[j][i]==0)
					{ef=ef+1;}
				tef=tef+ef;
			}
		}
		tef=(tef/2)*-1;		/*total energy*/
		printf ("total enery at start(tef) %f\n",tef);
	}

	printf("Working please wait..");
	for(kbT=5;kbT>0.05;kbT=kbT-0.05)
	{
		z=1;
		for(l=1;l<=15000;l++)
		{
			for(q=1;q<=n;q++)
			{
				for(p=1;p<=n;p++)
				{	if(p==1||p==n||q==1||q==n)
					{
						conf[q][0]=conf[q][n];	/*boundary conditions*/
						conf[q][n+1]=conf[q][1];
						conf[0][p]=conf[n][p];
						conf[n+1][p]=conf[1][p];
					}
					ef=0;/*energy contrib of point before flip*/
					if(conf[q-1][p]-conf[q][p]==0)
						{ef=1;}
					if(conf[q][p+1]-conf[q][p]==0)
						{ef=ef+1;}
					if(conf[q+1][p]-conf[q][p]==0)
						{ef=ef+1;}
					if(conf[q][p-1]-conf[q][p]==0)
						{ef=ef+1;}
					ef=ef*-1;
					confp=conf[q][p];/*saving prev conf*/
	
					r=random(32000);/*flipped*/
					if (r%2==0)
					{
						conf[q][p]=conf[q][p]+1;	
						if (conf[q][p]==4)
						{conf[q][p]=1;}
					}

					else
					{
						conf[q][p]=conf[q][p]-1;
						if(conf[q][p]==0)
						{conf[q][p]=3;}
					}

					
					fef=0;/*energy contrib of point after flip*/
					if(conf[q-1][p]-conf[q][p]==0)
						{fef=1;}
					if(conf[q][p+1]-conf[q][p]==0)
						{fef=fef+1;}
					if(conf[q+1][p]-conf[q][p]==0)
						{fef=fef+1;}
					if(conf[q][p-1]-conf[q][p]==0)
						{fef=fef+1;}
					fef=fef*-1;
					ftef=tef+fef-ef;/*energy after flip*/						
				
					
					/*basic initial check*/

					if (ftef<=tef)
					{
					 	/*if accepted changes*/
						tef=ftef;

					}
					else       /*second check compare exp wid random#*/
					{
						
						power=(double)(ftef-tef)/kbT;
						r1=exp(-1*power);
						/*r=random(32000);*/
						r2=(double)Random();

						if(r1>=r2)
							{
								tef=ftef;

							}
						else /*rejected no changes to config*/
							{
								conf[q][p]=confp ;

							}
					}
				}

			}
			if (l>=200&&l%5==0)
				{
					spin=0;spin1=0;spin2=0;spin3=0;
					for(b=1;b<=n;b++)
					{
						for(a=1;a<=n;a++)
						{	if(conf[b][a]==1){spin1=spin1+1;}
							else if(conf[b][a]==2){	spin2=spin2-0.5;/*REAL PART*/
										spin3=spin3+0.8660;/*IMAGINARY PART*/}
							else{spin2=spin2-0.5;/*REAL PART*/
								spin3=spin3-0.8660;/*IMAGINARY PART*/}
									
							
						}
					}
					spin=spin2+spin1;
					m1[z]=(float)spin/s;m2[z]=(float)spin3/s;
					e[z]=tef;
					z++;
				}
		}



		spin=0;
   		spin2=0;spin3=0;/*net spin*/
		energy=0;/*net energy*/
		chi=0;
		ham=0;/*sp heat*/
		z=z-1;

		for(a=1;a<=z;a++)
			{
				spin2=m1[a]+spin2;
				spin3=m2[a]+spin3;
				energy=e[a]+energy;
			}

		spin2=(float)spin2/z;/*spin2 is avg of m1[]*/
		spin3=(float)spin3/z;/*spin3 is avg of m2[]*/
		

		energy=(float)energy/z;/*energy is avg of e[]*/
		
		spin=spin2*spin2;/*magnitude*/
		spin3=spin3*spin3;
		

		spin=spin+spin3;
		spin=(float)sqrt(spin);
		
		for(a=1;a<=z;a++)
		{
			x[a]=m1[a]-spin2;
			
			h[a]=e[a]-energy;
			x[a]=x[a]*x[a];
			h[a]=h[a]*h[a];
			chi=x[a]+chi;
			ham=h[a]+ham;
		}

		

		chi=chi/(z*kbT);/*chi is avg of x[]^2*/
		ham=ham/(s*kbT*kbT*z);/*ham is avg of h[]^2*/
		chi=chi*chi;
					
		fprintf(obj,"%f, %.8lf, %.8lf, %.8lf, %.8lf\n",kbT,spin,energy,chi,ham);/*file */
	}

	printf(" done!!!\a");

	fclose(obj);
	getch();
}
      
static long seed[STREAMS] = {DEFAULT};  /* current state of each stream   */
static int  stream        = 0;          /* stream index, 0 is the default */
static int  initialized   = 0;          /* test for stream initialization */


   double Random(void)
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
        long t;

  t = MULTIPLIER * (seed[stream] % Q) - R * (seed[stream] / Q);
  if (t > 0) 
    seed[stream] = t;
  else 
    seed[stream] = t + MODULUS;
  return ((double) seed[stream] / MODULUS);
}
