// This code is available under GNU GPL 3 or later.
// By C. Hanuma Chaitanya and Angad S. Gill

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <math.h>
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
	int conf[52][52],e[3000];unsigned int l;
	int r,i,j,s,n,p,q,J=1;
	int a,b,z;
	float spin=0,kbT,m[3000],x[3000],h[3000];
	double r1,r2,energy,power,chi,ham,fef=0,ftef=0,ef=0,tef=0;

	FILE *obj;
	obj=fopen("ising.dat","w");

	clrscr ();
	randomize();
	printf ("enter the size of lattice:    ");
	scanf ("%d",&n);

	for (j=1;j<=n;j++)		/*random matrix*/
	{
	 	for (i=1;i<=n;i++)/*printing the first conf and calculating the first spin*/
		{

		 	r=rand();

		 	if (r%2==1)
	  		{
				printf (" 1");
				conf[j][i]=1;
				spin=spin+1;
			}

			else
			{
				printf ("-1");
				conf[j][i]=-1;
				spin=spin-1;
			}
		}
		printf ("\n");

	}
					/*SPIN CALCULATION*/
	s=n*n;
	printf("net spin of start config is %f\n",spin/s);

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
				ef=(conf[j-1][i]+conf[j][i+1]+conf[j+1][i]+conf[j][i-1])*conf[j][i] ;
				tef=tef+ef;
			}
		}
		tef=(tef/2)*-1*J;		/*total energy*/
		printf ("total enery at start(tef) %d(J)\n",tef);

	}

	printf("Working please wait..");

	for(kbT=5;kbT>0.05;kbT=kbT-0.05)
	{
		z=1;
		for(l=1;l<=10000;l++)
		{
			for(q=1;q<=n;q++)
			{
				for(p=1;p<=n;p++)
				{
					conf[q][p]=conf[q][p]*-1 ;	/*flipped*/

					if(p==1||p==n||q==1||q==n)
					{
						conf[q][0]=conf[q][n];	/*boundary conditions*/
						conf[q][n+1]=conf[q][1];
						conf[0][p]=conf[n][p];
						conf[n+1][p]=conf[1][p];
					}
					fef=(conf[q-1][p]+conf[q][p+1]+conf[q+1][p]+conf[q][p-1])*conf[q][p] ;
						/*energy cntrib of flipped pt*/

					fef=-2*fef;
					ftef=(tef+fef);		/*energy after flip*/

					/*basic initial check*/

					if (fef<=0)
					{
					 	/*if accepted changes*/
						tef=ftef;

					}
					else       /*second check compare exp wid random#*/
					{
						
						power=(double)fef/kbT;
						r1=exp(-1*power);
						/*r=random(32000);*/
						r2=(double)Random();

						if(r1>=r2)
							{
								tef=ftef;

							}
						else /*rejected no changes to config*/
							{
								conf[q][p]=conf[q][p]*-1 ;

							}
					}
				}

			}
			if (l>=200&&l%5==0)
				{
					spin=0;
					for(b=1;b<=n;b++)
					{
						for(a=1;a<=n;a++)
						{
							spin=conf[b][a]+spin;
						}
					}

					m[z]=(float)spin/s;
					e[z]=tef;
					z++;
				}
		}




   		spin=0;/*net spin*/
		energy=0;/*net energy*/
		chi=0;
		ham=0;/*sp heat*/
		z=z-1;

		for(a=1;a<=z;a++)
			{
				spin=m[a]+spin;
				energy=e[a]+energy;
			}

		spin=(float)spin/z;/*spin is avg of m[]*/
		energy=(float)energy/z;/*energy is avg of e[]*/

		for(a=1;a<=z;a++)
		{
			x[a]=m[a]-spin;
			h[a]=e[a]-energy;
			x[a]=x[a]*x[a];
			h[a]=h[a]*h[a];
			chi=x[a]+chi;
			ham=h[a]+ham;
		}

		

		chi=chi/(z*kbT);/*chi is avg of x[]^2*/
		ham=ham/(s*kbT*kbT*z);/*ham is avg of h[]^2*/

		fprintf(obj,"%f, %.8lf, %.8lf, %.8lf, %.8lf\n",kbT,spin,energy,chi,ham);/*file */



	}

	printf(" done!!!\a");

	fclose(obj);
	getch();
}
#include <stdio.h>
#include <time.h>

#define MODULUS    2147483647 /* DON'T CHANGE THIS VALUE                  */
#define MULTIPLIER 48271      /* DON'T CHANGE THIS VALUE                  */
#define CHECK      399268537  /* DON'T CHANGE THIS VALUE                  */
#define STREAMS    256        /* # of streams, DON'T CHANGE THIS VALUE    */
#define A256       22925      /* jump multiplier, DON'T CHANGE THIS VALUE */
#define DEFAULT    123456789  /* initial seed, use 0 < DEFAULT < MODULUS  */
      
static long seed[STREAMS] = {DEFAULT};  /* current state of each stream   */
static int  stream        = 0;          /* stream index, 0 is the default */
static int  initialized   = 0;          /* test for stream initialization */


   double Random(void)
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed 
 * between 0.0 and 1.0. 
 * ----------------------------------------------------------------
 */
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


