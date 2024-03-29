// This code is available under GNU GPL 3 or later
// This was writted by C. Hanuma Chaitanya and Y.Ravi Kiran

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
{	int conf[104][52],loop;unsigned int l;
	int r,i,j,s,n,p,q,J=1;
	int a,b,d,c,z,lp=10,l1,l2;
	float spin=0,kbT,m[3000],R[500],loopavg;
	double r1,r2,energy,power,chi,ham,fef=0,ftef=0,ef=0,tef=0,fef1,fef2;


	FILE *obj;
	obj=fopen("lattice.dat","w");

	clrscr ();



	randomize();



	printf ("enter no of boxes:    ");
	scanf ("%d",&n);
        c=n+1;
	d=2*n+2;
	for (j=1;j<=d;j++)		/*random matrix*/
	{
	 	for (i=1;i<=c;i++)/*printing the first conf and calculating the first spin*/
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
	s=(c)*(d);/*2n+2*/


	printf("net spin of start config is %f\n",spin/s);


	{
		for(j=1;j<=2*n+1;j=j+2)
		{
			for (i=1;i<=c;i++)
			{	if(i==c)

				{	/*boundary conditions*/
					conf[j+1][n+2]=conf[j+1][1];

				}
				if(j==2*n+1)

				{	/*boundary conditions*/
					conf[2*n+3][i]=conf[1][i];

				}
				ef=(conf[j][i]*conf[j+1][i]*conf[j+2][i]*conf[j+1][i+1]) ;
				tef=tef+ef;
				/*subtract n from this tef later contrib for ef=ef-`1*/
			}
		}



                tef=tef*-1;		/*total energy*/


		printf ("total enery at start(tef) %d(J)\n",tef);

	}

	printf("Working please wait..");

	for(kbT=2;kbT>0.05;kbT=kbT-0.05)
	{
		z=1;
		for(l=1;l<=10000;l++)
		{
			for(q=1;q<=d;q++)/*is it not 2n+2,we are not flipping the inbetween rows */
			{
				for(p=1;p<=c;p++)
				{
					conf[q][p]=conf[q][p]*-1 ;	
					
					if(p==c)

						{	/*boundary conditions*/
							conf[q+1][n+2]=conf[q+1][1];
					
						}

					if(q==2*n+1)/*each time it has to compute the value */


						{	/*boundary conditions*/
							conf[2*n+3][p]=conf[1][p];
						
						}
                                       
				
					
						if (q%2==0)
							{
								if (p!=1)
									{
										fef1=(conf[q-1][p-1]*conf[q][p-1]*conf[q+1][p-1]*conf[q][p]) ;
									 	fef2=(conf[q-1][p]*conf[q][p]*conf[q+1][p]*conf[q][p+1]) ; 
									}
									
								else
									{
										fef1=(conf[q-1][n+1]*conf[q][n+1]*conf[q+1][n+1]*conf[q][p]) ;
									 	fef2=(conf[q-1][p]*conf[q][p]*conf[q+1][p]*conf[q][p+1]) ;
					       				}
							}
						else
							{
								if (q!=1)
									{
										fef1=(conf[q][p]*conf[q+1][p]*conf[q+2][p]*conf[q+1][p+1]) ;
										fef2=(conf[q-2][p]*conf[q-1][p]*conf[q][p]*conf[q-1][p+1]) ;
									}
								else	{
										fef1=(conf[q][p]*conf[q+1][p]*conf[q+2][p]*conf[q+1][p+1]) ;
										fef2=(conf[2*n+1][p]*conf[d][p]*conf[1][p]*conf[d][p+1]) ;
									}
							}									
						/*energy cntrib of flipped pt*/

					fef=2*(fef1+fef2); /*chng in energy after flip*/

					ftef=(tef+fef);/*basic initial check*/

					if (fef<=0)
					{
					 	/*it accepted changes*/
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
			if (l>=2000&&l%100==0)
				{	/*lp=loop dimention also taking only five such loops*/
						loop=1;
						for (l1=1;l1<=lp;l1++)	
							{	
								loop=conf[2*lp+1][l1]*conf[1][l1]*loop;
							}	
						for (l2=2;l2<=2*lp;l2=l2+2)						
							{
								loop=conf[l2][1]*conf[l2][lp+1]*loop;
							}
						R[z]=loop;					
					/*printf(" chk pt 1 tef after 1 flip \n");*/
					spin=0;					

					for(b=1;b<=2*n+2;b++)/*2*n+2*/
					{
						for(a=1;a<=n+1;a++)/*n+1*/
						{
							spin=conf[b][a]+spin;
						}
					}
					
					m[z]=(float)spin/s;
					
					z++;
				}
		}




   		spin=0;/*net spin*/
		energy=0;/*net energy*/

		z=z-1;/*kk,so z=lmax /200+1*/

		
							
				for(a=1;a<=z;a++)
					{	
						loopavg=0;		
						loopavg=R[a]+loopavg;	
					}			
				
		
						
				loopavg=(float)loopavg/z; /*average spin on loop of size lp*/
							
		for(a=1;a<=z;a++)
			{
				spin=m[a]+spin;
				
			}

		spin=(float)spin/z;/*spin is avg of m[]*/
		

		fprintf(obj,"%f, %.8lf, %.8lf,\n",kbT,spin,loopavg);/*file */
		printf(".");


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


