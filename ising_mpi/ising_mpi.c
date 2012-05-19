// This code is available under GNU GPL 3 or later
// By: S.V.S. Aditya and C. Hanuma Chaitanya

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include "mpi.h"
#include <math.h>

int main (argc, argv)
int argc; 
char **argv;

{

int my_PE_num,np;  
 MPI_Init(&argc, &argv);  
 MPI_Comm_rank(MPI_COMM_WORLD, &my_PE_num);  
  MPI_Comm_size(MPI_COMM_WORLD, &np);
 printf("Hello from %d.\n", my_PE_num);   



	int conf[52][52],e[3000];unsigned int l;
	int r,i,j,s,n,y,p,q,J=1,arr;
	int a,b,z;
        
	float ind_spin[20],spin=0,kbT,m[3000],x[3000],h[3000];
	double r11,r12,r21,r22,ind_energy[20],ind_power[20],energy,ind_fef[20],ind_tef[20],ind_ftef[20],ind_ef[20];

	FILE *obj;
	obj=fopen("ising.dat","w");
        srand ( time(NULL) );
      

       
	printf ("enter the size of lattice:    ");
	scanf ("%d",&n);
       	
        y=n/np;
        
     

printf("%d\n",my_PE_num); 

if(my_PE_num==0)


{for (j=1;j<=n;j++)		/*random matrix*/
	{
	 	for (i=1;i<=n;i++)/*printing the first conf and calculating the first spin*/
		{

		 	r=rand();
                      

		 	if (r%2==1)
	  		{
				printf (" 1");
				conf[j][i]=1;
                               MPI_Bcast(&conf[j][i],1,MPI_INT,0,MPI_COMM_WORLD);
				
			}

			else
			{
				printf ("-1");
				conf[j][i]=-1;
                               MPI_Bcast(&conf[j][i],1,MPI_INT,0,MPI_COMM_WORLD);
				
			}
		}
		printf ("\n");

	}}





printf("%d is the size\n",np);





for (j=1;j<=n;j++)		/*random matrix*/
	{
	 	for (i=((my_PE_num)*(y))+1;i<=((my_PE_num+1)*y);i++)/*printing the first conf and calculating the first spin*/
		{

		 	r=rand();
                      

		 	if (conf[j][i]==1)
	  		{
				printf (" 1");
				
				ind_spin[my_PE_num]=ind_spin[my_PE_num]+1;
                              
                              
			}

			else
			{
				printf ("-1");
				conf[j][i]=-1;
				ind_spin[my_PE_num]=ind_spin[my_PE_num]-1;
			}
		}
		printf ("\n");

	}
					/*SPIN CALCULATION*/
	s=n*y;
	printf("net spin of start config is %f\n",ind_spin[my_PE_num]/s);

	{
		for (i=((my_PE_num)*(y))+1;i<((my_PE_num+1)*y)+1;i++)
		{
			for (j=1;j<=n;j++)
			{	if(i==1||i==n||j==1||j==n)

				{	conf[j][0]=conf[j][n];
                                       MPI_Bcast(&conf[j][0],1,MPI_INT,my_PE_num,MPI_COMM_WORLD);
                                            /*boundary conditions*/
					conf[j][n+1]=conf[j][1];
                                       MPI_Bcast(&conf[j][n+1],1,MPI_INT,my_PE_num,MPI_COMM_WORLD);     
					conf[0][i]=conf[n][i];
                                       MPI_Bcast(&conf[0][i],1,MPI_INT,my_PE_num,MPI_COMM_WORLD);     
					conf[n+1][i]=conf[1][i];
                                       MPI_Bcast(&conf[n+1][i],1,MPI_INT,my_PE_num,MPI_COMM_WORLD);  
				}
				ind_ef[my_PE_num]=(conf[j-1][i]+conf[j][i+1]+conf[j+1][i]+conf[j][i-1])*conf[j][i] ;
				ind_tef[my_PE_num]=ind_tef[my_PE_num]+ind_ef[my_PE_num];
			}
		}
		ind_tef[my_PE_num]=(ind_tef[my_PE_num]/2)*-1;		/*total energy*/
		printf ("total enery at start(tef) %d\n",ind_tef[my_PE_num]);

	}

	printf("Working please wait..");

	for(kbT=5;kbT>0.05;kbT=kbT-0.05)
	{
		z=1;
		for(l=1;l<=10000;l++)
		{
			for(q=1;q<=n;q++)
			{
				for(p=((my_PE_num)*(y))+1;p<=((my_PE_num+1)*y);p++)
				{
					conf[q][p]=conf[q][p]*-1 ;	/*flipped*/

					
					if(p==1||p==n||q==1||q==n)
					{
						conf[q][0]=conf[q][n];	
                                               MPI_Bcast(&conf[q][0],1,MPI_INT,my_PE_num,MPI_COMM_WORLD);/*boundary conditions*/
						conf[q][n+1]=conf[q][1];
					       MPI_Bcast(&conf[q][n+1],1,MPI_INT,my_PE_num,MPI_COMM_WORLD);	
                                                conf[0][p]=conf[n][p];
                                                MPI_Bcast(&conf[0][p],1,MPI_INT,my_PE_num,MPI_COMM_WORLD);  
						conf[n+1][p]=conf[1][p];
                                               MPI_Bcast(&conf[n+1][p],1,MPI_INT,my_PE_num,MPI_COMM_WORLD);
					}
					ind_fef[my_PE_num]=(conf[q-1][p]+conf[q][p+1]+conf[q+1][p]+conf[q][p-1])*conf[q][p] ;
						/*energy cntrib of flipped pt*/

					ind_fef[my_PE_num]=-2*ind_fef[my_PE_num];
					ind_ftef[my_PE_num]=(ind_tef[my_PE_num]+ind_fef[my_PE_num]);		/*energy after flip*/

					/*basic initial check*/

					if (ind_fef[my_PE_num]<=0)
					{
					 	/*if accepted changes*/
						ind_tef[my_PE_num]=ind_ftef[my_PE_num];
                                                MPI_Bcast(&conf[j][i],1,MPI_INT,my_PE_num,MPI_COMM_WORLD);
                                              

					}
					else       /*second check compare exp wid random#*/
					{
						
						ind_power[my_PE_num]=(double)ind_fef[my_PE_num]/kbT;
						r12=exp(-1*ind_power[my_PE_num]);
						r=rand();
                                                r22=(double)r/32000;


                                                   

						if(r12>=r22)
							{
								ind_tef[my_PE_num]=ind_ftef[my_PE_num];
                                                                MPI_Bcast(&conf[j][i],1,MPI_INT,my_PE_num,MPI_COMM_WORLD);

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
					ind_spin[my_PE_num]=0;
					for(b=1;b<=n;b++)
					{
						for(a=((my_PE_num)*(y))+1;a<=((my_PE_num+1)*y);a++)
						{
					       ind_spin[my_PE_num]=conf[b][a]+ind_spin[my_PE_num];
						}
					}

					m[z]=(float)ind_spin[my_PE_num]/s;
					e[z]=ind_tef[my_PE_num];
					z++;
				}
		}




   		ind_spin[my_PE_num]=0;/*net spin*/
		ind_energy[my_PE_num]=0;/*net energy*/
		
	        z=z-1;

		for(a=1;a<=z;a++)
			{
				ind_spin[my_PE_num]=m[a]+ind_spin[my_PE_num];
				ind_energy[my_PE_num]=e[a]+ind_energy[my_PE_num];
			}

		ind_spin[my_PE_num]=(float)ind_spin[my_PE_num]/z;/*spin is avg of m[]*/
		ind_energy[my_PE_num]=(float)ind_energy[my_PE_num]/z;/*energy is avg of e[]*/

       
                                       if(my_PE_num!=0){MPI_Bcast(&ind_spin[my_PE_num],1,MPI_INT,my_PE_num,MPI_COMM_WORLD);
                                                         MPI_Bcast(&ind_energy[my_PE_num],1,MPI_INT,my_PE_num,MPI_COMM_WORLD);}
                for(arr=1;arr<=np;arr++){ 
                                       ind_spin[0]=ind_spin[arr]+ind_spin[0];

                                        ind_energy[0]=ind_energy[0]+ind_energy[arr];
                                       }
                
               ind_spin[0]=ind_spin[0]/np;
                 		

		
          fprintf(obj,"%f, %.8lf, %.8lf, \n",kbT,spin,energy);/*file */
               
                printf("%f, %.8lf, %.8lf, \n",kbT,ind_spin[0],ind_energy[0]);
           



	}

	printf(" done!!!\a");

	fclose(obj);
	getch();
MPI_Finalize ();
return 0;
} 












