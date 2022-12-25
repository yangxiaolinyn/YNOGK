
#include <mpi.h>
#include <stdlib.h>
#include "myfunctions.h" 
#include "functions.h" 
#include "ConstantsUnits.h"


extern long I97, J97;

int main(void)
{
	/*~~~~~~~~~~~~~ MPI Settings~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	int myid, np;
	int namelen;
	char processor_name[ MPI_MAX_PROCESSOR_NAME ];
	int argc;
	char** argv;
	MPI_Status status;

	printf(" jjjjjjjjjjjjjj = ");

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np );
	MPI_Comm_rank(MPI_COMM_WORLD, &myid );
	MPI_Get_processor_name(processor_name, &namelen );

	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	printf( " There are %d processors and\n", np );
	printf( " My id is : %d \n", myid );
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

	long long mydutyphot, Total_Smple_Num;
	Total_Smple_Num = 5.0 * pow(10, 7);
	if( myid == np-1 ) 
		mydutyphot = Total_Smple_Num / np + ( Total_Smple_Num % np );
	else
		mydutyphot = Total_Smple_Num / np;

	printf(" My duty photon number is : %lld \n", mydutyphot);
	/*~~~~~~~~~~~~~ End MPI Settings~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
 
	initRandom( -1, 0 );


	//~~~ Set initial parameters and conditions ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	double Normal = {zero};
	int m1 = 20;
	int m2 = 10;
	int const nn = 500;
	double am[m1 + 1];
	double bm[m2 + 1];
	double prob[m1 + 1];
	double Pr[m1 + 1]; 
	double xi[nn+2], Px[nn+2], x_var[nn+2];


	double delta_x = twopi / nn;


	if( myid == np-1 ) { 
		for( int i = 1; i <= m1; i++) { 
			am[i] = RANMAR(); 
			Normal += am[i]; 
		};

		for( int i = 1; i <= m2; i++) { 
			bm[i] = RANMAR();
			Normal += bm[i]; 
		};

		Normal *= 1.10;
		am[0] = one;
		for( int i = 1; i <= m1; i++)
			am[i] /= Normal;

		for( int i = 1; i <= m2; i++)
			bm[i] /= Normal; 

		prob[0] = one;
		for( int i = 1; i <= m1; i++) { 
			prob[0] -= am[i];
			prob[i] = am[i];
		};

		double pro = zero;
		pro = prob[0]; 
		Pr[0] = pro;
		for( int i = 1; i <= m1; i++) {
			pro += prob[i];
			Pr[i] = pro; 
		} 
		//MPI_Barrier( MPI_COMM_WORLD ); 
		//cout << " ass = " << am[m1] << endl;  
	} 

   
	MPI_Bcast( am, m1 + 1, MPI_DOUBLE_PRECISION, np-1, MPI_COMM_WORLD );
	MPI_Bcast( bm, m2 + 1, MPI_DOUBLE_PRECISION, np-1, MPI_COMM_WORLD ); 
	MPI_Bcast( prob, m1 + 1, MPI_DOUBLE_PRECISION, np-1, MPI_COMM_WORLD ); 
	MPI_Bcast( Pr, m1 + 1, MPI_DOUBLE_PRECISION, np-1, MPI_COMM_WORLD );

	printf(" jjjjjjjjjjjjjj = "); 
 

	double Phi1; 
	long long i_count = 0;
	do
	{ 
		Phi1 = RANMAR() * twopi;
          
		double p_1 = zero;
		double temp_phi;
		p_1 = RANMAR();
		for( int i = 0; i <= m1; i++) { 
			if( p_1 <= Pr[i] ) {  
				temp_phi = cosm_phi_sampling( i, Phi1, RANMAR() );
				break;
			} 
		}
 
		double F_phi = zero, H_phi = zero, x; 
		for( int k = 0; k <= m1; k++ )
			F_phi += am[k] * cos(k * temp_phi);
		for( int k = 1; k <= m2; k++ )
			H_phi += bm[k] * sin(k * temp_phi); 
 

		if( (RANMAR() - one) <= (H_phi / F_phi) )
			x = temp_phi;
		else
			x = twopi - temp_phi; 
 
		int i;
		i = floor(x / delta_x) + 1; 

 
		if ( i > nn )
			continue; 

		xi[i] += 1.0; 

		i_count += 1;
 
		if( i_count >= mydutyphot ) 
			break;

		if( myid == np-1 && i_count % 500000 == 0 )
			printf("J is = %lld \n i_count = %lld \n", i_count, mydutyphot); 

	}while(1);
 
	/*collect data from all other MPI processors.*/

	if ( myid != np-1 ) {
		if ( np-2 >= 0 ) {
			long send_num, send_tag;
			send_num  = nn + 2;
			send_tag = 1;
			MPI_Send( xi, send_num, MPI_DOUBLE_PRECISION, np-1, \
				send_tag, MPI_COMM_WORLD );
				printf("MPI processor: %d has send xi array to Processor %d\n"
					, myid, np-1);
		} 
	} else {
		if (np-2 >= 0) { 
			for( int RECV_SOURCE = 0; RECV_SOURCE < np - 1; RECV_SOURCE++ ) {
				long recv_num, recv_tag;
				double xi_RECV[nn + 2]; 
				recv_num = nn + 2;
				recv_tag = 1;
				MPI_Recv( xi_RECV, recv_num, MPI_DOUBLE_PRECISION, \
				RECV_SOURCE, recv_tag, MPI_COMM_WORLD, &status );
				for( int i = 0; i < nn + 2; i++ )
					xi[i] += xi_RECV[i];
			}
		}
 
		FILE *fp=fopen("./plot/dataFF_xi.txt", "w"); 
		FILE *fp1=fopen("./plot/dataFF_Px.txt", "w"); 
		FILE *fp2=fopen("./plot/dataFF_x.txt", "w"); 
		if ( fp != NULL ) { 
			printf("File open successful");  
			printf("Finished writing to file, will close now"); 
			for( int i = 0; i < nn + 2; i++ ) {
				x_var[i] = i * delta_x;
				Px[i] = F_cos_mphi( m1, x_var[i], am ) + F_sin_mphi( m2, x_var[i], bm );
				//myFile << xi[i] << "      " << x_var << "      " << Px << endl;
			}
			fwrite(xi, sizeof xi[1], nn+2, fp);
			fwrite(Px, sizeof Px[1], nn+2, fp1);
			fwrite(x_var, sizeof x_var[1], nn+2, fp2);
			//myFile.close();
			fclose(fp);
			fclose(fp1);
			fclose(fp2);
		}

	}


	MPI_Finalize(); 
	return 0;
}












