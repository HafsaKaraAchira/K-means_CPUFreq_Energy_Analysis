/**
 * 
	Program to execute a CPU consuming phase 
	and a IO consuming phase on alternation
	
	executed with cmd : 
		clear && gcc -g program.c -o program -lm -D _GNU_SOURCE && ./program

	sudo su -
	cd /home/hafsa/Desktop/mini_tp
	(clear && gcc -g program.c -o program -lm -D _GNU_SOURCE && ./program) & (taskset -c 3 ./prog_script_cpu_io) & (taskset -c 2 watch -t -p -n 1 ./prog_script_freq)

 * 
 * **/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/resource.h>
#include <fcntl.h>
#include "libs/kmeans_utils.c"

#ifndef _GNU_SOURCE
#define	_GNU_SOURCE 1
#endif

#include <sched.h>

void save_phase_time(int loop, int phase , int iteration){
    char *cmd;
    asprintf(&cmd,"taskset -c 2 ./scripts/prog_script_phase %d %d %u",loop,phase,iteration);
	  system(cmd);
}

void generate_dataset_file(int dim_num, int point_num, int *seed , double * point , char *filename){
  r8mat_uniform_01 (dim_num,point_num,seed,point) ;
  r8mat_write (filename,dim_num,point_num,point) ;
  sleep(5);
  /*
  for ( int j = 0; j < point_num; j++ )
  {
    for ( int i = 0; i < dim_num; i++ ){printf ("  %24.16g", point[i+j*dim_num] );}
    printf ("\n" );
  }
  */
  //printf("Type any caractere to launch the clustering : ") ;
  //char c ;
  //scanf("%s",&c) ;
}




void kmeans_01 ( int dim_num, int point_num, int cluster_num, int it_max,
  int *it_num, double point[], int cluster[], double cluster_center[],
  int cluster_population[], double cluster_energy[] )

/******************************************************************************/
/*
  Purpose:

    KMEANS_01 applies the K-Means algorithm.

  Discussion:

    Given a matrix of POINT_NUM observations on DIM_NUM variables, the
    observations are to be allocated to CLUSTER_NUM clusters in such
    a way that the within-cluster sum of squares is minimized.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 October 2011

  Author:

    Original FORTRAN77 version by David Sparks.
    C++ version by John Burkardt.

  Reference:

    David Sparks,
    Algorithm AS 58:
    Euclidean Cluster Analysis,
    Applied Statistics,
    Volume 22, Number 1, 1973, pages 126-130.

  Parameters:

    Input, int DIM_NUM, the number of spatial dimensions.

    Input, int POINT_NUM, the number of points.

    Input, int CLUSTER_NUM, the number of clusters.

    Input, int IT_MAX, the maximum number of iterations.

    Output, int *IT_NUM, the number of iterations taken.

    Input, double POINT[DIM_NUM*POINT_NUM], the points.

    Output, int CLUSTER[POINT_NUM], indicates which cluster
    each point belongs to.

    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
    the cluster centers.

    Output, int CLUSTER_POPULATION[CLUSTER_NUM], the number
    of points in each cluster.

    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the
    cluster energies.
*/
{
  double dc;
  double de;
  double *f;
  int i;
  int il;
  int ir;
  int j;
  int j2;
  int k;
  double point_energy;
  double point_energy_min;
  int swap;

  *it_num = 0;
/*  Idiot checks.*/
  if ( cluster_num < 1 )
  {
    printf ( "\n" );
    printf ( "KMEANS_01 - Fatal error!\n" );
    printf ( "  CLUSTER_NUM < 1.\n" );
    exit ( 1 );
  }

  if ( dim_num < 1 )
  {
    printf ( "\n" );
    printf ( "KMEANS_01 - Fatal error!\n" );
    printf ( "  DIM_NUM < 1.\n" );
    exit ( 1 );
  }

  if ( point_num < 1 )
  {
    printf ( "\n" );
    printf ( "KMEANS_01 - Fatal error!\n" );
    printf ( "  POINT_NUM < 1.\n" );
    exit ( 1 );
  }
/*
  For each observation, calculate the distance from each cluster
  center, and assign to the nearest.
*/
  for ( j = 0; j < point_num; j++ ) // iterate on points
  {
    point_energy_min = HUGE_VAL;
    cluster[j] = -1;

    for ( k = 0; k < cluster_num; k++ ) // iterate on clusters
    {
      point_energy = 0.0;
      for ( i = 0; i < dim_num; i++ ) //iterate on dimensions
      {
        //calculate energy
        point_energy = point_energy +
          pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
      }

      if ( point_energy < point_energy_min )
      {
        point_energy_min = point_energy;  //update minimun point energy
        cluster[j] = k; //asign the point new cluster
      }
    }
  }
/*
  Determine the cluster population counts.
*/
  i4vec_zeros ( cluster_num, cluster_population );

  for ( j = 0; j < point_num; j++ ) //iterate on points
  {
    k = cluster[j];   //get the point cluter
    //increment population for retrieved cluster
    cluster_population[k] = cluster_population[k] + 1;
  }
/*
  Calculate the mean and sum of squares for each cluster.
*/
  r8vec_zeros ( dim_num * cluster_num, cluster_center );

  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    for ( i = 0; i < dim_num; i++ )
    {
      cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num]
        + point[i+j*dim_num];
    }
  }

  for ( k = 0; k < cluster_num; k++ )
  {
    if ( 0 < cluster_population[k] )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num] /
         ( double ) cluster_population[k];
      }
    }
  }
/*
  Set the point energies.
*/
  f = r8vec_zeros_new ( point_num );
  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    for ( i = 0; i < dim_num; i++ )
    {
      f[j] = f[j] + pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
    }
  }
/*
  Set the cluster energies.
*/
  r8vec_zeros ( cluster_num, cluster_energy );
  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    cluster_energy[k] = cluster_energy[k] + f[j];
  }
/*
  Adjust the point energies by a weight factor.
*/
  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];
    if ( 1 < cluster_population[k] )
    {
      f[j] = f[j] * ( double ) ( cluster_population[k] )
        / ( double ) ( cluster_population[k] - 1 );
    }
  }
/*
  Examine each observation in turn to see if it should be
  reassigned to a different cluster.
*/
  *it_num = 0;

  while ( *it_num < it_max )  //while not it_max reached
  {
    *it_num = *it_num + 1;

    swap = 0; //boolean if it needs swapping (non convergnece)

    for ( j = 0; j < point_num; j++ ) //iterate on point num
    {
      il = cluster[j];  //retrieve the point actual cluster
      ir = il;

      if ( cluster_population[il] <= 1 )  //if there is not point in cluster k , ignore
      {
        continue;
      }

      dc = f[j];  //retrieve the point-cluster energie (distance square)

      for ( k = 0; k < cluster_num; k++ ) //iterate on clusters
      {
        if ( k != il )  // for each cluter diffrent tha the point's actual cluster
        {
          de = 0.0;
          for ( i = 0; i < dim_num; i++ ) // iterate on dimensions
          {
            //calculate the distance square between the point and each other cluster
            de = de + pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
          }
          //weight the ditance with the other clusters population
          de = de * ( double ) cluster_population[k]
             / ( double ) ( cluster_population[k] + 1 );

          if ( de < dc )  // if the new distance is smaller thna the current one
          {
            //reassign the point curret energy
            dc = de;
            ir = k;
          }
        }
      }
/*
  If the lowest value was obtained by staying in the current cluster,
  then cycle.
*/
      if ( ir == il ){  continue;}
/*
  Reassign the point from cluster IL to cluster IR.
*/
      for ( i = 0; i < dim_num; i++ )
      {
        //update the current cluster center (remove the point)
        cluster_center[i+il*dim_num] = ( cluster_center[i+il*dim_num]
          * ( double ) ( cluster_population[il] ) - point[i+j*dim_num] )
          / ( double ) ( cluster_population[il] - 1 );

        //update the new cluster center (add the point)
        cluster_center[i+ir*dim_num] = ( cluster_center[i+ir*dim_num]
          * ( double ) ( cluster_population[ir] ) + point[i+j*dim_num] )
          / ( double ) ( cluster_population[ir] + 1 );
      }
      //update clusters (old and current) energies and population
      cluster_energy[il] = cluster_energy[il] - f[j];
      cluster_energy[ir] = cluster_energy[ir] + dc;
      cluster_population[ir] = cluster_population[ir] + 1;
      cluster_population[il] = cluster_population[il] - 1;
      
      //assign the cluster index to the point
      cluster[j] = ir;
/*
  Adjust the value of F for points in clusters IL and IR.
*/
      //iterate points
      for ( j2 = 0; j2 < point_num; j2++ )
      {
        //retirve the point clsuter
        k = cluster[j2];

        // if it is one of the udated clusters
        if ( k == il || k == ir )
        {
          f[j2] = 0.0;
          //update the point energy (after the clluster center change)
          for ( i = 0; i < dim_num; i++ )
          {
            f[j2] = f[j2] + pow ( point[i+j2*dim_num] - cluster_center[i+k*dim_num], 2 );
          }

          //multiply the energy by the new weight
          if ( 1 < cluster_population[k] )
          {
            f[j2] = f[j2] * ( double ) ( cluster_population[k] )
              / ( ( double ) ( cluster_population[k] - 1 ) );
          }
        }
      }
      //increment swap indicator
      swap = swap + 1;
    }
    /*
    Exit if no reassignments were made during this iteration.
    */
    if ( swap == 0 ){ break;}
  }

/*
  Compute the cluster energies.
*/
  r8vec_zeros ( cluster_num, cluster_energy );
  for ( j = 0; j < point_num; j++ )
  {
    k = cluster[j];

    point_energy = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      point_energy = point_energy +
        pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
    }
    cluster_energy[k] = cluster_energy[k] + point_energy;
  }

  free ( f );
  return;

}

/******************************************************************************/



/* Main program to test above function */
int main()
{

	cpu_set_t set;
	// clear cpu mask
	CPU_ZERO(&set);
	// set cpu 0
	CPU_SET(0, &set);
	// 0 is the calling process
	sched_setaffinity(0, sizeof(cpu_set_t), &set);
	//set priority
	setpriority(PRIO_PROCESS, 0, -20);

	// Initialization, should only be called once.
	srand(time(NULL));

  /************* PHASE ******************/

  int dim_num = 2;
  int point_num = file_row_count("./results/points.csv") ;
  int cluster_num = 10 ;
  int it_max = 24000000 , it_num = 0 ;
  double *point = ( double * ) malloc ( dim_num*point_num * sizeof ( double ) );
  int *cluster = ( int * ) malloc ( point_num * sizeof ( double ) ) ;
  double *cluster_center= ( double * ) malloc ( dim_num*point_num * sizeof ( double ) );
  int *cluster_population = ( int * ) malloc ( cluster_num * sizeof ( double ) ) ;
  double *cluster_energy = ( double * ) malloc ( cluster_num * sizeof ( double ) ) ;
  double * cluster_variance = ( double * ) malloc ( cluster_num * sizeof ( double ) ) ;
  
  int seed = 18 ;
  char *filename ="./results/points.csv" ;

  /*****************************************************/

  //generate_dataset_file(dim_num,point_num,&seed,point,filename) ;
  system("nohup ./scripts/prog_script_launch &")   ;

  /*****************************************************/

  point = r8mat_data_read (filename,dim_num,point_num) ;

  /*****************************************************/

  save_phase_time(1,1,0);
    cluster_initialize_1(dim_num,point_num,cluster_num,point) ;
  save_phase_time(1,1,1);

  /*****************************************************/

  save_phase_time(1,2,0);

    kmeans_01 (dim_num,point_num,cluster_num,it_max,&it_num,
            point,cluster,cluster_center,
            cluster_population,cluster_energy) ;
  
  save_phase_time(1,2,1);
  
  /*****************************************************/
  
  i4mat_write ("./results/clusters.csv",1,point_num,cluster) ;
  r8mat_write ("./results/centers.csv",dim_num,cluster_num,cluster_center) ;

  
  cluster_variance = cluster_variance_compute ( dim_num, point_num,cluster_num,
                      point,cluster, cluster_center) ;
          

  /*****************************************************/

  cluster_print_summary ( point_num, cluster_num,cluster_population, cluster_energy, cluster_variance ,it_num) ;

  //system("pkill watch");
	//system("killall -9 prog_script_cpu_io");

	return 0;
}

