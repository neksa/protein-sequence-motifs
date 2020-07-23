/**
 * Clustering  procedure
 *    Parallel implementation for use on HPC
 *    Coded in C programming language (C99 standard) and MPI
 *    Copyright (c) 2009 Alexandr Goncearenco
 *    Affiliation: CBU, BCCS, UNIFOB AS, UiB, Berezovsky Group.
 */

/*
INPUT:
    stored matrices in a text file: "output.matrix"
    counts are there and the initial parameters
    
    Parameters:
	rho	- this is threshold on K (number of matches in a profile
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <errno.h>
#include <error.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "mpi.h"
#include "PSSM.h"

#define MASTER 0

#include "clustering.h"

//#define DEBUG_DUMPALL 1

const char *source_matrices_filename = "output.matrix";
const char *unclustered_matrices_filename = "clustering.output/unclustered.matrix";
const char *clustered_matrices_filename = "clustering.output/clustered.matrix";
const char *iteration_clusters_filename = "clustering.output/iteration_clusters.tab";
const char *join_log_filename = "clustering.output/clustering_join_log.tab";

int main (int argc, char *argv[]) {
	int rank = 0;
	int size = 1;
	int rc = 0;
	
	rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		fprintf(stderr, "Error while initializing MPI\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		exit(1);
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int rho = 1;
	int a = 0;
	bool no_heuristics = false;
	bool reshuffle_input_matrices = false;

	// process command line parameters
	for (a=1; a<argc-1; a++) {
		if (argv[a] && !strcmp(argv[a], "-K")) {
			if (argv[a+1] && strlen(argv[a+1])) {
			    rho = atoi(argv[a+1]);
			}
		}

		if (argv[a] && !strcmp(argv[a], "-H")) {
			no_heuristics = true;
		}

		if (argv[a] && !strcmp(argv[a], "-r")) {
			reshuffle_input_matrices = true;
		}
	}


	#ifndef DNEBUG
	if (rank == MASTER) {
		fprintf(stderr, "Program started with %d CPU\n", size);
		fprintf(stderr, "OPTION: Rho (minimal K for matrix) = %d\n", rho);
		fprintf(stderr, "OPTION: no_heuristics = %d\n", (char)no_heuristics);
		fprintf(stderr, "OPTION: reshuffle input matrices = %d\n", (char)reshuffle_input_matrices);
	}
	#endif		

	int i=0, j=0, k=0, l=0, m=0;
	double (*matrices)[50][26] = NULL;
	double (**pmatrices)[50][26] = &matrices;
	double (*logodds)[50][26] = NULL;
	unsigned int *counts = NULL;
	unsigned int **pcounts = &counts;
	int N = 0;
	double composition[26];
	
	for (i=0; i<26; i++) {
		composition[i] = 0.0;
	}

	if (rank == MASTER) {
		mkdir("clustering.output", 0755);
		
		load_Composition(composition);
	
		N = load_Matrices(source_matrices_filename, pmatrices, pcounts, rho, reshuffle_input_matrices);
		#ifndef DNEBUG
		fprintf(stderr, "[rank %d] With rho=%d loaded %d matrices\n", rank, rho, N);
		#endif

		logodds = (double (*)[50][26])calloc(2*N, sizeof(double)*50*26);
		assert(N == 0 || logodds != NULL);

		// for each matrix
		for (m=0;m<N;m++) {
			calc_LogOdds_Matrix(composition, matrices[m], logodds[m]);
		}

		#ifndef DNEBUG
		fprintf(stderr, "[rank %d] LogOdds matrices calculated for %d matrices \n", rank, N);
		#endif
		
		if (N > 0) {
			write_Freq_Matrices(unclustered_matrices_filename, matrices, NULL, counts, N, false);
			#ifndef DNEBUG
			fprintf(stderr, "[rank %d] written unclustered matrices to file: %s\n", rank, unclustered_matrices_filename);
			#endif
		}
	}

	MPI_Bcast(&N, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	#ifndef DNEBUG
	fprintf(stderr, "[rank %d] Broadcasted N: %d\n", rank, N);
	#endif
	
	if (N == 0) {
		MPI_Abort(MPI_COMM_WORLD, rc);
		return 0;
	}

	MPI_Bcast(&composition, 26, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

	if (rank != MASTER) {
		matrices = (double (*)[50][26])ecalloc(2*N, sizeof(double)*50*26, "main space for matrices");
		assert(matrices != NULL);
		logodds = (double (*)[50][26])ecalloc(2*N, sizeof(double)*50*26, "main space for logodds matrices");
		assert(logodds != NULL);
		counts = (unsigned int *)ecalloc(2*N, sizeof(unsigned int), "counts");
		assert(counts != NULL);
	}
	
	MPI_Bcast(matrices, N*50*26, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	#ifndef DNEBUG
	fprintf(stderr, "[rank %d] Broadcasted matrices: %d\n", rank, N);
	#endif
	MPI_Bcast(logodds, N*50*26, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	#ifndef DNEBUG
	fprintf(stderr, "[rank %d] Broadcasted logodds: %d\n", rank, N);
	#endif
	MPI_Bcast(counts, N, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
	#ifndef DNEBUG
	fprintf(stderr, "[rank %d] Broadcasted counts: %d\n", rank, N);
	#endif

	/*
	if (rank != MASTER) {
		write_Freq_Matrices("unclustered.matrix2", matrices, counts, N);	
		#ifndef DNEBUG
		fprintf(stderr, "[rank %d] written matrices to file: %s\n", rank, unclustered_matrices_filename);
		#endif
	}
	*/
	
	#ifdef DEBUG_DUMPALL
	FILE *fd_dbg;
	#endif

	FILE *join_log = NULL;
	if (rank == MASTER) {
		join_log = fopen(join_log_filename, "w");
		fprintf(join_log, "iteration\taction\tc\ta\tb\tdistance\tshift\twindow_shift\tSI_c\tSI_a\tSI_b\tKL_c\tKL_a\tKL_b\tSI_complete_c\tSI_complete_a\tSI_complete_b\tKL_complete_c\tKL_complete_a\tKL_complete_b\n");
		fclose(join_log);		
	}
	
	
	double **DM;		// distances matrix:	 2N*2N = 4N^2 = 64N^2 bytes
	char **SM; 		// shifts matrix:	 2N*2N = 4N^2 = 4N^2 bytes
	char **WM; 		// window shifts matrix: 2N*2N = 4N^2 = 4N^2 bytes
	char **CM;     		// coloring matrix: iteration-clusters: N*2N = 2N^2 = 2N^2 bytes
	double *DITV;		// distance-iteration vector: N = 4N bytes
	double distance = 0.0;			

	
	// allocate memory:
	#ifndef DNEBUG
	fprintf(stderr, "[rank %d] Allocating memory for DM, SM, WM, CM matrices\n", rank);
	#endif
	DM = eppcallocf(2*N, 2*N, "DM pointer-to-pointer-to-double matrix");
	SM = eppcalloc(2*N, 2*N, "SM pointer-to-pointer-to-char matrix");
	WM = eppcalloc(2*N, 2*N, "WM pointer-to-pointer-to-char matrix");
	CM = eppcalloc(  N, 2*N, "CM pointer-to-pointer-to-char matrix");
	DITV = ecalloc(N, sizeof(double), "Distance on iteration pointer-to-double vector");
	
	fprintf(stderr, "[rank %d] Broadcasted counts: %d\n", rank, N);
	
	// Assume that MAX(iterations) = N

	for (j=0;j<N;j++) {
		DITV[j] = 0.0;
		CM[0][j] = 'x';
	}
	
	for (j=N;j<2*N;j++) {
		CM[0][j] = '.';
	}

	for (j=2*N;j--;) {
		for (i=1;i<N;i++) {
			CM[i][j] = '.';
		}
		for (i=2*N;i--;) {
			DM[i][j] = 0.0;
			SM[i][j] = 0;
			WM[i][j] = 0;
		}
	}	
	
	// calculate lower half-matrix
	for (k=rank; k<N; k+=size) {
		for (l=0; l<k; l++) {
			//distance = calc_Distance(logodds[k], logodds[l], &(DM[k][l]), &(SM[k][l]), &(WM[k][l]), composition);
			distance = calc_Distance(matrices[k], matrices[l], &(DM[k][l]), &(SM[k][l]), &(WM[k][l]), composition);

			#ifndef DNEBUG
			//fprintf(stderr, "[rank %d] Distance d(%d,%d) = {%f, %d, %d}\n", rank, k, l, DM[k][l], SM[k][l], WM[k][l]);
			#endif
		}
		#ifndef DNEBUG
		if (rank == MASTER) {
			fprintf(stderr, "[rank %d] Calculated distance d(%d,?)\n", rank, k);
		}
		#endif
	}

	// distribute lower half-matrices
	for (k = 1; k<N; k++) {
		#ifndef DNEBUG
		//fprintf(stderr, "[rank %d] k=%d, size=%d, broadcasting from rank k%%size=%d\n", rank, k, size, k%size);
		#endif
		MPI_Bcast(DM[k], k, MPI_DOUBLE, k%size, MPI_COMM_WORLD);
		MPI_Bcast(SM[k], k, MPI_CHAR, k%size, MPI_COMM_WORLD);
		MPI_Bcast(WM[k], k, MPI_CHAR, k%size, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	// each CPU does copy-transposed: copies lower half-matrix to higher half-matrix
	for (i=1; i<N; i++) {
		for (j=0; j<i; j++) {
			DM[j][i] = DM[i][j];
			SM[j][i] = 0 - SM[i][j];
			WM[j][i] = WM[i][j];			
		}
	}
	
	#ifndef DNEBUG
	fprintf(stderr, "[rank %d] Calculated and synchronized distance matrices\n", rank);
	#endif
	
	
	if (rank == MASTER) {
		write_Information_Distribution("clustering.output/information_distribution.initial.tab", composition, matrices, NULL, N);
		/*
		fd_dbg = fopen("clustering.0.debug", "w");
		for (i=0; i<N; i++) {
			for (j=0; j<N; j++) {
				fprintf(fd_dbg, "DM[%d][%d] = %f ; SM[%d][%d] = %d ; WM[%d][%d] = %d ; CM[%d][%d] = %c \n", i, j, DM[i][j], i,j, SM[i][j], i,j, WM[i][j], i,j, CM[i][j]);
			}
		}
		fclose(fd_dbg);
		*/
	}
	
	int iteration = 0;
	int N_joins = 0;
	int N_joins_now = 0;
	int N_discarded_now = 0;
	int done = 0;
	double mean = 0.0;
	double median = 0.0;
	double mean_of_squares = 0.0;
	double stddev = 0.0;
	double zscore = 0.0;
	double distance_threshold = 0.0;
	double min_min_distance = 0.0;
	double max_min_distance = 0.0;
	double min_max_distance = 0.0;
	double max_max_distance = 0.0;
	
	while (!done) {
		if (rank == MASTER) {
			// calculate mean of distances on current iteration
			// sort distances and calculate median on current iteration
			// calculate variance of distances
			// threshold would be distance with < 2 standard deviations from the mean
		
			char iteration_dir_name[50];
			snprintf(iteration_dir_name, 50, "clustering.output/%d", iteration);
			mkdir(iteration_dir_name, 0755);
		
			double *ordered_distances;
			ordered_distances = (double *)ecalloc(((N-1)*N)/2, sizeof(double), "ordered_distances vector");
			mean = 0.0;
			mean_of_squares = 0.0;
			stddev = 0.0;
			zscore = 0.0;
			//ordered_distances = calloc(N-iteration, sizeof(double));
			
			#ifdef DEBUG_DUMPALL
			char fname_dbg[50];
			snprintf(fname_dbg, 50, "clustering.output/%d/clustering.%d.CM.debug", iteration);
			fd_dbg = fopen(fname_dbg, "w");
			#endif
			i = 0;			
			for (k = 0; k< 2*N; k++) {
				if ('x' == CM[iteration][k]) {
					for (l = 0; l < k; l++) {
						if ('x' == CM[iteration][l]) {
							ordered_distances[i++] = DM[k][l];
							//printf("[rank %d] iteration=%d: i=%d d(%d, %d) = %f \n", rank, iteration, i, k, l, DM[k][l]);
							#ifdef DEBUG_DUMPALL
							fprintf(fd_dbg, "[rank %d] iteration=%d: i=%d d(%d, %d) = %f \n", rank, iteration, i-1, k, l, DM[k][l]);
							#endif
							mean += DM[k][l];
							mean_of_squares += DM[k][l] * DM[k][l];
						}
					}
				}
			}
			
			if (i < 1) {
			    done = 1;
			    break;
			}
			
			#ifdef DEBUG_DUMPALL
			for (k=0;k<2*N;k++){
			    fprintf(fd_dbg, "i=%d CM[%d][%d]=%d\n", i,  iteration, k, CM[iteration][k]);
			}
			fclose(fd_dbg);
			#endif

			mean /= i;
			mean_of_squares /= i;
			stddev = sqrt(mean_of_squares - mean*mean);
			
			qsort(ordered_distances, i, sizeof(double), compare_doubles);
			median = ordered_distances[i/2];
			free(ordered_distances);
			
			//for(j=0; j<i; j++) {
			//	fprintf(stderr, "[rank %d] Iteration=%d: ordered_distance[%d]=%f\n", rank, iteration, j, ordered_distances[j]);
			//}
			
			#ifndef DNEBUG
			fprintf(stderr, "[rank %d] Iteration = %d, mean=%f, median=%f, sd=%f \n", rank, iteration, mean, median, stddev);
			#endif
						
			char temp_fname[50];
			FILE *fd, *fd2;
			
			snprintf(temp_fname, 50, "clustering.output/%d/distances.%d.tab", iteration, iteration);
			fd = fopen(temp_fname, "w");
			assert(fd != NULL);
			fprintf(fd, "M1\tM2\tdistance\tdistance_zscore\n");

			// write half-matrix
			for (k = 0; k< 2*N; k++) {
				for (l = 0; l < k; l++) {
					//fprintf(stderr, "CM[%d]k[%d]=%d, CM[%d]l[%d]=%d\n", iteration, k, CM[iteration][k], iteration, l, CM[iteration][l]);
					if (CM[iteration][k] == 'x' && CM[iteration][l] == 'x') {
						fprintf(fd, "%d\t%d\t%.6f\t%.6f\n", k, l, DM[k][l], (DM[k][l] - mean) / stddev);
					}
				}
			}
			fclose(fd);
			
			// write min distance statistics:
			
			snprintf(temp_fname, 50, "clustering.output/%d/distances_min.%d.tab", iteration, iteration);
			fd = fopen(temp_fname, "w");
			assert(fd != NULL);
			fprintf(fd, "M1\tM2\tmin_distance\tmin_distance_zscore\n");

			snprintf(temp_fname, 50, "clustering.output/%d/distances_max.%d.tab", iteration, iteration);
			fd2 = fopen(temp_fname, "w");
			assert(fd2 != NULL);
			fprintf(fd2, "M1\tM2\tmax_distance\tmax_distance_zscore\n");

			// write min value for each row in full distance matrix
			double min_value, max_value;
			int min_k, min_l, max_k, max_l;
			min_min_distance = 1000000.0; //+inf
			max_min_distance = 0.0; // abs low
			max_max_distance = 0.0;
			min_max_distance = 1000000.0;

			for (k = 0; k< 2*N; k++) {
				min_value = 1000000.0; //+inf
				min_k=0;
				min_l=0;
				
				max_value = 0.0; //+inf
				max_k=0;
				max_l=0;
				
				for (l = 0; l < 2*N; l++) {
					if ( k != l && CM[iteration][k] == 'x' && CM[iteration][l] == 'x') {
						if (min_value > DM[k][l]) {
							min_value = DM[k][l];
							min_k = k;
							min_l = l;							
						}
						if (max_value < DM[k][l]) {
							max_value = DM[k][l];
							max_k = k;
							max_l = l;							
						}
					}
				}
				if (min_value != 1000000.0) {
					if (max_min_distance < min_value) {
						max_min_distance = min_value;
					}
					if (min_min_distance > min_value) {
						min_min_distance = min_value;
					}

					fprintf(fd, "%d\t%d\t%.6f\t%.6f\n", min_k, min_l, min_value, (min_value - mean) / stddev);
				}
				if (max_value != 0.0) {
					if (max_max_distance < max_value) {
						max_max_distance = max_value;
					}
					if (min_max_distance > max_value) {
						min_max_distance = max_value;
					}

					fprintf(fd2, "%d\t%d\t%.6f\t%.6f\n", max_k, max_l, max_value, (max_value - mean) / stddev);
				}
			}
			fclose(fd);
			fclose(fd2);

			#ifndef DNEBUG
			fprintf(stderr, "[rank %d] Calculated statistics on iteration %d!\n", rank, iteration);
			#endif
		} //Master
				
		MPI_Bcast(&mean, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
		MPI_Bcast(&median, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
		MPI_Bcast(&stddev, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

		MPI_Bcast(&min_min_distance, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
		MPI_Bcast(&max_min_distance, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

		MPI_Bcast(&min_max_distance, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
		MPI_Bcast(&max_max_distance, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
		
		if (no_heuristics) {
			distance_threshold = min_min_distance;
		} else {
			distance_threshold = min_min_distance + (max_min_distance - min_min_distance) / 20.0;
		}

		if (rank == MASTER) {
			#ifndef DNEBUG
			fprintf(stderr, "[rank %d] [iteration %d]  min_min=%f, max_min=%f, low 5%% distance threshold = %f.  min_max=%f, max_max=%f!\n", 
					rank, iteration, min_min_distance, max_min_distance, distance_threshold, min_max_distance, max_max_distance);
			#endif
		}
		
		// assume we keep what we have
		for (k = 0; k < 2*N; k++) {
			if (CM[iteration][k] == 'x') {
				CM[iteration+1][k] = 'x';
			}
		}
		
		N_joins_now = 0;
		N_discarded_now = 0;
		int min_k = -1;
		int min_l = -1;	
		
		
		// that's an openhearted :) (i.e. not greedy) algorithm to join the 5% segment of min_min distribution left tail
		// joins iteratively an absolutely minimal pair of matrices until the tail segment is completely processed
		if (rank == MASTER) {
			join_log = fopen(join_log_filename, "a");
		}
		
		while (1) {
			// not to use as an index without checking first!
			min_k = -1; 
			min_l = -1; 
			
			// find such (k, l) pair that DM[k][l] < distance_threshold and DM[k][l] is globally minimal
			for (k = 0; k < 2*N; k++) {
				if (CM[iteration][k] == 'x') {
					CM[iteration+1][k] = 'x';
					for (l = k+1; l < 2*N; l++) {
						if (CM[iteration][l] == 'x') {
							//zscore = (DM[k][l] - mean) / stddev;
							//zscore = (DM[k][l] - median) / stddev;
							//if (zscore < -0.55)
							if (DM[k][l] <= distance_threshold) {
								if (min_k < 0 || min_l < 0 || DM[k][l] < DM[min_k][min_l]) {
									min_k = k;
									min_l = l;
								}
							}
						}
					}
				}
			}
			
			// stop if global minimum not found
			// i.e. there exists no (k, l) pair that satisfies DM[k][l] < distance_threshold
			if (min_k < 0 || min_l < 0) {
				break;
			}
			
			
			int k_shift, l_shift;
			double SI_k, SI_l, SI_new, KL_k, KL_l, KL_new;
			double H_k_pos, H_l_pos, H_new_pos, new_freq;
			double SI_complete_k, SI_complete_l, SI_complete_new, KL_complete_k, KL_complete_l, KL_complete_new; 
			k_shift = (SM[min_k][min_l] > 0)?SM[min_k][min_l]:0;
			k_shift += WM[min_k][min_l];
			l_shift = (SM[min_k][min_l] > 0)?0:0-SM[min_k][min_l];
			l_shift += WM[min_k][min_l];
						
			SI_k = 0.0;
			SI_l = 0.0;
			SI_new = 0.0;
			//calculate average Kullback–Leibler divergence for K matrix
			KL_k = 0.0;
			KL_l = 0.0;
			KL_new = 0.0;
			for (i=0; i<30; i++) {
				H_k_pos = 0.0;
				H_l_pos = 0.0;
				H_new_pos = 0.0;
				for (j=0;j<26;j++) {
					if (composition[j] > 0.0) {
						if (matrices[min_k][i+k_shift][j] > 0.0) {
							KL_k += matrices[min_k][i+k_shift][j] * (log(matrices[min_k][i+k_shift][j]/composition[j])/log(2));
							H_k_pos += matrices[min_k][i+k_shift][j] * (log(matrices[min_k][i+k_shift][j])/log(2));
						}
						if (matrices[min_l][i+l_shift][j] > 0.0) {
							KL_l += matrices[min_l][i+l_shift][j] * (log(matrices[min_l][i+l_shift][j]/composition[j])/log(2));
							H_l_pos += matrices[min_l][i+l_shift][j] * (log(matrices[min_l][i+l_shift][j])/log(2));
						}
						new_freq = (matrices[min_k][i+k_shift][j] + matrices[min_l][i+l_shift][j]) / 2.0;
						if (new_freq > 0.0) {
							KL_new += new_freq * (log(new_freq/composition[j])/log(2));
							H_new_pos += new_freq * (log(new_freq)/log(2));
						}
					}
				}
				
				SI_k += log(20)/log(2) + H_k_pos;
				SI_l += log(20)/log(2) + H_l_pos;
				SI_new += log(20)/log(2) + H_new_pos;
			}
			SI_k /= 30.0;
			SI_l /= 30.0;
			SI_new /= 30.0;

			KL_k /= 30.0;
			KL_l /= 30.0;
			KL_new /= 30.0;

			if (rank == MASTER) {
				printf("[iteration %d] Kullback–Leibler divergence KL_k(%d)=%f, KL_l(%d)=%f, KL_new(%d)=%f. d(k,l) = %f. shift=%d, window_shift=%d, k_shift=%d, l_shift=%d\n",
					iteration, min_k, KL_k, min_l, KL_l, N + N_joins, KL_new, DM[min_k][min_l], SM[min_k][min_l], WM[min_k][min_l], k_shift, l_shift);
					
				//fprintf(join_log, "%d\t%d\t%d\t%d\t%f\t%d\t%d\n", iteration, N + N_joins, min_k, min_l, DM[min_k][min_l], SM[min_k][min_l], WM[min_k][min_l]);
			}
			
			// if new KL is less than the worst KL of l and k then discard k or l, whatever has lowest KL
			// else: close K, L, and leave NEW
			SI_complete_k = calculate_Shannon_information(composition, matrices[min_k]);
			SI_complete_l = calculate_Shannon_information(composition, matrices[min_l]);
			KL_complete_k = calculate_KL_divergence(composition, matrices[min_k]);
			KL_complete_l = calculate_KL_divergence(composition, matrices[min_l]);
			
			/* do not discard matrices
			
			if (KL_new < MIN(KL_k, KL_l)) {
				int discarded_profile = 0;
				int survived_profile = 0;
			
				if (KL_k < KL_l) {
					discarded_profile = min_k;
					survived_profile = min_l;
				} else {
					discarded_profile = min_l;
					survived_profile = min_k;
				}
				
				CM[iteration][discarded_profile] = '#'; 
				CM[iteration+1][discarded_profile] = '.';
				
				if (rank == MASTER) {
					printf("[iteration %d] d(%d,%d) = %f. shift=%d, window_shift=%d\n ",
						iteration, min_k, min_l, DM[min_k][min_l], SM[min_k][min_l], WM[min_k][min_l]);
					
						fprintf(join_log, "%d\t%s\t%d\t%d\t%d\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
						iteration, "D", survived_profile, min_k, min_l, DM[min_k][min_l], SM[min_k][min_l], WM[min_k][min_l],
						0.0, SI_k, SI_l, 0.0, KL_k, KL_l,
						0.0, SI_complete_k, SI_complete_l, 0.0, KL_complete_k, KL_complete_l);
				}		

				N_discarded_now++;
				continue;
			}
			*/
			
			// procede by joining k & l -> c
			
			join_Matrices(matrices[min_k], counts[min_k], matrices[min_l], counts[min_l], matrices[N + N_joins], SM[min_k][min_l], WM[min_k][min_l]);
			counts[N + N_joins] = MAX(counts[min_k], counts[min_l]);
			calc_LogOdds_Matrix(composition, matrices[N + N_joins], logodds[N + N_joins]);
			//MPI_Bcast(&(matrices[N + N_joins]), 50*26, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

			SI_complete_new = calculate_Shannon_information(composition, matrices[N+N_joins]);
			KL_complete_new = calculate_KL_divergence(composition, matrices[N+N_joins]);

			if (rank == MASTER) {
				printf("[iteration %d] Joining min d(%d,%d) = %f. shift=%d, window_shift=%d\n ",
					iteration, min_k, min_l, DM[min_k][min_l], SM[min_k][min_l], WM[min_k][min_l]);
					
					fprintf(join_log, "%d\t%s\t%d\t%d\t%d\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
					iteration, "J", N + N_joins, min_k, min_l, DM[min_k][min_l], SM[min_k][min_l], WM[min_k][min_l],
					SI_new, SI_k, SI_l, KL_new, KL_k, KL_l,
					SI_complete_new, SI_complete_k, SI_complete_l, KL_complete_new, KL_complete_k, KL_complete_l);
			}

			// mark pair as joined on current iteration:
			CM[iteration][min_k] = 'X'; // note capital X
			CM[iteration][min_l] = 'X'; // note capital X
			// remove joined matrices from the next iteration:
			CM[iteration+1][min_k] = '.';
			CM[iteration+1][min_l] = '.';
			// add a new resulting matrix to the next iteration:
			CM[iteration+1][N + N_joins] = 'x';

			N_joins_now++;
			N_joins++;
		}
		
		DITV[iteration] = distance_threshold;
		if (rank == MASTER) {
			
			fclose(join_log);
			
			#ifndef DNEBUG
			fprintf(stderr, "[rank %d] [iteration %d]  Joined %d matrices with distance threshold = %f!\n", rank, iteration, 2*N_joins_now, distance_threshold);
			#endif
		}

		MPI_Bcast(&N_joins_now, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (N_joins_now == 0 && N_discarded_now == 0) { // no joins on this iteration
			done = 1;
		} else {
		
			// recalculate lower half-matrix
			#ifndef DNEBUG
			if (rank == MASTER) {
				fprintf(stderr, "[rank %d] Recalculating lower half-matrix of distances with new matrices for the next iteration\n", rank);
			}
			#endif
			for (k=N+rank; k<2*N; k+=size) {
				if (CM[iteration+1][k] == 'x' && CM[iteration][k] == '.') {
					for (l=0; l<k; l++) {
						if (CM[iteration+1][l] == 'x') {
							distance = calc_Distance(logodds[k], logodds[l], &(DM[k][l]), &(SM[k][l]), &(WM[k][l]), composition);
						}
					}
				}
			}

			MPI_Barrier(MPI_COMM_WORLD);

			#ifndef DNEBUG
			if (rank == MASTER) {
				fprintf(stderr, "[rank %d] Distributing recalculated distances\n", rank);
			}
			#endif

			// distribute lower half-matrices
			for (k = N; k<2*N ; k++) {
				if (CM[iteration+1][k] == 'x' && CM[iteration][k] == '.') {
					#ifndef DNEBUG
					//fprintf(stderr, "[rank %d] k=%d, size=%d, broadcasting from rank k%%size=%d\n", rank, k, size, k%size);
					#endif
					MPI_Bcast(DM[k], k, MPI_DOUBLE, (k-N)%size, MPI_COMM_WORLD);
					MPI_Bcast(SM[k], k, MPI_CHAR, (k-N)%size, MPI_COMM_WORLD);
					MPI_Bcast(WM[k], k, MPI_CHAR, (k-N)%size, MPI_COMM_WORLD);
				}
			}

			MPI_Barrier(MPI_COMM_WORLD);
	
			// each CPU does copy-transposed: copies lower half-matrix to higher half-matrix
			for (i=1; i<2*N; i++) {
				for (j=0; j<i; j++) {
					DM[j][i] = DM[i][j];
					SM[j][i] = 0 - SM[i][j];
					WM[j][i] = WM[i][j];
				}
			}
		
			MPI_Barrier(MPI_COMM_WORLD);
			
			if (rank == MASTER) {
				#ifndef DNEBUG
				fprintf(stderr, "[rank %d] [iteration %d] Recalculated and distributed distances for new matrices.\n", 
						rank, iteration);
				#endif
			}				
		}
		
		if (MASTER == rank) {
			char fname_dbg[100];

			/*
			#ifdef DEBUG_DUMPALL
			snprintf(fname_dbg, 100, "clustering.output/%d/clustering.%d.enditeration.debug", iteration, iteration);

			fd_dbg = fopen(fname_dbg, "w");
			for (i=0; i<2*N; i++) {
				for (j=0; j<2*N; j++) {
					fprintf(fd_dbg, "DM[%d][%d] = %f \t\t SM[%d][%d] = %d \t\t WM[%d][%d] = %d \n", i, j, DM[i][j], i,j, SM[i][j], i,j, WM[i][j]);
				}
				fprintf(fd_dbg, "\n");				
			}
			fclose(fd_dbg);

			snprintf(fname_dbg, 100, "clustering.output/%d/clustering.%d.enditeration.counts.debug", iteration, iteration);

			fd_dbg = fopen(fname_dbg, "w");
			for (i=0; i<2*N; i++) {
				fprintf(fd_dbg, "counts[%d] = %d\n", i, counts[i]);
			}
			fclose(fd_dbg);
			#endif
			*/
			
			// take matrices with 'x' from iteration+1 and write to disk
			// counts should be: countsC = countsA + countsB
			
			//snprintf(fname_dbg, 100, "clustering.output/%d/clustered.%d.matrix", iteration, iteration);
			//write_Freq_Matrices(fname_dbg, matrices, CM[iteration+1], counts, N + N_joins, false);

			//snprintf(fname_dbg, 100, "clustering.output/%d/clustered.%d.logodds.matrix", iteration, iteration);
			//write_Freq_Matrices(fname_dbg, logodds, CM[iteration+1], counts, N + N_joins, true);

			snprintf(fname_dbg, 100, "clustering.output/%d/information_distribution.%d.tab", iteration, iteration);
			write_Information_Distribution(fname_dbg, composition, matrices, CM[iteration+1], N + N_joins);

			write_Iteration_Clusters(iteration_clusters_filename, CM, DITV, N, iteration);
		}
		
		iteration++;
	}
	#ifndef DNEBUG
	if (rank == MASTER) {
		fprintf(stderr, "End of clustering. Finalizing the program.\n");
		write_Freq_Matrices("clustering.output/clustered.all.matrix", matrices, (char *)NULL, counts, N + N_joins, false);
	}
	#endif
		
	eppfreef(DM, 2*N);
	eppfree(SM, 2*N);
	eppfree(WM, 2*N);
	eppfree(CM, N);
	free(DITV);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return(0);
}
