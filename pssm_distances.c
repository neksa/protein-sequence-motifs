/**
 * 	Calculate PSSM-PSSM matrix distances with gapless alignment shifts
 *	Author: Alexandr Goncearenco <ago064@uib.no>
 *	Feb 2009
 *	Affiliation: Berezovsky Group @ CBU, BCCS and University of Bergen
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../sqlite/sqlite3.h"

#define MASTER 0

typedef double Matrix [50][20];


int main (int argc, char *argv[]) {
	// MPI //////////////////////
	MPI_Datatype pssm_matrix_type;
	int blockcounts[2];
	MPI_Aint offsets[2], extent;
	MPI_Status status;
	MPI_Request req;
	int taskid = 0;
	int numtasks = 0;
	////////////////////////////
	// SQLITE //////////////////
	sqlite3 *db;
	int rc = 0;
	sqlite3_stmt *stmt, *select_matrix, *insert_score, *select_pattern;
	char *query;
	char *eee;
	char **results;
	int nrow, ncol;
	////////////////////////////
	int i, j, k, l, m, n, c;
	int matrix_counter = 0;
	int pattern_count = 0;
	int pattern_id = 0;
	int query_pattern = 0;
	double sum_over_aa_pairwise, sum_over_positions, d, distance;	
	Matrix matrix;
	Matrix *matrices;
	int *pattern_ids;
	
	FILE *dist_file;
	

	char amino_acids[20] = "ACDEFGHIKLMNPQRSTVWY";
	//		   A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
	int aa_to_index[26] = {0,-1, 1, 2, 3, 4, 5, 6, 7,-1, 8, 9,10,11,-1,12,13,14,15,16,-1,17,18,-1,19,-1};		
	
	// Check sqlite versions
	if (sqlite3_libversion_number() != SQLITE_VERSION_NUMBER) {
		fprintf(stderr, "Runtime and compiletime sqlite versions differ: %d and %d\n", sqlite3_libversion_number(), SQLITE_VERSION_NUMBER);
		exit(1);
	}
	if (SQLITE_VERSION_NUMBER < 3005009) {
		fprintf(stderr, "Upgrade sqlite to version 3.5.9 or greater\n");
		exit(1);	
	}
	
	// Obtain number of tasks and task ID
	rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		fprintf(stderr, "Error starting MPI program \n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		exit(1);
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	if (numtasks != 40) {
		fprintf(stderr, "Start this program with 40 tasks (40 CPUs), which is equal to the number of shifted matrix alignments\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		exit(1);	    
	}
	fprintf(stderr, "MPI task %d has started (total %d tasks)...\n", taskid, numtasks);
	
	// setup MPI type
	MPI_Type_contiguous(50*20, MPI_DOUBLE, &pssm_matrix_type);
	MPI_Type_commit(&pssm_matrix_type);

	if (taskid == MASTER) {
		rc = sqlite3_open_v2("merged.db", &db, SQLITE_OPEN_READONLY, 0);
		assert(rc == SQLITE_OK);

		sqlite3_extended_result_codes(db, 1);

		query = sqlite3_mprintf("select count(distinct id) from pattern where converged=1 and selected=1");
		sqlite3_get_table(db, query, &results, &nrow, &ncol, &eee);
		sqlite3_free(query);
		
		if (nrow > 0) {
			if (results && results[1]) {
				pattern_count = (int)strtol(results[1], NULL, 10);
			}
		}
		sqlite3_free_table(results);
		
		pattern_ids = (int *)calloc(pattern_count, sizeof(int));
	}

	rc = MPI_Bcast((void *)&pattern_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
	matrices = (Matrix *)calloc(pattern_count, sizeof(Matrix));

	MPI_Barrier(MPI_COMM_WORLD);
	fprintf(stderr, "[%d] Memory initialized for %d matrices.\n", taskid, pattern_count);

	if (taskid == MASTER) {
		// load matrices
		rc = sqlite3_prepare_v2(db, "SELECT DISTINCT id FROM pattern WHERE converged=1 AND selected>0 ORDER BY id", -1, &select_pattern, 0);
		assert(rc == SQLITE_OK);
	
	/*
	// to avoid border effect of unprocessed patterns in the last slice:
	if (taskid < numtasks - 1) {
	    sqlite3_bind_int(select_pattern, 1, slice);
	} else {
	    sqlite3_bind_int(select_pattern, 1, slice + pattern_count % (numtasks - 1));
	}
	sqlite3_bind_int(select_pattern, 2, (taskid - 1)*slice);
	*/
		
		rc = sqlite3_prepare_v2(db, 
		    "SELECT A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V ,W ,Y, id, position FROM pssm WHERE id=?", -1, &select_matrix, 0);
		assert(rc == SQLITE_OK);

		matrix_counter = 0;
		while (sqlite3_step(select_pattern) == SQLITE_ROW) {
			pattern_id = sqlite3_column_int(select_pattern, 0);
			pattern_ids[matrix_counter] = pattern_id;
			sqlite3_bind_int(select_matrix, 1, pattern_id);

			while (sqlite3_step(select_matrix) == SQLITE_ROW) {
				pattern_id = sqlite3_column_int(select_matrix, 20);
				i = sqlite3_column_int(select_matrix, 21) - 1; //position
				for (j = 0; j<20; j++) {
					matrices[matrix_counter][i][j] = sqlite3_column_double(select_matrix, j);
					//printf("M[%d][%d][%d] = %f\n", matrix_counter, i, j, matrices[matrix_counter][i][j]);
				}
			}
			sqlite3_reset(select_matrix);
			matrix_counter++;
		}
    		sqlite3_close(db);
		dist_file = fopen("selected_pssm_distances.csv", "w");
	}
	for (n = 0; n<pattern_count; n++) {
		rc = MPI_Bcast((void *)&matrices[n], 1, pssm_matrix_type, 0, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	d = 0.0; //local distance
	distance = 0.0; //global distance
	for (k = 0; k<pattern_count; k++) { // matrix 1
		for (l = 0; l<=k; l++) { // matrix 2
			MPI_Barrier(MPI_COMM_WORLD);
			if (l == k) {
				d = 0.0;
			} else {
				c = 0;
				while (1) { // just a simple block to use break. Do it one time.
					//printf("0 target = %d, current = %d\n", taskid, c);
					for (i =  0, j = 19; j >= 10; j--, c++) if (c == taskid) break;
					if (c == taskid) break;
					for (j = 10, i =  1; i <  10; i++, c++) if (c == taskid) break;
					if (c == taskid) break; 
					for (i = 10, j =  10; j >=  0; j--, c++) if (c == taskid) break;
					if (c == taskid) break;
					for (j =  0, i = 11; i <  20; i++, c++) if (c == taskid) break;
					//printf("4 target = %d, current = %d\n", taskid, c);
					break;
				}
				assert(c == taskid); // we should be now at our taskid with proper i and j selected
				// i is starting position in matrix 1 (k)
				// j is starting position in matrix 2 (l)
				sum_over_positions = 0.0;
				for (m = 0; m<30; m++) { //position offset
					sum_over_aa_pairwise = 0.0; // sum of squares of pairwise distances for each aminoacid on current position _m_
					for (n = 0; n<20; n++) {// amino acid
						sum_over_aa_pairwise += pow(matrices[k][i+m][n] - matrices[l][j+m][n], 2);
					}
					sum_over_positions += sqrt(sum_over_aa_pairwise);
				}
				d = sum_over_positions/30;
			}
			//printf("%d %d %d %f\n", taskid, k, l, d);
			
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Reduce(&d, &distance, 1, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
			
			if (taskid == MASTER) {
				printf("MASTER %d %d %f\n", k, l, distance);
				fprintf(dist_file, "%d\t%d\t%f\n", pattern_ids[k], pattern_ids[l], distance);
			}
		}
	}
	
	if (taskid == MASTER)
		fclose(dist_file);

	MPI_Finalize();
	return 0;
}
