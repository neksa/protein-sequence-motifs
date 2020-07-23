/**
 *    Compare sequence profiles based on their matches
 *    Coded in C programming language (C99 standard) and MPI
 *    Copyright (c) 2010 Alexandr Goncearenco
 *    
 *    Affiliations: 
 *	Department of Informatics, University of Bergen
 *	CBU, BCCS, Uni Research Berezovsky Group.
 */

/*
INPUT:
    search_matches.tab: matches produced by search PSSM procedure, should contain GI and positions of matches

OUTPUT:
    distance half-matrix for all profile-profile pairs

///////////////////////////

    writer node:
	reads profile matches into an array of records
	broadcast number of profiles
	distributes this array across compute nodes
	compute distance half-matrix using compute nodes (mod distributed)
	write distance half-matrix
    compute node:
	create array for profile matches
	read distributed values for profile matches for each profile
	compute a part of distance half-matrix
	broadcast - assemble distance half-matrix			

 DONE
///////////////////////////
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

//#include <sys/types.h>
//#include <regex.h>

#include "PSSM.h"
//#include "overlaps.h"
 
#include "mpi.h"
#define MASTER 0
#define DATA_TAG 1
#define BREAK_TAG 2
#define OVERLAP_TAG 3
#define PROFILE_TAG 4
#define PROFILE_MATCH_TAG 5

const char *matches_filename = "search_matches.tab";
const char *distances_filename = "profile-profile-overlaps.tab";

typedef struct ProfileMatch {
    int gi;
    int pos;
    struct ProfileMatch *next;
} ProfileMatch;

typedef struct Profile {
    int id;
    int match_counter;
    struct ProfileMatch *matches;
    struct Profile *next;
} Profile;


void free_profile_matches(ProfileMatch *m) {
    if (m == NULL) {
	return;
    }
    ProfileMatch *next = m->next;
    free(m);
    free_profile_matches(next);
}

void free_profile(Profile *p) {
    if (p == NULL) {
	return;
    }
    free_profile_matches(p->matches);
    free(p);
    p = NULL;
}

void free_profiles(Profile *p) {
    if (p == NULL) {
	return;
    }
    Profile *next = p->next;
    free_profile(p);
    free_profiles(next);
}

////////// GLOBAL VARIABLES ////////

int profile_counter = 0;
int max_profile_id = 1000;
Profile *profiles = NULL;
Profile **profile_by_id = NULL;

int rank = 0;
int size = 0;
///////////////////////////////////

void distribute_profile_matches(int from, int to)
{
	int i=0;
	Profile *profile = NULL;
	ProfileMatch *match = NULL;
    
	int c = 0;
	MPI_Status status;
	int message[3];

	//MPI_Bcast(number_of_profiles, MPI_INT, MASTER, MPI_COMM_WORLD);
	if (rank == MASTER) {
		for (i=from;i<=to;i++) {
			profile = profile_by_id[i];
			if (profile == NULL) continue;
			message[0] = i; // id
			message[1] = profile->match_counter;
			
			for (c=1; c<size; c++) MPI_Send(&message, 2, MPI_INT, c, PROFILE_TAG, MPI_COMM_WORLD);
		}
		
		for (c=1; c<size; c++) MPI_Send(0, 0, MPI_INT, c, BREAK_TAG, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);	
		
		for (i=from;i<=to;i++) {
			profile = profile_by_id[i];
			if (profile == NULL) continue;
			match = profile->matches;
			while (match != NULL) {
				message[0] = i;
				message[1] = match->gi;
				message[2] = match->pos;
				for (c=1; c<size; c++) MPI_Send(&message, 3, MPI_INT, c, PROFILE_MATCH_TAG, MPI_COMM_WORLD);
				match = match->next;
			}
		}
		
		for (c=1; c<size; c++) MPI_Send(0, 0, MPI_INT, c, BREAK_TAG, MPI_COMM_WORLD);
	
	} else {
		while (1) {
			MPI_Recv(&message, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			
			if (status.MPI_TAG == BREAK_TAG) break; 
			
			if (status.MPI_TAG == PROFILE_TAG) {
				profile = profile_by_id[message[0]];

				if (profile != NULL) {
					free_profile_matches(profile->matches);
				} else {
					profile = malloc(sizeof(Profile));
					profile->id = message[0];
					if (profile == NULL) {
						printf("memory allocation error for new profile\n");
					}
					profile->next = profiles;
					profiles = profile;
				}
				profile->matches = NULL;
				profile->match_counter = 0; // we'll increase it as matches arrive from MASTER
				profile_by_id[message[0]] = profile;
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);	
		
		while (1) {
			MPI_Recv(&message, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			
			if (status.MPI_TAG == BREAK_TAG) break; 
			
			if (status.MPI_TAG == PROFILE_MATCH_TAG) {
				profile = profile_by_id[message[0]];
				if (profile == NULL) continue; // actually, this should not happen

				match = NULL;
				match = malloc(sizeof(ProfileMatch));
				if (match == NULL) {
					printf("memory allocation error for new profile match\n");
				}
				match->gi = message[1];
				match->pos = message[2];
				match->next = NULL;
			
				if (profile->matches == NULL) {
					profile->matches = match;
				} else {
					// add new match to the head of the list
					match->next = profile->matches;
					profile->matches = match;
				}
				profile->match_counter++;
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);	
}

void load_search_matches_master(const char *search_matches_filename)
{
	int i = 0;

	FILE *fd = NULL;

	Profile *profile = NULL;
	ProfileMatch *match = NULL;	
	
	char *buf = calloc(256, sizeof(char));
	int profile_id, position, gi;
	double evalue, score;
	fd = fopen(search_matches_filename, "r");
	assert(fd != NULL);
	//4723    1.42    0.0     218     gi|14602183|ref|NP_146957.1|scop|c.42.1 hypothetical protein APE_0092.1
	while (fgets(buf, 256, fd)) {
		if (5 == sscanf(buf, "%d\t%lf\t%lf\t%d\tgi|%d", &profile_id, &score, &evalue, &position, &gi)) {
			//printf("scanned: %d\t%.2f\t%d\t%d\t%d\n", profile_id, score, (int)round(evalue), position, gi);
			
			if (profile_id < 0) {
			    printf("Error: profile_id is negative! skipping the line\n");
			    continue;
			}
			
			if (gi == 0) {
			    printf("Error: GI is zero! skipping the line\n");
			    continue;
			}
			
			/*
			// find if exists
			profile = profiles;
			while (profile != NULL) {
			    if (profile->id == profile_id) {
				//printf("Found profile id=%d\n", profile_id);
				break;
			    }
			    profile = profile->next;
			}*/
			profile = NULL;
			if (profile_id <= max_profile_id) {
			    profile = profile_by_id[profile_id]; // will be NULL if not found
			}
			
			// not found or empty
			if (profile == NULL) {
    				profile = malloc(sizeof(Profile));
    				profile_counter++;
				if (profile == NULL) {
				    printf("memory allocation error for new profile\n");
				}
				profile->id = profile_id;
				profile->match_counter = 0;
				profile->next = profiles;
				profile->matches = NULL;
				
				profiles = profile;
				
				if (profile_id > max_profile_id) {
				    i = max_profile_id;
				    max_profile_id = profile_id + 1000; // preallocation in chunks
				    profile_by_id = (Profile **)realloc(profile_by_id, (max_profile_id + 1) * sizeof(Profile *));
				    for (;i<max_profile_id;i++) {
					profile_by_id[i] = NULL;
				    }
				}
				profile_by_id[profile_id] = profile;
			}
			
			match = NULL;
			match = malloc(sizeof(ProfileMatch));
			if (match == NULL) {
			    printf("memory allocation error for new profile match\n");
			}
			match->gi = gi;
			match->pos = position;
			match->next = NULL;
			
			if (profile->matches == NULL) {
			    profile->matches = match;
			} else {
			    // add new match to the head of the list
			    match->next = profile->matches;
			    profile->matches = match;
			}
			profile->match_counter++;
			
		}
		else
		{
			printf("scan error. skipping line\n");
		}
	}
	free(buf);
	fclose(fd);

	printf("[MASTER] Profiles loaded: %d\n", profile_counter);
}

void calculate_pairwise_distance_overlap_master()
{	
	int message[10];
	int p1;
	int p2;
	int p1_card;
	int p2_card;
	int p1_inter_card;
	int p2_inter_card;
	double offset_mean;
	double offset_stddev;
	double offset_consensus_frequency;
	int offset_consensus;
	
	int break_counter = 0;
	double distance = 0.0;

	MPI_Status status;
	FILE *fd = NULL;

	fd = fopen(distances_filename, "w");
	assert(fd != NULL);
	fprintf(fd, "p1\tp2\tdistance\tp1_card\tp2_card\tp1_inter_card\tp2_inter_card\toffse_consensus\toffset_consensus_freq\toffset_mean\toffset_stddev\n");
	
	while (1) {
		MPI_Recv(&message, 10, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		p1 = message[0];
		p2 = message[1];
		p1_card = message[2];
		p2_card = message[3];
		p1_inter_card = message[4];
		p2_inter_card = message[5];
		offset_mean = message[6]/100000.0;
		offset_stddev = message[7]/100000.0;
		offset_consensus = message[8];
		offset_consensus_frequency = message[9]/100000.0;
	
		if (status.MPI_TAG == BREAK_TAG) {
			if (++break_counter == size - 1) {
				break;
			}
		}
		if (status.MPI_TAG == OVERLAP_TAG) {
		
		    if (p1_card + p2_card == 0) continue; // fool proof division by zero
		    // Jackard distance = 1 - intersect/union;
		    distance = 1.0 - (p1_inter_card + p2_inter_card) / (1.0 * (p1_card + p2_card));
		    #ifndef NDEBUG
		    //printf("[pair]: %d\t%d\t%d\t%d\t%d\t%d\t%.4lf\n",
		    //    *p1, *p2, *p1_card, *p2_card, *p1_inter_card, *p2_inter_card, distance);
		    #endif
		
		    fprintf(fd, "%d\t%d\t%.4lf\t%d\t%d\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\n",
		        p1, p2, distance, p1_card, p2_card, p1_inter_card, p2_inter_card,
		    	offset_consensus, offset_consensus_frequency, offset_mean, offset_stddev);
		}
	}
	
	fclose(fd);
}


void calculate_pairwise_distance_overlap()
{
	if (rank == MASTER) {
		calculate_pairwise_distance_overlap_master();
		return;
	}
		
	int message[10];

	int pp1, pp2;
	ProfileMatch *match1, *match2;
	int m1_index = 0;
	int m2_index = 0;
	bool m1_flag = false;
	int m1_counter = 0;
	int m2_counter = 0;
	
	int overlap_counter;
	char *offsets;
	int sum_of_offsets;
	int offset_consensus = 0;
	double offset_consensus_frequency;
	double offset_mean;
	double sum_of_squares;
	double offset_stddev;
	int offset_distribution[100]; // center it at 50
	
	unsigned char m2_bitmap[8129]; // 8 bit * 8129 bytes = 65536 bits or matches (should be as much as +inf for real-life cases)
	int i = 0;

	offsets = malloc(10000000 * sizeof(char)); // max 10M overlaps, should be sufficient
	
	for (pp1 = 0; pp1 < max_profile_id; pp1++) {
		// parallelization: comb split
		if (pp1 % (size-1) != (rank - 1))
		    continue;
	
		if (profile_by_id[pp1] == NULL)
		    continue;

		for (pp2 = 0; pp2 < pp1; pp2++) {
		    if (profile_by_id[pp2] == NULL)
			continue;
			
		    overlap_counter = 0;
		    
		    m1_index = 0;
		    m2_index = 0;
		    m1_counter = 0;
		    memset(&m2_bitmap, 0, sizeof(m2_bitmap));
		    //printf("comparing %d and %d with %d and %d matches\n", pp1, pp2, profile_by_id[pp1]->match_counter, profile_by_id[pp2]->match_counter);
		    match1 = profile_by_id[pp1]->matches;
		    while (match1 != NULL) {
			match2 = profile_by_id[pp2]->matches;
			m2_index = 0;
			while (match2 != NULL) {
			    if (match1->gi == match2->gi && (abs(match2->pos - match1->pos) < 20)) {
				    assert(overlap_counter < 10000000);
				    offsets[overlap_counter++] = match2->pos - match1->pos;

				    m1_flag = true;
				    m2_bitmap[m2_index / (8 * sizeof(unsigned char))] |= (1 << (m2_index % (8 * sizeof(unsigned char))));
			    }
			    match2 = match2->next;
			    m2_index++;
			}
			if (m1_flag) {
			    m1_counter++;
			    m1_flag = false;
			}
			match1 = match1->next;
			m1_index++;
		    }
		    m2_counter = 0;
		    //if (profile_by_id[pp2]->match_counter >= 8129 * 8 * sizeof(unsigned char)) {
		    //	printf("Error: profile %d match_counter = %d\n", pp2, profile_by_id[pp2]->match_counter);
		    //}
		    assert(profile_by_id[pp2]->match_counter < 8129 * 8 * sizeof(unsigned char));
		    for (m2_index=0; m2_index < profile_by_id[pp2]->match_counter; m2_index++) {
			m2_counter += (m2_bitmap[m2_index / (8 * sizeof(unsigned char))]
					& (1 << (m2_index % (8 * sizeof(unsigned char)))))
					    >> (m2_index % (8 * sizeof(unsigned char)));
		    }
		    
		    // offset consensus
		    for (i=0; i<100; i++) {
			    offset_distribution[i] = 0;
		    }
		    for (i=0; i<overlap_counter;i++) {
			    offset_distribution[offsets[i] + 50]++; //+50 for non-negative index
		    }
		    int max_offset_counter = 0;
		    offset_consensus = 0;
		    for (i=0; i<100; i++) {
			    if (max_offset_counter <= offset_distribution[i]) {
				    max_offset_counter = offset_distribution[i];
				    offset_consensus = i-50;
			    }
		    }
		    offset_consensus_frequency = (1.0 * max_offset_counter) / overlap_counter;
		    
		    // offset mean and std deviation
		    sum_of_offsets = 0;
		    offset_mean = 0;
		    for (i=0;i<overlap_counter;i++) {
			    sum_of_offsets += offsets[i];
		    }
		    if (overlap_counter > 1) {
			    offset_mean = (1.0 * sum_of_offsets) / overlap_counter;
			    sum_of_squares = 0.0;
			    for (i=0;i<overlap_counter;i++) {
				    sum_of_squares += (offsets[i] - offset_mean) * (offsets[i] - offset_mean);
			    }
			    offset_stddev = sqrt((1.0/(overlap_counter - 1.0)) * sum_of_squares);
		    } else {
			    offset_stddev = 0.0;
		    }
		    
		    // if intersect not empty
		    if (m2_counter + m1_counter > 0) {
			    /*
			    printf("overlaps: %d (%d) and %d (%d): Sum=%d Union=%d Intersect=%d (%d+%d)\n",
				pp1, profile_by_id[pp1]->match_counter, pp2, profile_by_id[pp2]->match_counter,
				(profile_by_id[pp1]->match_counter + profile_by_id[pp2]->match_counter),
				(profile_by_id[pp1]->match_counter + profile_by_id[pp2]->match_counter - m1_counter - m2_counter),
				m1_counter + m2_counter, m1_counter, m2_counter);
			    */
			
			    message[0] = pp1;
			    message[1] = pp2;
			    message[2] = profile_by_id[pp1]->match_counter;
			    message[3] = profile_by_id[pp2]->match_counter;
			    message[4] = m1_counter;
			    message[5] = m2_counter;
			    message[6] = (int)round(offset_mean*100000);
			    message[7] = (int)round(offset_stddev*100000);
			    message[8] = offset_consensus;
			    message[9] = (int)round(offset_consensus_frequency*100000);

			    MPI_Send(&message, 10, MPI_INT, MASTER, OVERLAP_TAG, MPI_COMM_WORLD);
		    }
		}
	}
	free(offsets);
	
	MPI_Send(0, 0, MPI_INT, MASTER, BREAK_TAG, MPI_COMM_WORLD);
}

int main (int argc, char *argv[])
{
	int i = 0;
	int rc = 0;
	
	rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		fprintf(stderr, "Error while initializing MPI\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		exit(1);
	}	
		
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	if (size < 2) {
		printf("More than one CPU is required\nExiting\n");
		MPI_Finalize();
		return(1);
	}
	assert(size > 1);
	
	#ifndef DNEBUG
	fprintf(stderr, "MPI rank %d started. Total %d CPU\n", rank, size);
	#endif

	if (rank == MASTER) {
		profile_by_id = calloc(max_profile_id + 1, sizeof(Profile *));
		for (i=0;i<max_profile_id;i++) {
			profile_by_id[i] = NULL;
		}
	        load_search_matches_master(matches_filename);
	}

	MPI_Bcast(&profile_counter, 1, MPI_INT, MASTER, MPI_COMM_WORLD);	
    	MPI_Bcast(&max_profile_id, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

	if (rank != MASTER) {
		profile_by_id = calloc(max_profile_id + 1, sizeof(Profile *));
		for (i=0;i<max_profile_id;i++) {
		    profile_by_id[i] = NULL;
		}
	}
	distribute_profile_matches(0, max_profile_id);
	printf("Profile matches are distributed\n");

	//debug: printf("[rank %d] Profile %d has %d matches\n", rank, profile_by_id[40336]->id, profile_by_id[40336]->match_counter);
	
	calculate_pairwise_distance_overlap();

	free(profile_by_id);
	free_profiles(profiles);

	printf("[rank %d] DONE.\n", rank);		

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return(0);
}

