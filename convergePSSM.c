/**
 * Converge initial_segment procedure
 *    Parallel implementation for use on HPC
 *    Coded in C programming language (C99 standard) and MPI
 *    Copyright (c) 2009 Alexandr Goncearenco
 *    Affiliation: CBU, BCCS, UNIFOB AS, UiB, Berezovsky Group.
 */

/*
INPUT:
    p:  global p-value
    proteome: protein sequences in fasta format
    

///////////////////////////

    +initial parameters: p=0.001, mu=0.05, eta=0.01, e=1, delta=5
    +read: one initial_segment from fasta into memory
    +read: proteome fasta sequences into memory
    +calculate composition from proteome sequences
    +make F count matrix from initial_segment, K=1
    +M = calcPSSM(F, K, composition)
    +iteration=0
    +write: M, p, S(p)=NULL, K=0, parent=NULL, iteration
    
 loop:
    read: max iteration
    read: M where iteration==max iteration
    +runP(M, permutation_vector, p-value) return S(p-value)
    +runS(M, S(p-value), &callback_match_found()) return: { count matrix F, number of matches K }
    +M' = calcPSSM(F, K, composition)
    +d = distance(M, M')
    +write: M', p, S(p), K, link to M, d, iteration
    +d < e => M' converged => DONE.
    +goto loop
    
 DONE
///////////////////////////
 
 then repeat procedure for another initial_segment
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <errno.h>
//#include <error.h>

#include  "PSSM.h"
#include  "convergePSSM.h"

#define MASTER 0

#define MAXLEN 250

// parameters
//
// with N = 20 M: p= 5^10-6 => E=10. rho = 10. epsilon=3
//
const int delta = 10;		// cutting step for initial segments, in residues
const double omega = 1.0;	// pseudocount scaling factor in PSSM
//const double p = 0.0000005;	// initial p-value, lower threshold
const double rho = 1;		// threshold ratio shuffled/natural number of matches where score > score(p)
const double mu = 0.05;		// upper threshold for p-value
const double eta = 0.01;	// increment of p-value in clustering
const double epsilon = 1;	// max distance(M1, M2) to consider M2 as converged
const int theta = 30;		// max # iterations

const char *initial_segments_filename = "initial.fasta";
const char *proteome_filename = "proteome.fasta";

int main (int argc, char *argv[]) {
	//Matrix *matrix = NULL;
	///////////////////////////////////
	
	int a;
	char initial_segments_filename[MAXLEN+1] = "initial.fasta";
	char proteome_filename[MAXLEN+1] = "proteome.fasta";
	char output_prefix[MAXLEN+1] = "output/output.matrix.%d";
	char composition_filename[MAXLEN+1] = "composition.csv";
	int number_of_randomizations = 10;
	double E_value = 1;
	int K_min = 1;
	int output_distributions = 0;
	int input_format = 0; // default: 0 = FASTA, 3 or 4 = Matrix
	bool origin_is_matrix = false;
	bool write_composition_flag = false;
	bool use_blosum = false;
	int BLOSUM[20][20];
	int length = 30; // desired resulting profile length
			
	// process command line parameters
	for (a=1; a<argc-1; a++) {
		//printf("OPTION %d\n", a);
		if (argv[a] && !strcmp(argv[a], "-i")) {
			if (argv[a+1] && strlen(argv[a+1]) > 0 && strlen(argv[a+1]) < MAXLEN) {
			    strncpy(initial_segments_filename, argv[a+1], MAXLEN);
			}
		}
		if (argv[a] && !strcmp(argv[a], "-p")) {
			if (argv[a+1] && strlen(argv[a+1]) > 0 && strlen(argv[a+1]) < MAXLEN) {
			    strncpy(proteome_filename, argv[a+1], MAXLEN);
			}
		}
		if (argv[a] && !strcmp(argv[a], "-o")) {
			if (argv[a+1] && strlen(argv[a+1]) > 0 && strlen(argv[a+1]) < MAXLEN) {
			    strncpy(output_prefix, argv[a+1], MAXLEN);
			}
		}
		if (argv[a] && !strcmp(argv[a], "-c")) {
			if (argv[a+1] && strlen(argv[a+1]) > 0 && strlen(argv[a+1]) < MAXLEN) {
			    strncpy(composition_filename, argv[a+1], MAXLEN);
			    write_composition_flag = true;
			}
		}
		if (argv[a] && !strcmp(argv[a], "-r")) {
			if (argv[a+1] && strlen(argv[a+1])) {
			    number_of_randomizations = atoi(argv[a+1]);
			    printf("OPTION: Randomizations = %d\n", number_of_randomizations);
			}
		}
		if (argv[a] && !strcmp(argv[a], "-E")) {
			if (argv[a+1] && strlen(argv[a+1])) {
			    E_value = atof(argv[a+1]);
			    printf("OPTION: E-value = %f\n", E_value);
			}
		}
		if (argv[a] && !strcmp(argv[a], "-K")) {
			if (argv[a+1] && strlen(argv[a+1])) {
			    K_min = atoi(argv[a+1]);
			    printf("OPTION: min K filtering = %d\n", K_min);
			}
		}
		if (argv[a] && !strcmp(argv[a], "-L")) {
			if (argv[a+1] && strlen(argv[a+1])) {
			    length = atoi(argv[a+1]);
			    printf("OPTION: desired profile length L = %d\n", length);
			    
			    if (length != 29 && length != 30 && length !=50) {
				printf("WARNING: incorrect profile length specified (only 29, 30 and 50 allowed)!\n");
			    }
			}
		}
		if (argv[a] && !strcmp(argv[a], "-f")) {
			if (argv[a+1] && strlen(argv[a+1])) {
			    input_format = atoi(argv[a+1]);
			    printf("OPTION: Input format = Matrix (%d)\n", input_format);
			    if (input_format != 3 && input_format !=4) {
				printf("WARNING: incorrect input format specified!\n");
			    } else {
				origin_is_matrix = true;
			    }
			}
		}
		if (argv[a] && !strcmp(argv[a], "-d")) {
			output_distributions = 1;
			printf("OPTION: output distributions = TRUE\n");
		}
		if (argv[a] && !strcmp(argv[a], "-B")) {
			use_blosum = true;
			printf("OPTION: use BLOSUM62 substitution table = TRUE\n");
		}
	}

	//////////////////////////////////
	
	int i=0, j=0, k=0, l=0, m=0, q=0, f=0;
	double aa_count[26]; for (j=26;j--;) aa_count[j] = 0.0;
	double composition[20]; for (j=20;j--;) composition[j] = 0.0;
	double F[MAX_PSSM_LENGTH][26]; 	for (i=MAX_PSSM_LENGTH;i--;) for (j=26;j--;) F[i][j] = 0.0;
	double PSSM[MAX_PSSM_LENGTH][26];	for (i=MAX_PSSM_LENGTH;i--;) for (j=26;j--;) PSSM[i][j] = 0.0;
	double M[MAX_PSSM_LENGTH][26];	for (i=MAX_PSSM_LENGTH;i--;) for (j=26;j--;) M[i][j] = 0.0;
	double M1[MAX_PSSM_LENGTH][26];	for (i=MAX_PSSM_LENGTH;i--;) for (j=26;j--;) M1[i][j] = 0.0;
	double Information[MAX_PSSM_LENGTH];	for (i=MAX_PSSM_LENGTH;i--;) Information[i] = 0.0;

	double S = 0.0;
		
	char ch = '\0';
	char *seq = NULL;
	char initial_segment[MAX_PSSM_LENGTH + 1] = "";
	char origin[61] = "";
	char *protein = NULL;
	
	Matrix *converged_matrices = NULL;
	int converged_matrices_count = 0;
	
	Sequence *initial_segments = NULL, **pinitial_segments = &initial_segments;
	int N_initial = 0;  	  // # initial sequences
	//int N_initial_slice = 0;  // # initial sequences
	Sequence *proteome = NULL, **pproteome = &proteome;
	int N_proteome = 0; // # proteome sequences
	int N = 0; 	    // # fragments in all proteomic sequences

	if (use_blosum) {
	    load_BLOSUM("BLOSUM62", BLOSUM);
	    printf ("BLOSUM OK");
	}
	
	N_proteome = load_fasta_sequences(proteome_filename, pproteome, false);
	
	unsigned long size_of_proteome = 0;
	unsigned long size_of_proteome_aminoacid = 0;
	for (k=0;k<N_proteome;k++) {
		for (i=0;(ch = (proteome + k)->sequence[i++]);) {
			aa_count[ch-'A']++;
			size_of_proteome++;
		}
		for (i=0; i+length <= strlen((proteome + k)->sequence); i++) {
			N++;
		}
	}
	
	for (j=20; j--;) {
		q = amino_acids[j]-'A';
		size_of_proteome_aminoacid += aa_count[q];
	}


	for (j=20; j--;) {
		q = amino_acids[j]-'A';
		composition[j] = aa_count[q] / size_of_proteome_aminoacid;
		#ifndef DNEBUG
		printf("%c %lf\n", amino_acids[j], composition[j]);
		#endif
	}

	if (write_composition_flag) {
		write_Composition(composition_filename, composition);
	}

	Matrix *matrix = NULL, *matrices = NULL, **pmatrices = &matrices;

	switch (input_format) {
	    case 3:
		N_initial = load_VariableMatricesCount(initial_segments_filename, pmatrices, composition);
		break;
	    case 4:
		N_initial = load_VariableMatricesFreq(initial_segments_filename, pmatrices, composition);
		break;
	    default:
		N_initial = load_fasta_sequences(initial_segments_filename, pinitial_segments, false);
		break;
	}
	
	double FP = E_value;
	double p = FP/N;
	
	List *list = NULL, *item = NULL, *curr = NULL, *prev = NULL;
	int list_length = 0;
	double maxS = 0.0;
	int from = 0;
	int to = N_initial;
	
	int last_informative_position = 0;
	
	//N_initial
	for (k = from; k < to; k++) {
		if (origin_is_matrix) {
			matrix = matrices + k;
			
			length = matrix->L;
			origin[60] = '\0';
			strncpy(origin, matrix->origin, strlen(matrix->origin));
			origin[strlen(matrix->origin)] = '\0';
		
			#ifndef NDEBUG
			printf("Input MATRIX %d (%s) K=%d L=%d\n", matrix->id, origin, matrix->K, matrix->L);
			#endif
		} else {
			seq = (initial_segments + k)->sequence;
			
			#ifndef NDEBUG
			printf("Sequence %d len=%d: %s\n", k, (int)strlen(seq), seq);
			#endif
		}
		
		for (l=0; origin_is_matrix || (l + length <= strlen(seq)); l+=delta) {

			if (!origin_is_matrix) {
				initial_segment[length] = '\0';
				strncpy(initial_segment, seq + l, length);
			
				bool flag1 = false;
    				for (i=length;i--;) {
					q = initial_segment[i];
					if (90==q || 74==q || 79==q || 85==q || 88==q || 66==q) {
						flag1 = true;
						break;
					}
				}
				if (flag1) continue;
				
				// new initial segment obtained here!
				for (i=length;i--;) for (j=26;j--;) {
					F[i][j] = 0;
					M[i][j] = 0;
					PSSM[i][j] = 0;
				}
			
				if (use_blosum) {
					int match = 0;
					for (i=0; i<length;i++) {
						q = initial_segment[i] - 'A';
					    	F[i][q] = 1; // initial number of amino acids
					    	M[i][q] = 1;
					    	for(j=20; j--;) {
					    		if (initial_segment[i] == amino_acids[j]) {
					    			match = j;
					    		}
					    	}
				    		//printf("F[%d][%c] = %d\t%f\n", i, q+'A', F[i][q], composition[q]);
						for (j=20; j--;) {
							q = amino_acids[j]-'A';
					    		PSSM[i][q] = BLOSUM[match][j];
				    		}
					}
				}
				else {
					double composition_sum_squares = 0.0;
	
					for (j=20;j--;) {
						composition_sum_squares += composition[j]*composition[j];
					}
					
					for (i=0; i<length;i++) {
						q = initial_segment[i] - 'A';
					    	F[i][q] = 1; // initial number of amino acids
				    		//printf("F[%d][%c] = %d\t%f\n", i, q+'A', F[i][q], composition[q]);
						for (j=20; j--;) {
							q = amino_acids[j]-'A';
					    		//M[i][j] = log( ((F[i][j] + omega)/(1 + 20*omega)) / composition[j] );
					    		//M[i][q] = log(((F[i][q] + 2*omega*composition[q])/(5 + 20*2*omega)) / composition[q]);
				    			
				    			M[i][q] = (F[i][q] + omega*composition[j]*composition[j]) / (1 + omega*composition_sum_squares);
					    		PSSM[i][q] = log( M[i][q] / composition[j] );
				    		}
					}
				}
				//print_PSSM(M, 0);
				#ifndef NDEBUG
				print_PSSM(PSSM, 0);
				#endif
				
			// origin is a matrix
			} else {
				for (i=length;i--;) for (j=26;j--;) {
					PSSM[i][j] = 0.0;
					M[i][j] = matrix->freq[i][j];
					//F[i][j] = (matrix->K * matrix->freq[i][j]);
				}
				for (i=length;i--;) {
					Information[i] = 0.0;
				}
			
				double sum_aa = 0;
				double fragment_freq = 0.0;
				double pc = 0.0;
				double H_sum = 0.0; // attention: it is in nats, not in bits
			
				//print_PSSM(matrix->freq, 0.0);
				
				last_informative_position = 0;
				for (i=0;i<length;i++) {
					H_sum = 0.0;
					for (j=20;j--;) {
						q = amino_acids[j] - 'A';
						if (M[i][q] > 0.0) {
							H_sum += M[i][q] * (log(M[i][q])/log(2));
						}
					}
					Information[i] = log(20)/log(2) + H_sum;
					if (Information[i] < 1.0) {
					    Information[i] = 0.0;
					} else {
					    last_informative_position = i;
					}
					//printf("%d %f %d\n", i, Information[i], last_informative_position);
				}
			
				for (j=20;j--;) {
					q = amino_acids[j] - 'A';
			    	
				    	sum_aa = 0;
				    	for (i=length;i--;) {
					    sum_aa += M[i][q] * matrix->K;
					}
					fragment_freq = sum_aa/(length * 1.0 * matrix->K);
				    	//compensating by observed frequencies in aligned fragment
					pc = 0.001;
					if (fragment_freq < composition[j]) {
						pc = composition[j] - fragment_freq;// /composition[j];
					}
					//printf("PC[%d][%c] = %f\n", i, q+'A', pc);
					for (i=length;i--;) {
						PSSM[i][q] =  Information[i] * 
								    log( ((M[i][q]*matrix->K + omega*pc) / (matrix->K + 20.0*omega*pc)) 
								    / composition[j] );
			    		}
			    		
			    		
			    		
			    		//M1[i][q] = (F[i][q] + omega*composition[j]*composition[j]) / (Kmatches + omega*composition_sum_squares);
					//PSSM[i][q] = log( M1[i][q] / composition[j] );
			    		
				}
				//print_PSSM(PSSM, 0);
			}

			/////////////////////////////// CONVERGENCE //////////////////////////////////

			maxS = 0.0;
			
			// continue processing of one segment here:
			bool converged = false;
			int iteration = 0;
			int randomization = 0;
			double sampled_maxS = -0xFFFF; //-inf
			double H_sum = 0.0;
			FILE *fd_scores = NULL;
			char scores_filename[50];
		
			while (!converged) {
			
			    for (i=length;i--;) {
				Information[i] = 0.0;
			    }
			    
			    H_sum = 0.0; 
			    
			    for (i=length;i--;) {
				H_sum = 0.0;
				for (j=20;j--;) {
					q = amino_acids[j] - 'A';
					if (M[i][q] > 0.0) {
						H_sum += M[i][q] * (log(M[i][q])/log(2));
					}
				}
				Information[i] = log(20)/log(2) + H_sum;
				if (Information[i] < 1.0) {Information[i] = 0.0;}
				//printf("Information[%d] = %f\n", i, Information[i]);			    
			    }
			    
			    // a hack
			    /*
			    for (i=length;i--;) {
				Information[i] *= 0.2;
			    }
			    Information[20] *= 5;
			    Information[21] *= 5;
			    Information[22] *= 5;
			    Information[23] *= 5;
			    */
			
			    maxS = 0.0;
			    sampled_maxS = -0xFFFF;

			    
			    /////////////////////////////// BACKGROUND DISTRIBUTIONS //////////////////////////////////
			
			    for (randomization=0; randomization < number_of_randomizations; randomization++) {
				if (output_distributions) {
				    	snprintf(scores_filename, 50, "converge_scores_shuffled_%d_%d_%d.tab", k, iteration, randomization);
					fd_scores = fopen(scores_filename, "w");
				}
			    
				for (f=0; f<N_proteome; f++) {
					protein = (proteome + f)->sequence;
					int protein_length = strlen(protein);
					for (m=0; m+length <= protein_length; m++) {
						S = 0.0;
						for (i=0;i<length;i++) {
							// compare 50 residues in one fragment from proteome with reshuffled matrix position
							q = protein[m+i];
							if ('X' == q) {
								//S = -0xFFFF;
								//goto next_segment;
								S += 0; // ATTN: X gives a zero score (which is not bad)
							} else {
								if (length == 30) {
								    S += Information[random_index_30[randomization][i]] * PSSM[random_index_30[randomization][i]][q-'A'];
								}
								if (length == 29) {
								    S += Information[random_index_29[randomization][i]] * PSSM[random_index_29[randomization][i]][q-'A'];
								}
								if (length == 50) {
								    S += Information[random_index_50[randomization][i]] * PSSM[random_index_50[randomization][i]][q-'A'];
								}
							}
						}
						S /= length;
					
						if (output_distributions) {
							fprintf(fd_scores, "%f\n", S);
						}
						
						if (list_length < FP || S > maxS) {
							//printf("S=%f, maxS=%f, list=%d, FP=%d\n", S, maxS, list_length, (int)FP);
							curr = list;
							prev = NULL;
							while (curr != NULL && curr->score < S) {
								prev = curr;
								curr = curr->next;
							}
							item = malloc(sizeof(List));
							item->score = S;
							item->next = curr;
							if (list == item->next) {
								list = item;
							}
							if (prev != NULL) {
								prev->next = item;
							}						
							list_length++;
													
							if (list_length > FP) {
								//truncate one from list head
								prev = list;
								list = list->next;
								free(prev);
								list_length--;
							}
						}
						maxS = list->score;
						continue;
					}
				}
				while (list != NULL) {
					prev = list;
					list = list->next;
					free(prev);
					list_length--;
				}
				assert(list_length == 0);
				// max score
				#ifndef NDEBUG
				printf("E(s > t) < %f. t=%f\n", FP, maxS);

				//file_name = snprintf("", 50, );
				//fprintf("S(p) = S(%f) = S(%d/%d) = %f\n", p, (int)FP, N, maxS);
				#endif
				
				if (sampled_maxS < maxS) {
					sampled_maxS = maxS;
				}
				
				if (output_distributions) {
					fclose(fd_scores);
				}
			} // randomization


			/////////////////////////////// SCORE DISTRIBUTION //////////////////////////////////

			maxS = sampled_maxS;
			printf("Sampled maxS = %f\n", maxS);

			int Kmatches = 0;
			
			for (i=length;i--;) for (j=26;j--;) {
				F[i][j] = 0.0;
			}
			
			// Use S(p) as threshold on natural matrix matches:
			if (output_distributions) {
				snprintf(scores_filename, 50, "converge_scores_%d_%d.tab", k, iteration);
				fd_scores = fopen(scores_filename, "w");
			}

			#ifndef NDEBUG
			printf("Bits: ");
			for (i=0; i<length; i++) {
				if (Information[i] > 1.5) {
					if (Information[i] > 2) {
						if (Information[i] > 3) {
							printf("*");
						} else {
							printf("+");
						}
					} else {
						printf(".");
					}
				} else {
					printf(" ");
				}
			}
			printf("\n");
			#endif

			for (f=0; f<N_proteome; f++) {
				protein = (proteome + f)->sequence;
				int protein_length = strlen(protein);
				for (m=0; m+length <= protein_length; m++) {
					S = 0.0;
					for (i=0;i<length;i++) {
						q = protein[m+i];
						if ('X' == q) {
							//S = -0xFFFF;
							//goto next_segment2;
							S += 0; // ATTN: X gives a zero score (which is not bad). This code should match the background score section
						} else {
							S += Information[i] * PSSM[i][q-'A'];
						}
					}
					S /= length;
					
					if (output_distributions && S > -0xFFFF) {
						fprintf(fd_scores, "%f\n", S);
					}

					if (S > maxS) {
						Kmatches++;
						#ifndef NDEBUG
						char fragment[MAX_PSSM_LENGTH + 1] = "";
						strncpy(fragment, protein + m, length);
						fragment[length] = '\0';
						
						printf("match[%s]{%2.2f}%-90.90s\n", fragment, S, (proteome + f)->description);
						#endif
						for (i=0;i<length;i++) {
							F[i][protein[m+i]-'A']++;
						}
					}
					continue;
				}
			}
			if (output_distributions) {
				fclose(fd_scores);
			}

			//if (iteration > 0 && 1.0*Kmatches/FP < rho) {
			//	#ifndef NDEBUG
			//	printf("BREAK: low rho=%f\n", 1.0*Kmatches/FP);
			//	#endif
			//	break;
			//}
			
			for (i=length;i--;) for (j=26;j--;) {
				PSSM[i][j] = 0.0;
				M1[i][j] = 0.0;
			}
			//double sum_aa = 0;
			//double fragment_freq = 0.0;
			//double pc = 0.0;
			double composition_sum_squares = 0.0;

			for (j=20;j--;) {
				composition_sum_squares += composition[j]*composition[j];
			}

			for (j=20;j--;) {
				q = amino_acids[j] - 'A';
			    	
			    	//sum_aa = 0;
			    	//for (i=length;i--;) {
				//    sum_aa += F[i][q];
				//}
				//fragment_freq = sum_aa/(length * 1.0 * Kmatches);
			    	//compensating by observed frequencies in aligned fragment
					
				//pc = 0.001;
				//if (fragment_freq < composition[j]) {
				//	pc = composition[j] - fragment_freq;// /composition[j];
				//}
				//printf("PC[%d][%c] = %f\n", i, q+'A', pc);
				
				for (i=length;i--;) {	
					//M1[i][q] = (F[i][q] + omega*pc) / (Kmatches + 20.0*omega*pc);
					M1[i][q] = (F[i][q] + omega*composition[j]*composition[j]) / (Kmatches + omega*composition_sum_squares);
					PSSM[i][q] = log( M1[i][q] / composition[j] );
			    	
			    		//printf("M1[%d][%c] = %f\n", i, q+'A', M1[i][q]);
			    	}
			}
			
			double sum = 0.0;
			double distance = 0.0;
			for (i=0;i<length;i++) {
				sum = 0.0;
				for (j=20;j--;) {
			    		q = amino_acids[j] - 'A';
			    		sum += (M1[i][q] - M[i][q]) * (M1[i][q] - M[i][q]);
			    	}
			    	distance += sqrt(sum);
			}
			
			#ifndef NDEBUG			
			printf("%s \titeration=%d \tK=%d \td(M1,M)=%f\n", initial_segment, iteration, Kmatches, distance);
			//print_PSSM(M1, 0);
			#endif
			printf("%d\t%d\t%f\t%f\n", k, iteration, 1.0*Kmatches/FP, distance);

			if (iteration++ >= theta) {
				#ifndef NDEBUG
				printf("FORCED CONVERGE: iteration >= %d\n", theta);
				#endif
				//break;
			}
			
			if (distance < epsilon || iteration >= theta) {
				converged = true;
				#ifndef NDEBUG
				printf("CONVERGED #%d!\n", converged_matrices_count);
				#endif

				if (Kmatches < K_min) {
					#ifndef NDEBUG
					printf("DISCARDED. LOW K: %d < %d!\n", Kmatches, K_min);
					#endif
					
					break;	
				}
				
				converged_matrices = realloc(converged_matrices, ++converged_matrices_count * sizeof(Matrix));
				matrix = converged_matrices + converged_matrices_count - 1;
				
				//strcpy(pmatrix->description, (initial_segments + k)->description);
				strncpy(matrix->initial_segment, initial_segment, MAX_PSSM_LENGTH + 1);
				strncpy(matrix->origin, origin, 61);
				matrix->K = Kmatches;
				matrix->L = length;
				matrix->N = N;
				matrix->p = p;
				matrix->id = k; //sequential number of the origin segment or matrix
				matrix->omega = omega;
				matrix->score = maxS;
				memcpy(&matrix->freq, &F, sizeof F);
			} else {
				memcpy(&M, &M1, sizeof M);
			}

			} // iterations. while not converged
			
			// live update of the derivation process
			write_Matrices("output.1.matrix", converged_matrices, converged_matrices_count);
			write_VariableMatrixFreq("output.4.matrix", converged_matrices, converged_matrices_count,  composition);
			
			if (origin_is_matrix) {
			    // do not iterate over segments if starting from matrix
			    break;
			}
		} // initial segment
	} // take another origin

	//printf("Final stage reached for rank %d.\n", rank);
	
	for (k=0; k < converged_matrices_count;k++) {
		matrix = converged_matrices + k;
		#ifndef NDEBUG
		printf("Matrix: %d %s %f\n", k, matrix->initial_segment, matrix->score);
		#endif
	}
	
	write_Matrices("output.1.matrix", converged_matrices, converged_matrices_count);
	write_VariableMatrixFreq("output.4.matrix", converged_matrices, converged_matrices_count, composition);

	#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	#endif
	return(0);
}
