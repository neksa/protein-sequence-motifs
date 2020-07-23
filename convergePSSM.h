

int load_fasta_sequences(const char *filename, Sequence **psequences, bool removex) {
	FILE *fd = NULL;
	char *seq = NULL;
	char desc[100];
	char *line = NULL;
	int line_len = 0;
	int seq_len = 0;
	int N = 0;
	int i = 0;

	fd = fopen(filename, "r");
	if (fd == NULL) {
		error(1, errno, "error while opening file %s", filename);
	}
	assert(fd != NULL);

	char *buf = (char *)calloc(4096, sizeof(char));

	for (;;) {
		line = fgets(buf, 4096, fd);
		//if (line[0] == 0) continue;
		
		if (line != NULL) {
			line_len = strlen(line);
			
			if (line[line_len-1] == '\n') {
				line[line_len-1] = 0x00;
				line_len--;
			}
			if (line[line_len-1] == 0x0D) {
				line[line_len-1] = 0x00;
				line_len--;
			}
		}
		
		if (seq) {
		    seq_len = strlen(seq);
		}
		
		// end of block
		if (line == NULL || line[0] == '>') {
			if (seq) {
				for (i=0; i < seq_len; i++) {
					seq[i] = (char)toupper((int)seq[i]);

					// replace all occurences non-natural amino acids with X:
				    	//    A CDEFGHI KLMN PQRST VW Y
				    	//     B       J    O     U  X Z
				    	//     66      74   79    85 8890
					if ('B' == seq[i] || 'J' == seq[i] || 'O' == seq[i] || 'U' == seq[i] || 'Z' == seq[i]) {
					     seq[i] = 'X';
					}
					if (seq[i] < 'A' || seq[i] > 'Z') {
					    seq[i] = 'X';
					}
									
					if (removex) {
						// remove all proteins with X:
						if ('X' == seq[i]) {
							free(seq);
							seq = NULL;
							break;
						}
					}
				}
				if (seq) {
					N++;
					*psequences = realloc(*psequences, N * sizeof(Sequence));
					strncpy((*psequences + N - 1)->description, desc, 100);
					(*psequences + N - 1)->description[99] = 0x00;
					(*psequences + N - 1)->sequence = seq;
					seq = NULL;
					seq_len = 0;
				}
			}
			if (line == NULL) {
				break; // EOF
			}
			//> fasta header
			strncpy(desc, line + 1, 100);
			desc[99] = 0x00;
			continue;
		}
		
		if (line_len == 0) {
			continue;
		}
						
		if (seq) {
			seq = realloc(seq, seq_len + line_len + 1);
			strcpy(seq + seq_len, line);
		} else {
			seq = calloc(line_len + 1, sizeof(char));
			strcpy(seq, line);
		}
	}
	free(buf);
	fclose(fd);
	return N;
}

void print_PSSM(double M[MAX_PSSM_LENGTH][26], double threshold) {
	int i=0, j=0, a=0;
	printf("PSSM");
	
	for (j=0; j<20; j++) {
		printf(" %c    ", amino_acids[j]);
	}
	printf("\n");
	for (i=0; i<MAX_PSSM_LENGTH; i++) {
		printf("%-2d ", i);
		for (j=0; j<20; j++) {
			a = amino_acids[j] - 'A';
			if (M[i][a] > threshold) {
				printf("% 1.2f ", M[i][a]);
			} else {
				printf("  .   ");
			}
		}
		printf("\n");
	}
	printf("\n");
}

void write_Matrices(const char *filename, Matrix *matrices, int n, int suffix) {
	FILE *fd = NULL;
	const int proto_format_version = 1;
	int m=0, i=0, j=0;
	Matrix *M=NULL;
	char full_filename[256] = "";

	if (n <= 0){
		return;
	}
	
	snprintf(full_filename, 256, "%s.%d", filename, suffix);

	fd = fopen(full_filename, "w");
	assert(fd != NULL);
	
	for (m=0; m<n; m++) {
		M = matrices + m;
		
		fprintf(fd, "PROTOTYPE %d\n", proto_format_version);
		fprintf(fd, "BEGIN\n");
		fprintf(fd, "SEGMENT %s\n", M->initial_segment);
		fprintf(fd, "MATRIX K=%d N=%d P=%1.8lf S=%f W=%f\n", M->K, M->N, M->p, M->score, M->omega);
	
		fprintf(fd, "50 ");
		for (j=0; j<strlen(amino_acids); j++) {
			fprintf(fd, "    %c ", amino_acids[j]);
		}
		fprintf(fd, "\n");
		for (i=0; i<MAX_PSSM_LENGTH; i++) {
			fprintf(fd, "%2d ", i);
			for (j=0; j<strlen(amino_acids); j++) {
				assert(M->freq[i][amino_acids[j] - 'A'] < 100000);
				fprintf(fd, "%5d ", (int)M->freq[i][amino_acids[j] - 'A']);
			}
			fprintf(fd, "\n");
		}
		fprintf(fd, "END\n");
	}
	fclose(fd);
}

// write composition to file
void write_Composition(const char *filename, double composition[20]) {
	FILE *fd = NULL;
	int j = 0;

	fd = fopen(filename, "w");
	assert(fd != NULL);
	
	for (j=20;j--;) {
	    fprintf(fd, "%c,%f\n", amino_acids[j], composition[j] * 100);
	}
	
	fclose(fd);
}

// returns number of loaded matrices
// matrix containing gaps, and counts
int load_VariableMatricesCount(const char *filename, Matrix **pmatrices, double composition[20]) {
	const int profile_format_version = 3;
	char *buf = calloc(256, sizeof(char));
	int m=0, j=0;
	int pos=0, tmp;
	int L;
	double b[21];
	FILE *fd = NULL;
	Matrix M;
	bool skip_profile = true;

	if (*pmatrices != NULL) {
		printf("Input error. load_Matrices expects an unallocated pointer to matrices (NULL)\n");
		return 0;
	}

	fd = fopen(filename, "r");
	assert(fd != NULL);

	while (fgets(buf, 256, fd)) {
		if (1 == sscanf(buf, "PROFILE %d", &tmp)) {
			if (tmp != profile_format_version) {
				printf("Format error. Wrong profile format version! Expecting version %d.\n", profile_format_version);
				return 0;
			}
		}
		if (strstr(buf, "BEGIN")) {
			for (pos=0; pos<MAX_PSSM_LENGTH; pos++) {
				for (j=20;j--;) {
					M.freq[pos][amino_acids[j] - 'A'] = composition[j];
				}
			}
			skip_profile = true;
		}
		if (strstr(buf, "END")) {
			//printf("Sizeof Matrix = %d\n", m*sizeof(Matrix));
			if (skip_profile) {
				printf("Matrix header not recognized, ignoring profile\n");
			} else {
				m++;
				*pmatrices = realloc(*pmatrices, m*sizeof(Matrix));
				memcpy(*pmatrices + m - 1, &M, sizeof(Matrix));
			}
		}
		if (2 == sscanf(buf, "MATRIX K=%d L=%d", &M.K, &L)) {
			// no conditions
			M.id = m;
			//printf("Matrix %d %d %d\n", m, M.K, L);
			skip_profile = false;
		}
		if (22 == sscanf(buf,
		    "%4d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		     &pos,&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7],&b[8],&b[9],&b[10],&b[11],&b[12],&b[13],&b[14],&b[15],&b[16],&b[17],&b[18],&b[19],&b[20])) {

			if (M.K <= 0) {
			    continue;
			}
		     
			if (pos<0 || pos>49) {
				printf("Format error: invalid position %d (borders 0-49)\n", pos);
				continue;
			}
			for (j=20; j--;) {
				if (b[j]<0) {
					printf("Format error: invalid position %d,%d = %lf\n", pos, j, b[j]);
				}
				M.freq[pos][amino_acids[j] - 'A'] = b[j] / (double)M.K;
			}
		}
	}
	//printf("Finished parsing file with matrices\n");
	fclose(fd);
	free(buf);
	return m;
}

/*
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
*/

// load BLOSUM62
void load_BLOSUM(const char *filename, int (*BLOSUM)[20]) {
	char *buf = calloc(256, sizeof(char));
	int i=0, j=0, k=0;
	int b[24];
	char A, B;
	bool skip = true;
	
	char listed_amino_acids[26];
	FILE *fd = NULL;

	fd = fopen(filename, "r");
	assert(fd != NULL);

	for (i=0; i< strlen(amino_acids); i++){
		for (j=i; j<strlen(amino_acids); j++) {
		    BLOSUM[i][j] = 0;
		}
	}
	while (fgets(buf, 256, fd)) {
		if (strstr(buf, "#") == buf) continue;

		if (strstr(buf, "  ") == buf) {
		    //printf("DEBUG: %s\n", buf);
		    i = 0;
		    for (j=0; j<strlen(buf); j++) {
			if (buf[j] != ' ') {
			    listed_amino_acids[i++] = buf[j];
			    //printf("DEBUG2: %c\n", buf[j]);

			}
		    }
		    
		    for (i=0; i<26; i++){
			//printf("DEBUG listed_aa[%d] = %c\n", i, listed_amino_acids[i]);
		    }
		    
		    continue;
		}
			
		if (25 == sscanf(buf,
		    "%c %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d", 
		     &A,&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7],&b[8],&b[9],&b[10],
		        &b[11],&b[12],&b[13],&b[14],&b[15],&b[16],&b[17],&b[18],&b[19],&b[20],&b[21],&b[22],&b[23])) {

			//printf("DEBUG 3! %c %2d %2d\n", A, b[0], b[1]);

			skip = true;
			for (i=0; i < strlen(amino_acids); i++){
			    if (amino_acids[i] == A) {
				skip = false;
				break;
			    }
			}
			if (skip) continue;
			
			for (k=0; k < 24; k++){
				skip = true;
				B = listed_amino_acids[k];
				for (j=0; j < strlen(amino_acids); j++){
				    if (amino_acids[j] == B) {
					skip = false;
					break;
				    }
				}
				if (skip) continue;				
				// we've got i,j here in ABCD coordinates
				//printf("A=%c B=%c ; BLOSUM[%d][%d] = %d\n", A, B, i, j, b[k]);
				BLOSUM[i][j] = b[k];
			}
		    }
	}
	//printf("Finished parsing BLOSUM\n");
	free(buf);

	/*
	printf(" ");
	for (i=0; i< strlen(amino_acids); i++){
		printf("  %c", amino_acids[i]);
	}
	printf("\n");
	for (i=0; i< strlen(amino_acids); i++){
		printf("%c ", amino_acids[i]);
		for (j=0; j<strlen(amino_acids); j++) {
		    printf("%2d ", BLOSUM[i][j]);
		}
		printf("\n");
	}
	*/
	fclose(fd);
}

// returns number of loaded matrices
int load_VariableMatricesFreq(const char *filename, Matrix **pmatrices, double composition[20]) {
	const int profile_format_version = 4;
	char *buf = calloc(256, sizeof(char));
	char *origin = NULL;
	int m=0, j=0;
	int pos=0, tmp;
	double b[20];
	FILE *fd = NULL;
	Matrix M;

	for (pos=0; pos<MAX_PSSM_LENGTH; pos++) {
		for (j=26;j--;) {
			M.freq[pos][j] = 0.0;
		}
	}

	if (*pmatrices != NULL) {
		printf("Input error. load_Matrices expects an unallocated pointer to matrices (NULL)\n");
		return 0;
	}

	fd = fopen(filename, "r");
	assert(fd != NULL);

	int check = 0;

	while (fgets(buf, 256, fd)) {
		if (1 == sscanf(buf, "PROFILE %d", &tmp)) {
			if (tmp != profile_format_version) {
				printf("Format error. Wrong profile format version! Expecting version %d.\n", profile_format_version);
				return 0;
			}
			check++;
		}
		if ((origin = strstr(buf, "ORIGIN ")) != NULL) {
		    origin += strlen("ORIGIN ");
		    char *eol = strstr(origin, "\n");
		    if (eol != NULL) {
			strncpy(M.origin, origin, MIN(eol-origin, 61));
			M.origin[60] = '\0';
			printf("ORIGIN: %s\n", M.origin);
		    }
		}
		if (strstr(buf, "BEGIN")) {

			for (pos=0; pos<MAX_PSSM_LENGTH; pos++) {
				for (j=20;j--;) {
					M.freq[pos][amino_acids[j] - 'A'] = composition[j];
				}
			}
			memset(M.initial_segment, '\0', MAX_PSSM_LENGTH + 1);
			memset(M.origin, '\0', 61);
			
			check++;
		}
		if (strstr(buf, "END")) {
			if (check < 3) {
				printf("PROFILE 4: load_VariableMatricesFreq detected an invalid profile in your input! only %d checks passed\n", check);
				return 0;
			}
			check = 0;
			
			//printf("Sizeof Matrix = %d\n", m*sizeof(Matrix));
			*pmatrices = realloc(*pmatrices, (m+1)*sizeof(Matrix));
			memcpy(*pmatrices + m, &M, sizeof(Matrix));
			m++;
		}
		if (3 == sscanf(buf, "MATRIX ID=%d K=%d L=%d", &M.id, &M.K, &M.L)) {
		    	//printf("PROFILE FORMAT 4. Read matrix %d, ID=%d K=%d L=%d\n", m, M.id, M.K, L);
			// no conditions
			if (M.L > 50) {
			    printf("Matrix length L=%d. Maximal allowed length is 50 residues!\n", M.L);
			    return 0;
			}

			check++;
		}
		if (21 == sscanf(buf,
		    "%2d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		     &pos,&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7],&b[8],&b[9],&b[10],&b[11],&b[12],&b[13],&b[14],&b[15],&b[16],&b[17],&b[18],&b[19])) {
		     
			if (pos<0 || pos>49 || pos>M.L-1) {
				printf("Format error: invalid position %d (hard borders 0-49) (matrix borders 0-%d) \n", pos, M.L);
				continue;
			}
			for (j=20; j--;) {
				if (b[j]<0 || b[j]>1) {
					printf("Format error: invalid position %d,%d = %lf\n", pos, j, b[j]);
				}
				M.freq[pos][amino_acids[j] - 'A'] = b[j];
			}
		}
	}
	//printf("Finished parsing file with matrices\n");
	fclose(fd);
	free(buf);
	return m;
}


void write_VariableMatrixFreq(const char *filename, Matrix *matrices, int n, int suffix, double composition[20]) {
	FILE *fd = NULL;
	int m=0, i=0, j=0;
	Matrix *M=NULL;
	char full_filename[256] = "";
	double summ = 0.0;
	double f = 0.0;
	int q;

	if (n <= 0){
		return;
	}
			
	snprintf(full_filename, 256, "%s.%d", filename, suffix);
	fd = fopen(full_filename, "w");
	assert(fd != NULL);
	
	for (m=0; m<n; m++) {
		M = matrices + m;
		
		fprintf(fd, "PROFILE 4\n");
		fprintf(fd, "BEGIN\n");
		fprintf(fd, "ORIGIN %s\n", M->origin);
		fprintf(fd, "MATRIX ID=%d K=%d L=%d\n", M->id, M->K, M->L);
	
		fprintf(fd, "%2d", M->L);
		for (j=0; j<strlen(amino_acids); j++) {
			fprintf(fd, "        %c", amino_acids[j]);
		}
		fprintf(fd, "\n");
		for (i=0; i < M->L; i++) {
			fprintf(fd, "%2d", i);

			summ = 0.0;

			for (j=strlen(amino_acids);j--;) {
			    q = amino_acids[j] - 'A';
			    summ += M->freq[i][q]; // it is a count, actually
			}
			
			for (j=0; j<strlen(amino_acids); j++) {
			    q = amino_acids[j] - 'A';
			    if (summ == 0) {
				f = composition[j];
			    } else {
				f = M->freq[i][q] / summ; // now it's a normalized frequency
			    }
			    //printf("DDDDDDDDDDDDD i=%d, j=%d, K=%d, SUMM=%lf, M->f=%lf, f=%lf\n", i, j, M->K, summ, M->freq[i][q], f);
			    assert(f <= 1.0); 
			    fprintf(fd, " %1.6lf", f); 
			}
			fprintf(fd, "\n");
		}
		fprintf(fd, "END\n");
	}
	fclose(fd);
}
