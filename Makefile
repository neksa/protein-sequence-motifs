#CC=mpicc
CC=cc
#CFLAGS=-tp barcelona-64 -fastsse -Mipa=fast
#CHOST=i686-pc-linux-gnu
#nocona
#CFLAGS=-Wall 
#CFLAGS=-Wall -march=pentium4 -O2 -pipe -fomit-frame-pointer -funroll-loops
#CFLAGS=-Wall -march=k8 -O2 -pipe -fforce-addr -fomit-frame-pointer -funroll-loops -std=c99 -pedantic
CFLAGS=-Wall -O2 -pipe -fforce-addr -fomit-frame-pointer -funroll-loops -std=c99 -pedantic
#CFLAGS=-Wall -march=barcelona -O2 -pipe -fomit-frame-pointer -funroll-loops
CXXFLAGS=${CFLAGS}

all: converge search

converge: convergePSSM.c
	${CC} ${CFLAGS} convergePSSM.c -lm -o converge

# -D NDEBUG

#mpiconverge: convergePSSM.c
#	${CC} -D MPI -D NDEBUG ${CFLAGS} convergePSSM.c -lm -o converge

#debug: clean convergePSSM.c
#	${CC} -g ${CFLAGS} convergePSSM.c -lm -o converge

#mpidebug: clean convergePSSM.c
#	${CC} -g ${CFLAGS} -D MPI  convergePSSM.c -lm -o converge

#clustering: clean clustering.c
#	${CC} -g ${CFLAGS} clustering.c -lm -o clustering

clean:
	rm -f converge
	rm -f search

search: searchPSSM.c
	rm -f search
	${CC} ${CFLAGS} searchPSSM.c -lm -o search

#searchdebug: searchPSSM.c
#	rm -f search
#	${CC} -D MPI -g ${CFLAGS} searchPSSM.c -lm -o search

#searchcontrol: searchPSSM_control.c
#	rm -f search_control
#	${CC} -g ${CFLAGS} searchPSSM_control.c -lm -o search_control

#reshuffle_matrix: reshuffle_matrix.c
#	rm -f reshuffle_matrix
#	${CC} -D NDEBUG ${CFLAGS} reshuffle_matrix.c -lm -o reshuffle_matrix

#convert_matrix_50_30: convert_matrix_50_30.c
#	rm -f convert_matrix_50_30
#	cc -D NDEBUG ${CFLAGS} convert_matrix_50_30.c -lm -o convert_matrix_50_30
    