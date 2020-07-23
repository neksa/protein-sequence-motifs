

find_permutations <- function(size=50, threshold_cross_matches=3, permutations=10) {
    for (k in 1:permutations) {
	s <- 0
	counter <- 100
	while (counter > threshold_cross_matches) {
		s <- sample(1:size - 1)
		#cat(s, "\n")
		counter <- 0
		for (i in 1:size) { 
			for (j in 1:size) {
				if ( (j>i) && (s[j] > s[i]) && (j-i == s[j]-s[i]) ) {
					#cat ("bad match: ", j-i, "=", j, "-", i, s[j], "-", s[i], "\n");
					counter <- counter+1;
				}
			}
		}
		#cat(counter, "\n")
	}

	cat("\nPermutation #", k, "\n")
	cat(counter, "\n")
	cat(s, "\n");

	for (i in 1:size) { 
		for (j in 1:size) {
			if ( (j>i) && (s[j] > s[i]) && (j-i == s[j]-s[i]) ) {
				cat ("bad match: ", j-i, "=", j, "-", i, s[j], "-", s[i], "\n");
			}
		}
	}
    }
}

find_permutations(30, 1, 10)
#find_permutations(50, 3, 10)
