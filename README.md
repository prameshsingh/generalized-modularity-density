# Generalized Modularity Density

Maximization of generalized modularity density ($Q_g$) using RenEEL algorithm


If you use this code, please cite the following preprint:

Resolution limit revisited: community detection using generalized modularity density.
J. Guo, P. Singh, K. E. Bassler, arXiv:2012.14543; https://arxiv.org/abs/2012.14543

Reduced Network Extremal Ensemble Learning (RenEEL) Algorithm
J. Guo, P. Singh, K. E. Bassler, Scientific Reports 9, 14234 (2019).

The example network used here for illustration is the Zachary’s Karate Club network (W. W. Zachary, Journal of Anthropological Research. 33 (4): 452–473 (1977)).


To use the code, follow the steps below:

1. Prepare data

	* 1.1 Use an edgelist file including 3 columns (separated by space or tab with no header) with the last one representing the weight. (If the network is unweighted, just put 1 to the whole column. The network should be undirected. Self-loops will be ignored. See example in karate.txt)
	* 1.2 Use bash script (work.sh) to generate the three files required by the program. 
	* Example:

		`sh work.sh karate.txt` 


2. Compile code

	* 2.1 Make sure the openmp library is installed.
	* 2.2 compile main.c, help.c and rg.c with required libraries (openmp, math).
	* Example:

		`gcc-9 main.c help.c rg.c -fopenmp -lm`

		This will generate file a.out.

3. Run
	run with 6 arguments, respectively as:
	* argument 1: Positive Integer, parameter for Randomized Greedy  (usually 2)
	* argument 2: Positive Integer, ensemble size divided by CPU numbers (It's like this because we use parallel programmming to generate the ensemble faster. Suppose you are using 4 CPU cores and 5 here, your ensemble size will be 4*5=20)
	* argument 3: Positive Integer, ensemble size divided by CPU numbers of partitions of the reduced network for iteration part in RenEEL (Suppose you are using 4 CPU cores and 2 here, your ensemble size during iteration of RenEEL will be 4*2=8)
	* argument 4: Positive Integer, seed for random number generator
	* argument 5: Non-negative Float, for parameter $\chi$ in Generalized Modularity Density formula ($\chi$ = 0 for Modularity)
	* argument 6: Input file name (edgelist), for example `karate.txt`. The file should have already been processed using `work.sh` as described in 1.

	* Example:

		`./a.out 2 5 2 12345 1.0 karate.txt`

4. Collect results

	* file 1: `partition_karate.txt` (or partition_[inputfile])
	* file 2: `results_karate.txt` (or results_[inputfile]), a copy will also be printed to stdout
