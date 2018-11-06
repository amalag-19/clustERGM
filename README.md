# clustERGM
Community Detection in Static Networks through ERGMs

#Static undirected density network example

#Documentation for simulate_network_wrapper function. Inputs are:
1) N: Number of nodes
2) K: Number of clusters
3) T_data (optional and to be given only for dynamic networks): Number of time points
4) dynamic_indicator: Indicator whether the network is dynamic
5) directed_indicator: Indicator whether the network is directed
6) OR_indicator: Indicator whether the in the directed network outgoing and reciprocity should be combined or not. This is meaningful only when directed_indicator is 1.
7) density_indicator: Indicator whether the network is clustered based on density parameter
8) stability_indicator: Indicator whether the network is clustered based on stability parameter
9) transitivity_indicator: Indicator whether the network is clustered based on transitivity parameter
10) alpha: vector of length K defining the mixture proportions of clusters in the network. For the dynamic(HMM) case this indicates the probability vector for initial time point cluster memberships. No value needs to be given for K=1 case
11) pi: transition probabilty array with dimensions K*K*(T_data-1) for the density parameter case and K*K*(T_data-2) for the stability/transitivity parameter case
12) theta_dens: vector of length K giving the true density parameters. Note that for stability and transitivity cases also, this is required for simulating the initial time point network
13) theta_stab(optional): vector of length K giving the true stability parameters
14) theta_trans(optional): vector of length K giving the true transitivity parameters

#Output is a list, for all K>=2 cases, with 
1) first element as the 2/3 dimensional simulated network in the form of adjacency matrices and 
2) second element being the true cluster memberships arranged as N*(T_data) matrix for dynamic networks and N dimensional vector for static networks.
For K=1 case, the output is just the simulated network array.

#########################################################################################################
# Documentation for estimation/clustering for static undirected function. All other cases are similar. Inputs are:
1) sim.net: Array of the simulated network adjacency matrices.
2) nclust: Number of clusers for which the estimation must be performed.
3) thres: convergence threshold for estimation, the default being 10^(-6).
4) theta_init: Initial value of the density/stability/transitivity parameters for which clustering is desired
5) sim_indicator: Indicator whether the network is simulated or coming from real data
6) If sim_indicator is 1 (i.e. it is a simulated network), then following additional inputs must be provided to calculate clustering and estimation performance.
6.1) theta_true: vector of length K giving the true value of theta parameter which can be density/stability/transitivity
6.2) K_true: a number indicating true number of clusters for which the network was simulated
6.3) cluster_ids_true: vector (matrix) of dimensions N*1 (N*T_data) giving the true cluster memberships of nodes for static (dynamic) network

#Output is a list with following elements:
1) First element is again a list containing converged parameters. The elements of this list in the sequential order are gamma, alpha, pi, Tau, theta. For static networks there is no pi and Tau.
2) Second element is a vector (matrix) of estimated cluster memberships for static (dynamic) case. This is absent for all K=1 cases.
3) Third element is the integreted classification likelihood value that can be used for model selection i.e. selecting the appropriate number of clusters.
4) The last three elements optional given sim_indicator is set as 1. These elements are 
4.1) Rand index: giving clustering performance
4.2) RASE theta: giving estimation performance of theta (only when K_true is same as nclust)
4.3) RASE pi: giving estimation performance of pi (only for dynamic network case)
