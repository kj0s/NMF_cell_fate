# Cosensus.py code: 
you look for factors in agreement across multiple runs. hence you rn a model w diff seed on the same data multiple times. 

A **consensus factor** is a factor that represents the common pattern found across many model runs.

this is trying to get stable consensus factors from a rinse and repeat model instead of just trying to fit smth to noisy data. 

## Class: ConsensusFactorModel 
### __init__()
k = no. of latent factors

n_reps = how many times to run the model

random_states = seeds for reproducibility

self.models here creates mulktiple independent factor models. many copies of same model but each one starts from a random seed, learns slightly different factors from the same data

### fit()

fits each model on the same data

each model has ts own m.factors[...]
- each model gives multiple outputs, for each run you ge multiple matrices, 

### aggregate_factors()

makes each factor unit length, used for clustering

self.factor_run_id keeps track of which model run each factor came from

stacks everything together using np.hstack, one big matrix of all factors from all runs. shape = [total features, k * n_runs]

### detect_outliers()

removes bad or noisy factors. method is:
1. build wnn using self.factor_graph. finds nearest neighbour for each factor
2. compute avg dist to neighbours, for each factor compute dist to every neighbour, take mean; self.mean_dist_to_k
3. self.is_outlier used to find outliers, categorised by if smth is in top q of distances, like q = 90%, its an outlier.

### get_consensus_factors()

first prepare data by running aggregate, outliers.

run k_means only on non-outlier factors.

each factor gets a cluster id (0 to k-1) or a -1 if its an outlier

regroup original factors using factors_raw, compute consensus using self.cons_factors.
- for each cluster, we take all factors in that cluster, compute the median across them, 

# factors.py code:
this is the model thats being repeated

## Class: LinearFactorModel 
### __init__()

1. k = no. of latent factors
2. params = disctionary of weights, how much to care about each datast, how strong is regularisation
3. iters = max no of traning iterations
4. tol = convergence threshold i.e. stop if improvement is tiny
5. reg = regulariation / extra structure like graphs and penalties
6. random_state = same as above, seed control for reproducibility
7. print_iter/check_convergence/verbose = how often to print, how often to check for stopping
8. normalise_activities = whether to normalise latent variables [check above code]
9. b_factor_type = non negative tensor factorisation [3d ver of nmf]
10. learn_factors = which factor matrices are allowed to update
11. OUTPUT: factors, activities
12. S = superdiagonal tensor, 3d I.

### loss model(self)
penalties for bad structure, and check how wrong the model is. the goal is to minimise this term

### fit()
training the model.
multiple datasets are jointly modelled, introduces U = activities.

activities / U = latent representation of samples. rows = samples, columns = latent factors.

U [activites] = how much each factor is active in each sample
V* [factors] = what each factor looks like in a data space

#### iterative updates: 
for i in range(self.iters) = repeat until full convergence

you update activities using a MULTIPLICATIVE UPDATE RULE to keep values nonnegative; increase U when it helps explain data, decrease U where it hurts

