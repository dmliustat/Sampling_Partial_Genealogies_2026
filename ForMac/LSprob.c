#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tables.h"
#include "header.h"
#include "gamma.h"
#include "common.h"
#include "wu_alg.h"


extern long M; /*number of genes observed at present*/
extern long L; /*number of loci*/
extern long NRUN;
extern double SEQLEN;
extern double THETA;
extern long NumLineage; /*The number of lineages to stop*/
extern double *positions;
extern long Ne;
extern long* ReTime; /*time (number of lineages) to resample*/
extern long Resample;


extern double calc_Wu(long* TYPE, long* ntype);

int Build(long* TYPE, long* type2node, long* ntype, const double THETA, double* logw, tsk_table_collection_t* tables, long* num_nodes, long* num_rec, long* num_mut, long* num_coal, long* num_nonrec, long* mut_by_site, long* mut_site, long* mut_state, long* mut_node, double* pHtau); /*success: return 0; otherwise, return 1*/

void FDupdate(long* TYPE, long* type2node, const double THETA, long k, long* ntype, long* nbranch, long* nanc, double* m, double* logw, tsk_table_collection_t* tables, int *rec_num, int *mut_num, int *coal_num, int* nonrecc, long* mut_by_site, long* mut_site, long* mut_state, long* mut_node, long* eligible);

void Resamplef(long* TYPE, long* nbranch, double* logw, tsk_table_collection_t* tables);

void SampleTypes();

double rexp(double rate);

int generateRandomValue(int min, int max);


void Resamplef(long* TYPE, long* nbranch, double* logw, tsk_table_collection_t* tables) {

	SampleTypes();
	
}


int Build(long* TYPE, long* type2node, long* ntype, const double THETA, double* logw, tsk_table_collection_t* tables, long* num_nodes, long* num_rec, long* num_mut, long* num_coal, long* num_nonrec, long* mut_by_site, long* mut_site, long* mut_state, long* mut_node, double* pHtau) {
    long *nanc, nbranch = M, *newseq, eligible, mut_elig; /*nanc: number of individuals ancestral at each locus */
    long *nanc_copy, *TYPE_tmp, ntype_copy, nlineage_copy, nbranch_copy; 
    long i, j, t, l, z, count, min, max, s1, s2, index_re=0;
    short int flag = 0; /*flag == 0, continue SIS; flag != 0, the stopping criteria is reached*/
    short int flag_MRCA = 0; /*flag_MRA == 0, the grand MRCA not found, continue loops*/
    //short int flag_H0 = 0; /*flag_H0 == 0, calculate the exponential rate for H_0*/
    double branchProb[TYPE_MAX], logH_tau=0.0, apprx;
	int u;
    double tot, temp, m[2];
    unsigned long seed_copy[2]; /* store initial seed; print out if error */
	//int rec_num = 0, mut_num = 0, coal_num = 0;
	long *shuffle_index, * shuffle_seq, *TYPE_flat, *elig_type;
	double logPAC, *logPAC_1, logPAC_avg, PAC_temp, logw_inner = 0.0;
	long nshuffle=1000;
	long* TYPE_st, ntype_st, nbranch_st, * nanc_st, * mut_by_site_st;
	double wu_prob;
	FILE* test_wt;
	char filename[200];

    nanc = long_vec_init(L);
    nanc_copy = long_vec_init(L);
    newseq = long_vec_init(L);


	TYPE_tmp = long_vec_init((L + 1) * TYPE_MAX);
	TYPE_flat = long_vec_init((L + 1) * NumLineage);
	shuffle_seq = long_vec_init((L + 1) * NumLineage);
	shuffle_index = long_vec_init(NumLineage);
	logPAC_1 = dou_vec_init(nshuffle);
	elig_type = long_vec_init(TYPE_MAX);


    for(i=0; i<L; i++) {
        *(nanc+i) = M;
    }
    /*for(i=0;i<L-1;i++) rec_no[i] = 0;*/

    *logw = 0.0; /*log SIS weight */
	*num_nodes = 0;
	*num_rec = 0;
	*num_mut = 0;
	*num_coal = 0;
	*num_nonrec = 0;
	*pHtau = 0.0;

    *seed_copy = *seed;
    *(seed_copy+1) = *(seed+1);

    while((!flag) && (!flag_MRCA)){
        /*stage 1: sample a lineage*/
		/*count the number of eligible sequences*/
		eligible = 0;
		t = 0;
        for(i=0; i<*ntype; i++){
			/*coalescence*/
			if (*(TYPE + (L + 1) * i + L) > 1) {
				eligible += *(TYPE + (L + 1) * i + L);
				*(elig_type + t) = i;
				t++;
				//fprintf(stderr, "Type %ld is eligible for coalescence, multiplicity = %ld.\n", i, *(TYPE + (L + 1) * i + L));
			}
			/*mutation*/
			else {
				/*multiplicity = 1*/
				if (*(TYPE + (L + 1) * i + L) != 1) {
					fprintf("ERROR! Multiplicity = %ld.\n", *(TYPE + (L + 1) * i + L));
					return(1);
				}
				else {
					mut_elig = 0;
					for (j = 0; j < L; j++) {
						if ((mut_by_site[j] == 1) && (*(TYPE + (L + 1) * i + j) == 1)) mut_elig += 1;
					}
					if (mut_elig > 0) {
						eligible += 1;
						//fprintf(stderr, "Type %ld is eligible, %ld mutation sites.\n", i, mut_elig);
						*(elig_type + t) = i;
						t++;
					}	
		
				}				
			}
        }
        if(eligible <= 0){
            fprintf(stdout, "Error. Number of eligible sequences = %ld.\n", eligible);
            return(1);
        }
		//fprintf(stdout, "%ld eligible sequences.\n", eligible);

        /*start sampling a branch (type) out of the current configuration*/
		
        u = generateRandomValue(1, eligible);
		//u = ceil(eligible * runif());
		if (u > eligible || u < 1) {
			fprintf(stderr, "Error generating random value!");
			return(1);
		}
		//fprintf(stderr, "eligible = %d, random value = %d\n", eligible, u);
        i = 0;
        while(u > 0){
			j = elig_type[i];
            u -= *(TYPE + (L + 1) * j + L);
            i++;
        }
        i--; 
		//fprintf(stderr, "The %ldth type is chosen. ", i);

		j = elig_type[i]; /*j is the chosen type*/
		//fprintf(stderr, "Type %ld is chosen. \n", j);

		/* Resample if the resampling time is reached*/
		if (nbranch == ReTime[index_re]) {
			/* resample the current configuration*/
			Resamplef(TYPE, &nbranch, logw, tables);

			index_re++;
			if (index_re > Resample) {
				fprintf(stderr, "Resampling error!");
				abort();
			}
		}
       
        if(nbranch == NumLineage){
			//fprintf(stderr, "Stopping time reached. logw = %lf\n", *logw);
			fprintf(stderr, "\nStopping criterion is met.\n");
			flag = 1;
        }
        else{
			/*Given the chosen chromosome j, select the genetic event. Update configuration and logw*/
			FDupdate(TYPE, type2node, THETA, j, ntype, &nbranch, nanc, m, logw, tables, num_rec, num_mut, num_coal, num_nonrec, mut_by_site, mut_site, mut_state, mut_node, &eligible);
			if (*logw == -INFINITY) {
				fprintf(stderr, "logw = -inf!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			}
		}

        /*check nbranch and ntype*/
        tot = 0;
        for(i=0; i<*ntype; i++) tot += *(TYPE+(L+1)*i+L);
		if ((long)tot != nbranch) {
			fprintf(stdout, "Error in the number of branches. nbranch = %ld, sum of # types = %ld\n", nbranch, (long)tot);
			return(1);
		}

        if(*ntype>TYPE_MAX){
            fprintf(stdout, "Number of distinct haplotypes exceeds the maximum\n");
            return(1);
        }
		*num_nodes = nbranch;
        /*flag_MRCA = 1 if the GMRCA is reached*/
        flag_MRCA = 1;
        for(j=0; j<L; j++){
            flag_MRCA *= (*(nanc+j)==1);
        }
    }

	fprintf(stderr, "num of nodes = %ld, recombinations = %ld; mutations = %ld; non-recurrent = %ld; coalescence = %ld\n", nbranch, *num_rec, *num_mut, *num_nonrec, *num_coal);
	fprintf(stderr, "Before calculate H_tau, the log(weight) = %lf\n", *logw);

	
	wu_prob = calc_Wu(TYPE, ntype);
	*logw += log(wu_prob);

	*pHtau = wu_prob;
	

    free(nanc);
    free(nanc_copy);
	free(newseq);
	free(elig_type);
    
    return(0);
}







void FDupdate(long* TYPE, long* type2node, const double THETA, long k, long* ntype, long* nbranch, long* nanc, double* m, double* logw, tsk_table_collection_t* tables, int *rec_num, int *mut_num, int* coal_num, int* nonrec, long* mut_by_site, long* mut_site, long* mut_state, long* mut_node, long* eligible)
{
    long i, j, a, e, *newseq, min, max, *TYPE_st, n_alpha_ip1;
    long onode, onode2, nnode, nnode2, edge, mut, *edge_arr, nodeTime, mut_n, *mut_s, multi;
    double tot, u, left, right; /* edge_sum = 0.0, rec_sum = 0.0 for edges*/
	int flag, flag_first = 0, flag_last = 0;
	char state[2];

    newseq = long_vec_init(L);
	edge_arr = long_vec_init(L);
	TYPE_st = long_vec_init((L + 1) * TYPE_MAX);
	mut_s = long_vec_init(L);


    /*onode is the index of the chosen/old node sampled from type k nodes*/
    onode = sample_node(TYPE, type2node, k);

    for(i=0; i<L; i++) {
        newseq[i] = *(TYPE+(L+1)*k+i);
    }

	for (e = 0; e < (L + 1) * *ntype; e++) TYPE_st[e] = *(TYPE + e);

    /*mutation events*/
	if (*(TYPE + (L + 1) * k + L) == 1) {
		fprintf(stderr, "Mutation occurs.\n");
		/*update samples*/
		/*onode is the index of the chosen/old node sampled from type k nodes*/
		onode = sample_node(TYPE, type2node, k);

		mut_n = 0;
		for (i = 0; i < L; i++) {
			if ((mut_by_site[i] == 1) && (*(TYPE + (L + 1) * k + i) == 1)) {
				mut_s[mut_n] = i;
				mut_n += 1;
			}
		}

		/*find the mutation site*/
		if (mut_n == 0) {
			fprintf(stderr, "Error! Not eligible sequence.\n");
			return;
		}
		if ((mut_n) == 1) {
			j = mut_s[0];
			newseq[j] = 0;
		}
		else {
			j = mut_s[(int)floor(runif() * mut_n)];
			newseq[j] = 0;			
		}
		*mut_num += 1;
		*(mut_by_site + j) -= 1;
		remove_old(TYPE, type2node, ntype, k, onode, *nbranch);
		*nbranch -= 1;
		i = add_new(TYPE, type2node, ntype, newseq, onode, *nbranch);
		*nbranch += 1;

		/*record the mutation site & allele*/
		mut_site[*mut_num - 1] = j;
		mut_state[*mut_num - 1] = 0;
		mut_node[*mut_num - 1] = onode;

		/*update weight*/
		n_alpha_ip1 = *(TYPE + (L + 1) * i + L);
		*logw += log(n_alpha_ip1 * THETA * *eligible / (*nbranch * (*nbranch - 1.0 + THETA)));
		if (*logw == -INFINITY || isnan(*logw) == 1) {
			fprintf(stderr, "Error!!! Mutation event logw = -inf or nan\n");
		}
		//fprintf(stderr, "After mutation, logw = %lf\n", *logw);
		//fprintf(stderr, "Mut: nbranch = %ld, ntype = %ld.\n", *nbranch, *ntype);

	}

    
    /*coalescence events*/
	else {
		fprintf(stderr, "Coalescence occurs.\n");
		multi = *(TYPE + (L + 1) * k + L);
		*coal_num += 1;
		nodeTime = *rec_num + *mut_num + *coal_num;
		n_alpha_ip1 = *(TYPE + (L + 1) * k + L) - 1;

		/*update samples*/
		/*For this moment, there is no recombination, so only the sequences of the same type can coalesce.*/
		/*However, with recombination events, different types can coalesce as long as they are identical on the ancestrol loci*/
		i = k; /*i is the type that coalesce with k*/
		for (j = 0; j < L; j++) {
			if ((*(TYPE + (L + 1) * i + j) > 0) && (*(TYPE + (L + 1) * k + j) > 0)) *(mut_by_site + j) -= 1;
		}
		
		/*sample the sequences that coalesce*/
		onode = sample_node(TYPE, type2node, k);
		onode2 = sample_node(TYPE, type2node, i);
		while (onode == onode2) {
			onode2 = sample_node(TYPE, type2node, i);
		}
		//fprintf(stderr, "Coalesce between node %ld (type %ld) and node %ld (type %ld)\n", onode, k, onode2, i);
		if (onode > onode2) {
			onode = onode + onode2;
			onode2 = onode - onode2;
			onode = onode - onode2;
			k = k + i;
			i = k - i;
			k = k - i;
		}
		nnode = tsk_node_table_add_row(&tables->nodes, 0, nodeTime, TSK_NULL, TSK_NULL, NULL, 0);
		/*edge_arr: 0 non-ancestral; 1 type k; 2 type i; 3: both*/
		for (j = 0; j < L; j++) {
			if (*(TYPE + (L + 1) * k + j) < 0 && *(TYPE + (L + 1) * i + j) >= 0) {
				newseq[j] = *(TYPE + (L + 1) * i + j);
				edge_arr[j] = 2;
			}
			else if (*(TYPE + (L + 1) * k + j) >= 0 && *(TYPE + (L + 1) * i + j) < 0) {
				newseq[j] = *(TYPE + (L + 1) * k + j);
				edge_arr[j] = 1;
			}
			else if (*(TYPE + (L + 1) * k + j) < 0 && *(TYPE + (L + 1) * i + j) < 0) {
				newseq[j] = -1;
				edge_arr[j] = 0;
			}
			else {
				/* *(TYPE+(L+1)*k+j)>=0 && *(TYPE+(L+1)*i+j)>=0 && they're identical on locus j */
				newseq[j] = *(TYPE + (L + 1) * k + j);
				edge_arr[j] = 3;
				nanc[j] -= 1;
			}
		}

		/*add to the edge table*/
		/*within each given parent, the edges must be sorted by child index and then left*/
		/*store nnode1 edges*/
		left = -1;
		flag = -1; /*flag=-1 haven't found the interval; flag=0 found left*/
		for (j = 0; j < L; j++) {
			if (flag == -1 && j < L - 1) {
				if (edge_arr[j] == 1 || edge_arr[j] == 3) {
					left = positions[j];
					flag = 0;
				}
			}
			else if (flag == 0) {
				if (j == L - 1) {
					right = positions[j];
					edge = tsk_edge_table_add_row(&tables->edges, left, right, nnode, onode, NULL, 0);
					if (edge < 0) {
						fprintf(stderr, "Adding edge failed.\n");
					}
				}
				else if (edge_arr[j] == 0 || edge_arr[j] == 2) {
					right = positions[j];
					edge = tsk_edge_table_add_row(&tables->edges, left, right, nnode, onode, NULL, 0);
					if (edge < 0) {
						fprintf(stderr, "Adding edge failed.\n");
					}
					flag = -1;
				}
			}
		}

		/*store nnode2 edges*/
		left = -1;
		flag = -1; /*flag=-1 haven't found the interval; flag=0 found left*/
		for (j = 0; j < L; j++) {
			if (flag == -1 && j < L - 1) {
				if (edge_arr[j] == 2 || edge_arr[j] == 3) {
					left = positions[j];
					flag = 0;
				}
			}
			else if (flag == 0) {
				if (j == L - 1) {
					right = positions[j];
					edge = tsk_edge_table_add_row(&tables->edges, left, right, nnode, onode2, NULL, 0);
					if (edge < 0) {
						fprintf(stderr, "Adding edge failed.\n");
					}
				}
				else if (edge_arr[j] == 0 || edge_arr[j] == 1) {
					right = positions[j];
					edge = tsk_edge_table_add_row(&tables->edges, left, right, nnode, onode2, NULL, 0);
					if (edge < 0) {
						fprintf(stderr, "Adding edge failed.\n");
					}
					flag = -1;
				}
			}
		}

		add_new(TYPE, type2node, ntype, newseq, nnode, *nbranch);
		*nbranch += 1;
		if (*(TYPE + (L + 1) * k + L) == 1 && k < i) i--;
		remove_old(TYPE, type2node, ntype, k, onode, *nbranch);
		*nbranch -= 1;
		remove_old(TYPE, type2node, ntype, i, onode2, *nbranch);
		*nbranch -= 1;

		/*update weight*/
		if (*nbranch > 1) {
			//*logw += log((*nbranch - 1.0) * *eligible * (multi - 1.0) / (2.0 * *nbranch * (*nbranch + THETA - 1.0)));
			*logw += log(n_alpha_ip1 * (*nbranch - 1.0) * *eligible / (*nbranch * (*nbranch - 1.0 + THETA) * (n_alpha_ip1 + 1.0)));
		}
	
		//fprintf(stderr, "After coalescence, logw = %lf\n", *logw);
		//fprintf(stderr, "Coal: nodes %ld and %ld, nbranch = %ld, ntype = %ld.\n", onode, onode2, *nbranch, *ntype);
		if (*logw == -INFINITY || isnan(*logw) == 1) {
			fprintf(stderr, "Error!!! Coalescence logw = -inf or nan\n");
		}

	}


    free(newseq);
	free(edge_arr);
	free(TYPE_st);
	free(mut_s);
}



double rexp(double rate){
    double r;
    r=runif();
    return (-log(r)/rate);
}


// Function to generate a random value between min and max (inclusive) 
int generateRandomValue(int min, int max) { 
	
	// Generate a random value between min and max (inclusive) 
	int randomValue = min + rand() % (max - min + 1); 
	return randomValue; 
}

