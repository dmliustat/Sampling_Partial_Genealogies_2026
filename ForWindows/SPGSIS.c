#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <direct.h>
#include "tables.h"
#include "common.h"
#include "header.h"
#include "gamma.h"
#include "wu_alg.h"
#include <time.h>

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        fprintf(stderr, "line %d: %s", __LINE__, tsk_strerror(val));                    \
        exit(EXIT_FAILURE);                                                             \
    }


/*Claim Functions defined in this file*/
void input1(FILE *infile);
void input2(FILE *infile);


extern int Build(long *TYPE, long *type2node, long *ntype, const double THETA, double *logw, tsk_table_collection_t *tables, long *num_nodes, long *num_rec, long *num_mut, long *num_coal, long* num_nonrec, long* mut_by_site, long* mut_site, long* mut_state, long* mut_node, double *pHtau);


extern void FDupdate(long* TYPE, long* type2node, const double THETA, long k, long* ntype, long* nbranch, long* n, double* m, double* logw, tsk_table_collection_t* tables, int* rec_num, int* mut_num, int* coal_num, int* nonrecc, long* mut_by_site, long* mut_site, long* mut_state, long* mut_node, long* eligible);

extern double rexp(double rate);



/*Constants from input*/
long NRUN;
long Ne; /*effective population size*/
long NumLineage; /*the number of lineages to stop */
double THETA; /*mutation rate*/
//long Resample; /*Number of resamplings*/
double SEQLEN; /*sequence length*/

long M; /*number of genes observed at present*/
long L; /*number of loci*/
//long *ReTime; /*time to resample*/
long *TYPE_copy; /*store observed haplotypes*/
long ntype_copy; /*number of distinct types*/
double *positions; /*genetic position for each site */

/*constants related to resampling*/
long resampling;
long *res_time; /*when to check the resampling threshold*/
double *ESS_threshold;
double *Entropy_threshold;


/*use a structure to store mutations*/
struct mutations
{
    long site_id;
    long state;
    long node_id;
};

static void state_to_string(long state, char *buffer, size_t buffer_len)
{
    snprintf(buffer, buffer_len, "%ld", state);
}
int cmp(const void* a, const void* b) {
    struct mutations* x = (struct mutations*)a;
    struct mutations* y = (struct mutations*)b;
    if ((*x).site_id < (*y).site_id) return 1;
    else return 0;
}


int main(int argc, char *argv[]){
    // Seed the random number generator with the current time 
    srand(time(NULL));
    /*declare the vatiables here*/
    FILE *in1, *in2, *out, *seed_file, *seed_file_w, *test;

    double logw, *wt, wt_sum=0.0, *logwt, *unnorm_wt, *PHtau, pHtau;
    long *TYPE, *type2node, *type2node_copy, *ntype, *mut_by_site;
    long *n_nodes, *n_rec, *n_mut, *n_nonrec, *n_coal, num_rec, num_mut, num_coal, num_nodes, num_nonrec, mut;
    char filename[200];
    char derived_state[32];

    long i,j, l, k, coal;
    int id_sort;
    int flag, ret;
    tsk_table_collection_t tables, tables_copy;
    tsk_table_collection_t *all_tables;
	//tsk_treeseq_t ts;
    long* mut_by_site_copy;
    long* mut_site, * mut_state, * mut_node;

    
    if(argc>1) fopen_s(&in1, argv[1],"r");
    else fopen_s(&in1, "C:\\Users\\blabl\\Desktop\\git repos\\Sampling_with_Resampling\\InputFiles\\infile1","r"); 
    
    input1(in1);
    fclose(in1);
    /*read input1*/
    
    
    if(argc>2) fopen_s(&in2, argv[2],"r");
    else fopen_s(&in2, "C:\\Users\\blabl\\Desktop\\git repos\\Sampling_with_Resampling\\InputFiles\\infile2","r");
    
    input2(in2);
    fclose(in2);
    /*read input2*/
    
	/*scaled mutation rate*/
    THETA = 4.0 * Ne * THETA * SEQLEN;


    TYPE = long_vec_init(((L+1)*TYPE_MAX)*NRUN);  
    type2node = long_vec_init(NODE_MAX);
    type2node_copy = long_vec_init(NODE_MAX);
    ntype = long_vec_init(NRUN);
    wt = dou_vec_init(NRUN);
    unnorm_wt = dou_vec_init(NRUN);
    PHtau = dou_vec_init(NRUN);
	n_nodes = long_vec_init(NRUN);
	n_rec = long_vec_init(NRUN);
	n_mut = long_vec_init(NRUN);
    n_nonrec = long_vec_init(NRUN);
	n_coal = long_vec_init(NRUN);
	logwt = dou_vec_init(NRUN);
    ret = tsk_table_collection_init(&tables, 0);
    ret += tsk_table_collection_init(&tables_copy, 0);
    if(ret<0) fprintf(stderr, "Initialise table collection failed.\n");
    all_tables = (tsk_table_collection_t *) calloc((size_t) NRUN, sizeof(tsk_table_collection_t));
    if (all_tables == NULL) {
        fprintf(stderr, "Unable to allocate tree collection storage.\n");
        return EXIT_FAILURE;
    }
    for (i = 0; i < NRUN; i++) {
        ret = tsk_table_collection_init(&all_tables[i], 0);
        if (ret < 0) {
            fprintf(stderr, "Initialise stored tree collection %ld failed.\n", i);
            return EXIT_FAILURE;
        }
    }
    mut_by_site_copy = long_vec_init(L);
    mut_by_site = long_vec_init(L);
    mut_site = long_vec_init(MUT_MAX);
    mut_state = long_vec_init(MUT_MAX);
    mut_node = long_vec_init(MUT_MAX);
    struct mutations mut_record[MUT_MAX];

    _mkdir("output");

   

    /* initialise sequence length and sites table */
    tables_copy.sequence_length = SEQLEN;
    for(i=0; i<L; i++){
        ret = tsk_site_table_add_row(&tables_copy.sites, positions[i], 0, 0, NULL, 0);
        if(ret<0) fprintf(stderr, "Initialise position %f failed.\n", positions[i]);
    }

    /* initialize sample nodes, flag=1 */
    for(i=0; i<M; i++){
        ret = tsk_node_table_add_row(&tables_copy.nodes, 1, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
        /*ret is the ID of the newly added node on success */
        if(ret<0) fprintf(stderr, "Initialise sample %ld failed.\n", i);
        type2node_copy[i] = ret;
    }

    /* initialize mutation number on each site*/
    for (i = 0; i < L; i++) {
        mut_by_site_copy[i] = 0;
        for (j = 0; j < ntype_copy; j++) {
            if(*(TYPE_copy + (L + 1) * j + i) == 1) mut_by_site_copy[i] += *(TYPE_copy + (L + 1) * j + L);
        }
    }

    /*Start the for loop from present to the stopping time*/
    for (coal = 0; coal < (M-NumLineage); coal++) {
        fprintf(stderr, "\n--------------------------------------------------------------");
        if (coal == 0) {
            fprintf(stderr, "\n From present to the 1st coalescence:\n");
        }
        else {
            fprintf(stderr, "\n From %ld coalescence to %ld coalescence:\n", coal, coal + 1);
        }
        /*for loop over the number of particles*/
        for (i = 0; i < NRUN; i++) {
            ntype[i] = ntype_copy;
            ret = tsk_table_collection_copy(&tables_copy, &tables, TSK_NO_INIT);
            if (ret < 0) fprintf(stderr, "Copy table collection failed.\n");

            for (j = 0; j < ((L + 1) * ntype[i]); j++) *(TYPE + j) = *(TYPE_copy + j);
            for (j = 0; j < M; j++) *(type2node + j) = *(type2node_copy + j);
            for (j = M; j < NODE_MAX; j++) *(type2node + j) = -1;
            for (j = 0; j < L; j++) *(mut_by_site + j) = *(mut_by_site_copy + j);


            flag = Build(TYPE, type2node, ntype, THETA, &logw, &tables, &num_nodes, &num_rec, &num_mut, &num_coal, &num_nonrec, mut_by_site, mut_site, mut_state, mut_node, &pHtau);           

            if (flag == 0) {
                /*sort mutation site id and record mutations.
                Mutations must be provided in non-decreasing site order and non-increasing time order within each site*/
                int* indices2 = sort_and_return_indices2((int*)mut_site, (int)num_mut);

                for (j = 0; j < num_mut; j++) {
                    id_sort = indices2[j];
                    mut_record[j].site_id = mut_site[j];
                    mut_record[j].state = mut_state[id_sort];
                    mut_record[j].node_id = mut_node[id_sort];
                }
                free(indices2);

                for (j = 0; j < num_mut; j++) {
                    state_to_string(mut_record[j].state, derived_state, sizeof(derived_state));
                    mut = tsk_mutation_table_add_row(&tables.mutations, mut_record[j].site_id, mut_record[j].node_id, -1, TSK_UNKNOWN_TIME, derived_state, 0, NULL, 0);
                    if (mut < 0) fprintf(stderr, "Adding mutation failed.\n");
                }

                /* Table collection must be indexed*/
                ret = tsk_table_collection_build_index(&tables, 0);
                if (ret < 0) fprintf(stderr, "Building indices for table collection failed.\n");

                /*.trees output */
                /*.trees are efficient and for analyzing; .txt are inefficient and for viewing only*/
                if(argc>3){
                    snprintf(filename, 200, "%s_%ld.trees", argv[3], i+1);
                    ret = tsk_table_collection_dump(&tables, filename, 0);
                    if (ret < 0) fprintf(stderr, "Dumping .trees failed!\n");
                }
                else{
                    snprintf(filename, 200, "output\\tree_%ld.trees", i+1);
                    ret = tsk_table_collection_dump(&tables, filename, 0);
                    if (ret < 0) fprintf(stderr, "Dumping .trees failed!\n");
                }

                if (logw < -312) {
                    fprintf(stderr, "Unnormalized weight too small. logw = %lf.\n", logw);
                }
                else {
                    fprintf(stderr, "Unnormalized weight logw = %lf.\nZ", logw);
                }
                logwt[i] = logw;
                wt[i] = exp(logw);
                n_nodes[i] = num_nodes;
                n_rec[i] = num_rec;
                n_mut[i] = num_mut;
                n_nonrec[i] = num_nonrec;
                n_coal[i] = num_coal;
                PHtau[i] = pHtau;

                fprintf(stderr, "Building ARG succeeds.\n");

                ret = tsk_table_collection_copy(&tables, &all_tables[i], TSK_NO_INIT);
                if (ret < 0) {
                    fprintf(stderr, "Copy tree collection for output %ld failed.\n", i);
                }

                //write the trees in .txt
                //snprintf(filename, 200, "C:\\Users\\blabl\\Dropbox\\output\\out_%d.txt", i + 1);
                //fopen_s(&test, filename, "w");
                //tsk_table_collection_print_state(&tables, test);
                //fclose(test);


            }
            else {
                printf("Number of distinct haplotypes too large, ignoring trees\n");
                i--;
            }

        }


    }
    for(i=0; i<NRUN; i++){
		*ntype = ntype_copy;
        ret = tsk_table_collection_copy(&tables_copy, &tables, TSK_NO_INIT);
        if(ret<0) fprintf(stderr, "Copy table collection failed.\n");

        for(j=0; j<((L+1)*(*ntype)); j++) *(TYPE+j) = *(TYPE_copy+j);
        for(j=0; j<M; j++) *(type2node+j) = *(type2node_copy+j);
        for(j=M; j<NODE_MAX; j++) *(type2node+j) = -1;
        for (j = 0; j < L; j++) *(mut_by_site + j) = *(mut_by_site_copy + j);

        flag = Build(TYPE, type2node, ntype, THETA, &logw, &tables, &num_nodes, &num_rec, &num_mut, &num_coal, &num_nonrec, mut_by_site, mut_site, mut_state, mut_node, &pHtau);
		
        if(flag==0){
            /*sort mutation site id and record mutations.
            Mutations must be provided in non-decreasing site order and non-increasing time order within each site*/            
            int* indices2 = sort_and_return_indices2((int*)mut_site, (int)num_mut);
 
            for (j = 0; j < num_mut; j++) {
                id_sort = indices2[j];
                mut_record[j].site_id = mut_site[j];
                mut_record[j].state = mut_state[id_sort];
                mut_record[j].node_id = mut_node[id_sort];
            }  
            free(indices2);

            for (j = 0; j < num_mut; j++) {
                state_to_string(mut_record[j].state, derived_state, sizeof(derived_state));
                mut = tsk_mutation_table_add_row(&tables.mutations, mut_record[j].site_id, mut_record[j].node_id, -1, TSK_UNKNOWN_TIME, derived_state, 0, NULL, 0);
                if (mut < 0) fprintf(stderr, "Adding mutation failed.\n");
            }

            /* Table collection must be indexed*/
            ret = tsk_table_collection_build_index(&tables, 0);
            if (ret < 0) fprintf(stderr, "Building indices for table collection failed.\n");

            /*.trees output */
			/*.trees are efficient and for analyzing; .txt are inefficient and for viewing only*/	
            if(argc>3){
                snprintf(filename, 200, "%s_%ld.trees", argv[3], i+1);
                ret = tsk_table_collection_dump(&tables, filename, 0);
                if (ret < 0) fprintf(stderr, "Dumping .trees failed!\n");
            } 
            else{
                snprintf(filename, 200, "output\\tree_%ld.trees", i+1);
                ret = tsk_table_collection_dump(&tables, filename, 0);
                if (ret < 0) fprintf(stderr, "Dumping .trees failed!\n");
            }

            if (logw < -312) {
                fprintf(stderr, "Unnormalized weight too small. logw = %lf.\n", logw);
            }
            else {
                fprintf(stderr, "Unnormalized weight logw = %lf.\n", logw);
            }
			logwt[i] = logw;
            wt[i] = exp(logw);
			n_nodes[i] = num_nodes;
			n_rec[i] = num_rec;
			n_mut[i] = num_mut;
            n_nonrec[i] = num_nonrec;
			n_coal[i] = num_coal;
            PHtau[i] = pHtau;
			
			fprintf(stderr, "Building ARG succeeds.\n");

			//write the trees in .txt
            //snprintf(filename, 200, "C:\\Users\\blabl\\Dropbox\\output\\out_%d.txt", i + 1);
			//fopen_s(&test, filename, "w");
			//tsk_table_collection_print_state(&tables, test);
			//fclose(test);
            

        } 
        else{
            printf("Number of distinct haplotypes too large, ignoring trees\n");
            i--;
        }		


        ret = tsk_table_collection_clear(&tables, TSK_CLEAR_METADATA_SCHEMAS);
        if(ret<0) fprintf(stderr, "Clear tree collection failed.\n");

        if(i%1000 ==999) fprintf(stderr, "End of run %d \n", (int)i+1);

    }

	
    /*calc the normalized SIS weights*/
	for (i = 0; i < NRUN; i++) {
        unnorm_wt[i] = wt[i];
		wt_sum += wt[i];
	}
	for (i = 0; i < NRUN; i++) {
		wt[i] /= wt_sum;
	}
	

    /*The stopping time is reached*/

    /*output the weights*/
    if(argc>3){
        snprintf(filename, 200, "%s_wt.txt", argv[3]);
        fopen_s(&out, filename, "w");
    } 
    else{
        snprintf(filename, 200, "output\\test_wt.txt");
        fopen_s(&out, filename, "w");
    }
    if(out){
		fprintf(out, "run \t normalized_weight \t log(weight) \t num_nodes \t num_recombinations \t num_mutations \t num_nonrecurrent \t num_coalescence \t unnormalized_weight \t P(Htau)\n");
		for (i = 0; i < NRUN; i++) {
			fprintf(out, "%ld \t %.15lf \t %lf \t %ld \t %ld \t %ld \t %ld  \t %ld \t %.15lf \t %.15lf \n", i + 1, wt[i], logwt[i], n_nodes[i], n_rec[i], n_mut[i], n_nonrec[i], n_coal[i], unnorm_wt[i], PHtau[i]);
		}
		
		fclose(out);
    }
	else {
		for (i = 0; i < NRUN; i++) fprintf(stderr, "Cannot open the weight file.");
	}

    /* Resample the trees with replacement using the original weights. */
    {
        long *resampled_index;
        double *cdf;
        double max_logwt;
        double sum_prob;
        double equal_weight;
        FILE *resampled_out;

        resampled_index = long_vec_init(NRUN);
        cdf = dou_vec_init(NRUN);
        if (resampled_index != NULL && cdf != NULL) {
            max_logwt = logwt[0];
            for (i = 1; i < NRUN; i++) {
                if (logwt[i] > max_logwt) {
                    max_logwt = logwt[i];
                }
            }

            sum_prob = 0.0;
            for (i = 0; i < NRUN; i++) {
                cdf[i] = exp(logwt[i] - max_logwt);
                sum_prob += cdf[i];
            }
            for (i = 0; i < NRUN; i++) {
                cdf[i] /= sum_prob;
                if (i > 0) {
                    cdf[i] += cdf[i - 1];
                }
            }

            snprintf(filename, 200, "output\\resampled_samples.txt");
            fopen_s(&resampled_out, filename, "w");
            if (resampled_out) {
                fprintf(resampled_out, "new_run\tselected_tree\toriginal_log_weight\tequal_weight\n");
            }

            equal_weight = 1.0 / (double) NRUN;
            for (i = 0; i < NRUN; i++) {
                double u = runif();
                long chosen = 0;

                while (chosen < NRUN - 1 && u > cdf[chosen]) {
                    chosen++;
                }

                resampled_index[i] = chosen;
                if (resampled_out) {
                    fprintf(resampled_out, "%ld\t%ld\t%.15lf\t%.15lf\n", i + 1, chosen + 1, logwt[chosen], equal_weight);
                }
            }

            if (resampled_out) {
                fclose(resampled_out);
            }

            /* Dump each resampled tree as a .trees file after sorting and indexing. */
            for (i = 0; i < NRUN; i++) {
                long chosen = resampled_index[i];
                tsk_table_collection_t resampled_tables;

                ret = tsk_table_collection_init(&resampled_tables, 0);
                if (ret < 0) {
                    fprintf(stderr, "Initialise resampled tree collection %ld failed.\n", i);
                    continue;
                }

                ret = tsk_table_collection_copy(&all_tables[chosen], &resampled_tables, TSK_NO_INIT);
                if (ret < 0) {
                    fprintf(stderr, "Copy selected tree %ld for resampling failed.\n", chosen);
                    tsk_table_collection_free(&resampled_tables);
                    continue;
                }

                ret = tsk_table_collection_sort(&resampled_tables, NULL, 0);
                if (ret < 0) {
                    fprintf(stderr, "Sort selected tree %ld for resampling failed.\n", chosen);
                    tsk_table_collection_free(&resampled_tables);
                    continue;
                }

                ret = tsk_table_collection_build_index(&resampled_tables, 0);
                if (ret < 0) {
                    fprintf(stderr, "Build index for selected tree %ld failed.\n", chosen);
                    tsk_table_collection_free(&resampled_tables);
                    continue;
                }

                snprintf(filename, 200, "output\\resampled_tree_%ld.trees", i + 1);
                ret = tsk_table_collection_dump(&resampled_tables, filename, 0);
                if (ret < 0) {
                    fprintf(stderr, "Dump resampled tree %ld failed.\n", i);
                }

                tsk_table_collection_free(&resampled_tables);
            }

            for (i = 0; i < NRUN; i++) {
                wt[i] = equal_weight;
                logwt[i] = log(equal_weight);
            }
        }

        free(resampled_index);
        free(cdf);
    }
    

    
    /*change .seed file*/
    fopen_s(&seed_file_w, ".seed", "w");
    seed_put(seed_file_w);
    fclose(seed_file_w);

    
    /*free the variables*/
    ret = tsk_table_collection_free(&tables);
    ret = tsk_table_collection_free(&tables_copy);
    free(TYPE);
    free(TYPE_copy);
    //free(ReTime);
    free(type2node);
    free(type2node_copy);
    free(wt);
    free(logwt);
    free(positions);

    free(n_nodes);
    free(n_rec);
    free(n_mut);
    free(n_nonrec);
    free(n_coal);
    free(mut_by_site_copy);
    free(mut_by_site);

    free(res_time);
    free(ESS_threshold);
    free(Entropy_threshold);
    for (i = 0; i < NRUN; i++) {
        tsk_table_collection_free(&all_tables[i]);
    }
    free(all_tables);
  

    return(0);

}

void input1(FILE* infile)
{
    char string[200];
    char iden[2];
    short lines, count;
    short flag = 1;
    short i, j;

    NRUN = -1; /*number of runs*/
    Ne = -1;
    //Resample = -1; /*number of resamplings*/
    THETA = -1; /*mutation rate*/
    SEQLEN = -1;
    NumLineage = -1;
    resampling = -1;
    lines = 0;
    count = 0;

    fprintf(stderr, "\n Data in file 1: \n");

    while (flag) {
        if (fgets(string, sizeof(string), infile) == NULL) {
            fprintf(stderr, "Error in file \n");
            abort();
        }
        else {
            i = 0;
            while (string[i] == ' ' || string[i] == '\t') i++;
            iden[0] = string[i];
            switch (iden[0]) {
            case EOF:
                flag = 0;
                break;
            case '#':
                flag = 0;
                break;
            case '\n':
                break;

            case 'N':
                if (NRUN != -1) {
                    fprintf(stderr, "number of runs redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf_s(string + i, "%ld", &NRUN);
                    fprintf(stderr, "number of runs %ld \n", NRUN);
                    break;
                }


            case 'E':
                if (Ne != -1) {
                    fprintf(stderr, "Effective population size redefined on line %d \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf_s(string + i, "%ld", &Ne);
                    fprintf(stderr, "effective population size %ld \n", Ne);
                    break;
                }


            case'M':
                if (THETA != -1) {
                    fprintf(stderr, "Mutation rate redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf_s(string + i, "%lf", &THETA);
                    fprintf(stderr, "mutation rate %3.9lf per generation per base pair\n", THETA);
                    break;
                }

            case'L':
                if (SEQLEN != -1) {
                    fprintf(stderr, "Sequence length redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf_s(string + i, "%lf", &SEQLEN);
                    fprintf(stderr, "sequence length %lf \n", SEQLEN);
                    break;
                }

            case'S':
                if (NumLineage != -1) {
                    fprintf(stderr, "Number of lineages redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf_s(string + i, "%ld", &NumLineage);
                    fprintf(stderr, "number of lineages %ld \n", NumLineage);
                    break;
                }


            case'R':
                if (resampling != -1) {
                    fprintf(stderr, "Number of resamplings redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf_s(string + i, "%ld", &resampling);
                    fprintf(stderr, "Number of resamplings is %ld \n", resampling);
                    break;
                }

            case'T':
                if (resampling == -1) {
                    fprintf(stderr, "Number of resamplings not defined yet. Ignore definition of resampling times on line %d\n", lines);
                    break;
                }
                else {
                    if ((count & 0x1) == 0x1) {
                        fprintf(stderr, "Resampling times redefined on line %d. Ignored \n", lines);
                        break;
                    }
                    res_time = long_vec_init(resampling);
                    fprintf(stderr, "Resampling times: ");
                    for (j = 0; j < resampling; j++) {
                        if (fscanf_s(infile, "%ld", &res_time[j]) != 1) {
                            fprintf(stderr, "Error in resampling times for %d th entry \n", j + 1);
                            abort();
                        }
                        else {
                            fprintf(stderr, "%ld ", res_time[j]);
                        }
                    }
                    fprintf(stderr, "\n");
                    count += 1;
                    break;
                }

            case '1':
                if (resampling == -1) {
                    fprintf(stderr, "Number of resamplings not defined yet. Ignore definition of resampling threshold on line %d\n", lines);
                    break;
                }
                else {
                    if ((count & 0x2) == 0x2) {
                        fprintf(stderr, "Resampling threshold for ESS redefined on line %d. Ignored \n", lines);
                        break;
                    }
                    ESS_threshold = dou_vec_init(resampling);
                    fprintf(stderr, "Resampling threshold for ESS:\n");
                    for (j = 0; j < resampling; j++) {
                        if (fscanf_s(infile, "%lf", &ESS_threshold[j]) != 1) {
                            fprintf(stderr, "Error in 1st resampling threshold for %d th entry \n", j + 1);
                            abort();
                        }
                        else {
                            fprintf(stderr, "%lf ", ESS_threshold[j]);
                        }
                    }
                    fprintf(stderr, "\n");
                    count += 2;
                    break;
                }

            case '2':
                if (resampling == -1) {
                    fprintf(stderr, "Number of resamplings not defined yet. Ignore definition of resampling threshold on line %d\n", lines);
                    break;
                }
                else {
                    if ((count & 0x4) == 0x4) {
                        fprintf(stderr, "Resampling threshold for entropy redefined on line %d. Ignored \n", lines);
                        break;
                    }
                    Entropy_threshold = dou_vec_init(resampling);
                    fprintf(stderr, "Resampling threshold for entropy:\n");
                    for (j = 0; j < resampling; j++) {
                        if (fscanf_s(infile, "%lf", &Entropy_threshold[j]) != 1) {
                            fprintf(stderr, "Error in 2nd resampling threshold for %d th entry \n", j + 1);
                            abort();
                        }
                        else {
                            fprintf(stderr, "%lf ", Entropy_threshold[j]);
                        }
                    }
                    fprintf(stderr, "\n");
                    count += 4;
                    break;
                }
                


            default:
                fprintf(stderr, "Non-standard input line %d in file 1. Ignored", lines);
                break;
            }
            lines++;
        }
    }

    flag = 1;

    if (NRUN <= 0) {
        fprintf(stderr, "Number of runs <= 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Number of runs = %ld checked\n", NRUN);

    if (Ne <= 0) {
        fprintf(stderr, "Effective population size <= 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Effective population size = %ld checked\n", Ne);

    if (THETA < 0) {
        fprintf(stderr, "Theta value < 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Theta value = %3.9lf checked\n", THETA);

    if (resampling < 0) {
        fprintf(stderr, "Numbjer of resamplings < 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Number of resamplings = %ld checked\n", resampling);

    if (SEQLEN < 0) {
        fprintf(stderr, "Sequence length < 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Sequence length = %lf checked\n", SEQLEN);

    if (NumLineage < 0) {
        fprintf(stderr, "Number of lineages < 0\n");
        flag = 0;
    }
    else fprintf(stderr, "Number of lineages = %ld to stop\n", NumLineage);
    if (flag == 0) abort();
}

void input2(FILE* infile)
{
    char string[200];
    char iden[2];
    short lines, flag, count;
    long i, j;
    long temp;
    double tot;

    L = 0; /*number of loci*/
    M = 0; /*number of genes*/
    flag = 1;
    lines = 0;
    count = 0;
    ntype_copy = 0;

    fprintf(stderr, "\n Data in file 2: \n");

    while (flag) {
        if (fgets(string, sizeof(string), infile) == NULL) {
            fprintf(stderr, "Error in file \n");
            abort();
        }
        else {
            i = 0;
            while (string[i] == ' ' || string[i] == '\t') i++;
            iden[0] = string[i];
            switch (iden[0]) {
            case EOF:
                flag = 0;
                break;
            case '#':
                flag = 0;
                break;
            case '\n':
                break;

            case 'L':
                if (L > 0) {
                    fprintf(stderr, "number of loci redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf_s(string + i, "%ld", &L);
                    fprintf(stderr, "number of loci %ld \n", L);
                    break;
                }

            case 'G':
                if (M > 0) {
                    fprintf(stderr, "number of genes redefined on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf_s(string + i, "%ld", &M);
                    fprintf(stderr, "number of genes %ld \n", M);
                    break;
                }

            case 'D':
                if (ntype_copy > 0) {
                    fprintf(stderr, "number of distinct haplotypes on line %d. Ignored \n", lines);
                    break;
                }
                else {
                    while (string[i] != '=') i++;
                    i++;
                    sscanf_s(string + i, "%ld", &ntype_copy);
                    fprintf(stderr, "number of distinct haplotypes %ld \n", ntype_copy);
                    break;
                }


            case 'H':
                if (L == 0 || ntype_copy == 0) {
                    fprintf(stderr, "Number of loci or types not defined. Ignored line %d\n", lines);
                    break;
                }
                else {
                    if ((count & 0x1) == 0x1) {
                        fprintf(stderr, "Haplotypes redefined on line %d. Ignored \n", lines);
                        break;
                    }
                    TYPE_copy = long_vec_init((L + 1) * ntype_copy);
                    fprintf(stderr, "Haplotypes:\n");
                    for (j = 0; j < (L + 1) * ntype_copy; j++) {
                        if (fscanf_s(infile, "%ld", &TYPE_copy[j]) != 1) {
                            fprintf(stderr, "Error in haplotypes for %ld th entry \n", j + 1);
                            abort();
                        }
                        else {
                            if (j % (L + 1) != L) fprintf(stderr, "%ld ", TYPE_copy[j]);
                            else fprintf(stderr, "%ld\n", TYPE_copy[j]);
                        }
                    }
                    fprintf(stderr, "\n");
                    count += 1;
                    break;
                }



            case 'P':
                if (L == 0) {
                    fprintf(stderr, "Genetic positions defined before number of loci on line %d - ignore \n", lines);
                    break;
                }
                if ((count & 0x8) == 0x8) {
                    fprintf(stderr, "Genetic positions redefined on line %d \n", lines);
                    break;
                }
                positions = dou_vec_init(L);
                for (i = 0; i < L; i++) {
                    if (fscanf_s(infile, "%lf", &positions[i]) != 1) {
                        fprintf(stderr, "Error in inputting genetic positions on locus %ld \n", i + 1);
                        abort();
                    }
                }
                count += 8;
                break;

            default:
                fprintf(stderr, "Non-standard input line %d in file 2. Ignored \n", lines);
                break;
            }
            lines++;
        }
    }
    flag = 1;
    if (L <= 1) {
        fprintf(stderr, "Number of loci <= 1 or undefined \n");
        flag = 0;
    }
    if (ntype_copy <= 0) {
        fprintf(stderr, "Number of haplotypes <= 0 or undefined \n");
        flag = 0;
    }
    if (M <= 0) {
        fprintf(stderr, "Number of genes <= 0 or undefined \n");
        flag = 0;
    }
    if (flag == 0) abort();


    /*check positions */
    if ((count & 0x8) != 0x8) {
        fprintf(stderr, "Genetic positions not defined \n");
        flag = 0;
    }
    else {
        for (i = 0; i < L; i++) {
            if (positions[i] < 0) {
                fprintf(stderr, "Genetic position %ld should be non-negative \n", i + 1);
                flag = 0;
            }
        }
    }

    if (flag == 0) abort();

    /*check TYPE*/
    temp = 0;
    for (i = 0; i < ntype_copy; i++) {
        for (j = 0; j < (L); j++) {
            if (TYPE_copy[(L + 1) * i + j]<0 || TYPE_copy[(L + 1) * i + j]>K - 1) {
                fprintf(stderr, "Type for haplotype %ld at site %ld is not an allele \n", j + 1, i + 1);
                flag = 0;
            }
        }
        if (TYPE_copy[(L + 1) * i + L] <= 0) {
            fprintf(stderr, "Number of %ld th haplotype is less than or equal to 0 \n", i);
            flag = 0;
        }
        temp += TYPE_copy[(L + 1) * i + L];
    }
    if (temp != M) {
        fprintf(stderr, "Number of genes from haplotypesdoes not equal the number of genes specified \n");
        flag = 0;
    }

    if (flag == 1) fprintf(stderr, "All checks done\n");
    if (flag == 0) abort();
}












