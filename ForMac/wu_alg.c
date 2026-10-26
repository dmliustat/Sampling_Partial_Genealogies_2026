#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define MAX_NUM_ROWS 50
#define MAX_NUM_COLUMNS 100
#define MUTATION_PATH_TERMINATOR 0

void print_types(int (*)[MAX_NUM_COLUMNS]);
int *sort_and_return_indices(int[], int);
void reorganize_array(int original[][MAX_NUM_COLUMNS], int order[], int rows, int cols, int new[][MAX_NUM_COLUMNS]);


extern double THETA;
extern long L;
extern long Ne;
extern double SEQLEN;
extern long NumLineage;

enum event_type
{
    coal = 0,
    mutn = 1
};

// External global variables
int m = 0, n = 0;
int observed_types[MAX_NUM_ROWS][MAX_NUM_COLUMNS];
double theta;

// Local global variables
int final_node_ac_mps[MAX_NUM_ROWS][MAX_NUM_COLUMNS];
int site_nums[MAX_NUM_COLUMNS];
int num_coal_table[MAX_NUM_COLUMNS];

struct AC
{
    int (*mps)[MAX_NUM_COLUMNS];
    int mult[MAX_NUM_ROWS];
    int num_mps;
};

struct AE
{
    enum event_type type;
    int parent;
    int num;
};

struct AES
{
    struct AE aes[MAX_NUM_ROWS + MAX_NUM_COLUMNS];
    int num_AEs;
};

struct node
{
    struct AC ac;
    struct AES aes;
    double p;
};

struct key
{
    char *keyStr;
    int length;
};

struct stage
{
    struct node *nodes;
    struct key *keys;
    int num_nodes;
};

struct child_mutns
{
    int child_mutns[MAX_NUM_ROWS];
    int num_child_mutns;
};

void removeColumn(int array[MAX_NUM_ROWS][MAX_NUM_COLUMNS], int colIndex)
{
    int i, j;

    for (i = 0; i < m; i++)
    {
        for (j = colIndex; j < n - 1; j++)
        {
            array[i][j] = array[i][j + 1];
        }
    }
}

int (*trim_and_sort_types())[MAX_NUM_COLUMNS]
{
    int colsum;
    for (int j = 0; j < n - 1; j++) // skip the last column which is the multiplicities
    {
        colsum = 0;
        for (int i = 0; i < m; i++)
        {
            colsum += observed_types[i][j];
        }

        if (colsum == 0 || colsum == m)
        {
            removeColumn(observed_types, j);
            j--;
            n--;
        }
    }

    int *ncopies = (int *)malloc((n - 1) * sizeof(int));
    if (ncopies == NULL)
        return NULL;
    else
    {
        for (int j = 0; j < n - 1; j++)
        {
            ncopies[j] = 0;
            for (int i = 0; i < m; i++)
            {
                ncopies[j] += observed_types[i][j];
            }
        }
    }

    int *indices = sort_and_return_indices(ncopies, n - 1);
    free(ncopies);
    if (indices == NULL)
    {
        return NULL;
    }

    for (int j = 0; j < n - 1; j++)
    {
        site_nums[j] = indices[j] + 1;
    }

    int(*ordered)[MAX_NUM_COLUMNS] = malloc(MAX_NUM_ROWS * sizeof(int[MAX_NUM_COLUMNS]));
    if (ordered == NULL)
    {
        free(indices);
        return NULL;
    }

    reorganize_array(observed_types, indices, m, n - 1, ordered);
    free(indices);

    return ordered;
}

struct node make_final_node()
{
    struct AC ac;
    ac.num_mps = m;

    for (int i = 0; i < m; i++)
    {
        ac.mult[i] = observed_types[i][n - 1];
    }

    // print_types(observed_types);
    int(*types)[MAX_NUM_COLUMNS] = trim_and_sort_types();
    if (types == NULL)
        exit(1);

    // print_types(types);

    int idx;
    for (int i = 0; i < m; i++)
    {
        idx = 0;
        int j;
        for (j = 0; j < n - 1; j++)
        {
            if (types[i][j] == 1)
            {
                final_node_ac_mps[i][idx++] = site_nums[j];
            }
        }
        final_node_ac_mps[i][idx] = MUTATION_PATH_TERMINATOR;
    }

    
    free(types);
    //print_types(final_node_ac_mps);
    ac.mps = final_node_ac_mps;

    struct AES aes;
    aes.num_AEs = 0;

    struct node n = {ac, aes, -1.0};
    return n;
}

int get_num_coal_mutn(int mutn, struct AC ac)
{
    bool set[MAX_NUM_COLUMNS] = {false};
    int cnt = 0;
    for (int i = 0; i < m; i++)
    {
        int ind = -1;
        for (int j = 0; j < n; j++)
        {
            if (ac.mps[i][j] == mutn)
            {
                ind = j;
                break;
            }
            else if (ac.mps[i][j] == MUTATION_PATH_TERMINATOR)
            {
                break;
            }
        }
        if (ind == -1) // not found
        {
            continue;
        }
        else if (ind == 0) // the first element in the array
        {
            if (!set[0])
            {
                set[0] = true;
                cnt++;
            }
        }
        else
        {
            if (!set[ac.mps[i][ind - 1]])
            {
                set[ac.mps[i][ind - 1]] = true;
                cnt++;
            }
        }
    }
    return cnt - 1;
}

void get_num_coal_table(struct node final_node)
{
    for (int j = 0; j < n; j++)
    {
        num_coal_table[j] = get_num_coal_mutn(j, final_node.ac);
        //printf("\n", num_coal_table[j]);
    }
    int *tip_mutns = (int *)malloc(m * sizeof(int));
    if (tip_mutns == NULL)
    {
        printf("tip_mutns malloc failed!\n");
        exit(1);
    }
    for (int i = 0; i < m; i++)
    {
        tip_mutns[i] = final_node.ac.mps[i][0];
    }
    for (int i = 0; i < m; i++)
    {
        num_coal_table[tip_mutns[i]] += final_node.ac.mult[i] - 1;
    }
    free(tip_mutns);
}

int get_num_coal(int mutn)
{
    return num_coal_table[mutn];
}

struct AES init_AEs()
{
    int ncoal = get_num_coal(0);
    struct AES aes;
    aes.aes[0].num = ncoal;
    aes.aes[0].parent = 0;
    aes.aes[0].type = 0;
    aes.num_AEs = 1;
    return aes;
}

struct AC init_AC()
{
    struct AC ac;
    int(*mps)[MAX_NUM_COLUMNS] = malloc(MAX_NUM_ROWS * sizeof(int[MAX_NUM_COLUMNS]));
    if (mps == NULL)
        exit(1);
    mps[0][0] = MUTATION_PATH_TERMINATOR;
    ac.mult[0] = 1;
    ac.num_mps = 1;
    ac.mps = mps;
    return ac;
}

struct node init_node()
{
    struct AC ac = init_AC();
    struct AES aes = init_AEs();
    double p = 1.0;
    struct node n = {ac, aes, p};
    return n;
}

struct stage *init_stage()
{
    struct stage *s = (struct stage *)malloc(sizeof(struct stage));
    if (s == NULL)
    {
        exit(1);
    }

    (*s).nodes = (struct node *)malloc(sizeof(struct node[1]));
    if ((*s).nodes == NULL)
        exit(1);
    (*s).nodes[0] = init_node();

    (*s).keys = (struct key *)malloc(sizeof(struct key[1]));
    if ((*s).keys == NULL)
    {
        free((*s).nodes[0].ac.mps);
        free((*s).nodes);
        exit(1);
    }

    (*s).keys[0].keyStr = (char *)malloc(sizeof(char[4]));
    if ((*s).keys[0].keyStr == NULL)
    {
        free((*s).nodes[0].ac.mps);
        free((*s).nodes);
        free((*s).keys);
        exit(1);
    }
    char str[] = "0:1";
    (*s).keys[0].keyStr[0] = '\0';
    strcat_s((*s).keys[0].keyStr, 4, str);
    (*s).keys[0].length = 3;
    (*s).num_nodes = 1;
    return s;
}

int find_path(int mutn, int (*mps)[MAX_NUM_COLUMNS], int num_mps)
{
    for (int i = 0; i < num_mps; i++)
    {
        if (mps[i][0] == mutn)
        {
            return i;
        }
    }
    // should always be able to find the path;
    exit(1);
}

struct AC make_desc_AC(struct AE ae, struct AC ac)
{
    int curr_path = find_path(ae.parent, ac.mps, ac.num_mps);
    // make a deep copy of ac
    struct AC new_ac;
    new_ac.mps = malloc(MAX_NUM_ROWS * sizeof(int[MAX_NUM_COLUMNS]));
    new_ac.num_mps = ac.num_mps;
    for (int i = 0; i < new_ac.num_mps; i++)
    {
        new_ac.mult[i] = ac.mult[i];
        for (int j = 0; j < n - 1; j++)
        {
            new_ac.mps[i][j] = ac.mps[i][j];
            if (new_ac.mps[i][j] == MUTATION_PATH_TERMINATOR)
                break;
        }
    }

    if (ae.type == 0) // coal
    {
        new_ac.mult[curr_path] += 1;
    }
    else // mutation
    {
        if (ac.mult[curr_path] == 1)
        {
            // shift all the items to the right by 1
            int tmp1 = ae.num, tmp2;
            for (int j = 0; j < n - 1; j++)
            {
                tmp2 = new_ac.mps[curr_path][j];
                new_ac.mps[curr_path][j] = tmp1;
                tmp1 = tmp2;
                if (tmp1 == MUTATION_PATH_TERMINATOR)
                {
                    new_ac.mps[curr_path][j + 1] = MUTATION_PATH_TERMINATOR;
                    break;
                }
            }
        }
        else
        {
            new_ac.mult[new_ac.num_mps] = 1;
            new_ac.mps[new_ac.num_mps][0] = ae.num;
            for (int j = 1; j < n - 1; j++)
            {
                new_ac.mps[new_ac.num_mps][j] = ac.mps[curr_path][j - 1];
                if (new_ac.mps[new_ac.num_mps][j] == MUTATION_PATH_TERMINATOR)
                    break;
            }
            new_ac.num_mps += 1;
            new_ac.mult[curr_path] -= 1;
        }
    }

    // print_types(new_ac.mps);
    return new_ac;
}

double weight(enum event_type type, int n, int mult, double theta)
{
    double wt;
    if (type == 0)
    {
        wt = mult / (n + theta);
    }
    else
    {
        wt = mult * theta / (n * ((double)n - 1 + theta));
    }
    return wt;
}

double event_weight(struct AE event, struct AC ac, double theta)
{
    int sum_multiplicities = 0;
    for (int i = 0; i < ac.num_mps; i++)
    {
        sum_multiplicities += ac.mult[i];
    }
    int mult = ac.mult[find_path(event.parent, ac.mps, ac.num_mps)];
    return weight(event.type, sum_multiplicities, mult, theta);
}

double calc_prog_contrib(int ae_index, struct node progen_node, double theta)
{
    double new_wt = event_weight(progen_node.aes.aes[ae_index], progen_node.ac, theta);
    double new_p = new_wt * progen_node.p;
    return new_p;
}

// Function to concatenate strings with a specified separator
char *str_c(const char *const *x, int n, const char *separator, int *ord)
{
    int total_length = 0;
    for (int i = 0; i < n; i++)
    {
        total_length += (int) strlen(x[i]);
    }
    total_length += (n - 1) * (int)strlen(separator) + 1;
    char *result = (char *)malloc(total_length * sizeof(char));
    if (result == NULL)
    {
        return NULL;
    }
    result[0] = '\0';
    for (int i = 0; i < n; i++)
    {
        int idx = i;
        if (ord != NULL)
        {
            idx = ord[i];
        }
        strcat_s(result, total_length, x[idx]);
        if (i < n - 1)
        {
            strcat_s(result, total_length, separator);
        }
    }
    return result;
}

// Function to convert an integer array to an array of strings
char **int_array_to_str_array(int *arr, int n)
{
    char **str_arr = (char **)malloc(n * sizeof(char *));
    if (str_arr == NULL)
        return NULL;
    for (int i = 0; i < n; i++)
    {
        str_arr[i] = (char *)malloc(6 * sizeof(char)); // Assuming int won't exceed 5 digits + 1 for null terminator
        if (str_arr[i] == NULL)
            return NULL;
        snprintf(str_arr[i], 6, "%d", arr[i]);
        if (arr[i] == MUTATION_PATH_TERMINATOR)
            break;
    }
    return str_arr;
}

int find_length(int *arr, int n, int terminator)
{
    for (int i = 0; i < n; i++)
    {
        if (arr[i] == terminator)
            return i + 1;
    }
    return n;
}

struct key make_key(struct AC ac)
{
    char *separator_comma = ",";
    char *separator_semicolon = ";";
    char *separator_colon = ":";

    int num_rows = ac.num_mps;
    // Convert integer arrays in mps to arrays of strings
    char **mpaths = (char **)malloc(num_rows * sizeof(char *));
    if (mpaths == NULL)
        exit(1);

    for (int i = 0; i < num_rows; i++)
    {
        int length = find_length(ac.mps[i], n, MUTATION_PATH_TERMINATOR);
        char **mps_str = int_array_to_str_array(ac.mps[i], length);
        mpaths[i] = str_c((const char **)mps_str, length, separator_comma, NULL);

        // Free dynamically allocated memory for mps_str
        for (int j = 0; j < length; j++)
        {
            free(mps_str[j]);
        }
        free(mps_str);
    }

    // Sort and concatenate mutation paths
    int *ord = (int *)malloc(num_rows * sizeof(int));
    if (ord == NULL)
    {
        exit(1);
    }
    for (int i = 0; i < num_rows; i++)
    {
        ord[i] = i;
    }
    for (int i = 0; i < num_rows - 1; i++)
    {
        for (int j = i + 1; j < num_rows; j++)
        {
            if (strcmp(mpaths[ord[i]], mpaths[ord[j]]) > 0)
            {
                int temp = ord[i];
                ord[i] = ord[j];
                ord[j] = temp;
            }
        }
    }

    char *sorted_mpaths = str_c((const char **)mpaths, num_rows, separator_semicolon, ord);

    // Concatenate multiplicities
    char **mult_str = int_array_to_str_array(ac.mult, num_rows);
    char *sorted_mults = str_c((const char **)mult_str, num_rows, separator_comma, ord);

    // Concatenate everything
    char *result = str_c((const char *[]){sorted_mpaths, sorted_mults}, 2, separator_colon, NULL);

    // Free dynamically allocated memory
    for (int i = 0; i < num_rows; i++)
    {
        free(mpaths[i]);
        free(mult_str[i]);
    }
    free(mpaths);
    free(mult_str);
    free(ord);
    free(sorted_mpaths);
    free(sorted_mults);

    struct key ret = {result, (int)strlen(result)};
    return ret;
}

int find_key(struct key new_key, struct stage next_stage)
{
    int ind = -1;
    for (int i = 0; i < next_stage.num_nodes; i++)
    {
        if (next_stage.keys[i].length != new_key.length)
            continue;
        int j;
        for (j = 0; j < new_key.length; j++)
        {
            if (new_key.keyStr[j] != next_stage.keys[i].keyStr[j])
                break;
        }
        if (j == new_key.length)
        {
            ind = i;
            break;
        }
    }
    return ind;
}

void update_prob(struct stage *next_stage, int ind, double new_p)
{
    (*next_stage).nodes[ind].p += new_p;
}

struct child_mutns get_child_mutns(int mutn)
{
    struct child_mutns cm;
    cm.num_child_mutns = 0;
    bool set[MAX_NUM_COLUMNS] = {false};
    for (int i = 0; i < m; i++)
    {
        int ind = -1;
        for (int j = 0; j < n; j++)
        {
            if (final_node_ac_mps[i][j] == mutn)
            {
                ind = j;
                break;
            }
            else if (final_node_ac_mps[i][j] == MUTATION_PATH_TERMINATOR)
            {
                break;
            }
        }
        if (ind > 0)
        {
            if (!set[final_node_ac_mps[i][ind - 1]])
            {
                set[final_node_ac_mps[i][ind - 1]] = true;
                cm.child_mutns[cm.num_child_mutns++] = final_node_ac_mps[i][ind - 1];
            }
        }
    }
    return cm;
}

void add_events_after_mutn(struct AE target_event, struct AES *desc_AEs)
{
    int ncoal = get_num_coal(target_event.num);
    if (ncoal > 0)
    {
        struct AE ae = {0, target_event.num, ncoal};
        (*desc_AEs).aes[(*desc_AEs).num_AEs++] = ae;
    }
    else
    {
        struct child_mutns cm = get_child_mutns(target_event.num);
        if (cm.num_child_mutns > 0)
        {
            struct AE ae = {1, target_event.num, cm.child_mutns[0]};
            (*desc_AEs).aes[(*desc_AEs).num_AEs++] = ae;
        }
    }
}

bool is_first_coal(struct AE target_event)
{
    return target_event.num == (get_num_coal(target_event.parent) - 1);
}

void add_events_after_coal(struct AE target_event, struct AES *desc_AEs)
{
    bool is_first = is_first_coal(target_event);
    if (is_first)
    {
        struct child_mutns cm = get_child_mutns(target_event.parent);
        for (int i = 0; i < cm.num_child_mutns; i++)
        {
            struct AE ae = {1, target_event.parent, cm.child_mutns[i]};
            (*desc_AEs).aes[(*desc_AEs).num_AEs++] = ae;
        }
    }
    if (target_event.num >= 1)
    {
        (*desc_AEs).aes[(*desc_AEs).num_AEs++] = target_event;
    }
}

void add_events(struct AE target_event, struct AES *desc_AEs)
{
    if (target_event.type == 0) // coal
    {
        target_event.num -= 1;
        add_events_after_coal(target_event, desc_AEs);
    }
    else // mutn
    {
        add_events_after_mutn(target_event, desc_AEs);
    }
}

struct AES make_desc_AEs(struct AES aes, int i)
{
    struct AE target_event = {aes.aes[i].type, aes.aes[i].parent, aes.aes[i].num};
    // make a deep copy of aes first
    struct AES desc_AEs;
    desc_AEs.num_AEs = aes.num_AEs - 1;
    int ind = 0;
    for (int j = 0; j < aes.num_AEs; j++)
    {
        if (j == i)
            continue;
        desc_AEs.aes[ind].type = aes.aes[j].type;
        desc_AEs.aes[ind].parent = aes.aes[j].parent;
        desc_AEs.aes[ind].num = aes.aes[j].num;
        ind++;
    }
    add_events(target_event, &desc_AEs);
    return desc_AEs;
}

struct node make_desc_node(struct node progen_node, int ae_index, struct AC new_AC, double new_p)
{
    struct AES new_AEs = make_desc_AEs(progen_node.aes, ae_index);
    struct node n = {new_AC, new_AEs, new_p};
    return n;
}

void add_desc_node(struct node new_node, struct key new_key, struct stage *new_stage)
{
    (*new_stage).nodes[(*new_stage).num_nodes] = new_node;
    (*new_stage).keys[(*new_stage).num_nodes++] = new_key;
}

void update_next_stage(int ae_index, struct node progen_node, struct stage *next_stage)
{
    struct AC new_AC = make_desc_AC(progen_node.aes.aes[ae_index], progen_node.ac);
    double new_p = calc_prog_contrib(ae_index, progen_node, theta);

    struct key new_key = make_key(new_AC);
    // printf("new_p: %f, Key: %s\n", new_p, new_key.keyStr);
    int ind = find_key(new_key, *next_stage);
    if (ind != -1) // this AC is already in the next stage
    {
        update_prob(next_stage, ind, new_p);
        free(new_AC.mps);
        free(new_key.keyStr);
    }
    else // this AC isn't yet in next stage, so create a new node and add it to next stage
    {
        struct node new_node = make_desc_node(progen_node, ae_index, new_AC, new_p);
        add_desc_node(new_node, new_key, next_stage);
    }
}

bool has_active_coal(struct AE ae, int *aparents, int parents_cnt)
{
    for (int i = 0; i < parents_cnt; i++)
    {
        if (ae.parent == aparents[i])
        {
            return true;
        }
    }
    return false;
}

bool has_one_mutn_free(struct AE ae, struct AC ac)
{
    for (int i = 0; i < ac.num_mps; i++)
    {
        if (ac.mps[i][0] == ae.parent && ac.mult[i] == 1)
        {
            return true;
        }
    }
    return false;
}

bool inelig(struct AE ae, int *aparents, int parents_cnt, struct AC ac)
{
    return has_active_coal(ae, aparents, parents_cnt) && has_one_mutn_free(ae, ac);
}

bool iselig(int ae_index, struct node n)
{
    if (n.aes.aes[ae_index].type == 0) // coal
    {
        return true;
    }
    else // mutn
    {
        int *aparents = (int *)malloc(n.aes.num_AEs * sizeof(int));
        if (aparents == NULL)
            exit(1);

        int parents_cnt = 0;
        for (int i = 0; i < n.aes.num_AEs; i++)
        {
            struct AE ae = n.aes.aes[i];
            if (ae.type == 0) // coal
            {
                aparents[parents_cnt++] = ae.parent;
            }
        }
        bool ret = !inelig(n.aes.aes[ae_index], aparents, parents_cnt, n.ac);
        free(aparents);
        return ret;
    }
}

void clean_up(struct stage *curr_stage)
{
    // free dynamically allocated memory
    for (int i = 0; i < (*curr_stage).num_nodes; i++)
    {
        free((*curr_stage).keys[i].keyStr);
        free((*curr_stage).nodes[i].ac.mps);
    }
    free((*curr_stage).nodes);
    free((*curr_stage).keys);
    free(curr_stage);
}

struct stage *construct_next_stage(struct stage *curr_stage)
{
    struct stage *next_stage = (struct stage *)malloc(sizeof(struct stage));
    if (next_stage == NULL)
    {
        clean_up(curr_stage);
        exit(1);
    }
    (*next_stage).num_nodes = 0;
    (*next_stage).nodes = malloc(sizeof(struct node) * (MAX_NUM_ROWS * MAX_NUM_COLUMNS));
    (*next_stage).keys = malloc(sizeof(struct key) * (MAX_NUM_ROWS * MAX_NUM_COLUMNS));
    for (int i = 0; i < (*curr_stage).num_nodes; i++)
    {
        struct node progen_node = (*curr_stage).nodes[i];
        for (int j = 0; j < progen_node.aes.num_AEs; j++)
        {
            if (iselig(j, progen_node))
            {
                update_next_stage(j, progen_node, next_stage);
            }
        }
    }

    clean_up(curr_stage);
    return next_stage;
}

double get_exact_prob(double theta2)
{
    int sum_multiplicities = 0;
    for (int i = 0; i < m; i++)
    {
        sum_multiplicities += observed_types[i][n - 1];
    }
    if (sum_multiplicities != NumLineage) printf("Error! Number of lineages at the stopping time!");

    struct node final_node = make_final_node();

    get_num_coal_table(final_node);
   
    theta = theta2;

    struct stage *curr_stage = init_stage();

    int num_mutn = n - 1;

    for (int i = 0; i < num_mutn + sum_multiplicities - 1; i++)
    {
        //printf("================Stage %d================\n", i);
        curr_stage = construct_next_stage(curr_stage);
    }

    double ret = (*curr_stage).nodes[0].p;
    clean_up(curr_stage);
    return ret;
}

// utilities
void print_types(int (*types)[MAX_NUM_COLUMNS])
{
    printf("print types\n");
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%d ", types[i][j]);
        }
        printf("\n");
    }
}

void swap(int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

int partition(int arr[], int index[], int low, int high)
{
    int pivot = arr[high];
    int i = low - 1;
    for (int j = low; j < high; j++)
    {
        if (arr[j] <= pivot)
        {
            i++;
            swap(&arr[i], &arr[j]);
            swap(&index[i], &index[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    swap(&index[i + 1], &index[high]);
    return i + 1;
}

void quick_sort(int arr[], int index[], int low, int high)
{
    if (low < high)
    {
        int pivot = partition(arr, index, low, high);
        quick_sort(arr, index, low, pivot - 1);
        quick_sort(arr, index, pivot + 1, high);
    }
}

int *sort_and_return_indices(int arr[], int size)
{
    int *indices = (int *)malloc(size * sizeof(int));
    if (indices == NULL)
        return NULL;
    // Initialize indices array
    for (int i = 0; i < size; i++)
    {
        indices[i] = i;
    }
    // Sort the array and indices simultaneously
    quick_sort(arr, indices, 0, size - 1);
    return indices;
}

void reorganize_array(int original[][MAX_NUM_COLUMNS], int order[], int rows, int cols, int new[][MAX_NUM_COLUMNS])
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            new[i][j] = original[i][order[j]];
        }
    }
}

double calc_Wu(long* TYPE, long* ntype)
{
    int i,j;
    double theta;

    theta =THETA;


    /*m = number of distinct haplotypes, n =number of sites + 1 (multiplicity) */
    m = *ntype;
    n = L + 1;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
        observed_types[i][j] = *(TYPE + (L + 1) * i + j);
    }
    

    /*FILE *fptr;
    char* next_token;
    fopen_s(&fptr, "./input.txt", "r");

    char line[1000];

    if (fptr != NULL)
    {
        while (fgets(line, 1000, fptr))
        {
            n = 0;

            char *pToken = strtok_s(line, &line, ",", &next_token);
            while (pToken != NULL)
            {
                observed_types[m][n++] = atoi(pToken);
                pToken = strtok_s(NULL, &line, ",", &next_token);
            }
            m++;
        }
    }*/


    print_types(observed_types);
    double p = get_exact_prob(theta);
    //double p = 0.01;
    if (p <= 0.0 || p >= 1.0) {
        fprintf(stderr, "Wu alg probability = %lf is wrong! \n", p);
        exit(1);
    }

    printf("\nWu alg result: %.50f\n", p);
    return p;
}