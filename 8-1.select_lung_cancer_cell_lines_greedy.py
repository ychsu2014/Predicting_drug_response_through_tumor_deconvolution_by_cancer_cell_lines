##### modify for multiple conditions
in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE/LungCancer_TCGA_CCLE_compare.txt"
in_oncoKB = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/LungCancer_TCGA_CCLE_compare_in_oncoKB_final.txt"
in_CGC = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/COSMIC/Census_allMon_Jan_30_21_08_13_2023.tsv"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/selected_cell_lines_for_lung_cancer.txt"

from itertools import chain, combinations

# CGC genes
f = open(in_CGC)
lines = f.readlines()
lines = lines[1:]
CGC_list = []
for line in lines:
    cols = line.strip("\n").split("\t")
    CGC_list.append(cols[0])
CGC_list = list(set(CGC_list))

# mutations that are targets of oncoKB drugs
f = open(in_oncoKB)
lines = f.readlines()
keep_oncoKB_mut_lines = []
for line in lines:
    cols = line.strip("\n").split("\t")
    drug = cols[len(cols)-1]
    ##### all oncoKB mutations
    #if drug != "":
    #mut_line = "\t".join(cols[:16]) + "\n"
    #keep_oncoKB_mut_lines.append(mut_line)
    ##### only actionable oncoKB mutations
    if drug != "":
        mut_line = "\t".join(cols[:16]) + "\n"
        keep_oncoKB_mut_lines.append(mut_line)
    
f.close()

# load file and get universe and subsets
f = open(in_file)
lines = f.readlines()
lines = lines[1:]
universe = set()
subsets = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    mut_line = "\t".join(cols[:16]) + "\n"
    gene = cols[6]
    ##### oncoKB + CGC
    #if mut_line in keep_oncoKB_mut_lines or gene in CGC_list:
    ##### only oncoKB
    if mut_line in keep_oncoKB_mut_lines:
        universe.add(mut_line)
        cell_lines = cols[16:]
        for cell_line in cell_lines:
            if cell_line in subsets.keys():
                mut_set = subsets[cell_line]
                mut_set.add(mut_line)
                subsets[cell_line] = mut_set
            else:
                subsets[cell_line] = {mut_line}

# cost dict
costs = {}
for cell_line in subsets.keys():
    costs[cell_line] = 1
f.close()

def greedy_set_cover(universe, subsets, costs):
    """Approximate greedy algorithm for set-covering. Can be used on large
    ¡inputs - though not an optimal solution.

    Args:
        universe (list): Universe of elements
        subsets (dict): Subsets of U {S1:elements,S2:elements}
        costs (dict): Costs of each subset in S - {S1:cost, S2:cost...}
    """
    elements = set(e for s in subsets.keys() for e in subsets[s])
    # elements don't cover universe -> invalid input for set cover
    if elements != universe:
        return None

    # track elements of universe covered
    covered = set()
    cover_sets = []

    while covered != universe:
        min_cost_elem_ratio = float("inf")
        min_set = None
        # find set with minimum cost:elements_added ratio
        for s, elements in subsets.items():
            new_elements = len(elements - covered)
            # set may have same elements as already covered -> new_elements = 0
            # check to avoid division by 0 error
            if new_elements != 0:
                cost_elem_ratio = costs[s] / new_elements
                if cost_elem_ratio < min_cost_elem_ratio:
                    min_cost_elem_ratio = cost_elem_ratio
                    min_set = s
        cover_sets.append(min_set)
        # union
        covered |= subsets[min_set]
    return cover_sets

greedy_cover = greedy_set_cover(universe, subsets, costs)
greedy_cost = sum(costs[s] for s in greedy_cover)

fout = open(out_file, "w")
for cell_line in greedy_cover:
    fout.write(cell_line + "\n")
fout.close()

print('Greedy Set Cover:')
print(greedy_cover)
print('Cost = %s' % greedy_cost)
