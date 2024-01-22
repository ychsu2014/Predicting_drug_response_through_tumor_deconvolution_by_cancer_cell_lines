##### modify for multiple conditions
in_mut_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE/"
in_oncoKB_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/"
in_CGC = "/home/CBBI/hsuy1/projects/drugResponse/data/COSMIC/Census_allMon_Jan_30_21_08_13_2023.tsv"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/greedy_algorithm_filtering_conditions_and_selected_cell_lines_for_all_cancers.txt"

from itertools import chain, combinations
import os

fout = open(out_file, "w")
# one cell line per line, one conditin multiple lines
fout.write("cancer_type\tfiltering_condition\tselected_cell_line_num\tcell_line\n")

# oncoKB files
all_files = os.listdir(in_oncoKB_path)
oncoKB_files = []
for afile in all_files:
    if "_in_oncoKB_final.txt" in afile:
        oncoKB_files.append(afile)

# mutation files
#all_files = os.listdir(in_mut_path)
#mut_files = []
#for afile in all_files:
#    if "_TCGA_CCLE_compare.txt" in afile:
#        mut_files.append(afile)

# CGC genes
f = open(in_CGC)
lines = f.readlines()
lines = lines[1:]
CGC_list = []
for line in lines:
    cols = line.strip("\n").split("\t")
    CGC_list.append(cols[0])
CGC_list = list(set(CGC_list))

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

def addToUniverseAndSubsets(inMutLine, inCellLines, inUniverse, inSubsets):
    inUniverse.add(inMutLine)
    for cellLine in inCellLines:
        if cellLine in inSubsets.keys():
            mutSet = inSubsets[cellLine]
            mutSet.add(inMutLine)
            inSubsets[cellLine] = mutSet
        else:
            inSubsets[cellLine] = {inMutLine}
    return([inUniverse, inSubsets])

def createCostDict(inSubsets):
    outCost = {}
    for cellLine in inSubsets.keys():
        outCost[cellLine] = 1
    return(outCost)

for oncoKB_file in oncoKB_files:
    cancer_type = oncoKB_file.split("_")[0]
    # mutations that are targets of oncoKB drugs
    f = open(in_oncoKB_path + oncoKB_file)
    lines = f.readlines()
    keep_oncoKB_mut_lines = []
    keep_oncoKB_actionable_mut_lines = []
    for line in lines:
        cols = line.strip("\n").split("\t")
        drug = cols[len(cols)-1]
        mut_line = "\t".join(cols[:16]) + "\n"
        ##### all oncoKB mutations
        keep_oncoKB_mut_lines.append(mut_line)
        ##### only actionable oncoKB mutations
        if drug != "":
            keep_oncoKB_actionable_mut_lines.append(mut_line)
    f.close()
    # load file and get universe and subsets
    mut_file = cancer_type + "_TCGA_CCLE_compare.txt"
    f = open(in_mut_path + mut_file)
    lines = f.readlines()
    lines = lines[1:]
    universe_oncoKB_CGC = set()
    subsets_oncoKB_CGC = {}
    universe_oncoKB = set()
    subsets_oncoKB = {}
    universe_oncoKB_actionable = set()
    subsets_oncoKB_actionable = {}
    for line in lines:
        cols = line.strip("\n").split("\t")
        mut_line = "\t".join(cols[:16]) + "\n"
        cell_lines = cols[16:]
        gene = cols[6]
        ##### oncoKB + CGC
        if mut_line in keep_oncoKB_mut_lines or gene in CGC_list:
            [universe_oncoKB_CGC, subsets_oncoKB_CGC] = addToUniverseAndSubsets(mut_line, cell_lines, universe_oncoKB_CGC, subsets_oncoKB_CGC)
        ##### only oncoKB
        if mut_line in keep_oncoKB_mut_lines:
            [universe_oncoKB, subsets_oncoKB] = addToUniverseAndSubsets(mut_line, cell_lines, universe_oncoKB, subsets_oncoKB)
        ##### only oncoKB actinable mutations
        if mut_line in keep_oncoKB_actionable_mut_lines:
            [universe_oncoKB_actionable, subsets_oncoKB_actionable] = addToUniverseAndSubsets(mut_line, cell_lines, universe_oncoKB_actionable, subsets_oncoKB_actionable)
    f.close()
    # cost dicts
    costs_oncoKB_CGC = createCostDict(subsets_oncoKB_CGC)
    costs_oncoKB = createCostDict(subsets_oncoKB)
    costs_oncoKB_actionable = createCostDict(subsets_oncoKB_actionable)
    # cell line selection (greedy algorithm) and write to file
    print(cancer_type)
    ##### oncoKB + CGC 
    greedy_cover_oncoKB_CGC = greedy_set_cover(universe_oncoKB_CGC, subsets_oncoKB_CGC, costs_oncoKB_CGC)
    greedy_cost_oncoKB_CGC = sum(costs_oncoKB_CGC[s] for s in greedy_cover_oncoKB_CGC)    
    fout.write(cancer_type + "\toncoKB + CGC\t" + str(len(greedy_cover_oncoKB_CGC)) + "\t" + str(greedy_cover_oncoKB_CGC) + "\n")
    print("oncoKB + CGC")
    print('Greedy Set Cover: ' + str(greedy_cover_oncoKB_CGC))
    print('Cost = %s' % greedy_cost_oncoKB_CGC)
    ##### only oncoKB
    greedy_cover_oncoKB = greedy_set_cover(universe_oncoKB, subsets_oncoKB, costs_oncoKB)
    greedy_cost_oncoKB = sum(costs_oncoKB[s] for s in greedy_cover_oncoKB)
    fout.write(cancer_type + "\toncoKB\t" + str(len(greedy_cover_oncoKB)) + "\t" + str(greedy_cover_oncoKB) + "\n")
    print("oncoKB")
    print('Greedy Set Cover: ' + str(greedy_cover_oncoKB))
    print('Cost = %s' % greedy_cost_oncoKB)
    ##### only oncoKB actinable mutations
    greedy_cover_oncoKB_actionable = greedy_set_cover(universe_oncoKB_actionable, subsets_oncoKB_actionable, costs_oncoKB_actionable)
    greedy_cost_oncoKB_actionable = sum(costs_oncoKB_actionable[s] for s in greedy_cover_oncoKB_actionable)
    fout.write(cancer_type + "\toncoKB actionable target\t" + str(len(greedy_cover_oncoKB_actionable)) + "\t" + str(greedy_cover_oncoKB_actionable) + "\n")
    print("oncoKB actionable")
    print('Greedy Set Cover: ' + str(greedy_cover_oncoKB_actionable))
    print('Cost = %s' % greedy_cost_oncoKB_actionable)   

fout.close()
