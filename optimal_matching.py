import score_helices
import numpy as np
from scipy.optimize import linear_sum_assignment
import itertools

#retrieved from
#https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
def levenshtein(s1, s2):
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            # j+1 instead of j since previous_row and current_row
            # are one character longer than s2
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]

#calculates the edit distances between all of the helix pairs between the two input structures
def edit_distance_matrix(structure_a,structure_b,helices,sequence):
    out_matrix = np.zeros((len(structure_a[1]),len(structure_b[1])))

    for (i, idx_a), (j, idx_b) in itertools.product(
            enumerate(structure_a[1]),enumerate(structure_b[1])):
        helix_a, helix_b = helices[idx_a - 1], helices[idx_b - 1]

        out_matrix[i, j] = levenshtein(
            list(helpers.helix_class_iterator(helix_a,sequence,sort=True)),
            list(helpers.helix_class_iterator(helix_b,sequence,sort=True)))

    return out_matrix

#calculates the absolute difference in energy between helices in each structure
def structure_energy_matrix(structure_a,structure_b,helices,sequence):
    sizes = len(structure_a[1]), len(structure_b[1])
    a_raw = np.array(score_helices.helix_class_energy([helices[idx - 1] for idx in structure_a[1]],sequence))
    b_raw = np.array(score_helices.helix_class_energy([helices[idx - 1] for idx in structure_b[1]],sequence))
    a_energies = np.tile(a_raw,(sizes[1],1))
    b_energies = np.tile(np.expand_dims(b_raw,axis=1),(1,sizes[0]))

    return np.abs(a_energies - b_energies)


#calculates the distance between the centers of the helices in each structure
def helix_distance_matrix(structure_a,structure_b,helices):
    out_matrix = np.zeros((len(structure_a[1]),len(structure_b[1])))

    for (i, idx_a), (j, idx_b) in itertools.product(
            enumerate(structure_a[1]),enumerate(structure_b[1])):
        helix_a, helix_b = helices[idx_a - 1][1], helices[idx_b - 1][1]

        out_matrix[i, j] = abs(helix_a[0] + helix_a[2]/2 - helix_b[0] - helix_b[2]/2) +\
                           abs(helix_a[1] - helix_a[2]/2 - helix_b[1] + helix_b[2]/2)

    return out_matrix
    

#returns the number of elements in the helix class as a default way to handle the cost of
#leaving a helix class unmatched
def default_deletion_cost(helix_class,sequence):
    return helix_class[1][2]

#cost of deleting a helix based on its energy
def energy_deletion_cost(helix_class,sequence):
    return score_helices.helix_class_energy(helix_class,sequence)[0]

#cost of deleting a helix based on the distance to other helices
#at the moment is just the length of the sequence to avoid leaving
#all helices unmatched, but could be made more sophisticated
def distance_deletion_cost(helix_class,sequence):
    return len(sequence)

#calculates the overall unrestricted cost of matching helices between structures including
#the possibilities of deleting/leaving structures unmatched
def matching_cost_matrix(
        structure_a,structure_b,helices,sequence,
        delete_cost_func=None,
        use_edit_distance=False):

    if delete_cost_func is None:
        if use_edit_distance:
            delete_cost_func = default_deletion_cost
        else:
            delete_cost_func =\
                lambda hel, seq : energy_deletion_cost(hel,seq) + distance_deletion_cost(hel,seq)

    #the out matrix consists of the costs of matching helices to other helices,
    #the costs of leaving them unmatched (matching with an index larger than the
    #size of the other structure) and matching non-helices (a cost of 0)
    sizes = len(structure_a[1]), len(structure_b[1])
    matrix_size = sizes[0] + sizes[1]
    out_matrix = np.zeros((matrix_size,matrix_size))

    if use_edit_distance:
        #set the portion of the matrix between helices
        out_matrix[:sizes[0],:sizes[1]] =\
            edit_distance_matrix(structure_a,structure_b,helices,sequence)
    else:
        out_matrix[:sizes[0],:sizes[1]] =\
            structure_energy_matrix(structure_a,structure_b,helices,sequence) +\
            helix_distance_matrix(structure_a,structure_b,helices)
    

    #set the portion of the matrix between helices and blank indices
    for (i,idx) in enumerate(structure_a[1]):
        out_matrix[i,sizes[1]:] =\
            delete_cost_func(helices[idx - 1],sequence)
    for (i,idx) in enumerate(structure_b[1]):
        out_matrix[sizes[0]:,i] =\
            delete_cost_func(helices[idx - 1],sequence)

    return out_matrix

#find a minimal matching between the helices of two structures
#use edit distance specifies whether to use the energy + distance or edit
#distance measures of cost
def match_structure_helices(
        structure_a,structure_b,helices,sequence,use_edit_distance=False):

    cost_matrix = matching_cost_matrix(
        structure_a,structure_b,helices,sequence,
        use_edit_distance=use_edit_distance)

    row_ind,col_ind = linear_sum_assignment(cost_matrix)
    score = cost_matrix[row_ind,col_ind].sum()

    pairs = []
    for x,y in zip(row_ind,col_ind):
        if x < len(structure_a[1]) and y < len(structure_b[1]):
            pairs.append((structure_a[1][x],structure_b[1][y]))

    return pairs, score

if __name__ == "__main__":
    #Code for testing matching functionality
    print("\nFirst 2 structure class descriptions")
    for struct in structures[:2]:
        print("Structure:",struct[1])

    cost_matrix = matching_cost_matrix(
        structures[0],structures[1],helices,sequence)

    print("\nCost matrix:")
    print(cost_matrix)

    pairs, score = match_structure_helices(
        structures[0],structures[1],helices,sequence)

    print("\nTotal matching score: " + str(score))
    print("Matched helices:")
    for pair in pairs:
        print(pair)

