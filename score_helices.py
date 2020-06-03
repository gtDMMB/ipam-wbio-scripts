
import helpers

#calculates the approximate energy of a base pair
def base_pair_energy(a, b):
    #convert to lowercase
    pair = (a.lower(), b.lower())

    if 'g' in pair and 'c' in pair:
        return 3
    if 'a' in pair and 'u' in pair:
        return 2
    if 'g' in pair and 'u' in pair:
        return 1

    return 0

#calculates the energy estimates for one or more helix classes. 
#Can either be passed a tuple representing a single helix class or
#a list of helix classes
def helix_class_energy(helices, sequence):
    #check whether input is a list
    if not isinstance(helices, list):
        helices = [helices]
    
    results = []
    for single_class in helices:
        #calculate a score for each base pair in the helix class
        scores = [base_pair_energy(*pair) 
                    for pair in helpers.helix_class_iterator(single_class, sequence)]
        results.append(sum(scores))

    return results

#calculates the energy approximation for one or more structures
def structure_energy(structures, helices, sequence):
    #check whether input is a list
    if not isinstance(structures, list):
        structures = [structures]

    results = []
    for single_structure in structures:
        #find the actual helices based on the indices in the structure
        structure_helices = [helices[idx - 1] for idx in single_structure[1]]

        #calculate the energies of the individual helix classes
        energies = helix_class_energy(structure_helices,sequence)
        results.append(sum(energies))

    return results 

if __name__ == "__main__":
    import sys
    import data
    
    if len(sys.argv) != 3:
        print("Usage: python score_helices.py structure.out verbose.txt")
        exit()

    structfile = sys.argv[1]
    verbosefile = sys.argv[2]

    structures = data.Read_Structures(structfile)
    helices = data.Read_HelixClasses(verbosefile)
    sequence = data.Read_Sequence(verbosefile)

    energies = helix_class_energy(helices, sequence)

    print("Helix class energies and descriptions")
    for i, helix_class, energy in zip(range(len(helices)),helices,energies):
        print("Helix class " + str(i+1) + " is:", helix_class[1])
        print("Base pairs:", list(helpers.helix_class_iterator(helix_class,sequence)))
        print("Energy:", energy)

    struct_energies = structure_energy(structures,helices,sequence)

    print("\nFirst 3 structure class energies and descriptions")
    for struct, energy in zip(structures[:3],struct_energies):
        print("Structure:",struct[1])
        print("Energy:",energy)

