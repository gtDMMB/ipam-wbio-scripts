import os
import helpers

#populate struct_list with tuples (index, list of helix classes) from data stored in inFile. 
def Read_Structures(inFile,profiles=False): 
    struct_list = []
    prof_list = []
    with open(inFile,'r') as f:
        for line in f:
            tokens = line.split()
            if tokens[0]=="Structure":
                index = int(tokens[1].rstrip(':'))
                struct_list.append((index,list(set((int(s) for s in tokens[2:] if s.isdigit())))))
            elif tokens[0]=="->":
                prof_list.append((index,[int(s) for s in tokens[1:] if s.isdigit()]))

    if profiles:
        return struct_list, prof_list

    return struct_list

#Read helix classes from verbose output of RNAprofile; stores ID in position 0, (i,j,k) in position 1 and freq in position 2
def Read_HelixClasses(inFile):
    hc_list = dict()
    with open(inFile, 'r') as f:
        for line in f:
            if line.startswith("Helix"):
                temp_list = [int(s) for s in line.split() if s.isdigit()]
                hc_list[temp_list[0]] = {"id":temp_list[0], "ijk":(temp_list[1], temp_list[2], temp_list[3]), "freq":temp_list[4]}
    return hc_list

#Read helix sequence from verbose output of RNAprofile; creates a list of the characters in the sequence
def Read_Sequence(inFile):
    with open(inFile, 'r') as f:
        first_line = f.readline()
        sequence = first_line.split()[4]
        return list(sequence)

def Calculate_Sequence_Data(sequenceFile,RNAStructureDirectory,seed):
    if RNAStructureDirectory == "":
        print("No RNAStructure directory provided and new seeds tested. Now exiting")
        exit()
    os.system('bash compute_rna_data.sh {} {} {}'.format(RNAStructureDirectory,sequenceFile,seed))    

def Get_File_Prefix(sequenceFile):
    sequence_name = os.path.splitext(os.path.basename(sequenceFile))[0]
    if sequence_name.endswith("_seq"):
        sequence_name = sequence_name[:-4]
    return sequence_name

def Load_All_Data(sequenceFile,RNAStructureDirectory="",
	recalculate_sample=False,seed=1234):
    sequence_name = Get_File_Prefix(sequenceFile)
    data_dir = "./test_files/{}/{}".format(sequence_name,seed)

    if not os.path.isdir(data_dir) or recalculate_sample:
        Calculate_Sequence_Data(sequenceFile,RNAStructureDirectory,seed)

    try:
        sequence = Read_Sequence("{}/output.txt".format(data_dir))
        helices = Read_HelixClasses("{}/output.txt".format(data_dir))
        structures, profiles = Read_Structures("{}/structure.out".format(data_dir),True)
        structures = helpers.Deduplicate_structures(structures)

        return sequence, helices, structures, profiles
    except Exception as e:
        print("Error loading data: {}".format(e))
        print("Recalculating sequence data")

        Calculate_Sequence_Data(sequenceFile,RNAStructureDirectory,seed)
    
    sequence = Read_Sequence("{}/output.txt".format(data_dir))
    helices = Read_HelixClasses("{}/output.txt".format(data_dir))
    structures, profiles = Read_Structures("{}/structure.out".format(data_dir),True)
    structures = helpers.Deduplicate_structures(structures)

    return sequence, helices, structures, profiles
    
def Save_Partition(outFile,partition,data,structure_partition=True,sequence=[],ch_index=None):
    total_counts = [sum(data[helix_id][0] for helix_id in community) for community in partition]
    partition = [x for _,x in reversed(sorted(zip(total_counts,partition)))]

    with open (outFile, 'w') as f:
        f.write("Sequence: {}\n".format("".join(sequence)))
        if ch_index is not None:
            f.write("CH-Index: {}\n".format(ch_index))
        f.write("\n")
        if structure_partition:
            f.write("partition id, multiplicity, structure helices\n")
        else:
            f.write("partition id, helix id, helix ijk\n")
        for i, community in enumerate(partition):
            total_count = sum(data[helix_id][0] for helix_id in community)
            all_helices = sorted(set(sum((list(data[helix_id][1]) for helix_id in community),[])))

            f.write("----------\n")
            f.write("Summary: {}, {}, {}\n".format(
                i, total_count,
                " ".join(str(x) for x in all_helices)))
            for helix_id in community:
                if structure_partition:
                    f.write("{}, {}, {}\n".format(
                        i, data[helix_id][0],
                        " ".join(str(x) for x in data[helix_id][1])))
                else:
                    f.write("{}, {}, {}\n".format(
                        i, data[helix_id]["id"], data[helix_id]["ijk"]))
        

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("Currently should be passed a file to read the sequece from and the RNAStructure directory")
        print("For example: python data.py [RNAStructure/directory/] [sequence_file.txt]")
        exit()

    sequence, helices, structures, profiles =\
        Load_All_Data(sys.argv[2], sys.argv[1])
    print(sequence[:10])
    print(helices[:10])
    print(structures[:10])

