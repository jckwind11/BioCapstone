# invocation: python3 display_connectomes.py [--all] [68 | 114 | 219 | 448 | 1000]
# 
#
# the default behavior of this code will go through each cvs and display a single subject's connectomes
# from each data set.
#
# when ran with the --all flag, this program will display every subject's connectomes from the given data sets.
#
# you can specifiy which file to analyze, by specifying the number of voxels [68 | 114 | 219 | 448 | 1000]

import sys
import os
import numpy as np
import networkx as nx 
import matplotlib.pyplot as plt

N_SUBJECTS = 70


def readFile(filename, num_of_subjects, specific_voxel_target):
    fileInfo = filename.split('.')[0].split('_')
    type = fileInfo[2]
    voxels = int(fileInfo[3])
    if (specific_voxel_target != -1) and (voxels != specific_voxel_target):
        return
    if type == "SC":
        print("DISPLAYING STRUCTURAL CONNECTIVITY MATRIX WITH " + str(voxels) + " nodes")
    elif type == "FC":
        print("DISPLAYING FUNCTIONAL CONNECTIVITY MATRIX WITH " + str(voxels) + " nodes")
    
    if (num_of_subjects == -1):
        connectome = nx.Graph()
        adj_matrix, edge_list = np.genfromtxt(filename, delimiter=',', max_rows=voxels, skip_header=0), []

        for rn, row in enumerate(adj_matrix):
            for cn, weight in enumerate(row):
                if weight != 0:
                    edge_list.append((rn, cn, weight))
        connectome.add_weighted_edges_from(edge_list)
        if type == "SC":
            plt.title("STRUCTURAL CONNECTIVITY MATRIX WITH " + str(voxels) + " NODES, SUBJECT: 0")
        elif type == "FC":
            plt.title("FUNCTIONAL CONNECTIVITY MATRIX WITH " + str(voxels) + " NODES, SUBJECT: 0")
        nx.draw(connectome)
        plt.show()
    else:
        connectomes = [nx.Graph()] * N_SUBJECTS
        for subject in range(N_SUBJECTS):
            adj_matrix, edge_list = np.genfromtxt(filename, delimiter=',', max_rows=voxels, skip_header=subject*voxels), []

            for rn, row in enumerate(adj_matrix):
                for cn, weight in enumerate(row):
                    if weight != 0:
                        edge_list.append((rn, cn, weight))
            connectomes[subject].add_weighted_edges_from(edge_list)
            if type == "SC":
                plt.title("STRUCTURAL CONNECTIVITY MATRIX WITH " + str(voxels) + " NODES, SUBJECT: " + str(subject))
            elif type == "FC":
                plt.title("FUNCTIONAL CONNECTIVITY MATRIX WITH " + str(voxels) + " NODES, SUBJECT: " + str(subject))
            #nx.draw(connectomes[subject])
            #plt.show()
    return connectomes

def main():
    voxel_target = -1
    num_of_subjects = -1
    if len(sys.argv) == 2:
        if sys.argv[1] == "--all":
            num_of_subjects = 0
        elif sys.argv[1] in {"68", "114", "219", "448"}:
            voxel_target = int(sys.argv[1])

    if len(sys.argv) == 3:
        if sys.argv[1] == "--all":
            num_of_subjects = 0
        if sys.argv[2] in {"68", "114", "219", "448"}:
            voxel_target = int(sys.argv[2])

    for filename in os.listdir(os.getcwd()):
        if filename.endswith(".csv"):
            connectomes = readFile(filename, num_of_subjects, voxel_target)
            
            continue
        else:
            continue

main()