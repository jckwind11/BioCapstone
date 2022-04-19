# invocation: python ./readConnectomes.py [SC | FC] [68 | 114 | 219 | 448 | 1000]
import sys
import numpy as np
import networkx as nx 
import matplotlib.pyplot as plt
import math
import time

TYPE, VOXELS, N_SUBJECTS, ETA = "SC", 68, 70, .999

def draw_figure_3d(m, subj):
    x, y = [], []
    for i in m:
        for j in m[i]:
            x.append(m[i][j][0])
            y.append(m[i][j][1])
    plt.scatter(x, y)
    plt.title("Subject "+str(subj)+" "+TYPE+", "+str(VOXELS)+" voxels")
    plt.ylabel("# random walkers")
    plt.xlabel("distortion")
    plt.show()

def draw_avg_degree_vs_gre(m):
    x, y = [], []
    for val in m:
        x.append(val[0])
        y.append(val[1])

    plt.scatter(x, y)
    plt.title("Avg Degree vs. Global Resource Efficiency for 70 Subjects, "+TYPE+" - "+str(VOXELS)+" voxels")
    plt.ylabel("Global Resource Efficiency")
    plt.xlabel("Avg. Degree")
    plt.show()

# returns highest probability of path from a to b
def getPathProb(g, p):
    probability = 1
    for i in range(len(p)-1):
        n, m = p[i], p[i+1]
        probability *= 2**(-g[n][m]["weight"])
    return probability

# read and validate cmd line args
if len(sys.argv) == 3:
    if sys.argv[1] in {"SC", "FC"}:
        TYPE = sys.argv[1]
    if sys.argv[2] in {"68", "114", "219", "448", "1000"}:
        VOXELS = int(sys.argv[2])

FILE_NAME = "Individual_Connectomes_"+TYPE+"_"+str(VOXELS)+".csv"
connectomes = [nx.Graph() for _ in range(N_SUBJECTS)]

for subject in range(N_SUBJECTS):
    adj_matrix, edge_list = np.genfromtxt(FILE_NAME, delimiter=',', max_rows=VOXELS, skip_header=subject*VOXELS), []
    for rn, row in enumerate(adj_matrix):
        for cn, weight in enumerate(row):
            if weight > 0 and rn != cn:
                edge_list.append((rn, cn, -math.log(weight)))

    connectomes[subject].add_weighted_edges_from(edge_list)
    # nx.draw(connectomes[subject])
    # plt.show()

# list of adj matrices
figure_3D_matrices, resource_efficiency_plt = [], []
start_time = time.time()

# calculate number of random walkers needed to ensure that the shortest path has probability Eta
for i, c in enumerate(connectomes):
    distortion, agg_resource_efficiency = dict(), 0

    for vox_a in c.nodes():
        distortion[vox_a] = dict()
        for vox_b in c.nodes():
            sp = list(nx.shortest_path(c, source=vox_a, target=vox_b, weight="weight"))
            prob = getPathProb(c, sp)
            if prob == 1:
                distortion[vox_a][vox_b] = (0, 0)
            else:
                distortion[vox_a][vox_b] = (1 - prob, math.log(1 - ETA) / math.log(1 - prob))
                agg_resource_efficiency += ( 1 / distortion[vox_a][vox_b][1])

    # save matricies for later processing
    figure_3D_matrices.append(distortion)
    # x - avg degreee, y - global resource efficiency
    resource_efficiency_plt.append((sum([n[1] for n in list(c.degree(weight="weight"))]) / VOXELS, agg_resource_efficiency / (VOXELS*VOXELS))) 

    # output 
    # draw_figure_3d(distortion, i)

draw_avg_degree_vs_gre(resource_efficiency_plt)
print("run time: ", (time.time() - start_time)/60, "mins")

# check if rij is symmetric 