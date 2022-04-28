# invocation: python ./readConnectomes.py [SC | FC] [68 | 114 | 219 | 448 | 1000]
import sys
import pickle
import numpy as np
import networkx as nx 
import matplotlib.pyplot as plt
import math
import random
import os

VOXELS, N_SUBJECTS, ETA,  = 68, 70, .999
# distortion vs random walkers
PLT_Y_LIM, MARKER_SIZE = 200, 2

PICKLE_PATH = os.getcwd() + "/pickles/"
CONNECTOMES_PATH = os.getcwd() + "/connectomes/"

# g1/g2 are lists of adjacency matrices
def draw_figure_3d(g1_mats, type1, g2_mats, type2):

    for n, mats in enumerate([g1_mats, g2_mats]):
        x, y = [], []
        for m in mats:
            for i in m:
                for j in m[i]:
                    x.append(m[i][j][0])
                    y.append(m[i][j][1])

        plt.scatter(x, y, c="blue" if n==0 else "red", s=MARKER_SIZE if n==0 else MARKER_SIZE*2, label=type1 if n==0 else type2, alpha=1 if n==0 else .01)

    if g2_mats == []:
        plt.title("{} {} Connectomes({} Voxels)".format(N_SUBJECTS, type1, VOXELS))

    elif type2 == "random":
        plt.title("{} {} Connectomes({} Voxels) vs. {} Erdos-Renyi Random Graphs({} Nodes)".format(N_SUBJECTS, type1, VOXELS, N_SUBJECTS, VOXELS))

    else:
        plt.title("{} {} Connectomes({} Voxels) vs. {} {} Connectomes({} Voxels)".format(N_SUBJECTS, type1, VOXELS, N_SUBJECTS, type2, VOXELS))


    plt.legend(loc='upper left')
    plt.ylabel("# random walkers")
    plt.xlabel("distortion")
    plt.ylim(0, PLT_Y_LIM)
    plt.show()


def draw_avg_degree_vs_clustering_vs_gre(m1, t1, m2, t2):

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel("Avg. Degree")
    ax.set_ylabel("Avg. Clustering Coefficient")
    ax.set_zlabel("Global Resource Efficiency")
    ax.legend(loc='upper left')

    for (m, t, cmap) in [(m1, t1, "greens"), (m2, t2, "reds")]:
        xs, ys, zs = [], [], []
        for (x, y, z) in m:
            xs.append(x)
            ys.append(y)
            zs.append(z)
            
        ax.scatter(xs, ys, zs, cmap=cmap, label=t)

    if m2 == []:
        ax.set_title("Avg Degree vs. Avg. CC vs. Global Resource Efficiency for 70 Subjects ("+t1+", "+str(VOXELS)+" voxels)")
    else:
        ax.set_title("70 Subjects ({}, {} voxels) vs 70 {}".format(t1, VOXELS, t2))

    plt.show()

# returns highest probability of path from a to b
def getPathProb(g, p, random):
    probability = 1
    for i in range(len(p)-1):
        n, m = p[i], p[i+1]

        if random:
            relative_prob = 1 / g.degree(n)
        else:
            relative_prob = (2**-g[n][m]['weight']) / sum([2**(-g[n][k]['weight']) for k in g[n]])

        probability *= relative_prob

    return probability


def parse_file(type):
    FILE_NAME = CONNECTOMES_PATH+"Individual_Connectomes_"+type+"_"+str(VOXELS)+".csv"
    connectomes = [nx.Graph() for _ in range(N_SUBJECTS)]

    for subject in range(N_SUBJECTS):
        adj_matrix, edge_list = np.genfromtxt(FILE_NAME, delimiter=',', max_rows=VOXELS, skip_header=subject*VOXELS), []
        for rn, row in enumerate(adj_matrix):
            for cn, weight in enumerate(row):
                if weight > 0 and rn != cn:
                    edge_list.append((rn, cn, -math.log(weight)))

        connectomes[subject].add_weighted_edges_from(edge_list)
    return connectomes


# calculate number of random walkers needed to ensure that the shortest path has probability Eta
def get_metrics(connectomes, random):

    figure_3D_matrices, resource_efficiency_plt, transitivity = [], [], []

    for i, c in enumerate(connectomes):
        distortion, agg_resource_efficiency = dict(), 0

        for vox_a in c.nodes():
            distortion[vox_a] = dict()
            for vox_b in c.nodes():

                if vox_a == vox_b: 
                    continue
                
                try:
                    sp = list(nx.shortest_path(c, source=vox_a, target=vox_b, weight=None if random else "weight"))

                except nx.NetworkXNoPath:
                    continue

                prob = getPathProb(c, sp, random)

                if prob == 1:
                    distortion[vox_a][vox_b] = (0, 0)

                else:
                    distortion[vox_a][vox_b] = (1 - prob, math.log(1 - ETA) / math.log(1 - prob))
                    agg_resource_efficiency += ( 1 / distortion[vox_a][vox_b][1])

        # save matricies for later processing
        figure_3D_matrices.append(distortion)

        # format: (avg degreee, global resource efficiency)
        if random:
            resource_efficiency_plt.append((sum([n[1] for n in list(c.degree())]) / VOXELS, nx.average_clustering(connectomes[i]), agg_resource_efficiency / (VOXELS*VOXELS)))
        else:
            resource_efficiency_plt.append((sum([n[1] for n in list(c.degree(weight="weight"))]) / VOXELS, nx.average_clustering(connectomes[i], weight="weight"), agg_resource_efficiency / (VOXELS*VOXELS)))

        print("connectome {} completed.".format(i))

    return (figure_3D_matrices, resource_efficiency_plt)



def main():

    type = None
    # read and validate cmd line args
    if len(sys.argv) == 3:

        if sys.argv[1] in {"SC", "FC", "BOTH"}:
            type = sys.argv[1]

        if sys.argv[2] in {"68", "114", "219", "448", "1000"}:
            VOXELS = int(sys.argv[2])
    else:
        sys.exit()

    # LOAD OR GENERATE CONNECTOMES
    for t in ["SC", "FC"]:
        if not os.path.isfile(PICKLE_PATH+t+"_"+str(VOXELS)+"_graphs.p"):
            graphs = parse_file(type)
            metrics = get_metrics(graphs, False)
            pickle.dump(graphs, open(PICKLE_PATH+t+"_"+str(VOXELS)+"_graphs.p", "wb"))
            pickle.dump(metrics, open(PICKLE_PATH+t+"_"+str(VOXELS)+"_metrics.p", "wb"))

    if not os.path.isfile(PICKLE_PATH+"random_"+str(VOXELS)+"_graphs.p"):
        RG_graphs = [nx.erdos_renyi_graph(VOXELS, random.random()) for _ in range(N_SUBJECTS)]
        RG_metrics = get_metrics(RG_graphs, True)
        pickle.dump(RG_graphs, open(PICKLE_PATH+"random_"+str(VOXELS)+"_graphs.p", "wb"))
        pickle.dump(RG_metrics, open(PICKLE_PATH+"random_"+str(VOXELS)+"_metrics.p", "wb"))


    # PLOT RESULTS
    SC_metrics = pickle.load(open(PICKLE_PATH+"SC_"+str(VOXELS)+"_metrics.p", "rb"))
    FC_metrics = pickle.load(open(PICKLE_PATH+"SC_"+str(VOXELS)+"_metrics.p", "rb"))
    RG_metrics = pickle.load(open(PICKLE_PATH+"random_"+str(VOXELS)+"_metrics.p", "rb"))

    #draw_avg_degree_vs_clustering_vs_gre(SC_metrics[1], "SC", RG_metrics[1], "random")
    draw_avg_degree_vs_clustering_vs_gre(SC_metrics[1], "SC", FC_metrics[1], "FC")
    draw_avg_degree_vs_clustering_vs_gre(SC_metrics[1], "SC", [], "")



if __name__ == "__main__":
    main()