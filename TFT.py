import time

import numpy as np
import pandas as pd
import networkx as nx
import pickle
import sys
import scipy.stats as stats

# sys.path.insert(1, 'D:/Documents/PycharmProjects/PPI/toolbox/')
sys.path.insert(1, '/home/data/dongyadong/Projects/PPI/toolbox/')
import wrappers
import network_utilities


class TFT(object):
    def __init__(self):

        # get the map between protein ID and gene name
        info = pd.read_csv('data/9606.protein.info.v12.0.txt', sep='\t')
        self.p2g_dict = {}
        self.g2p_dict = {}
        for index, row in info.iterrows():
            self.p2g_dict[row['string_protein_id']] = row['preferred_name']
            self.g2p_dict[row['preferred_name']] = row['string_protein_id']

        # get the pathway info
        pathways = open('data/c2.all.v7.5.1.symbols.gmt', 'r')
        self.pathway_dict = {}
        for pathway in pathways:
            pathway = pathway.strip().split('\t')
            self.pathway_dict[pathway[0]] = set(pathway[2:])

        # get the tf target info
        tfts = open('data/c3.tft.v2023.2.Hs.symbols.gmt', 'r')
        self.tft_dict = {}
        for tft in tfts:
            tft = tft.strip().split('\t')
            self.tft_dict[tft[0]] = set(tft[2:])

        self.graph = nx.Graph()
        for key in self.tft_dict:
            self.graph.add_edges_from([(key, target) for target in self.tft_dict[key]])

        # get the microRNA info
        miRNAs = open('data/c3.mir.v2023.2.Hs.symbols.gmt', 'r')
        self.miRNA_dict = {}
        for miRNA in miRNAs:
            miRNA = miRNA.strip().split('\t')
            self.miRNA_dict[miRNA[0]] = set(miRNA[2:])

        self.google_matrix = None
        self.nodes = set(self.graph.nodes())

        self.distance = None
        self.similarity = None

    def init_distance(self, file_name):
        # 'data/tft_shortest_path_length.pickle'
        with open(file_name, 'rb') as f:
            self.distance = pickle.load(f)

    def distance_zscore(self, genes1, genes2, n_random=1000, min_bin_size=100, seed=452456):

        genes1 = list(set(genes1).intersection(self.nodes))
        genes2 = list(set(genes2).intersection(self.nodes))
        dist = self.shortest_distance(genes1, genes2)
        bins = network_utilities.get_degree_binning(self.graph, min_bin_size)
        genes1_random = wrappers.get_random_nodes(genes1, self.graph, bins=bins, n_random=n_random,
                                                  min_bin_size=min_bin_size, seed=seed)
        genes2_random = wrappers.get_random_nodes(genes2, self.graph, bins=bins, n_random=n_random,
                                                  min_bin_size=min_bin_size, seed=seed)

        random_values_list = zip(genes1_random, genes2_random)
        values = []
        for i, values_random in enumerate(random_values_list):
            nodes1, nodes2 = values_random
            values.append(self.shortest_distance(nodes1, nodes2))

        pval = stats.mannwhitneyu(values, dist, alternative='greater').pvalue

        m, s = np.mean(values), np.std(values)
        if s == 0:
            z = 0.0
        else:
            z = (dist - m) / s
        return dist, z, m, s, pval

    def shortest_distance(self, nodes1, nodes2):
        distances = self.distance.loc[nodes1, nodes2].values
        # print('distance.shape:', distances.shape)
        return np.nanmean(distances.min(axis=1))

    def init_similarity(self, file_name):
        with open(file_name, 'rb') as f:
            self.similarity = pickle.load(f)

    def similarity_zscore(self, genes1, genes2, n_random=1000, min_bin_size=100, seed=452456):

        genes1 = list(set(genes1).intersection(self.nodes))
        genes2 = list(set(genes2).intersection(self.nodes))
        dist = self.maximum_similarity(genes1, genes2)
        bins = network_utilities.get_degree_binning(self.graph, min_bin_size)
        genes1_random = wrappers.get_random_nodes(genes1, self.graph, bins=bins, n_random=n_random,
                                                  min_bin_size=min_bin_size, seed=seed)
        genes2_random = wrappers.get_random_nodes(genes2, self.graph, bins=bins, n_random=n_random,
                                                  min_bin_size=min_bin_size, seed=seed)

        random_values_list = zip(genes1_random, genes2_random)
        values = []
        for i, values_random in enumerate(random_values_list):
            nodes1, nodes2 = values_random
            values.append(self.maximum_similarity(nodes1, nodes2))

        pval = stats.mannwhitneyu(values, dist, alternative='less').pvalue

        m, s = np.mean(values), np.std(values)
        if s == 0:
            z = 0.0
        else:
            z = (dist - m) / s
        return dist, z, m, s, pval

    def maximum_similarity(self, nodes1, nodes2):
        similarities = self.similarity.loc[nodes1, nodes2].values
        # if similarities.shape == (0, 0):
        #     return 0.0
        # print('similarities.shape ', similarities.shape)
        return np.nanmean(similarities.max(axis=1))

    def proximity_distance(self, nodes1, nodes2, n_random=1000, min_bin_size=100, seed=452456, lengths=None):
        return wrappers.calculate_proximity(self.graph, nodes1, nodes2, n_random=n_random,
                                            min_bin_size=min_bin_size, seed=seed, lengths=lengths)

    def init_google_matrix(self):
        nodes = list(self.nodes)
        start = time.time()
        self.google_matrix = pd.DataFrame(nx.algorithms.link_analysis.google_matrix(self.graph, nodelist=nodes),
                                          index=nodes, columns=nodes)
        print('It takes %.2f seconds to initialize google matrix' % (time.time() - start))

    def google_matrix_group_sim_zscore(self, genes1, genes2, n_random=1000, min_bin_size=100, seed=452456):

        genes1 = list(set(genes1).intersection(self.nodes))
        genes2 = list(set(genes2).intersection(self.nodes))
        sim = self.google_matrix_group_similarity(genes1, genes2)
        bins = network_utilities.get_degree_binning(self.graph, min_bin_size)
        genes1_random = wrappers.get_random_nodes(genes1, self.graph, bins=bins, n_random=n_random,
                                                  min_bin_size=min_bin_size, seed=seed)
        genes2_random = wrappers.get_random_nodes(genes2, self.graph, bins=bins, n_random=n_random,
                                                  min_bin_size=min_bin_size, seed=seed)

        random_values_list = zip(genes1_random, genes2_random)
        values = np.empty([n_random, 3])
        for i, values_random in enumerate(random_values_list):
            nodes1, nodes2 = values_random
            values[i] = self.google_matrix_group_similarity(nodes1, nodes2)

        pval_inner1 = stats.mannwhitneyu(values[:, 0], sim[0], alternative='less').pvalue
        pval_inner2 = stats.mannwhitneyu(values[:, 1], sim[1], alternative='less').pvalue
        pval_inter = stats.mannwhitneyu(values[:, 2], sim[2], alternative='less').pvalue

        m, s = np.mean(values, axis=0), np.std(values, axis=0)

        z = (np.array(sim) - m) / s
        return sim, z, (m, s), (pval_inner1, pval_inner2, pval_inter)

    def google_matrix_group_similarity(self, genes1, genes2):
        scores_group1 = self.google_matrix_inner_group_similarity(genes1)
        scores_group2 = self.google_matrix_inner_group_similarity(genes2)
        scores_inter_group = self.google_matrix_inter_group_similarity(genes1, genes2)
        return scores_group1, scores_group2, scores_inter_group

    def google_matrix_inner_group_similarity(self, genes):
        similarity = self.google_matrix.loc[genes, genes].values
        np.fill_diagonal(similarity, -np.inf)
        return np.mean(similarity.max(axis=1))

    def google_matrix_inter_group_similarity(self, genes1, genes2):
        similarity = self.google_matrix.loc[genes1, genes2].values
        return np.nanmax(similarity.max(axis=1))
