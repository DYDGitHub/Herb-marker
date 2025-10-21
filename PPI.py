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


class PPI(object):
    def __init__(self, g_name):
        self.g_name = g_name
        if g_name.split('.')[-1] == 'pickle':
            with open(g_name, 'rb') as f:
                self.graph = pickle.load(f)
        else:
            links = pd.read_csv(g_name, sep=' ')
            self.graph = nx.Graph()
            self.graph.add_weighted_edges_from(
                [(row['protein1'], row['protein2'], row['combined_score']) for index, row in links.iterrows()])

            with open('data/graph_with_combine_score.pickle', 'wb') as f:
                pickle.dump(self.graph, f)

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

        # get the microRNA info
        miRNAs = open('data/c3.mir.v2023.2.Hs.symbols.gmt', 'r')
        self.miRNA_dict = {}
        for miRNA in miRNAs:
            miRNA = miRNA.strip().split('\t')
            self.miRNA_dict[miRNA[0]] = set(miRNA[2:])

        self.expression = pd.read_csv('data/TCGA_Firehose_RNASeqV2_LIHC.csv', index_col=0, header=0)
        self.expression = self.expression[self.expression.std(axis=1) != 0]
        self.expression = self.expression[self.expression.columns[self.expression.columns.str[19:22] == '01A']]
        self.expression = np.log1p(self.expression)
        self.gene_set = set(self.expression.index)
        self.nodes = set(self.graph.nodes())
        self.express_corr = None
        self.google_matrix = None
        self.distance = None
        self.similarity = None

    def init_distance(self, file_name):
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
        return np.mean(distances.min(axis=1))

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
        return np.mean(similarities.max(axis=1))

    def init_express_corr(self, method='spearman'):
        self.express_corr = abs(self.expression.T.corr(method=method))

    def proximity_distance(self, nodes1, nodes2, n_random=1000, min_bin_size=100, seed=452456):
        return wrappers.calculate_proximity(self.graph, nodes1, nodes2, n_random=n_random, min_bin_size=min_bin_size,
                                            seed=seed)

    def path_similarity(self, pathway1, pathway2, function='average'):
        pathway_list = self.pathway_dict.values()
        pathway1 = set(pathway1)
        pathway2 = set(pathway2)
        combine = pathway1.union(pathway2)
        sim_list = []
        for pathway in pathway_list:
            if combine.intersection(pathway):
                sim1 = self.dice(pathway1, pathway)
                sim2 = self.dice(pathway2, pathway)
                sim_list.append(sim1 * sim2 / (sim1 + sim2))

        if sim_list:
            if function == 'maximum':
                return max(sim_list)
            elif function == 'average':
                return sum(sim_list) / len(sim_list)
            else:
                raise Exception('"function" should be chosen between "maximum" and "average".')
        else:
            return 0.0

    def tft_similarity(self, tft1, tft2, function='average'):
        tft_list = self.tft_dict.values()
        tft1 = set(tft1)
        tft2 = set(tft2)
        combine = tft1.union(tft2)
        sim_list = []
        for tft in tft_list:
            if combine.intersection(tft):
                sim1 = self.dice(tft1, tft)
                sim2 = self.dice(tft2, tft)
                sim_list.append(sim1 * sim2 / (sim1 + sim2))

        if sim_list:
            if function == 'maximum':
                return max(sim_list)
            elif function == 'average':
                return sum(sim_list) / len(sim_list)
            else:
                raise Exception('"function" should be chosen between "maximum" and "average".')
        else:
            return 0.0

    def miRNA_similarity(self, miRNA1, miRNA2, function='average'):
        miRNA_list = self.miRNA_dict.values()
        miRNA1 = set(miRNA1)
        miRNA2 = set(miRNA2)
        combine = miRNA1.union(miRNA2)
        sim_list = []
        for miRNA in miRNA_list:
            if combine.intersection(miRNA):
                sim1 = self.dice(miRNA1, miRNA)
                sim2 = self.dice(miRNA2, miRNA)
                sim_list.append(sim1 * sim2 / (sim1 + sim2))

        if sim_list:
            if function == 'maximum':
                return max(sim_list)
            elif function == 'average':
                return sum(sim_list) / len(sim_list)
            else:
                raise Exception('"function" should be chosen between "maximum" and "average".')
        else:
            return 0.0

    def path_sim_z_score(self, pathway1, pathway2, n_random=1000, min_bin_size=100, seed=452456):
        disease_targets_symbol = [self.p2g_dict[target] for target in pathway1]
        herb_targets_symbol = [self.p2g_dict[target] for target in pathway2]
        sim = self.path_similarity(disease_targets_symbol, herb_targets_symbol)
        bins = network_utilities.get_degree_binning(self.graph, min_bin_size)
        nodes_from_random = wrappers.get_random_nodes(pathway1, self.graph, bins=bins, n_random=n_random,
                                                      min_bin_size=min_bin_size, seed=seed)
        nodes_to_random = wrappers.get_random_nodes(pathway2, self.graph, bins=bins, n_random=n_random,
                                                    min_bin_size=min_bin_size, seed=seed)
        random_values_list = zip(nodes_from_random, nodes_to_random)
        values = np.empty(len(nodes_from_random))  # n_random
        for i, values_random in enumerate(random_values_list):
            nodes_from, nodes_to = values_random
            nodes_from_symbol = [self.p2g_dict[target] for target in nodes_from]
            nodes_to_symbol = [self.p2g_dict[target] for target in nodes_to]
            values[i] = self.path_similarity(nodes_from_symbol, nodes_to_symbol)
        pval = stats.mannwhitneyu(values, sim, alternative='less').pvalue
        m, s = np.mean(values), np.std(values)
        if s == 0:
            z = 0.0
        else:
            z = (sim - m) / s
        return sim, z, m, s, pval

    def tft_sim_z_score(self, tft1, tft2, n_random=1000, min_bin_size=100, seed=452456):
        disease_targets_symbol = [self.p2g_dict[target] for target in tft1]
        herb_targets_symbol = [self.p2g_dict[target] for target in tft2]
        sim = self.tft_similarity(disease_targets_symbol, herb_targets_symbol)
        bins = network_utilities.get_degree_binning(self.graph, min_bin_size)
        nodes_from_random = wrappers.get_random_nodes(tft1, self.graph, bins=bins, n_random=n_random,
                                                      min_bin_size=min_bin_size, seed=seed)
        nodes_to_random = wrappers.get_random_nodes(tft2, self.graph, bins=bins, n_random=n_random,
                                                    min_bin_size=min_bin_size, seed=seed)
        random_values_list = zip(nodes_from_random, nodes_to_random)
        values = np.empty(len(nodes_from_random))  # n_random
        for i, values_random in enumerate(random_values_list):
            nodes_from, nodes_to = values_random
            nodes_from_symbol = [self.p2g_dict[target] for target in nodes_from]
            nodes_to_symbol = [self.p2g_dict[target] for target in nodes_to]
            values[i] = self.tft_similarity(nodes_from_symbol, nodes_to_symbol)
        pval = stats.mannwhitneyu(values, sim, alternative='less').pvalue
        m, s = np.mean(values), np.std(values)
        if s == 0:
            z = 0.0
        else:
            z = (sim - m) / s
        return sim, z, m, s, pval

    def miRNA_sim_z_score(self, miRNA1, miRNA2, n_random=1000, min_bin_size=100, seed=452456):
        disease_targets_symbol = [self.p2g_dict[target] for target in miRNA1]
        herb_targets_symbol = [self.p2g_dict[target] for target in miRNA2]
        sim = self.miRNA_similarity(disease_targets_symbol, herb_targets_symbol)
        bins = network_utilities.get_degree_binning(self.graph, min_bin_size)
        nodes_from_random = wrappers.get_random_nodes(miRNA1, self.graph, bins=bins, n_random=n_random,
                                                      min_bin_size=min_bin_size, seed=seed)
        nodes_to_random = wrappers.get_random_nodes(miRNA2, self.graph, bins=bins, n_random=n_random,
                                                    min_bin_size=min_bin_size, seed=seed)
        random_values_list = zip(nodes_from_random, nodes_to_random)
        values = np.empty(len(nodes_from_random))  # n_random
        for i, values_random in enumerate(random_values_list):
            nodes_from, nodes_to = values_random
            nodes_from_symbol = [self.p2g_dict[target] for target in nodes_from]
            nodes_to_symbol = [self.p2g_dict[target] for target in nodes_to]
            values[i] = self.miRNA_similarity(nodes_from_symbol, nodes_to_symbol)
        pval = stats.mannwhitneyu(values, sim, alternative='less').pvalue
        m, s = np.mean(values), np.std(values)
        if s == 0:
            z = 0.0
        else:
            z = (d - m) / s
        return d, z, m, s, pval

    def path_tft_miRNA_sim_z_score(self, genes1, genes2, n_random=1000, min_bin_size=100, seed=452456):
        disease_targets_symbol = [self.p2g_dict[target] for target in genes1]
        herb_targets_symbol = [self.p2g_dict[target] for target in genes2]
        sim_path = self.path_similarity(disease_targets_symbol, herb_targets_symbol)
        sim_tft = self.tft_similarity(disease_targets_symbol, herb_targets_symbol)
        sim_miRNA = self.miRNA_similarity(disease_targets_symbol, herb_targets_symbol)

        bins = network_utilities.get_degree_binning(self.graph, min_bin_size)
        nodes_from_random = wrappers.get_random_nodes(genes1, self.graph, bins=bins, n_random=n_random,
                                                      min_bin_size=min_bin_size, seed=seed)
        nodes_to_random = wrappers.get_random_nodes(genes2, self.graph, bins=bins, n_random=n_random,
                                                    min_bin_size=min_bin_size, seed=seed)
        random_values_list = zip(nodes_from_random, nodes_to_random)
        values_path = np.empty(len(nodes_from_random))
        values_tft = np.empty(len(nodes_from_random))
        values_miRNA = np.empty(len(nodes_from_random))
        for i, values_random in enumerate(random_values_list):
            nodes_from, nodes_to = values_random
            nodes_from_symbol = [self.p2g_dict[target] for target in nodes_from]
            nodes_to_symbol = [self.p2g_dict[target] for target in nodes_to]
            values_path[i] = self.path_similarity(nodes_from_symbol, nodes_to_symbol)
            values_tft[i] = self.tft_similarity(nodes_from_symbol, nodes_to_symbol)
            values_miRNA[i] = self.miRNA_similarity(nodes_from_symbol, nodes_to_symbol)

        pval_path = stats.mannwhitneyu(values_path, sim_path, alternative='less').pvalue
        pval_tft = stats.mannwhitneyu(values_tft, sim_tft, alternative='less').pvalue
        pval_miRNA = stats.mannwhitneyu(values_miRNA, sim_miRNA, alternative='less').pvalue
        m = (np.mean(values_path), np.mean(values_tft), np.mean(values_miRNA))
        s = (np.std(values_path), np.std(values_tft), np.std(values_miRNA))

        z_path = (sim_path - m[0]) / s[0] if s[0] else 0.0
        z_tft = (sim_tft - m[1]) / s[1] if s[1] else 0.0
        z_miRNA = (sim_miRNA - m[2]) / s[2] if s[2] else 0.0

        return (sim_path, sim_tft, sim_miRNA), (z_path, z_tft, z_miRNA), (m, s), (pval_path, pval_tft, pval_miRNA)

    @staticmethod
    def dice(set1, set2):
        return len(set1 & set2) * 2.0 / (len(set1) + len(set2))

    def expression_similarity_matrix(self, genes1, genes2, method='pearson', min_bin_size=100, seed=452456):
        genes1_symbol = set([self.p2g_dict[target] for target in genes1]).intersection(self.gene_set)
        genes2_symbol = set([self.p2g_dict[target] for target in genes2]).intersection(self.gene_set)

        bins = network_utilities.get_degree_binning(self.graph, min_bin_size)
        genes1_random = wrappers.get_random_nodes(genes1, self.graph, bins=bins, n_random=1,
                                                  min_bin_size=min_bin_size, seed=seed)
        genes2_random = wrappers.get_random_nodes(genes2, self.graph, bins=bins, n_random=1,
                                                  min_bin_size=min_bin_size, seed=seed)

        genes1_random_symbol = set([self.p2g_dict[target] for target in genes1_random[0]]).intersection(self.gene_set)
        genes2_random_symbol = set([self.p2g_dict[target] for target in genes2_random[0]]).intersection(self.gene_set)

        gene_list = list(genes1_symbol - genes2_symbol) + list(genes1_symbol.intersection(genes2_symbol)) + list(
            genes2_symbol - genes1_symbol) + list(genes1_random_symbol - genes2_random_symbol) + list(
            genes1_random_symbol.intersection(genes2_random_symbol)) + list(
            genes2_random_symbol - genes1_random_symbol)

        len1 = len(genes1_symbol - genes2_symbol)
        len12 = len(genes1_symbol.intersection(genes2_symbol))
        len2 = len(genes2_symbol - genes1_symbol)
        len1_random = len(genes1_random_symbol - genes2_random_symbol)
        len12_random = len(genes1_random_symbol.intersection(genes2_random_symbol))
        len2_random = len(genes2_random_symbol - genes1_random_symbol)

        df = self.expression.T[gene_list]
        correlation = df.corr(method=method)
        return correlation, (len1, len12, len2, len1_random, len12_random, len2_random)

    def expression_group_sim_zscore(self, genes1, genes2, method, n_random=1000, min_bin_size=100, seed=452456):

        sim = self.expression_group_similarity_score(genes1, genes2, method)

        bins = network_utilities.get_degree_binning(self.graph, min_bin_size)
        genes1_random = wrappers.get_random_nodes(genes1, self.graph, bins=bins, n_random=n_random,
                                                  min_bin_size=min_bin_size, seed=seed)
        genes2_random = wrappers.get_random_nodes(genes2, self.graph, bins=bins, n_random=n_random,
                                                  min_bin_size=min_bin_size, seed=seed)

        random_values_list = zip(genes1_random, genes2_random)
        values = np.empty([n_random, 3])
        for i, values_random in enumerate(random_values_list):
            nodes1, nodes2 = values_random
            values[i] = self.expression_group_similarity_score(nodes1, nodes2, method=method)

        pval_inner1 = stats.mannwhitneyu(values[:, 0], sim[0], alternative='less').pvalue
        pval_inner2 = stats.mannwhitneyu(values[:, 1], sim[1], alternative='less').pvalue
        pval_inter = stats.mannwhitneyu(values[:, 2], sim[2], alternative='less').pvalue

        m, s = np.mean(values, axis=0), np.std(values, axis=0)

        z = (np.array(sim) - m) / s
        return sim, z, (m, s), (pval_inner1, pval_inner2, pval_inter)

    def expression_group_similarity_list(self, genes1, genes2, method):
        scores_group1 = self.expression_inner_group_similarity(genes1, method)
        scores_group2 = self.expression_inner_group_similarity(genes2, method)
        scores_inter_group = self.expression_inter_group_similarity(genes1, genes2, method)
        return scores_group1, scores_group2, scores_inter_group

    def expression_group_similarity_score(self, genes1, genes2, method):
        scores_group1 = self.expression_inner_group_similarity(genes1, method)
        scores_group2 = self.expression_inner_group_similarity(genes2, method)
        scores_inter_group = self.expression_inter_group_similarity(genes1, genes2, method)
        return np.nanmean(scores_group1), np.nanmean(scores_group2), np.nanmean(scores_inter_group)

    def expression_inner_group_similarity(self, genes, method):
        gene_symbols = list(set([self.p2g_dict[target] for target in genes]).intersection(self.gene_set))
        correlation = self.express_corr.loc[gene_symbols, gene_symbols].values

        return correlation[np.triu_indices(correlation.shape[0], 1)].flatten()

    def expression_inter_group_similarity(self, genes1, genes2, method):
        gene_symbols1 = list(set([self.p2g_dict[target] for target in genes1]).intersection(self.gene_set))
        gene_symbols2 = list(set([self.p2g_dict[target] for target in genes2]).intersection(self.gene_set))
        correlation = self.express_corr.loc[gene_symbols1, gene_symbols2].values

        return correlation

    def init_google_matrix(self):
        nodes = list(self.graph.nodes())
        gene_symbols = [self.p2g_dict[node] for node in nodes]

        self.google_matrix = pd.DataFrame(nx.algorithms.link_analysis.google_matrix(self.graph, nodelist=nodes),
                                          index=gene_symbols, columns=gene_symbols)

    def google_matrix_group_sim_zscore(self, genes1, genes2, n_random=1000, min_bin_size=100, seed=452456):

        sim = self.google_matrix_group_similarity_score(genes1, genes2)

        bins = network_utilities.get_degree_binning(self.graph, min_bin_size)
        genes1_random = wrappers.get_random_nodes(genes1, self.graph, bins=bins, n_random=n_random,
                                                  min_bin_size=min_bin_size, seed=seed)
        genes2_random = wrappers.get_random_nodes(genes2, self.graph, bins=bins, n_random=n_random,
                                                  min_bin_size=min_bin_size, seed=seed)

        random_values_list = zip(genes1_random, genes2_random)
        values = np.empty([n_random, 3])
        for i, values_random in enumerate(random_values_list):
            nodes1, nodes2 = values_random
            values[i] = self.google_matrix_group_similarity_score(nodes1, nodes2)

        pval_inner1 = stats.mannwhitneyu(values[:, 0], sim[0], alternative='less').pvalue
        pval_inner2 = stats.mannwhitneyu(values[:, 1], sim[1], alternative='less').pvalue
        pval_inter = stats.mannwhitneyu(values[:, 2], sim[2], alternative='less').pvalue

        m, s = np.mean(values, axis=0), np.std(values, axis=0)

        z = (np.array(sim) - m) / s
        return sim, z, (m, s), (pval_inner1, pval_inner2, pval_inter)

    def google_matrix_group_similarity_list(self, genes1, genes2):
        scores_group1 = self.google_matrix_inner_group_similarity(genes1)
        scores_group2 = self.google_matrix_inner_group_similarity(genes2)
        scores_inter_group = self.google_matrix_inter_group_similarity(genes1, genes2)
        return scores_group1, scores_group2, scores_inter_group

    def google_matrix_group_similarity_score(self, genes1, genes2):
        scores_group1 = self.google_matrix_inner_group_similarity(genes1)
        scores_group2 = self.google_matrix_inner_group_similarity(genes2)
        scores_inter_group = self.google_matrix_inter_group_similarity(genes1, genes2)
        return np.nanmean(scores_group1), np.nanmean(scores_group2), np.nanmax(scores_inter_group)

    def google_matrix_inner_group_similarity(self, genes):
        gene_symbols = [self.p2g_dict[target] for target in genes]
        similarity = self.google_matrix.loc[gene_symbols, gene_symbols].values
        return similarity[np.triu_indices(similarity.shape[0], 1)].flatten()

    def google_matrix_inter_group_similarity(self, genes1, genes2):
        gene_symbols1 = [self.p2g_dict[target] for target in genes1]
        gene_symbols2 = [self.p2g_dict[target] for target in genes2]
        similarity = self.google_matrix.loc[gene_symbols1, gene_symbols2].values
        return similarity


