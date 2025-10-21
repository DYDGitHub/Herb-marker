from TFT import TFT
import pickle
import pandas as pd
import networkx as nx
import time
import utils
import sys

i = int(sys.argv[1])

cancer_info = pd.read_excel('data/Cancer_target.xlsx', index_col=0)
cancer_name = cancer_info.index[i]
cancer_targets = set()
for UMLS_ID in cancer_info.loc[cancer_name].UMLS.split():
    cancer_targets.update(utils.get_disease_target_symbols(UMLS_ID))
cancer_targets = list(cancer_targets)

herb_targets = pd.read_csv('data/Data_S1_to_S13/S4-1. HIT_herb_target_data_0412dropna.csv')
herb_target_dict = {}
for index, row in herb_targets[['Herb ID', 'Gene Symbol']].iterrows():
    if row['Herb ID'] not in herb_target_dict:
        herb_target_dict[row['Herb ID']] = [row['Gene Symbol']]
    else:
        herb_target_dict[row['Herb ID']].append(row['Gene Symbol'])

start = time.time()
tft = TFT()
file_name = 'tft_jaccard_coefficient'
# file_name = 'tft_shortest_path_length'
tft.init_similarity('data/%s.pickle' % file_name)
# tft.init_distance('data/%s.pickle' % file_name)

print('It takes %.2f seconds to create the TFT object' % (time.time() - start))

predict_scores = []
for key in herb_target_dict:
    dts = cancer_targets
    hts = herb_target_dict[key]
    # dts = [ppi.g2p_dict[dt] for dt in dts if dt in tft.g2p_dict]
    # hts = [ppi.g2p_dict[ht] for ht in hts if ht in tft.g2p_dict]
    if hts:
        predict_scores.append([key, tft.similarity_zscore(dts, hts, 100)])


herb_ID = []
dist = []
z_score = []
m = []
s = []
pval = []
for score in predict_scores:
    herb_ID.append(score[0])
    dist.append(score[1][0])
    z_score.append(score[1][1])
    m.append(score[1][2])
    s.append(score[1][3])

df = pd.DataFrame(
    {'herb_ID': herb_ID, 'dist': dist, 'z_score': z_score, 'm': m, 's': s})

print('writing to file by pandas...')
df.to_csv('results/cancer_score/TFT_JC/%s_scores.csv' % cancer_name)
