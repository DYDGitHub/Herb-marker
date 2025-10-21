from PPI import PPI
import multiprocessing as mp
import pickle
import pandas as pd
import networkx as nx
import time


def shortest_distance_zscore(df, ppi, herb_target_dict, symptom_target_dict):
    predict_scores = []
    for index, row in df.iterrows():
        if index % 100 == 0:
            print(index)

        if row['SID'] in symptom_target_dict:
            dts = symptom_target_dict[row['SID']]
            hts = herb_target_dict[row['herb_UID']]
            dts = [ppi.g2p_dict[dt] for dt in dts if dt in ppi.g2p_dict]
            hts = [ppi.g2p_dict[ht] for ht in hts if ht in ppi.g2p_dict]
            predict_scores.append(
                [row['herb_UID'], row['SID'], ppi.distance_zscore(dts, hts, 10), row['clinical occurrence'],
                 row['effective'], row['PSMeffective']])

    return predict_scores


if __name__ == '__main__':

    herb_target = pd.read_csv('data/Data_S1_to_S13/S4-1. HIT_herb_target_data_0412dropna.csv')
    herb_target_dict = {}
    for index, row in herb_target[['tcm_id', 'Gene Symbol']].iterrows():
        if row['tcm_id'] not in herb_target_dict:
            herb_target_dict[row['tcm_id']] = [row['Gene Symbol']]
        else:
            herb_target_dict[row['tcm_id']].append(row['Gene Symbol'])

    symptom_target = pd.read_excel('data/Data_S1_to_S13/S6. Hospital_symptom_gene associations.xlsx')
    symptom_target_dict = {}
    for index, row in symptom_target.iterrows():
        symptom_target_dict[row['SID']] = row['gene_list'].strip().split('|')

    start = time.time()
    ppi = PPI('data/graph_with_combine_score.pickle')
    # file_name = 'ppi_jaccard_coefficient'
    file_name = 'ppi_shortest_path_length'
    # ppi.init_similarity('data/%s.pickle' % file_name)
    ppi.init_distance('data/%s.pickle' % file_name)
    print('It takes %.2f seconds to create the PPI object' % (time.time() - start))

    df = pd.read_csv('data/Data_S1_to_S13/S8. herb_SID_proximity_nofilterchemical.csv')
    df1 = df[['herb_UID', 'SID', 'clinical occurrence', 'effective', 'PSMeffective']]

    num_pairs = df1.shape[0]

    num_cores = 20
    split = num_pairs / 20
    pool = mp.Pool(20)
    results = [pool.apply_async(shortest_distance_zscore, args=(df1[int(i * split): int((i + 1) * split)], ppi,
                                                                herb_target_dict, symptom_target_dict)) for i in
               range(num_cores)]
    results = [p.get() for p in results]

    pool.close()
    pool.join()

    elapsed_sec = time.time() - start
    print("multiprocessing run: " + "{:.2f}".format(elapsed_sec) + " seconds")

    scores_list = []
    for result in results:
        scores_list.extend(result)


    herb_UID = []
    SID = []
    dist = []
    z_score = []
    clinical_occurrence = []
    effective = []
    PSMeffective = []
    for score in scores_list:
        herb_UID.append(score[0])
        SID.append(score[1])
        dist.append(score[2][0])
        z_score.append(score[2][1])
        clinical_occurrence.append(score[3])
        effective.append(score[4])
        PSMeffective.append(score[5])

    df = pd.DataFrame(
        {'herb_UID': herb_UID, 'SID': SID, 'dist': dist, 'z_score': z_score, 'clinical_occurrence': clinical_occurrence,
         'effective': effective, 'PSMeffective': PSMeffective})

    print('writing to file by pandas...')
    df.to_csv('results/predict_scores_%s_SID.csv' % file_name)
