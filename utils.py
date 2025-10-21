import pandas as pd


def get_disease_targets(disease_ID, g2p_dict):
    disease_data = pd.read_excel('data/disease/' + disease_ID + '_disease_gda_summary.xlsx')
    disease_targets = []
    for name in disease_data.Gene:
        if name in g2p_dict:
            disease_targets.append(g2p_dict[name])

    return disease_targets


# get herb targets
def get_herb_targets(herb_ID, g2p_dict):
    herb_data = pd.read_excel('data/herb/' + herb_ID + '-target.xlsx')
    herb_targets = []
    for name in herb_data['Gene Symbol']:
        if name in g2p_dict:
            herb_targets.append(g2p_dict[name])

    return herb_targets


def get_disease_target_symbols(disease_ID):
    disease_data = pd.read_excel('data/disease/' + disease_ID + '_disease_gda_summary.xlsx')
    return disease_data.Gene.tolist()


# get herb target symbols
def get_herb_target_symbols(herb_ID):
    herb_data = pd.read_excel('data/herb/' + herb_ID + '-target.xlsx')
    return herb_data['Gene Symbol'].tolist()
