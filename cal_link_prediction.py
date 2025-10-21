import multiprocessing as mp
from multiprocessing import Manager
import pickle
import pandas as pd
import networkx as nx
import time
import sys

method_list = [
    "resource_allocation_index",
    "jaccard_coefficient",
    "adamic_adar_index",
    "preferential_attachment",
    "cn_soundarajan_hopcroft",
    "ra_index_soundarajan_hopcroft",
    "within_inter_cluster",
    "common_neighbor_centrality",
]


def link_prediction(graph, ebunch, method):
    total = len(ebunch)
    function = 'nx.link_prediction.' + method
    generator = eval(function)(graph, ebunch=ebunch)
    result_dict = {}
    for i, value in enumerate(generator):
        if i % 10000 == 0:
            print(i, '/', total)
        if value[0] not in result_dict:
            result_dict[value[0]] = {value[1]: value[2]}
        else:
            result_dict[value[0]].update({value[1]: value[2]})
    return result_dict


if __name__ == '__main__':
    job_ID = int(sys.argv[1])
    method = method_list[job_ID]
    num_cores = 20
    num_split = 200
    start_t = time.time()
    with open('data/graph_with_combine_score.pickle', 'rb') as f:
        graph = pickle.load(f)

    nodes_list = list(graph.nodes())
    num_nodes = len(nodes_list)
    split = num_nodes / num_split

    pool = mp.Pool(num_cores)
    parameters = []
    for i in range(num_split):
        parameter = []
        for node in nodes_list[int(i * split): int((i + 1) * split)]:
            parameter += [(node, node2) for node2 in nodes_list]
        parameters.append(parameter)

    results = [pool.apply_async(link_prediction, args=(graph, parameters[i], method))
               for i in range(num_split)]

    pool.close()
    pool.join()
    print("Sub-process(es) done.")

    results = [p.get() for p in results]

    elapsed_sec = time.time() - start_t
    print("multiprocessing run: " + "{:.2f}".format(elapsed_sec) + " seconds")

    print('start transfer dict to dataframe...')
    lengths = {}
    for result in results:
        lengths.update(result)
    df = pd.DataFrame(lengths)
    print('writing to file by pickle...')
    with open('data/ppi_' + method + '.pickle', 'wb') as f2:
        pickle.dump(df, f2)

