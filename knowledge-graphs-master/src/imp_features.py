# script for imporatant features
# for averaged KGs

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import pandas as pd

from db_functions import DbHelper


# class ImpFeatures:
#     def __init__(self, phewas_code:str) -> None:
#         self.db_helper = DbHelper()
#         self.phe_code = phewas_code
    
#     def get_unique_patients(self) -> set:
#         icd_list = self.db_helper.get_features(self.phe_code)
#         unique_patients = set()
#         for icd in icd_list:
#             patients = self.db_helper.get_patient_list(icd)
#             unique_patients.update(patients)
#         return list(unique_patients)

# # DbHelper().get_labels(node)

#db_helper = DbHelper()

def unique_feature_patients(phe_code: str, db_helper) -> set:
    """ Gest list of unique patients with specific feature. """
    icd_list = db_helper.get_features(phe_code)
    unique_patients = set()
    for icd in icd_list:
        patients = db_helper.get_patient_list(icd)
        unique_patients.update(patients)
    return list(unique_patients)

def create_akg(node_list: pd.DataFrame, edge_list: pd.DataFrame, db_helper, labels: bool = True):
    """ Creates an averaged knowledge graph. """
    graph = nx.Graph()

    for _, node in node_list.iterrows():
        graph.add_node(node.node, color=node.color)
    
    for _, edge in edge_list.iterrows():
        graph.add_edge(edge.node1, edge.node2, color=edge.color, strength=edge.strength/10)
    
    # graph.remove_nodes_from(list(nx.isolates(graph)))

    nx.set_node_attributes(graph, dict(graph.degree), 'strength')

    edge_color = list(nx.get_edge_attributes(graph,'color').values())
    node_color = list(nx.get_node_attributes(graph,'color').values())

    edge_weight = list(nx.get_edge_attributes(graph,'strength').values())
    node_weight = list(nx.get_node_attributes(graph,'strength').values())
    node_weight = [node*200 for node in node_weight]

    fig, ax = plt.subplots(figsize=(15, 15))
    layout = nx.spring_layout(graph, k=0.5)
    
    if labels is True:
        node_labels = {n: db_helper.get_labels(node=n) for n in graph.nodes}
        nx.draw(graph, pos=layout, labels=node_labels, with_labels=True, node_size=node_weight, node_color=node_color, edge_color=edge_color, width=edge_weight)
    else:
        nx.draw(graph, pos=layout, with_labels=True, node_size=node_weight, node_color=node_color, edge_color=edge_color, width=edge_weight)

    mdas = mpatches.Patch(color='red', label='drug-diag')
    ddas = mpatches.Patch(color='blue', label='diag-diag')
    diag = mpatches.Patch(color='green', label='diag')
    drug = mpatches.Patch(color='orange', label='drug')

    ax.legend(handles=[drug, diag, mdas, ddas])

    return fig