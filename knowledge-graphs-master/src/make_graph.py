from datetime import timedelta
from typing import Dict
from copy import deepcopy

import matplotlib.pyplot as plt
import networkx as nx

from db_functions import DbHelper
from enum import Enum


class EdgeColor(Enum):
    """Enum class to set edge colors globally"""

    FOUND_DDA = "blue"
    FOUND_MDA = "red"
    DDA = "blue"
    MDA = "red"


class NodeColor(Enum):
    """Enum class to set node color globally"""

    DIAG = "green"
    DRUG = "orange"


class PKG:
    """ Class to make patient PKG """

    def __init__(self, 
        patient_drugs,
        patient_diagnosis,
        lit_drugs,
        mdas,
        lit_diagnosis,
        ddas,
        time_delta = 90,
        real_associations = True,
    ) -> None:
        self.patient_drugs = patient_drugs
        self.patient_diagnosis = patient_diagnosis
        self.lit_drugs = lit_drugs
        self.mdas = mdas
        self.lit_diagnosis = lit_diagnosis
        self.ddas = ddas
        self.time_delta = time_delta

        self.graph = self._create_pkg()
        if real_associations is True:
            self.found_drugs, self.found_diags, self.found_mdas, self.found_ddas = self._add_real_associations()
            self.found_edges = self.found_ddas.union(self.found_mdas)
    
    def _create_pkg(self):
        graph = nx.Graph()

        # Add drugs and diagnosis nodes
        for drug in self.patient_drugs.rx_cui.unique():
            graph.add_node(drug, color=NodeColor.DRUG.value)

        for diag in self.patient_diagnosis.icd_code.unique():
            graph.add_node(diag, color=NodeColor.DIAG.value)
        
        # Add drug-diag associations
        for _, (rxcui, icd) in self.mdas.iterrows():
            if rxcui in graph.nodes and icd in graph.nodes:
                graph.add_edge(rxcui, icd, color=EdgeColor.MDA.value, state="literature")
            
        # Add diag-diag associations
        for _, (disease1, disease2) in self.ddas.iterrows():
            if disease1 != disease2 and disease1 in graph.nodes and disease2 in graph.nodes:
                graph.add_edge(disease1, disease2, color=EdgeColor.DDA.value, state="literature")
        
        return graph
    

    def _add_real_associations(self):
        """ Add real associations based on patient data """

        found_drugs = set()
        found_diags = set()
        found_mdas = set()
        found_ddas = set()

        # drug-diag associations in real
        for rx, icd in self.graph.edges:
            if self.graph[rx][icd]["color"] == EdgeColor.MDA.value:
                diag_date = self.lit_diagnosis.loc[self.lit_diagnosis["icd_code"] == icd, "new_date"]
                drug_date = self.lit_drugs.loc[self.lit_drugs["rx_cui"] == rx, "new_date"]
                for visit in diag_date:
                    # time span in days between diagnosis and prescription date
                    end_date = visit + timedelta(days=self.time_delta)
                    for prescription in drug_date.unique():
                        if visit <= prescription <= end_date:
                            found_drugs.add(rx)
                            found_diags.add(icd)
                            found_mdas.add((rx, icd))
                            self.graph[rx][icd]["color"] = EdgeColor.FOUND_MDA.value
                            self.graph[rx][icd]["state"] = "found"

        # diag-diag associations in real
        for d1, d2 in self.graph.edges:
            if self.graph[d1][d2]["color"] == EdgeColor.DDA.value:
                for icd1, icd2 in [(d1, d2), (d2, d1)]:
                    diag1_date = self.lit_diagnosis.loc[self.lit_diagnosis["icd_code"] == icd1, "new_date"]
                    diag2_date = self.lit_diagnosis.loc[self.lit_diagnosis["icd_code"] == icd2, "new_date"]
                    for diag1_visit in diag1_date:
                        # time span in days between two diagnosis
                        end_date = diag1_visit + timedelta(days=self.time_delta)
                        for diag2_visit in diag2_date:
                            if diag1_visit <= diag2_visit <= end_date:
                                found_diags.update([icd1, icd2])
                                assos = tuple(sorted((icd1, icd2)))
                                found_ddas.add(assos)
                                self.graph[icd1][icd2]["color"] = EdgeColor.FOUND_DDA.value
                                self.graph[icd1][icd2]["state"] = "found"
        
        return found_drugs, found_diags, found_mdas, found_ddas


class PlotPKG(PKG):

    def __init__(self, patient_drugs, patient_diagnosis, lit_drugs, mdas, lit_diagnosis, ddas, time_delta=90, real_associations=True,) -> None:
        super().__init__(patient_drugs, patient_diagnosis, lit_drugs, mdas, lit_diagnosis, ddas, time_delta, real_associations)

        self.found_associations = real_associations
        self.layout = nx.spring_layout(self.graph, k=0.5)
    
    def _make_plot(self, graph, labels = False):

        fig, ax = plt.subplots(figsize=(15, 15))
        layout = self.layout
        attr = list(nx.get_node_attributes(graph, "color").values())  # node_color
        attr2 = list(nx.get_edge_attributes(graph, "color").values())  # edge_color

        if labels is True:
            node_labels = {n: DbHelper().get_labels(node=n) for n in graph.nodes}
            nx.draw(graph, pos=layout, labels=node_labels, with_labels=True, node_color=attr, edge_color=attr2, style='dashed', ax=ax,)
        else:
            nx.draw(graph, pos=layout, with_labels=True, node_color=attr, edge_color=attr2, style="dashed", ax=ax,)
        
        if self.found_associations == True:
            
            nx.draw_networkx_edges(graph, pos=layout, edgelist=list(self.found_mdas), width=4,
                alpha=1, edge_color=EdgeColor.FOUND_MDA.value, label="drug-diag", ax=ax,)
            
            nx.draw_networkx_edges(graph, pos=layout, edgelist=list(self.found_ddas), width=4, alpha=1,
                edge_color=EdgeColor.FOUND_DDA.value, label="diag-diag", ax=ax,)
            
            nx.draw_networkx_nodes(graph, pos=layout, nodelist=list(self.found_drugs), node_size=500,
                node_color=NodeColor.DRUG.value, label="drug", ax=ax, )
            
            nx.draw_networkx_nodes(graph, pos=layout, nodelist=list(self.found_diags), node_size=500,
                node_color=NodeColor.DIAG.value, label="diag", ax=ax, )
            
            ax.legend()

        return fig
    

    def complete_PKG(self, with_labels=False):
        return self._make_plot(graph=self.graph, labels=with_labels)
    

    def connected_PKG(self, with_labels=True):
        sub_graph = deepcopy(self.graph)
        sub_graph.remove_nodes_from(list(nx.isolates(sub_graph)))
        return self._make_plot(graph=sub_graph, labels=with_labels)
    

    def real_associations_PKG(self, with_labels=True):

        if self.found_associations is not True:
            raise AttributeError("'PlotPKG' object has no attribute 'real_associations'. Set parameter 'real_associations=True'.")

        sub_graph = deepcopy(self.graph)

        found_nodes = self.found_drugs.union(self.found_diags)
        remove_nodes = list(set(sub_graph.nodes) - found_nodes)
        remove_edges = set()
        for n1, n2, state in sub_graph.edges(data='state'):
            if state == 'literature':
                remove_edges.add((n1,n2))
        
        sub_graph.remove_nodes_from(remove_nodes)
        sub_graph.remove_edges_from(remove_edges)
        return self._make_plot(graph=sub_graph, labels=with_labels)
    
    # @staticmethod
    # def _get_labels(node):

    #     db_helper = DbHelper()
        
    #     label = db_helper.get_data(f""" SELECT drug_name FROM ka_medi_drugs WHERE rxcui='{node}'""").drug_name.values
    #     if len(label) != 0:
    #         return label[0] 
    #     else:
    #         label = db_helper.get_data(f""" SELECT diagnosis_name FROM ka_medi_diagnosis WHERE icd_code='{node}'""").diagnosis_name.values
    #         if len(label) != 0:
    #             return label[0]
    #         else:
    #             dnet_id = db_helper.get_data(f""" SELECT disease_id FROM ka_disgenet_mappings WHERE icd_code='{node}'""").disease_id.values
    #             if len(dnet_id) != 0:
    #                 label = db_helper.get_data(f""" SELECT disease_name FROM ka_disgenet_labels WHERE disease_id='{dnet_id[0]}'""").disease_name.values
    #                 return label[0] if len(label) != 0 else node

