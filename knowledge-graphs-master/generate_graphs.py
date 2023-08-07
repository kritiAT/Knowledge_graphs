"""Script to generate multiple PHKGs"""

# imports

from copy import deepcopy
from email import header
from email.policy import default
import os
from pathlib import Path
from typing import List

import click as click
import networkx as nx
from scipy import rand
from tqdm import tqdm
import matplotlib.pyplot as plt
from ast import literal_eval
import pandas as pd
import random

from src.db_functions import DbHelper
from src.make_graph import PKG, PlotPKG
from src.imp_features import unique_feature_patients, create_akg




@click.group()
def cli():
  pass


@cli.command('pkgs')
@click.argument("num_patients", type=int)
@click.argument("output_folder", type=str)
@click.argument("output_file", type=str)   # store stats table to file
@click.option('-p', '--plots', is_flag=True)
def personalised_kgs(num_patients, output_folder, output_file, plots):
    """
    Main program to generate patients

    Run with python - for instance - generate-patients.py 100 reports.
    """
    Path(output_folder).mkdir(exist_ok=True, parents=True)

    db_helper = DbHelper()

    medi_associations = db_helper.get_medi()
    icd_associations = db_helper.get_icd_associations()
    unique_icds = icd_associations['disease1'].append(icd_associations['disease2']).unique()

    # list of patients
    patient_ids = db_helper.get_data(
        f"SELECT * FROM ml_covid_joined_id LIMIT {num_patients}"
    ).explorys_patient_id.values

    # store stats
    ids: List[int] = [] 
    number_of_nodes: List[int] = []
    number_of_edges: List[int] = []
    number_of_found_edges: List[int] = []

    for i, patient_id in tqdm(
        enumerate(patient_ids),
        desc="Generate graphs",
        total=num_patients,
    ):

        drugs = db_helper.get_drugs(patient_id)
        diags = db_helper.get_diags(patient_id)

        # Commom drugs between patient data and medi data
        lit_drugs = drugs[(drugs['rx_cui'].isin(medi_associations.rxcui.unique()))]
        
        # Commom icd codes between patient data and literature (medi and disgenet)
        lit_diagnosis = diags[(diags['icd_code'].isin(medi_associations.icd_code.unique())) | 
            (diags['icd_code'].isin(unique_icds))]

        # drug_diagnosis realtions between common drugs and diagnosis
        mdas = medi_associations[(medi_associations['rxcui'].isin(lit_drugs.rx_cui.unique())) & 
            (medi_associations['icd_code'].isin(lit_diagnosis.icd_code.unique()))]
        
        # diagnosis_diagnosis realtions between common icds
        ddas = icd_associations[(icd_associations['disease1'].isin(lit_diagnosis.icd_code.unique())) & 
            (icd_associations['disease2'].isin(lit_diagnosis.icd_code.unique()))]

        ids.append(patient_id)    ## store patient ids

        # create KGs
        current_pkg = PlotPKG(patient_drugs=drugs, patient_diagnosis=diags,
            lit_drugs=lit_drugs, mdas=mdas, lit_diagnosis=lit_diagnosis, ddas=ddas, real_associations=True)
        
        number_of_nodes.append(len(current_pkg.graph.nodes))    ## add total nodes
        number_of_edges.append(len(current_pkg.graph.edges))    ## add total edges
        number_of_found_edges.append(len(current_pkg.found_edges))    ## add total found edges

        # store pickle file
        nx.write_gpickle(current_pkg.graph, os.path.join(output_folder, f"{patient_id}_graph.pkl"),)

        # plots
        if plots:

            pkg_complete = current_pkg.complete_PKG()
            pkg_connected = current_pkg.connected_PKG()
            pkg_real_associations = current_pkg.real_associations_PKG()

            # store plots
            pkg_complete.savefig(os.path.join(output_folder, f"{patient_id}_complete.png"))
            pkg_connected.savefig(os.path.join(output_folder, f"{patient_id}_without-isolates.png"))
            pkg_real_associations.savefig(os.path.join(output_folder, f"{patient_id}_only-found.png"))
            
            plt.close("all")
    
    # statistics table
    stats_table = pd.DataFrame(columns=['patient_ids', 'nodes', 'edges', 'found_edges'])

    stats_table['patient_ids'] = ids
    stats_table['nodes'] = number_of_nodes
    stats_table['edges'] = number_of_edges
    stats_table['found_edges'] = number_of_found_edges

    stats_table.to_csv(output_file, index=False, header=True)    ## store table in file
    
    print(f"Average Number of nodes: {sum(number_of_nodes)/num_patients}")
    print(f"Average Number of edges: {sum(number_of_edges)/num_patients}")
    print(f"Average Number of found edges: {sum(number_of_found_edges)/num_patients}")



@cli.command('akgs')
@click.argument("num_patients", type=int)
@click.argument("output_folder", type=str)
@click.argument("phecode", type=str, default='278.11') # morbid obesity
@click.option("-t", "--threshold", type=int, default=50)
def averaged_kgs(num_patients, output_folder, phecode, threshold):
    """
    Command to generate averaged graphs.

    try : python generate_graphs.py akgs  500 ..//..//averaged -t 25
    """

    Path(output_folder).mkdir(exist_ok=True, parents=True)

    db_helper = DbHelper()

    medi_associations = db_helper.get_medi()
    icd_associations = db_helper.get_icd_associations()
    unique_icds = icd_associations['disease1'].append(icd_associations['disease2']).unique()

    # list of patients
    patient_ids = unique_feature_patients(phecode, db_helper=db_helper)
    patient_ids = random.choices(patient_ids, k=num_patients)

    # store common edges and node attributes
    common_edges = {}
    all_nodes = {}

    # store stats
    ids: List[int] = [] 
    number_of_nodes: List[int] = []
    number_of_edges: List[int] = []
    number_of_found_edges: List[int] = []

    for i, patient_id in tqdm(
        enumerate(patient_ids),
        desc="Generate graphs",
        total=num_patients,
    ):

        drugs = db_helper.get_drugs(patient_id)
        diags = db_helper.get_diags(patient_id)

        # Commom drugs between patient data and medi data
        lit_drugs = drugs[(drugs['rx_cui'].isin(medi_associations.rxcui.unique()))]
        
        # Commom icd codes between patient data and literature (medi and disgenet)
        lit_diagnosis = diags[(diags['icd_code'].isin(medi_associations.icd_code.unique())) | 
            (diags['icd_code'].isin(unique_icds))]

        # drug_diagnosis realtions between common drugs and diagnosis
        mdas = medi_associations[(medi_associations['rxcui'].isin(lit_drugs.rx_cui.unique())) & 
            (medi_associations['icd_code'].isin(lit_diagnosis.icd_code.unique()))]
        
        # diagnosis_diagnosis realtions between common icds
        ddas = icd_associations[(icd_associations['disease1'].isin(lit_diagnosis.icd_code.unique())) & 
            (icd_associations['disease2'].isin(lit_diagnosis.icd_code.unique()))]

        ids.append(patient_id)    ## store patient ids

        # create KGs
        current_pkg = PKG(patient_drugs=drugs, patient_diagnosis=diags,
            lit_drugs=lit_drugs, mdas=mdas, lit_diagnosis=lit_diagnosis, ddas=ddas, real_associations=False)
        
        number_of_nodes.append(len(current_pkg.graph.nodes))    ## add total nodes
        number_of_edges.append(len(current_pkg.graph.edges))    ## add total edges

        # store edges and node attributes
        for edge1, edge2, data in current_pkg.graph.edges(data=True):
            edge = '_'.join(sorted((edge1, edge2)))
            if edge in common_edges:
                common_edges[edge]['count'] += 1
            else:
                common_edges[edge] = {'count' : 1, 'color' : data['color']}
            
            if edge1 not in all_nodes:
                all_nodes[edge1] = {'color' : current_pkg.graph.nodes[edge1]['color']}
            if edge2 not in all_nodes:
                all_nodes[edge2] = {'color' : current_pkg.graph.nodes[edge2]['color']}
        
        # # store pickle file
        # nx.write_gpickle(current_pkg.graph, os.path.join(output_folder, f"{patient_id}_graph.pkl"),)

    # averaged edges
    edges_df = pd.DataFrame(common_edges).T.rename_axis('edge').reset_index()
    nodes_df = pd.DataFrame(all_nodes).T.rename_axis('node').reset_index()

    # filter edges based on threshold
    cutoff = num_patients * ( threshold / 100 )
    averaged_edges = edges_df[edges_df['count'] >= cutoff]

    # process nodelist and edgelist
    # seperate edge str and add strength of edges and create nodelist
    averaged_edges[['node1', 'node2']] = averaged_edges['edge'].str.split('_', 1, expand=True)
    averaged_edges.drop(columns='edge', inplace=True)
    averaged_edges = averaged_edges[['node1', 'node2', 'count', 'color']]

    strength = averaged_edges['count'].apply(lambda x : round((x/num_patients) * 100, 2))
    averaged_edges.insert(averaged_edges.shape[1], column='strength', value=strength)

    averaged_nodes = nodes_df[(nodes_df['node'].isin(averaged_edges.node1.unique())) |
        (nodes_df['node'].isin(averaged_edges.node2.unique()))]
    
    # Create averaged plot and store
    akg_plot = create_akg(node_list=averaged_nodes, edge_list=averaged_edges, db_helper=db_helper)

    akg_plot.savefig(os.path.join(output_folder, f"{phecode}_averaged_plot.png"))

    plt.close("all")

    # store averaged edge and node list
    averaged_edges.to_csv(os.path.join(output_folder, f'{phecode}_edgelist.csv'), index=False, header=True)
    averaged_nodes.to_csv(os.path.join(output_folder, f'{phecode}_nodelist.csv'), index=False, header=True)


    print(f'Total edges : {edges_df.shape[0]}')



if __name__ == '__main__':
    cli()