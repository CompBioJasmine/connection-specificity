#!/usr/bin/python

import sys
sys.path.append("/Users/jasmine/python/python-graph/core")
from pygraph.classes.graph import graph
from filter2HybridData import *
from hypergeometric import *
from connectingNodesMCC import *

#takes a graph and graph/bait name, prints connection specificity and inputs
#to connection specificity to the screen or to file
def connectionSpecificityMain(gr, graphName):
    grA = makeOverlayGraphAnd(graphName)
    grO = makeOverlayGraphOr(graphName)
    ap_list = allAPDataNodes(graphName)
    ap_list_filtered = []
    for node in ap_list:
        if node in gr.nodes():
            ap_list_filtered.append(node)
    ap_list = ap_list_filtered
    #ap_list = grA.nodes()
    #ap_list = grO.nodes()

    print "\nBait\tNumber of Proteins in AP data and PIN (including ones that connect in)\tNumber of Proteins in PIN\tNumber of Proteins Not in AP Data that Interact\tSignificance Value"
    if len(connectingNodes(grA, grO)) != 0:
        significance_value = 0.01/(len(connectingNodes(grA, grO)))
    else:
        significance_value = "N/A"
    print str(graphName) + "\t" + str(len(ap_list)) + "\t" + str(len(gr.nodes())) + "\t" + str(len(connectingNodes(grA, grO))) + "\t" + str(significance_value)
    print "\nProtein\tNeighbors\tNeighbors in AP\tConnectionSpecificity"
    node_list = []
    for node in gr.nodes():
        if node not in ap_list:
            node_list.append(node)
    for node in node_list:
        intersecting_nodes = []
        for neighbor in gr.neighbors(node):
            if neighbor in ap_list:
                intersecting_nodes.append(node)
        
        connection_specificity = hypergeometric(len(intersecting_nodes), len(gr.neighbors(node)), len(ap_list), len(gr.nodes()))
        if connection_specificity < 1.0:
            print str(node) + "\t" + str(len(gr.neighbors(node))) + "\t" + str(intersecting_nodes) + "\t" + str(connection_specificity)

    return

gr = makeUniverseGraph()
for name in ["MED12C", "MED12L", "MED13", "MED13L", "CDK8", "CDK19", "BRD4", "CyclinC", "MED26C"]:
    connectionSpecificityMain(gr, name)

exit()

