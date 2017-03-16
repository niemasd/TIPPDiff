#! /usr/bin/env python3
'''
Niema Moshiri & Arya Kaul 2017

Differential microbiome profile analysis from TIPP output
The two trees MUST be identical
'''
import argparse,dendropy
from decimal import Decimal
from random import choice
from scipy.stats import expon
from statistics import variance

# parse user arguments
def parseArgs():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f1', '--tipp_placement_file1', required=True,  type=argparse.FileType('r'), default=None, help="TIPP Output Placement File 1")
    parser.add_argument('-f2', '--tipp_placement_file2', required=True,  type=argparse.FileType('r'), default=None, help="TIPP Output Placement File 2")
    parser.add_argument('-b',  '--bootstrap_replicates', required=False, type=int,                    default=100,  help="Number of Bootstrap replicates")
    parser.add_argument('-a',  '--alpha_threshold',      required=False, type=float,                  default=0.05, help="Significance p-value Threshold (alpha)")
    args = parser.parse_args()
    return args

# compute read proportions for all edges
def readProp(reads, pseudo, tree):
    prop = {}
    for edge in tree.postorder_edge_iter():
        prop[edge] = pseudo
    for read in reads:
        for edge,lwr in read:
            prop[edge] += lwr
    for node in tree.postorder_node_iter():
        prop[node.edge] += sum([prop[c.edge] for c in node.child_node_iter()])
    total = prop[tree.seed_node.edge]
    for node in tree.postorder_node_iter():
        prop[node.edge] /= total
    return prop

# TIPPDiff: Perform differential microbiome profile analysis
def TIPPDiff(tree, placements1, placements2, field2num):
    # preprocess placements
    edges = tree.edge_index
    for i in range(len(edges)):
        edges[i].edge_num = i
    tree.seed_node.edge.edge_num = len(edges)
    reads1 = [] # each read is a set of (edge, like_weight_ratio) tuples
    for read in placements1:
        curr = set()
        for p in read['p']:
            edge = edges[p[field2num['edge_num']]]
            curr.add((edge, Decimal(p[field2num['like_weight_ratio']])))
        reads1.append(curr)
    #from random import shuffle; shuffle(edges) # shuffle edges for testing
    reads2 = [] # each read is a set of (edge, like_weight_ratio) tuples
    for read in placements2:
        curr = set()
        for p in read['p']:
            edge = edges[p[field2num['edge_num']]]
            curr.add((edge, Decimal(p[field2num['like_weight_ratio']])))
        reads2.append(curr)

    # compute x, y, and score for all edges
    PSEUDO = Decimal(1.)/Decimal(len(edges))
    x = readProp(reads1, PSEUDO, tree)
    y = readProp(reads2, PSEUDO, tree)
    s = {e:abs(x[e].ln()-y[e].ln()) for e in x} # assuming ln(x/y) ~ Laplace(0,b), so |ln(x/y)| ~ Exponential(1/b)

    # estimate each edge's |X-Y| variance using bootstrapping
    var = {e:[s[e]] for e in s}
    for _ in range(BOOT):
        xP = readProp([choice(reads1) for __ in range(len(reads1))], PSEUDO, tree)
        yP = readProp([choice(reads2) for __ in range(len(reads2))], PSEUDO, tree)
        for e in var:
            var[e].append(abs(xP[e].ln()-yP[e].ln()))
    for e in var:
        var[e] = variance(var[e])
    
    # perform hypothesis test on each edge (with correction)
    p = [] # elements are (edge_num, p, significant) tuples
    significant = False
    for e in var:
        if var[e] == Decimal(0.):
            p.append((e.edge_num, 1, False))
        else:
            currP = expon.sf(float(s[e]), scale=float(var[e].sqrt()))
            sig = False
            if currP <= ALPHA/len(var): # bonferroni correction
                sig = True
                significant = True
            p.append((e.edge_num, currP, sig))
    return significant, p

# run TIPPDiff
if __name__== "__main__":
    args = parseArgs()
    global BOOT
    BOOT = args.bootstrap_replicates
    global ALPHA
    ALPHA = args.alpha_threshold
    tippout1 = eval(args.tipp_placement_file1.read())
    tippout2 = eval(args.tipp_placement_file2.read())
    tree = dendropy.Tree.get(data=tippout1['tree'].replace('[','{').replace(']','}'), schema='newick', is_parse_jplace_tokens=True)
    fields = tippout1['fields']
    field2num = {}
    for i in range(len(fields)):
        field2num[fields[i]] = i
    placements1 = tippout1['placements']
    placements2 = tippout2['placements']
    sig,p = TIPPDiff(tree, placements1, placements2, field2num)
    if sig:
        print("The two samples are significantly different")
    else:
        print("The two samples are not significantly different")
