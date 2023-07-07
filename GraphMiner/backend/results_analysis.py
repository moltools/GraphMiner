#!/usr/bin/env python3
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdFMCS, Draw
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=False
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage
# from sklearn.decomposition import PCA


def mol_to_fingerprint(mol: Chem.Mol, num_bits: int = 2048, radius: int = 3):
    """
    Creates Morgan fingerprint from RDKit Mol

    input:
    mol - molecule in Chem.Mol format
    num_bits - number of bits for making the fingerprint (int)
    radius - radius for making the fingerprint(int)

    output:
    bit_fingerprint - Fingerprint in a Numpy Array
    """
    # print(Chem.MolToSmiles(mol))
    bit_fingerprint = np.zeros((0,), dtype=int)
    morgan_bit_vector = AllChem.GetMorganFingerprintAsBitVect(mol, radius, num_bits)
    DataStructs.ConvertToNumpyArray(morgan_bit_vector, bit_fingerprint)
    return bit_fingerprint

def tanimoto_coefficient(first_fingerprint: np.array, second_fingerprint: np.array) -> float:
    """Calculates Tanimoto coefficient between two fingerprints

    input:
    first_fingerprint - Fingerprint in a Numpy Array
    second_fingerprint - Fingerprint in a Numpy Array

    output:
    tanimoto coefficient (float)
    """
    return (
        np.logical_and(first_fingerprint, second_fingerprint).sum() / (
            float(np.logical_or(first_fingerprint, second_fingerprint).sum())
        )
    )

def create_groups_dendrogram(dn):
    """
    Create the different groups which are shown in the dendrogram

    input:
    dn - information resulting from the dendrogram function

    output:
    vallist - dictionary with as key the groupnumber (str) and as value the
    substructures (list of str)
    """
    vallist = {}
    groupnum = 0
    vallist[str(groupnum)] = [dn['ivl'][0]]
    for index in range(1, len(dn['leaves_color_list'])):
        if dn['leaves_color_list'][index] != dn['leaves_color_list'][
            index - 1]:
            groupnum += 1
            vallist[str(groupnum)] = [dn['ivl'][index]]
        elif dn['leaves_color_list'][index] == dn['leaves_color_list'][
            index - 1]:
            vallist[str(groupnum)].append(dn['ivl'][index])
    return vallist

def draw_mol_fig(vallist, filepath):
    '''
    Draw the molecule as a structure

    input:
    vallist - dictionary with as key the groupnumber (str) and as value the
    substructures (list of str)
    filepath - pathway as to where to store the file
    '''
    for group in vallist:
        plt.title("Group " + group)
        mols = [Chem.MolFromSmiles(mol) for mol in vallist[group]]
        res = rdFMCS.FindMCS(mols)
        pattern = Chem.MolFromSmarts(res.smartsString)
        totalpath = filepath + '/mcs_group' + str(group) + '.png'
        Draw.MolToFile(pattern, totalpath)
        plt.title("Group " + group)
        lengths = [len(mol) for mol in vallist[group]]
        biggest = vallist[group][lengths.index(max(lengths))]
        bigmol = Chem.MolFromSmiles(biggest)
        totalpath = filepath + '/biggest_group' + str(group) + '.png'
        Draw.MolToFile(bigmol, totalpath)
    return


def plot_dendrogram(dist_matrix, substrsmiles, filename, args):
    '''
    Create and plot the dendrogram

    input:
    dist_matrix - matrix containing all pairwise distances
    substrsmiles - a list of all the substructures which are over/under enriched
    filename - pathway to where to store the output file
    args - arguments from command line

    output:
    dn - all information regarding the dendrogram
    '''
    X = squareform(dist_matrix)
    Z = linkage(X, "ward")
    fig = plt.figure(figsize=(25, 10))
    dn = dendrogram(Z, orientation="right", labels = substrsmiles, color_threshold=args.CutOffDendrogram)
    plt.savefig(filename)
    return dn


