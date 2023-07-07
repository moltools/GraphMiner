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


def mol_to_fingerprint(mol: Chem.Mol, num_bits: int = 2048, radius: int = 3) -> np.array:
    """Creates Morgan fingerprint from RDKit Mol."""
    # print(Chem.MolToSmiles(mol))
    bit_fingerprint = np.zeros((0,), dtype=int)
    morgan_bit_vector = AllChem.GetMorganFingerprintAsBitVect(mol, radius, num_bits)
    DataStructs.ConvertToNumpyArray(morgan_bit_vector, bit_fingerprint)
    return bit_fingerprint

def tanimoto_coefficient(first_fingerprint: np.array, second_fingerprint: np.array) -> float:
    """Calculates Tanimoto coefficient between two fingerprints."""
    return (
        np.logical_and(first_fingerprint, second_fingerprint).sum() / (
            float(np.logical_or(first_fingerprint, second_fingerprint).sum())
        )
    )

def create_groups_dendrogram(dn):
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

# def PCA():
#     pca = PCA(n_components=2)
#     pcs = pca.fit_transform(fps)
#     ev = pca.explained_variance_ratio_
#
#     x, y = pcs[:, 0], pcs[:, 1]
#     ev_x, ev_y = round(ev[0] * 100, 2), round(ev[1] * 100, 2)
#
#     # plt.scatter(x, y, s=1)
#     # plt.xlabel(f"PC1 ({ev_x}%)")
#     # plt.ylabel(f"PC2 ({ev_y}%)")
#     # plt.show()
#     return

def draw_mol_fig(vallist, filepath):
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
    X = squareform(dist_matrix)
    Z = linkage(X, "ward")
    fig = plt.figure(figsize=(25, 10))
    dn = dendrogram(Z, orientation="right", labels = substrsmiles, color_threshold=args.CutOffDendrogram)
    plt.savefig(filename)
    return dn


