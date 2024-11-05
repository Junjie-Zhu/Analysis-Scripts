import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np
import scipy
import tqdm
from scipy.spatial.distance import squareform


def read_pdb(file_path):
    """
    Read a PDB file and return a `biotite` `AtomArray` object.

    Parameters
    ----------
    file_path : str
        The path to the PDB file.

    Returns
    -------
    biotite.structure.AtomArray
        The `AtomArray` object representing the PDB file.
    """
    return strucio.load_structure(file_path)


class PdbParser():
    def __init__(self,
                 pdb_path,
                 traj_path=None):

        self.structure = read_pdb(pdb_path)

        self.mode = 'single' if self.structure.stack_depth() == 1 else 'multi'

    def split_models(self):
        """
        Split the structure into multiple models.

        Returns
        -------
        list
            A list of `AtomArray` objects, each representing a model.
        """
        assert self.mode == 'multi', 'The input structure is not a multi-model structure.'
        assert isinstance(self.structure, struc.AtomArrayStack), 'The structure is already split into models.'

        self.structure = [self.structure[model] for model in range(self.structure.stack_depth())]

    def select_atoms(self, atom_name):
        """
        Select atoms with the given atom name.

        Parameters
        ----------
        atom_name : str
            The name of the atom to select.

        Returns
        -------
        biotite.structure.AtomArray
            The `AtomArray` object containing the selected atoms.
        """
        if self.mode == 'single':
            return self.structure[self.structure.atom_name == atom_name]
        else:
            return [model[model.atom_name == atom_name] for model in self.structure]


def rmsd_first2all(structure):
    """
    Compute the RMSD between the first model and all other models.

    Returns
    -------
    numpy.ndarray
        The RMSD values.
    """
    assert structure.stack_depth() > 1, 'The input structure is not a multi-model structure.'

    trajectory, _ = struc.superimpose(structure[0], structure)
    rmsds = struc.rmsd(structure[0], trajectory)
    return rmsds


def rmsd_all2all(structure):
    """
    Compute the pairwise RMSD between all models.

    Returns
    -------
    numpy.ndarray
        The pairwise RMSD values.
    """
    assert structure.stack_depth() > 1, 'The input structure is not a multi-model structure.'

    rmsds = np.zeros((len(structure), len(structure)))
    for i in tqdm.tqdm(range(len(structure)), desc='Calculating pairwise RMSD'):
        trajectory, _ = struc.superimpose(structure[i], structure)
        rmsds = struc.rmsd(structure[i], trajectory)
        rmsds[i] = rmsds

    return rmsds


def rmsd_ref2all(structure, reference_structure):
    """
    Compute the RMSD between a reference structure and all other models.

    Parameters
    ----------
    reference_structure : biotite.structure.AtomArray
        The reference structure.

    Returns
    -------
    numpy.ndarray
        The RMSD values.
    """
    assert structure.stack_depth() > 1, 'The input structure is not a multi-model structure.'
    assert reference_structure.stack_depth() == 1, 'The reference structure is not a single-model structure.'

    trajectory, _ = struc.superimpose(reference_structure, structure)
    rmsds = struc.rmsd(reference_structure, trajectory)
    return rmsds


def cluster_by_rmsd(pairwise_rmsds, t=10, criterion='maxclust'):
    """
    Cluster structures based on pairwise RMSD values.
    """

    # Perform hierarchical clustering
    linkage = scipy.cluster.hierarchy.linkage(squareform(pairwise_rmsds), method='average')
    # dendrogram = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True)
    cluster_indices = scipy.cluster.hierarchy.fcluster(linkage, t=t, criterion=criterion)

    # get index of cluster centers
    cluster_center_index = []
    cluster_sizes_ = []

    cluster_indexes, cluster_sizes = np.unique(cluster_indices, return_counts=True)

    for cluster_index, cluster_size in zip(cluster_indexes, cluster_sizes):
        cluster = np.where(cluster_indices == cluster_index)[0]

        center_index = np.argmin(np.mean(pairwise_rmsds[cluster], axis=0))

        cluster_center_index.append(center_index)
        cluster_sizes_.append(cluster_size)

    return linkage, cluster_center_index, cluster_sizes_


def rmsf(structure):
    """
    Compute the root-mean-square fluctuation (RMSF) of each atom.

    Returns
    -------
    numpy.ndarray
        The RMSF values.
    """
    assert structure.stack_depth() > 1, 'The input structure is not a multi-model structure.'

    trajectory, _ = struc.superimpose(structure[0], structure)

    delta_coords = trajectory.coord - structure[0].coord
    rmsf = np.mean(delta_coords ** 2, axis=0) ** 0.5
    return rmsf


def rg(structure):
    """
    Compute the radius of gyration (Rg) of the structure.
    """
    return struc.gyration_radius(structure)


def e2e_distance(structure):
    """
    Compute the end-to-end distance of the structure.
    """

    if structure.stack_depth() > 1:
        e2e_dist = []
        for models in structure:
            e2e_dist.append(np.linalg.norm(models.coord[-1] - models.coord[0]))
        return np.stack(e2e_dist)
    else:
        return np.linalg.norm(structure.coord[-1] - structure.coord[0])


def bond_length(structure, atom_name1, atom_name2):
    """
    Compute the bond length between two atoms.

    Parameters
    ----------
    atom_name1 : str
        The name of the first atom.
    atom_name2 : str
        The name of the second atom.

    Returns
    -------
    numpy.ndarray
        The bond length values.
    """
    if structure.stack_depth() > 1:
        bond_lengths = []
        for model in structure:
            atom1 = model[model.atom_name == atom_name1]
            atom2 = model[model.atom_name == atom_name2]
            bond_lengths.append(np.linalg.norm(atom1.coord - atom2.coord, axis=-1))
        return np.stack(bond_lengths)
    else:
        atom1 = structure[structure.atom_name == atom_name1]
        atom2 = structure[structure.atom_name == atom_name2]
        return np.linalg.norm(atom1.coord - atom2.coord, axis=-1)


def bond_angle(structure, atom_name1, atom_name2, atom_name3):
    """
    Compute the bond angle between three atoms.

    Parameters
    ----------
    atom_name1 : str
        The name of the first atom.
    atom_name2 : str
        The name of the second atom.
    atom_name3 : str
        The name of the third atom.

    Returns
    -------
    numpy.ndarray
        The bond angle values.
    """
    if structure.stack_depth() > 1:
        bond_angles = []
        for model in structure:
            atom1 = model[model.atom_name == atom_name1]
            atom2 = model[model.atom_name == atom_name2]
            atom3 = model[model.atom_name == atom_name3]
            vec1 = atom1.coord - atom2.coord
            vec2 = atom3.coord - atom2.coord
            bond_angles.append(np.arccos(np.sum(vec1 * vec2, axis=-1) / (np.linalg.norm(vec1, axis=-1) * np.linalg.norm(vec2, axis=-1))))
        return np.stack(bond_angles)
    else:
        atom1 = structure[structure.atom_name == atom_name1]
        atom2 = structure[structure.atom_name == atom_name2]
        atom3 = structure[structure.atom_name == atom_name3]
        vec1 = atom1.coord - atom2.coord
        vec2 = atom3.coord - atom2.coord
        return np.arccos(np.sum(vec1 * vec2, axis=-1) / (np.linalg.norm(vec1, axis=-1) * np.linalg.norm(vec2, axis=-1)))


def phi_psi(structure):
    """
    Compute the phi and psi angles of the structure.
    """
    ph, ps, omg = struc.dihedral_backbone(structure)

    if len(ph.shape) > 1:
        ramachandran = {
            'phi': ph[:, 1:-1],
            'psi': ps[:, 1:-1],
        }
    else:
        ramachandran = {
            'phi': ph[1:-1],
            'psi': ps[1:-1],
        }
    return ramachandran

