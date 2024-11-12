import tree

import numpy as np

import residue_constants as rc

CA_IDX = rc.atom_order['CA']


def get_structure_length(protein_dict):
    """
    Get the length of the protein structure.
    """
    return len(protein_dict["aatype"])


def strip_ends(protein_dict):
    # Strip missing residues on both ends
    modeled_idx = np.where(protein_dict['aatype'] != 20)[0]
    min_idx, max_idx = np.min(modeled_idx), np.max(modeled_idx)
    protein_dict = tree.map_structure(
        lambda x: x[min_idx: (max_idx + 1)], protein_dict)
    return protein_dict


def recenter_and_scale_coords(protein_dict, coordinate_scale=1.0, eps=1e-8):
    # recenter and scale atom positions
    bb_pos = protein_dict['atom_positions'][:, CA_IDX]
    bb_center = np.sum(bb_pos, axis=0) / (np.sum(protein_dict['seq_mask']) + eps)
    centered_pos = protein_dict['atom_positions'] - bb_center[None, None, :]
    scaled_pos = centered_pos * coordinate_scale
    protein_dict['atom_positions'] = scaled_pos * protein_dict['atom_mask'][..., None]
    return protein_dict


