We intend to gather the scripts used for data analysis and plotting from our previous research in this project.

## Contents
- Processing PDB files
  - Read PDB files
  - Filter structure data
- Analyzing MD trajectories
  - RMSDs (first2all, all2all, ref2all)
  - Rgs, end-to-end distances
  - Bond lengths, angles and dihedrals
  - chemical shifts (N, CA, C, HN, HA, CB)
  - J couplings (HN-HA, HN-CA)  -----  to be updated

## Data Structure of Processed PDBs
- "{PDB_id}.{Model_id}.{Chain_id}.pdb"
  - "atom_positions": np.array(float32)  [N_res, 37, 3]
  - "atom_mask": np.array(float32)  [N_res, 37]
  - "aatype": np.array(int32)  [N_res]
  - "residue_index": np.array(int32)  [N_res]
  - "chain_ids": np.array(int32)  [N_res]

