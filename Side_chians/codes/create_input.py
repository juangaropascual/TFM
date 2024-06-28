import pandas as pd
import numpy as np
import os
from itertools import combinations
import mdtraj as md
from argparse import ArgumentParser


def args_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--SMILES",
        type=str,
        help="Canonical SMILES of the ligand.",
        required=True,
    )
    parser.add_argument(
        "--run_name",
        type=str,
        help="Name of the job.",
        required=True,
    )
    args = parser.parse_args()
    return args


def new_pdb(output_pdb, prot_a, prot_b, separate_chains=[10, 0, 0]):
    """
    This function concatenates multiple mdtraj objects into a single PDB file spacing them by a factor.

    Parameters
    ----------
    output_pdb : str. The path to the output PDB file.
    chains : list. A list of mdtraj objects.
    separate_chains : (list, optional) The distance (in Angstroms) to separate each chain, by default [10, 0, 0]

    Returns
    -------
    None. Generates a pdb file with the merged proteins properly spaced.

    """

    prot_a = md.load(prot_a)
    prot_b = md.load(prot_b)

    top_a = prot_a.top  # extract protein topology
    top_b = prot_b.top

    xyz_a = prot_a.xyz  # extract the protein coordinates
    xyz_b = prot_b.xyz

    xyz_b += separate_chains  # apply separation desired for every chain

    merged_top = top_a.join(top_b)  # merge into a single pdb
    merged_xyz = np.concatenate((xyz_a, xyz_b), axis=1)

    merged_pdb = md.Trajectory(xyz=merged_xyz, topology=merged_top)

    merged_pdb.save(f"{output_pdb}")
    print(f"File {output_pdb} generated successfully")


def generate_input(structure_dir, smiles, run_name):

    # Get all the structures
    folders = ["AF", "DF", "OF"]
    for f in folders:
        print(f)
        try:
            files = os.listdir(f"{structure_dir}/{f}")
            files.sort()
        except:
            continue

        try:
            os.mkdir(f"protein_structures/{run_name}")
        except:
            pass

        n_structures = len(files)

        # # Generate all possible pairs
        ind = np.arange(n_structures)
        # pair_ind = list(combinations(ind, 2))

        # print(
        #     f"--> Using {n_structures} structures to generate {len(pair_ind)} different possible pairs."
        # )
        # Generate the input file
        input = pd.DataFrame(
            columns=[
                "complex_name",
                "ligand_description",
                "protein_path",
                "protein_sequence",
            ]
        )
        ids = []
        s = []
        p = []
        kd_list = []

        # Open the folding file and search if this file contains any affinities.
        fold = pd.read_csv(f"inputs/folding_{run_name}_input.csv")
        
        no_aff = False
        try:
            aff = fold['Affinity']
        except:
            no_aff = True

        for name in ind:
            prot = files[name]
            # prot_b = files[pair[1]]
            id = f"{prot.split('.')[0]}"

            if not no_aff:
                # aff_a = fold[fold['id']==prot_a.split('.')[0]]['Affinity'].values[0]
                # aff_b = fold[fold['id']==prot_b.split('.')[0]]['Affinity'].values[0]
                aff_prot = fold[fold['id']==prot.split('.')[0]]['Affinity'].values[0]
                kd_list.append(aff_prot)
            if not os.path.exists(
                f"protein_structures/{run_name}//{f}"
            ):
                os.makedirs(f"protein_structures/{run_name}/{f}")

            ids.append(id)
            s.append(smiles)
            p.append(
                os.path.abspath(
                    f"protein_structures/{run_name}/{f}/{id}.pdb"
                )
            )

        input["complex_name"] = ids
        input["ligand_description"] = s
        input["protein_path"] = p

        input.to_csv(f"inputs/docking_{run_name}_{f}_input.csv", index=False)


    most_affine = pd.DataFrame(columns=["ID", "Kd"])
    most_affine["ID"] = ids

    if no_aff:    
        most_affine["Kd"] = [1 for i in range(len(ids))]
        most_affine.to_csv(f"results/most_affine_{run_name}_dummy.csv", index=False)
    else:
        most_affine["Kd"] = kd_list
        most_affine.to_csv(f"results/most_affine_{run_name}.csv", index=False)        


if __name__ == "__main__":
    args = args_parser()
    print("")
    print(f"Generating docking input for {args.run_name}")

    generate_input(f"protein_structures/{args.run_name}", args.SMILES, args.run_name)
