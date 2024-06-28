import pandas as pd
import numpy as np
import os
from argparse import ArgumentParser

def args_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--ds",
        type=str,
        help="list containing the index of the first and last aa of the docking site of the protein.",
        required=True,
    )
    parser.add_argument(
        "--n_mutations",
        type=int,
        help="Number of different random mutations to be performed at the docking site.",
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

def mutate(run_name, prot_name, seq, ds, n_mutations=1):

    try:
        os.mkdir(f'protein_structures/{run_name}/mutations')
    except:
        pass

    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
          'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    seq = list(seq) # list all the aa of the sequence

    ds = [int(i)-1 for i in ds.split(',')] # Change aa position for its index

    mutated_sequences = []
    ids = []

    if n_mutations > 20*len(ds):
        print(f'the maximum number of mutations possible is {20*len(ds)}')
        return

    for i in range(n_mutations):

        check = False

        while check == False:
            aux = seq.copy()

            r1 = np.random.choice(ds)
            
            r2 = np.random.randint(20)

            prev_aa = aux[r1]

            aux[r1] = aa[r2]

            
            new_seq = ''.join(aux)
            
            if new_seq not in mutated_sequences:
                mutated_sequences.append(new_seq)
                ids.append(f'{prev_aa}_{r1+1}_{aux[r1]}')
                check = True
    
    df = pd.DataFrame()
    df['id'] = ids
    df['sequence'] = mutated_sequences
    df.to_csv(f'protein_structures/{run_name}/mutations/mutations_{prot_name}.csv', index=False)

    return df

def mutate_all(seq_file,run_name,ds,n_mutations=3):
    df = pd.read_csv(seq_file)

    for i in range(len(df)):
        prot_name = df.iloc[i]['id']
        seq = df.iloc[i]['sequence']

        mutate(run_name,prot_name,seq,ds,n_mutations=n_mutations)  


if __name__ == '__main__':
    args = args_parser()
    print('')
    print(f'Mutating proteins from {args.run_name}')
    # mutate(args.run_name,args.seq,args.ds,args.n_mutations)

    ds = args.ds # 
    mutate_all(f'inputs/folding_{args.run_name}_input.csv',args.run_name, args.ds,args.n_mutations)

