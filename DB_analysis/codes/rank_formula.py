import numpy as np
import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(
    "--results_file",
    type=str,
    help="Path to the file containing the DD results. Must be csv",
    required=True,
)
parser.add_argument(
    "--output_file",
    type=str,
    help="File that will be saved in output directory and will contain the predictions of the DD run.",
    required=True,
)
parser.add_argument(
    "--data_a",
    type=str,
    help="Path to the file containing the proteins and ligands with their affinities of proteins A (more affine). Must be csv",
    required=True,
)
parser.add_argument(
    "--data_b",
    type=str,
    help="Path to the file containing the proteins and ligands with their affinities of proteins B (less affine). Must be csv",
    required=True,
)
args = parser.parse_args()


def prediction(df, out_file):
    """
    This function takes a pandas dataframe as input and predicts the chain in which the ligand has attached.
    The prediction is based on the confidence scores of each prediction chain.

    Args:
        df (pd.DataFrame): A pandas dataframe containing the prediction results and customer information.
                            The columns of the dataframe must include:
                                'ID', 'Confidence', 'Chain(A=0)(B=1)'

    Returns:
        pd.DataFrame: A pandas dataframe containing the chain ID, prediction, and Reliability of the prediction.
    """

    pred_df = pd.DataFrame(columns=["ID", "Prediction", "Reliability"])

    for i in np.unique(df["ID"]):

        slice = df[df["ID"] == i].copy()
        slice = slice[slice["Confidence"] != -1000]  # exclude failed attempts
        samples = len(slice)
        slice = slice.reset_index(drop=True)

        slice.loc[slice["Chain(A=0)(B=1)"] == 0, "Chain(A=0)(B=1)"] = -1
        # slice.loc[slice['Chain(A=0)(B=1)'] == 1, 'Chain(A=0)(B=1)'] = 0
        # slice.loc[slice['Chain(A=0)(B=1)'] == -1, 'Chain(A=0)(B=1)'] = 1
        # slice.loc[slice['Chain(A=0)(B=1)'] == 0, 'Chain(A=0)(B=1)'] = -1

        # print(slice)
        ci = np.array(slice["Confidence"])
        # ci[ci > -2] *= 2
        # print(ci)
        c_max = ci[0]
        c_min = ci[len(slice) - 1]

        numerator = ci - c_min
        denominator = c_max - c_min

        norm_conf = np.array(numerator) / denominator

        slice["Normal conf"] = norm_conf

        pred_1 = np.sum(
            np.array(slice["Chain(A=0)(B=1)"]) * np.array(slice["Normal conf"])
        )
        pred_1 = pred_1 / samples
        # if slice.shape[0] > samples*2/3:
        
        if pred_1 > 0.0:
            pred_df = pd.concat(
                [
                    pred_df if not pred_df.empty else None,
                    pd.DataFrame(
                        [
                            {
                                "ID": i,
                                "Prediction": 1,
                                "Reliability": round(abs(pred_1) * 100, 2),
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )
            # print(f'The prediction is chain 1 with certainty of {pred_1*100:.2f}%')
        else:
            pred_df = pd.concat(
                [
                    pred_df if not pred_df.empty else None,
                    pd.DataFrame(
                        [
                            {
                                "ID": i,
                                "Prediction": 0,
                                "Reliability": round(abs(pred_1) * 100, 2),
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )
            # print(f'The prediction is chain 0 with certainty of {(1-pred_1)*100:.2f}%')
        # print(pred_1)

        pred_df.to_csv("output/" + out_file)
    return pred_df


def plot_rel_rat(dfa, dfb, pred_df, out):
    import matplotlib.pyplot as plt

    # if DD has failed there will be different values at the dfa and dfb so we must avoid them
    for i in range(len(dfa["Prot ID"])):
        id = f"{dfa['Prot ID'].loc[i]}_{dfb['Prot ID'].loc[i]}_{int(dfa['PubChem CID'].loc[i])}"
        # print(id)
        if id in pred_df["ID"].values:
            pass
        else:
            dfa = dfa.drop(i)
            dfb = dfb.drop(i)

    # pred_df = pred_df[pred_df['Prediction'] == 1]
    x = np.array(pred_df["Reliability"])
    y = np.log10(np.array(dfa["Ki (nM)"]) / np.array(dfb["Ki (nM)"]))

    col = np.array(pred_df["Prediction"], dtype=float)
    col[col == 0] = 0.20
    correct = len(col[col == 0.2])
    failed = len(col) - correct

    if len(y) != len(col):
        print(
            f'{out.split(".")[0]}_rel.pdf will not be generated since the prediction file and the data_A and data_B are not corresponding: {dfa.shape}, {len(col)}'
        )
        return

    plt.figure(figsize=(20, 15))
    a = plt.scatter(x, y, c=col, cmap="PiYG_r", s=35, clim=(0, 1))
    plt.xlabel("Reliability", fontsize=35)
    plt.ylabel(rf"$log(K_i^A/K_i^B)$", fontsize=35)
    plt.legend(
        fontsize=35,
        markerscale=4.0,
        handles=a.legend_elements()[0],
        labels=[f"Guessed: {correct}", f"Failed: {failed}"],
        ncol=2,
        bbox_to_anchor=(0.5, 1.05),
        fancybox=True,
        shadow=True,
        loc="upper left",
    )

    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)
    plt.tight_layout()
    plt.savefig(f"output/figures/{out.split('.')[0]}_rel.pdf")
    # plt.show()


if __name__ == "__main__":
    print(f"Ranking {args.results_file.split('/')[-1]} values...")
    df = pd.read_csv(args.results_file)
    dfa = pd.read_csv(args.data_a)
    dfb = pd.read_csv(args.data_b)

    df = prediction(df, args.output_file)
    plot_rel_rat(dfa, dfb, df, args.output_file)
