import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
from argparse import ArgumentParser
import os

parser = ArgumentParser()
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
parser.add_argument(
    "--threshold",
    type=float,
    help="Threshold from which the data has been studied. Default is 0.1",
    default=0.1,
    required=False,
)
parser.add_argument(
    "--predictions",
    type=str,
    help="Path to file containing the DD predictions.",
    required=True,
)
parser.add_argument(
    "--counts",
    type=bool,
    help="Whether to plot or not the counts graph. Default is False",
    required=False,
    default=False,
)
parser.add_argument(
    "--failed_file",
    type=str,
    help="File containing if the DiffDock run seems to fave failed when only running one protein with one ligand. Default is None",
    required=False,
    default=None,
)
parser.add_argument(
    "--plot_th",
    type=int,
    help="Whether or not to only plot the points below the threshold. Default is False.",
    required=False,
    default=0,
)
args = parser.parse_args()


def plot_ab_count(dfa, dfb, threshold=0.1):
    """
    This function plots the dataframes dfa and dfb, where dfa contains the affinities of proteins A with the same ligand, and dfb contains the affinities of proteins
    B with the same ligand.

    Args:
        dfa (pd.DataFrame): The dataframe containing the affinities of proteins A with the same ligand.
        dfb (pd.DataFrame): The dataframe containing the affinities of proteins B with the same ligand.
        threshold (float, optional): The threshold used for sorting. Defaults to 0.1.
        counts_graph (bool, optional): If True, a histogram of the ligand counts is plotted alongside the error bars. Defaults to False.
        filename (str, optional): The name of the file containing the DD guess. Defaults to None.

    Returns:
        None

    """

    un_lig, counts = np.unique(dfa["SMILES"], return_counts=True)
    col = []

    for lig in dfa["SMILES"]:
        col.append(counts[list(un_lig).index(lig)])

    c_list = cm.viridis(np.linspace(0, 1, 4))

    x = np.array(dfb["Kd (nM)"])
    y = np.array(dfa["Kd (nM)"])

    err_x = np.abs(np.array(dfb["kd SEM"]) / x)
    err_y = np.abs(np.array(dfa["kd SEM"]) / y)

    mask_x = err_x != 0
    mask_y = err_y != 0

    mask = mask_x | mask_y
    err_x = err_x[mask]
    err_y = err_y[mask]

    x_err = x[mask]
    y_err = y[mask]

    x = np.log10(x)
    y = np.log10(y)

    x_err = np.log10(x_err)
    y_err = np.log10(y_err)

    x_values = np.linspace(min(min(x), min(y)) - 5, max(max(x), max(y) + 5), 100)

    plt.figure(figsize=(20, 15))

    plt.errorbar(
        x_err,
        y_err,
        yerr=err_y,
        xerr=err_x,
        alpha=0.3,
        ecolor="r",
        capsize=5,
        fmt="none",
        zorder=0,
    )
    plt.scatter(x, y, c=col, cmap="viridis", s=35)

    plt.plot(x_values, x_values, label="$log(K_D^A) = log(K_D^B)$", color=c_list[0])
    plt.plot(
        x_values,
        x_values + np.log10(threshold),
        color=c_list[1],
        label=f"$log(K_D^A/K_D^B) = log({threshold}$)",
    )

    plt.xlabel(rf"$log(K_D^B)$", fontsize=35)
    plt.ylabel(rf"$log(K_D^A)$", fontsize=35)
    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=35)
    cbar.set_label("Ligand counts", fontsize=35)
    plt.legend(fontsize=35)

    try:
        min_value = min(min(min(x), min(x_err)), min(min(y), min(y_err))) - 0.2
        max_value = max(max(max(x), max(x_err)), max(max(y), max(y_err))) + 0.2
    except ValueError:
        min_value = min(min(x), min(y)) - 0.2
        max_value = max(max(x), max(y)) + 0.2

    plt.xlim([min_value, max_value])
    plt.ylim([min_value, max_value])
    plt.tight_layout()

    plt.savefig("output/figures/ka_kb_count.png")
    # plt.show()


def plot_bool_res(dfa, dfb, threshold=0.1, filename=None, plot_th=0):
    """
    This function plots the dataframes dfa and dfb, where dfa contains the affinities of proteins A with the same ligand, and dfb contains the affinities of proteins
    B with the same ligand.

    Args:
        dfa (pd.DataFrame): The dataframe containing the affinities of proteins A with the same ligand.
        dfb (pd.DataFrame): The dataframe containing the affinities of proteins B with the same ligand.
        threshold (float, optional): The threshold used for sorting. Defaults to 0.1.
        counts_graph (bool, optional): If True, a histogram of the ligand counts is plotted alongside the error bars. Defaults to False.
        filename (str, optional): The name of the file containing the DD guess. Defaults to None.

    Returns:
        None

    """

    col = pd.read_csv(filename)
    c_list = cm.RdYlGn(np.linspace(0, 1, 2))

    # if DD has failed there will be different values at the dfa and dfb so we must avoid them
    for i in range(len(dfa["Prot ID"])):
        id = f"{dfa['Prot ID'].loc[i]}_{dfb['Prot ID'].loc[i]}_{int(dfa['PubChem CID'].loc[i])}"

        if id in col["ID"].values:
            pass
        else:
            dfa = dfa.drop(i)
            dfb = dfb.drop(i)
    dfa = dfa.reset_index(drop=True)
    dfb = dfb.reset_index(drop=True)

    if threshold != 1:
        mask = dfa[dfa["Kd (nM)"] / dfb["Kd (nM)"] >= threshold]
        col_th = col.drop(mask.index)

        col_th = np.array(col_th["Prediction"], dtype=float)
        col_th[col_th == 0] = 0.2
        correct_th = len(col_th[col_th == 0.2])
        failed_th = len(col_th) - correct_th

        if plot_th == 1:
            dfa = dfa.drop(mask.index)
            dfb = dfb.drop(mask.index)
            col = col.drop(mask.index)

    col = np.array(col["Prediction"], dtype=float)
    col[col == 0] = 0.2
    correct = len(col[col == 0.2])
    failed = len(col) - correct

    x = np.array(dfb["Kd (nM)"])
    y = np.array(dfa["Kd (nM)"])

    if len(x) != len(col):
        print(
            f'{filename.split("/")[-1].split(".")[0]}_kakb_pred.png will not be generated since the prediction file and the data_A and data_B are not corresponding:{dfa.shape}, {len(col)}'
        )
        return
    err_x = np.array(dfb["kd SEM"]) / x
    err_y = np.array(dfa["kd SEM"]) / y

    mask_x = err_x != 0
    mask_y = err_y != 0

    mask = mask_x | mask_y
    err_x = err_x[mask]
    err_y = err_y[mask]

    x_err = x[mask]
    y_err = y[mask]

    x = np.log10(x)
    y = np.log10(y)

    x_err = np.log10(x_err)
    y_err = np.log10(y_err)

    x_values = np.linspace(min(min(x), min(y)) - 5, max(max(x), max(y) + 5), 100)

    plt.figure(figsize=(15, 15))

    plt.errorbar(
        x_err, y_err, yerr=err_y, xerr=err_x, alpha=0.3, capsize=5, fmt="none", ecolor='gray', zorder=0
    )
    a = plt.scatter(x, y, c=col, cmap="RdBu", s=35, clim=(0, 1))  # 'RdYlGn'

    plt.plot(x_values, x_values, label="$log(K_D^A) = log(K_D^B)$", color='black')
    plt.plot(
        x_values,
        x_values + np.log10(threshold),
        color='gray',
        label=f"$log(K_D^A/K_D^B) = log({threshold}$)",
    )

    plt.xlabel(rf"$log(K_D^B)$", fontsize=35)
    plt.ylabel(rf"$log(K_D^A)$", fontsize=35)
    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)
    # cbar = plt.colorbar()
    # cbar.ax.tick_params(labelsize=35)

    if plot_th == 1:
        legend = plt.legend(
            fontsize=35,
            markerscale=4.0,
            handles=a.legend_elements()[0],
            labels=[f"Guessed: {correct}", f"Failed: {failed}"],
        )
        plt.gca().add_artist(legend)
        legend_l =plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=2, fontsize=35)
    else:
        legend = plt.legend(
            fontsize=35,
            markerscale=4.0,
            handles=a.legend_elements()[0],
            labels=[f"Guessed: {correct}", f"Failed: {failed}"],
        )
        plt.gca().add_artist(legend)
        if threshold != 1:
            legend2 = plt.legend(
                fontsize=35,
                markerscale=4.0,
                handles=a.legend_elements()[0],
                labels=[f"Guessed (th): {correct_th}", f"Failed (th): {failed_th}"],
                bbox_to_anchor=(0.528, 0.85),
            )
            plt.gca().add_artist(legend2)
            legend_l =plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=2, fontsize=35)
    # cbar.set_label("DiffDock guess", fontsize=35)
    # cbar.set_ticks([0, 1])
    # cbar.set_ticklabels(["A", "B"])

    try:
        min_value = min(min(min(x), min(x_err)), min(min(y), min(y_err))) - 0.2
        max_value = max(max(max(x), max(x_err)), max(max(y), max(y_err))) + 0.2
    except ValueError:
        min_value = min(min(x), min(y)) - 0.2
        max_value = max(max(x), max(y)) + 0.2

    plt.xlim([min_value, max_value])
    plt.ylim([min_value, max_value])

    plt.tight_layout()
    plt.savefig(f'output/figures/{filename.split("/")[-1].split(".")[0]}_kakb_pred.pdf')
    # plt.show()


def plot_rel_res(dfa, dfb, threshold=0.1, filename=None, plot_th=0):
    """
    This function plots the dataframes dfa and dfb, where dfa contains the affinities of proteins A with the same ligand, and dfb contains the affinities of proteins
    B with the same ligand.

    Args:
        dfa (pd.DataFrame): The dataframe containing the affinities of proteins A with the same ligand.
        dfb (pd.DataFrame): The dataframe containing the affinities of proteins B with the same ligand.
        threshold (float, optional): The threshold used for sorting. Defaults to 0.1.
        counts_graph (bool, optional): If True, a histogram of the ligand counts is plotted alongside the error bars. Defaults to False.
        filename (str, optional): The name of the file containing the DD guess. Defaults to None.

    Returns:
        None

    """
    col = pd.read_csv(filename)
    c_list = cm.cividis_r(np.linspace(0, 1, 2))

    # if DD has failed there will be different values at the dfa and dfb so we must avoid them
    for i in range(len(dfa["Prot ID"])):
        id = f"{dfa['Prot ID'].loc[i]}_{dfb['Prot ID'].loc[i]}_{int(dfa['PubChem CID'].loc[i])}"

        if id in col["ID"].values:
            pass
        else:
            dfa = dfa.drop(i)
            dfb = dfb.drop(i)

    dfa = dfa.reset_index(drop=True)
    dfb = dfb.reset_index(drop=True)

    if threshold != 1:
        mask = dfa[dfa["Kd (nM)"] / dfb["Kd (nM)"] >= threshold]
        col_th = col.drop(mask.index)

        col_th = np.array(col_th["Reliability"], dtype=float)
        correct_th = len(col_th[col_th >= 60])
        failed_th = len(col_th) - correct_th

        if plot_th == 1:
            dfa = dfa.drop(mask.index)
            dfb = dfb.drop(mask.index)
            col = col.drop(mask.index)

    col = np.array(col["Reliability"])
    correct = len(col[col >= 60])
    failed = len(col) - correct

    x = np.array(dfb["Kd (nM)"])
    y = np.array(dfa["Kd (nM)"])

    if len(x) != len(col):
        print(
            f'{filename.split("/")[-1].split(".")[0]}_kakb_rel.png will not be generated since the prediction file and the data_A and data_B are not corresponding: {dfa.shape}, {len(col)}'
        )
        return

    err_x = np.array(dfb["kd SEM"]) / x
    err_y = np.array(dfa["kd SEM"]) / y

    mask_x = err_x != 0
    mask_y = err_y != 0

    mask = mask_x | mask_y
    err_x = err_x[mask]
    err_y = err_y[mask]

    x_err = x[mask]
    y_err = y[mask]

    x = np.log10(x)
    y = np.log10(y)

    x_err = np.log10(x_err)
    y_err = np.log10(y_err)

    x_values = np.linspace(min(min(x), min(y)) - 5, max(max(x), max(y) + 5), 100)

    plt.figure(figsize=(20, 15))

    plt.errorbar(
        x_err, y_err, yerr=err_y, xerr=err_x, alpha=0.3, capsize=5, fmt="none", zorder=0
    )
    a = plt.scatter(
        x, y, c=col, cmap="cividis_r", s=35, clim=(min(col), max(col))
    )  # 'RdYlGn'

    plt.plot(x_values, x_values, label="$log(K_D^A) = log(K_D^B)$", color=c_list[0])
    plt.plot(
        x_values,
        x_values + np.log10(threshold),
        color=c_list[1],
        label=f"$log(K_D^A/K_D^B) = log({threshold})$",
    )

    plt.xlabel(rf"$log(K_D^B)$", fontsize=35)
    plt.ylabel(rf"$log(K_D^A)$", fontsize=35)
    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=35)

    if plot_th == 1:
        legend = plt.legend(
            fontsize=35,
            markerscale=4.0,
            handles=[a.legend_elements()[0][-1], a.legend_elements()[0][0]],
            labels=[f"Reliability > 60%: {correct}", f"Reliability < 60%: {failed}"],
        )
        plt.gca().add_artist(legend)
        legend_l =plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=2, fontsize=35)
        # plt.gca().add_artist(legend_l)
    else:
        legend = plt.legend(
            fontsize=35,
            markerscale=4.0,
            handles=[a.legend_elements()[0][-1], a.legend_elements()[0][0]],
            labels=[f"Reliability > 60%: {correct}", f"Reliability < 60%: {failed}"],
        )
        plt.gca().add_artist(legend)
        if threshold != 1:
            legend2 = plt.legend(
                fontsize=35,
                markerscale=4.0,
                handles=[a.legend_elements()[0][-1], a.legend_elements()[0][0]],
                labels=[
                    f"Reliability > 60% (th): {correct_th}",
                    f"Reliability < 60% (th): {failed_th}",
                ],
                bbox_to_anchor=(0.62, 0.85),
            )
            plt.gca().add_artist(legend2)
            legend_l =plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=2, fontsize=35)

    cbar.set_label("DiffDock guess reliability", fontsize=35)
    # cbar.set_ticks([50, 75, 100])
    # cbar.set_ticklabels(["50%", "75%", "100%"])

    try:
        min_value = min(min(min(x), min(x_err)), min(min(y), min(y_err))) - 0.2
        max_value = max(max(max(x), max(x_err)), max(max(y), max(y_err))) + 0.2
    except ValueError:
        min_value = min(min(x), min(y)) - 0.2
        max_value = max(max(x), max(y)) + 0.2

    plt.xlim([min_value, max_value])
    plt.ylim([min_value, max_value])

    plt.tight_layout()

    plt.savefig(f'output/figures/{filename.split("/")[-1].split(".")[0]}_kakb_rel.png')
    # plt.show()


def plot_failed(dfa, dfb, failed, threshold=0.1, counts=False, predictions=None, plot_th=0):

    try:
        fail_df = pd.read_csv(failed)
    except:
        print(f"Wrong file was inputed: {failed}")
        return

    pred = pd.read_csv(predictions)
    pred = pred[["ID", "Prediction", "Reliability"]]

    # if DD has failed there will be different values at the dfa and dfb so we must avoid them
    order = []
    for i in range(len(dfa["Prot ID"])):
        id = f"{dfa['Prot ID'].loc[i]}_{dfb['Prot ID'].loc[i]}_{int(dfa['PubChem CID'].loc[i])}"
        order.append(id)
        if id in pred["ID"].values:
            pass
        else:
            dfa = dfa.drop(i)
            dfb = dfb.drop(i)

    pred = pred.sort_values(by=["ID"])

    # Sort the predictions to match the dfa and dfb
    pred["ID"] = pd.Categorical(pred["ID"], categories=order, ordered=True)
    pred = pred.sort_values(by="ID")

    pred = pred.reset_index()

    for i in range(len(fail_df["Fail"])):
        f = fail_df.loc[i, "Fail"]
        prot = fail_df.loc[i, "Protein"]
        cid = fail_df.loc[i, "PubChem CID"]

        if f == 1:
            mask_a = dfa["Prot ID"].eq(prot) & dfa["PubChem CID"].eq(cid)
            idx_to_drop = dfa.index[mask_a]
            dfa = dfa.drop(idx_to_drop)
            dfb = dfb.drop(idx_to_drop)
            pred = pred.drop(idx_to_drop)

            mask_b = dfb["Prot ID"].eq(prot) & dfb["PubChem CID"].eq(cid)
            idx_to_drop = dfa.index[mask_b]
            dfa = dfa.drop(idx_to_drop)
            dfb = dfb.drop(idx_to_drop)
            pred = pred.drop(idx_to_drop)

            dfa = dfa.reset_index(drop=True)
            dfb = dfb.reset_index(drop=True)
            pred = pred.reset_index(drop=True)

    pred.to_csv(predictions.split(".")[0] + "_no_fails.csv")

    predictions = predictions.split(".")[0] + "_no_fails.csv"
    if counts:
        plot_ab_count(dfa, dfb, threshold=threshold)
    plot_bool_res(dfa, dfb, threshold=threshold, filename=predictions, plot_th=0)
    plot_rel_res(dfa, dfb, threshold=threshold, filename=predictions, plot_th=0)


def plot_conf_per_sample(n):
    results = os.listdir("DD")

    for r in results:
        df = pd.read_csv("DD/" + r)
        samples = len(df[df["ID"] == df["ID"].values[0]])
        samp_list = np.arange(1, samples + 1)

        conf = df["Confidence"].values
        c_list = cm.viridis(np.linspace(0, 1, 8))
        plt.figure()
        for i in range(n):
            y = conf[i * samples : (i + 1) * samples]
            m = max(y)
            y_best = y[y > m - abs(m / 2)]
            x_best = samp_list[: len(y_best)]
            if m >= 0:
                c = 0
            else:
                c = 4
            plt.plot(samp_list, y, color=c_list[c], alpha=0.7)
            plt.plot(x_best, y_best, color="r")

        plt.xlabel("Sample")
        plt.ylabel("Confidence")
        plt.title(f"DD runs: {n}")
        plt.savefig(
            f'output/figures/{r.split("/")[-1].split(".")[0]}_conf_per_sample_{n}.png'
        )
        # plt.show()


if __name__ == "__main__":
    print(f"Plotting {args.predictions}")
    dfa = pd.read_csv(args.data_a)
    dfb = pd.read_csv(args.data_b)
    # dfa = dfa.loc[[i for i in range(51)]]
    # dfb = dfb.loc[[i for i in range(51)]]

    if args.failed_file != "None":
        print(
            f"Using {args.failed_file} to remove pairs that have individually failed with their ligand..."
        )
        print("")
        plot_failed(
            dfa,
            dfb,
            args.failed_file,
            threshold=args.threshold,
            counts=args.counts,
            predictions=args.predictions,
            plot_th=args.plot_th
        )

    else:
        if args.counts:
            plot_ab_count(dfa, dfb, threshold=args.threshold)
            plot_conf_per_sample(10)
            plot_conf_per_sample(20)
            plot_conf_per_sample(30)
        plot_bool_res(dfa, dfb, threshold=args.threshold, filename=args.predictions, plot_th=args.plot_th)
        # plot_rel_res(dfa, dfb, threshold=args.threshold, filename=args.predictions, plot_th=args.plot_th)
