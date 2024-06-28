import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
from argparse import ArgumentParser
import seaborn as sns
import os
import matplotlib.patches as mpatches


def args_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--run_name",
        type=str,
        help="Name of the job.",
        required=True,
    )
    parser.add_argument(
        "--AF",
        type=int,
        help="whether AF was used or not. Default is 1 (true)",
        required=False,
        default=1,
    )
    parser.add_argument(
        "--DF",
        type=int,
        help="whether DF was used or not. Default is 1 (true)",
        required=False,
        default=1,
    )
    parser.add_argument(
        "--OF",
        type=int,
        help="whether OF was used or not. Default is 1 (true)",
        required=False,
        default=1,
    )
    args = parser.parse_args()
    return args


def plot_best_chain(run_name, fold="all"):

    df = pd.read_csv(f"results/ranked_{run_name}_{fold}.csv")
    ids = df["ID"].values

    prot1 = []
    prot2 = []

    for i in ids:
        id = i.split("_")
        prot1.append(id[0])
        prot2.append(id[1])

    df["Prot 1"] = prot1
    df["Prot 2"] = prot2

    c_list = cm.coolwarm(np.linspace(0, 1, 2))

    data = df.pivot(index="Prot 2", columns="Prot 1", values="Rank")
    data = data.iloc[
        data.isnull().sum(axis=1).mul(-1).argsort()
    ]  # make sure the data is sorted

    plot_num = np.sum(df["Rank"].values * np.where(df["Affine"].values == 0, 1, 0)) / (
        (data.shape[0] + 1) * data.shape[0] / 2
    )

    plt.figure()
    hm = sns.heatmap(
        data,
        cmap="coolwarm",
        cbar_kws={"label": "Rank"},
        linewidth=1.5,
        annot=True,
        vmin=0,
        vmax=1,
    )
    bbox = dict(boxstyle="round", ec="black", fc=c_list[0], alpha=0.7)
    plt.setp(hm.get_xticklabels(), bbox=bbox)
    bbox = dict(boxstyle="round", ec="black", fc=c_list[1], alpha=0.7)
    plt.setp(hm.get_yticklabels(), bbox=bbox)

    plt.title(f"{run_name} ({fold}) {plot_num:.3f}")

    plt.tight_layout()
    plt.savefig(f"results/figures/rank_{run_name}_{fold}.png")
    # plt.show()


def plot_best_chain_guess(run_name, fold="all"):

    df = pd.read_csv(f"results/ranked_{run_name}_{fold}.csv")
    ids = df["ID"].values

    prot1 = []
    prot2 = []

    for i in ids:
        id = i.split("_")
        prot1.append(id[0])
        prot2.append(id[1])

    df["Prot 1"] = prot1
    df["Prot 2"] = prot2

    most_affine = pd.read_csv(f"results/most_affine_{run_name}.csv")
    df["Affine"] = np.where(
        most_affine["Affine"].values - np.round(df["Rank"].values) == 0, 0, 1
    )
    c_list = cm.coolwarm(np.linspace(0, 1, 2))
    c_list2 = cm.PiYG(np.linspace(0, 1, 2))

    data = df.pivot(index="Prot 2", columns="Prot 1", values="Rank")
    data = data.iloc[
        data.isnull().sum(axis=1).mul(-1).argsort()
    ]  # make sure the data is sorted

    plot_num = np.sum(df["Rank"].values * np.where(df["Affine"].values == 0, 1, 0)) / (
        (data.shape[0] + 1) * data.shape[0] / 2
    )

    data_aff = df.pivot(index="Prot 2", columns="Prot 1", values="Affine")
    data_aff = data_aff.iloc[data_aff.isnull().sum(axis=1).mul(-1).argsort()]

    cmap1 = plt.get_cmap("coolwarm").copy()
    cmap1.set_under("none")
    cmap2 = plt.get_cmap("PiYG").copy()
    cmap2.set_over("none")

    width = data.shape[0]
    height = data.shape[1]
    xpos, ypos = np.meshgrid(np.arange(width) + 0.5, np.arange(height) + 0.5)
    xpos = np.where(~np.isnan(data), xpos, np.nan)
    ypos = np.where(~np.isnan(data), ypos, np.nan)

    plt.figure()
    hm = sns.heatmap(
        data,
        cmap=cmap1,
        cbar_kws={"label": "Rank"},
        linewidth=1.5,
        annot=True,
        vmin=0,
        vmax=1,
    )

    if (data.shape[1] - 3) >= 0:
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if not np.isnan(xpos[i][j]):
                    if data_aff.iloc[i, j] == 0:
                        plt.text(
                            xpos[i][j] + 0.3,
                            ypos[i][j] - 0.3,
                            "   ",
                            size=10 - 2 * (data.shape[1] - 3),
                            bbox=dict(
                                boxstyle="round",
                                ec="black",
                                fc=c_list2[1],
                            ),
                        )
                    else:
                        plt.text(
                            xpos[i][j] + 0.3,
                            ypos[i][j] - 0.3,
                            "   ",
                            size=10 - 2 * (data.shape[1] - 3),
                            bbox=dict(
                                boxstyle="round",
                                ec="black",
                                fc=c_list2[0],
                            ),
                        )
    else:
        if data.shape[1] == 2:
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    if not np.isnan(xpos[i][j]):
                        if data_aff.iloc[i, j] == 0:
                            plt.text(
                                xpos[i][j] + 0.3,
                                ypos[i][j] - 0.3,
                                "   ",
                                size=10,
                                bbox=dict(
                                    boxstyle="round",
                                    ec="black",
                                    fc=c_list2[1],
                                ),
                            )
                        else:
                            plt.text(
                                xpos[i][j] + 0.3,
                                ypos[i][j] - 0.3,
                                "   ",
                                size=10,
                                bbox=dict(
                                    boxstyle="round",
                                    ec="black",
                                    fc=c_list2[0],
                                ),
                            )

    plt.text(
        xpos[data.shape[1] - 1][data.shape[1] - 1] + 0.2,
        xpos[data.shape[1] - 1][data.shape[1] - 1] - data.shape[1] + 1 - 0.3,
        f"{len(df[df['Affine']==1])}",
        size=15,
        bbox=dict(
            boxstyle="round",
            ec="black",
            fc=c_list2[0],
        ),
        color="white",
    )
    plt.text(
        xpos[data.shape[1] - 1][data.shape[1] - 1] + 0.2,
        xpos[data.shape[1] - 1][data.shape[1] - 1]
        - data.shape[1]
        + 1
        + 0.1 * (data.shape[1] - 3),
        f"{len(df[df['Affine']==0])}",
        size=15,
        bbox=dict(
            boxstyle="round",
            ec="black",
            fc=c_list2[1],
        ),
        color="white",
    )

    bbox = dict(boxstyle="round", ec="black", fc=c_list[0], alpha=0.7)
    plt.setp(hm.get_xticklabels(), bbox=bbox)
    bbox = dict(boxstyle="round", ec="black", fc=c_list[1], alpha=0.7)
    plt.setp(hm.get_yticklabels(), bbox=bbox)

    plt.title(f"{run_name} ({fold}) {plot_num:.3f}")

    plt.tight_layout()
    plt.savefig(f"results/figures/rank_{run_name}_{fold}.png")
    # plt.show()


def plot_logk(run_name, fold):

    c_list = cm.coolwarm(np.linspace(0, 1, 2))

    df = pd.read_csv(f"results/ranked_{run_name}_{fold}.csv")
    ids = df["ID"].values
    most_affine = pd.read_csv(f"results/most_affine_{run_name}.csv")

    ratios = np.log10(most_affine["K2/K1"].values)
    rank = df["Rank"].values
    rank = np.where(rank > 0.5, rank, 1 - rank)

    df["Affine"] = np.where(
        most_affine["Affine"].values - np.round(df["Rank"].values) == 0, 0, 1
    )

    plt.figure()
    plt.scatter(rank, ratios, cmap="PiYG_r", c=df["Affine"].values)
    plt.hlines(
        0.05, min(rank) - 0.2, max(rank) + 0.2, linestyles="--", colors=c_list[0]
    )
    plt.hlines(
        -0.05, min(rank) - 0.2, max(rank) + 0.2, linestyles="--", colors=c_list[0]
    )

    plt.ylabel("log(K2/K1)")
    plt.xlabel("Rank")
    plt.title(f"{run_name} ({fold})")

    plt.xlim([min(rank) - 0.05, max(rank) + 0.05])
    plt.ylim([min(ratios) - 0.2, max(ratios) + 0.2])

    plt.tight_layout()
    plt.savefig(f"results/figures/logk_{run_name}_{fold}.png")

def plot_single_affinity(run_name, fold):
    df_results = pd.read_csv(f"results/docking_{run_name}_{fold}.csv")
    confidence = df_results.groupby('ID')['Confidence'].apply(lambda x: x.head(15).mean()).reset_index()
    df_affinity = pd.read_csv(f"results/most_affine_{run_name}.csv")
    # print(df_affinity)
    affinity = df_affinity['Kd']

    df_plot = pd.concat([confidence, affinity], axis=1)
    df_plot = df_plot.sort_values(by='Kd')
    df_plot = df_plot.dropna()
    df_plot = df_plot.loc[df_plot['Kd'] != 1].loc[df_plot['Confidence'] > -10]
    df_plot = df_plot.iloc[:]
    print(df_plot.shape)
    # print(df_plot)
    x = df_plot['Kd']
    y = df_plot['Confidence']
    coefficients = np.polyfit(x, y, 1)
    trendline = np.poly1d(coefficients)
    # print(x,y)
    # Calcula el R^2
    residuals = y - trendline(x)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y-np.mean(y))**2)
    r2 = 1 - (ss_res / ss_tot)

    plt.figure()
    plt.plot((x), y, 'o')
    plt.plot((x), trendline(x), '-')
    plt.xlabel('Kd')
    plt.ylabel('DD Confidence')
    plt.title(f'Confidence correlation - {run_name} {fold}')
    plt.text(0.82, 1.05, r'$R^2$'f' = {r2:.3f}', transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))    
    plt.tight_layout()
    plt.savefig(f'results/figures/confidence_{run_name}_{fold}.png')




if __name__ == "__main__":
    args = args_parser()
    print(f"Plotting ranked predictions from {args.run_name}")

    if args.AF == 1:
        if os.path.exists(f"results/most_affine_{args.run_name}.csv"):
            plot_single_affinity(args.run_name, "AF")

    if args.DF == 1:
        if os.path.exists(f"results/most_affine_{args.run_name}.csv"):
            plot_single_affinity(args.run_name, "DF")

    if args.OF == 1:
        if os.path.exists(f"results/most_affine_{args.run_name}.csv"):
            plot_single_affinity(args.run_name, "OF")

    # if os.path.exists(f"results/most_affine_{args.run_name}.csv"):
    #     plot_single_affinity(args.run_name, "all")
