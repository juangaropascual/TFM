import numpy as np
import pandas as pd
from argparse import ArgumentParser
import os


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
    parser.add_argument(
        "--n_samples",
        type=int,
        help="How many samples are being considered for ranking",
        required=False,
        default=1,
    )
    args = parser.parse_args()
    return args


def rank_docking_all(run_name, AF, DF, OF, n_samples):
    results = [
        f"docking_{run_name}_AF.csv" * AF,
        f"docking_{run_name}_DF.csv" * DF,
        f"docking_{run_name}_OF.csv" * OF,
    ]

    results = [i for i in results if i != ""]  # extract all results that don't exist

    rank_total = np.array([])
    fold = []
    for i, r in enumerate(results):
        rank = []

        res = pd.read_csv(f"results/{r}")

        # samples = len(res[res['ID'] == res['ID'].iloc[0]])

        ids = np.unique(res["ID"].values)
        if len(rank_total) == 0:
            rank_total = np.empty(len(ids))

        for id in ids:
            rank.append(
                np.sum(res[res["ID"] == id]["Chain(A=0)(B=1)"].values[:n_samples])
                / n_samples
            )
            fold.append(r.split(".")[0].split("_")[-1])
        # for i in range(len(rank)):
        #     print(ids[i],rank[i])

        rank_total += np.array(rank)

    ranked = pd.DataFrame()

    ranked["ID"] = ids
    ranked["Rank"] = np.round(rank_total / len(results), 4)
    ranked["Folding"] = ["all" for i in range(len(ids))]

    ranked.to_csv(f"results/ranked_{run_name}_all.csv", index=False)


def rank_one(run_name, n_samples, fold="AF"):
    results = f"results/docking_{run_name}_{fold}.csv"

    res = pd.read_csv(results)
    # samples = len(res[res['ID'] == res['ID'].iloc[0]])

    ids = np.unique(res["ID"].values)

    rank = []
    for id in ids:
        rank.append(
            np.sum(res[res["ID"] == id]["Chain(A=0)(B=1)"].values[:n_samples])
            / n_samples
        )

    # for i in range(len(rank)):
    #     print(ids[i],rank[i])

    ranked = pd.DataFrame()

    ranked["ID"] = ids
    ranked["Rank"] = rank
    ranked["Folding"] = [fold for i in range(len(ids))]

    ranked.to_csv(f"results/ranked_{run_name}_{fold}.csv", index=False)


if __name__ == "__main__":
    args = args_parser()
    print(f" Ranking docking results from {args.run_name}...")

    if args.AF == 1:
        rank_one(args.run_name, args.n_samples)
    if args.DF == 1:
        rank_one(args.run_name, args.n_samples, fold="DF")
    if args.OF == 1:
        rank_one(args.run_name, args.n_samples, fold="OF")

    rank_docking_all(args.run_name, args.AF, args.DF, args.OF, args.n_samples)
