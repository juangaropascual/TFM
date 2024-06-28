import numpy as np
import pandas as pd
from argparse import ArgumentParser
from joblib import dump, load
from sklearn.preprocessing import StandardScaler,PolynomialFeatures
# import detect_DD_fail as ddfail

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
parser.add_argument(
    "--model_name",
    type=str,
    help="Name of the model that will be saved in a .joblib file.",
    required=True,
)
parser.add_argument(
    "--n_samples",
    type=int,
    help="The number of different samples from which the features will be extracted.",
    required=True,
)
parser.add_argument(
    "--validation_data",
    type=str,
    help="Name of file containing the set for validation.",
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
        slice = slice.reset_index(drop=True)
        # slice = slice.iloc[:80]
        samples = len(slice)

        slice.loc[slice["Chain(A=0)(B=1)"] == 0, "Chain(A=0)(B=1)"] = -1
        slice.loc[slice['Chain(A=0)(B=1)'] == 1, 'Chain(A=0)(B=1)'] = 0
        slice.loc[slice['Chain(A=0)(B=1)'] == -1, 'Chain(A=0)(B=1)'] = 1
        # slice.loc[slice['Chain(A=0)(B=1)'] == 0, 'Chain(A=0)(B=1)'] = -1  

        ci = np.array(slice["Confidence"])
        c_max = ci[0]
        c_min = ci[samples - 1]

        numerator = ci - c_min
        denominator = c_max - c_min

        norm_conf = np.array(numerator) / denominator

        slice["Normal conf"] = norm_conf
        d_0 = slice.loc[slice["Chain(A=0)(B=1)"] == 0, "Normal conf"] 
        c_0 = np.sum(np.array(d_0))/samples
   
        d_1 = slice.loc[slice["Chain(A=0)(B=1)"] == 1, "Normal conf"]
        c_1 = np.sum(np.array(d_1))/samples

        # ni/N
        # w_0 = len(d_0)/samples
        # w_1 = len(d_1)/samples

        # sum(ind/samples) / w_max 
        w_max = np.sum(1/np.arange(1,samples+1))
        w_0 = np.sum(1/(np.array(d_0.index) + 1))/w_max
        w_1 = np.sum(1/(np.array(d_1.index) + 1))/w_max

        # simply 1
        # w_0 = 1
        # w_1 = 1

        # change weight every 10 samples
        # w_max = np.sum([(i+1) for i in range(int(samples/10))])
        # ind_0 = np.array(d_0.index) + 1
        # ind_1 = np.array(d_1.index) + 1

        # i_0 = []
        # i_1 = []
        # factor = 10
        # for j in (ind_0):
        #     if j < factor:
        #         i_0.append(1/factor)
        #     else:
        #         factor += 10
        #         i_0.append(1/factor)

        # factor = 10
        # for j in ind_1:
        #     if j < factor:
        #         i_1.append(1/factor)
        #     else:
        #         factor += 10
        #         i_1.append(1/factor) 

        # w_0 = np.sum((np.array(i_0)))/w_max
        # w_1 = np.sum((np.array(i_1)))/w_max
        
        pred_0 = c_0*w_0
        pred_1 = c_1*w_1

        # print(pred_0,pred_1) 
        # print(w_0,w_1, c_0, c_1)   

        if pred_1 > pred_0:
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
                                "Reliability": round(abs(pred_0) * 100, 2),
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

def prediction_rank1(df, out_file):
    pred_df = pd.DataFrame(columns=["ID", "Prediction", "Reliability"])

    for i in np.unique(df["ID"]):

        slice = df[df["ID"] == i].copy()
        slice = slice[slice["Confidence"] != -1000]  # exclude failed attempts
        slice = slice.reset_index(drop=True)
        # slice = slice.iloc[:80]
        samples = len(slice)

        slice.loc[slice["Chain(A=0)(B=1)"] == 0, "Chain(A=0)(B=1)"] = -1
        slice.loc[slice['Chain(A=0)(B=1)'] == 1, 'Chain(A=0)(B=1)'] = 0
        slice.loc[slice['Chain(A=0)(B=1)'] == -1, 'Chain(A=0)(B=1)'] = 1
        # slice.loc[slice['Chain(A=0)(B=1)'] == 0, 'Chain(A=0)(B=1)'] = -1  

        pred = slice["Chain(A=0)(B=1)"].loc[0]
        c = slice['Confidence'].loc[0] 

        if pred == 1:
            pred_df = pd.concat(
                [
                    pred_df if not pred_df.empty else None,
                    pd.DataFrame(
                        [
                            {
                                "ID": i,
                                "Prediction": 1,
                                "Reliability": round(abs(c), 2),
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
                                "Reliability": round(abs(c), 2),
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )
        pred_df.to_csv("output/" + out_file)
    return pred_df

def rank_by_forest_model(model,validation_data,out_file):
    # load the model from the saved file
    clf = load(f'model/{model}') 

    # load the validation data
    val_df = pd.read_csv(f'model/{validation_data}')

    ids = val_df['ID']
    val_df = val_df.drop('ID',axis=1)
    exp_pred = val_df['Experimental Pred']
    val_df = val_df.drop('Experimental Pred',axis=1)
    prot_a = []
    prot_b = []
    ligs = []
    
    y_pred = clf.predict(val_df.values)

    pred_df = pd.DataFrame()

    for i,id in enumerate(ids):
        id = id.split('_')
        if (exp_pred[i] == 'A' or exp_pred[i] == 'C'):
            prot_a.append(id[0])
            prot_b.append(id[1])
        elif exp_pred[i] == 'B':
            prot_a.append(id[1])
            prot_b.append(id[0])
        ligs.append(id[2])

    pred_df['Prot A ID'] = prot_a
    pred_df['Prot B ID'] = prot_b
    pred_df['Ligand ID'] = ligs

    pred_df['Experimental Pred'] = exp_pred
    pred_df['Prediction'] = y_pred
    # Change 0 and 1 for A and B on the selected chains.
    pred_df.loc[pred_df["Prediction"] == 'A', "Prediction"] = 0
    pred_df.loc[pred_df["Prediction"] == 'B', "Prediction"] = 1
    pred_df.loc[pred_df["Prediction"] == 'C', "Prediction"] = 0.7

    pred_df.loc[pred_df["Experimental Pred"] == 'A', "Experimental Pred"] = 0
    pred_df.loc[pred_df["Experimental Pred"] == 'B', "Experimental Pred"] = 1
    pred_df.loc[pred_df["Experimental Pred"] == 'C', "Experimental Pred"] = 0.7

    ac = clf.score(val_df.values, exp_pred)
    print(f'Score of the model: {ac}')
    pred_df.to_csv(f'output/{out_file}')

    return pred_df

def rank_by_svc_model(model,validation_data,out_file):
    # load the model from the saved file
    clf = load(f'model/{model}') 

    # load the validation data
    val_df = pd.read_csv(f'model/{validation_data}')

    ids = val_df['ID']
    val_df = val_df.drop('ID',axis=1)
    exp_pred = val_df['Experimental Pred']
    val_df = val_df.drop('Experimental Pred',axis=1)
    prot_a = []
    prot_b = []
    ligs = []
    X_val = val_df.values
    # Make data separable by using polynomial behaviour
    # poly = PolynomialFeatures(3)
    # X_val = poly.fit_transform(val_df.values)

    # Scale the data 
    feature_scaler = StandardScaler()
    X_val = feature_scaler.fit_transform(X_val)
    y_pred = clf.predict(X_val)

    pred_df = pd.DataFrame()

    for i,id in enumerate(ids):
        id = id.split('_')
        if (exp_pred[i] == 'A' or exp_pred[i] == 'C'):
            prot_a.append(id[0])
            prot_b.append(id[1])
        elif exp_pred[i] == 'B':
            prot_a.append(id[1])
            prot_b.append(id[0])
        ligs.append(id[2])

    pred_df['Prot A ID'] = prot_a
    pred_df['Prot B ID'] = prot_b
    pred_df['Ligand ID'] = ligs

    pred_df['Experimental Pred'] = exp_pred
    pred_df['Prediction'] = y_pred

    # Change 0 and 1 for A and B on the selected chains.
    pred_df.loc[pred_df["Prediction"] == 'A', "Prediction"] = 0
    pred_df.loc[pred_df["Prediction"] == 'B', "Prediction"] = 1
    pred_df.loc[pred_df["Prediction"] == 'C', "Prediction"] = 0.7

    pred_df.loc[pred_df["Experimental Pred"] == 'A', "Experimental Pred"] = 0
    pred_df.loc[pred_df["Experimental Pred"] == 'B', "Experimental Pred"] = 1
    pred_df.loc[pred_df["Experimental Pred"] == 'C', "Experimental Pred"] = 0.7

    ac = clf.score(X_val, exp_pred)
    print(f'Score of the model: {ac:.3f}')
    pred_df.to_csv(f'output/{out_file}')

    return pred_df

    




if __name__ == "__main__":
    print(f"Ranking {args.results_file.split('/')[-1]} values...")
    df = pd.read_csv(args.results_file)
    dfa = pd.read_csv(args.data_a)
    dfb = pd.read_csv(args.data_b)

    # df = prediction_rank1(df, args.output_file)
    # df = rank_by_forest_model(args.model_name,'validation_data.csv',args.output_file)
    df =rank_by_svc_model(args.model_name,'validation_data.csv',args.output_file)