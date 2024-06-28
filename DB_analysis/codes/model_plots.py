import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
from argparse import ArgumentParser
import os
import seaborn as sns
from sklearn.metrics import confusion_matrix

def violin_confidence(results,samples):
    # create a dataframe with the predictions
    
    res = pd.read_csv(results)
    df = pd.DataFrame(columns=["Correct", "Prediction", "Truth", 'Confidence'])
    df["Prediction"] = res['Chain(A=0)(B=1)'].values[::samples].astype(str)
    df["Truth"] = res['Experimental Pred'].values[::samples].astype(str)
    df['Confidence'] = res['Confidence'].values[::samples]

    # Change 0 and 1 for A and B on the selected chains.
    df.loc[df["Prediction"] == '0', "Prediction"] = 'A'
    df.loc[df["Prediction"] == '1', "Prediction"] = 'B'
    df.loc[df['Truth'] == 'C', 'Truth'] = 'A'


    # order = ["A", "B"]
    # df["Prediction"] = pd.Categorical(df["Prediction"], categories=order, ordered=True)
    # df["Truth"] = pd.Categorical(df["Truth"], categories=order, ordered=True)
    df['Correct'] = df['Truth'] == df['Prediction']

    colors = list(plt.cm.PiYG(np.linspace(0.05, 1, 2)))
    g = sns.catplot(data=df, x="Correct",hue="Correct", y="Confidence", kind="violin", inner=None,orient='v',palette=colors,legend=False)
    sns.swarmplot(data=df, x="Correct", y="Confidence", color="black", size=3, ax=g.ax)

    plt.savefig("output/figures/violin.png")
    print('Violin plot generated successfully.')

def plot_cm(y_test,y_pred,out_name):
    # Computing the confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    ab = cm[0][0] + cm[1][1] + cm[0][1] + cm[0][2] + cm[1][2]
    c = cm[1][0] + cm[2][0] + cm[2][1] + cm[2][2]
    cm_compact = np.array(
        [
            [(cm[0][0] + cm[1][1]), (cm[0][1] + cm[0][2] + cm[1][2])],
            [(cm[1][0] + cm[2][0] + cm[2][1]), (cm[2][2])],
        ]
    )

    cm_normalized = cm_compact.astype('float') / cm_compact.sum(axis=1)[:, np.newaxis]*100

    # Plotting the confusion matrix using Matplotlib
    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(cm_normalized, interpolation="nearest", cmap=plt.cm.PuRd)
    ax.figure.colorbar(im, ax=ax)
    ax.set(
        xticks=np.arange(cm_normalized.shape[1]),
        yticks=np.arange(cm_normalized.shape[0]),
        xticklabels=["AB", "C"],
        yticklabels=["AB", "C"],
        title="Confusion Matrix",
        ylabel="True label",
        xlabel="Predicted label",
    )
    plt.setp(ax.get_xticklabels(), ha="right")
    thresh = cm_normalized.max() / 2.0
    for i in range(cm_normalized.shape[0]):
        for j in range(cm_normalized.shape[1]):
            ax.text(
                j,
                i,
                f'{cm_normalized[i, j]:.2f}%\n\n({cm_compact[i, j]})',
                ha="center",
                va="center",
                color="white" if cm_normalized[i, j] > thresh else "black",
            )
    fig.tight_layout()
    plt.savefig(f"output/figures/{out_name}.png")
    print('Confusion matrix of the model generated successfully.')

def plot_feature_importance(estimator,n_features,n_trials,X_train,y_train,X_test,y_test,out_name):

    # Make n_trials to obtain an average of the feature weights
    feature_weights = []
    feat=np.arange(n_features//2)
    columns=[]
    for i in feat:
        columns.append(f'Pred {i//2+1}')
        columns.append(f'Conf {i//2+1}')

    columns[-2:] = ['RMSD 0', 'RMSD 1']     

    ac = []
    df = pd.DataFrame(columns=columns)
    for i in range(n_trials):
        estimator = estimator.fit(X_train, y_train)
        feature_weights.append(estimator.feature_importances_)
        std = np.std(estimator.feature_importances_, axis=0)
        # df_aux = pd.DataFrame(estimator.feature_importances_.T, columns=columns)
        df.loc[len(df)] = estimator.feature_importances_
        ac.append(estimator.score(X_test,y_test))


    colors = list(plt.cm.PiYG(np.linspace(0.05, 1, 2)))
    plt.figure()
    feat = [np.mean(df[col]) for col in df]
    plt.bar(columns,feat,color=colors[0])
    plt.xticks(rotation=90, ha='center')
    plt.title(f'Mean accuracy: {np.mean(ac):.3f}')
    plt.tight_layout()
    plt.savefig(f'output/figures/{out_name}.png')
    print('Feature importance plot generated successfully.')


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--training_data",
        type=str,
        help="Path to the file containing the training data. Must be csv",
        required=True,
    )
    parser.add_argument(
        "--samples",
        type=int,
        help="The number of samples computed by the docking program.",
        required=True,
    )
    parser.add_argument(
        "--y_pred",
        type=int,
        help="Values predicted by the model.",
        required=True,
    )
    parser.add_argument(
        "--y_test",
        type=int,
        help="Real values of the labels.",
        required=True,
    )
    args = parser.parse_args()
    violin_confidence(args.training_data,args.samples,'confusion_matrix')
    plot_cm(args.y_test,args.y_pred,'features')