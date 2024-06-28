import numpy as np
import pandas as pd
from argparse import ArgumentParser
from sklearn.model_selection import ShuffleSplit
import os
import model_plots as mp


parser = ArgumentParser()
parser.add_argument(
    "--training_data",
    type=str,
    help="Path to the file containing the training data. Must be csv",
    required=True,
)
parser.add_argument(
    "--n_samples",
    type=int,
    help="The number of different samples from which the features will be extracted.",
    required=True,
)
parser.add_argument(
    "--labels",
    type=str,
    help="Kind of labels: options are: 'ABC' or 'logk'. Default is 'logk'",
    default='logk',
    required=False,
)
args = parser.parse_args()


# Load training data

df = pd.read_csv(args.training_data)
df_rmsd = pd.read_csv("output/RMSD_DF.csv")
id = df["ID"].loc[0]
samples = len(df[df["ID"] == id])

res_dir = args.training_data.split('.')[0]

if args.labels == 'ABC':
    mp.violin_confidence(args.training_data,samples)
else:
    print('Use --labels ABC to plot the violin of the results.')


try:
    os.mkdir(res_dir)
except:
    pass


# Add experimental prediction to test 
experimental = df["Experimental Pred"].iloc[::samples].values

# Normalize by withdrawing rank1 to every inferior rank
rank1 = df["Confidence"].iloc[::samples].values #rank1 values only
for i in range(len(rank1)):
    df.loc[i * samples : (i + 1) * samples, "Confidence"] -= rank1[i]

chain = df["Chain(A=0)(B=1)"].values
confidence = df["Confidence"].values
assert len(chain) == len(
    confidence
), "Length of chain and confidence arrays do not match."

# Add chain and confidence of each sample as features
interleaved_array = np.ravel(np.column_stack((chain, confidence)))
features = interleaved_array.reshape(len(np.unique(df["ID"])), samples * 2)

# Add RMSD as a feature 
rmsd_0 = df_rmsd["RMSD_0"].values
rmsd_1 = df_rmsd["RMSD_1"].values

# Select how many of the sample features should be taken into account (n_samples)
features = features[:, : args.n_samples * 2]
X = np.column_stack((features, rmsd_0, rmsd_1))
y = experimental


# Create training, test and validation sets.
ShuffleSplit(n_splits=1, test_size=0.2, random_state=1).get_n_splits(X,y)
train_index, val_index = next(ShuffleSplit(n_splits=1, test_size=0.2,random_state=1).split(X, y)) 
X_train, X_val = X[train_index], X[val_index] 
y_train, y_val = y[train_index], y[val_index]

ShuffleSplit(n_splits=1, test_size=0.25,random_state=1).get_n_splits(X_train,y_train)
train_index, test_index = next(ShuffleSplit(n_splits=1, test_size=0.25,random_state=1).split(X_train, y_train)) 
X_train, X_test = X[train_index], X[test_index] 
y_train, y_test = y[train_index], y[test_index]

# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)
# X_train, X_val, y_train, y_val = train_test_split(
#     X_train, y_train, test_size=0.25, random_state=1
# )


columns=[]
for i in range(int(len(X[0])/2)):
    columns.append(f'Pred {i+1}')
    columns.append(f'Conf {i+1}')

columns[-2:] = ['RMSD 0', 'RMSD 1']   

train_data = pd.DataFrame(columns=columns)
for i,col in enumerate(columns):
    train_data[col] = X_train[:,i]
train_data['y_train'] = y_train

train_data.to_csv(f'{res_dir}/training_data.csv', index=False)

test_data = pd.DataFrame(columns=columns)
for i,col in enumerate(columns):
    test_data[col] = X_test[:,i]
test_data['y_test'] = y_test

test_data.to_csv(f'{res_dir}/test_data.csv', index=False)

val_df = pd.DataFrame()
val_df['ID'] = df['ID'].iloc[val_index*samples] 
val_df['Experimental Pred'] = df['Experimental Pred'].iloc[val_index*samples]

for i,col in enumerate(columns):
    val_df[col] = X_val[:,i]

val_df.to_csv(f'{res_dir}/validation_data.csv', index=False)