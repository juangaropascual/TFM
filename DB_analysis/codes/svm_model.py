import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from argparse import ArgumentParser
from sklearn import metrics
from sklearn.preprocessing import StandardScaler,PolynomialFeatures,MinMaxScaler
from sklearn.model_selection import cross_val_predict,GridSearchCV
import model_plots as mp
from sklearn import svm
from joblib import dump, load
import datetime
import time
from sklearn.model_selection import ShuffleSplit
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

parser = ArgumentParser()
parser.add_argument(
    "--data",
    type=str,
    help="Path to the directory containing the test and training data.",
    required=True,
)
parser.add_argument(
    "--n_samples",
    type=int,
    help="The number of different samples from which the features will be extracted.",
    required=True,
)
parser.add_argument(
    "--model_name",
    type=str,
    help="Name of the model that will be saved in a .joblib file.",
    required=True,
)
args = parser.parse_args()

# Load training data

train = pd.read_csv(f'{args.data}/training_data.csv')
test = pd.read_csv(f'{args.data}/test_data.csv')

X_train = train.iloc[:, :-1]
y_train = train['y_train']

X_test = test.iloc[:, :-1]
y_test = test['y_test']

# Make data separable by using polynomial behaviour
# poly = PolynomialFeatures(3)
# X_train = poly.fit_transform(X_train)
# X_test = poly.transform(X_test)


# Scale the data 
feature_scaler = StandardScaler()
# feature_scaler = MinMaxScaler()
X_train = feature_scaler.fit_transform(X_train)
X_test = feature_scaler.transform(X_test)

clf = svm.SVC(C=5, kernel='poly', degree=3,coef0=1).fit(X_train, y_train)

y_pred = clf.predict(X_test)

precision_ini = metrics.precision_score(y_test, y_pred, average="macro",labels=np.unique(y_pred))
recall_ini = metrics.recall_score(y_test, y_pred, average="macro",labels=np.unique(y_pred))
score = clf.score(X_test,y_test)

print(f' Score: {score:.3f}')
print(f' Precision: {precision_ini:.3f}')
print(f' Recall: {recall_ini:.3f}')
print(f" F1: {2/(1/precision_ini+1/recall_ini):.3f}")

grid_param = {
    'C': [1,2,3,4,5, 6, 7, 8, 9, 10, 11, 12],
    'degree': [2,3,4,5],
    'coef0': [1,2,3,4,5],
    'max_iter': [1000,1500,2000,2500]
}

t1 = time.time()
print(f'----> Current time:{datetime.datetime.now()}')

clf = svm.SVC(kernel='poly')
gd_sr = GridSearchCV(estimator=clf,
                    param_grid=grid_param,
                    scoring='balanced_accuracy',
                    cv=5,
                    verbose=1).fit(X_train, y_train)
t2 = time.time()

# Collect the best model acording to gd_sr
best_parameters = gd_sr.best_params_
print(best_parameters)
best_estimator = gd_sr.best_estimator_

# Warn the user if the best parameters are the last ones tried (suggest to extend the grid search)
for i in best_parameters:
    if best_parameters[i] == grid_param[i][-1]:
        try:
            float(best_parameters[i])
            print(f'WARNING: Best {i} is {best_parameters[i]} which is the LAST parameter tried!!!')
        except ValueError:
            pass

    if best_parameters[i] == grid_param[i][0]:
        try:
            float(best_parameters[i])
            print(f'WARNING: Best {i} is {best_parameters[i]} which is the FIRST parameter tried!!!')
        except ValueError:
            pass
best_estimator.fit(X_train, y_train)
ac = best_estimator.score(X_test, y_test)
print(ac)




grid_param = {
    'C': [0.001,0.01,0.1,0.5,1,2,3,4,5],
    'gamma': [0.001,0.01,0.1,0.5,1,2,3,4,5,6,7,8],
}


t1 = time.time()
print(f'----> Current time:{datetime.datetime.now()}')

clf = svm.SVC(kernel='rbf')
gd_sr = GridSearchCV(estimator=clf,
                    param_grid=grid_param,
                    scoring='balanced_accuracy',
                    cv=5,
                    verbose=1).fit(X_train, y_train)
t2 = time.time()

# Collect the best model acording to gd_sr
best_parameters = gd_sr.best_params_
print(best_parameters)
best_estimator = gd_sr.best_estimator_

# Warn the user if the best parameters are the last ones tried (suggest to extend the grid search)
for i in best_parameters:
    if best_parameters[i] == grid_param[i][-1]:
        try:
            float(best_parameters[i])
            print(f'WARNING: Best {i} is {best_parameters[i]} which is the LAST parameter tried!!!')
        except ValueError:
            pass

    if best_parameters[i] == grid_param[i][0]:
        try:
            float(best_parameters[i])
            print(f'WARNING: Best {i} is {best_parameters[i]} which is the FIRST parameter tried!!!')
        except ValueError:
            pass
best_estimator.fit(X_train, y_train)
ac = best_estimator.score(X_test, y_test)
print(ac)


dump(best_estimator, f'model/{args.model_name}') # To load the model: clf = load('model/model.joblib') 
