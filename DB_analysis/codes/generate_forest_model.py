import numpy as np
import pandas as pd
from sklearn.tree import DecisionTreeClassifier,DecisionTreeRegressor
from sklearn.model_selection import train_test_split
from argparse import ArgumentParser
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, GridSearchCV
import model_plots as mp
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from joblib import dump, load
import datetime
import time
import pickle


def load_features(data):
    # Load training data

    train = pd.read_csv(f"{data}/training_data.csv")
    test = pd.read_csv(f"{data}/test_data.csv")

    X_train = train.iloc[:, :-1]
    y_train = train["y_train"]

    X_test = test.iloc[:, :-1]
    y_test = test["y_test"]

    return X_train, y_train, X_test, y_test


def scale(X_train, X_test):
    # Scale the data
    feature_scaler = StandardScaler()
    X_train = feature_scaler.fit_transform(X_train)
    X_test = feature_scaler.transform(X_test)
    return X_train, X_test


def tree(data,model):

    red = "\033[31m"
    green = "\033[92m"
    end = "\033[0m"

    X_train, y_train, X_test, y_test = load_features(data)
    X_train, X_test = scale(X_train, X_test)

    # Make a first tree with default values
    clf = DecisionTreeClassifier(max_depth=3, splitter="best").fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    # Evaluating the model
    precision_ini = metrics.precision_score(
        y_test, y_pred, average="macro", labels=np.unique(y_pred)
    )
    recall_ini = metrics.recall_score(
        y_test, y_pred, average="macro", labels=np.unique(y_pred)
    )
    all_accuracies = cross_val_score(estimator=clf, X=X_train, y=y_train, cv=cv)

    if len(set(y_test) - set(y_pred)) > 0:
        print("")
        print(f"************* {green}TREE{end} *************")
        print(
            f"This labels were {red}never predicted{end}: {set(y_test)-set(y_pred)}, may be ill-defined"
        )
        print("********************************")
        print("")

    print("")
    print(f"~~~~~~~ {green}Tree metrics{end} ~~~~~~~")
    print(f" Mean Accuracy: {all_accuracies.mean()} +- {all_accuracies.std()}")
    print(
        f" Initial Precision: {precision_ini:.3f}"
    )  # (TA/(TA+FA)+TB/(TB+FB)+TC/(TC+FC))/3
    print(
        f" Initial Recall: {recall_ini:.3f}"
    )  # (TA/(TA+FB+FC)+TB/(TB+FA+FC)+TC/(TC+FB+FA))/3
    print(f" F1: {2/(1/precision_ini+1/recall_ini):.3f}")
    return clf


def grid_tree(
    cv,
    clf,
    data,
    max_depth=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    splitter=["random", "best"],
):
    red = "\033[31m"
    green = "\033[92m"
    end = "\033[0m"

    X_train, y_train, X_test, y_test = load_features(data)
    X_train, X_test = scale(X_train, X_test)

    grid_param = {
        "max_depth": max_depth,
        "splitter": splitter,
    }

    # Retrain to check different parameter values for the tree estimator
    gd_sr = GridSearchCV(
        estimator=clf, param_grid=grid_param, scoring="accuracy", cv=cv
    ).fit(X_train, y_train)

    best_parameters = gd_sr.best_params_
    best_estimator = gd_sr.best_estimator_

    # Predict the values using the best estimator
    y_pred = best_estimator.predict(X_test)

    print("")
    print(f" Best parameters: {best_parameters}")
    print(f" Best estimator: {best_estimator}")

    accuracy_fin = gd_sr.best_score_
    precision_fin = metrics.precision_score(
        y_test, y_pred, average="macro", labels=np.unique(y_pred)
    )
    recall_fin = metrics.recall_score(
        y_test, y_pred, average="macro", labels=np.unique(y_pred)
    )

    if len(set(y_test) - set(y_pred)) > 0:
        print("")
        print(f"********** {green}GRID TREE{end} **********")
        print(
            f"This labels were {red}never predicted{end}: {set(y_test)-set(y_pred)}, may be ill-defined"
        )
        print("************************************")
        print("")

    print("")
    print(f" Accuracy: {accuracy_fin}")
    print(
        f" Final Precision: {precision_fin:.3f}"
    )  # (TA/(TA+FA)+TB/(TB+FB)+TC/(TC+FC))/3
    print(
        f" Final Recall: {recall_fin:.3f}"
    )  # (TA/(TA+FB+FC)+TB/(TB+FA+FC)+TC/(TC+FB+FA))/3
    print(f" F1: {2/(1/precision_fin+1/recall_fin):.3f}")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("")

    # Make the plots of the model
    mp.plot_cm(y_test=y_test, y_pred=y_pred, out_name="confusion_matrix_tree")

    mp.plot_feature_importance(
        best_estimator,
        n_features=X_train.shape[1],
        n_trials=100,
        X_train=X_train,
        y_train=y_train,
        X_test=X_test,
        y_test=y_test,
        out_name="features_tree",
    )
    return best_estimator


# FOREST
def grid_forest(
    data,
    cv,
    model_name,
    max_depth=[5, 6, 7, 8, 9, 10, 11, 12],
    criterion=["gini", "log_loss", "entropy"],
    n_estimators=[200, 300, 400, 500],
    bootstrap=[True, False],
    min_samples_split=[2, 3, 4, 5],
    warm_start=[True, False],
):

    red = "\033[31m"
    green = "\033[42m"
    end = "\033[0m"

    X_train, y_train, X_test, y_test = load_features(data)
    X_train, X_test = scale(X_train, X_test)

    grid_param = {
        "max_depth": max_depth,
        "criterion": criterion,
        "n_estimators": n_estimators,
        "bootstrap": bootstrap,
        "min_samples_split": min_samples_split,
        "warm_start": warm_start,
    }

    est = False
    try:
        with open(".grid_params.pkl", "rb") as f:
            loaded_dict = pickle.load(f)

        if loaded_dict == grid_param:
            with open(".et.txt", "r") as file:
                t = file.readlines()
            est = True
        else:
            with open(".grid_params.pkl", "wb") as f:
                pickle.dump(grid_param, f)

    except:
        with open(".grid_params.pkl", "wb") as f:
            pickle.dump(grid_param, f)

    t1 = time.time()
    if est:
        print(f"----> Estimated time: {t} (Current time:{datetime.datetime.now()})")
    else:
        print(f"----> Current time:{datetime.datetime.now()}")

    clf = RandomForestClassifier(n_jobs=-1)
    gd_sr = GridSearchCV(
        estimator=clf,
        param_grid=grid_param,
        scoring="accuracy",
        cv=cv,
        verbose=1,
        error_score="raise",
        refit=True,
    )

    gd_sr.fit(X_train, y_train)
    t2 = time.time()

    print(f"----> Elapsed time: {str(datetime.timedelta(seconds=round(t2-t1)))}")

    with open(".et.txt", "w") as f:
        f.write(str(datetime.timedelta(seconds=round(t2 - t1))))

    # Collect the best model acording to gd_sr
    best_parameters = gd_sr.best_params_
    best_estimator = gd_sr.best_estimator_

    # Warn the user if the best parameters are the last ones tried (suggest to extend the grid search)
    for i in best_parameters:
        if best_parameters[i] == grid_param[i][-1]:
            try:
                if best_parameters[i] not in [True, False]:
                    float(best_parameters[i])
                    print(
                        f"{red}WARNING{end}: Best {i} is {best_parameters[i]} which is the {red}LAST{end} parameter tried!!!"
                    )
            except ValueError:
                pass

        if best_parameters[i] == grid_param[i][0]:
            try:
                if best_parameters[i] not in [True, False]:
                    float(best_parameters[i])
                    print(
                        f"{red}WARNING{end}: Best {i} is {best_parameters[i]} which is the {red}FIRST{end} parameter tried!!!"
                    )
            except ValueError:
                pass

    # test the best estimator
    clf = best_estimator
    clf.fit(X_train, y_train)
    ac = clf.score(X_test, y_test)
    y_pred = clf.predict(X_test)

    precision_fin = metrics.precision_score(
        y_test, y_pred, average="macro", labels=np.unique(y_pred)
    )
    recall_fin = metrics.recall_score(
        y_test, y_pred, average="macro", labels=np.unique(y_pred)
    )

    print("")
    print(f"~~~~~~~ {green}Forest metrics{end} ~~~~~~~")
    print(f" Using {best_parameters['n_estimators']} trees")
    print(f" Criterion: {best_parameters['criterion']}")
    print(f" Max Depth: {best_parameters['max_depth']}")
    print(f" Bootstrap: {best_parameters['bootstrap']}")
    print(f" Accuracy: {ac:.3f}")
    print(f" Precision: {precision_fin:.3f}")  # (TA/(TA+FA)+TB/(TB+FB)+TC/(TC+FC))/3
    print(f" Recall: {recall_fin:.3f}")  # (TA/(TA+FB+FC)+TB/(TB+FA+FC)+TC/(TC+FB+FA))/3
    print(f" F1: {2/(1/precision_fin+1/recall_fin):.3f}")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    # Save pickle of the model
    dump(
        clf, f"model/{model_name}"
    )  # To load the model: clf = load('model/model.joblib')

    # Make the plots of the forest
    mp.plot_cm(y_test=y_test, y_pred=y_pred, out_name="confusion_matrix_forest")

    mp.plot_feature_importance(
        best_estimator,
        n_features=X_train.shape[1],
        n_trials=100,
        X_train=X_train,
        y_train=y_train,
        X_test=X_test,
        y_test=y_test,
        out_name="features_forest",
    )
    return best_estimator


def regression_tree(data):
    red = "\033[31m"
    green = "\033[92m"
    end = "\033[0m"

    X_train, y_train, X_test, y_test = load_features(data)
    X_train, X_test = scale(X_train, X_test)

    # Make a first tree with default values
    clf = DecisionTreeRegressor(max_depth=3, splitter="best").fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    # Evaluating the model
    
    all_accuracies = cross_val_score(estimator=clf, X=X_train, y=y_train, cv=cv)


    print("")
    print(f"~~~~~~~ {green}Tree metrics{end} ~~~~~~~")
    print(f" Mean Error: {all_accuracies.mean()} +- {all_accuracies.std()}")
    return clf


def grid_regression_tree(
    cv,
    clf,
    data,
    max_depth=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    splitter=["random", "best"],
):
    red = "\033[31m"
    green = "\033[92m"
    end = "\033[0m"

    X_train, y_train, X_test, y_test = load_features(data)
    X_train, X_test = scale(X_train, X_test)

    grid_param = {
        "max_depth": max_depth,
        "splitter": splitter,
    }

    # Retrain to check different parameter values for the tree estimator
    gd_sr = GridSearchCV(
        estimator=clf, param_grid=grid_param, scoring="max_error", cv=cv
    ).fit(X_train, y_train)

    best_parameters = gd_sr.best_params_
    best_estimator = gd_sr.best_estimator_

    # Predict the values using the best estimator
    y_pred = best_estimator.predict(X_test)

    print("")
    print(f" Best parameters: {best_parameters}")
    print(f" Best estimator: {best_estimator}")

    accuracy_fin = gd_sr.best_score_

    print("")
    print(f" Max Error: {accuracy_fin}")



    # Make the plots of the model
    # mp.plot_cm(y_test=y_test, y_pred=y_pred, out_name="confusion_matrix_tree")

    mp.plot_feature_importance(
        best_estimator,
        n_features=X_train.shape[1],
        n_trials=100,
        X_train=X_train,
        y_train=y_train,
        X_test=X_test,
        y_test=y_test,
        out_name="features_tree",
    )
    return best_estimator


def grid_regression_forest(
    data,
    cv,
    model_name,
    max_depth=[5, 6, 7, 8, 9, 10, 11, 12],
    criterion=["squared_error", "friedman_mse", "absolute_error", "poisson"],
    n_estimators=[200, 300, 400, 500],
    bootstrap=[True, False],
    min_samples_split=[2, 3, 4, 5],
    warm_start=[True, False]
):
    red = "\033[31m"
    green = "\033[42m"
    end = "\033[0m"

    X_train, y_train, X_test, y_test = load_features(data)
    X_train, X_test = scale(X_train, X_test)

    grid_param = {
        "max_depth": max_depth,
        "criterion": criterion,
        "n_estimators": n_estimators,
        "bootstrap": bootstrap,
        "min_samples_split": min_samples_split,
        "warm_start": warm_start,
    }

    est = False
    try:
        with open(".grid_params.pkl", "rb") as f:
            loaded_dict = pickle.load(f)

        if loaded_dict == grid_param:
            with open(".et.txt", "r") as file:
                t = file.readlines()
            est = True
        else:
            with open(".grid_params.pkl", "wb") as f:
                pickle.dump(grid_param, f)

    except:
        with open(".grid_params.pkl", "wb") as f:
            pickle.dump(grid_param, f)

    t1 = time.time()
    if est:
        print(f"----> Estimated time: {t} (Current time:{datetime.datetime.now()})")
    else:
        print(f"----> Current time:{datetime.datetime.now()}")

    clf = RandomForestRegressor()
    gd_sr = GridSearchCV(
        estimator=clf,
        param_grid=grid_param,
        scoring='max_error',
        cv=cv,
        verbose=1,
        error_score="raise",
        refit=True,
    )

    gd_sr.fit(X_train, y_train)
    t2 = time.time()

    print(f"----> Elapsed time: {str(datetime.timedelta(seconds=round(t2-t1)))}")

    with open(".et.txt", "w") as f:
        f.write(str(datetime.timedelta(seconds=round(t2 - t1))))

        # Collect the best model acording to gd_sr
    best_parameters = gd_sr.best_params_
    best_estimator = gd_sr.best_estimator_

    # Warn the user if the best parameters are the last ones tried (suggest to extend the grid search)
    for i in best_parameters:
        if best_parameters[i] == grid_param[i][-1]:
            try:
                if best_parameters[i] not in [True, False]:
                    float(best_parameters[i])
                    print(
                        f"{red}WARNING{end}: Best {i} is {best_parameters[i]} which is the {red}LAST{end} parameter tried!!!"
                    )
            except ValueError:
                pass

        if best_parameters[i] == grid_param[i][0]:
            try:
                if best_parameters[i] not in [True, False]:
                    float(best_parameters[i])
                    print(
                        f"{red}WARNING{end}: Best {i} is {best_parameters[i]} which is the {red}FIRST{end} parameter tried!!!"
                    )
            except ValueError:
                pass

    # test the best estimator
    clf = best_estimator
    clf.fit(X_train, y_train)
    ac = clf.score(X_test, y_test)
    y_pred = clf.predict(X_test)

    precision_fin = metrics.precision_score(
        y_test, y_pred, average="macro", labels=np.unique(y_pred)
    )
    recall_fin = metrics.recall_score(
        y_test, y_pred, average="macro", labels=np.unique(y_pred)
    )

    print("")
    print(f"~~~~~~~ {green}Regression metrics{end} ~~~~~~~")
    print(f" Using {best_parameters['n_estimators']} trees")
    print(f" Criterion: {best_parameters['criterion']}")
    print(f" Max Depth: {best_parameters['max_depth']}")
    print(f" Bootstrap: {best_parameters['bootstrap']}")
    print(f" Accuracy: {ac:.3f}")
    print(f" Precision: {precision_fin:.3f}")  # (TA/(TA+FA)+TB/(TB+FB)+TC/(TC+FC))/3
    print(f" Recall: {recall_fin:.3f}")  # (TA/(TA+FB+FC)+TB/(TB+FA+FC)+TC/(TC+FB+FA))/3
    print(f" F1: {2/(1/precision_fin+1/recall_fin):.3f}")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    # Save pickle of the model
    dump(
        clf, f"model/{model_name}"
    )  # To load the model: clf = load('model/model.joblib')

    # Make the plots of the forest
    mp.plot_cm(y_test=y_test, y_pred=y_pred, out_name="confusion_matrix_regression")

    mp.plot_feature_importance(
        best_estimator,
        n_features=X_train.shape[1],
        n_trials=100,
        X_train=X_train,
        y_train=y_train,
        X_test=X_test,
        y_test=y_test,
        out_name="features_regression",
    )
    return best_estimator


if __name__ == "__main__":
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
    parser.add_argument(
        "--labels",
        type=str,
        help="Kind of labels: options are: 'ABC' or 'logk'. Default is 'logk'",
        default="logk",
        required=False,
    )
    args = parser.parse_args()

    cv = 5

    if args.labels == "ABC":
        clf_tree = tree(args.data)
        best_tree = grid_tree(cv, clf_tree, args.data)
        print("")
        best_forest = grid_forest(args.data, cv, args.model_name)
    elif args.labels == "logk":
        clf_tree = regression_tree(args.data)
        best_tree = grid_regression_tree(cv, clf_tree, args.data)
        print('')
        clf_reg = grid_regression_forest(args.data,cv,args.model_name)
    else:
        print(f"Error: {args.labels} is not a valid label. Use 'logk' or 'ABC'")

