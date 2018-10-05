
#################################################################################
# FUNCTIONS                                                                     #
#################################################################################


#################################################################################
# GLOBALS                                                                       #
#################################################################################


#################################################################################
# RULES                                                                         #
#################################################################################

rule all:
    input: "data/interim/gbt/tio_model_random_stage4.pkl",

rule tio:
    input: kmers="data/interim/kmers/kmer_matrix.npz",
        mics="data/interim/mic_class_dataframe2.pkl",
        classes="data/interim/mic_class_order_dict2.pkl"
    output: data="data/interim/gbt/tio_stage1.npz",
    run:
        import numpy as np
        import pandas as pd
        from sklearn.externals import joblib
        with np.load(input.kmers) as data:
            kmers = data['kmers']
            kmer_order = data['kmer_order']
            genome_order = data['genome_order']

        micsdf = joblib.load(input.mics)
        class_orders = joblib.load(input.classes)

        tio_labels = class_orders['TIO']
        tio_label_index = { k: v for v, k in enumerate(tio_labels) }
        y_tio = np.array([ tio_label_index[m] if not pd.isna(m) else m for m in micsdf.loc[genome_order, 'TIO'] ])
        labels, counts = np.unique(y_tio, return_counts=True)
        ok = labels[counts >= 5]

        mask = np.in1d(y_tio, ok) # Since Nan is not a label, this also filters invalid MICs
        y_tio = y_tio[mask]
        X_tio = kmers[mask,:]
        tio_samples = genome_order[mask]

        np.savez(output.data, y=y_tio, X=X_tio, samples=tio_samples, labels=tio_labels)


rule clonalfilter:
    input: data="data/interim/gbt/tio_stage1.npz",
        groups="data/interim/mash_population_groups.csv"
    output: data="data/interim/gbt/tio_filtered_clonal_groups_stage2.npz"
    params:
        min=5
    run:
        import numpy as np
        import pandas as pd

        clonaldf = pd.read_csv(input.groups, index_col=0)
        with np.load(input.data) as data:
            X = data['X']
            y = data['y']
            samples = data['samples']

            f = lambda x: np.unique(clonaldf.loc[samples[x > 0]]).size >= params.min
            mask = np.apply_along_axis(f, 0, X)

            print(mask.shape)
            print(mask)
            print(sum(mask))

            X1 = X[:,mask]

            np.savez(output.data, mask=mask, X=X1)

rule clonalsignificance:
    input: xdata="data/interim/gbt/tio_filtered_clonal_groups_stage2.npz",
        ydata="data/interim/gbt/tio_stage1.npz",
        groups="data/interim/mash_population_groups.csv"
    output:
        data="data/interim/gbt/tio_significant_clonal_groups_stage3.npz"
    params:
        min=5
    run:
        import numpy as np
        import pandas as pd
        from sklearn.feature_selection import chi2, f_classif

        clonaldf = pd.read_csv(input.groups, index_col=0)
        with np.load(input.xdata) as xdata, np.load(input.ydata) as ydata:
            X = xdata['X']
            y = ydata['y']
            samples = ydata['samples']

            clonaldf = clonaldf.loc[samples]

            n = X.shape[1]
            print(n)
            clusters = clonaldf['cluster'].value_counts()
            significant = np.zeros(n)
            pcutoff = 0.05/n
            print(pcutoff)
            for c in clusters.index:
                mask = clonaldf.cluster == c
                dx = X[mask,:]
                dy = y[mask]
                fvals, pvals = f_classif(dx,dy)

                sig = pvals < pcutoff
                print(sum(sig))

                significant[sig] += 1


            passed = significant > params.min

            print("final")
            print(sum(passed))

            np.savez(output.data, mask=passed, X=X[:,passed])

rule random_model:
    input: data="data/interim/gbt/tio_stage1.npz",
    output:
        report="data/interim/gbt/tio_performance_report_random_stage4.txt",
        model="data/interim/gbt/tio_model_random_stage4.pkl",
    params:
        features=270
    run:
        import numpy as np
        import pandas as pd
        from xgboost import XGBClassifier

        from sklearn import metrics
        from sklearn.externals import joblib
        from sklearn.model_selection import train_test_split, RandomizedSearchCV, GridSearchCV
        from sklearn.feature_selection import SelectKBest, f_classif
        from sklearn.metrics import matthews_corrcoef, classification_report, precision_recall_fscore_support, confusion_matrix

        from hpsklearn import HyperoptEstimator, xgboost_classification
        from hyperopt import tpe

        # Load
        with np.load(input.data) as data:
            X = data['X']
            y = data['y']

            # Split test / train
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=36, stratify=y)

            # Feature selection
            fsel = SelectKBest(f_classif, k=params.features)
            X_train_fs = fsel.fit_transform(X_train, y_train)
            X_test_fs = fsel.transform(X_test)

            # Train / find best params:
            model = HyperoptEstimator( classifier=xgboost_classification('xbc'), preprocessing=[], algo=tpe.suggest, trial_timeout=2000 )
            print(type(X_train_fs))
            print(type(y_train))
            model.fit( X_train_fs, y_train )
            best_model = model.best_model()['learner']

            # Compute performance
            # Accuracy
            prediction = best_model.predict(X_test_fs)
            acc = sum(prediction == y_test)/len(y_test)
            print(acc)

        with open(output.report, 'w') as outfh:
            tmpfh = open(output.model, 'a')
            print(output.report)
            outfh.write('Accuracy: {}\n'.format(acc))
            print(acc)
            tmpfh.write('HHHHHH {}'.format(type(best_model)))
