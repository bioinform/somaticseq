#!/usr/bin/env python3

import argparse
import logging
from copy import copy

import pandas as pd
import somaticseq.ntchange_type as ntchange
import xgboost as xgb
from somaticseq._version import __version__

FORMAT = "%(levelname)s %(asctime)-15s %(name)-20s %(message)s"
logger = logging.getLogger("Somatic_Xgboost")
logger.setLevel(logging.DEBUG)
logging.basicConfig(level=logging.DEBUG, format=FORMAT)


DEFAULT_PARAM = {
    "max_depth": 8,
    "nthread": 1,
    "objective": "binary:logistic",
    "seed": 0,
    "tree_method": "hist",
    "grow_policy": "lossguide",
}
NON_FEATURE = [
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "Strelka_QSS",
    "Strelka_TQSS",
    "if_COSMIC",
    "COSMIC_CNT",
    "TrueVariant_or_False",
]

DEFAULT_XGB_BOOST_ROUNDS = 500
DEFAULT_NUM_TREES_PREDICT = 100


def param_list_to_dict(param_list, existing_param_dict=DEFAULT_PARAM):
    """
    param_list is what will be passed from the CLI, e.g.,
    ['scale_pos_weight:0.8', 'seed:42'] If the value is integer, float, bool,
    etc., it will be eval'ed as such. Otherwise, it'll remain as a string.
    """

    updated_param_dict = copy(existing_param_dict)
    for param_string in param_list:
        param_i, value_i = param_string.split(":")
        try:
            value_i = eval(value_i)
        except NameError:
            pass

        updated_param_dict[param_i] = value_i

    return updated_param_dict


def save_feature_importance_to_file(xgb_model, filename):
    feature_gain = xgb_model.get_score(importance_type="gain")
    feature_weight = xgb_model.get_score(importance_type="weight")
    feature_cover = xgb_model.get_score(importance_type="cover")
    feature_total_gain = xgb_model.get_score(importance_type="total_gain")
    feature_total_cover = xgb_model.get_score(importance_type="total_cover")

    line_i = "{}\t{}\t{}\t{}\t{}\t{}\n".format(
        "FEATURE", "GAIN", "WEIGHT", "COVER", "TOTAL_GAIN", "TOTAL_COVER"
    )

    with open(filename, "w") as fout:
        fout.write(line_i)

        for feature_i in sorted(feature_gain):
            line_i = "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                feature_i,
                feature_gain[feature_i],
                feature_weight[feature_i],
                feature_cover[feature_i],
                feature_total_gain[feature_i],
                feature_total_cover[feature_i],
            )
            fout.write(line_i)

    return True


def builder(
    input_tsvs,
    param=DEFAULT_PARAM,
    non_feature=NON_FEATURE,
    num_rounds=DEFAULT_XGB_BOOST_ROUNDS,
    model=None,
):
    logger = logging.getLogger("xgboost_" + builder.__name__)
    logger.info("TRAINING {} for XGBOOST".format(",".join(input_tsvs)))
    logger.info("Columns removed before training: {}".format(", ".join(non_feature)))
    logger.info("Number of boosting rounds = {}".format(num_rounds))
    logger.info(
        "Hyperparameters: " + ", ".join(["{}={}".format(i, param[i]) for i in param])
    )

    if not model:
        model = input_tsvs[0] + ".xgb.v{}.classifier".format(__version__)

    input_data = pd.concat(
        [
            pd.read_csv(input_tsv_i, sep="\t", low_memory=False)
            for input_tsv_i in input_tsvs
        ]
    )

    train_data = ntchange.ntchange(input_data)
    for non_feature_i in non_feature:
        if non_feature_i in train_data:
            train_data.drop(
                [
                    non_feature_i,
                ],
                axis=1,
                inplace=True,
            )

    train_label = input_data["TrueVariant_or_False"]

    dtrain = xgb.DMatrix(train_data, label=train_label)
    bst = xgb.train(param, dtrain, num_boost_round=num_rounds)

    bst.save_model(model)
    bst.dump_model(model + ".txt")
    save_feature_importance_to_file(bst, model + ".feature_importance.txt")

    return model


def predictor(
    model,
    input_tsv,
    output_tsv,
    non_feature=NON_FEATURE,
    iterations=DEFAULT_NUM_TREES_PREDICT,
):
    logger = logging.getLogger("xgboost_" + predictor.__name__)
    logger.info("Columns removed for prediction: {}".format(",".join(non_feature)))
    logger.info("Number of trees to use = {}".format(iterations))

    xgb_model = xgb.Booster()
    xgb_model.load_model(model)

    # Read in chunks, start with write/Header, and append/NoHeader afterwards.
    chunksize = 10000
    writeMode, writeHeader = "w", True

    for input_data in pd.read_csv(
        input_tsv, sep="\t", chunksize=chunksize, low_memory=False
    ):
        test_data = ntchange.ntchange(input_data)
        for non_feature_i in non_feature:
            if non_feature_i in test_data:
                test_data.drop(non_feature_i, axis=1, inplace=True)

        dtest = xgb.DMatrix(test_data)
        scores = xgb_model.predict(dtest, iteration_range=(0, iterations))
        predicted = input_data.assign(SCORE=scores)

        predicted.to_csv(
            output_tsv,
            sep="\t",
            index=False,
            mode=writeMode,
            header=writeHeader,
            na_rep="nan",
        )

        writeMode, writeHeader = "a", False

    return output_tsv


################################################################################################
# Execute:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run XGBoost",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    sample_parsers = parser.add_subparsers(title="mode")

    # TRAINING mode
    parser_train = sample_parsers.add_parser("train")
    parser_train.add_argument(
        "-tsvs",
        "--tsvs-in",
        type=str,
        nargs="+",
        help="labeled tsv file(s)",
        required=True,
    )
    parser_train.add_argument(
        "-out", "--model-out", type=str, help="output model file name"
    )
    parser_train.add_argument(
        "-threads", "--num-threads", type=int, help="num threads."
    )
    parser_train.add_argument(
        "-depth", "--max-depth", type=int, help="tree max depth. default=8"
    )
    parser_train.add_argument(
        "-seed", "--seed", type=int, help="random seed. default=0"
    )
    parser_train.add_argument(
        "-method", "--tree-method", type=str, help="tree method. default=hist"
    )
    parser_train.add_argument(
        "-iter",
        "--num-boost-rounds",
        type=int,
        help="num boosting rounds, i.e., number of trees",
        default=1000,
    )
    parser_train.add_argument(
        "--extra-params",
        nargs="*",
        type=str,
        help="extra xgboost training parameters in format of PARAM_1:VALUE_1 PARAM_2:VALUE_2. Will overwrite defaults and other options.",
    )
    parser_train.add_argument(
        "--features-excluded",
        nargs="*",
        type=str,
        help="features to exclude for xgboost training. Must be same for train/predict.",
        default=[],
    )
    parser_train.set_defaults(which="train")

    # PREDICTION mode
    parser_predict = sample_parsers.add_parser("predict")
    parser_predict.add_argument(
        "-model", "--model", type=str, help="xgboost model", required=True
    )
    parser_predict.add_argument(
        "-tsv", "--tsv-in", type=str, help="tsv file in", required=True
    )
    parser_predict.add_argument(
        "-out", "--predicted-tsv", type=str, help="tsv file out", required=True
    )
    parser_predict.add_argument(
        "-ntrees",
        "--num-trees",
        type=int,
        help="only use this many trees to classify",
        default=100,
    )
    parser_predict.add_argument(
        "--features-excluded",
        nargs="*",
        type=str,
        help="features to exclude for xgboost training. Must be same for train/predict.",
        default=[],
    )
    parser_predict.set_defaults(which="predict")

    args = parser.parse_args()

    if args.which == "train":
        PARAM = copy(DEFAULT_PARAM)

        if args.num_threads:
            PARAM["nthread"] = args.num_threads

        if args.max_depth:
            PARAM["max_depth"] = args.max_depth

        if args.tree_method:
            PARAM["seed"] = args.seed

        if args.tree_method:
            PARAM["tree_method"] = args.tree_method

        # If they're integers, floats, bools, etc., they will be eval'ed to them.
        # If they're strings, they're remain strings.
        if args.extra_params:
            PARAM = param_list_to_dict(args.extra_params, PARAM)

        for feature_i in args.features_excluded:
            NON_FEATURE.append(feature_i)

        _ = builder(
            args.tsvs_in,
            param=PARAM,
            non_feature=NON_FEATURE,
            num_rounds=args.num_boost_rounds,
            model=args.model_out,
        )

    elif args.which == "predict":
        for feature_i in args.features_excluded:
            NON_FEATURE.append(feature_i)

        _ = predictor(
            args.model,
            args.tsv_in,
            args.predicted_tsv,
            non_feature=NON_FEATURE,
            iterations=args.num_trees,
        )
