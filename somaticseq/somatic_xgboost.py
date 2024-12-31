#!/usr/bin/env python3

import argparse
import logging
from copy import copy
from typing import Any

import pandas as pd
import xgboost as xgb

import somaticseq.ntchange_type as ntchange
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


def param_list_to_dict(
    param_list: list[str], existing_param_dict: dict[str, Any] = DEFAULT_PARAM
) -> dict[str, Any]:
    """
    Args:
        param_list: this is what will be passed from the CLI, e.g.,
            ["scale_pos_weight:0.8", "seed:42", "grow_policy:lossguide"]. If the
            value is integer, float, bool, etc., it will be eval'ed as such.
            Otherwise, it'll remain as a string.
        existing_param_dict: a pre-existing set of params and their values
            before param_list adds and/or modifies it.

    Returns:
        An updated params dict
    """
    updated_param_dict = copy(existing_param_dict)
    for param_string in param_list:
        param_i, value_i = param_string.split(":")
        try:
            value_i = eval(value_i)
        except NameError:
            # value_i stays a string if cannot be eval'ed
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
    input_tsvs: list[str],
    param: dict[str, Any] = DEFAULT_PARAM,
    non_feature: list[str] = NON_FEATURE,
    num_rounds: int = DEFAULT_XGB_BOOST_ROUNDS,
    model: str | None = None,
) -> str:
    """
    Build SomaticSeq's somatic mutation classifiers

    Args:
        input_tsvs: a list of labeled training data sets
        params: hyperparameters for model training
        non_features: list of columns in the tsvs not to be used as training.
            Those typically include chromosome coordinates and the data label.
        num_rounds: number of boosting rounds
        model: the output classifier file name. If None, will use input file
            name as basename.
    Returns:
        The classifier file path
    """
    logger = logging.getLogger("xgboost_" + builder.__name__)
    logger.info("TRAINING {} for XGBOOST".format(",".join(input_tsvs)))
    logger.info("Columns removed before training: {}".format(", ".join(non_feature)))
    logger.info(f"Number of boosting rounds = {num_rounds}")
    logger.info("Hyperparameters: " + ", ".join([f"{i}={param[i]}" for i in param]))

    if model is None:
        model = input_tsvs[0] + f".xgb.v{__version__}.classifier"

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
    bst.dump_model(f"{model}.json", with_stats=True, dump_format="json")
    bst.dump_model(f"{model}.txt", with_stats=True, dump_format="text")
    save_feature_importance_to_file(bst, f"{model}.feature_importance.txt")

    return model


def predictor(
    model: str,
    input_tsv: str,
    output_tsv: str,
    non_feature: list[str] = NON_FEATURE,
    iterations: int = DEFAULT_NUM_TREES_PREDICT,
) -> str:
    """
    Uses an existing SomaticSeq classifier to predict somatic mutations

    Args:
        model: path to the existing SomaticSeq classifier
        input_tsv: the SomaticSeq tsv file for which to classify somatic
            mutation status
        output_tsv: adds predicted label to the input_tsv above
        non_feature: features to exclude from input_tsv
        iterations: number of trees to use from the model

    Returns:
        output_tsv file path
    """
    logger = logging.getLogger("xgboost_" + predictor.__name__)
    logger.info("Columns removed for prediction: {}".format(",".join(non_feature)))
    logger.info(f"Number of trees to use = {iterations}")

    xgb_model = xgb.Booster()
    xgb_model.load_model(model)

    # Read in chunks, start with write/Header, and append/NoHeader afterwards.
    chunksize = 10000
    write_or_append, write_header_or_not = "w", True

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
            mode=write_or_append,
            header=write_header_or_not,
            na_rep="nan",
        )
        write_or_append, write_header_or_not = "a", False

    return output_tsv


# Execute:
def main() -> None:
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
        help=(
            "extra xgboost training parameters in format of "
            "PARAM_1:VALUE_1 PARAM_2:VALUE_2. "
            "Will overwrite defaults and other options."
        ),
    )
    parser_train.add_argument(
        "--features-excluded",
        nargs="*",
        type=str,
        help=(
            "features to exclude for xgboost training. Must be same for train/predict."
        ),
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
        help=(
            "features to exclude for xgboost training. Must be same for train/predict."
        ),
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

        builder(
            args.tsvs_in,
            param=PARAM,
            non_feature=NON_FEATURE,
            num_rounds=args.num_boost_rounds,
            model=args.model_out,
        )
    elif args.which == "predict":
        for feature_i in args.features_excluded:
            NON_FEATURE.append(feature_i)

        predictor(
            args.model,
            args.tsv_in,
            args.predicted_tsv,
            non_feature=NON_FEATURE,
            iterations=args.num_trees,
        )


if __name__ == "__main__":
    main()
