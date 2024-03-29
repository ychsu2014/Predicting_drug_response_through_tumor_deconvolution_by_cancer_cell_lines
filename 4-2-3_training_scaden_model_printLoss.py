# Code below are modified from the souce codes of the Scaden paper: https://www.science.org/doi/10.1126/sciadv.aba2619 
import logging
from scaden.model.architectures import architectures
#from scaden.model.scaden import Scaden

logger = logging.getLogger(__name__)

import os
import logging
import sys
import gc
import tensorflow as tf
import numpy as np
import pandas as pd
from anndata import read_h5ad
import collections
#from .functions import sample_scaling
from rich.progress import Progress, BarColumn
from sklearn import preprocessing as pp

in_data = sys.argv[1]
model_folder = sys.argv[2]
out_file = sys.argv[3]

# for testing
#in_data = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/processed_simulation_data/training_Bile_Duct_Cancer_processed.h5ad"
#model_folder = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/processed_simulation_data/Bile_Duct_Cancer/model/"
#out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/processed_simulation_data/Bile_Duct_Cancer/model/Bile_Duct_Cancer_loss_record.txt"

logger = logging.getLogger(__name__)

def sample_scaling(x, scaling_option):
    """
    Apply scaling of data
    :param x:
    :param scaling_option:
    :return:
    """

    if scaling_option == "log_min_max":
        # Bring in log space
        x = np.log2(x + 1)

        # Normalize data
        mms = pp.MinMaxScaler(feature_range=(0, 1), copy=True)

        # it scales features so transpose is needed
        x = mms.fit_transform(x.T).T

    return x

class Scaden(object):
    """
    scaden class
    """

    def __init__(
        self,
        model_dir,
        model_name,
        batch_size=128,
        learning_rate=0.0001,
        num_steps=1000,
        seed=0,
        hidden_units=[256, 128, 64, 32],
        do_rates=[0, 0, 0, 0],
    ):

        self.model_dir = model_dir
        self.batch_size = batch_size
        self.model_name = model_name
        self.beta1 = 0.9
        self.beta2 = 0.999
        self.learning_rate = learning_rate
        self.data = None
        self.n_classes = None
        self.labels = None
        self.x = None
        self.y = None
        self.num_steps = num_steps
        self.scaling = "log_min_max"
        self.sig_genes = None
        self.sample_names = None
        self.hidden_units = hidden_units
        self.do_rates = do_rates

        # Set seeds for reproducibility
        tf.random.set_seed(seed)
        os.environ["TF_DETERMINISTIC_OPS"] = "1"
        np.random.seed(seed)

    def scaden_model(self, n_classes):
        """Create the Scaden model"""

        model = tf.keras.Sequential()
        model.add(tf.keras.layers.Dense(self.hidden_units[0], activation=tf.nn.relu))
        model.add(tf.keras.layers.Dropout(self.do_rates[0]))
        model.add(tf.keras.layers.Dense(self.hidden_units[1], activation=tf.nn.relu))
        model.add(tf.keras.layers.Dropout(self.do_rates[1]))
        model.add(tf.keras.layers.Dense(self.hidden_units[2], activation=tf.nn.relu))
        model.add(tf.keras.layers.Dropout(self.do_rates[2]))
        model.add(tf.keras.layers.Dense(self.hidden_units[3], activation=tf.nn.relu))
        model.add(tf.keras.layers.Dropout(self.do_rates[3]))
        model.add(tf.keras.layers.Dense(n_classes, activation=tf.nn.softmax))

        return model

    def compute_loss(self, logits, targets):
        """
        Compute L1 loss
        :param logits:
        :param targets:
        :return: L1 loss
        """
        loss = tf.reduce_mean(input_tensor=tf.math.square(logits - targets))
        return loss

    def compute_accuracy(self, logits, targets, pct_cut=0.05):
        """
        Compute prediction accuracy
        :param targets:
        :param pct_cut:
        :return:
        """
        equality = tf.less_equal(
            tf.math.abs(tf.math.subtract(logits, targets)), pct_cut
        )
        accuracy = tf.reduce_mean(input_tensor=tf.cast(equality, tf.float32))
        return accuracy

    def correlation_coefficient(self, logits, targets):
        """
        Calculate the pearson correlation coefficient
        :param logits:
        :param targets:
        :return:
        """
        mx = tf.reduce_mean(input_tensor=logits)
        my = tf.reduce_mean(input_tensor=targets)
        xm, ym = logits - mx, targets - my
        r_num = tf.reduce_sum(input_tensor=tf.multiply(xm, ym))
        r_den = tf.sqrt(
            tf.multiply(
                tf.reduce_sum(input_tensor=tf.square(xm)),
                tf.reduce_sum(input_tensor=tf.square(ym)),
            )
        )
        r = tf.divide(r_num, r_den)
        r = tf.maximum(tf.minimum(r, 1.0), -1.0)
        return r

    def visualization(self, logits, targets, classes):
        """
        Create evaluation metrics
        :param targets:
        :param classes:
        :return:
        """
        # add evaluation metrics
        rmse = tf.compat.v1.metrics.root_mean_squared_error(logits, targets)[1]
        pcor = self.correlation_coefficient(logits, targets)
        eval_metrics = {"rmse": rmse, "pcor": pcor}

        for i in range(logits.shape[1]):
            eval_metrics[
                "mre_" + str(classes[i])
            ] = tf.compat.v1.metrics.mean_relative_error(
                targets[:, i], logits[:, i], targets[:, i]
            )[
                0
            ]
            eval_metrics[
                "mae_" + str(classes[i])
            ] = tf.compat.v1.metrics.mean_absolute_error(
                targets[:, i], logits[:, i], targets[:, i]
            )[
                0
            ]
            eval_metrics["pcor_" + str(classes[i])] = self.correlation_coefficient(
                targets[:, i], logits[:, i]
            )

        eval_metrics["mre_total"] = tf.compat.v1.metrics.mean_relative_error(
            targets, logits, targets
        )[1]

        eval_metrics["mae_total"] = tf.compat.v1.metrics.mean_relative_error(
            targets, logits, targets
        )[1]

        eval_metrics["accuracy01"] = self.compute_accuracy(
            logits, targets, pct_cut=0.01
        )
        eval_metrics["accuracy05"] = self.compute_accuracy(
            logits, targets, pct_cut=0.05
        )
        eval_metrics["accuracy1"] = self.compute_accuracy(logits, targets, pct_cut=0.1)

        # Create summary scalars
        for key, value in eval_metrics.items():
            tf.compat.v1.summary.scalar(key, value)

        tf.compat.v1.summary.scalar("loss", self.loss)

        merged_summary_op = tf.compat.v1.summary.merge_all()

        return merged_summary_op

    def load_h5ad_file(self, input_path, batch_size, datasets):
        """
        Load input data from a h5ad file and divide into training and test set
        :param input_path: path to h5ad file
        :param batch_size: batch size to use for training
        :param datasets: a list of datasets to extract from the file
        :return: Dataset object
        """
        try:
            print(input_path)
            raw_input = read_h5ad(input_path)
        except:
            logger.error(
                "Could not load training data file! Is it a .h5ad file generated with `scaden process`?"
            )
            sys.exit()

        # Subset dataset if --train_datasets is given
        if len(datasets) > 0:
            all_ds = collections.Counter(raw_input.obs["ds"])

            # Check that given datasets are all actually available
            for ds in datasets:
                if not ds in all_ds:
                    logger.warn(
                        f"The dataset '[cyan]{ds}[/cyan]' could not be found in the training data! Is the name correct?"
                    )

            for ds in all_ds:
                if ds not in datasets:
                    raw_input = raw_input[raw_input.obs["ds"] != ds].copy()

        # Create training dataset
        ratios = [raw_input.obs[ctype] for ctype in raw_input.uns["cell_types"]]
        self.x_data = raw_input.X.astype(np.float32)
        self.y_data = np.array(ratios, dtype=np.float32).transpose()
        self.data = tf.data.Dataset.from_tensor_slices((self.x_data, self.y_data))
        self.data = self.data.shuffle(1000).repeat().batch(batch_size=batch_size)
        self.data_iter = iter(self.data)

        # Extract celltype and feature info
        self.labels = raw_input.uns["cell_types"]
        self.sig_genes = list(raw_input.var_names)

    def load_prediction_file(self, input_path, sig_genes, labels, scaling=None):
        """
        Load a file to perform prediction on it
        :param input_path: path to input file
        :param sig_genes: the signature genes to use
        :param scaling: which scaling to perform
        :return: Dataset object
        """
        # Load data
        data = pd.read_table(input_path, sep="\t", index_col=0)
        sample_names = list(data.columns)

        # check for duplicates
        data_index = list(data.index)
        if not (len(data_index) == len(set(data_index))):
            logger.warn(
                "Scaden Warning: Your mixture file conatins duplicate genes! The first occuring gene will be used for every duplicate."
            )
            data = data.loc[~data.index.duplicated(keep="first")]

        data = data.loc[sig_genes]
        data = data.T

        # Scaling
        if scaling:
            data = sample_scaling(data, scaling_option=scaling)

        self.data = data

        return sample_names

    def build_model(self, input_path, train_datasets, mode="train"):
        """
        Build the model graph
        :param reuse:
        :return:
        """
        self.global_step = tf.Variable(0, name="global_step", trainable=False)

        # Load training data
        if mode == "train":
            self.load_h5ad_file(
                input_path=input_path,
                batch_size=self.batch_size,
                datasets=train_datasets,
            )

        # Load prediction data
        if mode == "predict":
            self.sample_names = self.load_prediction_file(
                input_path=input_path,
                sig_genes=self.sig_genes,
                labels=self.labels,
                scaling=self.scaling,
            )

        # Build the model or load if available
        self.n_classes = len(self.labels)

        try:
            self.model = tf.keras.models.load_model(self.model_dir, compile=False)
            logger.info(f"Loaded pre-trained model: [cyan]{self.model_name}")
        except:
            self.model = self.scaden_model(n_classes=self.n_classes)

    def train(self, input_path, train_datasets):
        """
        Train the model
        :param num_steps:
        :return:
        """
        ##### get the loss for the model training        
        loss_record_dict = {}

        # Define the optimizer
        optimizer = tf.keras.optimizers.Adam(learning_rate=self.learning_rate)

        # Build model graph
        self.build_model(
            input_path=input_path, train_datasets=train_datasets, mode="train"
        )

        # Training loop
        progress_bar = Progress(
            "[bold blue]{task.description}",
            "[bold cyan]Step: {task.fields[step]}, Loss: {task.fields[loss]}",
            BarColumn(bar_width=None),
        )

        training_progress = progress_bar.add_task(
            self.model_name, total=self.num_steps, step=0, loss=1
        )
        count = 0
        with progress_bar:

            for step in range(self.num_steps):

                x, y = self.data_iter.get_next()

                with tf.GradientTape() as tape:
                    self.logits = self.model(x, training=True)
                    loss = self.compute_loss(self.logits, y)
                    #####
                    count += 1
                    loss_record_dict[(count, step)] = loss
                grads = tape.gradient(loss, self.model.trainable_weights)

                optimizer.apply_gradients(zip(grads, self.model.trainable_weights))

                progress_bar.update(
                    training_progress, advance=1, step=step, loss=f"{loss:.4f}"
                )

                # Collect garbage after 100 steps - otherwise runs out of memory
                if step % 100 == 0:
                    gc.collect()

        # Save the trained model
        self.model.save(self.model_dir)
        pd.DataFrame(self.labels).to_csv(
            os.path.join(self.model_dir, "celltypes.txt"), sep="\t"
        )
        pd.DataFrame(self.sig_genes).to_csv(
            os.path.join(self.model_dir, "genes.txt"), sep="\t"
        )
        #####
        return(loss_record_dict)

    def predict(self, input_path):
        """
        Perform prediction with a pre-trained model
        :param input_path: prediction data path
        :return:
        """
        # Load signature genes and celltype labels
        sig_genes = pd.read_table(self.model_dir + "/genes.txt", index_col=0)
        self.sig_genes = list(sig_genes["0"])
        labels = pd.read_table(self.model_dir + "/celltypes.txt", index_col=0)
        self.labels = list(labels["0"])

        # Build model graph
        self.build_model(input_path=input_path, train_datasets=[], mode="predict")

        predictions = self.model.predict(self.data)

        pred_df = pd.DataFrame(
            predictions, columns=self.labels, index=self.sample_names
        )
        return pred_df


#PARAMETERS

# ==========================================#

# Extract architectures
M256_HIDDEN_UNITS = architectures["m256"][0]
M512_HIDDEN_UNITS = architectures["m512"][0]
M1024_HIDDEN_UNITS = architectures["m1024"][0]
M256_DO_RATES = architectures["m256"][1]
M512_DO_RATES = architectures["m512"][1]
M1024_DO_RATES = architectures["m1024"][1]

# ==========================================#

def training(
    data_path, train_datasets, model_dir, batch_size = 128, learning_rate = 0.0001, num_steps = 5000, seed=0
):
    # Convert training datasets
    if train_datasets == "":
        train_datasets = []
    else:
        train_datasets = train_datasets.split(",")
        logger.info(f"Training on: [cyan]{train_datasets}")

    # Training of M256 model
    logger.info("[cyan]Training M256 Model ... [/]")
    cdn256 = Scaden(
        model_dir=model_dir + "/m256",
        model_name="m256",
        batch_size=batch_size,
        learning_rate=learning_rate,
        num_steps=num_steps,
        seed=seed,
        hidden_units=M256_HIDDEN_UNITS,
        do_rates=M256_DO_RATES,
    )
    #####
    train_loss_256 = cdn256.train(input_path=data_path, train_datasets=train_datasets)
    del cdn256

    # Training of M512 model
    logger.info("[cyan]Training M512 Model ... [/]")
    cdn512 = Scaden(
        model_dir=model_dir + "/m512",
        model_name="m512",
        batch_size=batch_size,
        learning_rate=learning_rate,
        num_steps=num_steps,
        seed=seed,
        hidden_units=M512_HIDDEN_UNITS,
        do_rates=M512_DO_RATES,
    )
    #####
    train_loss_512 = cdn512.train(input_path=data_path, train_datasets=train_datasets)
    del cdn512

    # Training of M1024 model
    logger.info("[cyan]Training M1024 Model ... [/]")
    cdn1024 = Scaden(
        model_dir=model_dir + "/m1024",
        model_name="m1024",
        batch_size=batch_size,
        learning_rate=learning_rate,
        num_steps=num_steps,
        seed=seed,
        hidden_units=M1024_HIDDEN_UNITS,
        do_rates=M1024_DO_RATES,
    )
    #####
    train_loss_1024 = cdn1024.train(input_path=data_path, train_datasets=train_datasets)
    del cdn1024

    logger.info("[green]Training finished.")
    return([train_loss_256, train_loss_512, train_loss_1024])

[loss_256, loss_512, loss_1024] = training(in_data, "", model_folder)
fout = open(out_file, "w")
def writeLossToFile(inFile, inDict, inModelName):
    for akey in inDict.keys():
        (count, step) = akey
        loss = inDict[akey].numpy()
        fout.write("\t".join([inModelName, str(count), str(step), str(loss)]) + "\n")
writeLossToFile(fout, loss_256, "256")
writeLossToFile(fout, loss_512, "512")
writeLossToFile(fout, loss_1024, "1024")
fout.close()
