# Retrosynthesis prediction enhanced by in-silico reaction data augmentation


[[arxiv](https://arxiv.org/abs/2402.00086)]

The directory contains source code of the unpublished article:

Zhang et al's Retrosynthesis prediction enhanced by in-silico reaction data augmentation.

## Data

USPTO-50K: [Dai's version of USPTO dataset](https://www.dropbox.com/sh/6ideflxcakrak10/AADTbFBC0F8ax55-z-EDgrIza)

USPTO-MIT: https://github.com/wengong-jin/nips17-rexgen/blob/master/USPTO/data.zip

USPTO-FULL: [Dai's version of USPTO dataset](https://www.dropbox.com/sh/6ideflxcakrak10/AADTbFBC0F8ax55-z-EDgrIza)

Extra data used for USPTO-50K:

https://figshare.com/articles/dataset/RetroWISE-data_zip/25272451

## Environment Preparation

Please make sure you have installed anaconda or miniconda. The version about `pytorch` and `cudatoolkit` should be depended on your machine. The version of `pytorch` should not be smaller than 1.6 according to the OpenNMT-py.

```shell
conda create -n retrowise python=3.7 \
conda activate retrowise \
pip3 install -r requirements.txt
```

## Overview of the workflow

We follow the OpenNMT architecture to train the Transformer. The workflow is

* generate the dataset with SMILES augmentation;
* build the vocab based on the dataset;
* Pretrain the model;
* train a base model with the dataset
* augment the dataset with in-silico reactions generated from the base model
* start training a more powerful model;
* averaging the checkpoints to get the final checkpoint.
* choose a checkpoint to make model predictions.
* score the prediction to get the accuracy.

All the config files are placed in the `pretrain_finetune` folder. Using OpenNMT commands to run the codes and modifiing them according to the needs.

## Data preprocessing

Download the datasets and put them in the `dataset` directory like this:

```shell
USPTO-50K: dataset/USPTO_50K/raw_train.csv
USPTO-MIT: dataset/USPTO-MIT/train.txt
USPTO_full: dataset/USPTO_full/raw_train.csv
```

Generate the pretrain data with SMILES augmentation:

```
python preprocessing/generate_pretrain_data.py -mode product -augmentation 5
python preprocessing/generate_masked_pretrain_data.py -dataset USPTO_full_pretrain_aug5_product
```

Generate the dataset with R-SMILES augmentation:

```
python preprocessing/generate_PtoR_data.py -dataset USPTO_50K -augmentation 20 -processes 8
python preprocessing/generate_PtoR_data.py -dataset USPTO-MIT -augmentation 5 -processes 8
python preprocessing/generate_PtoR_data.py -dataset USPTO_full -augmentation 5 -processes 8
```

Similar operations (need to modify some codes of generate_PtoR_data.py) for Rare-subsets:

```
Rare2: dataset/low-frequency/Rare2.csv
Rare5: dataset/low-frequency/Rare5.csv
Rare10: dataset/low-frequency/Rare10.csv
```


## Pretrain 

​(**Generate all the datasets in advance to build a full vocab before pretraining the model.**)

Just run the prepared shell command to start pretraining according to your pretrain dataset.

  ```shell
  bash shell_cmds/product_pretrain.sh
  bash shell_cmds/reactant_pretrain.sh
  bash shell_cmds/both_pretrain.sh
  ```
You can get the pretrained model from [here](https://figshare.com/articles/dataset/RetroWISE-data_zip/25272451).

## Train a base model 

Run the OpenNMT train command with the prepared config. 

  ```shell
  onmt_train -config pretrain_finetune/finetune/RtoP/RtoP-50K-aug20-config.yml
  bash pretrain_finetune/finetune/RtoP/RtoP-50K-aug20-average.sh
  ```

## Generate in-silico data

Run the bash script to generate in silico reactions.

 ```shell
 onmt_translate -config pretrain_finetune/backT/RtoP-t/RtoP-extra-50K-aug5-translate.yml
 ```
 
## Filtering

(step 1) Template matching

```
cd templates
bash match-template.sh
```

(step 2) Molecular Similarity Comparison

```
cd filter-post
bash compare.sh
```
Alternatively, you can get the generated in-silico data from [here](https://figshare.com/articles/dataset/RetroWISE-data_zip/25272451).

## Retrain the model

Place the generated data in the correct directory. 

Run the bash script to get the final prediction score.

```shell
bash btbash/50k.sh
```
## Acknowledgement

OpenNMT-py: https://github.com/OpenNMT/OpenNMT-py

R-SMILES: https://github.com/otori-bird/retrosynthesis