# TCRclub
a tool to clustering T cells by integrating sc-RNA seq and sc-TCR seq on local harmony

## Introduction
TCRclub is a novel approach that identifies the functional relevance of T cells at single-cell resolution. TCRclub receives scRNA expression and the numeric embeddings of the CDR3β sequences as inputs. It aims to bridge the gap between scRNA expression and TCRs by focusing on the inner relationship between T cells with similar functions. To achieve this, TCRclub models the relationship between pair-wise TCR embedding and pair-wise expression distances according to the local harmony. Local harmony means the nearby homogeneity existing in the local neighbours of any cell, since the neighbouring cells in the distance space are more likely to have similar characteristics and belong to the same category. By emphasizing local harmony, TCRclub reduces noise and increases the robustness of integration. Considering the built-in cell structure, TCRclub builds the T-cell hierarchy based on the distances updated by the integration and extracts the T-cell clubs. Finally, TCRclub repeats multiple times to obtain the consensus results of the clubs as the final output. 

## Prerequisite
TCRclub is implemented in Python and requires a GPU for the acceleration. 

## Environment Setup
For optimal performance, we recommend using CUDA versions 11.2 or 11.6 along with cudnn8.1.0. Ensure that your CUDA environment supports TensorFlow versions 2.6 to 2.9.

To avoid conflicts with your existing environment, we suggest creating a new Anaconda environment. There are several ways to set up the environment:

1. ### Use Provided Conda Environment
You can directly utilize the provided conda environment using the following instructions.

```
# Navigate to your anaconda3 directory 
cd /home/XX/anaconda3/envs/
 
mkdir -p TCRclub
tar zxvf /<the path you store the provided .tar.gz>/TCRclubConda.tar.gz -C TCRclub
source /home/XX/anaconda3/envs/TCRclub/bin/activate
conda-unpack

# Check the new environment
conda env list
## Output
conda environments:
base    */home/XX/anaconda3
TCRclub /home/XX/anaconda3/envs/TCRclub
```

2. ### Create Conda Environment
* Clone the repository.
* Create a conda environment with python3.8 or python3.9, for example `conda create -n <Environment Name> python=3.8`
* Activate the conda environment you just created. Navigate to the TCRclub directory and execute [install.sh](./TCRclub/install.sh).
3. ### Use Your Existing Environment
You have the option to install the required Python packages manually in your existing environment. However, for successful execution of TCRclub, we recommend installing the following Python packages in sequential order along with their respective versions:

* pyseat (version 0.0.1.4)
* tensorflow (version 2.6.0~2.9.0)
* torch (version 1.13.0)
* numba (version 0.59.0)
* networkx (version 2.8.8)
  
[pySEAT](https://github.com/deepomicslab/SEAT) can be installed by the following command:

`pip install pyseat`

## Running TCRclub
TCRclub includes a VAE model for converting CDR3β sequences into embeddings. To ensure the proper functioning of TCRclub, please download the "autoencoder" folder and the Python script TCRclub.py if you do not want to clone the whole repository.

### Input data
TCRclub requires two files to identify functional-similar clubs for T cells. The first column of both files should be named "barcode", with each element being unique. Each row represents a single T cell, and the order of T cells should be the same in both files.

For the sc-TCR sequences file, two columns are necessary: "barcode" and "cdr3". If the file contains cells from multiple samples, you should include a column labeled "sample" to indicate the source of each T cell. Because of the structure of our autoencoder, CDR3 longer than 30 should be removed before using TCRclub.
![Image text](https://github.com/deepomicslab/TCRclub/blob/e48f9a7903811dd043a2e26f4402704a68c69bb1/img/required_tcr_file.png)

**Fig.1** An example of required scTCR file in .csv format.

For the sc-RNA expression file, each row corresponds to a T cell, and the columns (except the first column) correspond to genes. We suggest using the top 10% highly expressed genes extracted from the original sc-RNA expression file. You can select the input genes according to your own criteria. Normalization and log-transformation are recommended. Be cautious with selecting too many genes, as it may cause GPU memory issues.
![Image text](https://github.com/deepomicslab/TCRclub/blob/e48f9a7903811dd043a2e26f4402704a68c69bb1/img/required_rna_file.png)

**Fig.2** An example of required scRNA expression file in .csv format.

### Parameters
TCRclub accepts several parameters, as listed in the table below:
| Parameters  | Description |
| ------------- | ------------- |
| tcr_file | **(Required)** The path of the .csv file contains the scTCR sequences (see **Input data**). |
| rna_file | **(Required)** The path of the .csv file contains the scRNA expression (see **Input data**). |
| k | Number of selected neareast neighbours. Default: 10.|
| repeat_times  | Repeat times for obtaining the consensus results. Default: 50. Set to 1 if using *fixed_initialization*.|
| beta  | L2 regularization parameter. Default: 1e-7. |
| single_cutoff  | Cut-off parameter to split the cell hierarchy in a single run. It can be increased(decreased) for a higher(lower) clustering coverage of the individual result. Default: 1e-4. |
| con_cutoff  | Cut-off parameter to split the cell hierarchy based on the consensus matrix. It can be increased(decreased) for a higher(lower) clustering coverage of the consensus result. Default: 5e-4. |
| con_topk  | Parameter to choose the number of results with the smallest loss from the repeat_times results to produce the consensus matrix. Default: 15. |
| out | Output directory name. |
| multiple_sample  | A binary value indicating whether the input T cells are derived from different samples. If this parameter is selected, the input TCR file should contain a column specified as "sample". Default: False|
| fixed_initialization | A binary value indicating whether the initialization of TCRclub starts in the default way (randomness). If this parameter is selected, the initialization of matrix C in TCRclub will be fixed. In this case, we suggest the parameter *repeat_times* should be set as 1. Default: False.|

### Identifying T-cell clubs
To identify T-cell clubs, follow the instructions below based on your specific scenario:

If the cells come from a single sample/patient, run the Python script [TCRclub.py](TCRclub.py) using the following command:

`python3 TCRclub.py --tcr_file ./example_data/processed_tcr.csv --rna_file ./example_data/processed_rna.csv`

If the cells come from multiple samples/patients, run python script [TCRclub.py](TCRclub.py) with with the additional *--multiple_sample* parameter:

`python3 TCRclub.py --tcr_file ./example_data/processed_tcr.csv --rna_file ./example_data/processed_rna.csv --multiple_sample`

If you choose the fixed_initialization option, run python script [TCRclub.py](TCRclub.py) with the following command.

`python3 TCRclub.py --tcr_file ./example_data/processed_tcr.csv --rna_file ./example_data/processed_rna.csv --fixed_initialization --repeat_times 1`

The T-cell clubs will be saved in the "consensus_result.csv" file within the output directory specified by the *out* parameter. The output file will include a new column named "club" in the input TCR file, where T cells with the same club ID are considered to belong to the same club.
![Image text](https://github.com/deepomicslab/TCRclub/blob/e48f9a7903811dd043a2e26f4402704a68c69bb1/img/example_result_file.png)
**Fig.3** An example of the produced result in .csv format.

