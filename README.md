# TCRclub
a tool to clustering T cells by integrating sc-RNA seq and sc-TCR seq on local harmony

## Introduction
TCRclub is a novel approach that identifies the functional relevance of T cells at single-cell resolution. TCRclub receives scRNA expression and the numeric embeddings of the CDR3β sequences as inputs. It aims to bridge the gap between scRNA expression and TCRs by focusing on the inner relationship between T cells with similar functions. To achieve this, TCRclub models the relationship between pair-wise TCR embedding and pair-wise expression distances according to the local harmony. Local harmony means the nearby homogeneity existing in the local neighbours of any cell, since the neighbouring cells in the distance space are more likely to have similar characteristics and belong to the same category. By emphasizing local harmony, TCRclub reduces noise and increases the robustness of integration. Considering the built-in cell structure, TCRclub builds the T-cell hierarchy based on the distances updated by the integration and extracts the T-cell clubs. Finally, TCRclub repeats multiple times to obtain the consensus results of the clubs as the final output. 

## Prerequisite
TCRclub is implemented in Python and requires a GPU. Please install the GPU version of PyTorch and TensorFlow before proceeding. We recommend executing the Python scripts using Linux shell commands.

### Python Packges
To run TCRclub successfully, we recommend installing the following Python packages with their respective versions:

* tensorflow (version 2.5.0)
* torch (version 1.12.0)
* pandas (version 1.3.5)
* scanpy (version 1.9.1)
* numba (version 0.54.1)
* pyseat (version 0.0.1.3)

[pySEAT](https://github.com/deepomicslab/SEAT) can be installed by the following command:

`pip install pyseat`

## Running TCRclub
TCRclub includes a VAE model for converting CDR3β sequences into embeddings. To ensure the proper functioning of TCRclub, please download the "autoencoder" folder and the Python script TCRclub.py if you do not want to clone the whole repository.

### Input data
TCRclub requires two files to identify functional-similar clubs for T cells. The first column of both files should be named "barcode", with each element being unique. Each row represents a single T cell, and the order of T cells should be the same in both files.

For the sc-TCR sequences file, two columns are necessary: "barcode" and "cdr3". If the file contains cells from multiple samples, you should include a column labeled "sample" to indicate the source of each T cell.
![Image text](https://github.com/deepomicslab/TCRclub/blob/ab0b9ad869b7b0a47ab6fc95d23e34e3b3960bb4/img/required_tcr_file.png)

**Fig.1** An example of required scTCR file in .csv format.

For the sc-RNA expression file, each row corresponds to a T cell, and the columns (except the first column) correspond to genes. We suggest using the top 10% highly expressed genes extracted from the original sc-RNA expression file. You can select the input genes according to your own criteria. Normalization and log-transformation are recommended. Be cautious with selecting too many genes, as it may cause GPU memory issues.
![Image text](https://github.com/deepomicslab/TCRclub/blob/ab0b9ad869b7b0a47ab6fc95d23e34e3b3960bb4/img/required_rna_file.png)

**Fig.2** An example of required scRNA expression file in .csv format.

### Parameters
TCRclub accepts several parameters, as listed in the table below:
| Parameters  | Description |
| ------------- | ------------- |
| tcr_file | **(Required)** The path of the .csv file contains the scTCR sequences (see **Input data**). |
| rna_file | **(Required)** The path of the .csv file contains the scRNA expression (see **Input data**). |
| k | Number of selected neareast neighbours. Default: 5.|
| repeat_times  | Repeat times for obtaining the consensus results. Default: 50. Set to 1 if using *fixed_initialization*.|
| alpha  | L2 regularization parameter. Default: 1e-7. |
| beta  | L2 regularization parameter. Default: 1e-7. |
| single_cutoff  | Cut-off parameter to split the cell hierarchy in a single run. Default: 1e-4. |
| con_cutoff  | Cut-off parameter to split the cell hierarchy based on the consensus matrix . Default: 2e-3. |
| con_topk  | Parameter to choose the number of results with the smallest loss from the repeat_times results to produce the consensus matrix. Default: 25. |
| out | Output directory name. |
| multiple_sample  | A binary value indicating whether the input T cells are derived from different samples. If this parameter is selected, the input TCR file should contain a column specified as "sample". Default: False|
| fixed_initialization | A binary value indicating whether the initialization of TCRclub starts in the default way (randomness). If this parameter is selected, the initialization of matrix C in TCRclub will be fixed. In this case, we suggest the parameter *repeat_times* should be set as 1. Default: False.|

### Identifying T-cell clubs
To identify T-cell clubs, follow the instructions below based on your specific scenario:

If the cells come from a single sample/patient, run the Python script [TCRclub.py](TCRclub.py) using the following command:

`python3 TCRclub.py --tcr_file ./example_data/example_tcr.csv --rna_file ./example_data/example_rna.csv`

If the cells come from multiple samples/patients, run python script [TCRclub.py](TCRclub.py) with with the additional *--multiple_sample* parameter:

`python3 TCRclub.py --tcr_file ./example_data/example_tcr.csv --rna_file ./example_data/example_rna.csv --multiple_sample`

If you choose the fixed_initialization option, run python script [TCRclub.py](TCRclub.py) with the following command.

`python3 TCRclub.py --tcr_file ./example_data/example_tcr.csv --rna_file ./example_data/example_rna.csv --fixed_initialization --repeat_times 1`

The T-cell clubs will be saved in the "consensus_result.csv" file within the output directory specified by the *out* parameter. The output file will include a new column named "club" in the input TCR file, where T cells with the same club ID are considered to belong to the same club. 
