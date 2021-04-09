# Inference and analysis of allele-specific CNAs, LOH, and WGD

Allele-specific CNAs, LOHs, and WGDs have been inferred using [HATCHet](https://github.com/raphael-group/hatchet) and the inferred results for all analyzed samples are reported here in the folder `hatchet`, with a folder for each tumour.
Note that all samples from the same tumour are analyzed jointly and therefore the corresponding results are in the same files.
This repository includes a Jupyter notebook [alleleCNAs_LOH_WGD](./alleleCNAs_LOH_WGD.ipynb) to reproduce all the related analyses and plots.

## Requirements

The Jupyter notebook requires `python2.7` and the following packages: `numpy`, `pandas`, `matplotlib`, `seaborn`, and `scipy`. Moreover, [jupyter](https://jupyter.org/) is required to execute the notebook.
For the sake of space limitations, all the HATCHet processed results that are included in this repository have been zipped using `gzip`.
The execution of the notebook thus requires the user to unpack all files with the following command (executing from this directory):

```shell
find hatchet/ -type f -name 'hatchet*.gz' -exec gzip -d {} +
```

## Installation

The following commands are sufficient to install all the requirements within a [conda](https://docs.conda.io/en/latest/) environment without any requirement in any *nix system (if you are in OSx please substitute the downloading link as appropriate from [miniconda](https://docs.conda.io/en/latest/miniconda.html)) when executed from this directory:

```shell
curl -L https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh > miniconda.sh
rm -rf ./conda/
bash miniconda.sh -b -f -p ./conda/
conda/bin/conda create -n hatchet-pdx python=2.7 jupyter numpy pandas matplotlib seaborn scipy -y
source conda/bin/activate hatchet-pdx
```

Before any re-execution in a new session, please only run the following command from this folder:

```shell
source conda/bin/activate hatchet-pdx
```

## Results

The resulting plots are generated in this directory as PNG images after running the notebook but they can also simply be visualized by clicking on [alleleCNAs_LOH_WGD](./alleleCNAs_LOH_WGD.ipynb) without executing the script.

## Contact

Author: Dr Simone Zaccaria\
Old affilliation: Princeton University, NJ (USA)\
New affilliation: UCL Cancer Institute, London (UK)\
Correspondence: s.zaccaria@ucl.ac.uk\
Website: [www.ucl.ac.uk/cancer/zaccaria-lab](https://www.ucl.ac.uk/cancer/zaccaria-lab)
