#!/usr/bin/env bash
set -e # Exit the script if anything fails
# This script builds the basic hbdx enviroment

# CONFIGURATION
#------------------------------------------------------------------------------
# The path to the users bash config file
bash_rc_file=$HOME/.bashrc

# Name of the enviroment to be installed
env_name=hbdx-scgen

# Python version to be installed
python_version="3.7"
#------------------------------------------------------------------------------

add_to_bashrc () {
    # Add a line to the users bash_rc file if the line is not present yet
    grep -qF -- "$1" "$bash_rc_file" || echo "$1" >> "$bash_rc_file"
}

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=MacOSX;;
esac

for cdir in anaconda3 miniconda3 conda; do
    [ -f $HOME/$cdir/bin/conda ] && conda_bin=$HOME/$cdir/bin/conda
    [ -f $HOME/opt/$cdir/bin/conda ] && conda_bin=$HOME/opt/$cdir/bin/conda
done

# If conda is not installed, install miniconda!
if [[ $conda_bin == "" ]]; then
    echo "Conda is not installed, will install miniconda3"
    target="Miniconda3-latest-${machine}-x86_64.sh"
    wget "https://repo.anaconda.com/miniconda/${target}"
    bash $target -b -p $HOME/conda
    rm $target

    # Add init code to bashrc
    $HOME/conda/bin/conda init
fi

echo "Anaconda install is ${conda_bin}"
eval "$(${conda_bin} shell.bash hook)"

# Remove the hbdx enviroment if already installed
if conda env list | grep -q $env_name; then
    conda deactivate
    conda env remove --name $env_name
fi

# Create a new enviroment called $env_name
conda create -n $env_name python="${python_version}" -y
eval "$(${conda_bin} shell.bash hook)"
conda activate $env_name

# Update bash.rc...................................................
add_to_bashrc "conda activate ${env_name}"
add_to_bashrc "export PYTHONBREAKPOINT=ipdb.set_trace"

conda env list

# Install the mamba package manager
conda install mamba -c conda-forge -y

# Installing pakages by repository source
# PyTorch..........................................................
if [[ $machine == "MacOSX" ]]; then
        mamba install https://conda.anaconda.org/pytorch/osx-64/pytorch-1.7.0-py3.8_0.tar.bz2
else
    if hash nvidia-smi &>/dev/null ; then
        mamba install -c pytorch -y pytorch==1.7.0 cudatoolkit=10.2
        # Hard-fixing for now as the repodata.json was offline
        #mamba install https://conda.anaconda.org/pytorch/linux-64/pytorch-1.7.0-py3.8_cuda10.2.89_cudnn7.6.5_0.tar.bz2
    else
        mamba install -c pytorch -y pytorch==1.7.0 cpuonly
    fi
fi

#mamba install scvi-tools -c conda-forge -c bioconda

# Conda-Forge......................................................
mamba install -c conda-forge -y python"=3.7" h5py"=2.10" rpy2 mkl conda-build colorama ipdb isort pandas pdoc3 rich scikit-learn scipy seaborn umap-learn tqdm pylint pytest black anndata">=0.7.5" dill hydra-core jupyter rope mpmath louvain python-xxhash python-annoy

# BioConda.........................................................
#mamba install -c bioconda -y harmonypy scanpy bbknn

# PIP..............................................................
pip install combat harmonypy scanorama scanpy bbknn tensorflow-gpu==1.15 keras flit biopython
#rm -rf scGen
#git clone https://github.com/theislab/scGen
#cd scGen
#flit install
#cd ..

pip install https://github.com/theislab/scGen


# Set the local hummingbird repo into development mode.............
if [[ -d .git ]] && [[ -d hummingbird ]]; then
    conda develop "${PWD}"
fi
