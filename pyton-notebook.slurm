#!/bin/bash
#SBATCH --partition=plgrid-gpu-a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=15:00:00
#SBATCH --job-name=jupyter-notebook-gpu
#SBATCH --output=jupyter-log-%J.txt
#SBATCH --error=jupyter-err-%J.txt
#SBATCH --gpus=2

set -e  # stop on first error

## -------------------------
## CONFIG
VENV_DIR=/net/tscratch/people/plgkoba/venvs/deepmreye-gpu
REQ_FILE=/net/tscratch/people/plgkoba/requirements.txt
PYTHON_MODULE=Python/3.10.8
CUDA_MODULE=CUDA/12.2
CUDNN_MODULE=cuDNN/8.9.2.26-CUDA-12.2.0

## -------------------------
## SSH tunneling info
XDG_RUNTIME_DIR=""
ipnport=$(shuf -i8000-9999 -n1)
ipnip=$(hostname -i)
user=$USER

echo -e "
Copy/Paste this in your local terminal to SSH tunnel with remote
-----------------------------------------------------------------
ssh -o ServerAliveInterval=300 -N -L $ipnport:$ipnip:$ipnport ${user}@athena.cyfronet.pl
-----------------------------------------------------------------

Then open a browser on your local machine to the following address
------------------------------------------------------------------
http://localhost:$ipnport
------------------------------------------------------------------
"

## -------------------------
## LOAD PYTHON + GPU MODULES
module purge
module load GCCcore/12.2.0
module load $PYTHON_MODULE
module load $CUDA_MODULE
module load $CUDNN_MODULE

# Fix for libdevice not found (XLA GPU)
export XLA_FLAGS=--xla_gpu_cuda_data_dir=/net/software/v1/software/CUDA/12.2.0

## -------------------------
## CREATE OR ACTIVATE VENV
create_venv=0
if [ -d "$VENV_DIR" ]; then
    VENV_PYTHON="$VENV_DIR/bin/python3"
    if [ -x "$VENV_PYTHON" ]; then
        VENV_VER=$($VENV_PYTHON -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
        MODULE_VER=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
        if [ "$VENV_VER" != "$MODULE_VER" ]; then
            echo ">>> Python version mismatch. Recreating venv."
            create_venv=1
            rm -rf $VENV_DIR
        fi
    else
        create_venv=1
    fi
else
    create_venv=1
fi

if [ $create_venv -eq 1 ]; then
    echo ">>> Creating virtualenv at $VENV_DIR"
    python3 -m venv $VENV_DIR
fi

echo ">>> Activating virtualenv"
source $VENV_DIR/bin/activate
pip install --upgrade pip

## -------------------------
## INSTALL REQUIRED PACKAGES
pip install tensorflow==2.15.0 deepmreye
if [ -f "$REQ_FILE" ]; then
    pip install --no-cache-dir -r $REQ_FILE
fi

# Ensure Jupyter + widgets are installed
pip install notebook jupyterlab ipywidgets anywidget ipykernel

## -------------------------
## INSTALL KERNEL
python -m ipykernel install --user --name=deepmreye-gpu --display-name "Python (deepmreye-gpu)"

## -------------------------
## START JUPYTER NOTEBOOK
jupyter-notebook --no-browser --port=$ipnport --ip=$ipnip
