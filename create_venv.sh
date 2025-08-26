#!/bin/bash
#SBATCH --partition=plgrid-gpu-a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --job-name=create-venv
#SBATCH --output=create-venv-%J.out
#SBATCH --error=create-venv-%J.err

set -e  # stop on first error

## -------------------------
## CONFIG
VENV_DIR=/net/tscratch/people/plgkoba/venvs/deepmreye
REQ_FILE=/net/tscratch/people/plgkoba/requirements.txt
LOG_FILE=/net/tscratch/people/plgkoba/logs/create_venv.log

## -------------------------
## PREPARE DIRECTORIES
mkdir -p $(dirname $VENV_DIR)   # parent dir for venv
mkdir -p $(dirname $LOG_FILE)   # parent dir for logs

echo ">>> Starting venv creation at $(date)" | tee $LOG_FILE

## -------------------------
## LOAD PYTHON MODULE
module purge
module load GCCcore/12.2.0
module load Python/3.10.8

## -------------------------
## REMOVE OLD VENV IF EXISTS
if [ -d "$VENV_DIR" ]; then
    echo ">>> Removing old virtualenv at $VENV_DIR" | tee -a $LOG_FILE
    rm -rf $VENV_DIR
fi

## -------------------------
## CREATE NEW VENV
echo ">>> Creating new virtualenv at $VENV_DIR" | tee -a $LOG_FILE
python3 -m venv $VENV_DIR 2>&1 | tee -a $LOG_FILE

## -------------------------
## ACTIVATE VENV
source $VENV_DIR/bin/activate

## -------------------------
## UPGRADE PIP
echo ">>> Upgrading pip" | tee -a $LOG_FILE
pip install --upgrade pip 2>&1 | tee -a $LOG_FILE

## -------------------------
## INSTALL REQUIREMENTS
if [ -f "$REQ_FILE" ]; then
    echo ">>> Installing packages from $REQ_FILE" | tee -a $LOG_FILE
    pip install --no-cache-dir -r $REQ_FILE 2>&1 | tee -a $LOG_FILE
else
    echo ">>> ERROR: requirements.txt not found at $REQ_FILE" | tee -a $LOG_FILE
    exit 1
fi

## -------------------------
## ENSURE JUPYTER + IPYKERNEL
echo ">>> Installing Jupyter + ipykernel" | tee -a $LOG_FILE
pip install notebook ipykernel ipywidgets 2>&1 | tee -a $LOG_FILE

## -------------------------
## REGISTER JUPYTER KERNEL
echo ">>> Registering Jupyter kernel" | tee -a $LOG_FILE
python -m ipykernel install --user --name=deepmreye --display-name "Python (deepmreye)" 2>&1 | tee -a $LOG_FILE

echo ">>> Virtual environment setup complete at $(date)" | tee -a $LOG_FILE
