#!/bin/bash

# ------------------------------
# CONFIG
# ------------------------------
VENV_DIR=~/venvs/deepmreye-gpu
server_address="CemalKoba@XXXXXX"
XDG_RUNTIME_DIR=""

# ------------------------------
# 1️Create or activate venv
# ------------------------------
if [ -d "$VENV_DIR" ]; then
    echo "Activating existing venv at $VENV_DIR"
else
    echo "Creating venv at $VENV_DIR"
    mkdir -p ~/venvs
    python3 -m venv $VENV_DIR
fi

source $VENV_DIR/bin/activate

# ------------------------------
# 2️ Upgrade pip
# ------------------------------
pip install --upgrade pip

# ------------------------------
# 3️ Install required packages
# ------------------------------
REQUIRED_PACKAGES=(
    tensorflow==2.15.0
    deepmreye==0.3
    numpy==1.26.4
    scipy==1.13.1
    pandas==2.3.1
    scikit-learn
    plotly
    matplotlib
    ipykernel
    ipywidgets
    notebook
)

for pkg in "${REQUIRED_PACKAGES[@]}"; do
    pip show "$pkg" >/dev/null 2>&1 || pip install "$pkg"
done

# ------------------------------
# 4️ Register kernel if not present
# ------------------------------
if ! jupyter kernelspec list | grep -q "deepmreye-gpu"; then
    python -m ipykernel install --user --name=deepmreye-gpu --display-name "Python (deepmreye-gpu)"
fi

# ------------------------------
# 5️ Pick a random port
# ------------------------------
ipnport=$(shuf -i8000-9999 -n1)
ipnip=$(hostname -i)
user=$USER

echo -e "
Copy/Paste this in your local terminal to SSH tunnel with remote
-----------------------------------------------------------------
ssh -o ServerAliveInterval=300 -N -L $ipnport:$ipnip:$ipnport ${server_address}
-----------------------------------------------------------------

Then open a browser on your local machine to the following address
------------------------------------------------------------------
http://localhost:$ipnport
------------------------------------------------------------------
"
# ------------------------------
# 6️ Start Jupyter Notebook with cleanup
# ------------------------------
# ------------------------------
# 6️ Start Jupyter Notebook in /mnt/compneuro
# ------------------------------
NOTEBOOK_DIR=/mnt/compneuro/
trap "echo 'Stopping Jupyter Notebook'; pkill -f 'jupyter-notebook'; exit" SIGINT SIGTERM

jupyter-notebook --notebook-dir=$NOTEBOOK_DIR --no-browser --port=$ipnport --ip=$ipnip



