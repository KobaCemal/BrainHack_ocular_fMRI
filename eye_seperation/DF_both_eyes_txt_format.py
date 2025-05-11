# ----------------------------------------------
# This code is for files that is in txt format.
# NOT FOR Niftii files. 
# ----------------------------------------------
# This file splits the data on patients 
# control and pathological group and left/right eye.
# ----------------------------------------------

import nibabel as nib
import numpy as np
import pandas as pd
from dipy.align.imaffine import (MutualInformationMetric, AffineRegistration)
from dipy.align.transforms import RigidTransform3D
from dipy.viz import regtools
from scipy.spatial.transform import Rotation as R
import os, glob

# file path
data_path = "folder_path"

# eye deviation
left_eye = "left"
right_eye = "right"
eyes = [left_eye, right_eye]

# file search engine
right_files = glob.glob(os.path.join(data_path, "sub-*_run-0_eyes_only_mc_mc_right.txt"))
subject_ids = sorted({ os.path.basename(f).split("_")[0] for f in right_files })

# result list
all_motion_params = []

for sub in subject_ids:
        for eye in eyes:
                # specifying file path
                file_path = os.path.join(data_path, f"{sub}_run-0_eyes_only_mc_mc_{eye}.txt")
                if not os.path.exists(file_path):
                        print(f"file not exist: {file_path}")
                        continue
                
                print(f"Processing: {file_path}")
                
                # Load data from txt extension files
                data = pd.read_csv(file_path, delim_whitespace=True, header=None, )
                data.columns = ["tx", "ty", "tz", "rx", "ry", "rz"]
                data["time"] = data.index # naming the time step
                data["eye"] = eye # specifying the eye
                data["subject"] = sub # subject number
                # sorting the columns
                data = data[["subject", "eye", "time", "tx", "ty", "tz", "rx", "ry", "rz"]]

                # DataFrame output
                out_file = os.path.join(data_path, f"group_motion_params_{sub}_{eye}.csv")
                all_motion_params.append(data)

df_all = pd.concat(all_motion_params, ignore_index=True)
df_all.to_csv(os.path.join(data_path, "motion_params_all_eyes.csv"), index=False)

print("Motion correction done. Motion parameters saved to motion_parameters_DF.csv")
