# --- Libreries
import nibabel as nib
import numpy as np
import pandas as pd
from dipy.align.imaffine import (MutualInformationMetric, AffineRegistration)
from dipy.align.transforms import RigidTransform3D
from dipy.viz import regtools
from scipy.spatial.transform import Rotation as R

# eye deviation
left_eye = "left"
right_eye = "right"
eyes = [left_eye, right_eye]

# list
all_motion_params = []

for eye in eyes:
        file_path = nib.load(f"/Users/mikolajturczyniak/Desktop/connected/Visual_Studio/python/BrainHack/example data/{eye}_eye_time.nii.gz")

        # Load fMRI data (assumes 4D NIfTI file)
        fmri_data = file_path.get_fdata()

        # Extract the first time point as reference
        reference_img_data = fmri_data[:, :, :, 0]  

        # Set up the registration
        n_volumes = fmri_data.shape[3]
        corrected_data = np.zeros_like(fmri_data)
        motion_params_list = {"time": [], "tx": [], "ty": [], "tz": [], "rx": [], "ry": [], "rz": []}
        metric = MutualInformationMetric(nbins=32)  # Mutual information metric
        reg = AffineRegistration()  # Register images
        rigid_transform = RigidTransform3D()  # Rigid transformation (6 DOF)
        params0 = np.array([0, 0, 0, 0, 0, 0]) # Starting point of eye analysis

        for t in range(n_volumes):
                moving_img_data = fmri_data[:, :, :, t]
                # transformation
                transform = reg.optimize(static=reference_img_data,
                        moving=moving_img_data,
                         transform=rigid_transform,
                         params0=params0)

                # transformation on volume
                corrected_volume = transform.transform(moving_img_data)
                corrected_data[:, :, :, t] = corrected_volume
                affine_matrix = transform.affine # matrix 4x4
        
                # translation
                translation = affine_matrix[:3, 3]
        
                # rotation extraction 3x3
                rotation_matrix = affine_matrix[:3, :3]
                r = R.from_matrix(rotation_matrix)
                euler_angles = r.as_euler('xyz', degrees=True)
        
                # translation and rotation in one vector
                motion_params_list["time"].append(t)
                motion_params_list["tx"].append(translation[0])
                motion_params_list["ty"].append(translation[1])
                motion_params_list["tz"].append(translation[2])
                motion_params_list["rx"].append(euler_angles[0])
                motion_params_list["ry"].append(euler_angles[1])
                motion_params_list["rz"].append(euler_angles[2])

        # DataFrame
        df = pd.DataFrame(motion_params_list)
        df["eye"] = eye
        df = df[["eye", "time", "tx", "ty", "tz", "rx", "ry", "rz"]]
        df.to_csv(f"motion_params_{eye}.csv", index=False)
        all_motion_params.append(df)

df.to_csv("motion_params_DF.csv", index=False)
df_all = pd.concat(all_motion_params, ignore_index=True)
df_all.to_csv("motion_params_all_eyes.csv", index=False)

print("Motion correction done. Motion parameters saved to motion_parameters_DF.csv")
