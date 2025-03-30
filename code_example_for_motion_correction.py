import nibabel as nib
import numpy as np
import pandas as pd
from dipy.align.imaffine import (MutualInformationMetric, AffineRegistration)
from dipy.align.transforms import RigidTransform3D
from dipy.viz import regtools

# Load fMRI data (assumes 4D NIfTI file)
fmri_img = nib.load("fmri.nii.gz")
fmri_data = fmri_img.get_fdata()

# Extract the first time point as reference
reference_img_data = fmri_data[:,:,:,0]  # First time point (volume)

# Set up the registration
metric = MutualInformationMetric(nbins=32)  # Mutual information metric
reg = AffineRegistration()  # Register images
rigid_transform = RigidTransform3D()  # Rigid transformation (6 DOF)

# Perform rigid-body registration
params0 = None  # Initial parameters (identity transformation)
transform = reg.optimize(static=reference_img_data, moving=fmri_data, transform=rigid_transform, params0=params0)

# Extract motion parameters (translation and rotation)
motion_params = {
    "tx": transform.params[0],  # Translation in X
    "ty": transform.params[1],  # Translation in Y
    "tz": transform.params[2],  # Translation in Z
    "rx": transform.params[3],  # Rotation around X-axis
    "ry": transform.params[4],  # Rotation around Y-axis
    "rz": transform.params[5]   # Rotation around Z-axis
}

# Convert motion parameters to a DataFrame
df = pd.DataFrame([motion_params])

# Save motion parameters to CSV
df.to_csv("motion_parameters.csv", index=False)

# Apply the transformation to all volumes
corrected_data = transform.transform(fmri_data)

# Save the corrected fMRI image
corrected_img = nib.Nifti1Image(corrected_data, fmri_img.affine)
nib.save(corrected_img, "fmri_corrected.nii.gz")

print("Motion correction done. Motion parameters saved to motion_parameters.csv")
