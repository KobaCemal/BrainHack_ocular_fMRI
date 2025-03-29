import ants
import numpy as np
import pandas as pd

# Load fMRI image
fmri_img = ants.image_read("fmri.nii.gz")

# Extract the first time point (volume) as the reference image
reference_img = fmri_img[:,:,:,0]  # First 3D volume

# Perform rigid-body registration
motion_corrected = ants.registration(fixed=reference_img, moving=fmri_img, type_of_transform="Rigid")

# Save the corrected fMRI image
ants.image_write(motion_corrected['warpedmovout'], "fmri_corrected.nii.gz")

# Extract motion parameters (affine transformation matrix)
affine_matrix = motion_corrected['fwdtransforms'][0]  # Forward transformation file

# Load the affine transformation matrix
affine_values = ants.read_transform(affine_matrix).tolist()

# Convert to readable format (translation + rotation)
motion_params = {
    "tx": affine_values[0][-1],  # Translation X
    "ty": affine_values[1][-1],  # Translation Y
    "tz": affine_values[2][-1],  # Translation Z
    "rx": affine_values[0][:3],  # Rotation row 1
    "ry": affine_values[1][:3],  # Rotation row 2
    "rz": affine_values[2][:3]   # Rotation row 3
}

# Convert to DataFrame for easy export
df = pd.DataFrame([motion_params])

# Save motion parameters as CSV
df.to_csv("motion_parameters.csv", index=False)

print("Motion correction done. Motion parameters saved to motion_parameters.csv")
