# Libreries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_all = pd.read_csv("motion_params_all_eyes.csv")

# Left and right eye
df_left = df_all[df_all["eye"] == "left"].set_index("time")[["rx", "ry", "rz"]]
df_right = df_all[df_all["eye"] == "right"].set_index("time")[["rx", "ry", "rz"]]

# data for the same time point
df_merged = df_left.join(df_right, lsuffix="_left", rsuffix="_right")

# Euclidean distance
df_merged["euclidean_distance"] = np.sqrt(
    (df_merged["rx_left"] - df_merged["rx_right"])**2 +
    (df_merged["ry_left"] - df_merged["ry_right"])**2 +
    (df_merged["rz_left"] - df_merged["rz_right"])**2
)

# CSV file
df_merged.to_csv("euclidean_distance_rotation.csv")

# Graph
plt.figure(figsize=(10, 6))
plt.plot(df_merged.index, df_merged["euclidean_distance"], label="Euclidean Distance", color="red")
plt.plot()
plt.xlabel("Time (volume number)")
plt.ylabel("Euclidean Distance")
plt.title("Euclidean Distance between Left and Right Eye Rotations")
plt.legend()
plt.grid(True)
plt.show()
