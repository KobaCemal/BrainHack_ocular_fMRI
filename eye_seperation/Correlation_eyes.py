# Libreries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df_all = pd.read_csv("motion_params_all_eyes.csv")

# separate data
df_left = df_all[df_all["eye"] == "left"].set_index("time")[["rx"]]
df_right = df_all[df_all["eye"] == "right"].set_index("time")[["rx"]]

# data merging
df_merged = df_left.join(df_right, lsuffix="_left", rsuffix="_right")

# Correlation
Correlation = df_merged["rx_left"].corr(df_merged["rx_right"])

# Graph
plt.figure(figsize=(10, 6))
plt.plot(df_merged.index, df_merged["rx_left"], label="Left Eye - rx", color="blue")
plt.plot(df_merged.index, df_merged["rx_right"], label="Right Eye - rx", color="red")
plt.xlabel("Time (volume number)")
plt.ylabel("Rotation X (degrees)")
plt.title(f"Comparison of X Rotation (Correlation: {Correlation:.2f})")
plt.legend()
plt.grid(True)
plt.show()
