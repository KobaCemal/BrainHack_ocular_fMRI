# --------------------------------------------
# Graph for correlation, eclidean distance
# and mean rotation for every patient
# --------------------------------------------

# Libreris
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Data file path
df = pd.read_csv("file_path", sep=",")

# Aggregation - mean values for patient and time step
# Mean for numerical columns only
df_avg = df.groupby(["subject", "time"])[["tx", "ty", "tz", "rx", "ry", "rz"]].mean().reset_index()

# Group column
df_avg["group"] = df_avg["subject"].apply(
    lambda x: "control" if "CON" in x.upper() else ("pathological" if "PAT" in x.upper() else "unknown")
)

print("Example mean: \n")
print(df_avg.head())

print("""
--------------------------
1. X axis mean values graph (rx) in time for both subject groups
--------------------------
""") 

# Mean value group only for numerical columns
group_time_mean = df_avg.groupby(["group", "time"])[["rx"]].mean().reset_index()

plt.figure(figsize=(10, 6))
for grp in group_time_mean["group"].unique():
    grp_data = group_time_mean[group_time_mean["group"] == grp]
    plt.plot(grp_data["time"], grp_data["rx"], marker='o', label=f"{grp} (rx)")
plt.xlabel("Time (volume number)")
plt.ylabel("Average rx (degrees)")
plt.title("Average X Rotation over Time by Group")
plt.legend()
plt.grid(True)
plt.show()

print("""
--------------------------
2. Euclidean distance value for mean values in rotation
--------------------------
""")

# Mean values for subject group and time step
group_params = df_avg.groupby(["group", "time"])[["rx", "ry", "rz"]].mean().reset_index()

# Pivot
pivot_rx = group_params.pivot(index="time", columns="group", values="rx")
pivot_ry = group_params.pivot(index="time", columns="group", values="ry")
pivot_rz = group_params.pivot(index="time", columns="group", values="rz")

# [rx, ry, rz] mean values in euclidean distance 
# for 'control' and 'pathological' group
# Euclidean distance function
def euclidean_distance(df):
    return np.sqrt(df["rx"]**2 + df["ry"]**2 + df["rz"]**2)

df_avg["euclidean_distance"] = euclidean_distance(df_avg)
group_euclidean = df_avg.groupby(["group", "time"])["euclidean_distance"].mean().reset_index()

plt.figure(figsize=(10, 6))
for grp in group_euclidean["group"].unique():
    grp_data = group_euclidean[group_euclidean["group"] == grp]
    plt.plot(grp_data["time"], grp_data["euclidean_distance"], marker='o', label=f"{grp}")
plt.xlabel("Time (volume number)")
plt.ylabel("Euclidean Distance")
plt.title("Average Euclidean Distance of Rotation for Each Group")
plt.legend()
plt.grid(True)
plt.show()

print("""
--------------------------
3. correlation matrix for mean rotation values
--------------------------
""")

# correlation matrix for numerical columns in control group only
df_control = df_avg[df_avg["group"] == "control"]
corr_matrix_control = df_control[["tx", "ty", "tz", "rx", "ry", "rz"]].corr()

plt.figure(figsize=(8, 6))
sns.heatmap(corr_matrix_control, annot=True, cmap="coolwarm", fmt=".2f")
plt.title("Correlation Matrix of Motion Parameters (Averaged over Eyes)")
plt.show()
