# ------------------------------------------
# T-test graphs for grouped subjects 
# by stroke and left/right eye
# ------------------------------------------

import scipy.stats as stats

# Mean values for every subject
df_subject_mean = df_avg.groupby("subject")[["tx", "ty", "tz", "rx", "ry", "rz"]].mean().reset_index()

# control/pathological group
df_subject_mean["group"] = df_subject_mean["subject"].apply(
    lambda x: "control" if "CON" in x.upper() else "pathological"
)

# group deviation
control_group = df_subject_mean[df_subject_mean["group"] == "control"]
pathological_group = df_subject_mean[df_subject_mean["group"] == "pathological"]

# T-test for every numerical column
t_results = {}
for param in ["tx", "ty", "tz", "rx", "ry", "rz"]:
    t_stat, p_value = stats.ttest_ind(control_group[param], pathological_group[param], equal_var=False)
    t_results[param] = (t_stat, p_value)

# result graph for t-test
params = list(t_results.keys())
t_values = [t_results[p][0] for p in params]
p_values = [t_results[p][1] for p in params]

plt.figure(figsize=(10, 5))
plt.bar(params, t_values, color=["blue" if p > 0.05 else "red" for p in p_values])
plt.axhline(y=0, color="black", linestyle="--")
plt.xlabel("Motion Parameter")
plt.ylabel("t-Statistic")
plt.title("T-Student Test for Motion Parameters (Control vs Pathological)")
plt.grid(axis="y")
plt.show()

# result presentation
print("T-test rotation parameters:")
for param, (t_stat, p_value) in t_results.items():
    print(f"{param}: t = {t_stat:.2f}, p = {p_value:.4f} {'(Statistically Significant)' if p_value < 0.05 else ''}")

# Mean values for rotation for every patient
df_subject_mean = df_avg.groupby("subject")[["rx", "ry", "rz"]].mean().reset_index()

# Euclidean value for patients
df_subject_mean["euclidean_distance"] = np.sqrt(
    df_subject_mean["rx"]**2 + df_subject_mean["ry"]**2 + df_subject_mean["rz"]**2
)

# control/pathological columns
df_subject_mean["group"] = df_subject_mean["subject"].apply(
    lambda x: "control" if "CON" in x.upper() else "pathological"
)

# group separation
control_group = df_subject_mean[df_subject_mean["group"] == "control"]["euclidean_distance"]
pathological_group = df_subject_mean[df_subject_mean["group"] == "pathological"]["euclidean_distance"]

# T-test
t_stat, p_value = stats.ttest_ind(control_group, pathological_group, equal_var=False)

# Euclidean Distance for every group
plt.figure(figsize=(8, 6))
plt.boxplot([control_group, pathological_group], labels=["Control", "Pathological"], patch_artist=True,
            boxprops=dict(facecolor="lightblue"), medianprops=dict(color="red"))
plt.ylabel("Euclidean Distance")
plt.title("Comparison of Euclidean Distance (Control vs Pathological)")
plt.grid(axis="y")

plt.text(1.5, max(df_subject_mean["euclidean_distance"]) * 0.9, f"t = {t_stat:.2f}, p = {p_value:.4f}", ha="center")

plt.show()

# t-Test result
print(f"T-Student test:\n t = {t_stat:.2f}, p = {p_value:.4f}")
if p_value < 0.05:
    print("There is a significant difference.")
else:
    print("There is NO significant difference.")
