import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.patches as patches

# Data Hall C
Eg_values_hall_c = np.array([9.175, 9.375, 9.475, 9.625, 9.775, 9.925, 10.075, 10.225, 10.375, 10.525])
eg_errors_hall_c = np.array([0.025, 0.009, 0.017, 0.018, 0.043, 0.038, 0.024, 0.027, 0.033, 0.038])
mA_values_hall_c = np.array([1.69, 2.30, 1.59, 1.52, 1.18, 1.26, 1.43, 1.31, 1.33, 1.32])
mA_errors_hall_c = np.array([0.45, 0.47, 0.22, 0.17, 0.16, 0.17, 0.16, 0.14, 0.19, 0.20])
mass_radius_values_hall_c = np.array([0.404, 0.297, 0.430, 0.450, 0.579, 0.544, 0.477, 0.520, 0.515, 0.517])
mass_radius_errors_hall_c = np.array([0.107, 0.061, 0.059, 0.052, 0.079, 0.075, 0.053, 0.057, 0.074, 0.079])
color_hall_c = '#154734'
color_rec_hall_c='#79863c'

# Data CLAS12
Eg_values_CLAS12 = np.array([9.00,9.86])
eg_errors_up_CLAS12 = np.array([0.28,0.50])
eg_errors_down_CLAS12 = np.array([0.80, 0.58])
mA_values_CLAS12 = np.array([1.43, 1.27])
mA_errors_CLAS12 = np.array([0.23, 0.09])
mass_radius_values_CLAS12 = np.array([0.48, 0.54])
mass_radius_errors_CLAS12 = np.array([0.08, 0.04])
color_CLAS12 = '#00b4d8'#'ffb81c'

# Data GlueX
Eg_values_GlueX = np.array([8.93, 9.86, 10.82])
eg_errors_up_GlueX = np.array([0.35, 0.5, 0.62])
eg_errors_down_GlueX = np.array([0.73, 0.58, 0.46])
mA_values_GlueX = np.array([1.105, 1.472, 1.313])
mA_errors_GlueX = np.array([0.168, 0.075, 0.049])
mass_radius_values_GlueX = np.array([0.619, 0.464, 0.521])
mass_radius_errors_GlueX = np.array([0.094, 0.024, 0.020])
color_glueX = '#ce0037'

# Calculate mean and standard deviation for the first four data points of the blue dataset
mean_hall_c_1 = np.mean(mass_radius_values_hall_c[:4])
std_dev_hall_c_1 = np.std(mass_radius_values_hall_c[:4])

mean_hall_c_2 = np.mean(mass_radius_values_hall_c[4:])
std_dev_hall_c_2 = np.std(mass_radius_values_hall_c[4:])

# Create a square grid with two columns for the two plots
fig = plt.figure(figsize=(8, 8))
gs = GridSpec(1, 2, width_ratios=[1, 1])

# Plot mA as a function of Eg with rotated graph
ax0 = plt.subplot(gs[0])
ax0.errorbar(mA_values_hall_c, Eg_values_hall_c, xerr=mA_errors_hall_c,  marker='o', linestyle='', color=color_hall_c, capsize=5, label='Hall C-007 (2022)')
ax0.errorbar(mA_values_GlueX, Eg_values_GlueX, xerr=mA_errors_GlueX, yerr=[eg_errors_down_GlueX, eg_errors_up_GlueX], marker='o', linestyle='', color=color_glueX, capsize=5, label='GlueX (2023)')
ax0.errorbar(mA_values_CLAS12, Eg_values_CLAS12, xerr=mA_errors_CLAS12, yerr=[eg_errors_down_CLAS12, eg_errors_up_CLAS12], marker='o', linestyle='', color=color_CLAS12, capsize=5, label='CLAS12')
ax0.set_ylabel('Eg (GeV)')
ax0.set_xlabel('mA (GeV)')

# Plot mass radius as a function of Eg with rotated graph
ax1 = plt.subplot(gs[1], sharey=ax0)
ax1.errorbar(mass_radius_values_hall_c, Eg_values_hall_c, xerr=mass_radius_errors_hall_c, marker='o', linestyle='', color=color_hall_c, capsize=5, label='Hall C-007')
ax1.errorbar(mass_radius_values_CLAS12, Eg_values_CLAS12, xerr=mass_radius_errors_CLAS12, yerr=[eg_errors_down_CLAS12, eg_errors_up_CLAS12], marker='o', linestyle='', color=color_CLAS12, capsize=5, label='CLAS12')
ax1.errorbar(mass_radius_values_GlueX, Eg_values_GlueX, xerr=mass_radius_errors_GlueX, yerr=[eg_errors_down_GlueX, eg_errors_up_GlueX], marker='o', linestyle='', color=color_glueX, capsize=5, label='GlueX')
ax1.set_xlabel('Mass Radius (fm)')
ax1.tick_params(axis='y', labelleft=False)  # Remove y-axis label on the right plot

# Add rectangle for the mean and standard deviation of the first four data points of the blue dataset
#rect = patches.Rectangle((min(Eg_values[:4]), mean_blue - std_dev_blue / 2), max(Eg_values[:4]) - min(Eg_values[:4]), std_dev_blue, linewidth=1, edgecolor='none', facecolor='blue')
rect_1 = patches.Rectangle((mean_hall_c_1, Eg_values_hall_c[0]),std_dev_hall_c_1, Eg_values_hall_c[3]-Eg_values_hall_c[0], linewidth=1, edgecolor='none', facecolor=color_rec_hall_c, alpha=0.5)
ax1.add_patch(rect_1)

rect_2 = patches.Rectangle((mean_hall_c_2, Eg_values_hall_c[4]),std_dev_hall_c_2, Eg_values_hall_c[9]-Eg_values_hall_c[4], linewidth=1, edgecolor='none', facecolor=color_rec_hall_c, alpha=0.5)
ax1.add_patch(rect_2)

# Add legend to the right plot without box
legend = ax0.legend(loc='upper right', ncol=1, fancybox=True, shadow=False, fontsize=13)
frame = legend.get_frame()
frame.set_linewidth(0)  # Remove the box around the legend

# Remove top and right border lines
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)

# Remove top and right border lines
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Force the y-axis to finish at 10.6
ax0.set_ylim(ax0.get_ylim()[0], 11.5)

# Force the x-axis to finish at 0.7
ax1.set_xlim(ax1.get_xlim()[0], 0.8)

ax0.set_xlabel('m$_s$ (GeV)')  # Replace 'mA' with 'ms'
ax0.set_ylabel('E$γ$ (GeV)')  # Replace 'Eg' with 'Eγ'


plt.text(0.5, 0.5, 'Preliminary', fontsize=100, rotation=45, ha='center', va='center', alpha=0.2, transform=fig.transFigure)

# Adjust layout to prevent overlap and remove space between graphs
plt.subplots_adjust(wspace=0)

# Save the figure as a PNG file
plt.savefig('M_S_and_r_m.png')

# Show the plots
plt.show()
