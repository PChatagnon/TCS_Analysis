import numpy as np
import math 
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define your custom function
def cross_section(input_var, C_0, m_A, m_C):
    mt, Epho = input_var
    
    t = -1 * mt
    
    A_0 = 0.414
    
    Mv = 3.1 # Mass of the JPsi
    Mn = 0.938 # Mass of the proton
    
    Conv_factor = 0.3894*(10^6) # GeV-2 to nb
    
    alpha_EM = 1/137.
    alpha_S = 0.3
    e_q = (2/3)
    Phi_V_2 = 1.0952/(4. * math.pi)
    
    W = (Mn**2 + 2*Mn*Epho)
    tau = Mv**2 / ( W - Mn**2 )
    xi = tau /(2 - tau)
    #xi = (Mv**2)/(4*Mn*Epho - Mv**2) # Not sure of this one
    
    A_t = A_0 / ( 1- (t/(m_A**2)) )**2
    C_t = C_0 / ( 1- (t/(m_C**2)) )**2
    
    H2 = A_t + ((2*xi)**2) * C_t
    E2 = -1.0 * ((2*xi)**2) * C_t
    G_2 = (1.0/xi**4) * ( (1-(t/(4*Mn**2)))*E2**2 - 2*E2*(H2+E2) +(1-xi**2)*(H2+E2)**2 )
    
    sigma = Conv_factor*( (alpha_EM * e_q**2 ) / ( 4*(W*W - Mn**2)**2 ) ) * ( (16.0*math.pi*alpha_S)**2 / (3*Mv**3) ) * Phi_V_2 * G_2
    
    #sigma = C_t*C_t
    
    return sigma

with open("/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/CS_Extraction_combine_t_05_rad_new_8.20_9.28_Latex_Table.txt", "r") as file:
    lines = file.readlines()
    data_1 = np.array([[float(value) for value in line.split("&")] for line in lines])
    
# Assign constant y value for the first file
constant_y_value_file1 = 9.00
y_data_file1 = np.full_like(data_1[1], constant_y_value_file1)
    
with open("/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/CS_Extraction_combine_t_05_rad_new_2_9.28_10.36_Latex_Table.txt", "r") as file:
    lines = file.readlines()
    data_2 = np.array([[float(value) for value in line.split("&")] for line in lines])

# Assign constant y value for the second file
constant_y_value_file2 = 9.86
y_data_file2 = np.full_like(data_2[1], constant_y_value_file2)

# Concatenate data from both files
x_data = np.concatenate((data_1[0], data_2[0]))
y_data = np.concatenate((y_data_file1, y_data_file2))
z_data = np.concatenate((data_1[2], data_2[2]))

sigma_z = np.concatenate((data_1[3], data_2[3]))

print(x_data)

print(z_data)

# Perform the fit
#initial_guess = [1.0, 1.5, 1.0]  # Initial guess for the parameters
initial_guess = [-10000.2, 1.0, 1.25]  # Initial guess for the parameters
params, covariance = curve_fit(cross_section, (x_data, y_data), z_data, p0=initial_guess, sigma=sigma_z, absolute_sigma=True)

# Extract the errors from the covariance matrix
errors = np.sqrt(np.diag(covariance))

# Print the fitted parameters and their errors
print("Fitted parameters:")
print("C_0 =", params[0], "±", errors[0])
print("m_A =", params[1], "±", errors[1])
print("m_C =", params[2], "±", errors[2])

# Generate meshgrid for 3D plot
x_min, x_max = np.min(x_data), np.max(x_data)
y_min, y_max = np.min(y_data), np.max(y_data)
x_range = np.linspace(x_min, x_max, 100)
y_range = np.linspace(y_min, y_max, 100)
X, Y = np.meshgrid(x_range, y_range)
Z = cross_section((X, Y), *params)
#Z = cross_section((X, Y), -10000.2, 10002.7, 1.25)


print(params)

# Create 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_data, y_data, z_data, color='blue', label='Data')
ax.errorbar(x_data, y_data, z_data, zerr=sigma_z, fmt='o', color='blue', label='Data')
ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.5, label='Fit')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
#plt.legend()
plt.savefig("3d_plot.png")
plt.show()