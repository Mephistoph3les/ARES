# execute with      python3 scattering_calibration_plot_burj.py 1       for      Dataset1
# or                python3 scattering_calibration_plot_burj.py 2       for      Dataset2

import pandas as pd 
import sys
import matplotlib.pyplot as plt 
import math
import numpy as np



dataset = int(sys.argv[1])
print('the chosen dataset is: Dataset',dataset)

def from_mean_to_AAD(file_name, reference_file_name):
    
    #reading in the data
    data = pd.read_csv(file_name, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])
    #reading in reference beam data to subtract from actual data for beam widening
    reference = pd.read_csv(reference_file_name, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])

    ref_mean = reference.loc[:, 'Width'].mean()
    
    d = 75*10**(-3) # distance burj <-> detector in m
    data.loc[:, ["Width"]] = data.loc[:, ["Width"]] - ref_mean
    data_mean = data.loc[:, 'Width'].mean()
    data_length = data.shape[0]
    sum = float(0)
    for i in range (data_length):
        sum += abs(data.iloc[i, 2] - data_mean)
    AAD = (1/data_length) * sum     #right now AAD still in Pixels
    AAD = (AAD/2)*(55*10**(-6))
    AAD = np.arctan(AAD/d)*10**(+3)
    return AAD**2

def highland_AAD(material_budget):
    
    beam_energy_MeV = 154.5 #Beam energy in MeV
    electron_mass = 0.510998950 # in MeV
    momentum = ((beam_energy_MeV)**2 - (electron_mass)**2)**(0.5)
    Theta = (13.6/(momentum))*((material_budget)**(0.5))*(1+0.038*np.log(material_budget))
    Theta = Theta*10**(+3)
    return Theta**2

if __name__ == '__main__':
    if dataset == 1:
        files_list=(['Dataset1/output_measurement_191.dat', 'Dataset1/output_measurement_192.dat', 'Dataset1/output_measurement_193.dat', 'Dataset1/output_measurement_194.dat', 'Dataset1/output_measurement_195.dat', 'Dataset1/output_measurement_196.dat', 'Dataset1/output_measurement_197.dat', 'Dataset1/output_measurement_198.dat'])
        reference_file_name = "Dataset1/output_measurement_190.dat"
    elif dataset == 2:
        files_list=(['Dataset2/output_measurement_223.dat', 'Dataset2/output_measurement_224.dat', 'Dataset2/output_measurement_225.dat', 'Dataset2/output_measurement_226.dat', 'Dataset2/output_measurement_227.dat', 'Dataset2/output_measurement_228.dat', 'Dataset2/output_measurement_229.dat', 'Dataset2/output_measurement_230.dat'])
        reference_file_name = "Dataset2/output_measurement_222.dat"

    n = len(files_list)
    AAD_list=([])

    for i in range (n):
        AAD_list.append(from_mean_to_AAD(files_list[i], reference_file_name))
    print("and the list of AAD's: ")
    print(AAD_list)


    radiation_length_PEEK = 319.6   #in mm
    error_in_radiation_length = 6.4 # ~2% error
    list_of_material_thicknesses = ([1, 6, 11, 16, 21, 26, 31, 36]) #in mm
    error_in_material_thickness = 0.1
    
    m = len(list_of_material_thicknesses)

    material_budget_list = ([])
    highland_prediction_list = ([])
    material_budget_error = ([])
    for i in range (m):
        material_budget_list.append(list_of_material_thicknesses[i]/radiation_length_PEEK)
        material_budget_error = ((error_in_material_thickness/radiation_length_PEEK)**2+(error_in_radiation_length*list_of_material_thicknesses[i]/(radiation_length_PEEK)**2)**2)**(0.5)
        highland_prediction_list.append(highland_AAD(material_budget_list[i]))
    print("Thats the list of material budgets")
    print(material_budget_list)
    print("Thats the list of highlands predictions")
    print(highland_prediction_list)



    fig, ax = plt.subplots(figsize=(10,10), layout='constrained')
    ax.scatter(AAD_list, material_budget_list, label = 'Burj (PEEK)') #Plotting data onto the axes
    ax.errorbar(AAD_list, material_budget_list, yerr = material_budget_error, fmt="o")
    ax.plot(highland_prediction_list, material_budget_list, label = 'Highland')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$AAD^2$ [$mrad^2$]', style='normal')
    ax.set_ylabel('Material Budget')
    ax.set_title('Calibration Plot')
    ax.legend()
    plt.show()