# execute with      python3 scattering_calibration_plot_burj.py 1       for      Dataset1
# or                python3 scattering_calibration_plot_burj.py 2       for      Dataset2

import pandas as pd 
import sys
import matplotlib.pyplot as plt 
import math
import numpy as np



dataset = int(sys.argv[1])
print('the chosen dataset is: Dataset',dataset)
print("")

def get_mean(file_name, reference_file_name):
    
    #reading in the data
    data = pd.read_csv(file_name, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])
    #reading in reference beam data to subtract from actual data for beam widening
    reference = pd.read_csv(reference_file_name, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])

    ref_mean = reference.loc[:, 'Width'].mean()
    
    d = 75*10**(-3) # distance burj <-> detector in m
    data.loc[:, ["Width"]] = np.sqrt((data.loc[:, ["Width"]])**2 - (ref_mean)**2)      #background deduction
    data_mean = data.loc[:, 'Width'].mean() *(55*10**(-6))
    mean = np.arctan(data_mean/d)*10**(+3)
    
    return mean**2

def highland(material_budget):
    
    beam_energy_MeV = 154.5                                                                     #beam energy in MeV
    electron_mass = 0.510998950                                                                 #electron mass in MeV
    momentum = ((beam_energy_MeV)**2 - (electron_mass)**2)**(0.5)                               #relativistic momentum of electron beam
    lorentz_factor = np.sqrt(1+(momentum/electron_mass))    
    beta = np.sqrt(1-(1/lorentz_factor)**2)
    Theta = (13.6/(beta*momentum))*((material_budget)**(0.5))*(1+0.038*np.log(material_budget))  #Highland formula with c=1
    Theta = Theta*10**(+3)                                                                      #from rad to mrad
    return Theta**2, beta, momentum

if __name__ == '__main__':
    if dataset == 1:
        files_list=(['Dataset1/output_measurement_191.dat', 'Dataset1/output_measurement_192.dat', 'Dataset1/output_measurement_193.dat', 'Dataset1/output_measurement_194.dat', 'Dataset1/output_measurement_195.dat', 'Dataset1/output_measurement_196.dat', 'Dataset1/output_measurement_197.dat', 'Dataset1/output_measurement_198.dat'])
        reference_file_name = "Dataset1/output_measurement_190.dat"
    elif dataset == 2:
        files_list=(['Dataset2/output_measurement_223.dat', 'Dataset2/output_measurement_224.dat', 'Dataset2/output_measurement_225.dat', 'Dataset2/output_measurement_226.dat', 'Dataset2/output_measurement_227.dat', 'Dataset2/output_measurement_228.dat', 'Dataset2/output_measurement_229.dat', 'Dataset2/output_measurement_230.dat'])
        reference_file_name = "Dataset2/output_measurement_222.dat"

    n = len(files_list)
    means_list=([])

    for i in range (n):
        means_list.append(get_mean(files_list[i], reference_file_name))
    print("Thats the list of means of the data width: ")
    print(means_list)


    radiation_length_PEEK = 319.6                                       #in mm
    error_in_radiation_length = 6.4                                     # ~2% error
    list_of_material_thicknesses = ([1, 6, 11, 16, 21, 26, 31, 36])     #in mm
    error_in_material_thickness = 0.1                                   #error in mm
    
    m = len(list_of_material_thicknesses)

    material_budget_list = ([])
    highland_prediction_list = ([])
    material_budget_error_list = ([])
    for i in range (m):
        material_budget_list.append(list_of_material_thicknesses[i]/radiation_length_PEEK)
        material_budget_error = np.sqrt(((error_in_material_thickness/radiation_length_PEEK)**2+(error_in_radiation_length*list_of_material_thicknesses[i]/(radiation_length_PEEK)**2)**2))
        material_budget_error_list.append(material_budget_error)
        highland_prediction_list.append(highland(material_budget_list[i])[0])
    print("Thats the list of material budgets: ")
    print(material_budget_list)
    print("Thats the list of material budget errors: ")
    print(material_budget_error_list)
    print("Thats the list of highlands predictions: ")
    print(highland_prediction_list)
    print("")
    print("beta is: ", highland(material_budget_list[0])[1])
    print("electron momentum is: ", highland(material_budget_list[0])[2], "MeV")



    fig, ax = plt.subplots(figsize=(10,10), layout='constrained')
    ax.scatter(means_list, material_budget_list, label = 'Burj (PEEK)') #Plotting data onto the axes
    ax.errorbar(means_list, material_budget_list, yerr = material_budget_error_list, fmt="o")
    ax.plot(highland_prediction_list, material_budget_list, label = 'Highland')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$mean^2$ [$mrad^2$]', style='normal')
    ax.set_ylabel('Material Budget')
    ax.set_title('Calibration Plot')
    ax.legend()
    plt.show()