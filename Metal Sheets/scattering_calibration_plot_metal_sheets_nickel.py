# execute with      python3 scattering_calibration_plot_metal_sheets.py

import pandas as pd 
import sys
import matplotlib.pyplot as plt 
import math
import numpy as np

def time_cut_of_data(file_name, logbook_name):
    data = pd.read_csv(file_name, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])
    logbook = pd.read_csv(logbook_name, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])
    

def get_mean(file_name, reference_file_name):
    
    #reading in the data
    data = pd.read_csv(file_name, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])
    #reading in reference beam data to subtract from actual data for beam widening
    reference = pd.read_csv(reference_file_name, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])

    ref_mean = reference.loc[:, 'Width'].mean()
    
    d = 75*10**(-3) # distance burj <-> detector in m
    data.loc[:, ["Width"]] = (np.sqrt((data.loc[:, ["Width"]])**2 - (ref_mean)**2)) *(55*10**(-6))     #background deduction
    data.loc[:, ["Width"]] = (np.arctan(data.loc[:, ["Width"]]/d)*10**(+3))
    data_mean = data.loc[:, 'Width'].mean() 
    data_std = data.loc[:, 'Width'].std() 
    N = data.shape[0]

    return data_mean**2, data_std, N

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
    # files_list=(['Dataset1/output_measurement_191.dat', 'Dataset1/output_measurement_192.dat', 'Dataset1/output_measurement_193.dat', 'Dataset1/output_measurement_194.dat', 'Dataset1/output_measurement_195.dat', 'Dataset1/output_measurement_196.dat', 'Dataset1/output_measurement_197.dat', 'Dataset1/output_measurement_198.dat'])
    # reference_file_name = "Dataset1/output_measurement_190.dat"         
    n = len(files_list)
    means_list=([])
    std_list=([])
    N_list=([])
    for i in range (n):
        means_list.append(get_mean(files_list[i], reference_file_name)[0])
        std_list.append(get_mean(files_list[i], reference_file_name)[1])
        N_list.append(get_mean(files_list[i], reference_file_name)[2])
    print("Thats the list of means of the data width: ")
    print(means_list)
    print("The list of standard deviations of those: ")
    print(std_list)
    print("The list of sample N's: ")
    print(N_list)


    radiation_length_nickel= 14.24                                      #in mm
    error_in_radiation_length = 0.28                                    # ~2% error
    list_of_material_thicknesses_overnight = ([0.025, 0.075, 0.05, 0.15, 0.1, 0.35, 0.25, 0.75, 0.5, 1.5, 1, 3, 2, 3, 1, 1.5, 0.5, 0.75, 0.25, 0.35, 0.1, 0.15, 0.05, 0.075, 0.025])     #in mm
#    list_of_material_thicknesses_shortscan = ([2, 3, 1, 1.5, 0.5, 0.75, 0.25, 0.35, 0.1, 0.15, 0.05, 0.075, 0.025, 0])     #in mm
    error_in_material_thickness = 0.0005                                   #error in mm
    
    m = len(list_of_material_thicknesses_overnight)

    material_budget_list = ([])
    highland_prediction_list = ([])
    material_budget_error_list = ([])
    for i in range (m):
        material_budget_list.append(list_of_material_thicknesses_overnight[i]/radiation_length_nickel)
        material_budget_error = np.sqrt(((error_in_material_thickness/radiation_length_nickel)**2+(error_in_radiation_length*list_of_material_thicknesses_overnight[i]/(radiation_length_nickel)**2)**2))
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
    ax.scatter(means_list, material_budget_list, label = 'Nickel') #Plotting data onto the axes
    ax.errorbar(means_list, material_budget_list, yerr = material_budget_error_list , xerr = std_list/(np.sqrt(N_list)) , fmt="o")
    ax.plot(highland_prediction_list, material_budget_list, label = 'Highland')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$mean^2$ [$mrad^2$]', style='normal')
    ax.set_ylabel('Material Budget')
    ax.set_title('Calibration Plot')
    ax.legend()
    plt.show()