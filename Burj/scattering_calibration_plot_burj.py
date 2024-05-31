# execute with      python3 scattering_calibration_plot_burj.py 1       for      Dataset1
# or                python3 scattering_calibration_plot_burj.py 2       for      Dataset2

import pandas as pd 
import sys
import matplotlib.pyplot as plt 
import math
import numpy as np


def get_mean(file_name, reference_file_name):
    
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

    return data_mean**2, data_std**2, N

def get_mean_both(file_name, file_name_2, reference_file_name, reference_file_name_2):
    
    #reading in the data
    data = pd.read_csv(file_name, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])
    #reading in reference beam data to subtract from actual data for beam widening
    reference = pd.read_csv(reference_file_name, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])
    
    data_2 = pd.read_csv(file_name_2, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])
    #reading in reference beam data to subtract from actual data for beam widening
    reference_2 = pd.read_csv(reference_file_name_2, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])

    ref_mean = reference.loc[:, 'Width'].mean()
    ref_mean_2 = reference_2.loc[:, 'Width'].mean()
    
    d = 75*10**(-3) # distance burj <-> detector in m
    data.loc[:, ["Width"]] = (np.sqrt((data.loc[:, ["Width"]])**2 - (ref_mean)**2)) *(55*10**(-6))     #background deduction
    data.loc[:, ["Width"]] = (np.arctan(data.loc[:, ["Width"]]/d)*10**(+3))
    
    data_2.loc[:, ["Width"]] = (np.sqrt((data_2.loc[:, ["Width"]])**2 - (ref_mean_2)**2)) *(55*10**(-6))     #background deduction
    data_2.loc[:, ["Width"]] = (np.arctan(data_2.loc[:, ["Width"]]/d)*10**(+3))
    
    frames = [data, data_2]
    result = pd.concat(frames)
    data_mean = result.loc[:, 'Width'].mean() 
    data_std = result.loc[:, 'Width'].std() 
    N = data.shape[0]

    return data_mean**2, data_std**2, N

def highland(material_budget):
    
    beam_energy_MeV = 154.5                                                                     #beam energy in MeV
    electron_mass = 0.510998950                                                                 #electron mass in MeV
    momentum = ((beam_energy_MeV)**2 - (electron_mass)**2)**(0.5)                               #relativistic momentum of electron beam
    lorentz_factor = np.sqrt(1+(momentum/electron_mass))    
    beta = np.sqrt(1-(1/lorentz_factor)**2)
    Theta = (13.6/(beta*momentum))*((material_budget)**(0.5))*(1+0.038*np.log(material_budget))  #Highland formula with c=1
    Theta = Theta*10**(+3)                                                                      #from rad to mrad
    return Theta**2, beta, momentum

def ladderplot_both(file_names_1, file_names_2, reference_file_name_1, reference_file_name_2):
    list_of_material_thicknesses = ([0, 1, 6, 11, 16, 21, 26, 31, 36])
    material = "Peek"
    widths=[]
    materials = []
    ladder_filename = f"ladderplot_{material.lower()}_merged"

    data_list_1 =([0, 0, 0, 0, 0, 0, 0, 0])
    data_list_2 =([0, 0, 0, 0, 0, 0, 0, 0])
    for i in range(len(data_list_1)):
        data_list_1[i] = pd.read_csv(file_names_1[i], sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])
        data_list_2[i] = pd.read_csv(file_names_2[i], sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])

    reference_1 = pd.read_csv(reference_file_name_1, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])
    reference_2 = pd.read_csv(reference_file_name_2, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp']) 

    for i in range(len(list_of_material_thicknesses)):
        if i == 0:
            data_1 = reference_1.loc[:, "Width"]
            data_2 = reference_2.loc[:, "Width"]
        else:
            data_1 = data_list_1[i-1].loc[:, "Width"]
            data_2 = data_list_2[i-1].loc[:, "Width"]

        frames = [data_1, data_2]
        result = pd.concat(frames)
        shape = result.shape[0]
        material_thicknesses=[]
        for _ in range(shape):
            material_thicknesses.append(list_of_material_thicknesses[i])
        arr = np.array(material_thicknesses)
        materials.append(arr)
        result = result.to_numpy()
        widths.append(result)


    fig, ax = plt.subplots(figsize=(10,10), layout='constrained')
    for i in range(len(list_of_material_thicknesses)):
        ax.scatter(widths[i], materials[i], label = ('Material Thickness (in mm): ', list_of_material_thicknesses[i]))
    ax.set_xlabel('Width (in pixels)', style='normal')
    ax.set_ylabel('Material Thickness (in mm)')
    ax.set_title(f'Ladderplot {material} merged')
    ax.legend()
    ax.grid()
    ax.set_axisbelow(True)
    plt.xlim([0, 50])
    plt.savefig(ladder_filename)
    print("Ladderplot saved as " + ladder_filename)
    
def ladderplot(file_names, reference_file_name, dataset):
    list_of_material_thicknesses = ([0, 1, 6, 11, 16, 21, 26, 31, 36])
    material = "Peek"
    widths=[]
    materials = []
    
    ladder_filename = f"ladderplot_{material.lower()}_dataset_{str(dataset)}"
        
    data_list =([0, 0, 0, 0, 0, 0, 0, 0])
    
    for i in range(len(data_list)):
        data_list[i] = pd.read_csv(file_names[i], sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])

    reference = pd.read_csv(reference_file_name, sep='\s+', names= ['Event', 'Macropulse', 'Width', 'Timestamp'])    

        
    for i in range(len(list_of_material_thicknesses)):
        if i == 0:
            data = reference.loc[:, "Width"]
        else:
            data = data_list[i-1].loc[:, "Width"]
        shape = data.shape[0]
        material_thicknesses=[]
        for _ in range(shape):
            material_thicknesses.append(list_of_material_thicknesses[i])
        arr = np.array(material_thicknesses)
        materials.append(arr)
        result = data.to_numpy()
        widths.append(result)

    fig, ax = plt.subplots(figsize=(10,10), layout='constrained')
    for i in range(len(list_of_material_thicknesses)):
        ax.scatter(widths[i], materials[i], label = ('Material Thickness (in mm): ', list_of_material_thicknesses[i]))
    ax.set_xlabel('Width (in pixels)', style='normal')
    ax.set_ylabel('Material Thickness (in mm)')
    ax.set_title(f'Ladderplot {material} of dataset {str(dataset)}')
    ax.legend()
    ax.grid()
    ax.set_axisbelow(True)
    plt.xlim([0, 50])
    plt.savefig(ladder_filename)
    print("Ladderplot saved as " + ladder_filename)
    
    
if __name__ == '__main__':

    radiation_length_PEEK = 319.6                                       #in mm
    error_in_radiation_length = 0.5                                     # ~2% error
    list_of_material_thicknesses = ([1, 6, 11, 16, 21, 26, 31, 36])     #in mm
    error_in_material_thickness = 0.05                                   #error in mm
    
    m = len(list_of_material_thicknesses)

    material_budget_list = ([])
    highland_prediction_list = ([])
    material_budget_error_list = ([])
    for i in range (m):
        material_budget_list.append(list_of_material_thicknesses[i]/radiation_length_PEEK)
        material_budget_error = np.sqrt(((error_in_material_thickness/radiation_length_PEEK)**2+(error_in_radiation_length*list_of_material_thicknesses[i]/(radiation_length_PEEK)**2)**2))
        material_budget_error_list.append(material_budget_error)
        highland_prediction_list.append(highland(material_budget_list[i])[0])


    files_list_1=(['Dataset1/output_measurement_191.dat', 'Dataset1/output_measurement_192.dat', 'Dataset1/output_measurement_193.dat', 'Dataset1/output_measurement_194.dat', 'Dataset1/output_measurement_195.dat', 'Dataset1/output_measurement_196.dat', 'Dataset1/output_measurement_197.dat', 'Dataset1/output_measurement_198.dat'])
    reference_file_name_1 = "Dataset1/output_measurement_190.dat"
    files_list_2=(['Dataset2/output_measurement_223.dat', 'Dataset2/output_measurement_224.dat', 'Dataset2/output_measurement_225.dat', 'Dataset2/output_measurement_226.dat', 'Dataset2/output_measurement_227.dat', 'Dataset2/output_measurement_228.dat', 'Dataset2/output_measurement_229.dat', 'Dataset2/output_measurement_230.dat'])
    reference_file_name_2 = "Dataset2/output_measurement_222.dat"

    n = len(files_list_1)
    means_list=([])
    std_list=([])
    N_list=([])
    

    ask = input("Do you want to analyse both datasets (y) or just one (n)? ")
    if ask.lower() in ["y","yes"]:
        dataset = 0
        ladderplot_both(files_list_1, files_list_2, reference_file_name_1, reference_file_name_2)
        for i in range (n):
            help_list = get_mean_both(files_list_1[i], files_list_2[i], reference_file_name_1, reference_file_name_2)
            means_list.append(help_list[0])
            std_list.append(help_list[1])
            N_list.append(help_list[2])
        xerror = std_list/(np.sqrt(N_list))
        
        peek = {'Material Budget Peek': material_budget_list, ' Material Budget error': material_budget_error_list, ' mean-squared deviation angle from reference beam': means_list, ' xerror': xerror}
        plot_data_peek = pd.DataFrame(data=peek)
        plot_data_peek.to_csv('Calibration_data/calibration_data_peek.csv', sep=',', index=False)
    else:
        dataset = int(input("Do you want to analyse dataset 1 or dataset 2 ? "))
        if dataset == 1:
            ladderplot(files_list_1, reference_file_name_1, dataset)
            for i in range (n):
                help_list = get_mean(files_list_1[i], reference_file_name_1)
                means_list.append(help_list[0])
                std_list.append(help_list[1])
                N_list.append(help_list[2])
            xerror = std_list/(np.sqrt(N_list))
            peek = {'Material Budget Peek': material_budget_list, ' Material Budget error': material_budget_error_list, ' mean-squared deviation angle from reference beam': means_list, ' xerror': xerror}
            plot_data_peek = pd.DataFrame(data=peek)
            plot_data_peek.to_csv('Calibration_data/calibration_data_peek_dataset1.csv', sep=',', index=False)
        elif dataset == 2:
            ladderplot(files_list_2, reference_file_name_2, dataset)
            for i in range (n):
                help_list = get_mean(files_list_2[i], reference_file_name_2)
                means_list.append(help_list[0])
                std_list.append(help_list[1])
                N_list.append(help_list[2])
            xerror = std_list/(np.sqrt(N_list))
            peek = {'Material Budget Peek': material_budget_list, ' Material Budget error': material_budget_error_list, ' mean-squared deviation angle from reference beam': means_list, ' xerror': xerror}
            plot_data_peek = pd.DataFrame(data=peek)
            plot_data_peek.to_csv('Calibration_data/calibration_data_peek_dataset2.csv', sep=',', index=False)
        else:
            print("invalid dataset")
            sys.exit()

    fig, ax = plt.subplots(figsize=(10,10), layout='constrained')
    ax.scatter(means_list, material_budget_list, label = 'Burj (PEEK)') #Plotting data onto the axes
    ax.errorbar(means_list, material_budget_list, yerr = material_budget_error_list , xerr = std_list/(np.sqrt(N_list)) , fmt="o")
    ax.plot(highland_prediction_list, material_budget_list, label = 'Highland')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$mean^2$ [$mrad^2$]', style='normal')
    ax.set_ylabel('Material Budget')
    plt.xlim([10, 1000])
    ax.legend()
    if dataset == 1:
        ax.set_title('Calibration Plot Dataset 1')
        plt.savefig('Burj_calibration_plot_dataset1')
        print("Calibration plot saved as Burj_calibration_plot_dataset1")
    elif dataset ==2:
        ax.set_title('Calibration Plot Dataset 2')
        plt.savefig('Burj_calibration_plot_dataset2')
        print("Calibration plot saved as Burj_calibration_plot_dataset2")
    else:  
        ax.set_title('Calibration Plot both datasets merged')
        plt.savefig('Burj_calibration_plot_merged')
        print("Calibration plot saved as Burj_calibration_plot_merged")
    plt.show()