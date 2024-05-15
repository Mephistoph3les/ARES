# execute with      python3 scattering_calibration_plot_metal_sheets.py

from numba import jit, cuda 
import pandas as pd 
import matplotlib.pyplot as plt 
import math
import numpy as np
import time



def time_cut_of_data(file_name, logbook_name, log_starting_point, file_number, list_of_x_positions, list_of_material_thicknesses):
    data = pd.read_csv(file_name, sep='\s+', names= ['Event', 'Timestamp', 'Width', 'Intensity'])
    logbook = pd.read_csv(logbook_name, sep='\s+', skiprows=(16+log_starting_point-1) , names= ['Event', 'X Position (mm)', 'Y Index', '[ADD 1] Position (deg)', '[ADD 2] Position (mm)' , '[ADD 3] Count ()', '[ADD 4] Charge (pC)', 'UNIX Time Stamp (seconds)'])
    logbook['Event'] = logbook['Event'].astype(int)
    logbook_list = ([])
    help_list = ([])

    for i in range(len(list_of_x_positions)):
        if i==0:
            help_list.extend((logbook.loc[logbook['X Position (mm)'] == list_of_x_positions[i]], logbook.loc[logbook['X Position (mm)'] == 266]))
            help = pd.concat(help_list)
            logbook_list.append(help)
        else:
            logbook_list.append(logbook.loc[logbook['X Position (mm)'] == list_of_x_positions[i]])

    list_of_boundaries = [get_boundaries(logbook) for logbook in logbook_list]
    
    file_details = {
        0: (1712160150, 'Nickel_data/output_measurement_203_amended.csv'),
        1: (1712224363, 'Aluminum_data/output_measurement_211_amended.csv'),
        2: (1712235885.5, 'Aluminum_data/output_measurement_221_amended.csv')}
    timestamp_offset, output_filename = file_details[file_number]

    data_timestamps = data.loc[:, "Timestamp"]
    length_of_dataset = data_timestamps.shape[0]
    data_timestamp_array = data_timestamps.to_numpy()

    print(list_of_x_positions)

    scan = timestamp_scan(list_of_x_positions, list_of_material_thicknesses, list_of_boundaries, timestamp_offset, length_of_dataset, data_timestamp_array)
    thickness_region_allocation = scan[0]
    thickness_allocation = scan[1]
    rows_to_drop = data.shape[0]- len(thickness_region_allocation)
    data.drop(data.tail(rows_to_drop).index, inplace = True)
    print("dropped ", rows_to_drop, " rows from data")
    data.insert(4, "Stage x Position", thickness_region_allocation, True)
    data.insert(5, "Material Thickness", thickness_allocation, True)
    print(data)
    data.to_csv(output_filename, sep=',', index=False) 

    return data
@jit
def timestamp_scan(list_of_x_positions, list_of_material_thicknesses, list_of_boundaries, timestamp_offset, length_of_dataset, data_timestamp_array):
    thickness_region_allocation = []
    thickness_allocation = []
    
    for i in range(length_of_dataset):                                  #iteration over every datapoint of the measurement
        found = False
        for m in range(len(list_of_boundaries)):                        #iteration over every x position of the stage in logfile corresponding to material thickness 
            if found == True:
                break
            for l in range(len(list_of_boundaries[m][0])):              #iteration over every start and end point of measurement repetition
                start = list_of_boundaries[m][0][l] - timestamp_offset
                end = list_of_boundaries[m][1][l] - timestamp_offset
                print("Comparing:", start, "<= ", data_timestamp_array[i], "<=", end)
                if start <= data_timestamp_array[i] <= end:
                    thickness_region_allocation.append(list_of_x_positions[m])
                    thickness_allocation.append(list_of_material_thicknesses[m])
                    print("added", list_of_x_positions[m], "to slot: ", i)
                    found = True
                    break
                else:
                    print("skipped width for: ", i)
        if not found:
            thickness_region_allocation.append(0)
            thickness_allocation.append(0)

    return(thickness_region_allocation, thickness_allocation)

def get_boundaries(logbook_list):    
    log = logbook_list.loc[:, "Event"]              #Take a slice with only the Event numbers
    print(log)
    array = log.to_numpy()                          #convert event numbers into numpy array to compare consecutive entries if there is a skip to determine boundaries
    print(array)
    i=1
    list_of_start_event_boundaries=([array[0]])
    list_of_end_event_boundaries=([])
    while i< len(array)-1:
        if array[i] > array[i-1]+1:
            list_of_start_event_boundaries.append(array[i])
            list_of_end_event_boundaries.append(array[i-1])
        i+=1
    list_of_end_event_boundaries.append(array[(len(array)-1)])
    
    list_of_start_time_boundaries=([])                  #When the corresponding events are found, those rows are searched and the timestamp will be given
    list_of_end_time_boundaries=([])
    for i in range(len(list_of_start_event_boundaries)):
        event_start_frame = logbook_list.loc[logbook_list['Event'] == list_of_start_event_boundaries[i]]
        event_end_frame = logbook_list.loc[logbook_list['Event'] == list_of_end_event_boundaries[i]]

        event_start_frame = event_start_frame.iloc[0, 7]
        event_end_frame = event_end_frame.iloc[0, 7]
        list_of_start_time_boundaries.append(event_start_frame)
        list_of_end_time_boundaries.append(event_end_frame)

    return list_of_start_time_boundaries, list_of_end_time_boundaries

def get_mean(file_name, reference_file_name):
    
    #reading in the data
    data = pd.read_csv(file_name, sep='\s+', names= ['Event', 'Timestamp', 'Width', 'Intensity'])
    #reading in reference beam data to subtract from actual data for beam widening
    reference = pd.read_csv(reference_file_name, sep='\s+', names= ['Event', 'Timestamp', 'Width', 'Intensity'])

    ref_mean = reference.loc[:, 'Width'].mean()
    
    d = 75*10**(-3) # distance burj <-> detector in m
    data.loc[:, ["Width"]] = (np.sqrt((data.loc[:, ["Width"]])**2 - (ref_mean)**2)) *(55*10**(-6))     #background deduction
    data.loc[:, ["Width"]] = (np.arctan(data.loc[:, ["Width"]]/d)*10**(+3))
    data_mean = data.loc[:, 'Width'].mean() 
    data_std = data.loc[:, 'Width'].std() 
    N = data.shape[0]

    return data_mean**2, data_std, N

def separate_data(data, list_of_x_positions):
    data_separated_list = ([])
    for i in range(len(list_of_x_positions)):
        data_separated_list.append(data.loc[data['Material Thickness'] == list_of_x_positions[i]])
    return data_separated_list

def ladderplot(data_list, dataset, list_of_x_positions):
    material = "Nickel" if dataset==0 else "Aluminum"
    widths=[]
    material_thicknesses=[]
    help=[]
    for i in range(len(list_of_x_positions)-1):
        help = data_list[i].loc[:, "Width"]
        help = help.to_numpy()
        widths.append(help)
        
        help2 = data_list[i].loc[:, "Material Thickness"]
        help2 = help2.to_numpy()
        material_thicknesses.append(help2)

    fig, ax = plt.subplots(figsize=(10,10), layout='constrained')
    for i in range(len(list_of_x_positions)-1):
        ax.scatter(widths[i], material_thicknesses[i], label = ('Stage x position: ', list_of_x_positions[i]))
    ax.set_xlabel('Width', style='normal')
    ax.set_ylabel('Material Thickness')
    ax.set_title(('Ladderplot', material, ' of dataset ', dataset ))
    ax.legend()
    plt.show()

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
    start_time = time.time()
    list_of_x_positions = ([35.0, 51.5 , 68.0, 84.5, 101.0, 117.5, 134.0, 150.5, 167.0, 183.5, 200.0, 216.5, 233.0, 249.5])       #266 is skipped because that is the second zero/reference measurement
    radiation_length_nickel = 35.3                                      #in mm
    error_in_radiation_length = 0.05                                    # ~2% error
    radiation_length_aluminum = 8.9                                     #in mm
    list_of_material_thicknesses = ([0.0 ,0.025, 0.075, 0.05, 0.15, 0.1, 0.35, 0.25, 0.75, 0.5, 1.5, 1.0, 3.0, 2.0])     #in mm (CAREFUL! ALUMINUM HAS DIFFERENT THICKNESSES WHICH ARE NOT IMPLEMENTED RN)
    files_list=(['Nickel_data/output_measurement_203.dat', 'Aluminum_data/output_measurement_211.dat', 'Aluminum_data/output_measurement_221.dat'])
    ask = input("Do you want to do the time slicing of measurement data (no skips to analysis)?  ")
    if ask.lower() in ["y","yes"]:
        print("0 for nickel overnight \n1 for aluminum scan1 \n2 for aluminum scan2 ")
        dataset = int(input("Which dataset do you want to slice? "))
        answer = input("This computation will take a long time to do, continue? ")
        if answer.lower() in ["y","yes"]:
            logbook_files_list = (['Nickel_data/eCTLogger_ni_overnight.txt', 'Aluminum_data/eCTLogger_al_shortscan1.txt', 'Aluminum_data/eCTLogger_al_shortscan2.txt'])
            list_of_log_starting_points = ([17, 654, 492])
            data = time_cut_of_data(files_list[dataset],logbook_files_list[dataset],list_of_log_starting_points[dataset], dataset, list_of_x_positions, list_of_material_thicknesses)
    else:
        print("0 for nickel overnight \n1 for aluminum scan1 \n2 for aluminum scan2 ")
        dataset = int(input("Which dataset do you want to analyse? "))
        files_list_amen=(['Nickel_data/output_measurement_203_amended.csv', 'Aluminum_data/output_measurement_211_amended.csv', 'Aluminum_data/output_measurement_221_amended.csv'])
        data = pd.read_csv(files_list_amen[dataset], sep=',', skiprows=(1), names= ['Event', 'Timestamp', 'Width', 'Intensity', 'Stage x Position', 'Material Thickness'])
        data_list = separate_data(data, list_of_x_positions)        
        ask = input("Do you want to create a ladderplot? ")
        if ask.lower() in ["y","yes"]:
            ladderplot(data_list, dataset, list_of_x_positions)


    
#    list_of_material_thicknesses_shortscan = ([2, 3, 1, 1.5, 0.5, 0.75, 0.25, 0.35, 0.1, 0.15, 0.05, 0.075, 0.025, 0])     #in mm
    error_in_material_thickness = 0.0005                                   #error in mm
    
    m = len(list_of_material_thicknesses)

    material_budget_list = ([])
    highland_prediction_list = ([])
    material_budget_error_list = ([])
    for i in range (1, m-1):
        material_budget_list.append(list_of_material_thicknesses[i]/radiation_length_nickel)
        material_budget_error = np.sqrt(((error_in_material_thickness/radiation_length_nickel)**2+(error_in_radiation_length*list_of_material_thicknesses[i]/(radiation_length_nickel)**2)**2))
        material_budget_error_list.append(material_budget_error)
        highland_prediction_list.append(highland(material_budget_list[i-1])[0])
    #print("Thats the list of material budgets: ")
    #print(material_budget_list)
    #print("Thats the list of material budget errors: ")
    #print(material_budget_error_list)
    #print("Thats the list of highlands predictions: ")
    #print(highland_prediction_list)
    #print("")
    #print("beta is: ", highland(material_budget_list[0])[1])
    #print("electron momentum is: ", highland(material_budget_list[0])[2], "MeV")



    #fig, ax = plt.subplots(figsize=(10,10), layout='constrained')
    #ax.hist( material_budget_list, label = 'Nickel') #Plotting data onto the axes
    #ax.errorbar( material_budget_list, yerr = material_budget_error_list , fmt="o")
    #ax.plot(highland_prediction_list, material_budget_list, label = 'Highland')
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #ax.set_xlabel('$mean^2$ [$mrad^2$]', style='normal')
    #ax.set_ylabel('Material Budget')
    #ax.set_title('Calibration Plot')
    #ax.legend()
    #plt.show()
    end_time = time.time()
    print("Computing time was:  ", round(end_time - start_time, 2), " s")
    print("Which in minutes is: ", round((end_time - start_time)/60, 2) , " min")