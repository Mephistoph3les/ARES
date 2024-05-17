<<<<<<< HEAD
# execute with      python3 scattering_calibration_plot_metal_sheets.py

from numba import jit, cuda 
import pandas as pd 
import matplotlib.pyplot as plt 
import math
import numpy as np
import time
import csv



def time_cut_of_data(file_name, logbook_name, log_starting_point, file_number, list_of_x_positions, list_of_material_thicknesses):
    data = pd.read_csv(file_name, sep='\s+', names= ['Event', 'Timestamp', 'Width', 'Intensity'])
    logbook = pd.read_csv(logbook_name, sep='\s+', skiprows=(16+log_starting_point-1) , names= ['Event', 'X Position (mm)', 'Y Index', '[ADD 1] Position (deg)', '[ADD 2] Position (mm)' , '[ADD 3] Count ()', '[ADD 4] Charge (pC)', 'UNIX Time Stamp (seconds)'])
    logbook['Event'] = logbook['Event'].astype(int)
    logbook_list = ([])

    for i in range(len(list_of_x_positions)):
        if i==1:
            help_list=([logbook.loc[logbook['X Position (mm)'] == list_of_x_positions[i]], logbook.loc[logbook['X Position (mm)'] == 51.5]])
            help = pd.concat(help_list)             #Measurements on stage position 51.5 and 68.0 are actually the same thickness, here these two are merged
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

    scan = timestamp_scan(list_of_x_positions, list_of_material_thicknesses, list_of_boundaries, timestamp_offset, length_of_dataset, data_timestamp_array)
    thickness_region_allocation = scan[0]
    thickness_allocation = scan[1]
    rows_to_drop = data.shape[0]- len(thickness_region_allocation)
    data.drop(data.tail(rows_to_drop).index, inplace = True)
    print("dropped ", rows_to_drop, " rows from data to be of same length with analysis")
    data.insert(4, "Stage x Position", thickness_region_allocation, True)
    data.insert(5, "Material Thickness", thickness_allocation, True)
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
    array = log.to_numpy()                          #convert event numbers into numpy array to compare consecutive entries if there is a skip to determine boundaries
    list_of_start_event_boundaries=([array[0]])
    list_of_end_event_boundaries=([])
    for i in range(1, len(array)):
        if array[i] > array[i-1]+1:
            list_of_start_event_boundaries.append(array[i])
            list_of_end_event_boundaries.append(array[i-1])
    list_of_end_event_boundaries.append(array[(len(array)-1)])

    start_time_boundaries=([])                  #When the corresponding events are found, those rows are searched and the timestamp will be given
    end_time_boundaries=([])
    for i in range(len(list_of_start_event_boundaries)):
        event_start_frame = logbook_list.loc[logbook_list['Event'] == list_of_start_event_boundaries[i]]
        event_end_frame = logbook_list.loc[logbook_list['Event'] == list_of_end_event_boundaries[i]]

        event_start_frame = event_start_frame.iloc[0, 7]
        event_end_frame = event_end_frame.iloc[0, 7]
        start_time_boundaries.append(event_start_frame)
        end_time_boundaries.append(event_end_frame)

    return start_time_boundaries, end_time_boundaries

def drop_unphysical_data(data_list):
    new_data = ([])
    for i in range(len(data_list)):
        help = data_list[i]
        help = help[help['Width'] < 100.0]
        new_data.append(help)
    return new_data

def get_mean(data, reference_mean):
    
    d = 75*10**(-3) # distance metal sheets <-> detector in m
    data.loc[:, ["Width"]] = (np.sqrt((data.loc[:, ["Width"]])**2 - (reference_mean[0])**2)) *(55*10**(-6))     #background deduction
    data.loc[:, ["Width"]] = (np.arctan(data.loc[:, ["Width"]]/d)*10**(+3))
    data_mean = data.loc[:, 'Width'].mean() 
    data_std = data.loc[:, 'Width'].std() 
    N = data.shape[0]

    return data_mean**2, data_std, N

def separate_data(data, list_of_x_positions):
    return [
        data.loc[data['Stage x Position'] == list_of_x_positions[i]]
        for i in range(len(list_of_x_positions))
    ]

def ladderplot(data_list, dataset, list_of_x_positions, list_of_material_thicknesses):
    material = "Nickel" if dataset==0 else "Aluminum"
    widths=[]
    material_thicknesses=[]
    help=[]
    ladder_filename = f"ladderplot_{material.lower()}_dataset_{str(dataset)}"

    for i in range(len(list_of_x_positions)):
        help = data_list[i].loc[:, "Width"]
        help = help.to_numpy()
        widths.append(help)

        help2 = data_list[i].loc[:, "Stage x Position"]
        help2 = help2.to_numpy()
        material_thicknesses.append(help2)

    fig, ax = plt.subplots(figsize=(10,10), layout='constrained')
    for i in range(len(list_of_x_positions)):
        ax.scatter(widths[i], material_thicknesses[i], label = ('Material Thickness (in mm): ', list_of_material_thicknesses[i]))
    ax.set_xlabel('Width (in pixels)', style='normal')
    ax.set_ylabel('Stage x Position (in mm)')
    ax.set_title(f'Ladderplot {material} of dataset {str(dataset)}')
    ax.legend()
    ax.grid()
    ax.set_axisbelow(True)
    if dataset == 1 or dataset == 2:
        plt.xlim([0, 35])
    else:  
        plt.xlim([0, 90])
    plt.savefig(ladder_filename)
    print("Ladderplot saved as " + ladder_filename)
    plt.show()

def calibrationplot(material_nickel, material_aluminum, means_list_nickel, means_list_aluminum, std_list_nickel, std_list_aluminum, N_list_nickel, N_list_aluminum):   
    radiation_length_nickel = 14.24                       #in mm                                   
    radiation_length_aluminum = 88.97                     #in mm   
    error_in_radiation_length = 0.002                     #considering the decimal value given by pdg 
    error_in_material_thickness = 0.0002                  #error in mm

    m = len(list_of_material_thicknesses)

    material_budget_nickel = ([])
    material_budget_aluminum = ([])
    highland_prediction_list = ([])
    material_budget_for_highland = ([])
    highland_prediction_list_nickel = ([])
    highland_prediction_list_aluminum = ([])
    material_budget_error_nickel = ([])
    material_budget_error_aluminum = ([])
    for i in range (1, m):
        material_budget_nickel.append(material_nickel[i]/radiation_length_nickel)
        material_budget_aluminum.append(material_aluminum[i]/radiation_length_aluminum)
        material_budget_error_nickel_i = np.sqrt(((error_in_material_thickness/radiation_length_nickel)**2+(error_in_radiation_length*material_nickel[i]/(radiation_length_nickel)**2)**2))
        material_budget_error_aluminum_i = np.sqrt(((error_in_material_thickness/radiation_length_aluminum)**2+(error_in_radiation_length*material_aluminum[i]/(radiation_length_aluminum)**2)**2))
        material_budget_error_nickel.append(material_budget_error_nickel_i)
        material_budget_error_aluminum.append(material_budget_error_aluminum_i)
        highland_prediction_list_nickel.append(highland(material_budget_nickel[i-1])[0])       #retrieves Theta**2 from Highland prediction
        highland_prediction_list_aluminum.append(highland(material_budget_aluminum[i-1])[0])       #retrieves Theta**2 from Highland prediction

    highland_prediction_list.extend(
        (highland_prediction_list_nickel, highland_prediction_list_aluminum))
    material_budget_for_highland.extend(
        (material_budget_nickel, material_budget_aluminum))
    nickel = {'Material Budget Nickel': material_budget_nickel, ' Material Budget error': material_budget_error_nickel, ' mean-squared deviation angle from reference beam': means_list_nickel, ' xerror': std_list_nickel/(np.sqrt(N_list_nickel))}
    plot_data_nickel = pd.DataFrame(data=nickel)
    aluminum = {'Material Budget Aluminum': material_budget_aluminum, ' Material Budget error': material_budget_error_aluminum, ' mean-squared deviation angle from reference beam': means_list_aluminum, ' xerror': std_list_aluminum/(np.sqrt(N_list_aluminum))}
    plot_data_aluminum = pd.DataFrame(data=aluminum)
    plot_data_nickel.to_csv('Calibration_data/calibration_data_nickel.csv', sep=',', index=False)
    plot_data_aluminum.to_csv('Calibration_data/calibration_data_aluminum.csv', sep=',', index=False) 


    fig, ax = plt.subplots(figsize=(10,10), layout='constrained')
    ax.scatter( means_list_nickel, material_budget_nickel, label = 'Nickel') #Plotting data onto the axes
    ax.errorbar( means_list_nickel, material_budget_nickel, yerr = material_budget_error_nickel , xerr = std_list_nickel/(np.sqrt(N_list_nickel)), fmt="o")
    ax.scatter( means_list_aluminum, material_budget_aluminum, label = 'Aluminum') #Plotting data onto the axes
    ax.errorbar( means_list_aluminum, material_budget_aluminum, yerr = material_budget_error_aluminum , xerr = std_list_aluminum/(np.sqrt(N_list_aluminum)), fmt="o")
    ax.plot(highland_prediction_list, material_budget_for_highland, label = 'Highland Prediction')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$mean^2$ deviation angle from reference beam [$mrad^2$]', style='normal')            #think of a new name for mean... it's technically not a mean
    ax.set_ylabel('Material Budget')
    ax.set_title('Calibration Plot')
    ax.legend()
    ax.grid()
    ax.set_axisbelow(True)
    #plt.xlim([0, 35])
    plt.savefig('metal_sheet_calibration_plot')
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
    list_of_x_positions = ([35.0, 68.0, 84.5, 101.0, 117.5, 134.0, 150.5, 167.0, 183.5, 200.0, 216.5, 233.0, 249.5, 266.0])       #266 is skipped because that is the second zero/reference measurement
    list_of_material_thicknesses_nickel = ([0.0, 0.025, 0.075, 0.05, 0.15, 0.1, 0.35, 0.25, 0.75, 0.5, 1.5, 1.0, 3.0, 2.0])     #in mm 
    list_of_material_thicknesses_aluminum = ([0.0, 0.025, 0.095, 0.07, 0.27, 0.2, 0.6, 0.4, 1.4, 1.0, 3.0, 2.0, 6.0, 4.0])
    files_list=(['Nickel_data/output_measurement_203.dat', 'Aluminum_data/output_measurement_211.dat', 'Aluminum_data/output_measurement_221.dat'])
    ask = input("Do you want to do the time slicing of measurement data (no skips to analysis)?  ")
    if ask.lower() in ["y","yes"]:
        print("0 for nickel overnight \n1 for aluminum scan1 \n2 for aluminum scan2 ")
        dataset = int(input("Which dataset do you want to slice? "))
        if dataset == 0:
            list_of_material_thicknesses = list_of_material_thicknesses_nickel
        else:
            list_of_material_thicknesses = list_of_material_thicknesses_aluminum
        answer = input("This computation will take a long time to do, continue? ")
        if answer.lower() in ["y","yes"]:
            logbook_files_list = (['Nickel_data/eCTLogger_ni_overnight.txt', 'Aluminum_data/eCTLogger_al_shortscan1.txt', 'Aluminum_data/eCTLogger_al_shortscan2.txt'])
            list_of_log_starting_points = ([17, 654, 1])
            data = time_cut_of_data(files_list[dataset], logbook_files_list[dataset], list_of_log_starting_points[dataset], dataset, list_of_x_positions, list_of_material_thicknesses)
    else:
        print("0 for nickel overnight \n1 for aluminum scan1 \n2 for aluminum scan2 ")
        dataset = int(input("Which dataset do you want to analyse? "))
        if dataset == 0:
            list_of_material_thicknesses = list_of_material_thicknesses_nickel
        else:
            list_of_material_thicknesses = list_of_material_thicknesses_aluminum
        files_list_amen=(['Nickel_data/output_measurement_203_amended.csv', 'Aluminum_data/output_measurement_211_amended.csv', 'Aluminum_data/output_measurement_221_amended.csv'])
        data = pd.read_csv(files_list_amen[dataset], sep=',', skiprows=(1), names= ['Event', 'Timestamp', 'Width', 'Intensity', 'Stage x Position', 'Material Thickness'])
        data_list = separate_data(data, list_of_x_positions) 
        
        ask = input("Do you want to drop unphysical data? ")    #all data points with width > 100 will be dropped
        if ask.lower() in ["y","yes"]:
            data_list = drop_unphysical_data(data_list)
        
        ask = input("Do you want to create a ladderplot? ")
        if ask.lower() in ["y","yes"]:
            ladderplot(data_list, dataset, list_of_x_positions, list_of_material_thicknesses)

        ask = input("Do you want to produce a calibration plot? ")
        if ask.lower() in ["y","yes"]:
            data_nickel = pd.read_csv(files_list_amen[0], sep=',', skiprows=(1), names= ['Event', 'Timestamp', 'Width', 'Intensity', 'Stage x Position', 'Material Thickness'])
            data_list_nickel = separate_data(data_nickel, list_of_x_positions)
            data_list_nickel = drop_unphysical_data(data_list_nickel)
            data_aluminum = pd.read_csv(files_list_amen[2], sep=',', skiprows=(1), names= ['Event', 'Timestamp', 'Width', 'Intensity', 'Stage x Position', 'Material Thickness'])
            data_list_aluminum = separate_data(data_aluminum, list_of_x_positions)
            data_list_aluminum = drop_unphysical_data(data_list_aluminum)
            
            means_list_nickel = ([])
            std_list_nickel = ([])
            N_list_nickel = ([])
            means_list_aluminum = ([])
            std_list_aluminum = ([])
            N_list_aluminum = ([])
            
            reference_nickel = ([])
            reference_aluminum = ([])
            data_ref_nickel = data_list_nickel[0]
            reference_nickel.append(data_ref_nickel.loc[:, 'Width'].mean()) 
            reference_nickel.append(data_ref_nickel.loc[:, 'Width'].std())
            reference_nickel.append(data_ref_nickel.shape[0])
            data_ref_aluminum = data_list_aluminum[0]
            reference_aluminum.append(data_ref_aluminum.loc[:, 'Width'].mean())
            reference_aluminum.append(data_ref_aluminum.loc[:, 'Width'].std())
            reference_aluminum.append(data_ref_aluminum.shape[0])

            
            for i in range (1, len(data_list)):
                means_list_nickel.append(get_mean(data_list_nickel[i], reference_nickel)[0])
                std_list_nickel.append(get_mean(data_list_nickel[i], reference_nickel)[1])
                N_list_nickel.append(get_mean(data_list_nickel[i], reference_nickel)[2])
                means_list_aluminum.append(get_mean(data_list_aluminum[i], reference_aluminum)[0])
                std_list_aluminum.append(get_mean(data_list_aluminum[i], reference_aluminum)[1])
                N_list_aluminum.append(get_mean(data_list_aluminum[i], reference_aluminum)[2])
                
            calibrationplot(list_of_material_thicknesses_nickel, list_of_material_thicknesses_aluminum, means_list_nickel, means_list_aluminum, std_list_nickel, std_list_aluminum, N_list_nickel, N_list_aluminum)
    end_time = time.time()
    print("Computing time was:  ", round(end_time - start_time, 2), " s")
    print("Which in minutes is: ", round((end_time - start_time)/60, 2) , " min")
=======
# execute with      python3 scattering_calibration_plot_metal_sheets.py

from numba import jit, cuda 
import pandas as pd 
import matplotlib.pyplot as plt 
import math
import numpy as np
import time
import csv



def time_cut_of_data(file_name, logbook_name, log_starting_point, file_number, list_of_x_positions, list_of_material_thicknesses):
    data = pd.read_csv(file_name, sep='\s+', names= ['Event', 'Timestamp', 'Width', 'Intensity'])
    logbook = pd.read_csv(logbook_name, sep='\s+', skiprows=(16+log_starting_point-1) , names= ['Event', 'X Position (mm)', 'Y Index', '[ADD 1] Position (deg)', '[ADD 2] Position (mm)' , '[ADD 3] Count ()', '[ADD 4] Charge (pC)', 'UNIX Time Stamp (seconds)'])
    logbook['Event'] = logbook['Event'].astype(int)
    logbook_list = ([])

    for i in range(len(list_of_x_positions)):
        if i==1:
            help_list=([logbook.loc[logbook['X Position (mm)'] == list_of_x_positions[i]], logbook.loc[logbook['X Position (mm)'] == 51.5]])
            help = pd.concat(help_list)             #Measurements on stage position 51.5 and 68.0 are actually the same thickness, here these two are merged
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

    scan = timestamp_scan(list_of_x_positions, list_of_material_thicknesses, list_of_boundaries, timestamp_offset, length_of_dataset, data_timestamp_array)
    thickness_region_allocation = scan[0]
    thickness_allocation = scan[1]
    rows_to_drop = data.shape[0]- len(thickness_region_allocation)
    data.drop(data.tail(rows_to_drop).index, inplace = True)
    print("dropped ", rows_to_drop, " rows from data to be of same length with analysis")
    data.insert(4, "Stage x Position", thickness_region_allocation, True)
    data.insert(5, "Material Thickness", thickness_allocation, True)
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
    array = log.to_numpy()                          #convert event numbers into numpy array to compare consecutive entries if there is a skip to determine boundaries
    list_of_start_event_boundaries=([array[0]])
    list_of_end_event_boundaries=([])
    for i in range(1, len(array)):
        if array[i] > array[i-1]+1:
            list_of_start_event_boundaries.append(array[i])
            list_of_end_event_boundaries.append(array[i-1])
    list_of_end_event_boundaries.append(array[(len(array)-1)])

    start_time_boundaries=([])                  #When the corresponding events are found, those rows are searched and the timestamp will be given
    end_time_boundaries=([])
    for i in range(len(list_of_start_event_boundaries)):
        event_start_frame = logbook_list.loc[logbook_list['Event'] == list_of_start_event_boundaries[i]]
        event_end_frame = logbook_list.loc[logbook_list['Event'] == list_of_end_event_boundaries[i]]

        event_start_frame = event_start_frame.iloc[0, 7]
        event_end_frame = event_end_frame.iloc[0, 7]
        start_time_boundaries.append(event_start_frame)
        end_time_boundaries.append(event_end_frame)

    return start_time_boundaries, end_time_boundaries

def drop_unphysical_data(data_list):
    new_data = ([])
    for i in range(len(data_list)):
        help = data_list[i]
        help = help[help['Width'] < 100.0]
        new_data.append(help)
    return new_data

def get_mean(data, reference_mean):
    
    d = 75*10**(-3) # distance metal sheets <-> detector in m
    data.loc[:, ["Width"]] = (np.sqrt((data.loc[:, ["Width"]])**2 - (reference_mean[0])**2)) *(55*10**(-6))     #background deduction
    data.loc[:, ["Width"]] = (np.arctan(data.loc[:, ["Width"]]/d)*10**(+3))
    data_mean = data.loc[:, 'Width'].mean() 
    data_std = data.loc[:, 'Width'].std() 
    N = data.shape[0]

    return data_mean**2, data_std, N

def separate_data(data, list_of_x_positions):
    return [
        data.loc[data['Stage x Position'] == list_of_x_positions[i]]
        for i in range(len(list_of_x_positions))
    ]

def ladderplot(data_list, dataset, list_of_x_positions, list_of_material_thicknesses):
    material = "Nickel" if dataset==0 else "Aluminum"
    widths=[]
    material_thicknesses=[]
    help=[]
    ladder_filename = f"ladderplot_{material.lower()}_dataset_{str(dataset)}"

    for i in range(len(list_of_x_positions)):
        help = data_list[i].loc[:, "Width"]
        help = help.to_numpy()
        widths.append(help)

        help2 = data_list[i].loc[:, "Stage x Position"]
        help2 = help2.to_numpy()
        material_thicknesses.append(help2)

    fig, ax = plt.subplots(figsize=(10,10), layout='constrained')
    for i in range(len(list_of_x_positions)):
        ax.scatter(widths[i], material_thicknesses[i], label = ('Material Thickness (in mm): ', list_of_material_thicknesses[i]))
    ax.set_xlabel('Width (in pixels)', style='normal')
    ax.set_ylabel('Stage x Position (in mm)')
    ax.set_title(f'Ladderplot {material} of dataset {str(dataset)}')
    ax.legend()
    ax.grid()
    ax.set_axisbelow(True)
    if dataset == 1 or dataset == 2:
        plt.xlim([0, 35])
    else:  
        plt.xlim([0, 90])
    plt.savefig(ladder_filename)
    print("Ladderplot saved as " + ladder_filename)
    plt.show()

def calibrationplot(material_nickel, material_aluminum, means_list_nickel, means_list_aluminum, std_list_nickel, std_list_aluminum, N_list_nickel, N_list_aluminum):   
    radiation_length_nickel = 14.24                       #in mm                                   
    radiation_length_aluminum = 88.97                     #in mm   
    error_in_radiation_length = 0.002                     #considering the decimal value given by pdg 
    error_in_material_thickness = 0.0002                  #error in mm

    m = len(list_of_material_thicknesses)

    material_budget_nickel = ([])
    material_budget_aluminum = ([])
    highland_prediction_list = ([])
    material_budget_for_highland = ([])
    highland_prediction_list_nickel = ([])
    highland_prediction_list_aluminum = ([])
    material_budget_error_nickel = ([])
    material_budget_error_aluminum = ([])
    for i in range (1, m):
        material_budget_nickel.append(material_nickel[i]/radiation_length_nickel)
        material_budget_aluminum.append(material_aluminum[i]/radiation_length_aluminum)
        material_budget_error_nickel_i = np.sqrt(((error_in_material_thickness/radiation_length_nickel)**2+(error_in_radiation_length*material_nickel[i]/(radiation_length_nickel)**2)**2))
        material_budget_error_aluminum_i = np.sqrt(((error_in_material_thickness/radiation_length_aluminum)**2+(error_in_radiation_length*material_aluminum[i]/(radiation_length_aluminum)**2)**2))
        material_budget_error_nickel.append(material_budget_error_nickel_i)
        material_budget_error_aluminum.append(material_budget_error_aluminum_i)
        highland_prediction_list_nickel.append(highland(material_budget_nickel[i-1])[0])       #retrieves Theta**2 from Highland prediction
        highland_prediction_list_aluminum.append(highland(material_budget_aluminum[i-1])[0])       #retrieves Theta**2 from Highland prediction

    highland_prediction_list.extend(
        (highland_prediction_list_nickel, highland_prediction_list_aluminum))
    material_budget_for_highland.extend(
        (material_budget_nickel, material_budget_aluminum))
    nickel = {'Material Budget Nickel': material_budget_nickel, ' Material Budget error': material_budget_error_nickel, ' mean-squared deviation angle from reference beam': means_list_nickel, ' xerror': std_list_nickel/(np.sqrt(N_list_nickel))}
    plot_data_nickel = pd.DataFrame(data=nickel)
    aluminum = {'Material Budget Aluminum': material_budget_aluminum, ' Material Budget error': material_budget_error_aluminum, ' mean-squared deviation angle from reference beam': means_list_aluminum, ' xerror': std_list_aluminum/(np.sqrt(N_list_aluminum))}
    plot_data_aluminum = pd.DataFrame(data=aluminum)
    plot_data_nickel.to_csv('Calibration_data/calibration_data_nickel.csv', sep=',', index=False)
    plot_data_aluminum.to_csv('Calibration_data/calibration_data_aluminum.csv', sep=',', index=False) 


    fig, ax = plt.subplots(figsize=(10,10), layout='constrained')
    ax.scatter( means_list_nickel, material_budget_nickel, label = 'Nickel') #Plotting data onto the axes
    ax.errorbar( means_list_nickel, material_budget_nickel, yerr = material_budget_error_nickel , xerr = std_list_nickel/(np.sqrt(N_list_nickel)), fmt="o")
    ax.scatter( means_list_aluminum, material_budget_aluminum, label = 'Aluminum') #Plotting data onto the axes
    ax.errorbar( means_list_aluminum, material_budget_aluminum, yerr = material_budget_error_aluminum , xerr = std_list_aluminum/(np.sqrt(N_list_aluminum)), fmt="o")
    ax.plot(highland_prediction_list, material_budget_for_highland, label = 'Highland Prediction')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$mean^2$ deviation angle from reference beam [$mrad^2$]', style='normal')            #think of a new name for mean... it's technically not a mean
    ax.set_ylabel('Material Budget')
    ax.set_title('Calibration Plot')
    ax.legend()
    ax.grid()
    ax.set_axisbelow(True)
    #plt.xlim([0, 35])
    plt.savefig('metal_sheet_calibration_plot')
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
    list_of_x_positions = ([35.0, 68.0, 84.5, 101.0, 117.5, 134.0, 150.5, 167.0, 183.5, 200.0, 216.5, 233.0, 249.5, 266.0])       #266 is skipped because that is the second zero/reference measurement
    list_of_material_thicknesses_nickel = ([0.0, 0.025, 0.075, 0.05, 0.15, 0.1, 0.35, 0.25, 0.75, 0.5, 1.5, 1.0, 3.0, 2.0])     #in mm 
    list_of_material_thicknesses_aluminum = ([0.0, 0.025, 0.095, 0.07, 0.27, 0.2, 0.6, 0.4, 1.4, 1.0, 3.0, 2.0, 6.0, 4.0])
    files_list=(['Nickel_data/output_measurement_203.dat', 'Aluminum_data/output_measurement_211.dat', 'Aluminum_data/output_measurement_221.dat'])
    ask = input("Do you want to do the time slicing of measurement data (no skips to analysis)?  ")
    if ask.lower() in ["y","yes"]:
        print("0 for nickel overnight \n1 for aluminum scan1 \n2 for aluminum scan2 ")
        dataset = int(input("Which dataset do you want to slice? "))
        if dataset == 0:
            list_of_material_thicknesses = list_of_material_thicknesses_nickel
        else:
            list_of_material_thicknesses = list_of_material_thicknesses_aluminum
        answer = input("This computation will take a long time to do, continue? ")
        if answer.lower() in ["y","yes"]:
            logbook_files_list = (['Nickel_data/eCTLogger_ni_overnight.txt', 'Aluminum_data/eCTLogger_al_shortscan1.txt', 'Aluminum_data/eCTLogger_al_shortscan2.txt'])
            list_of_log_starting_points = ([17, 654, 1])
            data = time_cut_of_data(files_list[dataset], logbook_files_list[dataset], list_of_log_starting_points[dataset], dataset, list_of_x_positions, list_of_material_thicknesses)
    else:
        print("0 for nickel overnight \n1 for aluminum scan1 \n2 for aluminum scan2 ")
        dataset = int(input("Which dataset do you want to analyse? "))
        if dataset == 0:
            list_of_material_thicknesses = list_of_material_thicknesses_nickel
        else:
            list_of_material_thicknesses = list_of_material_thicknesses_aluminum
        files_list_amen=(['Nickel_data/output_measurement_203_amended.csv', 'Aluminum_data/output_measurement_211_amended.csv', 'Aluminum_data/output_measurement_221_amended.csv'])
        data = pd.read_csv(files_list_amen[dataset], sep=',', skiprows=(1), names= ['Event', 'Timestamp', 'Width', 'Intensity', 'Stage x Position', 'Material Thickness'])
        data_list = separate_data(data, list_of_x_positions) 
        
        ask = input("Do you want to drop unphysical data? ")    #all data points with width > 100 will be dropped
        if ask.lower() in ["y","yes"]:
            data_list = drop_unphysical_data(data_list)
        
        ask = input("Do you want to create a ladderplot? ")
        if ask.lower() in ["y","yes"]:
            ladderplot(data_list, dataset, list_of_x_positions, list_of_material_thicknesses)

        ask = input("Do you want to produce a calibration plot? ")
        if ask.lower() in ["y","yes"]:
            data_nickel = pd.read_csv(files_list_amen[0], sep=',', skiprows=(1), names= ['Event', 'Timestamp', 'Width', 'Intensity', 'Stage x Position', 'Material Thickness'])
            data_list_nickel = separate_data(data_nickel, list_of_x_positions)
            data_list_nickel = drop_unphysical_data(data_list_nickel)
            data_aluminum = pd.read_csv(files_list_amen[1], sep=',', skiprows=(1), names= ['Event', 'Timestamp', 'Width', 'Intensity', 'Stage x Position', 'Material Thickness'])
            data_list_aluminum = separate_data(data_aluminum, list_of_x_positions)
            data_list_aluminum = drop_unphysical_data(data_list_aluminum)
            
            means_list_nickel = ([])
            std_list_nickel = ([])
            N_list_nickel = ([])
            means_list_aluminum = ([])
            std_list_aluminum = ([])
            N_list_aluminum = ([])
            
            reference_nickel = ([])
            reference_aluminum = ([])
            data_ref_nickel = data_list_nickel[0]
            reference_nickel.append(data_ref_nickel.loc[:, 'Width'].mean()) 
            reference_nickel.append(data_ref_nickel.loc[:, 'Width'].std())
            reference_nickel.append(data_ref_nickel.shape[0])
            data_ref_aluminum = data_list_aluminum[0]
            reference_aluminum.append(data_ref_aluminum.loc[:, 'Width'].mean())
            reference_aluminum.append(data_ref_aluminum.loc[:, 'Width'].std())
            reference_aluminum.append(data_ref_aluminum.shape[0])

            
            for i in range (1, len(data_list)):
                means_list_nickel.append(get_mean(data_list_nickel[i], reference_nickel)[0])
                std_list_nickel.append(get_mean(data_list_nickel[i], reference_nickel)[1])
                N_list_nickel.append(get_mean(data_list_nickel[i], reference_nickel)[2])
                means_list_aluminum.append(get_mean(data_list_aluminum[i], reference_aluminum)[0])
                std_list_aluminum.append(get_mean(data_list_aluminum[i], reference_aluminum)[1])
                N_list_aluminum.append(get_mean(data_list_aluminum[i], reference_aluminum)[2])
                
            calibrationplot(list_of_material_thicknesses_nickel, list_of_material_thicknesses_aluminum, means_list_nickel, means_list_aluminum, std_list_nickel, std_list_aluminum, N_list_nickel, N_list_aluminum)
    end_time = time.time()
    print("Computing time was:  ", round(end_time - start_time, 2), " s")
    print("Which in minutes is: ", round((end_time - start_time)/60, 2) , " min")
>>>>>>> bd4d72e87d9c6423d7038ddd563f8b148b211a04
