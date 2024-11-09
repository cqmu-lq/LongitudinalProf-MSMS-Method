import os, subprocess, pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyopenms import MSExperiment, MzMLFile, PeakPickerHiRes
import pymzml

# This file contains common methods, divided into the following modules: file reading/writing, data checking, data generation, data calculation, and data extraction

"""=============================
File Reading/Writing
"""

# File copy
import shutil
def copyfile(infile, outfile):
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        makedir(outdir)
        while True:
            if os.path.exists(outdir):
                break
    shutil.copy(infile, outfile)
    
# Create folder: create if it does not exist
def makedir(path):
    if not os.path.exists(path):
        subprocess.Popen('mkdir "{}"'.format(path), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

def make_folder(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

# Check if the file path exists: if not, create it, and return True after creation
def dir_is_exist(filepath):
    if os.path.exists(filepath):
        return True
    else:
        os.makedirs(filepath)
        while True:
            if os.path.exists(filepath):
                break
        return True

# Write to a pickle file
def to_pkl(data, w_path):
    with open(w_path, 'wb') as f:
        pickle.dump(data, f)

# Read a pickle file
def read_pkl(r_path):
    with open(r_path, 'rb') as f:
        data = pickle.load(f)
    return data
        
# Read mzML data
def readmzML(filepath):
    exp = MSExperiment()
    MzMLFile().load(filepath, exp)
    return exp

# Parse the mzML file
def get_TIC_from_File(file_path, tic_pkl):
    '''
    return [ {'rt', 'scanID',,,,},    ]
    '''
    if os.path.exists(tic_pkl):
        TIC = read_pkl(tic_pkl)
    else:
        # Load the file
        if isinstance(file_path, str):
            exp = MSExperiment()
            MzMLFile().load(file_path, exp)
        elif isinstance(file_path, MSExperiment):
            exp = file_path
        else:
            return

        profiled_exp = MSExperiment()
        temp_cen_list = []
        for spec in exp:
            if spec.getType() == 2:
                profiled_exp.addSpectrum(spec)
            else:
                temp_cen_list.append(spec)
        
        # Normalize data
        centroided_exp = MSExperiment()
        if profiled_exp.size() > 0:
            PeakPickerHiRes().pickExperiment(exp, centroided_exp)

        if len(temp_cen_list) > 0:
            for spec in temp_cen_list:
                centroided_exp.addSpectrum(spec)
        centroided_exp.sortSpectra()
        
        # Read with pymzml to get isolation window target m/z
        TIC2 = pymzml.run.Reader(file_path)

        TIC = []
        for idx, spec in enumerate(centroided_exp):
            spectrum = {}
            spectrum['ID'] = idx + 1
            spectrum['RT'] = spec.getRT()
            spectrum['Polarity'] = spec.getInstrumentSettings().getPolarity()
            spectrum['MSlevel'] = spec.getMSLevel()
            cur_id = idx if spectrum['MSlevel'] == 1 else cur_id
            cur_RT = spectrum['RT'] if spectrum['MSlevel'] == 1 else cur_RT
            spectrum['MS1_id'] = cur_id if spectrum['MSlevel'] == 2 else None
            spectrum['MS1_RT'] = cur_RT if spectrum['MSlevel'] == 2 else None
            
            # Modify this line
            spectrum['pecursor'] = None if spectrum['MSlevel'] == 1 else TIC2[idx + 1]['MS:1000827'] 
            
            spectrum['peaks'] = list(zip(*spec.get_peaks()))
            spectrum['TIC'] = np.sum(np.array(spectrum['peaks'])[:, 1])
            TIC.append(spectrum)
        
        to_pkl(TIC, tic_pkl)
    return TIC

# Merge all CSV results in the folder
def marge_data_fun(folder_path, qianzhui, begin_0, end_0, step_0, outfile):
    all_files_list = ["{}/{}_{}_{}.csv".format(folder_path, qianzhui, i, i + step_0) for i in range(begin_0, end_0, step_0)]
    merged_data_df = pd.DataFrame()
    for file in all_files_list:
        df = pd.read_csv(file, index_col=0)
        merged_data_df = pd.concat([merged_data_df, df], ignore_index=True)
    print(merged_data_df.shape, outfile)
    merged_data_df.to_excel(outfile, encoding='utf8')

"""
=============================
Data Checking
"""

# Data validation: 1. Count the number of MS1 and MS2 spectra; 2. Check if each MS1 spectrum is followed by 4 MS2 spectra
def check_data(tic_list):
    # Count the number of MS1 and MS2 spectra
    MS1_count, MS2_count = 0, 0
    for tic_ in tic_list:
        if tic_['MSlevel'] == 1:
            MS1_count += 1
        else:
            MS2_count += 1
    count_all = len(tic_list)
    print('1. Total count, MS1 count, and MS2 count are:', count_all, MS1_count, MS2_count)
    i = 0
    id_list = []
    for one_tic in tic_list:
        if one_tic['MSlevel'] == 1:
            id = one_tic['ID']
            end = id + 4 if id + 4 < count_all else count_all

            if not [tic_list[id_]['MSlevel'] for id_ in range(id, end)] == [2, 2, 2, 2]:
                id_list.append(id - 1)
    print("2. Output indices that are not followed by 4 MS2 spectra (already decremented by 1):", id_list)

"""
=============================
Data Generation
"""

# Generate file list in the format ‘filename-energy-pos/neg.mzML’ for 26386 -- (only used for the first time), the second time filename is “{}-NCE{}-{}”
def get_filename_list(filepre_list):
    filename_list = []
    ev_list = ['1', '10', '100', '15', '20', '25', '30', '35', '40', '50', '60', '70', '80']
    for filepre in filepre_list:
        for ev in ev_list:
            for tail in ['pos', 'neg']:
                filename_list.append('{}-{}-{}'.format(filepre, ev, tail))
    return filename_list

# Description: Generate MS1 data files
# Input: temp_df - list of compounds provided by Xiao Yin; tic_MS1_df is the list of MS1 spectra from tic_list
# pos_neg indicates positive or negative ionization; MS1_2dict_file is a two-dimensional dictionary for MS1 spectrum data, where the key is CAS and the content is rt and intensity
def make_MS1_data(tic_list, temp_df, tic_MS1_df, pos_neg, MS1_2dict_file):
    if os.path.exists(MS1_2dict_file):
        MS1_2dict = read_pkl(MS1_2dict_file)
    else:
        MS1_2dict = {}
        rt_list = []
        # First, build the dictionary         
        for temp_index in temp_df.index:
            cas = temp_df.loc[temp_index, 'CAS No.']
            MS1_2dict[cas] = {'rt': rt_list, 'intensity': []}
        # Add data to the dictionary
        for tic_MS1_index in tic_MS1_df.index:
            tic_one_dict = tic_list[tic_MS1_index]
            d_df = pd.DataFrame(tic_one_dict['peaks'], columns=['mass', 'intensity'])
            rt_list.append(tic_one_dict['RT'])
            for i in temp_df.index:
                cas = temp_df.loc[i, 'CAS No.']
                # Old version
                # mz = temp_df.loc[i, 'm/z-' + pos_neg]
                # Modified on 240409
                mz = temp_df.loc[i, 'm/z']
                
                filter_df = d_df[d_df['mass'].apply(lambda x: less_5ppm_true(abs(cal_ppm(mz, x))))]
                if filter_df.shape[0] != 0:
                    MS1_2dict[cas]['intensity'].append(filter_df.sum()[1])
                else:
                    MS1_2dict[cas]['intensity'].append(0)
        to_pkl(MS1_2dict, MS1_2dict_file)
    return MS1_2dict

"""
=============================
Data Calculation
Change 5ppm to 10ppm here, modify these 5 numbers
"""

# Calculate ppm
def cal_ppm(expected_value, real_value):
    return (expected_value - real_value) / expected_value * 1000000

# Is it within ±5ppm?
# Description: The first-level calculation uses 5ppm, while the second-level uses 10ppm. To minimize code changes, the less_5ppm_true method is added for first-level 5ppm.
def less_5ppm_true(ppm_value):
    return True if abs(ppm_value) <= 5 else False

def less_5ppm(ppm_value):
    return True if abs(ppm_value) <= 10 else False

def less_10ppm(ppm_value):
    return True if abs(ppm_value) <= 10 else False

# Calculate theoretical mz upper limit
def get_mz_uplimit(expected_value, diff=10):
    return expected_value + diff * expected_value / 1000000

def get_mz_lowlimit(expected_value, diff=10):
    return expected_value - diff * expected_value / 1000000

def get_mz_range(expected_value, diff=10):
    mz_up = get_mz_uplimit(expected_value, diff)
    mz_low = get_mz_lowlimit(expected_value, diff)
    return mz_low, mz_up

# Matrix acceleration
def get_mz_range_matrix(expected_values, diff=10):
    # Convert expected_values to a numpy array
    expected_values = np.array(expected_values)
    
    # Calculate the mz range using matrix operations
    mz_low = expected_values - diff * expected_values / 1000000
    mz_up = expected_values + diff * expected_values / 1000000
    
    # Combine the lower and upper limits into a list of tuples
    mz_range_list = np.vstack((mz_low, mz_up)).T
    return mz_range_list.tolist()

# Check if mz is greater than the limits, return a boolean array, used to replace apply for acceleration
def mz_is_less5(data_Series, theo_mz):
    mz_low, mz_up = get_mz_range(theo_mz)
    return (data_Series.values >= mz_low) & (data_Series.values <= mz_up)

# Check if a specific column in df is within the limits
def between_low_up(data_Series, low_, up_):
    return (data_Series.values >= low_) & (data_Series.values <= up_)

# Check if greater than the upper limit
def mz_bigger_uplimit(expected_value, real_value):
    mz_uplimit = get_mz_uplimit(expected_value)
    return real_value <= mz_uplimit


"""
=============================
Data Extraction
"""

# Description: Get upper and lower limits
# First generate a threshold; if the minimum value < 1%, take 1%; if the minimum value > 1%, take (min/max + 1%)
def get_threshold(data_df):
    min_ratio = (data_df['intensities_list'].min() - 100) / data_df['intensities_list'].max()
    if min_ratio < 0.01:
        return 0.01
    else:
        return min_ratio + 0.01

# Get the range where the maximum intensity is located
def get_max_rt_limit(limit_2list, max_intensity_rt):
    for low_, up_ in limit_2list:
        if (low_ <= max_intensity_rt <= up_):
            return [low_, up_]
            break
            
def get_rt_limit(rt_list, intensities_list, refer_rt=-1, name=""):
    # MS1 spectrum
    mz_RT_df = pd.DataFrame({'rt_list': rt_list, 'intensities_list': intensities_list})

    # If all intensities are zero, it means this raw file did not find this compound
    no0_df = mz_RT_df[mz_RT_df['intensities_list'] != 0]
    if no0_df.shape[0] == 0:
        explain = 'All intensities are 0'
        return [-1, -1], explain
    
    no0_df = no0_df.reset_index(drop=True)
    max_intensity = no0_df['intensities_list'].max()
    max_intensity_rt = no0_df[no0_df['intensities_list'] == max_intensity].iloc[0, 0]

    # Get threshold
    threshold = get_threshold(no0_df)
    threshold_intensity = max_intensity * threshold
    
    # Create a segmented dataframe
    threshold_df = mz_RT_df['intensities_list'] >= threshold_intensity
    threshold_diff = threshold_df.diff()
    threshold_diff.loc[0] = True  # Include head and tail
    threshold_diff.loc[threshold_diff.shape[0] - 1] = True
    mz_threshold_diff_df = mz_RT_df[threshold_diff]
    
    # Segment data and get limits, merge [(1,2),(3,4)] and [(2,3)]
    limit_2list = list(zip(mz_threshold_diff_df.iloc[::2, 0], mz_threshold_diff_df.iloc[1::2, 0]))
    limit_2list.extend(list(zip(mz_threshold_diff_df.iloc[1::2, 0], mz_threshold_diff_df.iloc[2::2, 0])))
    limit_2list = sorted(limit_2list, key=lambda x: x[0])

    lower_limit_rt, upper_limit_rt = [-1, -1]
    explain = ''
    max_rt_low, max_rt_up = get_max_rt_limit(limit_2list, max_intensity_rt)
    for low_, up_ in limit_2list:
        
        # Data within the range
        range_df = mz_RT_df[(mz_RT_df['rt_list'] >= low_) & (mz_RT_df['rt_list'] <= up_)]
        # Check if the range is all below the minimum threshold; skip current range
        if (range_df.iloc[:-1, 1] < threshold_intensity).all():
            continue
        
        # With reference rt        
        if (refer_rt != -1):
            # Peak width exceeds 60:
            if (up_ - low_) > 60:
                # If refer_rt is within ±20 of max intensity, take max ±20; otherwise take refer_rt ±20
                if (max_intensity_rt - 20 <= refer_rt <= max_intensity_rt + 20):
                    lower_limit_rt, upper_limit_rt = max_intensity_rt - 20, max_intensity_rt + 20
                    explain = 'Peak width > 60, refer_rt is within max ±20'
                    break
                else:
                    lower_limit_rt, upper_limit_rt = refer_rt - 20, refer_rt + 20
                    explain = 'Peak width > 60, refer_rt is outside max, only select refer_rt ±20'
                    break
            else:
                # Peak width < 60, find the range that includes both refer_rt and max_intensity_rt; then find the range that contains refer_rt; finally find the range with max intensity ±12
                if (low_ <= refer_rt <= up_) & (low_ <= max_intensity_rt <= up_):
                    lower_limit_rt, upper_limit_rt = low_ - 6, up_ + 6
                    explain = 'Peak width < 60, both refer_rt and max are in the range, range ±6'
                    break
                elif (low_ <= refer_rt <= up_):
                    # If max is within rt ±6, use the range of max intensity ±6
                    if (low_ - 6 <= max_intensity_rt <= up_ + 6):
                        lower_limit_rt, upper_limit_rt = max_rt_low - 6, max_rt_up + 6
                        explain = 'Peak width < 60, only refer_rt in the range, max within ±6, take max range ±6'
                        break
                    # If the current range's maximum intensity is greater than 5e5
                    range_max_intensity = range_df['intensities_list'].max()
                    if range_max_intensity >= 5e5:
                        lower_limit_rt, upper_limit_rt = low_ - 6, up_ + 6
                        explain = 'Peak width < 60, only refer_rt in the range, max not within ±6, range_max > 5e5, take current range ±6'
                        break
                elif (low_ <= max_intensity_rt <= up_):
                    if (low_ - 6 < refer_rt < up_ + 6):
                        lower_limit_rt, upper_limit_rt = low_ - 6, up_ + 6
                        explain = 'Peak width < 60, max filtered, refer_rt is within ±6, range ±6'
                        break
                    
        # If no reference rt, directly look between 6 and 30
        if (refer_rt == -1):
            lower_limit_rt, upper_limit_rt = 6, 30
            explain = 'No refer_rt, directly look for the middle between 6-30'
            break
                
    if ([lower_limit_rt, upper_limit_rt] == [-1, -1]):
        # Additional rule added on 3.23
        if max_intensity >= 1e6:
            lower_limit_rt, upper_limit_rt = max_rt_low - 6, max_rt_up + 6
            explain = 'No range found, max > 1e6, directly take max range ±6'
        else:
            explain = 'Completely not found'
            return [-1, -1], explain

    return [lower_limit_rt, upper_limit_rt], explain

#=================
#=== Optimal Spectrum Judgment ===
#=================

# Description: Get the parent ion peak from the current mass spectrometry data; return [-1,-1] if none
# Input: mass spectrometry data, basic molecular information
# Output: [parent ion mass, intensity]
def get_real_val(intensity_df, theo_mz):
    temp_df = intensity_df[mz_is_less5(intensity_df['Mass'], theo_mz)]
    
    if temp_df.empty:
        return [-1, -1]
    else:
        temp_max = temp_df['Intensity'].max()
        temp_max_df = temp_df[temp_df['Intensity'] == temp_max]
        return [temp_max_df.iloc[0, 0], temp_max_df.iloc[0, 1]]

# Description: Get the base peak. Returns [-1,-1] if data is empty after removing upper limit peaks
# Input: mass spectrometry data, basic molecular information
# Output: [base peak mass, intensity]
def get_basepeak(data_df, theo_mz=""):
    try:
        if theo_mz != "":
            # Calculate the upper limit of the parent ion peak
            mz_uplimit = get_mz_uplimit(theo_mz)
            # Remove peaks exceeding the upper limit
            data_df = data_df[data_df['Mass'].values <= mz_uplimit]
        temp_max = data_df['Intensity'].max()
        temp_df = data_df[data_df['Intensity'] == temp_max]
        return [temp_df.iloc[0, 0], temp_df.iloc[0, 1]]
    except Exception as e:
        print('get_basepeak error:', e)
        return [-1, -1]

# Description: Check if the base peak mass is reasonable
# The difference between the base peak and parent ion peak should be in the range of 4-13; 19-25 is unreasonable. Updated on 3.11 to (4-13,20-25). Updated on 12.13 to [3,13],[21,24]
# Output: [is reasonable, difference]
def basepeak_is_reasonable(molpeak_mass, basepeak_mass):
    mass_dif = molpeak_mass - basepeak_mass
    if (mass_dif < -1) | (mass_dif <= 0.99999) | (3 <= mass_dif <= 13) | (21 <= mass_dif <= 24):
        return [False, mass_dif]
    else:
        return [True, mass_dif]

# Description: Check if the molecular ion peak is greater than the base peak
# Input: Molecular peak intensity, base peak intensity, threshold 5%
# Output: [Is greater than 5%, percentage]
def greater_base_5per(molpeak_intensity, basepeak_intensity, Threshold):
    percent = round(molpeak_intensity / basepeak_intensity, 4) * 100
    if percent >= Threshold:
        return True, percent
    else:
        return False, percent

# Description: Calculate the number of fragment peaks with intensity greater than a specified percentage of the base peak
# Input: Read each mass spectrometry file (two columns, mass and Intensity), threshold for fragment peaks
# Output: DataFrame of peaks greater than the specified threshold
def get_fragmentpeak_df(intensity_df, basepeak_intensity, Threshold, theo_mz):
    fragmentpeak_df = intensity_df[(intensity_df['Intensity'] / basepeak_intensity) * 100 >= Threshold]
    # Remove unreasonable peaks using basepeak_is_reasonable method
    fragmentpeak_df = fragmentpeak_df[fragmentpeak_df['Mass'].apply(lambda x: basepeak_is_reasonable(theo_mz, x)).apply(lambda x: pd.Series(x))[0]]
    return fragmentpeak_df

# Description: Check if the number of fragment peaks is greater than or equal to 1
def fragmentpeak_num_greater1(fragmentpeak_df):
    return True if fragmentpeak_df.shape[0] > 1 else False

# Description: Get the filenames of MS2 files in the input folder
def get_ms2_file(foldername):
    # First check if the folder exists
    if not os.path.exists(foldername):
        return -1
    ms2_filename_list = []
    for i in os.listdir(foldername):
        if '_MS2_' in i:
            ms2_filename_list.append(i)
    if ms2_filename_list == []:
        return -1
    else:
        return ms2_filename_list

# Get the intensity of the base peak in the current MS2 graph, note to remove peaks greater than the molecular ion upper limit
# Input file path
def get_ms2_base(filepath, theo_mz):
    if not os.path.exists(filepath):
        return 0
    try:
        ms2_data = pd.read_csv(filepath, index_col=0)
        # Calculate the upper limit of the molecular ion peak
        mz_uplimit = get_mz_uplimit(theo_mz)
        # Remove peaks that exceed the upper limit
        ms2_data = ms2_data[ms2_data['Mass'] <= mz_uplimit]
    except Exception:
        print(filepath)
        return 0
    return get_basepeak(ms2_data)[1]

"""=============================
Similarity Calculation
"""
from sklearn.metrics.pairwise import cosine_similarity

# Description: Round DataFrame values
def round_df(input_df, rate=10):
    if input_df.shape[0] > 0:
        input_df['Mass'] = round(input_df['Mass'] * rate)
        max_indexes = input_df.groupby(by=['Mass'])['Intensity'].idxmax()
        input_df = input_df.loc[max_indexes]
        input_df['Mass'] = input_df['Mass'].astype(int)
        return input_df
    else:
        input_df['Mass'] = input_df['Mass'].astype(int)
        return input_df

# Description: Layered rounding of DataFrame values
def layer_round_df(input_df):
    input_low_df = input_df[input_df['Mass'] <= 100]
    input_high_df = input_df[input_df['Mass'] > 100]

    input_low_df = round_df(input_low_df, 1000)
    input_high_df = round_df(input_high_df, 100)
    input_high_df = round_df(input_high_df, 10)

    input_df = pd.concat([input_low_df, input_high_df])
    return input_df

# Description: Read data from file with layered rounding
def layer_read_data(filepath, mz, rate=10):
    data_df = pd.read_csv(filepath, usecols=["Mass", "Intensity"])
    # Calculate the upper limit of the molecular ion peak and remove peaks that exceed it
    mz_uplimit = get_mz_uplimit(mz)
    data_df = data_df[data_df['Mass'] <= mz_uplimit]
    
    data_df = layer_round_df(data_df)
    data_df.set_index(['Mass'], inplace=True)
    return data_df

# Description: Read data from file with rounding
def read_data(filepath, mz, rate=100):
    data_df = pd.read_csv(filepath, usecols=["Mass", "Intensity"])
    # Calculate the upper limit of the molecular ion peak and remove peaks that exceed it
    mz_uplimit = get_mz_uplimit(mz)
    data_df = data_df[data_df['Mass'] <= mz_uplimit]
    
    data_df = round_df(data_df, rate)
    data_df.set_index(['Mass'], inplace=True)
    return data_df

# Description: Change data by applying upper limit and rounding
def change_data(data_df, mz):
    # Calculate the upper limit of the molecular ion peak and remove peaks that exceed it
    mz_uplimit = get_mz_uplimit(mz)
    data_df = data_df[data_df['Mass'] <= mz_uplimit]
    
    data_df = layer_round_df(data_df)
    data_df.set_index(['Mass'], inplace=True)
    return data_df

# ====== Similarity Calculation ======
def fu_jaccard_similarity_matrix(X, Y=None):
    # If Y is None, use X for both arguments (i.e., compare all pairs within X)
    if Y is None:
        Y = X

    # Ensure input is numpy arrays
    X = np.atleast_2d(X)
    Y = np.atleast_2d(Y)

    # Compute dot products
    dot_product = np.dot(X, Y.T)
    
    # Compute the sum of squares for each vector
    sum_of_squares_X = np.sum(X ** 2, axis=1)[:, np.newaxis]  # Make column vector
    sum_of_squares_Y = np.sum(Y ** 2, axis=1)  # Make row vector
    
    # Compute the Jaccard similarity matrix
    denominator = sum_of_squares_X + sum_of_squares_Y - dot_product
    jaccard_index = dot_product / denominator
    
    return jaccard_index


import numpy as np
from sklearn.utils import check_array
def entropy_similarity(x, y=None):
    """
    Calculate the entropy similarity between two vectors or matrices.

    Parameters:
    - x: array-like, shape (n_samples, n_features) or (n_features,)
        First input array or vector.
    - y: array-like, shape (n_samples, n_features) or (n_features,), optional
        Second input array or vector. If None, the similarity will be computed
        between the rows of x.

    Returns:
    - Similarity matrix or vector, depending on the shape of the inputs.
    """
    def cal_Ipie(x_one, s_x):
        if s_x >= 3:
            return x_one
        else:
            w = 0.25 + s_x * 0.25
            return np.power(x_one, w)

    def compute_entropy(x):
        x = np.clip(x, 1e-10, None) # 将0都替换为1e-10
        x_sum = np.sum(x)
        x_norm = x / x_sum
        s_x = -np.sum([i * np.log(i) for i in x_norm])
        # print(x_norm, s_x)
        return x_norm, s_x

    def compute_entropy_pie(x, s_x):
        x = np.clip(x, 1e-10, None)
        x_pie = np.array([cal_Ipie(i, s_x) for i in x])
        x_pie_sum = np.sum(x_pie)
        x_pie_norm = x_pie / x_pie_sum
        s_x_pie = -np.sum([i * np.log(i) for i in x_pie_norm])
        return x_pie_norm, s_x_pie

    x = check_array(x, accept_sparse=False, ensure_2d=False)
    if y is not None:
        y = check_array(y, accept_sparse=False, ensure_2d=False)

    if x.ndim == 1:
        # Handle case where x is a vector
        x_norm, s_x = compute_entropy(x)
        if y is None:
            # Calculate similarity with itself
            y_norm, s_y = compute_entropy(x)
            x_pie_norm, s_x_pie = compute_entropy_pie(x, s_x)
            y_pie_norm, s_y_pie = compute_entropy_pie(y, s_y)
            ab_pie_norm = (x_pie_norm + y_pie_norm) / 2
            s_ab_pie = -np.sum([i * np.log(i) for i in ab_pie_norm])
            entropy_sim = 1 - ((2 * s_ab_pie - s_x_pie - s_y_pie) / np.log(4))
            return entropy_sim
        else:
            # Calculate similarity between x and y
            y_norm, s_y = compute_entropy(y)
            x_pie_norm, s_x_pie = compute_entropy_pie(x, s_x)
            y_pie_norm, s_y_pie = compute_entropy_pie(y, s_y)
            ab_pie_norm = (x_pie_norm + y_pie_norm) / 2
            s_ab_pie = -np.sum([i * np.log(i) for i in ab_pie_norm])
            entropy_sim = 1 - ((2 * s_ab_pie - s_x_pie - s_y_pie) / np.log(4))
            return entropy_sim

    elif x.ndim == 2:
        # Handle case where x is a matrix
        n_samples_x = x.shape[0]
        n_samples_y = y.shape[0] if y is not None else n_samples_x

        if y is None:
            similarity_matrix = np.zeros((n_samples_x, n_samples_x))
            for i in range(n_samples_x):
                for j in range(i, n_samples_x):
                    entropy_sim = entropy_similarity(x[i], x[j])
                    similarity_matrix[i, j] = entropy_sim
                    similarity_matrix[j, i] = entropy_sim
            return similarity_matrix
        else:
            similarity_matrix = np.zeros((n_samples_x, n_samples_y))
            for i in range(n_samples_x):
                for j in range(n_samples_y):
                    entropy_sim = entropy_similarity(x[i], y[j])
                    similarity_matrix[i, j] = entropy_sim
            return similarity_matrix
# The following is for the -1 column as an unknown substance #
# Modification 240802: Added entropy similarity calculation
def cal_pos_sim(mz_filter_df, sim_method="cosine"):
    pos_sim_df = mz_filter_df[mz_filter_df[-1] != 0]
    if pos_sim_df.shape[0] == 0:
        return [0 for i in range(pos_sim_df.shape[1] - 1)]
    if sim_method == "cosine":
        pos_matrix = cosine_similarity(pos_sim_df.T)
    elif sim_method == "jaccard":
        pos_matrix = fu_jaccard_similarity_matrix(pos_sim_df.T)
    elif sim_method == "entropy":
        pos_matrix = entropy_similarity(pos_sim_df.T)
    pos_sim_list = list(np.round(pos_matrix[0][1:], 4))
    return pos_sim_list

# Reverse mode calculation
def cal_neg_sim(standard_ms2_df, mz_filter_df, sim_method="cosine"):
    standard_columns_list = list(standard_ms2_df.columns)
    neg_sim_list = []
    for temp_index in standard_columns_list:
        temp_df = mz_filter_df.loc[:, [-1, temp_index]]
        temp_df = temp_df[temp_df[temp_index] != 0]
        # Some MS2 graphs are empty, similarity is set to -1
        if temp_df.shape[0] == 0:
            neg_sim_list.append(-1)
        else:
            if sim_method == "cosine":
                neg_sim_list.append(cosine_similarity(temp_df.T)[0][1])
            elif sim_method == "jaccard":
                neg_sim_list.append(fu_jaccard_similarity_matrix(temp_df.T)[0][1])
            elif sim_method == "entropy":
                neg_sim_list.append(entropy_similarity(temp_df.T)[0][1])
    neg_sim_list = list(np.round(neg_sim_list, 4))
    
    return neg_sim_list

# Bilateral calculation
def cal_both_sim(mz_filter_df, sim_method="cosine"):
    both_sim_df = mz_filter_df
    if sim_method == "cosine":
        both_matrix = cosine_similarity(both_sim_df.T)
    elif sim_method == "jaccard":
        both_matrix = fu_jaccard_similarity_matrix(both_sim_df.T)
    elif sim_method == "entropy":
        both_matrix = entropy_similarity(both_sim_df.T)
    both_sim_list = list(np.round(both_matrix[0][1:], 4))
    return both_sim_list

# Convert to relative abundance
def to_relative_abundance(data_df):
    # Find the maximum value of each column
    max_values = data_df.max()
    # Divide each column by its maximum value
    result = data_df.div(max_values)
    return result

"""=============================
Plotting
"""

# Draw Total Ion Chromatogram (TIC)
def draw_tic(exp):
    tic = exp.calculateTIC()
    retention_times, intensities = tic.get_peaks()
    plt.plot(retention_times, intensities)
    plt.title('TIC')
    plt.xlabel('time (s)')
    plt.ylabel('intensity (cps)')
    plt.show()

# Plotting
# Input: x-axis, y-axis, title, plot type (affects axis labels), x label range, vertical red lines
def plot_pic(x_list, y_list, title, pic_type, label_lim=None, rt_line_list=None):
    title = str(title)
    plt.plot(x_list, y_list)

    plt.title(title)
    if pic_type == 'rt':
        xlabel, ylabel = 'time (s)', 'intensity (cps)'
    elif pic_type == 'mz':
        xlabel, ylabel = 'mass', 'intensity'   
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    if label_lim is not None:
        plt.xlim(*label_lim)
    if rt_line_list is not None:
        for rt_line in rt_line_list:
            plt.axvline(rt_line, color='red')
    plt.show()
