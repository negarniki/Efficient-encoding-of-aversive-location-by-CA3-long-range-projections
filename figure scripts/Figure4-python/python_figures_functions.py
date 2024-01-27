# -*- coding: utf-8 -*-
"""
Functions used in the computation and plots for the paper "Efficient encoding of aversive location by CA3 long-range projections"

Created by Albert Miguel López

We used the decoder codes made available by the Kording Lab: https://github.com/KordingLab/Neural_Decoding/tree/master

"""





#Numpy-related imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import h5py


#Specific imports
from sklearn import decomposition
from sklearn.model_selection import KFold
from Neural_Decoding.decoders import WienerFilterRegression, WienerCascadeRegression, KalmanFilterRegression, SVRegression, NaiveBayesRegression, XGBoostRegression
from sklearn.metrics import r2_score as r2_score_fun
from scipy.ndimage import gaussian_filter
from mpl_toolkits.mplot3d.art3d import Line3DCollection

#Global parameters
import python_figures_parameters as figparams



############# GENERAL FUNCTIONS #############

#Preprocessing functions

def load_place_cell_boolean(mnum, snum, criteria='dombeck'):
    from scipy.io import loadmat
    ''' Loads the place cell arrays for a particular animal and session 
        WARNING: the path is hardcoded!
    ''' 
    
    if criteria == 'dombeck':     
        place_cell_path = figparams.PLACE_CELL_PATH_DOMBECK   
    elif criteria == 'losonczy':
        place_cell_path = figparams.PLACE_CELL_PATH_LOSONCZY   

    place_cell_bool_dataset = loadmat(place_cell_path)['place_cell_bool']
    place_cell_bool = np.squeeze(place_cell_bool_dataset[snum, mnum])
    
    return place_cell_bool




def smoothing(data, std, axis=0):
    '''
    Applies a gaussian filter of given size and std along the specified axis
    '''
    
    sigma = np.zeros(data.ndim)
    sigma[axis] = std
    return gaussian_filter(data, sigma=sigma)



def sum_array_by_chunks(array, bin_size, axis=-1):
    ''' Sums every X elements of array along specified axis, with X being chunk_size
    
        Parameters:
        ----------
        array: numpy array, any dimensions
        chunk_size: number of elements from array that will be summed
        axis: axis of array along which the sum will be performed
        pro
        Returns:
        ------------
        array_out: array of dimensions same as input array, except for "axis" which will be its original size divided by chunk_size
        
        eg: 
            array = [[0 0 1 0 0 1],
                     [1 1 0 1 1 1]]
            axis = -1 (or, equivalently, 1)
            chunk_size = 2
            
            array_out = [[0, 1, 1],
                         [2, 1, 2]]
        
    '''
    
    if 0 > axis:
        axis = axis + array.ndim
    
    shape = array.shape
    
    # Pad the axis "axis" with zeros to make it even with the chunk size
    axis_residual = shape[axis] % bin_size
    if axis_residual > 0:
        pad_width = np.zeros((array.ndim, 2))
        pad_width[axis][1] = bin_size - axis_residual
        array = np.pad(array, pad_width.astype(int), mode='constant',constant_values = 0)

    array_reshaped = array.reshape(shape[:axis] + (-1, bin_size) + shape[axis+1:])
    array_out = array_reshaped.sum(axis=axis+1)
    return array_out


def normalize_data(data, axis=1):
    ''' Elements along given axis have their average subtracted and are divided by their standard deviation '''
    data_mean = np.mean(data, axis=axis, keepdims=True)
    data_std = np.std(data, axis=axis, keepdims=True)
    data_std[data_std==0] = 1 #For cases where there's no variation in the data to avoid dividing by zero
    return (data-data_mean)/data_std, data_mean, data_std







def compute_velocity_and_eliminate_zeros(position, data, pos_max = figparams.MAX_POS):
    ''' Given a 1D array "position" of size "timepoints", return the velocity at each point. "position" is assumed periodic in the range [0, pos_max]
        Additionally eliminate all the points where it is zero, by collapsing the position elements (e.g. 1,3,5,5,5,6 becomes 1,3,5,5,6)
        "Data" is a related 2D matrix of size "features X timepoints", the collapsed points are averaged accordingly so no information is ignored
    '''
    
    v = get_periodic_difference(position[2:], position[:-2], figparams.MAX_POS)/2 #Centered difference
    v0 = get_periodic_difference([position[1]], [position[0]], figparams.MAX_POS) #Forward difference
    vend = get_periodic_difference([position[-1]], [position[-2]], figparams.MAX_POS) #Backward difference
    v = list(v0) + list(v) + list(vend)
    # v = compute_velocity(position, pos_max=figparams.MAX_POS)
    zero_v_bool = np.array(v) < 1e-10
    non_zero_v_bool = np.invert(zero_v_bool)

    if np.any(zero_v_bool): #for ever-increasing velocity this can only happen when consecutive positions remain the same (3 minimum at the center, 2 min at the tails of the array)              
    
        diff_for_zeros = np.diff((zero_v_bool==0).astype(int)) #Will be 1 before 0, -1 at the point where a 0 stops
        zero_vel_interval_starts = np.where(diff_for_zeros==-1)[0]
        zero_vel_interval_ends = np.where(diff_for_zeros==1)[0]+1
        
        starts_num = len(zero_vel_interval_starts); ends_num = len(zero_vel_interval_ends)
        #Is the first interval at the start? Calculate the mean neuronal activity of this interval, assign it to the last element of the interval (which has v != 0)
        if ends_num > 0 and ((starts_num == 0 and ends_num == 1) or zero_vel_interval_starts[0] > zero_vel_interval_ends[0]):
            end = zero_vel_interval_ends[0]
            data[:, end] = np.mean(data[:, :end+1], axis=1)
            zero_vel_interval_ends = zero_vel_interval_ends[1:]
            ends_num = ends_num-1
        
        
        
        #Is the last interval at the end? Calculate the mean neuronal activity of this interval, assign it to the last element of the interval (which has v != 0)
        if starts_num > 0 and ((starts_num == 1 and ends_num == 0) or zero_vel_interval_starts[-1] > zero_vel_interval_ends[-1]):
            start = zero_vel_interval_starts[-1]
            data[:, start] = np.mean(data[:, start:], axis=1)
            zero_vel_interval_starts = zero_vel_interval_starts[:-1]
            starts_num = starts_num-1
    
        #Are there any zero vel intervals in the middle? Average neuronal activities within the zero vel intervals, substitute value at its edges            
        if starts_num > 0:
            for start,end in zip(zero_vel_interval_starts, zero_vel_interval_ends):
                mid = int(np.ceil(start + (end-start)/2))
                data[:, start] = np.mean(data[:, start:mid], axis=1)
                data[:, end] = np.mean(data[:, mid:end+1], axis=1)
            
        position = position[non_zero_v_bool]
        data = data[:, non_zero_v_bool]
        v = list(np.array(v)[non_zero_v_bool])

        
    return position, data, v, non_zero_v_bool

def get_data_from_datadict(data_dict, data_used, distance_limits=None, pos_max=figparams.MAX_POS):
    ''' Convenience function that returns the relevant data from the data_dict
        data_dict: output of "read_CAIM" function
        data_used: spikes, amplitudes, or scaled spikes
        running: if True, only running datapoints are used
        eliminate_zero_v: if True, all the points were v=0 (using centered differences) are averaged out, so that distance increases almost everywhere
            NOTE: this is different from "running". "running" is a boolean from the original dataset which sometimes includes points with no detected distance increase, v=0 is applied on top of this
        
        distance_limits: if not None, must be a two element array with the minimum and maximum distance to analyze
        
    
    '''
    
    spikes = data_dict['spikes']; spikes_binned_normalized = data_dict['spikes_binned_normalized']
    amplitudes = data_dict['amplitudes']; amplitudes_binned_normalized = data_dict['amplitudes_binned_normalized']
    distance = data_dict['distance']
    times = data_dict['times']
    num_neurons = spikes.shape[0]        
    
    if data_used == 'spikes':
        output_data = np.copy(spikes_binned_normalized)
    elif data_used == 'amplitudes':
        output_data = np.copy(amplitudes_binned_normalized)
    elif data_used == 'scaled spikes':
        output_data = np.copy(amplitudes_binned_normalized*spikes_binned_normalized)

    if distance_limits is not None:
        dmin, dmax = distance_limits
        dlim_bool = np.bitwise_and(dmin < distance, distance < dmax)
        output_data = output_data[:, dlim_bool]
        distance = distance[dlim_bool]
        times = times[dlim_bool]

    return output_data, distance, times

def read_CAIM(data_cluster, mouse_num, session_num, trim_data_selection=None):
    '''
    Returns relevant variables from the data cluster
    mouse_num: between 0 and 7, 0-3 are V-D and 4-7 are D-D
    session_num: between 0 and 17
        - 0-2 are baseline trials (B1-B3)
        - 3-6 are airpuff trials at 50cm (T1, T2, Tn-1, Tn)
        - 7-8 are extinction trials (P1-P2)
        - 9-11 are 2nd round of baseline trials (B1-B2)
        - 12-13 are airpuff trials at 100cm (T1, T2)
        - 14 is extinction trial (P1)
        
    trim_data_selection: if not None, it should be a two-element tuple with the start and end of the timepoints to select. Must be floats between 0 and 1.
                        Example: for 5000 timepoints, (0.1,0.2) would select the data from timepint 500 to 1000        
    '''
    
    spikes_ref = data_cluster['CAIM']['S'][mouse_num, session_num] #CAIM for data, S for spikes, first axis is rat number (0-7), second is session (0-17)
    spikes = np.transpose(data_cluster[spikes_ref])
        
    amplitudes = data_cluster[data_cluster['CAIM']['SRaw'][mouse_num, session_num]]
    amplitudes = np.transpose(amplitudes)

    times = data_cluster[data_cluster['CAIM']['behave'][mouse_num,session_num]]['tsscn'][0]
    
    distance = data_cluster[data_cluster['CAIM']['behave'][mouse_num, session_num]]['distance'][0]

    airpuff_ref = data_cluster['CAIM']['AP'][mouse_num, session_num]
    airpuff_bool = data_cluster[airpuff_ref][0]
    
    running_bool = data_cluster[data_cluster['CAIM']['behave'][mouse_num, session_num]]['running'][0]
    
    skaggs_info = data_cluster[data_cluster['CAIM']['cclust'][mouse_num, session_num]][15,:]
    
    num_neurons, num_timepoints = spikes.shape
    

    ######## Mouse 2, session B1 (indexes 1,1) is messed up after t=10030, we deal with it separately
    if (int(mouse_num), int(session_num)) == (1,1):
        thr = 10030
        spikes = spikes[:,:thr]
        
        amplitudes = amplitudes[:, :thr]
        amplitudes = amplitudes[:, :thr]
        times = times[:thr]        
        distance = distance[:thr]
        airpuff_bool = airpuff_bool[:thr]
        running_bool = running_bool[:thr]
        num_neurons, num_timepoints = spikes.shape
        
    #Trim data in a specified way, if specified
    if trim_data_selection is not None:
        start_prop, end_prop = trim_data_selection
        start = int(start_prop * num_timepoints)
        end = int(end_prop * num_timepoints)
        spikes = spikes[:, start:end]
        amplitudes = amplitudes[:, start:end]
        times = times[start:end]
        distance = distance[start:end]
        airpuff_bool = airpuff_bool[start:end]
        running_bool = running_bool[start:end]
        num_timepoints = len(times)
        

    
    
    data_dict = {'spikes':spikes, 'amplitudes': amplitudes, 'times':times, 'distance':distance, 'AP':airpuff_bool, 
                 'running':running_bool, 'num_neurons':num_neurons, 'num_timepoints':num_timepoints,
                 'mouse_num':mouse_num, 'session_num':session_num, 'session_name':figparams.SESSION_NAMES[session_num],
                 'skaggs':skaggs_info}
    
    # data_dict['running'] = running_bool
        
    return data_dict


def read_and_preprocess_data(data_cluster, mouse_num, session_num, gaussian_size, time_bin_size, distance_bin_size, trim_data_selection=None, 
                             only_running=True, eliminate_v_zeros=False, distance_limits=None, pos_max=figparams.MAX_POS, **kwargs):
    ''' Apply read_CAIM function and a bunch of preprocessing, convenience function
        data_cluster: the CAIM .mat file, opened with h5py
        mouse_num: int, should be 0-7
        session_num: int, should be 0-17
        gaussian_size: size (in bin elements) of the gaussian filter
        time_bin_size: size (in bin elements) that are averaged together for analysis
        distance_bin_size: size (* in mm *) of the distance bin. 
            e.g. if "10", all elements 5-15 will become 10, 15-25 to 20, etc.
        trim_data_selection: if not None, it should be a two-element tuple with the start and end of the timepoints to select
        **kwargs are the keyword arguments for the "read_CAIM" function
        
        
        
        '''

    data_dict = read_CAIM(data_cluster, mouse_num, session_num, trim_data_selection=trim_data_selection, **kwargs)
    
    #Spikes
    spikes = data_dict['spikes']
    output_data = [['spikes', spikes]]
    
    if 'amplitudes' in data_dict:
        amplitudes = data_dict['amplitudes']
        output_data.append(['amplitudes', amplitudes])
        
    #Spike pre-processing (smoothing, binning)
    for idx, [data_name, data] in enumerate(output_data):
        data = smoothing(data.astype(float), gaussian_size, axis=1)
        data = sum_array_by_chunks(data, time_bin_size, 1)/time_bin_size # num neurons X num time bins
        output_data[idx][1] = data
    
    #Times
    times = data_dict['times']
    dt = np.average(times[1:] - times[:-1]) #Time interval in miliseconds
    
    #Position
    distance = data_dict['distance']
    distance = distance % 1500 #Eliminates negative distances, which I've seen at least in animal 5 session 2 for some reason
    distance = sum_array_by_chunks(distance, time_bin_size, 0)/time_bin_size
    distance = (distance // distance_bin_size) * distance_bin_size
        
    if 'running' in data_dict:
        running_bool = sum_array_by_chunks(data_dict['running'].astype(float), time_bin_size, 0)/time_bin_size
        running_bool = np.around(running_bool, decimals=0).astype(bool)
        data_dict['running'] = running_bool
        
    if only_running == True:
        for idx, [_, data] in enumerate(output_data):
            data_cut = data[:, running_bool]
            output_data[idx][1] = data_cut
        distance = distance[running_bool]
        times = times[running_bool]
        
    if eliminate_v_zeros == True:
        distance_original = np.copy(distance)
        for idx, [_, data] in enumerate(output_data):
            distance, data_cut, _, _ = compute_velocity_and_eliminate_zeros(distance_original, data, pos_max=figparams.MAX_POS)
            output_data[idx][1] = data_cut
            
        times = np.arange(len(distance)) * dt

    #After all other pre-processing is done, normalize
    for data_name, data in output_data:
        data_normalized, data_mean, data_std = normalize_data(data, axis=1)
        data_dict[data_name + '_binned_normalized'] = data_normalized*1
        data_dict[data_name + '_mean'] = data_mean
        data_dict[data_name + '_std'] = data_std
        # data_dict[data_name + '_binned_normalized'] = data ##!"·"!·## Uncomment for the "old" normalizing behavior ##!"·"!·##
        
    #Position-related values    
    dist_diff = np.diff(distance)
    overround = [0] + list(np.where(dist_diff < -1000)[0]+1)
    num_trials = len(overround)-1
    distance_by_trial = [distance[overround[i] : overround[i+1]] for i in range(num_trials)]           
        
    data_dict['dt'] = dt
    data_dict['distance'] = distance    
    data_dict['overround'] = overround
    data_dict['num_trials'] = num_trials
    data_dict['distance_by_trial'] = distance_by_trial
    
    return data_dict


#PCA related functions

def save_pca_data(mouse_list, session_list):
    ''' Convenience function to pre-compute PCA information '''
    
    ## Get PCA weights ##
    components_dict = {} #(mnum, snum) = (pca_dims, num_neurons)
    variance_explained_dict = {} #(mnum, snum) = (pca_dims, num_neurons)
    PCA_dict = {} #(mnum, snum) = (position, pca_data) with size(position) = timepoints and size(pca_data) = (pca_dims X timepoints)
    input_data_dict = {} #(mnum, snum) = (position, preprocessed_spikes) with size(position) = timepoints and size(preprocessed_spikes) = (num_neurons X timepoints)

    time_bin_size = 1  # Number of elements to average over, each dt should be ~65ms
    distance_bin_size = 1  # mm, track is 1500mm, data is in mm
    gaussian_size = 25  # Why not
    data_used = 'amplitudes'
    running = True
    eliminate_v_zeros = False
    
    with h5py.File(figparams.FAT_CLUSTER_PATH, 'r') as fat_cluster:
        
        for mnum in mouse_list:        
            
            for snum in session_list:
                
                
                #Get data
                data_dict = read_and_preprocess_data(fat_cluster, mnum, snum, gaussian_size, time_bin_size, distance_bin_size, 
                                                        only_running=running, eliminate_v_zeros=eliminate_v_zeros, pos_max=figparams.MAX_POS)
                pca_input_data, position, times = get_data_from_datadict(data_dict, data_used)    
                num_neurons = pca_input_data.shape[0]
                
                if eliminate_v_zeros == True:
                    position, pca_input_data, _, non_zero_bool = compute_velocity_and_eliminate_zeros(position, pca_input_data, pos_max = figparams.MAX_POS)
                
                
            
                #PCA
                pca = decomposition.PCA(n_components=num_neurons)
                pca.fit(pca_input_data.T)
                
                #PCA over time
                pca_data = project_spikes_PCA(pca_input_data, pca_instance = pca, num_components = num_neurons)
                
                # if eliminate_v_zeros == True:
                #     position, pca_data, _, non_zero_bool = compute_velocity_and_eliminate_zeros(position, pca_data, pos_max = figparams.MAX_POS)
                # pca_input_data = pca_input_data[:, non_zero_bool]

                input_data_dict[mnum, snum] = (position, pca_input_data)            
                PCA_dict[mnum, snum] = (position, pca_data)
                
                #Components
                components = pca.components_
                components_dict[mnum, snum] = components
                
                #Variance explained
                variance_explained = pca.explained_variance_ratio_
                variance_explained_dict[mnum, snum] = variance_explained
                
        np.save(figparams.OUTPUT_PATH + "input_data_dict.npy", input_data_dict)
        np.save(figparams.OUTPUT_PATH + "pca_components_dict.npy", components_dict)
        np.save(figparams.OUTPUT_PATH + "variance_explained_dict.npy", variance_explained_dict)
        np.save(figparams.OUTPUT_PATH + "PCA_dict.npy", PCA_dict)

def dimensions_to_explain_variance(variance_explained, variance_to_explain):
    ''' 
        Returns number of dimensions needed to explain a given pca variance.
    
        variance_explained should be the output of pca_instance.explained_variance_ratio_ (from sklearn)
        variance_to_explain must be a float between 0 and 1, representing % of variance to explain
    '''
    variance_explained_cum = np.cumsum(variance_explained)
    dimensions_for_x = np.argmax(variance_explained_cum > variance_to_explain) + 1
    return dimensions_for_x

def project_spikes_PCA(pca_input_data, pca_instance = None, num_components = None, return_pca_instance = False):
    ''' Projects spikes to PCA space.
        pca_input_data: spikes, shape "num features (e.g. neurons)" X "num samples (e.g. timepoints)"
        pca_instance: if None, a pca instance will be created and trained. If not None, it must be a trained instace of sklearn's PCA. 
        num_components: how many pca dimensions to take.
            *if None, take all
            *if int, take that many
            *if 'X%', where "X" is a scalar, take dimensions that explain X% of the variance
            *if 'X', where "X" is a float between 0 and 1, interpret it as taking X*100% of the variance
        
            
    '''
    num_features = pca_input_data.shape[0]
    if pca_instance is None:
        pca_instance = decomposition.PCA(n_components=num_features)
        pca_instance.fit(pca_input_data.T) 
            
    if num_components is None:
        num_components = num_features
        
    elif type(num_components) in [int, float, np.int32]:
        num_components = int(num_components)
        
    elif type(num_components) in [str]:
        variance_explained = pca_instance.explained_variance_ratio_        
        if '%' in num_components: #Dimensions to explain 'X%' of the variance
            variance_to_explain = float(num_components[:num_components.find('%')])
            num_components = dimensions_to_explain_variance(variance_explained, variance_to_explain/100)
        else: #num_components is also variance to explain, but as a float between 0 and 1
            variance_to_explain = float(num_components)        
            num_components = dimensions_to_explain_variance(variance_explained, variance_to_explain)

        # variance_explained_cum = np.cumsum(variance_explained)
        # dimensions_for_x = np.argmax(variance_explained_cum > variance_to_explain) + 1

        
    else:
        raise(TypeError('Unknown data type for "num_components"'))
        
    transform_m = pca_instance.components_[:num_components, :]
    spikes_projected = transform_m @ pca_input_data
    if return_pca_instance == False:
        return spikes_projected
    else:
        return spikes_projected, pca_instance

def compute_pca(pca_input_data, num_components=None):
    num_neurons = pca_input_data.shape[0]
    
    pca = decomposition.PCA(n_components=num_neurons)
    pca.fit(pca_input_data.T)
    pca_data = project_spikes_PCA(pca_input_data, pca_instance = pca, num_components = num_components)
    return pca, pca_data


#### POSITION PREDICTION FUNCTIONS

def from_linear_to_circular_position(position, pmin, pmax):
    
    ''' Takes 1D position data and transforms to circular.
        pmin and pmax are the minimum and maximum position, respectively.
        Returns data in the form "samples" X 2
    '''
    
    angle_d = (position - pmin)/(pmax-pmin) * 2 * np.pi - np.pi
    sin_d = np.sin(angle_d)
    cos_d = np.cos(angle_d)
    angle_pos = np.vstack((sin_d, cos_d)).T
    
    return angle_pos

def initiate_predictor(predictor_name): #Initiates predictor object
    if predictor_name == 'Wiener':              
        predictor = WienerFilterRegression()
        
    elif predictor_name == 'Wiener Cascade':
        predictor = WienerCascadeRegression(degree=4)
        
    elif predictor_name == 'Kalman':
        predictor = KalmanFilterRegression(C=1)
        
    elif predictor_name == 'SVR':
        predictor = SVRegression(C=3, max_iter=-1) #-1 means no iteration limit
        
    elif predictor_name == 'Naive Bayes':
        predictor = NaiveBayesRegression(encoding_model = 'quadratic', res=100)
        
    elif predictor_name == 'XGBoost':
        predictor = XGBoostRegression(max_depth=3, num_round=300, eta=0.3, gpu=-1)
        
        
    else:
        raise NameError('Predictor name for position predictor is not allowed')
    return predictor


def get_prediction(predictor, X, Y): #Get function to handle Kalman filter's different function call
    ''' Convenience function that handles the Kalman filter case differently, X is data set to predict, Y is the actual labels.
        Note that Y is only ever called in the Kalman prediction but only the array's shape is used
        "predictor" must be scikit type trained object (or from "decoders.py")
        X of shape "samples X features"
        Y of shape "samples"
        '''
        
        
    #Predict
    if type(predictor) in [type(KalmanFilterRegression()), type(NaiveBayesRegression())]:
        pred = predictor.predict(X, Y)
    else:
        pred = predictor.predict(X)
    return pred


def get_periodic_difference(prediction, objective, period, get_sign = False):
    ''' Computes the difference between a prediction and its objective but taking periodicity into account
        Assumes they are all 1D array of same size
        If "get_sign" is True, the sign of the difference is kept.
            Example: for period=1500, 20-50 = -30; 1250-50 = -300; 50-20 = 30, 50-1250 = 300
            Explanation: the shortest path from 50 to 1250 is to move 300 units left. The shortest path from 1250 to 50 is 300 units to the right.
    '''
    
    if type(prediction) in [int, float, np.float64]:
        prediction = [prediction]
    if type(objective) in [int, float, np.float64]:
        objective = [objective]

    prediction = np.array(prediction); objective = np.array(objective)
    
    assert prediction.ndim==1 and objective.ndim==1
    extended_objective = np.stack((objective-period, objective, objective+period))
    raw_diff = extended_objective - prediction
    extended_diff_arg= np.argmin(np.abs(raw_diff), axis=0)
    extended_diff = (raw_diff)[extended_diff_arg, range(len(prediction))]
    if get_sign == False:
        extended_diff = np.abs(extended_diff)
    return extended_diff

def get_sse(prediction, objective, period=None):
    ''' Computes the average squared sum of error. prediction and objective are assumed to be 1d and same size '''
    if period is None:
        diff = prediction-objective
        
    else:
        diff = get_periodic_difference(prediction, objective, period)
    
    sse = np.sum((diff)**2)/(prediction.size)
    # sse = np.trapz(np.abs(diff))
    return sse

def get_error_dict(position, position_pred, period=1500, crosscorr_align = False, crosscorr_error='sse'):
    ''' Returns error dict for objective "distance" and predicted "distance pred" '''

    sse = get_sse(position_pred, position, period=period)
    sse = np.sqrt(sse)/10 #From mm to cm
    best_delay = 0
    
    if period is not None:
        diff = get_periodic_difference(position_pred, position, period=period)/10 #The 10 is the conversion to cm!
    else:
        diff = (position_pred - position)/10
    diff_std = np.std(diff)
    diff_std = diff_std/2 #Done to match Martin Pofahl's code
    diff_avg = np.mean(np.abs(diff))
    
    r2 = r2_score_fun(position, position_pred)
    
    error_dict = {'sse':sse, 'diff':diff, 'diff_std':diff_std, 'diff_avg':diff_avg, 'r2':r2,
                  'best_delay':best_delay}
            

    return error_dict, position_pred
    
def predict_position_CV(data, position, n_splits=5, shuffle=False, periodic=True, pmin=0, pmax=1500,
                        predictor_name='Wiener', predictor_default=None, return_error='sse'):
    ''' 
    1D position data is converted to a 2D circle so the quantity to predict is periodic
    
    data: shape "features" X "samples"
    position: 1D array of size "samples" with real values to predict
    n_splits: number of CV folds. If 1 or less no cross-validation is performed (so no test sets)
    shuffle: if True, data is shuffled [DOESN'T WORK AS OF NOW]
    periodic: if True, data is converted                                        
    dmin, dmax: min and max position to do the periodicity
    predictor name: 'Wiener', 'Kalman', 'SVR'
    predictor_default: if None, a type of predictor specified by "predictor_name" is trained on the data. If one is given, that is used instead (we assume it has been trained)
    return_error: must be a error name. 
    '''
    
    
    #Convert to periodic 2D circle
    X = data.T
    if periodic:    
        Y = from_linear_to_circular_position(position, pmin, pmax) # num samples X num features (which are two)
    else:
        Y = np.vstack(position)

  
    #Predict position
    if predictor_default is None:
        #If no previous predictor is given, train it from scratch
        if n_splits > 1:
            kf = KFold(n_splits=n_splits, shuffle=shuffle)
            Y_pred = np.zeros(Y.shape)
            
            for train_index, test_index in kf.split(X):
                # Y[train_index]
                Y[test_index]
                X_train, X_test = X[train_index], X[test_index]
                Y_train, Y_test = Y[train_index], Y[test_index]
                
                predictor = initiate_predictor(predictor_name)
                predictor.fit(X_train, Y_train)
                Y_pred[test_index] = get_prediction(predictor, X_test, Y_test)
                
        else: #No validation set
            predictor = initiate_predictor(predictor_name)
            predictor.fit(X,Y)
            Y_pred = get_prediction(predictor, X, Y)

            
    else:
        #If a trained predictor is given, test
        predictor = predictor_default
        Y_pred = get_prediction(predictor, X, Y)


    if periodic:
        angle_d_pred = np.arctan2(Y_pred[:,0], Y_pred[:,1])
        position_pred = ((angle_d_pred+np.pi)/(2*np.pi)) * (pmax-pmin) + pmin
        period = pmax
    else:
        period = None
        position_pred = np.squeeze(Y_pred)
    
    error_dict, position_pred = get_error_dict(position, position_pred, period=period)
    error = error_dict[return_error]
    return position_pred, error, predictor


#Plot functions

def generate_title_sessions_and_mice(session_list, mouse_list):
    '''
        Given a list of sessions and mice, returns a title that describes them.
        Takes into account conditions, e.g. if the sessions are 0,1,2 then it calls them baseline
    '''
    
    both_exist = True
    
    if session_list is not None and len(session_list) > 0:
        sessions_title = 'sessions %s' %session_list
        
        if list(session_list)  == [3,4,5,6]:
            label1 = 'Airpuff'
            sessions_title = 'aversive sessions'
        elif list(session_list)  == [0,1,2]:
            sessions_title = 'baseline sessions'
        elif list(session_list) == [7,8]:
            sessions_title = 'probe sessions'
        elif list(session_list) == [0,1,2,3,4,5,6,7,8]:
            sessions_title = 'all sessions'
        elif list(session_list)  == [9,10]:
            sessions_title = 'extinction sessions'
        elif list(session_list)  == [11,12]:
            sessions_title = '2nd baseline sessions'
        elif list(session_list) == [13,14]:
            sessions_title = '2nd aversive sessions'
        elif list(session_list) == [15,16]:
            sessions_title = '2nd probe sessions'        
        elif list(session_list) == [9,10,11,12,13,14,15,16]:
            sessions_title = 'all 2nd round sessions'
        else:
            sessions_title = ' '.join([figparams.SESSION_NAMES[s] for s in session_list])
    else:
        sessions_title = ''
        both_exist = False
        
    if mouse_list is not None and len(mouse_list) > 0:
        
        mouse_title = 'mice %s' %str(mouse_list)[1:-1]
        
        if list(mouse_list) == [0,1,2,3]:
            mouse_title = 'i-D mice'
        elif list(mouse_list) == [4,5,6,7]:
            mouse_title = 'D-D mice'
        elif list(mouse_list) == [0,1,2,3,4,5,6,7]:
            mouse_title = 'all mice'
        else:
            mouse_title = ' '.join(['M%d'%m for m in mouse_list])
    else:
        mouse_title = ''
        both_exist = False
    
    if both_exist == True:
        final_title = sessions_title + ', ' + mouse_title
    else:
        final_title = sessions_title + mouse_title
    
    return final_title

def plot_colored_line_3d(x, y, z, color_scalar, cmap_name='viridis', fig=1, ax=None, lw=2, cbar=True, 
                         xlim=None, ylim=None, zlim=None, color_norm_limits = None):
    ''' Efficiently plots a 2d colored line using LineCollection
        x: dimension N
        y: dimension N
        z: dimension N
        color_scalar: dimension N, each element is a scalar that defines its color
        cmap: string, indicates matplotlib cmap to use for the colorbar
        figure: which figure. Can be int for its number or a figure instance from matplotlib
        ax: if not None, takes precedence over "figure"
        lw: linewidth, thickness of drawing       
        color_norm: if None, the color will use the min and max of color scalar as limits.
                    If two-element list of the sort [min, max], it will use those instead
        
        Outputs the figure, axis, and colorbar objects
        
    '''
    
    cmap = plt.get_cmap(cmap_name)
    color = cmap(color_scalar)
    
    # We are going to create a line for every pair of points, so we reshape our data into a
    # N x line_length x dimensions array. Line length is 2 and dimension 3 (x, y, z).
    #   e.g. element [0,:,:] is [[x0,x1],[y0,y1]]
    points = np.array([x,y,z]).T
    segments = np.stack((points[:-1], points[1:]), axis=-1).transpose(0,2,1)


    if ax is not None: #Axis is given, get figure
        plt.sca(ax)
        fig = plt.gcf()
        
    elif type(fig) is int: #Only fig num is given, create fig and axis
        fig = plt.figure(fig)
        ax = fig.add_subplot(111, projection='3d')
        
    else: #Figure is given, get axis
        # fig = plt.figure(fig)
        ax = fig.gca()
        
    if color_norm_limits is None:
        norm = plt.Normalize(np.min(color_scalar), np.max(color_scalar))
    else:
        norm = plt.Normalize(color_norm_limits[0], color_norm_limits[1])
    lc = Line3DCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(color_scalar)
    lc.set_linewidth(lw)
    line = ax.add_collection(lc)
    
    if cbar:
        cbar = fig.colorbar(line)
        
    if xlim is None:
        ax.set_xlim([np.min(x), np.max(x)])
    else:
        ax.set_xlim(xlim)
        
    if ylim is None:
        ax.set_ylim([np.min(y), np.max(y)])
    else:
        ax.set_ylim(ylim)
        
    if zlim is None:
        ax.set_zlim([np.min(z), np.max(z)])
    else:
        ax.set_zlim(zlim)
    
    return fig, ax, cbar

    
def add_distance_cbar(fig, cmap, vmin = 0, vmax = 1500, fs=15, cbar_label = 'Position (mm)'):
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), orientation='vertical', fraction=0.04)
    cbar.ax.set_title(cbar_label, fontsize=fs)
    cbar.ax.tick_params(axis='both', which='major', labelsize=fs-3)

def plot_pca_with_position(pca, position, ax=None, max_pos = 1500, cmap_name = figparams.PCA_CMAP, fs=15, scatter=True, cbar=False, cbar_label='Position (mm)',
                           alpha=1, angle=50, angle_azim=None, axis = 'off', show_axis_labels=True, axis_label=None, ms = 10, lw=3):
    '''
    Plots the first three dimensions of a PCA.
    
    Parameters
    ----------
    pca : array of size "num features" X "num samples"
    position : array of size "num samples"
        Indicates the position of the animal per sample.
    ax : TYPE, optional
        DESCRIPTION. The default is None.
    max_pos : int, optional
        Maximum value of the position value. The default is 1500.
    cmap_name : string, optional
        Name of the colormap to use. The default is 'hsv'.
    fs : int, optional
        Fontsize. The default is 15.
    scatter : bool, optional
        If True, PCA is plotted using the scatter function. Otherwise lines are drawn. The default is True.
    cbar : bool, optional
        if True, a colorbar is plotted
    alpha : float, optional
        Transparency of the plot. The default is 1.
    angle : int, optional
        Main 3d angle. The default is 50.
    angle_azim : int, optional
        Azimuthal angle. The default is None.
    axis : string, optional
        'on' or 'off'. The default is 'off'.
    show_axis_labels : bool, optional
        Whether or not to label the axis. The default is True.
    axis_label : str, optional
        If None, each axis is labelled as a PCA axis. Otherwise this label is used. The default is None.
    ms : int, optional
        Marker size for the scatter plots. The default is 10.
    lw : int, optional
        Line width for the average PCA plot. The default is 3.

    Returns
    -------
    ax : matplotlib axis object. If None, a new one is created.

    '''
    
    if ax is None:
        ax = plt.subplot(1,1,1,projection = '3d')
    
    fig = plt.gcf()        
    cmap = plt.get_cmap(cmap_name)        
    xlim = ax.get_xlim(); ylim = ax.get_ylim(); zlim = ax.get_zlim()
    
    x, y, z = pca[:3]
    if scatter==True:
        ax.scatter(x, y, z, color=cmap(position/max_pos), s=ms, marker='o', alpha=alpha)        
    else:
        fig, ax, _ = plot_colored_line_3d(x, y, z, position, cmap_name=cmap_name, ax=ax, lw=lw, cbar=False, 
                                          xlim=xlim, ylim=ylim, zlim=zlim, color_norm_limits = [0,max_pos])
    
    if axis_label is None:
        axis_label = 'PCA'

        
    if show_axis_labels == True:    
        ax.set_xlabel(axis_label + ' D1', fontsize=fs)
        ax.set_ylabel(axis_label + ' D2', fontsize=fs)
        ax.set_zlabel(axis_label + ' D3', fontsize=fs) 
    ax.view_init(elev=angle, azim=angle_azim)
    if cbar == True:
        add_distance_cbar(fig, cmap, vmin = 0, vmax = max_pos, fs=fs, cbar_label=cbar_label)
        
    ax.set_xlim([np.minimum(np.min(x), xlim[0]), np.maximum(np.max(x), xlim[1])])
    ax.set_ylim([np.minimum(np.min(y), ylim[0]), np.maximum(np.max(y), ylim[1])])
    ax.set_zlim([np.minimum(np.min(z), zlim[0]), np.maximum(np.max(z), zlim[1])])
    
    if axis == 'off':
        ax.set_axis_off()
        
    elif axis == 'on':
        ax.grid(False)
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) 
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) 
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) 

    return ax

def get_significance_label(pval, thresholds = [0.001, 0.0001], asterisk=False, ns=False):
    
    thresholds = np.sort(thresholds)[::-1]
    signum = np.sum(pval < np.array(thresholds)) #How many thresholds is pval below?
    
    if signum == 0:
        if ns == True:
            return 'ns'
        else:
            return ''
    else:
        if asterisk == False:
            return "p<" + ("%f" %thresholds[signum-1]).rstrip('0').rstrip('.')
        else:
            return "*" * signum

def plot_pca_and_average(pca, position, fig_num, figsize=(6,6), dims = None, max_pos=figparams.MAX_POS, pca_distance_bin_size=20, **pca_kwargs):
    '''
    Plots pca values colored by position vector, next to its trial average

    Parameters
    ----------
    pca : array of size "num features" X "num samples"
    position : array of size "num samples"
        Indicates the position of the animal per sample.
    figsize : tuple
    fig_num : int or matplotlib fig object
    dims : None or list, optional ----------------- TO IMPLEMENT ------------------
        If None, first three dimensions are taken. If list, the indicated dimensions are plotted instead
    max_pos : int, optional
        Maximum value of the position value. The default is 1500.
    pca_distance_bin_size : int, optional
        Bin size (in terms of the position values) for the average PCA plot. The default is 20.
    pca_kwargs: keyword arguments for the "plot_pca_with_position" function

    Returns
    -------
    fig : matplotlib fig object
    fig_num : next empty figure number

    '''
    if type(dims) is list:
        assert len(dims) == 3
        pca = pca[dims]
    
    if type(fig_num) is int:
        figsize = (figsize[0], figsize[1]*2)
        fig = plt.figure(num=fig_num, figsize=figsize); fig_num += 1
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        
    else: #Figure already exists, must have two axes
        fig = fig_num
        print(fig.number)
        fig_num = fig.number
        axs = fig.axes
        ax1 = axs[0]
        ax2 = axs[1]
        
    pca_kwargs['max_pos'] = max_pos
    
    pca_kwargs_nonaverage = copy_dict(pca_kwargs)
    pca_kwargs_nonaverage['cbar'] = False
    plot_pca_with_position(pca, position, ax1, **pca_kwargs_nonaverage)
        
    position_unique, pca_average, pca_std = compute_average_data_by_position(pca, position, position_bin_size=pca_distance_bin_size, max_pos=max_pos)
    
    pca_kwargs_average = copy_dict(pca_kwargs)
    pca_kwargs_average['scatter'] = False
    if 'cbar' not in pca_kwargs:
        pca_kwargs_average['cbar'] = True
    plot_pca_with_position(pca_average, position_unique, ax=ax2, **pca_kwargs_average)


    fig.tight_layout()
    
    return fig, fig_num


#Overall functions

def get_round_endtimes(position, diff_thresh = 1000):
    ''' Get round endtimes. Assumes position goes from its maximum value to 0 after a round end.
        Accounts for reverse movement 
            e.g. 1495, 1499, 3, 1497, 6, 10... counts as 1 single round change
        position is 1-D array
        diff_thresh is how much the position value must change in a single step to consider it a round end.
            Recommend to take 1/2 or 2/3 of maximum position
    '''
    pos_diff = np.diff(position)
    overround = list(np.where(pos_diff < -diff_thresh)[0]+1)
    any_round_back_jumps = np.any(pos_diff > diff_thresh) #E.g. 0, 3, 5, 1448, 1449...

    if any_round_back_jumps:
        overround_back = list(np.where(pos_diff > diff_thresh)[0]+1)
        num_timepoints = len(position)
        #Round counter jumps +1 when going from 1500 to 0, and jumps -1 when going 0 to 1500
        round_counter = np.zeros(num_timepoints, dtype = int)
        for round_end in overround:
            round_counter[round_end:] += 1
        for jump_back in overround_back:
            round_counter[jump_back:] -= 1
        round_counter = list(round_counter)
        num_rounds = np.max(round_counter)
        overround = [round_counter.index(round_idx) for round_idx in range(1,num_rounds+1)]

    return overround

def compute_average_data_by_position(data, position, position_bin_size=None, max_pos=figparams.MAX_POS):
    ''' Given a dataset, calculate its average value for each observed position (using the given bin size)
        Input:
            data: matrix of size "pca dim" X "timepoints"
            position: array of size "timepoints"
            position_bin_size: if None the values from "position" are used to average. If "int", the position values are approximated to the nearest mcm with the bin size
            max_pos: maximum value of the position
            
        Returns:
            position_bins: array of size "num of position bins", contains the position values used to average the data
                e.g.: if bin_size = 20 and the first element is 0, the first data_average element will contain the average of the PCA of all the times the position was between 0 and 20
            data_average: array of size "pca dim" X "num of unique positions", contains the average PCA at each corresponding position bin
        
    '''
    num_dimensions, num_timepoints = data.shape
    position = position % (max_pos+1)
    
    if position_bin_size is not None:
        position = (position // position_bin_size) * position_bin_size
    position_bins = np.unique(position)
    num_bins = len(position_bins)
        
    overround = get_round_endtimes(position, diff_thresh=max_pos * 0.66)
    num_rounds = len(overround)
    
    
    data_average = np.zeros((num_dimensions, num_bins))
    data_std = np.zeros((num_dimensions, num_bins))
    for idx, d in enumerate(position_bins):
        data_filtered = data[:,position==d]
        data_average[:,idx] = np.mean(data_filtered, axis=1)
        data_std[:, idx] = np.std(data_filtered, axis=1)/np.sqrt(num_rounds)
    return position_bins, data_average, data_std

def copy_dict(d):
    ''' Returns independent copy of a dictionary d '''
    return {k:v for k,v in d.items()}

def do_pca_and_predict(position, pca_input_data, pos_max = figparams.MAX_POS, num_components=None, cv_folds=5, predictor_name='Wiener', error_type='sse'):
    
    pca, pca_data = compute_pca(pca_input_data, num_components = num_components)
    
    position_pred, error, _ = predict_position_CV(pca_data, position, n_splits=cv_folds, shuffle=False,
                                                              periodic=True, pmin=0, pmax=pos_max, predictor_name=predictor_name, predictor_default=None,
                                                              return_error=error_type)

    
    return pca_data, position_pred, error


if __name__ == '__main__':
    
    pass