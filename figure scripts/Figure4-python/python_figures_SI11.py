# -*- coding: utf-8 -*-
"""

Functions to create the Supporting Information Figure 11 for the paper "Efficient encoding of aversive location by CA3 long-range projections"

Created by Albert Miguel López

"""




#Global project variables

import time
import csv
from itertools import combinations

import scipy
import scipy.io
from scipy.special import comb

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import h5py

# import processing_functions as pf
import python_figures_functions as figfuns
import python_figures_parameters as figparams
from python_figures_parameters import RESULTS_PATH, MAX_POS


from functools import partial


            
def main():
    
    ''' 
    
        - Running times -
            
        save_pca_data: around 10s
        
        compute_and_plot_SI11_figures:
            Using Wiener: 20 seconds
    
    
    '''
    
    # # #### RUN THIS FIRST TO PROCESS RAW DATA INTO PCA ####
    figfuns.save_pca_data(np.arange(8), np.arange(9))
    # # #### RUN THIS FIRST TO PROCESS RAW DATA INTO PCA ####
    
    ### Plot all figure S11 subfigures ###
    compute_and_plot_SI11_figures()


    ### Uncomment to compute and plot figures independently ###
    # figSI11_a_b()    
    # figSI11_c()
    # figSI11_d_e()
    # figSI11_d_and_e_to_excel()
    
    ### Uncomment to just plot the figures (assumes you have done the computations before) ###
    # plot_precomputed_SI11_subfigures()
    
    
def compute_and_plot_SI11_figures():
    
    ''' Using all the predictors takes around 30 minutes for place cell figure, 8 minutes for the fig4h
    
        Wiener takes around 40s for everything
    
        For XGBoost: around 1h and 30 minutes

            
    '''
    
    # error_type = ['sse']
    error_type = 'diff_std'

    # predictor_list = ['Wiener', 'Wiener Cascade', 'Kalman', 'SVR', 'XGBoost']
    predictor_list = ['Wiener']    
    # predictor_list = ['XGBoost'] 
    
    
    # # #### RUN THIS FIRST TO PROCESS RAW DATA INTO PCA (can uncomment later) ####
    figfuns.save_pca_data(np.arange(8), np.arange(9))
    # # #### RUN THIS FIRST TO PROCESS RAW DATA INTO PCA (can uncomment later) ####

    
    for predictor in predictor_list:
        #SI figures
        figSI11_a_b(predictor, error_type)    
        figSI11_c(predictor, error_type)
        figSI11_d_e(predictor, error_type)
        figSI11_d_and_e_to_excel(predictor, error_type)        
        
        
def plot_precomputed_SI11_subfigures():
    
    ''' 
    Assuming that the "SI11" functions have been run, plot the data.
    Convenience function to play with the plots without having to recalculate everything
    
    '''
    
    # error_type = ['sse']
    error_type = 'diff_std'

    # predictor_list = ['Wiener', 'Wiener Cascade', 'Kalman', 'SVR', 'XGBoost']
    predictor_list = ['Wiener']    
    # predictor_list = ['XGBoost'] 
    
    
    # # #### RUN THIS FIRST TO PROCESS RAW DATA INTO PCA (can uncomment later) ####
    # figfuns.save_pca_data(np.arange(8), np.arange(9))
    # # #### RUN THIS FIRST TO PROCESS RAW DATA INTO PCA (can uncomment later) ####

    
    for predictor in predictor_list:
        #SI figures
        figSI11_a_b(predictor, error_type)    
        figSI11_c_plots(predictor, error_type)
        figSI11_d_e_plots(predictor, error_type)

def figSI11_a_b(predictor_name = None, error_type=None):
    ''' Requires saved data from running "fig4i" '''
    
    global RESULTS_PATH
    
    if predictor_name == None:
        predictor_name = 'XGBoost'

    if error_type is None:
        error_type = 'diff_std' 
        
    fig_num = plt.gcf().number+1
    condition_labels = ['i-D', 'D-D']
    cell_type_labels = ['Place', 'Non-Place']
    
    pval_thresholds = [0.05, 0.005, 0.0005]
    
    num_cell_array_path = RESULTS_PATH + "cell_num_array"
    num_cell_array = scipy.io.loadmat(num_cell_array_path)['cell_num_array']
    
    
    mat_path = RESULTS_PATH + "Position_prediction_by_session_and_cell_type_all_sessions_" + predictor_name + '_' + error_type
    error_array = scipy.io.loadmat(mat_path)['error_array']
    
    
    
    mouse_list_by_condition = [np.arange(4), np.arange(4, 8)]
    
    ## Figure 1: how much is the error related to number of cells?
    
    for plot_idx in range(5): #All cells, only iD, only DD, only P, only NP
        figsize = (6,6); fs = 20
        fig = plt.figure(fig_num, figsize=figsize); fig_num += 1
        ax = plt.gca()
        cell_type_color_list = ['royalblue', 'indianred']
        condition_type_marker_list = ['o', 'x']
        color_counter = -1
        
        error_list_flattened = [[],[]] # 0 is PC, 1 is NPC
        num_cell_flattened = [[],[]]
        
        for cond_idx, mlist in enumerate(mouse_list_by_condition):
            
            if plot_idx == 1 and cond_idx == 1:
                continue
            
            if plot_idx == 2 and cond_idx == 0:
                continue
            
            
            for ctype_idx in range(2):
                
                if plot_idx == 3 and ctype_idx == 1:
                    continue
                if plot_idx == 4 and ctype_idx == 0:
                    continue
                
                if plot_idx == 0:
                    marker = condition_type_marker_list[cond_idx]
                else:
                    marker = 'o'
                
                color_counter += 1
                label = condition_labels[cond_idx] + ' ' + cell_type_labels[ctype_idx]
                err_list = error_array[mlist, :, ctype_idx]
                not_nan = ~np.isnan(err_list)
                err_list = err_list[not_nan]
                num_list = num_cell_array[mlist, :, ctype_idx]
                num_list = num_list[not_nan]
                
                num_cell_flattened[ctype_idx].extend(num_list)
                error_list_flattened[ctype_idx].extend(err_list)
    
                ax.scatter(num_list, err_list, marker = marker, zorder=((ctype_idx+1)%2), s=25, label=label, 
                                     color=cell_type_color_list[ctype_idx], alpha=1)
                            
        if plot_idx == 0:
            ylim_max = ax.get_ylim()[1] #Save it for other plots
            ax.set_ylim((0, ylim_max))
            unique_label = ''
            
        if plot_idx == 1:
            ax.set_ylim((0, ylim_max))
            unique_label = 'id_'
            
        elif plot_idx == 2:
            ax.set_ylim((0, ylim_max))
            unique_label = 'dd_'
            
        ax.legend(fontsize=fs-8, loc='lower right')
        ax.set_title('Error vs cell number', fontsize=fs)
        ax.set_xlabel('Num Cells', fontsize=fs)    
        ax.set_ylabel('%s error (cm)'%error_type, fontsize=fs)
        ax.tick_params(axis='x', labelsize=fs-3)
        ax.tick_params(axis='y', labelsize=fs-3)

        

            
        fig.tight_layout()                
        label_save = "Error_vs_cell_num_" + unique_label + predictor_name + '_' + error_type
        plt.savefig(RESULTS_PATH + label_save + ".svg")
        
        
        continue
            


def figSI11_c(predictor_name=None, error_type=None):
    ''' Compare errors using the same value of place and non-place cells
        - We select the number of cells according to the minimum of place/non-place cells
        - Whichever group has more cells, randomize a selection of them (same amount of shuffles as the "random" comparison)
        
        Wiener takes 30s with 25 shuffles
        XGBoost takes ~17 minutes with 2 shuffles, expected 6 minutes per extra shuffle. 70 minutes for 10, 120 for 20
    '''
    
    np.random.seed(10)
    
    mouse_list = range(8)
    
    session_list = range(9)


    mouse_num = len(mouse_list)
    session_num = len(session_list)
    
    #PCA params
    num_components = '0.9'
    # num_components = 5
    

    
    #Predictor parameters
    cv_folds = 5
    if predictor_name is None:
        predictor_name = 'Wiener'
        
    if error_type is None:
        error_type = 'diff_std'
        
    #Dropout parameters
    shuffling_samples = 25       
    min_cell_num = 3
    
    PCA_dict = np.load(figparams.OUTPUT_PATH + "PCA_dict.npy", allow_pickle=True)[()]
    input_data_dict = np.load(figparams.OUTPUT_PATH + "input_data_dict.npy", allow_pickle=True)[()]

    mouse_num = len(mouse_list)
    session_num = len(session_list)
    
    do_pca_and_predict_partial = partial(figfuns.do_pca_and_predict, pos_max=MAX_POS, num_components=num_components,
                                          cv_folds=cv_folds, predictor_name=predictor_name, error_type=error_type)
    
    error_array_all = np.zeros((mouse_num, session_num, 2)) #Last dimensions are: place cell, non-place cell
    
 
    
    for midx, mnum in enumerate(mouse_list):
        for sidx, snum in enumerate(session_list):
            #Get neuronal data
            position, pca_input_data = input_data_dict[mnum, snum]
            
            #Get PCA data (of full dataset)
            position, pca_data = PCA_dict[mnum, snum]

            #Get place cell data
            place_cell_bool = figfuns.load_place_cell_boolean(mnum, snum, criteria='dombeck').astype(bool)
            place_cell_idxs = np.where(place_cell_bool)[0]
            nplace_cell_idxs = np.where(np.invert(place_cell_bool))[0]
            
            num_place = len(place_cell_idxs)
            num_nplace = len(nplace_cell_idxs)
            num_cells_min = np.minimum(num_place, num_nplace)

            
            if num_cells_min < min_cell_num:
                error_array_all[midx, sidx, 0] = error_array_all[midx, sidx, 1] = np.nan                
                continue
            
            cell_type_with_least = [num_place, num_nplace].index(num_cells_min)
            cell_type_with_most = (1 + cell_type_with_least)%2
            cell_idxs_least = [place_cell_idxs, nplace_cell_idxs][cell_type_with_least]
            cell_idxs_most = [place_cell_idxs, nplace_cell_idxs][cell_type_with_most]

            #Do position prediction with type with least cells
            pca_input_cell_data = pca_input_data[cell_idxs_least]
            pca_data_cell, position_pred, error = do_pca_and_predict_partial(position, pca_input_cell_data)
            error_array_all[midx, sidx, cell_type_with_least] = error
            
            #Then the one with most cells, try different iterations
            error_avg = 0
            for it in range(shuffling_samples):
                selected_neurons = np.random.choice(cell_idxs_most, size=num_cells_min, replace=False)
                pca_input_cell_data = pca_input_data[selected_neurons]
                pca_data_cell, position_pred, error = do_pca_and_predict_partial(position, pca_input_cell_data)
                error_avg += error
                
            error_avg = error_avg/shuffling_samples
            error_array_all[midx, sidx, cell_type_with_most] = error_avg        
            
    label_save = "Place_cell_dropout_same_number_" + predictor_name + '_' + error_type
    scipy.io.savemat(RESULTS_PATH + label_save + ".mat", {"error_array":error_array_all})
    
    figSI11_c_plots(predictor_name, error_type)
    
    return         

def figSI11_c_plots(predictor_name, error_type):
    ''' 
    
    Plots data from "figSI11_c", that function must be run first.
    Allows to play around with the plots without having to re-run the analysis
    
    '''
    
    label_save = "Place_cell_dropout_same_number_" + predictor_name + '_' + error_type
    error_array_path = RESULTS_PATH + label_save + ".mat"
    error_array = scipy.io.loadmat(error_array_path)['error_array']

    session_list = np.arange(9)
    session_num = len(session_list)
    
    mouse_list = range(8)
    mouse_list_by_condition = [np.arange(4), np.arange(4, 8)]
    
    #Ploting parameters
    fig_num = plt.gcf().number+1
    fs = 15
    sample_title = figfuns.generate_title_sessions_and_mice(session_list, mouse_list)
    bar_width = 0.7
    condition_labels = ['i-D', 'D-D']    
    condition_color = ['forestgreen', 'khaki']
    cell_labels = ['Place', 'Non-Place']
    cell_colors = ['royalblue', 'indianred']
       
    ##Figures
    
    #Overall statistics and figures
    fs = 20
    fig, ax = plt.subplots(figsize=(5,8)); fig_num += 1
    # ax = plt.gca()
    
    # for axon_type_idx, mouse_axon_list in enumerate([mouse_id, mouse_dd]):
    for axon_type_idx, mouse_axon_list in enumerate(mouse_list_by_condition):

        prev_errors = None
        xx_pos = axon_type_idx
        
        for cell_type_idx in range(2):
            errors = error_array[mouse_axon_list,:,cell_type_idx].ravel()
            errors = errors[~np.isnan(errors)]
            tot_samples = len(errors)
                        
            if axon_type_idx == 0:
                label = cell_labels[cell_type_idx]
            else:
                label = None
            
            ax.scatter([xx_pos]*tot_samples, errors, color = cell_colors[cell_type_idx], zorder = 1, label=label)
            
            # ax.errorbar([xx_pos], np.average(errors), yerr=np.std(errors), fmt='_', color=cell_colors[cell_type_idx], markersize=20, markeredgewidth=3, elinewidth = 3, zorder=2)
    
            parts = ax.violinplot(errors, [xx_pos], widths=0.5, showmeans=True, showextrema=False, showmedians=False, bw_method='scott')
            for pc in parts['bodies']:
                pc.set_facecolor(cell_colors[cell_type_idx])
                pc.set_alpha(.4)
    
            parts['cmeans'].set_edgecolor(cell_colors[cell_type_idx])
                            
            if cell_type_idx == 0:
                prev_errors = errors
            
            
        #t-test: PC vs NPC
        tstat, pval = scipy.stats.ttest_ind(prev_errors, errors, equal_var=True, permutations=None, alternative='two-sided')
                
        pval_label = figfuns.get_significance_label(pval, thresholds = [0.05, 0.005, 0.0005], asterisk=True)
        dx = 0.25
        x0 = xx_pos + dx; x1 = x0+dx/2
        y0 = np.average(errors); y1 = np.average(prev_errors)
        ax.plot([x0, x1, x1, x0],[y0, y0, y1, y1], 'k', lw=2)
        ax.text(x1+dx/3, (y0+y1)/2, pval_label, fontsize = fs, style='italic')    
    
        
    ax.legend(fontsize=fs-1)
    
    ax.set_xticks([0,1])
    ax.set_xticklabels(condition_labels, fontsize=fs+4)
    ax.set_xlim([-0.5, 1.5])
    ax.set_ylim([0, ax.get_ylim()[1]])
    ax.tick_params(axis='y', labelsize=fs+3) 
    ax.set_ylabel('Prediction error %s (cm)'%error_type, fontsize=fs)
    ax.set_title('Same number dropout, ' + sample_title, fontsize=fs-8)    
    fig.tight_layout()
    
    plt.savefig(RESULTS_PATH + label_save + ".svg")    
    
    

    #Figure 2: just repeating 4g, plotting average error per session and condition, but for all sessions.
    mouse_list_by_condition = [np.arange(4), np.arange(4, 8)]
    # for axon_type_idx, mouse_axon_list in enumerate([mouse_id, mouse_dd]):
        
    figsize = (9,7); fs = 20
    fig = plt.figure(fig_num, figsize=figsize); fig_num += 1
    ax = plt.gca()
            
    for sidx, snum in enumerate(session_list):
        
        for cond_idx, mlist in enumerate(mouse_list_by_condition):
            for ctype_idx in range(2):
                
                if sidx == 0:
                    label = condition_labels[cond_idx] + ' ' + cell_labels[ctype_idx]
                else:
                    label = None
                    
                if ctype_idx == 1:
                    alpha = 0.6
                else:
                    alpha = 1
                    
                sidx_original = list(session_list).index(snum)
                err_list = error_array[mlist, sidx_original, ctype_idx]
                not_nan = ~np.isnan(err_list)
                err_list = err_list[not_nan]

                avg = np.average(err_list)
                std = np.std(err_list)
                
                pos = sidx + (len(session_list) + 1)*cond_idx + (bar_width/4) * (2*ctype_idx - 1)
                
                # plt.scatter([pos]*len(err_list), err_list, zorder=2, s=12, color='dimgray')
                # ax.bar(pos, avg, yerr=std, width=bar_width/2, zorder=1, color=condition_color[cond_idx], label=label, alpha=alpha)
                
                ax.bar(pos, avg, width=bar_width/2, zorder=1, color=condition_color[cond_idx], label=label, alpha=alpha)
                _, caplines, _ = ax.errorbar(pos, avg, yerr=[[0],[std]], lolims=True, capsize = 0, ls='None', color='k')
                caplines[0].set_marker('_')
                caplines[0].set_markersize(8 * (9/len(session_list)))                
        
    slabels = [figparams.SESSION_NAMES[snum] for snum in session_list]
    xlabels = slabels + [''] + slabels
    xpos = np.arange(2*session_num+1)

    ax.set_title('Prediction error, same number dropout', fontsize=fs)
    ax.set_ylabel('%s error (cm)'%error_type, fontsize=fs)
    ax.set_xlabel('Session', fontsize=fs)
    ax.set_xticks(xpos, xlabels, fontsize=fs - 6)
    ax.legend(fontsize=fs-10, loc='upper center')
    ax.tick_params(axis='y', labelsize=fs-3)

    label_save = "Place_cell_same_number_dropout_by_session_and_cell_type_" + predictor_name + '_' + error_type
    plt.savefig(RESULTS_PATH + label_save + ".svg")


            
            
            
            
def figSI11_d_e(predictor_name=None, error_type=None):
    ''' Compare errors using the same value of a cell type across all mice for each session
        e.g.: if you have this amount of place cells for session B3: 4, 10, 20, 7 (id) and 23, 32, 25, 40 (dd)
                You select 4 for all of them, randomizing the chosen ones the mice with larger amounts
        
        Wiener takes 30s with 25 shuffles
        XGBoost takes 4h minutes with 25 shuffles
    '''
    
    np.random.seed(10)
    
    mouse_list = range(8)
    # mouse_list = range(4)
    # mouse_list = range(4,8)
    # mouse_list = [2, 4]
    
    session_list = range(9)
    # session_list = [0,1,2]
    # session_list = range(7,9)
    # session_list = [2,3,6,7]

    mouse_num = len(mouse_list)
    session_num = len(session_list)
    
    #PCA params
    num_components = '0.9'
    # num_components = 5
    

    
    #Predictor parameters
    cv_folds = 5
    if predictor_name is None:
        predictor_name = 'Wiener'
        
    if error_type is None:
        error_type = 'diff_std'
        
    #Dropout parameters
    shuffling_samples = 25
    dropout_cell_threshold = 4
    
    PCA_dict = np.load(figparams.OUTPUT_PATH + "PCA_dict.npy", allow_pickle=True)[()]
    input_data_dict = np.load(figparams.OUTPUT_PATH + "input_data_dict.npy", allow_pickle=True)[()]
    
    num_cell_array_path = RESULTS_PATH + "cell_num_array"
    num_cell_array = scipy.io.loadmat(num_cell_array_path)['cell_num_array']

    mouse_num = len(mouse_list)
    session_num = len(session_list)
    
    do_pca_and_predict_partial = partial(figfuns.do_pca_and_predict, pos_max=MAX_POS, num_components=num_components,
                                          cv_folds=cv_folds, predictor_name=predictor_name, error_type=error_type)
        
    error_array_all = np.zeros((mouse_num, session_num, 2)) #Last dimensions are: place cell, non-place cell
    
    num_cell_array_filtered = ma.masked_less(num_cell_array.astype(int), dropout_cell_threshold)
    minimum_cell_num_per_session = np.min(num_cell_array_filtered, axis=0) #Shape "num sessions X 2"
    minimum_cell_num_per_session[:,:] = np.min(minimum_cell_num_per_session) #Force the same number for all mice and sessions
        
    for sidx, snum in enumerate(session_list):
                
        for midx, mnum in enumerate(mouse_list):
            
            print('Session number: %d, Mouse number: %d' %(mnum, snum))

        
            #Get neuronal data
            position, pca_input_data = input_data_dict[mnum, snum]
            
            #Get PCA data (of full dataset)
            position, pca_data = PCA_dict[mnum, snum]

            #Get place cell data
            place_cell_bool = figfuns.load_place_cell_boolean(mnum, snum, criteria='dombeck').astype(bool)
            place_cell_idxs = np.where(place_cell_bool)[0]
            nplace_cell_idxs = np.where(np.invert(place_cell_bool))[0]
            
            for cell_type_idx, cell_idxs in enumerate([place_cell_idxs, nplace_cell_idxs]):
                num_cells = len(cell_idxs)
                num_cells_min = minimum_cell_num_per_session[sidx, cell_type_idx]
                
                if num_cells < dropout_cell_threshold:
                    error_array_all[midx, sidx, cell_type_idx] = np.nan
                    continue
                
                #Try different iterations choosing that amount of cells for that session
                
                ## Efficency: if there are less possible combinations than the selected "shuffling", choose those instead
                max_combinations = comb(len(cell_idxs), num_cells_min) #Possible combinations
                if max_combinations < shuffling_samples:
                    selected_neurons_list = list(combinations(cell_idxs, num_cells_min))
                    selected_neurons_list = [list(e) for e in selected_neurons_list]
                else:
                    selected_neurons_list = [np.random.choice(cell_idxs, size=num_cells_min, replace=False) for i in range(shuffling_samples)]

                # selected_neurons_list = [np.random.choice(cell_idxs, size=num_cells_min, replace=False) for i in range(shuffling_samples)]
        
                num_iterations = len(selected_neurons_list)
                error_avg = 0
                for it, selected_neurons in enumerate(selected_neurons_list):
                    # selected_neurons = np.random.choice(cell_idxs, size=num_cells_min, replace=False)
                    pca_input_cell_data = pca_input_data[selected_neurons]
                    pca_data_cell, position_pred, error = do_pca_and_predict_partial(position, pca_input_cell_data)
                    error_avg += error

                
                    
                error_avg = error_avg/num_iterations
                error_array_all[midx, sidx, cell_type_idx] = error_avg                 
      
            
      
       
    label_save = "Place_cell_dropout_same_number_across_sessions_" + predictor_name + '_' + error_type
    scipy.io.savemat(RESULTS_PATH + label_save + ".mat", {"error_array":error_array_all})
    
    figSI11_d_e_plots(predictor_name, error_type)
                
    return         

def figSI11_d_e_plots(predictor_name, error_type):    
    
    ''' Plots data from "figSI11_d_e", that function must be run first.
        Allows to play around with the plots without having to re-run the analysis
    '''
    label_save = "Place_cell_dropout_same_number_across_sessions_" + predictor_name + '_' + error_type
    error_array_path = RESULTS_PATH + label_save + ".mat"
    error_array = scipy.io.loadmat(error_array_path)['error_array']

    session_list = np.arange(9)
    
    mouse_list = range(8)
    mouse_list_by_condition = [np.arange(4), np.arange(4, 8)]   
    
    #Ploting parameters
    sessions_to_plot = [2,3,6,7]
    session_num_to_plot = len(sessions_to_plot)
    fig_num = plt.gcf().number+1
    fs = 15
    sample_title = figfuns.generate_title_sessions_and_mice(session_list, mouse_list)
    bar_width = 0.7
    condition_labels = ['i-D', 'D-D']    
    condition_color = ['forestgreen', 'khaki']
    cell_labels = ['Place', 'Non-Place']
    cell_colors = ['royalblue', 'indianred']    
    
    #Overall statistics and figures
    fs = 20
    fig, ax = plt.subplots(figsize=(5,8)); fig_num += 1
    # ax = plt.gca()
    
    for axon_type_idx, mouse_axon_list in enumerate(mouse_list_by_condition):

        prev_errors = None
        xx_pos = axon_type_idx
        
        for cell_type_idx in range(2):
            errors = error_array[mouse_axon_list,:,cell_type_idx].ravel()
            errors = errors[~np.isnan(errors)]
            tot_samples = len(errors)
                        
            if axon_type_idx == 0:
                label = cell_labels[cell_type_idx]
            else:
                label = None
            
            ax.scatter([xx_pos]*tot_samples, errors, color = cell_colors[cell_type_idx], zorder = 1, label=label)
            
            # ax.errorbar([xx_pos], np.average(errors), yerr=np.std(errors), fmt='_', color=cell_colors[cell_type_idx], markersize=20, markeredgewidth=3, elinewidth = 3, zorder=2)
    
            parts = ax.violinplot(errors, [xx_pos], widths=0.5, showmeans=True, showextrema=False, showmedians=False, bw_method='scott')
            for pc in parts['bodies']:
                pc.set_facecolor(cell_colors[cell_type_idx])
                pc.set_alpha(.4)
    
            parts['cmeans'].set_edgecolor(cell_colors[cell_type_idx])
                            
            if cell_type_idx == 0:
                prev_errors = errors
            
            
        #t-test: PC vs NPC
        tstat, pval = scipy.stats.ttest_ind(prev_errors, errors, equal_var=True, permutations=None, alternative='two-sided')
                
        pval_label = figfuns.get_significance_label(pval, thresholds = [0.05, 0.005, 0.0005], asterisk=True)
        dx = 0.25
        x0 = xx_pos + dx; x1 = x0+dx/2
        y0 = np.average(errors); y1 = np.average(prev_errors)
        ax.plot([x0, x1, x1, x0],[y0, y0, y1, y1], 'k', lw=2)
        ax.text(x1+dx/3, (y0 + y1)/2, pval_label, fontsize = fs, style='italic')    
    
        
    ax.legend(fontsize=fs-1)
    
    ax.set_xticks([0,1])
    ax.set_xticklabels(condition_labels, fontsize=fs+4)
    ax.set_xlim([-0.5, 1.5])
    ax.set_ylim([0, ax.get_ylim()[1]])
    ax.tick_params(axis='y', labelsize=fs+3) 
    ax.set_ylabel('Prediction error %s (cm)'%error_type, fontsize=fs)
    ax.set_title('Same number across sessions dropout, ' + sample_title, fontsize=fs-8)    
    fig.tight_layout()
    
    plt.savefig(RESULTS_PATH + label_save + ".svg")    
    
    
    
    #Figure 2: just repeating 4g, plotting average error per session and condition, but for all sessions.        
    figsize = (9,7); fs = 20
    fig = plt.figure(fig_num, figsize=figsize); fig_num += 1
    ax = plt.gca()
            
    for sidx, snum in enumerate(sessions_to_plot):
        
        for cond_idx, mlist in enumerate(mouse_list_by_condition):
            for ctype_idx in range(2):
                
                if sidx == 0:
                    label = condition_labels[cond_idx] + ' ' + cell_labels[ctype_idx]
                else:
                    label = None
                    
                if ctype_idx == 1:
                    alpha = 0.6
                else:
                    alpha = 1
                    
                sidx_original = list(session_list).index(snum)
                err_list = error_array[mlist, sidx_original, ctype_idx]
                not_nan = ~np.isnan(err_list)
                err_list = err_list[not_nan]
    
                avg = np.average(err_list)
                std = np.std(err_list)
                
                pos = sidx + (len(sessions_to_plot) + 1)*cond_idx + (bar_width/4) * (2*ctype_idx - 1)
                
                #Basic error bars and show points
                plt.scatter([pos]*len(err_list), err_list, zorder=2, s=12, color='dimgray')
                ax.bar(pos, avg, yerr=std, width=bar_width/2, zorder=1, color=condition_color[cond_idx], label=label, alpha=alpha)
                
                #Half error bars, no individual points
                # ax.bar(pos, avg, width=bar_width/2, zorder=1, color=condition_color[cond_idx], label=label, alpha=alpha)
                # _, caplines, _ = ax.errorbar(pos, avg, yerr=[[0],[std]], lolims=True, capsize = 0, ls='None', color='k')
                # caplines[0].set_marker('_')
                # caplines[0].set_markersize(8 * (9/len(sessions_to_plot)))                 
        
    slabels = [figparams.SESSION_NAMES[snum] for snum in sessions_to_plot]
    xlabels = slabels + [''] + slabels
    xpos = np.arange(2*session_num_to_plot+1)
    
    ax.set_title('Prediction error, same number across sessions dropout', fontsize=fs)
    ax.set_ylabel('%s error (cm)'%error_type, fontsize=fs)
    ax.set_xlabel('Session', fontsize=fs)
    ax.set_xticks(xpos, xlabels, fontsize=fs - 6)
    ax.legend(fontsize=fs-6, loc='upper center')
    ax.tick_params(axis='y', labelsize=fs-3)
    
    label_save = "Place_cell_dropout_same_number_across_sessions_by_session_and_cell_type_" + predictor_name + '_' + error_type
    plt.savefig(RESULTS_PATH + label_save + ".svg")
    
    

            
            
def figSI11_d_and_e_to_excel(predictor='XGBoost', error_type='Wiener'):
    ''' Exports data used in figSI11_c to excel file for further processing '''

    slist = [2,3,6,7]

    mouse_list_by_condition = [np.arange(4), np.arange(4, 8)]

     
    mat_filename = "Place_cell_dropout_same_number_across_sessions_" + predictor + "_" + error_type
    mat_path = RESULTS_PATH + mat_filename + ".mat"
    error_array = scipy.io.loadmat(mat_path)['error_array']
        
    axon_type_label = ['iD', 'DD']
    session_label = ['Day 1', 'Day 2', 'Day 3', 'Day 4']
    cell_type_label = ['Place', 'Non-Place']
    
    # print(error_array.shape)
    # return
    
    #Create csv to analyze contradictions later
    with open(RESULTS_PATH + mat_filename + ".csv", mode='w', newline='') as results_csv:
        csvw = csv.writer(results_csv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csvw.writerow(['Cell Type', 'Axon Type', 'Day', 'Error (' + error_type + ')'])
        csvw.writerow([])
        
        
        for cell_type_idx in range(2):
        
            for cond_idx, mlist in enumerate(mouse_list_by_condition):
                error_array_by_axon_type = error_array[mlist]
            
                for sidx, snum in enumerate(slist):
                    errors = error_array_by_axon_type[:, snum, cell_type_idx] #### Assumes that snum is the same as sidx for the results array!!
                    csvw.writerow([cell_type_label[cell_type_idx], axon_type_label[cond_idx], session_label[sidx]] + list(errors))
        
            csvw.writerow([])
        csvw.writerow([])
        

        

        



with h5py.File(figparams.FAT_CLUSTER_PATH, 'r') as fat_cluster:

    if __name__ == '__main__':
        
        # project_parameters
        
        tt = time.time()
        main()
        np.random.seed(1)
        print('Time Ellapsed: %.1f' % (time.time() - tt))
        plt.show()