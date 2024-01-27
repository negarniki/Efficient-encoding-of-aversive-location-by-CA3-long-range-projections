# -*- coding: utf-8 -*-
"""

Functions to create the low-dimensional analysis subplots in Figure 4 for the paper "Efficient encoding of aversive location by CA3 long-range projections"

Created by Albert Miguel López

                                      
"""




#Global project variables

import os
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

from scipy.stats import linregress

from functools import partial


            
def main():
    
    
    ''' 
    
        - Running times -
            
        save_pca_data: around 10s
        
        compute_and_plot_fig4_subfigures:
            Using Wiener decoder: around 5s
            Using XGBoost decoder: around 8 minutes
    
    
    '''
    
    # # # #### RUN THIS FIRST TO PROCESS RAW DATA INTO PCA (can uncomment later) ####
    # figfuns.save_pca_data(np.arange(8), np.arange(9))
    # # # #### RUN THIS FIRST TO PROCESS RAW DATA INTO PCA (can uncomment later) ####

    # ### Plot all figure 4 subfigures ###
    compute_and_plot_fig4_subfigures()


    ### Uncomment to compute and plot figures independently ###
    # fig4e()
    # fig4f()
    # fig4g()
    # fig4h()
    # fig4i()
    # fig4h_to_excel()
    # fig4h_by_cell_type_to_excel()
    
    ### Uncomment to just plot the figures (assumes you have done the computations before) ###
    # plot_precomputed_fig4_subfigures()
    
    
def compute_and_plot_fig4_subfigures():
    ''' Using all the predictors takes around 30 minutes for place cell figure, 8 minutes for the fig4h
    
        Wiener takes around 40s for everything
    
        For XGBoost:
            - 4f 9mins
            - 4g place cell 12mins
            - 4h place cell randomization 70min with 10 shuffles (~6-7 per shuffle)
            - 4i 50min with 10 shuffles (~5 mins per shuffle)
            
            - For 20 shuffles in each, 4h total.
            
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
        ### Main figures ###
        fig4e()
        fig4f()
        fig4g(predictor, error_type=error_type)
        fig4h(predictor, error_type=error_type)
        fig4i(predictor, error_type=error_type)
        fig4h_to_excel(predictor, error_type=error_type)
        
        
def plot_precomputed_fig4_subfigures():
    ''' 
    Assuming that the "fig4" functions have been run, plot the data.
    Convenience function to play with the plots without having to recalculate everything
    
    '''
    # error_type = ['sse']
    error_type = 'diff_std'

    # predictor_list = ['Wiener', 'Wiener Cascade', 'Kalman', 'SVR', 'XGBoost']
    predictor_list = ['Wiener']    
    # predictor_list = ['XGBoost'] 
    
    for predictor in predictor_list:
        ### Main figures ###
        fig4e()
        fig4f()
        fig4g_plots(predictor, error_type)        
        fig4h_plots(predictor, error_type=error_type)
        fig4i_plots(predictor, error_type=error_type)

def fig4e():
    
    
    mnum = 2
    snum = 0
    
    #PCA params
    num_components = 3
    num_components = '.9'
    
    #Ploting parameters
    fig_num = plt.gcf().number+1
    fs = 15
    figsize = (8,8)

    PCA_dict = np.load(figparams.OUTPUT_PATH +"PCA_dict.npy", allow_pickle=True)[()]
    position, pca_data = PCA_dict[mnum, snum]
    fig, fig_num = figfuns.plot_pca_and_average(pca_data[:3, :], position, fig_num, figsize, cbar=True, 
                                           pca_distance_bin_size=25, fs=fs, axis='on', ms=6, lw=1)

def fig4f():
    
    
    mnum = 6
    snum = 3
    
    #PCA params
    num_components = 3
    num_components = '.9'
    
    #Ploting parameters
    fig_num = plt.gcf().number+1
    fs = 15
    figsize = (8,8)

    PCA_dict = np.load(figparams.OUTPUT_PATH +"PCA_dict.npy", allow_pickle=True)[()]
    position, pca_data = PCA_dict[mnum, snum]
    fig, fig_num = figfuns.plot_pca_and_average(pca_data[:3, :], position, fig_num, figsize, cbar=True, 
                                           pca_distance_bin_size=25, fs=fs, axis='on', ms=6, lw=1)
        

def fig4g(predictor_name = None, error_type=None):
    ''' Plots std error of position prediction against shuffled data
        Data is shuffled by cycling the time dimension
        
        Computation times:
        Wiener: 1.5s
        Wiener Cascade: 2s 
        Kalman: 35s
        XGBoost: 15min with 2 shuffles, 50min with 10 (so 5 per shuffle + initial 5)
        SVR: 2min
    '''
            
    #Set seed for reproducibility
    np.random.seed(5)
    
    #Take only baseline sessions
    session_list = [0,1,2]    
    
    #Take all mice
    mouse_list = np.arange(8)
    
    

    
    #Preprocessing and pca parameters
    num_components = 3
    num_components = '.9'
    
    #Predictor parameters
    if predictor_name is None:
        predictor_name = 'Wiener' #'Wiener', 'Wiener Cascade', 'Kalman', 'XGBoost', 'SVR'
    if error_type is None:
        error_type = 'diff_std' 
    cv_folds = 5
    number_of_random_shuffles = 1

    
    session_num = len(session_list)
    mouse_num = len(mouse_list)

        
    input_data_dict = np.load(figparams.OUTPUT_PATH + "input_data_dict.npy", allow_pickle=True)[()]
    
    
    do_pca_and_predict_partial = partial(figfuns.do_pca_and_predict, pos_max=MAX_POS, num_components=num_components,
                                         cv_folds=cv_folds, predictor_name=predictor_name, error_type=error_type)
    
    error_array = np.zeros((mouse_num, session_num, number_of_random_shuffles + 1))
    
    for midx, mnum in enumerate(mouse_list):
        for sidx, snum in enumerate(session_list):
            print('Mouse %d, session number %d'%(midx, sidx))
            #Get neuronal data
            position, pca_input_data = input_data_dict[mnum, snum]
            timepoints = len(position)
            
            #Unshuffled data
            pca_data, position_pred, error = do_pca_and_predict_partial(position, pca_input_data)
            error_array[midx, sidx, -1] = error
            
            #Shuffled data
            shift_amounts = np.random.randint(-timepoints, timepoints, size=number_of_random_shuffles)
            for shift_idx, shift in enumerate(shift_amounts):
                idxs = list(range(len(position)))
                idxs = idxs[shift:] + idxs[:shift]
                # position_shuffled = position[idxs]
                pca_input_data_shifted = pca_input_data[:, idxs]
                
                pca_data, position_pred, error = do_pca_and_predict_partial(position, pca_input_data_shifted)
                error_array[midx, sidx, shift_idx] = error
        

    label_save = "Position_prediction_shuffled" + "_" + predictor_name + '_' + error_type
    scipy.io.savemat(RESULTS_PATH + label_save + ".mat", {"error_array":error_array})
    
    fig4g_plots(predictor_name, error_type)

    
    return

def fig4g_plots(predictor_name, error_type):
    
    ''' 
    
    Plots data from "fig4g", that function must be run first.
    Allows to play around with the plots without having to re-run the analysis
    
    '''
    
    
    label_save = "Position_prediction_shuffled" + "_" + predictor_name + '_' + error_type
    error_array_path = RESULTS_PATH + label_save + ".mat"
    error_array = scipy.io.loadmat(error_array_path)['error_array']
    mouse_list_by_condition = [np.arange(4), np.arange(4, 8)]

    #Ploting parameters
    fig_num = plt.gcf().number+1
    condition_color = ['forestgreen', 'khaki']
    pval_thresholds = [0.05, 0.005, 0.0005]
    
    #Figure 1: just repeating 4e, plotting average error per condition vs shuffled
    figsize = (6,12); fs = 25
    fig = plt.figure(fig_num, figsize=figsize); fig_num += 1
    ax = plt.gca()
    
    err_list_both = []
    
    for cond_idx, mlist in enumerate(mouse_list_by_condition):
        err_list_original = error_array[mlist, :, -1].ravel()
        err_list_shuffled = error_array[mlist, :, :-1].ravel()
        
        for err_list_idx, err_list in enumerate([err_list_original, err_list_shuffled]):
            
            if err_list_idx == 1:
                color = 'gray'
                alpha = 0.8
                
            else:
                alpha = 1
                color = condition_color[cond_idx]
            
            avg = np.average(err_list)
            std = np.std(err_list)
            
            xpos = err_list_idx + 3*cond_idx
            # plt.scatter([xpos]*len(err_list), err_list, zorder=2, s=12, color='dimgray')
            # ax.bar(xpos, avg, yerr=std, width=.7, zorder=1, color=color, alpha = alpha)

            ax.bar(xpos, avg, width=.7, zorder=1, color=color, alpha = alpha)            
            _, caplines, _ = ax.errorbar(xpos, avg, yerr=[[0],[std]], lolims=True, capsize = 0, ls='None', color='k')
            caplines[0].set_marker('_')
            caplines[0].set_markersize(25)

        equal_var = True
        if error_type == 'diff_std':
            equal_var = False
            
            
        tstat, pval = scipy.stats.ttest_ind(err_list_original, err_list_shuffled, equal_var=equal_var , permutations=None, alternative='two-sided')
        pval_label = figfuns.get_significance_label(pval, thresholds = pval_thresholds, asterisk=True, ns=True)
            
        dy = .5
        x0 = xpos - 1; x1 = xpos
        y0 = np.max(np.hstack((err_list_original, err_list_shuffled))) + 1; y1 = y0 + dy
        ax.plot([x0, x0, x1, x1],[y0, y1, y1, y0], 'k', lw=2)
        ax.text((x1 + x0)/2 - 0.15, y1+dy, pval_label, fontsize = fs, style='italic')            

        err_list_both.append(np.copy(err_list_original))
        
    tstat, pval = scipy.stats.ttest_ind(err_list_both[0], err_list_both[1], equal_var=True, permutations=None, alternative='two-sided')
    pval_label = figfuns.get_significance_label(pval, thresholds = pval_thresholds, asterisk=True, ns=True)
        
    dy = .5
    x0 = 0; x1 = 3
    y0 = (y1+dy)*1.1; y1 = y0 + dy
    ax.plot([x0, x0, x1, x1],[y0, y1, y1, y0], 'k', lw=2)
    ax.text((x1 + x0)/2 - 0.15, y1+dy, pval_label, fontsize = fs, style='italic')            
    
    
    xlabels = ['i-D', 'Shuffle', 'D-D', 'Shuffle']
    xpos = [0,1,3,4]
    ax.set_xticks(xpos, xlabels, fontsize=fs-2)
    ax.set_ylim([0, y1 * 1.15])

    ax.set_title('Normal vs shuffled prediction', fontsize=fs)
    ax.set_ylabel('%s error (cm)'%error_type, fontsize=fs)
    ax.tick_params(axis='y', labelsize=fs-3)
    fig.tight_layout()
    
    plt.savefig(RESULTS_PATH + label_save + ".svg")    




def fig4h(predictor_name = None, error_type = None):
    ''' Plots std error of position prediction by session, averaging over mice
        Wiener: 1.5s
        Wiener Cascade: 2s 
        Kalman: 35s
        XGBoost: 9min
        SVR: 2min
    '''
        
    session_list = np.arange(9)
    # session_list = [2,3,6,7]

    session_list_to_plot = [2,3,6,7]
    # session_list = range(7,9)
    # session_list = [2,3]
    
    session_num = len(session_list)
    
    #PCA parameters
    num_components = '.9'
    
    #Predictor parameters
    cv_folds = 5
    if predictor_name is None:
        predictor_name = 'Wiener'
    if error_type is None:
        error_type = 'diff_std'         
    
    input_data_dict = np.load(figparams.OUTPUT_PATH + "input_data_dict.npy", allow_pickle=True)[()]
    
    
    do_pca_and_predict_partial = partial(figfuns.do_pca_and_predict, pos_max=MAX_POS, num_components=num_components,
                                          cv_folds=cv_folds, predictor_name=predictor_name, error_type=error_type)
    
    error_array_all = np.zeros((8, session_num))
        
    
    for sidx, snum in enumerate(session_list):
        for midx, mnum in enumerate(range(8)):        
            #Get neuronal data
            position, pca_input_data = input_data_dict[mnum, snum]    
            pca_data, position_pred, err = do_pca_and_predict_partial(position, pca_input_data)
            error_array_all[midx, sidx] = err

    

    #Figure 1: just repeating 4e, plotting average error per session and condition
        
    
    for slist_idx, slist in enumerate([session_list_to_plot, session_list]):

        if slist_idx == 0:
            label_save = "Position_prediction_key_sessions_" + predictor_name + '_' + error_type
            
        if slist_idx == 1:
            label_save = "Position_prediction_all_sessions_" + predictor_name + '_' + error_type
            
        scipy.io.savemat(RESULTS_PATH + label_save + ".mat", {"error_array":error_array_all})               
    fig4h_plots(predictor_name, error_type)

    return


def fig4h_plots(predictor_name, error_type):
    
    ''' 
    
    Plots data from "fig4h", that function must be run first.
    Allows to play around with the plots without having to re-run the analysis
    
    '''
    
    
    #Ploting parameters
    session_list = np.arange(9)
    session_list_to_plot = [2,3,6,7] 
    mouse_list_by_condition = [np.arange(4), np.arange(4, 8)]      
    
    fig_num = plt.gcf().number+1
    fs = 15
    condition_labels = ['i-D', 'D-D']
    condition_color = ['forestgreen', 'khaki']
    pval_thresholds = [0.05, 0.005, 0.0005]
    
    for slist_idx, slist in enumerate([session_list_to_plot, session_list]):      
        
        if slist_idx == 0:
            label_save = "Position_prediction_key_sessions_" + predictor_name + '_' + error_type
            
        if slist_idx == 1:
            label_save = "Position_prediction_all_sessions_" + predictor_name + '_' + error_type
        
        error_array_path = RESULTS_PATH + label_save + ".mat"
        error_array = scipy.io.loadmat(error_array_path)['error_array']
        
        
        slist_idxs = [list(session_list).index(s) for s in slist]
        figsize = (6,8); fs = 20
        fig = plt.figure(fig_num, figsize=figsize); fig_num += 1
        ax = plt.gca()
        
        for sidx_loopidx, sidx in enumerate(slist_idxs):
            for cond_idx, mlist in enumerate(mouse_list_by_condition):
                
                if sidx_loopidx == 0:
                    label = condition_labels[cond_idx]
                else:
                    label = None
                err_list = error_array[mlist, sidx]
                avg = np.average(err_list)
                std = np.std(err_list)
                
                pos = sidx_loopidx + (len(slist) + 1)*cond_idx
                
                ax.bar(pos, avg, width=.7, zorder=1, color=condition_color[cond_idx], label=label)
                _, caplines, _ = ax.errorbar(pos, avg, yerr=[[0],[std]], lolims=True, capsize = 0, ls='None', color='k')
                caplines[0].set_marker('_')
                caplines[0].set_markersize(10 * (9/len(slist_idxs)))
                
                if sidx_loopidx != len(slist)-1 and slist_idx == 0:
                    tstat, pval = scipy.stats.ttest_ind(error_array[mlist, sidx], error_array[mlist, sidx+1], equal_var=True, permutations=None, alternative='two-sided')
                    pval_label = figfuns.get_significance_label(pval, thresholds = pval_thresholds)
    
    
        slabels = [figparams.SESSION_NAMES[snum] for snum in slist]
        xlabels = slabels + [''] + slabels
        xpos = np.arange(2*len(slist)+1)
            
        ax.set_title('Error of position prediction', fontsize=fs)
        ax.set_ylabel('%s error (cm)'%error_type, fontsize=fs)
        ax.set_xlabel('Session', fontsize=fs)
        xticks_fs = fs + 1.5*(4 - len(slist))
        ax.set_xticks(xpos, xlabels, fontsize=xticks_fs)
        ax.legend(fontsize=fs-3, loc='upper center')
        ax.tick_params(axis='y', labelsize=fs-3)
        fig.tight_layout()
        plt.savefig(RESULTS_PATH + label_save + ".svg")    




def fig4i(predictor_name=None, fig_num=1, error_type = None):
    ''' Plots sse error of position prediction by session, averaging over mice, for place and non-place cells 
    
        Wiener takes around 12 seconds
        Kalman takes around 2:30 minutes
        XGBoost takes around 12 minutes
        SVR takes around 5 minutes
        '''
        
        
    session_list = [2,3,6,7]
    # session_list = [0,1,2]
    # session_list = range(7,9)
    # session_list = [2,3]
    session_list = np.arange(9)
    
    session_num = len(session_list)
    
    #PCA params
    num_components = 3
    num_components = '.9'

    
    #Predictor parameters
    cv_folds = 5    
    if predictor_name is None:
        predictor_name = 'Wiener'
    if error_type is None:
        error_type = 'diff_std' 
        
    #Place cell parameters
    min_cell_num = 3
    
    input_data_dict = np.load(figparams.OUTPUT_PATH + "input_data_dict.npy", allow_pickle=True)[()]
    PCA_dict = np.load(figparams.OUTPUT_PATH + "PCA_dict.npy", allow_pickle=True)[()]
    
    
    do_pca_and_predict_partial = partial(figfuns.do_pca_and_predict, pos_max=MAX_POS, num_components=num_components,
                                         cv_folds=cv_folds, predictor_name=predictor_name, error_type=error_type)
    
    error_array_all = np.zeros((8, session_num, 2)) #Third axis is cell type
    num_cell_array = np.zeros((8, session_num, 2)) #Third axis is cell type
        
    
    for midx, mnum in enumerate(range(8)):
        for sidx, snum in enumerate(session_list):
            #Get neuronal data            
            position, pca_input_data = input_data_dict[mnum, snum]            
                
            #Get PCA data (of full dataset)
            position, pca_data_original = PCA_dict[mnum, snum]
            

            #Get place cell data
            place_cell_bool = figfuns.load_place_cell_boolean(mnum, snum, criteria='dombeck').astype(bool)
            place_cell_idxs = np.where(place_cell_bool)[0]
            nplace_cell_idxs = np.where(np.invert(place_cell_bool))[0]
            
            for cell_type_idx, cell_idxs in enumerate([place_cell_idxs, nplace_cell_idxs]):
                # print(sidx, snum, cell_type_idx)

                num_cells = len(cell_idxs)
                num_cell_array[midx, sidx, cell_type_idx] = num_cells
                
                if num_cells < min_cell_num:
                    error_array_all[midx, sidx, cell_type_idx] = np.nan    
                    continue
                pca_input_cell_data = pca_input_data[cell_idxs]
                # print(position.shape, pca_input_cell_data.shape)
                pca_data_cell, position_pred, error = do_pca_and_predict_partial(position, pca_input_cell_data)
                error_array_all[midx, sidx, cell_type_idx] = error
        
        
    scipy.io.savemat(RESULTS_PATH + "cell_num_array" + ".mat", {"cell_num_array":num_cell_array})

    #Figure 1: just repeating 4g, plotting average error per session and condition        
    label_save = "Position_prediction_by_session_and_cell_type_all_sessions_" + predictor_name + '_' + error_type            
    scipy.io.savemat(RESULTS_PATH + label_save + ".mat", {"error_array":error_array_all})
    
    fig4i_plots(predictor_name, error_type)

            
    return


def fig4i_plots(predictor_name, error_type):
    ''' 
    
    Plots data from "fig4i", that function must be run first.
    Allows to play around with the plots without having to re-run the analysis
    
    '''
    
    label_save = "Position_prediction_by_session_and_cell_type_all_sessions_" + predictor_name + '_' + error_type
    error_array_path = RESULTS_PATH + label_save + ".mat"
    error_array = scipy.io.loadmat(error_array_path)['error_array']
        
    #Ploting parameters
    fig_num = plt.gcf().number+1
    fs = 15
    figsize = (8,8)
    condition_labels = ['i-D', 'D-D']
    cell_type_labels = ['place', 'non-place']
    condition_color = ['forestgreen', 'khaki']
    bar_width = 0.7
    
    session_list = np.arange(9)
    session_lists_to_plot = ([2,3,6,7], session_list)
    mouse_list_by_condition = [np.arange(4), np.arange(4, 8)]

    for session_list_idx, session_list_to_plot in enumerate(session_lists_to_plot):
        
        figsize = (9,7); fs = 20
        fig = plt.figure(fig_num, figsize=figsize); fig_num += 1
        ax = plt.gca()
        
        session_num_to_plot = len(session_list_to_plot)
        
        for sidx, snum in enumerate(session_list_to_plot):
            
            for cond_idx, mlist in enumerate(mouse_list_by_condition):
                for ctype_idx in range(2):
                    
                    if sidx == 0:
                        label = condition_labels[cond_idx] + ' ' + cell_type_labels[ctype_idx]
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
                    
                    pos = sidx + (len(session_list_to_plot) + 1)*cond_idx + (bar_width/4) * (2*ctype_idx - 1)
                    
                    plt.scatter([pos]*len(err_list), err_list, zorder=2, s=12, color='dimgray')
                    ax.bar(pos, avg, yerr=std, width=bar_width/2, zorder=1, color=condition_color[cond_idx], label=label, alpha=alpha)
                    
            
        slabels = [figparams.SESSION_NAMES[snum] for snum in session_list_to_plot]
        xlabels = slabels + [''] + slabels
        xpos = np.arange(2*session_num_to_plot+1)
    
        ax.set_title('Error of position prediction', fontsize=fs)
        ax.set_ylabel('%s error (cm)'%error_type, fontsize=fs)
        ax.set_xlabel('Session', fontsize=fs)
        ax.set_xticks(xpos, xlabels, fontsize=fs - 6*session_list_idx)
        ax.legend(fontsize=fs-6, loc='best')
        ax.tick_params(axis='y', labelsize=fs-3)
        fig.tight_layout()
        if session_list_idx == 0:
            session_list_idx_label = 'key_sessions'
        else:
            session_list_idx_label = 'all_sessions'
        label_save = "Position_prediction_by_session_and_cell_type_" + session_list_idx_label + "_" + predictor_name + '_' + error_type
        plt.savefig(RESULTS_PATH + label_save + ".svg")
                
        
    fs = 20
    fig, ax = plt.subplots(figsize=(5,8)); fig_num += 1
    # ax = plt.gca()
    cell_colors = ['royalblue', 'indianred']
    
    for axon_type_idx, mouse_axon_list in enumerate(mouse_list_by_condition):
    
        prev_errors = None
        xx_pos = axon_type_idx
        
        for cell_type_idx in range(2):
            errors = error_array[mouse_axon_list,:,cell_type_idx].ravel()
            errors = errors[~np.isnan(errors)]
            tot_samples = len(errors)
                        
            if axon_type_idx == 0:
                label = cell_type_labels[cell_type_idx]
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
        ax.text(x1+dx/3, (y0+y1)/2.05, pval_label, fontsize = fs, style='italic')    
    
        
    ax.legend(fontsize=fs-1)
    
    ax.set_xticks([0,1])
    ax.set_xticklabels(condition_labels, fontsize=fs+4)
    ax.set_xlim([-0.5, 1.5])
    
    ax.set_ylim([0, ax.get_ylim()[1]])
    ax.tick_params(axis='y', labelsize=fs+3) 
    ax.set_ylabel('Prediction error %s (cm)'%error_type, fontsize=fs)
    
    ax.set_title('Place vs Non-place cells', fontsize=fs-8)    
    fig.tight_layout()
    
    label_save = "Position_prediction_by_cell_type_" + predictor_name + '_' + error_type
    plt.savefig(RESULTS_PATH + label_save + ".svg") 
    
    
    
def fig4g_to_excel(predictor='XGBoost', error_type='diff_std'):
    
    ''' Exports data used in fig4g to excel file for further processing '''
    
    
    slist = [2,3,6,7]

    mouse_list_by_condition = [np.arange(4), np.arange(4, 8)]

     
    mat_filename = "Position_prediction_key_sessions_" + predictor + "_" + error_type
    
    mat_path = RESULTS_PATH + mat_filename + ".mat"
    error_array = scipy.io.loadmat(mat_path)['error_array']
    
    print(np.around(error_array, decimals=1))
    
    axon_type_label = ['iD', 'DD']
    session_label = ['Day 1', 'Day 2', 'Day 3', 'Day 4']
    
    #Create csv to analyze contradictions later
    with open(RESULTS_PATH + mat_filename + ".csv", mode='w', newline='') as results_csv:
        csvw = csv.writer(results_csv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csvw.writerow(['Axon Type', 'Day', 'Error (' + error_type + ')'])
        csvw.writerow([])
        
        for cond_idx, mlist in enumerate(mouse_list_by_condition):
            error_array_by_axon_type = error_array[mlist]
            
            for sidx, snum in enumerate(slist):
                errors = error_array_by_axon_type[:, snum] #### Assumes that snum is the same as sidx for the results array!!
                csvw.writerow([axon_type_label[cond_idx], session_label[sidx]] + list(errors))
        
            csvw.writerow([])            
            
            
def fig4h_to_excel(predictor='XGBoost', error_type='diff_std'):
    ''' Exports data used in fig4h to excel file for further processing '''
    
    slist = [2,3,6,7]

    mouse_list_by_condition = [np.arange(4), np.arange(4, 8)]

     
    mat_filename = "Position_prediction_by_session_and_cell_type_all_sessions_" + predictor + "_" + error_type
    mat_path = RESULTS_PATH + mat_filename + ".mat"
    error_array = scipy.io.loadmat(mat_path)['error_array']
    
    a = error_array[:, slist]
    # print(np.around(error_array[:, :, 0], decimals =1))
        
    # return    
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