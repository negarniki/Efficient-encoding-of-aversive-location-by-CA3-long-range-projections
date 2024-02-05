#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 16:21:44 2020

@author:
"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import os, fnmatch
import subprocess
from scipy import stats as st
from PIL import Image
import cv2 as cv
from scipy.signal import find_peaks
import glob, re


cv.destroyAllWindows()
MiceList=['M278']#['M259','M261','M262','M263','M270','M271','M272','M278']
MainPath='/Users/invisibleobserver/IEECR/Files_txt/GrandAverageNew/TimePerRound/Third_Sep20/TimeToImage/2nd/M270_Eval/'#/2nd/Train/
cv.destroyAllWindows()

listSession= ['Train']#['PreTrain_d1', 'PreTrain_d2', 'Train', 'PostTrain_d1' , 'PostTrain_d2'] #'PreTrain_d2']#,    
#SAVEPath=MainPath+'OneThresh/Train/'
SAVEPath=MainPath
for s, session in enumerate(listSession):
    fileTotWrite=SAVEPath+session+'_all_stats_automated.txt'
    subprocess.call(['touch', fileTotWrite]) # creat the file to write the fit results in
    f_tot=open (fileTotWrite, 'r+')
    f_tot.write('# mouse, session, top line location (cm)), accum_thresh\n')
    subprocess.call(['touch', SAVEPath+session+'_noise_stats.txt'])
    f1=open(SAVEPath+session+'_noise_stats.txt', 'r+')
    f1.write('# mouse, session, ROI_mean (No AP, No Rew), ROI_std(No AP, No Rew), ROI_AP_mean, ROI_AP_std, precentage of pixels with values greater than 200\n')
    for mouse in MiceList:
        print (mouse)
       # plt.figure()
        #ESM=mouse+'/%s/'%session
        # Start analysing the photo using OpenCV
        orig=cv.imread(MainPath+mouse+'_'+session+'.png', cv.IMREAD_GRAYSCALE)
        NoRew_img=cv.imread(MainPath+mouse+'_'+session+'_NoRew.png', cv.IMREAD_GRAYSCALE)
        #cv.imshow('Original No Reward', NoRew_img)
        ROI1, ROI2, ROI3= NoRew_img[0::, 8:28], NoRew_img[0::, 60:80], NoRew_img[0::, 100:120]
        ROIs=np.concatenate((ROI1,ROI2,ROI3), axis=0)
        ROI_mean, ROI_std=np.nanmean(ROIs), np.nanstd(ROIs)
        ROI_AP=NoRew_img[0::, 40:60]
        #AP_mean, AP_std=np.nanmean(NoRew_img[np.where(ROI_AP)]), np.nanstd(NoRew_img[np.where(ROI_AP)])
        AP_mean, AP_std=np.nanmean(ROI_AP), np.nanstd(ROI_AP)
        #np.int(np.percentile(cv.bitwise_not(NoRew_img), 85))
        # This part will make sure that if the mouse is running faster in the AP zone, we will detect that as well.
        if AP_mean > (3*ROI_mean):
        #if mouse=='M262' and session=='PostTrain_d2': # If the mouse is running faster in the AP zone, we do not need to create a bitwise image
            thr_val=np.percentile(np.ndarray.flatten(NoRew_img), 87)    
            ret, thresh = cv.threshold(np.uint8(NoRew_img), thr_val, 255, cv.THRESH_BINARY)
            #print(ret)
            point_num1=sum(map(lambda x: x> 190, np.ndarray.flatten(NoRew_img)))/len(np.ndarray.flatten(NoRew_img))
            point_num2=sum(map(lambda x: x> 198, np.ndarray.flatten(NoRew_img)))/len(np.ndarray.flatten(NoRew_img))
            point_num3=sum(map(lambda x: x> 200, np.ndarray.flatten(NoRew_img)))/len(np.ndarray.flatten(NoRew_img))
            point_num4=sum(map(lambda x: x> 220, np.ndarray.flatten(NoRew_img)))/len(np.ndarray.flatten(NoRew_img))
            prec1=np.percentile(np.ndarray.flatten(NoRew_img), 85)
            # plt.hist(np.ndarray.flatten(NoRew_img), bins=100)
            #plt.title(mouse+' '+session)
            #plt.xlabel('Pixel Values')
            #plt.ylabel('Number')
            #plt.savefig(SAVEPath+mouse+'_'+session+'Hist.png')
            print('Higher speed in AP zone')
            cv.imshow('Tresh', thresh)
        else:
            bit=cv.bitwise_not(NoRew_img)
            thr_val=np.percentile(np.ndarray.flatten(bit), 85)
            ret, thresh = cv.threshold(np.uint8(bit), thr_val, 255, cv.THRESH_BINARY)
            #thr2 = cv.adaptiveThreshold(bit,255,cv.ADAPTIVE_THRESH_MEAN_C,cv.THRESH_BINARY,103, 10)
            #thr3 = cv.adaptiveThreshold(bit,255,cv.ADAPTIVE_THRESH_GAUSSIAN_C,cv.THRESH_BINARY,103,10)
            #print(ret)
            #cv.imshow('Tresh 2', thr2)
            #cv.imshow('Tresh 3', thr3)
            #cv.imshow('Tresh', thresh)
            point_num1=sum(map(lambda x: x> 190, np.ndarray.flatten(bit)))/len(np.ndarray.flatten(bit))
            point_num2=sum(map(lambda x: x> 198, np.ndarray.flatten(bit)))/len(np.ndarray.flatten(bit))
            point_num3=sum(map(lambda x: x> 200, np.ndarray.flatten(bit)))/len(np.ndarray.flatten(bit))
            point_num4=sum(map(lambda x: x> 220, np.ndarray.flatten(bit)))/len(np.ndarray.flatten(bit))
            prec1=np.percentile(np.ndarray.flatten(bit), 85)
            #plt.hist(np.ndarray.flatten(bit), bins=100)
            #plt.title(mouse+' '+session)
            #plt.xlabel('Pixel Values')
            #plt.ylabel('Number')
            #plt.savefig(SAVEPath+mouse+'_'+session+'Hist.png')
            
        fileWrite=SAVEPath+mouse+'_'+session+'_thresh_top1.txt'
        subprocess.call(['rm','-r', fileWrite]) # delete the file if it exists
        subprocess.call(['touch', fileWrite]) # creat the file to write the fit results in
       
        # Find vertical lines (thetha=90)     img.shape[0] is the number of rounds 
        accum_thresh=0
        while accum_thresh>-1:
            accum_thresh=accum_thresh+1
            lines_hogh = cv.HoughLinesP(thresh,rho = 1,theta = 90*np.pi/90,threshold = accum_thresh, maxLineGap=NoRew_img.shape[0]-1, minLineLength =NoRew_img.shape[0]/6)
            if lines_hogh is None:
                print('NONE')
                break
        lines_hogh = cv.HoughLinesP(thresh,rho = 1,theta = 90*np.pi/90,threshold = accum_thresh-1, maxLineGap=NoRew_img.shape[0]-1, minLineLength =NoRew_img.shape[0]/6)
        
        # result = cv.cvtColor(NoRew_img, cv.COLOR_GRAY2BGR)
        f=open (fileWrite, 'r+')
        
        for i in range(len(lines_hogh[:,:])):
            x1,y1,x2,y2=lines_hogh[i,0]
            print(x1+10,y1,x2+10,y2, accum_thresh-1)      
           # cv.line(result,(x1+10,y1),(x2+10,y2),(255,86,180),1)
            f.write('%s,%i,%i,%i,%i,%i \n'%(session, x1+10, y1, x2+10, y2, accum_thresh-1))
            f_tot.write('%s,%s,%i,%i\n'%(mouse, session, x1+10, accum_thresh-1))
        f.close()
        #cv.imshow('result', result)
        #cv.imshow('thresh', thresh)
        #cv.imshow('bit', bit)
        #cv.imwrite(SAVEPath+mouse+'_'+session+'houghlines_top3.jpg',result)
        # Use standard deviation in three different winsows for noise estimation
       
        #f1=open(SAVEPath+'noise_stats.txt', 'r+')
        f1.write('%s, %s, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, , %.2f, %0.2f\n'%(mouse, session, ROI_mean, ROI_std, AP_mean, AP_std,point_num1, point_num2, point_num3, point_num4, prec1, np.median(scipy.stats.mode(bit)[0])))
        # len(np.ndarray.flatten(bit)[np.where(np.ndarray.flatten(bit)>220)])/len(np.ndarray.flatten(bit))*100, np.int(np.mean(bit[np.where(ROI_AP)]))
        #f1.close()
    f_tot.close()
    f1.close()
        #cv.destroyAllWindows()
        #%%
#MainPath='/Users/invisibleobserver/IEECR/Files_txt/GrandAverageNew/TimePerRound/Third_Sep20/TimeToImage/2nd/Post/'
MiceList=['M278']#['M259','M261','M262','M263','M270','M271','M272']
#mouse=np.loadtxt(MainPath+'all_stats_Tr.txt', dtype='str',comments='#', delimiter=',', unpack=True, usecols=[0])
#x=np.loadtxt(MainPath+'all_stats_Tr.txt', comments='#', delimiter=',', unpack=True, usecols=[2])
SAVEPath='/Users/invisibleobserver/IEECR/Files_txt/GrandAverageNew/TimePerRound/Third_Sep20/TimeToImage/2nd/FirstTrainSess/Automated/'

listSession= ['FirstTrainSess']# , 'PostTrain_d2']
for s, session in enumerate(listSession):
    for mouse in MiceList:
        print (mouse)
       
        # Figures with 3 axes:
        csfont = {'fontname':'Open Sans'}
        # define name of session for plotting: 
        if session=='PreTrain_d1':
            prop_ses='Habituation day$\,$1'
        elif session=='PreTrain_d2':
            prop_ses='Habituation day$\,$2'
        elif session=='Train':
            prop_ses='Air-Puff day'
        elif session=='PostTrain_d1':
            prop_ses='Retrieval day$\,$1'
        elif session=='PostTrain_d2':
            prop_ses='Retrieval day$\,$2 '
        elif session=='FirstTrainSess':
            prop_ses='Train First Day'
        print(prop_ses)
        
       # FigLineList=glob.glob(SAVEPath+'*%s*%shoughlines*.jpg'%(mouse, session))
        #img_orig=cv.imread(MainPath+mouse+'_'+session+'.png', cv.IMREAD_GRAYSCALE)
        TextList=glob.glob(MainPath+mouse+'_'+session+'_*.txt')
        #TextList=glob.glob(MainPath+'Automated/'+mouse+'_'+session+'_*.txt')
        for t in enumerate(TextList):
            fig, axes = plt.subplots(3,1, sharex=False)
            num_lines = sum(1 for line in open(TextList[t[0]]))

            x1, y1, x2, y2, accum_thresh =np.loadtxt(TextList[t[0]], comments='#', delimiter=',', unpack=True, usecols=[1,2,3,4,5])
            #for p in enumerate(FigLineList):
            top=re.findall("\d+", TextList[t[0]])[-1] # Find the number of lines (e.g. top3 lines)
            #SaveName='%s_%s_top_report_%s'%(mouse, session, top)
            img_orig=cv.imread(MainPath+mouse+'_'+session+'.png', cv.IMREAD_GRAYSCALE)
            NoRew_img=cv.imread(MainPath+mouse+'_'+session+'_NoRew.png', cv.IMREAD_GRAYSCALE)
            bit2=cv.bitwise_not(img_orig)
            thr_val=np.percentile(np.ndarray.flatten(bit2), 85)
            _, thresh = cv.threshold(np.uint8(bit2), thr_val, 255, cv.THRESH_BINARY)
            #img_line=cv.imread()
            #img_line = cv.cvtColor(img_line, cv.COLOR_BGR2RGB)
                
            plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9, wspace=0.1)
            axes[0].imshow(img_orig, aspect = "auto",cmap=plt.cm.gray)
            axes[1].imshow(thresh, aspect = "auto", cmap=plt.cm.gray)#, extent=(np.rad2deg(angles[-1]), np.rad2deg(angles[0]), dists[-1], dists[0]))
            axes[2].imshow(img_orig, aspect = "auto", cmap=plt.cm.gray)
            # Title and axis labels
            axes[0].set_title('%s: %s'%(mouse, prop_ses), size=25, color='k', **csfont)
            axes[0].set_xticks(np.linspace(0, 150, 16))
            axes[0].set_ylabel('Rounds', size=13, **csfont)
            
            axes[1].set_xticks(np.linspace(0, 150, 16))
            axes[1].set_ylabel('Rounds', size=13, **csfont)
            
            axes[2].set_xlabel('Distance on Belt [cm]', size=13, **csfont)
            axes[2].set_xticks(np.linspace(0, 150, 16))
            axes[2].set_ylabel('Rounds', size=13, **csfont)
            # Plot the detected lines from OpenCv Hough line detection algorithm
            for i in range(num_lines):
                axes[2].plot((x1,x2),(y1,y2), color='#ffcc00', lw=2)
                # Show the airpuff zone
            axes[0].axvline(x=40,color='#ed0e02',linestyle='--', alpha=0.8)
            axes[0].axvline(x=60,color='#ed0e02', linestyle='--', alpha=0.8)
            axes[1].axvline(x=40,color='#ed0e02',linestyle='--', alpha=0.8)
            axes[1].axvline(x=60,color='#ed0e02', linestyle='--', alpha=0.8)
            axes[2].axvline(x=40,color='#ed0e02',linestyle='--', alpha=0.8)
            axes[2].axvline(x=60,color='#ed0e02', linestyle='--', alpha=0.8)
            # Show the arrow for air-puff zone marking
            bbox_props = dict(boxstyle="darrow,pad=0.3", fc="#ed0e02", ec="k", alpha=0.3, lw=2)
            axes[0].text(42, img_orig.shape[0]+img_orig.shape[0]/2, "Air-Puff Zone", ha="left", va="bottom", rotation=0,size=7,bbox=bbox_props, **csfont)
            
            #plt.tight_layout()
            plt.show()
            #plt.savefig(SAVEPath+'%s'+'.pdf'%SaveName, dpi=100)
            plt.savefig(SAVEPath+mouse+session+'.pdf', dpi=100)
            #plt.close(fig)
#%%
        alpha = 4.2 # Simple contrast control
        beta = 50    # Simple brightness control
        gamma=0.1        
        alpha_image = np.zeros(NoRew_img.shape, NoRew_img.dtype)
        #try:
         #   alpha = float(input('* Enter the alpha value [1.0-3.0]: '))
          #  beta = int(input('* Enter the beta value [0-100]: '))
        #except ValueError:
         #       print('Error, not a number')
        alpha_image = cv.convertScaleAbs(NoRew_img, alpha=alpha, beta=beta)
        cv.imshow('Alpha Beta Image', alpha_image)
        #gamma
        lookUpTable = np.empty((1,256), np.uint8)
        for i in range(256):
            lookUpTable[0,i] = np.clip(pow(i / 255.0, gamma) * 255.0, 0, 255)
        gamma_image = cv.LUT(NoRew_img, lookUpTable)
        cv.imshow('Gamma', gamma_image )
        #Change the colorscale so thresholding can work (we need to put values above a threshold to max)
        bit=cv.bitwise_not(NoRew_img)
        cv.imshow('BitWise', bit)
        # Thresholding the image # Estimate a threshold level: 5*std
      #  _, thresh1 = cv.threshold(np.uint8(bit), 220, 245, cv.THRESH_TRUNC)
      #  cv.imshow('First Thresholding Image',thresh1)
        _, thresh = cv.threshold(np.uint8(bit),35, 255, cv.THRESH_BINARY)
        cv.imshow('Thresholding Image',thresh)
