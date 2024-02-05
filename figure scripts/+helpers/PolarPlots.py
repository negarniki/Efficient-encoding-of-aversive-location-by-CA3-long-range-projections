#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:09:20 2020

@author: invisibleobserver
"""
# Polar plots
import numpy as np
import matplotlib.pyplot as plt
#%%
MainPath='/Users/invisibleobserver/IEECR/Files_txt/GrandAverageNew/TimePerRound/Third_Sep20/TimeToImage/2nd/Train/'
mouse=np.loadtxt(MainPath+'all_stats_Tr.txt', dtype='str',comments='#', delimiter=',', unpack=True, usecols=[0])
x=np.loadtxt(MainPath+'all_stats_Tr.txt', comments='#', delimiter=',', unpack=True, usecols=[2])
ax = plt.axes(polar=True)
unit=2*np.pi/150

colors=['r', 'g', 'b']
hatches=['.', '*', '/']
bar1 = plt.bar(25*unit, 1, width=50*unit, bottom=0.0, alpha=0.3, fill=False, hatch='.', linewidth=0)
bar2 = plt.bar(75*unit, 1, width=50*unit, bottom=0.0, alpha=0.3, fill=False, hatch='+', linewidth=0)
bar3 = plt.bar(125*unit, 1, width=50*unit, bottom=0.0, alpha=0.3, fill=False, hatch='/', linewidth=0)
bar4 = plt.bar(50*unit, 1, width=20*unit, bottom=0.0, color='#a80a2f',alpha=0.5, fill=True, linewidth=0)
bar5 = plt.bar(0*unit, 1, width=20*unit, bottom=0.0, color='#129304',alpha=0.5, fill=True, linewidth=0)
for i in range(len(x)):
    ax.quiver(0,0,x[i]*unit,1.0, color='black', angles="xy", scale_units='xy', scale=1. , width=0.007, zorder=6)

ax.set_theta_direction(1)
ax.set_title('Top voted lines detected by Hough transform \n Airpuff Day', size=20)
#ax.set_thetagrids()
ax.set_theta_zero_location('E')
#ax.set_yticklabels([])
ax.set_yticks([])
ax.set_xticks([50*unit, 100*unit, 150*unit])
ax.set_xticklabels(['  50$\,$cm', '  100$\,$cm', '  0$\,$cm \n     150$\,$cm'])
plt.tight_layout()
#plt.savefig(MainPath+'polar_Tr_test.png')
plt.show()
#%% ALl TOGETHR:
MainPath='/Users/invisibleobserver/IEECR/Files_txt/GrandAverageNew/TimePerRound/Third_Sep20/TimeToImage/2nd/Present/'
mouse=np.loadtxt(MainPath+'Train_all_stats_BestCase.txt', dtype='str',comments='#', delimiter=',', unpack=True, usecols=[0])
x1=np.loadtxt(MainPath+'Train_all_stats_BestCase.txt', comments='#', delimiter=',', unpack=True, usecols=[2])
#x2=np.loadtxt(MainPath+'Train_all_stats_BestCase.txt', comments='#', delimiter=',', unpack=True, usecols=[2])
fig = plt.figure(figsize=(8, 10))
ax = fig.add_axes([0.1, 0.04, 0.78, 0.78],projection='polar')
jittered_x1 = x1 - 1.2 * np.random.rand(len(x1)) #for train
#jittered_x1 = x1 + 0.5 * np.random.rand(len(x1))+0.2 #for Post
#jittered_x1 = x1 + 0.5 * np.random.rand(len(x1))+0.5 #for pre
#jittered_x2 = x2 #+ 0.1 * np.random.rand(len(x2)) 

#ax = plt.axes(polar=True)
#ax = plt.axes(polar=True)
unit=2*np.pi/150

colors=['r', 'g', 'b']
hatches=['.', '*', '/']
bar1 = plt.bar(25*unit, 1, width=50*unit, bottom=0.0, alpha=0.3, fill=False, hatch='.', linewidth=0)
bar2 = plt.bar(75*unit, 1, width=50*unit, bottom=0.0, alpha=0.3, fill=False, hatch='+', linewidth=0)
bar3 = plt.bar(125*unit, 1, width=50*unit, bottom=0.0, alpha=0.3, fill=False, hatch='/', linewidth=0)
bar4 = plt.bar(50*unit, 1, width=20*unit, bottom=0.0, color='#a80a2f',alpha=0.5, fill=True, linewidth=0)
bar5 = plt.bar(0*unit, 1, width=20*unit, bottom=0.0, color='#129304',alpha=0.5, fill=True, linewidth=0)
#c='#003f7f'
#c='k'
c='#560303'#red
for i in range(len(x1)):
    if i==1:
        ax.quiver(0,0,jittered_x1[i]*unit,1.0, color=c, angles="xy", scale_units='xy', scale=1. , width=0.005, zorder=6,label='Day 1')
        #ax.quiver(0,0,x2[i]*unit,1.0, color=c, angles="xy", scale_units='xy', scale=1. , width=0.005, zorder=7, label='Day 2')
    else:
        ax.quiver(0,0,jittered_x1[i]*unit,1.0, color=c, angles="xy",scale_units='xy', scale=1. , width=0.005, zorder=6)
       # ax.quiver(0,0,x2[i]*unit,1.0, color=c, angles="xy", scale_units='xy', scale=1. , width=0.005, zorder=7)
#handles, labels = ax.get_legend_handles_labels()

ax.set_theta_direction(1)
ax.set_title('Top voted lines detected by Hough transform \n Train', size=20)
#ax.set_thetagrids()
ax.set_theta_zero_location('E')
#ax.set_yticklabels([])
ax.set_yticks([])
ax.set_ylim([0, 1.005])
ax.set_xticks([50*unit, 100*unit, 150*unit])
ax.set_xticklabels(['  50$\,$cm', '  100$\,$cm', '  0$\,$cm \n     150$\,$cm'])
#plt.legend()
#fig.legend(loc='lower left', ncol=1)
#plt.tight_layout()
#plt.savefig(MainPath+'polar_Train_all_BestCase_Jit100.pdf', dpi=100)
plt.show()
#%% ALl TOGETHR CLOCKWISE
MainPath='/Users/invisibleobserver/IEECR/Files_txt/GrandAverageNew/TimePerRound/Third_Sep20/TimeToImage/2nd/FirstTrainSess/'
mouse=np.loadtxt(MainPath+'FirstTrainSess_all_stats_automated.txt', dtype='str',comments='#', delimiter=',', unpack=True, usecols=[0])
x1=np.loadtxt(MainPath+'FirstTrainSess_all_stats_automated.txt', comments='#', delimiter=',', unpack=True, usecols=[2])
#x1=np.loadtxt(MainPath+'PreTrain_d1_BestCase.txt', comments='#', delimiter=',', unpack=True, usecols=[2])
#x2=np.loadtxt(MainPath+'PostTrain_d2_BestCase.txt', comments='#', delimiter=',', unpack=True, usecols=[2])
fig = plt.figure(figsize=(8, 10))
ax = fig.add_axes([0.1, 0.04, 0.78, 0.78],projection='polar')

jittered_x1 = x1 - 1.2 * np.random.rand(len(x1)) #for train
#jittered_x1 = x1 + 0.5 * np.random.rand(len(x1))+0.2 #for Post
#jittered_x1 = x1 + 0.5 * np.random.rand(len(x1))+0.5 #for pre
#jittered_x2 = x2 #+ 0.1 * np.random.rand(len(x2)) 

#ax = plt.axes(polar=True)
#ax = plt.axes(polar=True)
unit=2*np.pi/150 

colors=['r', 'g', 'b']
hatches=['.', '*', '/']
bar1 = plt.bar(25*unit, 1, width=50*unit, bottom=0.0, alpha=0.3, fill=False, hatch='.', linewidth=0)
bar2 = plt.bar(75*unit, 1, width=50*unit, bottom=0.0, alpha=0.3, fill=False, hatch='+', linewidth=0)
bar3 = plt.bar(125*unit, 1, width=50*unit, bottom=0.0, alpha=0.3, fill=False, hatch='/', linewidth=0)
bar4 = plt.bar(50*unit, 1, width=20*unit, bottom=0.0, color='#a80a2f',alpha=0.5, fill=True, linewidth=0)
bar5 = plt.bar(0*unit, 1, width=20*unit, bottom=0.0, color='#129304',alpha=0.5, fill=True, linewidth=0)
#c='#003f7f'
c='k'
#c='#560303'#red
for i in range(len(x1)):
    if i==1:
        ax.quiver(0,0,jittered_x1[i]*unit,1.0, color=c, angles="xy", scale_units='xy', scale=1. , alpha=0.99, width=0.005, zorder=6,label='Day 1')
        #ax.quiver(0,0,x2[i]*unit,1.0, color=c, angles="xy", scale_units='xy', scale=1. , width=0.005, zorder=7, label='Day 2')
    else:
        ax.quiver(0,0,jittered_x1[i]*unit,1.0, color=c, angles="xy",scale_units='xy', scale=1 , alpha=0.99, width=0.005, zorder=6)
        #ax.quiver(0,0,x2[i]*unit,1.0, color=c, angles="xy", scale_units='xy', scale=1. , width=0.005, zorder=7)
#handles, labels = ax.get_legend_handles_labels()

ax.set_theta_direction(1)
ax.set_title('Top voted lines detected by Hough transform \n First Train Day', size=20)
#ax.set_thetagrids()
ax.set_theta_zero_location('W')
ax.set_theta_direction(-1)# minus sign make everything clockwise
#ax.set_yticklabels([])
ax.set_yticks([])
ax.set_ylim([0, 1.005])
ax.set_xticks([50*unit, 100*unit, 150*unit])
ax.set_xticklabels([ '  50$\,$cm', '  100$\,$cm','0$\,$cm         \n150$\,$cm          '], fontsize=10)
#plt.legend()
#fig.legend(loc='lower left', ncol=1)
#plt.tight_layout()
plt.savefig(MainPath+'polar_FirstTrainSess_Jit100_clock.pdf', dpi=100)
plt.show()
