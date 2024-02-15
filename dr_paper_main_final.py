# -*- coding: utf-8 -*-
"""
For "Multiplicative Joint Coding in Preparatory Activity for Reaching Sequence 
     in Macaque Motor Cortex".

@author: yc, Cui Lab
"""

#%% Setup

import os
import operator
import torch
from torch import nn
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import tqdm
from scipy import optimize
# from scipy import integrate
import seaborn as sns
from sklearn.linear_model import LinearRegression
from scipy.optimize import leastsq
from sklearn.decomposition import PCA
import scipy.stats as stats
from scipy.spatial.distance import pdist, squareform
import similarity as sim
import random
from matplotlib.colors import ListedColormap


PATH_SAVE = 'F://Aac//Output//Script//M1_DR//toNC//code//final_results//'
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


#%% Define RNN cell

class rRNNCell(nn.Module):
    
    global TAU, g
    TAU = 50
    g = 1.5
    
    def __init__(self, input_size, hidden_size, output_size):

        super(rRNNCell, self).__init__()

        self.w_ih = nn.Linear(input_size, hidden_size)
        self.w_hh = nn.Linear(hidden_size, hidden_size)
        self.w_ho = nn.Linear(hidden_size, output_size)

        nn.init.normal_(self.w_hh.weight, mean=0, std=g/math.sqrt(hidden_size))
        nn.init.uniform_(self.w_ih.weight, a=-g/math.sqrt(hidden_size)/1000,
                        b=g/math.sqrt(hidden_size)/1000)
        nn.init.uniform_(self.w_ho.weight, a=-g/math.sqrt(hidden_size)/1000,
                        b=g/math.sqrt(hidden_size)/1000)

    def forward(self, input_now, hidden_last):

        relu_filter = nn.ReLU()
        r_last = relu_filter(torch.tanh(hidden_last))

        # dynamic evolving rule
        h_now = (1 - 1/TAU)*hidden_last + (self.w_ih(input_now) + self.w_hh(r_last))/TAU

        r_now = relu_filter(torch.tanh(h_now))
        output_now = self.w_ho(r_now)

        return h_now, r_now, output_now


class rRNN(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):

        super(rRNN, self).__init__()
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.output_size = output_size

        self.cell = rRNNCell(input_size, hidden_size, output_size)


    def forward(self, input_series, initial_states):

        time_steps = input_series.size(1)
        #(batch, time, neuron)

        h = initial_states

        h_series = []
        r_series = []
        z_series = []

        for t in range(time_steps):
            h, r, z = self.cell(input_series[:, t, :], h)
            h_series.append(h)
            r_series.append(r)
            z_series.append(z)

        return (torch.stack(h_series, dim=1),
                torch.stack(r_series, dim=1),
                torch.stack(z_series, dim=1))
    
    
#%% Define task-relevant variables

def gaussian_f(mu, sigma, x):
    # gaussian function
    return (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-(x-mu)**2/2/sigma**2)


# calculate points in a Gaussian distribution from idt to idt2, idt at peak=1
gaussian_p = lambda t, idt, idt2, sigma: gaussian_f(
        t[idt], sigma, t[min([idt, idt2]):max([idt, idt2])])/max(
        gaussian_f(t[idt], sigma, t[min([idt, idt2]):max([idt, idt2])]))


def calculate_best_sigma_for_gaussian(ts, tp1, tp2):
    """
    Calculate the sigma of a Gaussian distribution covering the given period.

    Parameters
    ----------
    ts: array, time points.
    tp1: int, time point 1 (start point)
    tp2: int, time point 2 (end point)

    Returns
    -------
    best_sigma: float, sigma optimized for the given period

    """

    cost_function = lambda v: abs(min(gaussian_p(ts, tp1, tp2, v)) - 1e-3)

    res = optimize.minimize(cost_function, (20), method='Powell')
    best_sigma = res.x

    return best_sigma    


def gen_input(target1, target2, form, factor=1):
    """
    Generate a condition DataFrame including 2nd target
    in visual coordinates (Dir2) or in movement coordinates (MovDir2).

    Parameters
    ----------
    target1: list, target direction of 1st reach in degrees.
    target2: list, target direction of 2nd reach in degrees (absolutely).
    form: str, coordinates for 2nd reach, 'actual' or 'relative'
    factor: float, optional, amplitude of inputs for 2nd reach.

    Returns
    -------
    x: tensor, 5 inputs including two 2d target coordinates and a Go-pulse.
    cond_df: pd.DataFrame, a detailed condition table

    """

    nTrial = len(target1)

    # signal duration
    h1 = np.zeros((1, nTime))
    h1[0, TO:GO] = 1  # assuming signals maintained by other brain areas

    # get condition dataframe
    cond_df = pd.DataFrame(columns=['Dir1', 'Dir2', 'TypeIdx', 'Type', 'CondName'])
    for iTrial in range(nTrial):
        cond_df.loc[iTrial, 'Dir1'] = target1[iTrial]
        if target2[iTrial] == -1:
            cond_df.loc[iTrial, 'Dir2'] = 0
            cond_df.loc[iTrial, 'TypeIdx'] = 0
            cond_df.loc[iTrial, 'Type'] = 'SR'
            cond_df.loc[iTrial, 'CondName'] = 'SR'
        else:
            cond_df.loc[iTrial, 'Dir2'] = target2[iTrial]
            if np.mod(target2[iTrial] - target1[iTrial], 360) < 180:
                cond_df.loc[iTrial, 'TypeIdx'] = 1
                cond_df.loc[iTrial, 'Type'] = 'CCW'
                cond_df.loc[iTrial, 'CondName'] = 'CCW%d' % np.mod(target2[iTrial] - target1[iTrial], 360) 
            elif np.mod(target2[iTrial] - target1[iTrial], 360) > 180:
                cond_df.loc[iTrial, 'TypeIdx'] = 2
                cond_df.loc[iTrial, 'Type'] = 'CW'
                cond_df.loc[iTrial, 'CondName'] = 'CW%d' % np.mod(target1[iTrial] - target2[iTrial], 360) 
            else:
                cond_df.loc[iTrial, 'TypeIdx'] = 1
                cond_df.loc[iTrial, 'Type'] = 'CCW'
                cond_df.loc[iTrial, 'CondName'] = str(np.mod(target1[iTrial] - target2[iTrial], 360))

    # calculate relative 2nd reach direction (based on 1st reach's endpoint)
    cond_df['MovDir2'] = -1
    for iTrial in range(nTrial):
        dir1 = cond_df.loc[iTrial, 'Dir1']
        dir2 = cond_df.loc[iTrial, 'Dir2']
        if dir2 - cond_df.loc[iTrial, 'Dir1'] < 0:
            dir2 = dir2 + 360
        if cond_df.loc[iTrial, 'Type'] != 'SR':
            cond_df.loc[iTrial, 'MovDir2'] = (dir1+180-(180-(dir2-dir1))/2)

    cond_df.loc[cond_df['TypeIdx']>0, 'MovDir2'] = np.mod(
                        cond_df.loc[cond_df['TypeIdx']>0, 'MovDir2'], 360)
    # print(cond_df)

    # assemble inputs
    if form == 'actual':
        targ1_cos = np.cos(np.array(cond_df['Dir1']).astype(float)/180*np.pi)
        targ1_sin = np.sin(np.array(cond_df['Dir1']).astype(float)/180*np.pi)
        targ2_cos = np.cos(np.array(cond_df['Dir2']).astype(float)/180*np.pi)
        targ2_sin = np.sin(np.array(cond_df['Dir2']).astype(float)/180*np.pi)

    elif form == 'relative':
        targ1_cos = np.cos(np.array(cond_df['Dir1']).astype(float)/180*np.pi)
        targ1_sin = np.sin(np.array(cond_df['Dir1']).astype(float)/180*np.pi)
        targ2a_cos = np.cos(np.array(cond_df['Dir2']).astype(float)/180*np.pi)
        targ2a_sin = np.sin(np.array(cond_df['Dir2']).astype(float)/180*np.pi)

        dis = np.sqrt((targ2a_cos-targ1_cos)**2+(targ2a_sin-targ1_sin)**2)
        targ2_cos = np.cos(np.array(cond_df['MovDir2']).astype(float)/180*np.pi)*dis
        targ2_sin = np.sin(np.array(cond_df['MovDir2']).astype(float)/180*np.pi)*dis

    # GO pulse
    sigma_for_go = calculate_best_sigma_for_gaussian(timepoints, GO, GO+10)
    go_pulse = np.zeros((1, nTime))
    go_pulse[0, GO:(GO+10)] = gaussian_p(timepoints, GO, GO+10, sigma_for_go)
    go_pulse[0, (GO-10):GO] = gaussian_p(timepoints, GO, GO-10, sigma_for_go)
    
    ip1 = np.kron(h1.T, targ1_cos)
    ip2 = np.kron(h1.T, targ1_sin)
    ip3 = np.kron(h1.T, targ2_cos)*factor
    ip4 = np.kron(h1.T, targ2_sin)*factor
    ip5 = np.tile(go_pulse.T, (1, ip1.shape[1]))

    ip3[:, cond_df['Type']=='SR'] = 0
    ip4[:, cond_df['Type']=='SR'] = 0

    x_time = np.stack([ip1, ip2, ip3, ip4, ip5]).T

    # numpy to torch.tensor
    x = torch.from_numpy(np.array(x_time, dtype=np.float32))

    return x, cond_df


def gen_PV_as_target(target1_list, target21_list, overlap1=0, overlap2=0):
    """
    Generate desired population vectors.

    Parameters
    ----------
    target1_list: array, list of 1st reach directions, in degrees.
    target21_list: array, list of 2nd reach directions, in degrees.
    overlap1: int, PV1 plan ends with respect to TT, default=0
    overlap2: int, PV2 plan starts with respect to MO2, default=0

    Returns
    -------
    y_time: array, assembled desired 2d population vectors
    """

    # Gaussian segments to assemble PV
    sigma11 = calculate_best_sigma_for_gaussian(timepoints, GO, MO)
    sigma12 = calculate_best_sigma_for_gaussian(timepoints, MO, TT+overlap1)
    sigma21 = calculate_best_sigma_for_gaussian(timepoints, MO2-50-overlap2, MO2)
    sigma22 = calculate_best_sigma_for_gaussian(timepoints, MO2, TT2)

    # PV1
    pv1 = np.zeros((1, nTime))
    pv1[0, GO:MO] = gaussian_p(timepoints, MO, GO, sigma11)
    pv1[0, MO:(TT+overlap1)] = gaussian_p(timepoints, MO, TT+overlap1, sigma12)

    theta1_cos = np.cos(target1_list/180*np.pi)
    theta1_sin = np.sin(target1_list/180*np.pi)

    pv1_cos = np.kron(pv1.T, theta1_cos)
    pv1_sin = np.kron(pv1.T, theta1_sin)

    # PV2
    pv2 = np.zeros((1, nTime))
    pv2[0, (MO2-50-overlap2):MO2] = gaussian_p(timepoints, MO2,
                                               MO2-50-overlap2, sigma21)
    pv2[0, MO2:TT2] = gaussian_p(timepoints, MO2, TT2, sigma22)

    theta2_cos = np.cos(target21_list/180*np.pi)
    theta2_sin = np.sin(target21_list/180*np.pi)

    theta2_cos[target21_list==-1] = 0
    theta2_sin[target21_list==-1] = 0

    pv2_cos = np.kron(pv2.T, theta2_cos)
    pv2_sin = np.kron(pv2.T, theta2_sin)

    # sum PV1 and PV2
    cos_sum = pv1_cos + pv2_cos
    sin_sum = pv1_sin + pv2_sin

    y_time = np.stack([cos_sum, sin_sum]).T*5 #amplified to match the magnitude

    return y_time


def gen_target(cond_df, *args):
    """
    Generate desired outputs.

    Parameters
    ----------
    cond_df: pd.DataFrame, condition DataFrame.

    Returns
    -------
    y_time: array, assembled desired 2d population vectors
    """
    targ1_theta = np.array(cond_df['Dir1']).astype(float)
    targ2_theta = np.array(cond_df['MovDir2']).astype(float)

    if len(args)>0:
      overlap1, overlap2 = args
      y_time = gen_PV_as_target(targ1_theta, targ2_theta, overlap1, overlap2)
    else:
      y_time = gen_PV_as_target(targ1_theta, targ2_theta)
    # torch
    y = torch.from_numpy(np.array(y_time, dtype=np.float32))

    return y


def gen_input_target_specific(targ1, targ2, input_form, **kws):
    # kws = {'factor':1, 'overlaps':(0, 0)}

    x_tensor, cond_df = gen_input(targ1, targ2, input_form, factor=kws['factor'])
    y_tensor = gen_target(cond_df, kws['overlap1'], kws['overlap2'])

    return x_tensor, y_tensor, cond_df


def gen_2tinput(target1, target2, form, factor=1):
    """
    Generate a condition DataFrame including 2nd target
    in visual coordinates (Dir2) or in movement coordinates (MovDir2),
    but with 2 triggers (GO-pulses).
    This differs with gen_input only in Go-pulse.

    Parameters
    ----------
    target1: list, target direction of 1st reach in degrees.
    target2: list, target direction of 2nd reach in degrees (absolutely).
    form: str, coordinates for 2nd reach, 'actual' or 'relative'
    factor: float, optional, amplitude of inputs for 2nd reach.

    Returns
    -------
    x: tensor, 5 inputs including two 2d target coordinates and Go pulse.
    cond_df: pd.DataFrame, a detailed condition table

    """

    nTrial = len(target1)

    # signal duration
    h1 = np.zeros((1, nTime))
    h1[0, TO:GO] = 1  # assuming signals maintained by other brain areas

    # get condition dataframe
    cond_df = pd.DataFrame(columns=['Dir1', 'Dir2', 'TypeIdx', 'Type', 'CondName'])
    for iTrial in range(nTrial):
        cond_df.loc[iTrial, 'Dir1'] = target1[iTrial]
        if target2[iTrial] == -1:
            cond_df.loc[iTrial, 'Dir2'] = 0
            cond_df.loc[iTrial, 'TypeIdx'] = 0
            cond_df.loc[iTrial, 'Type'] = 'SR'
            cond_df.loc[iTrial, 'CondName'] = 'SR'
        else:
            cond_df.loc[iTrial, 'Dir2'] = target2[iTrial]
            if np.mod(target2[iTrial] - target1[iTrial], 360) < 180:
                cond_df.loc[iTrial, 'TypeIdx'] = 1
                cond_df.loc[iTrial, 'Type'] = 'CCW'
                cond_df.loc[iTrial, 'CondName'] = 'CCW%d' % np.mod(target2[iTrial] - target1[iTrial], 360) 
            elif np.mod(target2[iTrial] - target1[iTrial], 360) > 180:
                cond_df.loc[iTrial, 'TypeIdx'] = 2
                cond_df.loc[iTrial, 'Type'] = 'CW'
                cond_df.loc[iTrial, 'CondName'] = 'CW%d' % np.mod(target1[iTrial] - target2[iTrial], 360) 
            else:
                cond_df.loc[iTrial, 'TypeIdx'] = 1
                cond_df.loc[iTrial, 'Type'] = 'CCW'
                cond_df.loc[iTrial, 'CondName'] = str(np.mod(target1[iTrial] - target2[iTrial], 360))

    # calculate relative 2nd reach direction (based on 1st reach's endpoint)
    cond_df['MovDir2'] = -1
    for iTrial in range(nTrial):
        dir1 = cond_df.loc[iTrial, 'Dir1']
        dir2 = cond_df.loc[iTrial, 'Dir2']
        if dir2 - cond_df.loc[iTrial, 'Dir1'] < 0:
            dir2 = dir2 + 360
        if cond_df.loc[iTrial, 'Type'] != 'SR':
            cond_df.loc[iTrial, 'MovDir2'] = (dir1+180-(180-(dir2-dir1))/2)

    cond_df.loc[cond_df['TypeIdx']>0, 'MovDir2'] = np.mod(
                        cond_df.loc[cond_df['TypeIdx']>0, 'MovDir2'], 360)
    # print(cond_df)

    # assemble inputs
    if form == 'actual':
        targ1_cos = np.cos(np.array(cond_df['Dir1']).astype(float)/180*np.pi)
        targ1_sin = np.sin(np.array(cond_df['Dir1']).astype(float)/180*np.pi)
        targ2_cos = np.cos(np.array(cond_df['Dir2']).astype(float)/180*np.pi)
        targ2_sin = np.sin(np.array(cond_df['Dir2']).astype(float)/180*np.pi)

    elif form == 'relative':
        targ1_cos = np.cos(np.array(cond_df['Dir1']).astype(float)/180*np.pi)
        targ1_sin = np.sin(np.array(cond_df['Dir1']).astype(float)/180*np.pi)
        targ2a_cos = np.cos(np.array(cond_df['Dir2']).astype(float)/180*np.pi)
        targ2a_sin = np.sin(np.array(cond_df['Dir2']).astype(float)/180*np.pi)

        dis = np.sqrt((targ2a_cos-targ1_cos)**2+(targ2a_sin-targ1_sin)**2)
        targ2_cos = np.cos(np.array(cond_df['MovDir2']).astype(float)/180*np.pi)*dis
        targ2_sin = np.sin(np.array(cond_df['MovDir2']).astype(float)/180*np.pi)*dis

    # GO pulse
    sigma_for_go = calculate_best_sigma_for_gaussian(timepoints, GO, GO+10)
    go_pulse = np.zeros((1, nTime))
    go_pulse[0, GO:(GO+10)] = gaussian_p(timepoints, GO, GO+10, sigma_for_go)
    go_pulse[0, (GO-10):GO] = gaussian_p(timepoints, GO, GO-10, sigma_for_go)
   
    go_pulse2 = np.zeros((1, nTime))
    go_pulse2[0, (MO2-50):(MO2-50+10)] = gaussian_p(timepoints, MO2-50, MO2-50+10,
                                                    sigma_for_go)
    go_pulse2[0, (MO2-50-10):(MO2-50)] = gaussian_p(timepoints, MO2-50, MO2-50-10,
                                                    sigma_for_go)
    
    ip1 = np.kron(h1.T, targ1_cos)
    ip2 = np.kron(h1.T, targ1_sin)
    ip3 = np.kron(h1.T, targ2_cos)*factor
    ip4 = np.kron(h1.T, targ2_sin)*factor
    ip5 = np.tile(go_pulse.T+go_pulse2.T, (1, ip1.shape[1]))

    ip3[:, cond_df['Type']=='SR'] = 0
    ip4[:, cond_df['Type']=='SR'] = 0

    x_time = np.stack([ip1, ip2, ip3, ip4, ip5]).T

    # numpy to torch.tensor
    x = torch.from_numpy(np.array(x_time, dtype=np.float32))

    return x, cond_df


def gen_2tinput_target_specific(targ1, targ2, input_form, **kws):
    # kws = {'factor':1, 'overlaps':(0, 0)}

    x_tensor, cond_df = gen_2tinput(targ1, targ2, input_form, factor=kws['factor'])
    y_tensor = gen_target(cond_df, kws['overlap1'], kws['overlap2'])

    return x_tensor, y_tensor, cond_df


#%% Define RNN training

def training_rnn_regularity(model, trainingset, nsteps_ntrain,
                            input_form, ini_state, save_name, **kws):

    # set training set
    train_targ1 = trainingset['target1']
    train_targ2 = trainingset['target2']
    nCond = len(train_targ1)
    x_tensor_o, y_tensor_o, train_cond_df_o = gen_input_target_specific(
                train_targ1, train_targ2, input_form, **kws)

    # set training parameters
    nsteps, ntrain = nsteps_ntrain
    learning_rate = 0.001
    optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad,
                                        model.parameters()),
                                 lr=learning_rate)
    msef = nn.MSELoss()
    loss_list = []

    # training
    with tqdm.tqdm(total=nsteps) as pbar:
        for step in range(nsteps):

            optimizer.zero_grad()

            # make ntrain trials have random conditions
            random_condid = np.random.choice((nCond), size=(ntrain))
            train_cond_df = train_cond_df_o.loc[random_condid].rename_axis('CondIdx').reset_index()

            x_tensor = x_tensor_o[random_condid, :, :].to(device)
            y_tensor = y_tensor_o[random_condid, :, :].to(device)

            h, r, pred = model(x_tensor, ini_state)

            r2 = torch.square(r).sum()/model.hidden_size/nTime/nCond
            # nTime is a global variable and defined in main()

            E_r = msef(y_tensor, pred) + r2 

            E_r.requires_grad_(True)
            E_r.backward()
            optimizer.step()

            loss_list.append(E_r.cpu().detach().numpy())
            pbar.update(1)

    torch.save(model.state_dict(), save_name+'.pth')

    plt.figure()
    plt.plot(loss_list)
    plt.title(save_name)
    plt.show()

    return {'nSteps': nsteps, 'nTrain':ntrain,
            'Trainingset': (train_cond_df_o, train_cond_df,
                            x_tensor.cpu().detach().numpy(),
                            y_tensor.cpu().detach().numpy()),
            'learningrate': learning_rate, 'losslist': loss_list}


def train():
    ini_state = torch.zeros((1, hidden_size)).to(device)
    # hidden_size is a global variable and defined in main()

    # sequence pool
    seq_pool = {'target1': [0, 60, 120, 180, 240, 300,
                            0, 60, 120, 180, 240, 300,
                            0, 60, 120, 180, 240, 300,
                            0, 60, 120, 180, 240, 300,
                            0, 60, 120, 180, 240, 300,
                            0, 60, 120, 180, 240, 300], # in degrees
                'target2': [120, 180, 240, 300, 0,   60,
                            240, 300, 0,   60,  120, 180,
                            180, 240, 300, 0,   60,  120,
                            60,  120, 180, 240, 300, 0,
                            300, 0,   60,  120, 180, 240,
                            -1, -1, -1, -1, -1, -1]} # in degrees

    # training parameters
    steps_params = (500, 100) #nsteps, ntrain (trial number for each step)

    # train 100 models
    for i in range(100):
        # initialize network
        rrnn = rRNN(input_size, hidden_size, output_size).to(device)
        
        ft = training_rnn_regularity(rrnn, seq_pool, steps_params,
                                      'relative', ini_state, 
                                      PATH_SAVE + 'dr_rnn_%d' % i,
                                      factor=factor21, 
                                      overlap1=overlap1m, 
                                      overlap2=overlap2m)
        np.save(PATH_SAVE + 'dr_rnn_%d' % i, ft)
        print(i)
        print("")
        

#%% Analysis functions

def cal_pd(meanfr, coshere, sinhere):

    # calculate neuronal preferred direction
    # fr = a*coshere + bsinhere + c
    # pd = arctan2(b, a)
    # depth = sqrt(a^2 + b^2)

    numtime = len(np.arange(100, End+1-100, 20))
    numneuron = meanfr.shape[2]

    pd_fit_np = np.zeros((numtime, numneuron))
    pd_baseline_np = np.zeros_like(pd_fit_np)
    pd_depth_np = np.zeros_like(pd_fit_np)
    pd_r2_np = np.zeros_like(pd_fit_np)

    with tqdm.tqdm(total=int(numtime*numneuron)) as pbar:
        for iloc, iTime in enumerate(np.arange(100, End+1-100, 20)):
            for iNeuron in range(numneuron):
                reg_pd = LinearRegression()
                reg_pd.fit(np.vstack([coshere, sinhere]).T,
                           np.mean(meanfr[:, iTime-100:iTime+100, iNeuron], axis=1))

                pd_coef = reg_pd.coef_
                pd_fit_np[iloc, iNeuron] = np.arctan2(pd_coef[1], pd_coef[0])
                pd_baseline_np[iloc, iNeuron] = reg_pd.intercept_
                pd_depth_np[iloc, iNeuron] = np.linalg.norm(pd_coef)
                pd_r2_np[iloc, iNeuron] = reg_pd.score(np.vstack([coshere, sinhere]).T,
                                                       np.mean(meanfr[:, iTime-100:iTime+100, iNeuron], axis=1))
                pbar.update(1)

    pd_dict = {'fit': pd_fit_np, 'baseline': pd_baseline_np,
               'depth': pd_depth_np, 'r2':pd_r2_np}

    return pd_dict 

def fitting_PD(r_test_np, test_cd, load_name, load_path):

    loadornot = 0
    for root, dirs, files in os.walk(load_path):
        for file_name in files:
            if file_name.split('.')[0] == load_name: loadornot=1

    if loadornot == 0:
        # cs_label = ['Dir1', 'Dir2', 'MovDir2']

        numN, numT = r_test_np.shape[2], len(np.arange(100, End+1-100, 20))

        ccw_index = test_cd[test_cd['Type']=='CCW'].index.tolist()
        cw_index = test_cd[test_cd['Type']=='CW'].index.tolist()
        all_index = ccw_index + cw_index
        index_dict = {'ccw': ccw_index, 'cw': cw_index, 'all': all_index}

        # fitted pd with different directions/coordinate systems
        pd_dict = {}
        for idxkey in index_dict.keys():

            cos_d1 = np.cos(
                test_cd.loc[index_dict[idxkey], 'Dir1'].values.astype(np.float64)/180*np.pi)
            sin_d1 = np.sin(
                test_cd.loc[index_dict[idxkey], 'Dir1'].values.astype(np.float64)/180*np.pi)

            cos_d2 = np.cos(
                test_cd.loc[index_dict[idxkey], 'Dir2'].values.astype(np.float64)/180*np.pi)
            sin_d2 = np.sin(
                test_cd.loc[index_dict[idxkey], 'Dir2'].values.astype(np.float64)/180*np.pi)

            cos_d21 = np.cos(
                test_cd.loc[index_dict[idxkey], 'MovDir2'].values.astype(np.float64)/180*np.pi)
            sin_d21 = np.sin(
                test_cd.loc[index_dict[idxkey], 'MovDir2'].values.astype(np.float64)/180*np.pi)

            pd_1 = cal_pd(r_test_np[index_dict[idxkey], :, :], cos_d1, sin_d1)
            pd_2 = cal_pd(r_test_np[index_dict[idxkey], :, :], cos_d2, sin_d2)
            pd_3 = cal_pd(r_test_np[index_dict[idxkey], :, :], cos_d21, sin_d21)
            pd_list = [pd_1, pd_2, pd_3]
            pd_dict[idxkey] = pd_list

        pd_fit_np = np.zeros((numT, numN))
        pd_baseline_np = np.zeros_like(pd_fit_np)
        pd_depth_np = np.zeros_like(pd_fit_np)
        pd_r2_np = np.zeros_like(pd_fit_np)

        for iloc in range(numT):
            for iNeuron in range(numN):
                # here idxkey = 'all'
                pd_candidate = np.array([pd_1['r2'][iloc, iNeuron],
                                         pd_2['r2'][iloc, iNeuron],
                                         pd_3['r2'][iloc, iNeuron]])
                pickid = np.argmax(pd_candidate)

                pd_fit_np[iloc, iNeuron] = pd_list[pickid]['fit'][iloc, iNeuron]
                pd_baseline_np[iloc, iNeuron] = pd_list[pickid]['baseline'][iloc, iNeuron]
                pd_depth_np[iloc, iNeuron] = pd_list[pickid]['depth'][iloc, iNeuron]
                pd_r2_np[iloc, iNeuron] = pd_list[pickid]['r2'][iloc, iNeuron]

        pdfit_contents = {'pd': pd_fit_np, 'baseline': pd_baseline_np,
                          'depth': pd_depth_np, 'r2': pd_r2_np}
        pdfit_comparison = pd_dict

        np.save(load_path + load_name + '.pd_fit_contents02.npy',
                pdfit_contents, allow_pickle=True)
        np.save(load_path + load_name + '.pd_fit_comparison02.npy',
                pdfit_comparison, allow_pickle=True)

        print('Already saved.')

    else:
        pdfit_contents = np.load(load_path + load_name + '.pd_fit_contents02.npy',
                                 allow_pickle=True).item()
        pdfit_comparison = np.load(load_path + load_name + '.pd_fit_comparison02.npy',
                                   allow_pickle=True).item()
        print('Already loaded')

    return pdfit_contents, pdfit_comparison

#@title @fitting_full_models

def fullfit(x, p):
    a, c, e, pd_, f = p
    return (a*np.cos(x[:,0]-pd_)
            # + b*np.cos(x[:, 1]-pd)
            + c*np.cos(x[:, 2]-pd_)
            # + d*np.cos(x[:,0]-pd)*np.cos(x[:, 1]-pd)
            + e*np.cos(x[:, 0]-pd_)*np.cos(x[:, 2]-pd_) + f)


def residuals(p, y, x): 
    return y - fullfit(x, p)


def cal_r2(yy1, yy2):
    return 1 - (np.sum((yy2-yy1)**2)/np.sum((yy1-np.mean(yy1))**2)).astype(np.float64)


def fitting_full_models(r_test_np, pdfit_contents, test_cd, load_name, load_path):
    loadornot = 0
    for root, dirs, files in os.walk(load_path):
        for file_name in files:
            if file_name.split('.')[0] == load_name: loadornot=1

    if loadornot == 0:
        # numN, numT = r_test_np.shape[2], len(np.arange(100, End+1-100, 20))
        numN = r_test_np.shape[2]

        FrMean = r_test_np.copy()[test_cd['TypeIdx']>0, :, :]

        full_collection = []

        print('\n' + load_name + '\n')
        with tqdm.tqdm(total=numN) as pbar:
            for iNeuron in range(numN):

                # for each neuron
                full_coef_list, full_bl_list, full_r2_list = [], [], []

                for iloc, iTime in enumerate(np.arange(100, End+1-100, 20)):

                    frtofit = np.mean(FrMean[:, iTime-100:iTime+100, iNeuron], axis=1)
                    plsq = leastsq(residuals,
                                    [pdfit_contents['depth'][iloc, iNeuron],
                                    0, 0,
                                    pdfit_contents['pd'][iloc, iNeuron],
                                    pdfit_contents['baseline'][iloc, iNeuron]],
                                    args=(frtofit,
                                          np.array([test_cd.loc[test_cd['TypeIdx']>0, 'Dir1'].values.astype(np.float64)/180*np.pi,
                                                    test_cd.loc[test_cd['TypeIdx']>0, 'Dir2'].values.astype(np.float64)/180*np.pi,
                                                    test_cd.loc[test_cd['TypeIdx']>0, 'MovDir2'].values.astype(np.float64)/180*np.pi]).T))

                    predfr = fullfit(np.array([test_cd.loc[test_cd['TypeIdx']>0, 'Dir1'].values.astype(np.float64)/180*np.pi,
                                              test_cd.loc[test_cd['TypeIdx']>0, 'Dir2'].values.astype(np.float64)/180*np.pi,
                                              test_cd.loc[test_cd['TypeIdx']>0,'MovDir2'].values.astype(np.float64)/180*np.pi]).T,
                                      plsq[0])

                    full_coef = plsq[0][0:-1]
                    full_baseline = plsq[0][-1]
                    full_r2 = cal_r2(frtofit, predfr)

                    full_coef_list.append(full_coef)
                    full_bl_list.append(full_baseline)
                    full_r2_list.append(full_r2)

                full_coef_np = np.array(full_coef_list) # time x coefficients

                full_dict = {'coef_list': full_coef_list,
                            'coef_np': full_coef_np,
                            'baseline': full_bl_list,
                            'r2': full_r2_list}

                full_collection.append(full_dict) # add each neuron
                pbar.update(1)

            full_dr_models = {'full_collection': full_collection}

        np.save(load_path + load_name + '.fit_models02.npy', full_dr_models, allow_pickle=True)
        print('Already saved.')

    else:

        full_dr_models = np.load(load_path + load_name + '.fit_models02.npy', allow_pickle=True).item()
        print('Already loaded.')

    return full_dr_models

#@title @getting_neural_states_withPCA

def getting_neural_states_withPCA(r_test_np, t_left, t_right):
    r_c_nt = np.zeros((r_test_np.shape[0], r_test_np.shape[2]))
    for ic in range(r_test_np.shape[0]):
        r_c_nt[ic, :] = np.mean(r_test_np[ic, t_left:t_right, :], axis=0) #time-averaged

    pca = PCA(n_components=3)
    pcs = pca.fit_transform(r_c_nt)
    evr = pca.explained_variance_ratio_

    return pcs, evr

#@title @getting_neuronal_PC

def getting_neuronal_PC(r_test_np):
    r_n_ct = np.zeros((r_test_np.shape[0]*r_test_np.shape[1], r_test_np.shape[2]))
    for iN in range(r_test_np.shape[2]):
        r_n_ct[:, iN] = r_test_np[:, :, iN].flatten()

    pca = PCA(n_components=3)
    pcs = pca.fit_transform(r_n_ct)
    pcs_np = np.zeros((r_test_np.shape[0], r_test_np.shape[1], 3))
    for iPC in range(3):
        pcs_np[:, :, iPC] = pcs[:, iPC].reshape(r_test_np.shape[0], r_test_np.shape[1])
    evr = pca.explained_variance_ratio_
    return pcs_np, evr


#%% Perturbation functions

def getting_inactivation(ic_idx, candi_kd_idx, model, nsteps, 
                         x_tensor, ini_state, y_test_np):

    msebs_list = []
    msebs_1 = []
    msebs_2 = []
    
    with tqdm.tqdm(total=nsteps) as pbar: 
        for bootsteps in range(nsteps):
    
            rnn_kd = rRNN(input_size, hidden_size, output_size)
            rnn_kd.load_state_dict(model.state_dict())
    
            kd_idx = random.sample(candi_kd_idx, 10)
            rnn_kd.state_dict()['cell.w_ih.weight'][kd_idx, :] = 0
            rnn_kd.state_dict()['cell.w_hh.weight'][kd_idx, :] = 0
            rnn_kd.state_dict()['cell.w_hh.weight'][:, kd_idx] = 0
            rnn_kd.state_dict()['cell.w_ho.weight'][:, kd_idx] = 0
    
            hkd, rkd, predkd = rnn_kd(x_tensor, ini_state)
            p_kd_np = predkd.detach().numpy()
    
            ## Behavioral accuracy
            mse_list = []
            temp1 = []
            temp2 = []
            for ic in ic_idx: 
                
                mse_c = np.sum((y_test_np[ic, MO:TT2, :] - p_kd_np[ic, MO:TT2, :])**2, axis=0)/y_test_np.shape[1]
                mse_list.append(mse_c)
                
                mse_c = np.sum((y_test_np[ic, MO:TT, :] - p_kd_np[ic, MO:TT, :])**2, axis=0)/y_test_np.shape[1]
                temp1.append(mse_c)
                
                mse_c = np.sum((y_test_np[ic, MO2:TT2, :] - p_kd_np[ic, MO2:TT2, :])**2, axis=0)/y_test_np.shape[1]
                temp2.append(mse_c)
   
            mse_np = np.array(mse_list)
            mse_mean = np.mean(mse_np)
            
            msebs_list.append(mse_mean)
            msebs_1.append(np.mean(np.array(temp1)))
            msebs_2.append(np.mean(np.array(temp2)))
            pbar.update(1)

    return (msebs_list, msebs_1, msebs_2)


def getting_highlow_indices(full_collection_dict, r_np, idxforcoef):
  mul_coef_in = []
  nozero_idx = []
  nNeuron = r_np.shape[2]
  for iNeuron in range(nNeuron):
      mul_coef_in.append(full_collection_dict[iNeuron]['coef_np'][:, idxforcoef])
      if r_np[:, :, iNeuron].max()>0:
          nozero_idx.append(iNeuron)
  mul_coef_np = abs(np.array(mul_coef_in))
  nozero_idx = np.array(nozero_idx)
  nz_mul_coef_np = mul_coef_np[nozero_idx, :]

  # time_period = np.arange(TO, TT2)
  # goindex = np.where((np.arange(100, End+1-100, 20)<GO)&(np.arange(100, End+1-100, 20)>=GO-600))
  high_idxfromzi = np.argwhere(
              np.mean(nz_mul_coef_np, axis=1)>
              np.percentile(np.mean(nz_mul_coef_np, axis=1), 90)).flatten().tolist()
  high_idx = nozero_idx[high_idxfromzi].tolist()

  low_idxfromzi = np.argwhere(
              np.mean(nz_mul_coef_np, axis=1)<
              np.percentile(np.mean(nz_mul_coef_np, axis=1), 10)).flatten().tolist()
  low_idx = nozero_idx[low_idxfromzi].tolist()

  return high_idx, low_idx, nozero_idx.tolist()


#%% 2 triggers
def training_2trnn_regularity(model, trainingset, nsteps_ntrain,
                            input_form, ini_state, save_name, **kws):

    # set training set
    train_targ1 = trainingset['target1']
    train_targ2 = trainingset['target2']
    nCond = len(train_targ1)
    x_tensor_o, y_tensor_o, train_cond_df_o = gen_2tinput_target_specific(
                train_targ1, train_targ2, input_form, **kws)

    # set training parameters
    nsteps, ntrain = nsteps_ntrain
    learning_rate = 0.001
    optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad,
                                        model.parameters()),
                                        lr=learning_rate)
    msef = nn.MSELoss()
    loss_list = []

    # training
    with tqdm.tqdm(total=nsteps) as pbar:
        for step in range(nsteps):

            optimizer.zero_grad()

            random_condid = np.random.choice((nCond), size=(ntrain))
            train_cond_df = train_cond_df_o.loc[random_condid].rename_axis('CondIdx').reset_index()

            x_tensor = x_tensor_o[random_condid, :, :].to(device)
            y_tensor = y_tensor_o[random_condid, :, :].to(device)

            h, r, pred = model(x_tensor, ini_state)

            r2 = torch.square(r).sum()/model.hidden_size/nTime/nCond
            # nTime is a global variable and defined elsewhere

            E_r = msef(y_tensor, pred) + r2

            E_r.requires_grad_(True)
            E_r.backward()
            optimizer.step()

            loss_list.append(E_r.cpu().detach().numpy())
            pbar.update(1)

    torch.save(model.state_dict(), save_name+'.pth')

    plt.figure()
    plt.plot(loss_list)
    plt.title(save_name)
    plt.show()

    return {'nSteps': nsteps, 'nTrain':ntrain,
            'Trainingset': (train_cond_df_o, train_cond_df,
                            x_tensor.cpu().detach().numpy(),
                            y_tensor.cpu().detach().numpy()),
            'learningrate': learning_rate, 'losslist': loss_list}


def train_2t():
    ini_state = torch.zeros((1, hidden_size)).to(device)

    # sequence pool
    seq_pool = {'target1': [0, 60, 120, 180, 240, 300,
                            0, 60, 120, 180, 240, 300,
                            0, 60, 120, 180, 240, 300,
                            0, 60, 120, 180, 240, 300,
                            0, 60, 120, 180, 240, 300,
                            0, 60, 120, 180, 240, 300], # in degrees
                'target2': [120, 180, 240, 300, 0,   60,
                            240, 300, 0,   60,  120, 180,
                            180, 240, 300, 0,   60,  120,
                            60,  120, 180, 240, 300, 0,
                            300, 0,   60,  120, 180, 240,
                            -1, -1, -1, -1, -1, -1]} # in degrees


    # training parameters
    steps_params = (500, 100)

    # train finely
    for i in range(1):
        # initialize network
        rrnn = rRNN(input_size, hidden_size, output_size).to(device)
        
        ft3 = training_2trnn_regularity(rrnn, seq_pool, steps_params,
                                      'relative', ini_state, 
                                      PATH_SAVE + 'dr_2trnn_%d' % i,
                                      factor=factor21, 
                                      overlap1=overlap1m, 
                                      overlap2=overlap2m)
        np.save(PATH_SAVE + 'dr_2trnn_%d' % i, ft3)
        print(i)
        print("")


#%% Plot functions

def plot_training_process(loss_list):
    plt.figure(figsize=(6, 4))
    plt.plot(loss_list)
    plt.xlabel('Step')
    plt.ylabel('Cost function')
    plt.title('Training process', pad=10)
    plt.subplots_adjust(top=0.8, left=0.15, right=0.9, bottom=0.2)
    plt.show()
    # fig = plt.gcf()
    # fig.savefig(figsavepath + 'training_process_%s.png' % model_name, dpi=300)


def plot_inputs(test_cd, x_test_np):
    plt.figure(figsize=(16,4))
    for isub in range(6):
        plt.subplot(3, 6, isub+1)
        plt.plot(x_test_np[isub, :, 0])
        plt.plot(x_test_np[isub, :, 1])
        plt.ylim(-1.5, 1.5)
        plt.axis('off')
        plt.title(test_cd.loc[isub, 'CondName'])

        plt.subplot(3, 6, isub+7)
        plt.plot(x_test_np[isub, :, 2])
        plt.plot(x_test_np[isub, :, 3])
        plt.ylim(-2, 2)
        plt.axis('off')

        plt.subplot(3, 6, isub+13)
        plt.plot(x_test_np[isub, :, 4])
        plt.ylim(-0.3, 1.7)
        plt.axis('off')
    plt.suptitle('Present input')
    plt.subplots_adjust(top=0.9, wspace=0.3)
    plt.show()


def plot_output_compared_with_target(test_cd, p_test_np, y_test_np):
    plot_idx = np.arange(p_test_np.shape[0]) #np.random.choice(x_test_np.shape[0], 18)
    plt.figure(figsize = (16, 3*int(p_test_np.shape[0]/6)))
    for isub, itrial in enumerate(plot_idx):
        plt.subplot(int(p_test_np.shape[0]/6), 6, isub+1)
        plt.plot(np.arange(p_test_np.shape[1]), y_test_np[itrial, :, :], ls='--', label='Target')
        plt.plot(np.arange(p_test_np.shape[1]), p_test_np[itrial, :, :], label='Model')

        plt.title(test_cd.loc[isub, 'CondName'])
        plt.ylim(1.2*p_test_np.min(), 1.2*p_test_np.max())
        plt.xlim(TO-5, TT2+5)
        plt.title(itrial)
    plt.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.1,
                        wspace=0.6, hspace=0.8)
    # fig3 = plt.gcf()
    # fig3.savefig(figsavepath + 'behavior_%s.png' % model_name)
    plt.show()


def plot_neuronal_PSTH(r_test_np, test_cd, condtoplot, figsavepath=None):

    nN = r_test_np.shape[2]
    for iN in range(nN):
        plt.figure(figsize=(5, 4))
        # iN = 0
        cm = plt.cm.get_cmap('copper', 
                             len(test_cd[test_cd['CondName']==condtoplot].index.tolist()))
        for ic in test_cd[test_cd['CondName']==condtoplot].index.tolist():
            plt.plot(np.arange(r_test_np.shape[1]), r_test_np[ic, :, iN], color=cm(ic))
        plt.axvline(x=MO, ls='--', color='k', lw=0.5)
        plt.axvline(x=MO2, ls='--', color='k', lw=0.5)
        plt.ylim(0, 1)
        plt.xticks([TO, GO, MO, TT, MO2, TT2], labels=[])
        plt.yticks([])
        plt.ylabel('Activity', labelpad=15)
        plt.title(condtoplot)

        if figsavepath is not None:
          plt.savefig(figsavepath + 'Node_{}.png'.format(iN))
        plt.show()
        plt.close()


def plot_special_conds(r_test_np):
    plt.figure(figsize=(5, 3))
    # iN = 0
    for iN in range(r_test_np.shape[2]):
      plt.figure(figsize=(5, 3))
      cm = ['#42407B', '#962E44', '#2A768E', '#F3985B']

      for icolor, ic in enumerate([0, 6, 1, 24]):
          plt.plot(timepoints, r_test_np[ic, :, iN].T,
                    color=cm[icolor], lw=3.5)
      plt.ylim(0, 1)
      plt.title('Node %03d' % (iN+1), pad=15, fontsize=20)
      plt.xticks([TO, 300, 600, GO, MO, MO2, TT2],
                  labels = ['TO', '0.3', '0.6', 'GO', 'MO', 'MO2', 'TT2'],
                  fontsize=12)
      plt.axvline(x=MO, ls='--', color='k', lw=1.5)
      plt.axvline(x=MO2, ls='--', color='k', lw=1.5)
      plt.yticks([])
      plt.xlabel('Time (s)', labelpad=10, fontsize=12)
      plt.ylabel('Firing rate', labelpad=10)
      ax = plt.gca()
      ax.spines['bottom'].set_linewidth(2)
      ax.spines['top'].set_linewidth(0)
      ax.spines['left'].set_linewidth(0)
      ax.spines['right'].set_linewidth(0)
      plt.tight_layout()
      plt.axis('off')
      plt.subplots_adjust(right=0.95, left=0.05, top=0.8, bottom=0.1)
      # plt.savefig(figsavepath + 'Node_{}.png'.format(iN), dpi=600)
      plt.show()


def plot_PD_comparison(pdfit_comparison, load_name):
    cs_label = ['Dir1', 'Dir2', 'MovDir2']
    for idxkey in pdfit_comparison.keys():
            # here idxkey = 'all'
        plt.figure(figsize=(12, 2))
        for idf in range(3):
            plt.plot(np.arange(100, End+1-100, 20),
                    np.mean(pdfit_comparison[idxkey][idf]['r2'], axis=1),
                    label=cs_label[idf])
        plt.axvline(x=0, color='k', ls='--')
        plt.xlabel('Time (ms)')
        plt.axhline(y=0.5, color='k', ls='--')
        plt.ylim(-0.05, 1.01)
        plt.yticks([0, 0.5, 1.0], labels=[0, 0.5, 1.0])
        plt.xticks([TO, 200, 400, 600, 800, GO, MO, MO2, TT2],
                    labels = ['TO', 200, 400, 600, 800, 'GO', 'MO', 'MO2', 'TT2'])
        plt.axvline(x=MO, ls='--', color='k', lw=0.5)
        plt.axvline(x=MO2, ls='--', color='k', lw=0.5)
        plt.ylabel('R-square')
        plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        plt.suptitle('R-square for different coordinate system in %s' % load_name)
        plt.subplots_adjust(wspace=0.1, hspace=0.5, top=0.8)
        # plt.savefig(PATH_SAVE + 'PDfit//' + 'PDfit_%s_%s.png' % (idxkey, load_name))
        plt.show()


def plot_fitting_meanR2(full_dr_models, load_name):

    full_collection_dict = full_dr_models['full_collection']
    numN = len(full_collection_dict)
    numT = len(full_collection_dict[0]['r2'])

    r2_np = np.zeros((numT, numN))
    for iNeuron in range(numN):
        r2_np[:, iNeuron] = full_collection_dict[iNeuron]['r2']

    r2_mean = np.nanmean(r2_np, axis=1)
    r2_sem = np.nanstd(r2_np, axis=1)/np.sqrt(numN)

    plt.figure(figsize=(8, 3))
    plt.plot(np.arange(100, End+1-100, 20), r2_mean, lw=3)
    plt.fill_between(np.arange(100, End+1-100, 20),
                    r2_mean - 2*r2_sem,
                    r2_mean + 2*r2_sem, alpha=0.3)
    # plt.axvline(x=0, color='k', ls='--')
    plt.xlabel('Time (s)', fontsize=16, labelpad=15)
    # plt.yticks(np.arange(0, 1.05, 0.2), labels=[])
    # plt.ylim(-0.05, 0.55)
    # plt.yticks(np.arange(0, 1.05, 0.2),
    #            labels=np.arange(0, 1.05, 0.2))
    # plt.xticks([TO, 300, 600, GO, MO, MO2, TT2],
    #             labels = ['TO', '0.3', '0.6', 'GO', 'MO', 'MO2', 'TT2'],
    #             fontsize=16)
    plt.yticks([0.7, 0.9], fontsize=16)
    plt.axvline(x=MO, ls='--', color='k', lw=1.5)
    plt.axvline(x=MO2, ls='--', color='k', lw=1.5)
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    # plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    plt.suptitle('%s Mean R-square' % load_name)
    plt.subplots_adjust(top=0.8, bottom=0.3, left=0.1, right=0.9)
    # plt.savefig(PATH_SAVE + 'fit//' + 'mean_r2_%s_%s.png' % (idxkey, load_name))
    plt.show()


def plot_fitting_coefficients(coef_mean, coef_std, nNeuron, load_name):

    # plot
    label_list = ['1st reach weight', 'Add. weight', 'Mul. weight']
    plt.figure(figsize=(4, 2), dpi=300)
    for idl, labelh in enumerate(label_list):
        plt.plot(np.arange(100, End+1-100, 20)-MO, coef_mean[:, idl],
                lw=3, label=labelh)
        # plt.fill_between(timepoints,
        #                  coef_mean[:, idl] - coef_std[:, idl],
        #                  coef_mean[:, idl] + coef_std[:, idl], alpha=0.1)
        plt.fill_between(np.arange(100, End+1-100, 20)-MO,
                          coef_mean[:, idl] - 2*coef_std[:, idl]/np.sqrt(nNeuron),
                          coef_mean[:, idl] + 2*coef_std[:, idl]/np.sqrt(nNeuron), alpha=0.3)
    plt.xlim(-800, 400)
    # plt.xticks([TO, 200, 400, 600, 800, GO, MO, MO2, TT2],
                # labels = ['TO', 200, 400, 600, 800, 'GO', 'MO', 'MO2', 'TT2'])
    # plt.axvline(x=MO, ls='--', color='k', lw=1.5)
    # plt.axvline(x=MO2, ls='--', color='k', lw=1.5)
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", fontsize=20)
    # plt.legend(loc='upper left', fontsize=20)
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(1.5)
    # ax.spines['top'].set_linewidth('1.5')
    ax.spines['left'].set_linewidth(1.5)
    # ax.spines['right'].set_linewidth('1.5')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.title('%s mean(abs(coef))' % (load_name), pad=15)
    plt.xticks(np.array([-400, 0, 400]), fontsize=16)
    plt.yticks([0, round(coef_mean.max()*1.1, 2)], fontsize=16)
    plt.xlabel('Time to MO', fontsize=16, labelpad=15)
    plt.subplots_adjust(top=0.8, bottom=0.1, left=0.1, right=0.9)
    plt.title(load_name)
    # plt.savefig(PATH_SAVE + 'fit//' + 'coeffs_%s_%s.png' % (idxkey, load_name))
    plt.show()


def plot_neural_states(pcs, evr, title):

    plt.figure(figsize=(10, 3), dpi=300)
    # cmap = plt.get_cmap('copper', 6)
    shapes = ['^', 's', 'o', 'd', 'x']
    plt.subplot(131)
    for iplot in range(pcs.shape[0]):
      plt.scatter(pcs[iplot, 0], pcs[iplot, 1],
                  facecolors=cmap(iplot%6), marker=shapes[iplot//6])
    plt.xlabel('PC1 (%.2f' % (evr[0]*100) +'%)')
    plt.ylabel('PC2 (%.2f' % (evr[1]*100) +'%)')
    plt.xlim(-pcs.max()*1.2, pcs.max()*1.2)
    plt.ylim(-pcs.max()*1.2, pcs.max()*1.2)
    ax = plt.gca()
    ax.set_aspect(1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xticks([])
    plt.yticks([])
    # plt.axis('off')

    plt.subplot(132)
    for iplot in range(pcs.shape[0]):
      plt.scatter(pcs[iplot, 0], pcs[iplot, 2],
                  facecolors=cmap(iplot%6), marker=shapes[iplot//6])
    plt.xlabel('PC1 (%.2f' % (evr[0]*100) +'%)')
    plt.ylabel('PC3 (%.2f' % (evr[2]*100) +'%)')
    plt.xlim(-pcs.max()*1.2, pcs.max()*1.2)
    plt.ylim(-pcs.max()*1.2, pcs.max()*1.2)
    ax = plt.gca()
    ax.set_aspect(1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xticks([])
    plt.yticks([])

    plt.subplot(133)
    for iplot in range(pcs.shape[0]):
      plt.scatter(pcs[iplot, 1], pcs[iplot, 2],
                  facecolors=cmap(iplot%6), marker=shapes[iplot//6])
    plt.xlabel('PC2 (%.2f' % (evr[1]*100) +'%)')
    plt.ylabel('PC3 (%.2f' % (evr[2]*100) +'%)')
    plt.xlim(-pcs.max()*1.2, pcs.max()*1.2)
    plt.ylim(-pcs.max()*1.2, pcs.max()*1.2)
    ax = plt.gca()
    ax.set_aspect(1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xticks([])
    plt.yticks([])
    # plt.legend()
    plt.suptitle(title)
    plt.subplots_adjust(wspace=0.4)
    plt.show()


def plot_neural_states_together(pcs1, pcs2, pcs3, evr1, evr2, evr3):

  plt.figure(figsize=(10, 3), dpi=300)
  # cmap = plt.get_cmap('copper', 6)
  shapes = ['^', 's', 'o', 'd', 'x']
  plt.subplot(131)
  for iplot in range(pcs1.shape[0]):
    plt.scatter(pcs1[iplot, 0], pcs1[iplot, 1],
                facecolors=cmap(iplot%6), marker=shapes[iplot//6])
  plt.xlabel('PC1 (%.1f' % (evr1[0]*100) +'%)')
  plt.ylabel('PC2 (%.1f' % (evr1[1]*100) +'%)')
  plt.xlim(-pcs1.max()*1.2, pcs1.max()*1.2)
  plt.ylim(-pcs1.max()*1.2, pcs1.max()*1.2)
  ax = plt.gca()
  ax.set_aspect(1)
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  plt.xticks([])
  plt.yticks([])
  # plt.axis('off')
  plt.title('[GO-600, GO]')

  plt.subplot(132)
  for iplot in range(pcs2.shape[0]):
    plt.scatter(pcs2[iplot, 0], pcs2[iplot, 1],
                facecolors=cmap(iplot%6), marker=shapes[iplot//6])
  plt.xlabel('PC1 (%.1f' % (evr2[0]*100) +'%)')
  plt.ylabel('PC2 (%.1f' % (evr2[1]*100) +'%)')
  plt.xlim(-pcs2.max()*1.2, pcs2.max()*1.2)
  plt.ylim(-pcs2.max()*1.2, pcs2.max()*1.2)
  ax = plt.gca()
  ax.set_aspect(1)
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  plt.xticks([])
  plt.yticks([])
  plt.title('[MO-100, MO+100')

  plt.subplot(133)
  for iplot in range(pcs3.shape[0]):
    plt.scatter(pcs3[iplot, 0], pcs3[iplot, 1],
                facecolors=cmap(iplot%6), marker=shapes[iplot//6])
  plt.xlabel('PC1 (%.1f' % (evr3[0]*100) +'%)')
  plt.ylabel('PC2 (%.1f' % (evr3[1]*100) +'%)')
  plt.xlim(-pcs3.max()*1.2, pcs3.max()*1.2)
  plt.ylim(-pcs3.max()*1.2, pcs3.max()*1.2)
  ax = plt.gca()
  ax.set_aspect(1)
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  plt.xticks([])
  plt.yticks([])
  plt.title('[MO2-100, MO2+100]')
  # plt.legend()
  plt.subplots_adjust(wspace=0.4)
  plt.show()


def plot_timeline():
    ## Visualize
    plt.figure(figsize=(16, 3))
    plt.plot(timepoints, np.zeros_like(timepoints))
    # plot time markers
    plt.scatter([TO, TOFF, GO, MO, TT, MO2, TT2], [-1, -1, -1, -1, -1, -1, -1])
    for v in ['TO', 'TOFF', 'GO', 'MO', 'TT', 'MO2', 'TT2']:
      plt.text(eval(v), -1.5, v+'\n'+'%d ms'%eval(v), ha='center', wrap=True)
    plt.ylim(-1.5, 1.5)
    plt.axis('off')
    plt.show()


def plot_example_targets():
    test_targ1 = [0, 60, 120, 180, 240, 300,
                  0, 60, 120, 180, 240, 300,
                  0, 60, 120, 180, 240, 300,
                  0, 60, 120, 180, 240, 300,
                  0, 60, 120, 180, 240, 300,
                  0, 60, 120, 180, 240, 300] # in deg
    test_targ2 = [120, 180, 240, 300, 0,   60,
                  240, 300, 0,   60,  120, 180,
                  180, 240, 300, 0,   60,  120,
                  60,  120, 180, 240, 300, 0,
                  300, 0,   60,  120, 180, 240,
                  -1, -1, -1, -1, -1, -1] # in deg
    
    x_tensor, _ = gen_input(test_targ1, test_targ2, 'actual')
    # print(cond_df)
    
    x_tensor_r, yp_tensor, cd = gen_input_target_specific(test_targ1, test_targ2,
                                                          'relative',
                                                          factor=1,
                                                          overlap1=10, overlap2=10)
    random_condid = np.arange(36) #np.random.choice((36), size=(36))
    test_cd = cd.loc[random_condid].rename_axis('CondIdx').reset_index()
    print(test_cd)
    x_tensor_r = x_tensor_r[random_condid, :, :]
    yp_tensor = yp_tensor[random_condid, :, :]
    
    icond = 7
    # Distribution of targets (inputs: actual/relative)
    plt.figure(figsize=(1.5, 1.5), dpi=300)
    plt.plot(np.cos(np.linspace(0, 2*np.pi, 100)),
              np.sin(np.linspace(0, 2*np.pi, 100)), '--')
    plt.plot(x_tensor.data.numpy()[icond, TO, 0],
              x_tensor.data.numpy()[icond, TO, 1], color='g', marker='s')
    
    plt.plot(x_tensor.data.numpy()[icond, TO, 2],
              x_tensor.data.numpy()[icond, TO, 3], color='g', marker='^')
    
    plt.plot(x_tensor_r.data.numpy()[icond, TO, 2],
              x_tensor_r.data.numpy()[icond, TO, 3], color='b', marker='^')
    plt.xlim(-2.5, 2.5)
    plt.ylim(-2.5, 2.5)
    plt.plot(0, 0, color='k', marker='.')
    ax = plt.gca()
    ax.set_aspect(1)
    plt.axis('off')
    #plt.axis('equal')
    plt.title('%s %d' % (test_cd.loc[icond, 'Type'], test_cd.loc[icond, 'Dir1']))
    plt.show()


    ## time line
    plt.figure(figsize=(1, 1.5), dpi=300)
    
    plt.subplot(4, 1, 1)
    plt.plot(x_tensor_r.data.numpy()[icond, :, 0])
    plt.plot(x_tensor_r.data.numpy()[icond, :, 1])
    plt.ylim(-1.5, 1.5)
    plt.axis('off')
    
    plt.subplot(4, 1, 2)
    plt.plot(x_tensor_r.data.numpy()[icond, :, 2])
    plt.plot(x_tensor_r.data.numpy()[icond, :, 3])
    plt.ylim(-2, 2)
    plt.axis('off')
    
    plt.subplot(4, 1, 3)
    plt.plot(x_tensor_r.data.numpy()[icond, :, 4])
    plt.ylim(-0.3, 1.7)
    plt.axis('off')
    
    plt.subplot(4, 1, 4)
    plt.plot(yp_tensor.data.numpy()[icond, :, 0])
    plt.plot(yp_tensor.data.numpy()[icond, :, 1])
    plt.ylim(-12, 12)
    plt.axis('off')
    
    # plt.suptitle('%s %d' % (test_cd.loc[icond, 'Type'], test_cd.loc[icond, 'Dir1']))
    plt.subplots_adjust(top=0.9, wspace=0.3)
    plt.show()

def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "ns"


def plot_pca_distance(pca_distance_result):
    epochlist = pca_distance_result.keys()
    for idp, ide in enumerate(epochlist):
        
        r2_dis_pd = pca_distance_result[ide]['in1Group']
        r1_dis_pd = pca_distance_result[ide]['in21Group']
        
        plt.figure(figsize=(6, 3), dpi=300)
        ax1 = plt.subplot(1, 2, 1)
        bplot1 = plt.boxplot(r2_dis_pd, patch_artist=True, labels=r2_dis_pd.columns,
                             medianprops=dict(linewidth=1.5))
        plt.ylim(-0.1, 1.6)
        plt.yticks([0.0, 0.5, 1.0])
        plt.title('For 1st-reach group %s' % ide)
        for patch in bplot1['boxes']:
            patch.set_facecolor('black')
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        
        ax2 = plt.subplot(1, 2, 2)
        bplot2 = plt.boxplot(r1_dis_pd, patch_artist=True, labels=r1_dis_pd.columns,
                             medianprops=dict(linewidth=1.5))
        plt.ylim(-0.1, 1.6)
        plt.yticks([0.0, 0.5, 1.0])
        plt.title('For DR-condition group %s' % ide)
        for patch in bplot2['boxes']:
            patch.set_facecolor('black')
            
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        plt.show()
        
        
def plot_kd_performance(kd_df):
    for rtype in ['MSE-total', 'MSE-1', 'MSE-2']:
        plt.figure(figsize=(4, 5), dpi=300)
        sns.boxplot(data=kd_df, x='kdidx', y=rtype, hue='cftype')
        plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
        plt.title(rtype)
        if rtype == 'MSE-total':            
            plt.ylim(0, 1.8)
        else:
            plt.ylim(0, 1)
        # plt.yticks([0.6, 0.7, 0.8, 0.9, 1.0])
    
        s, p = stats.mannwhitneyu(kd_df.query("kdidx=='high_idx' and cftype=='mul'")[rtype].values,
                                  kd_df.query("kdidx=='high_idx' and cftype=='add'")[rtype].values,
                                  alternative='greater')
        print('p(high mul > high add) = %.4f' %p)
    
        # plt.plot([0, 0, 0.27, 0.27], [1, 1.02, 1.02, 1], lw=1, c="k")
        # # plt.plot([0, 0, 0.27, 0.27], [0.7, 0.75, 0.75, 0.7], lw=1, c="k")
        # plt.text((0.27)*.5, 1.02, convert_pvalue_to_asterisks(p),
        #         ha='center', va='bottom', color="k")
    
        s, p = stats.mannwhitneyu(kd_df.query("kdidx=='high_idx' and cftype=='mul'")[rtype].values,
                                  kd_df.query("kdidx=='high_idx' and cftype=='1st'")[rtype].values,
                                  alternative='greater')
        print('p(high mul > high 1st) = %.4f' %p)
    
        # plt.plot([-0.27, -0.27, 0.27, 0.27], [1.05, 1.07, 1.07, 1.05], lw=1, c="k")
        # # plt.plot([-0.27, -0.27, 0.27, 0.27], [0.8, 0.85, 0.85, 0.8], lw=1, c="k")
        # plt.text(0, 1.07, convert_pvalue_to_asterisks(p),
        #         ha='center', va='bottom', color="k")
    
        s, p = stats.mannwhitneyu(kd_df.query("kdidx=='high_idx' and cftype=='mul'")[rtype].values,
                                  kd_df.query("kdidx=='low_idx' and cftype=='mul'")[rtype].values,
                                  alternative='greater')
        print('p(high mul > low mul) = %.4f' %p)
    
        # plt.plot([0.27, 0.27, 1.27, 1.27], [0.9, 0.95, 0.95, 0.9], lw=1, c="k")
        # plt.plot([0.27, 0.27, 1.27, 1.27], [1.1, 1.12, 1.12, 1.1], lw=1, c="k")
        # plt.text((0.27+1.27)/2, 1.12, convert_pvalue_to_asterisks(p),
        #         ha='center', va='bottom', color="k")
    
        s, p = stats.mannwhitneyu(kd_df.query("kdidx=='high_idx' and cftype=='mul'")[rtype].values,
                                  kd_df.query("kdidx=='random_idx' and cftype=='mul'")[rtype].values,
                                  alternative='greater')
        print('p(high mul > random mul) = %.4f' %p)
    
        # plt.plot([0.27, 0.27, 2.27, 2.27], [1.13, 1.15, 1.15, 1.13], lw=1, c="k")
        # # plt.plot([0.27, 0.27, 2.27, 2.27], [1.0, 1.05, 1.05, 1.0], lw=1, c="k")
        # plt.text((0.27+2.27)/2, 1.15, convert_pvalue_to_asterisks(p),
        #         ha='center', va='bottom', color="k")
        
        s, p = stats.mannwhitneyu(kd_df.query("kdidx=='low_idx' and cftype=='mul'")[rtype].values,
                                  kd_df.query("kdidx=='low_idx' and cftype=='1st'")[rtype].values,
                                  alternative='greater')
        print('p(low mul > low 1st) = %.4f' %p)
    
        # plt.plot([0.73, 0.73, 1.27, 1.27], [0.83, 0.85, 0.85, 0.83], lw=1, c="k")
        # # plt.plot([0.27, 0.27, 2.27, 2.27], [1.0, 1.05, 1.05, 1.0], lw=1, c="k")
        # plt.text((0.73+1.27)/2, 0.85, convert_pvalue_to_asterisks(p),
        #         ha='center', va='bottom', color="k")
        
        
        s, p = stats.mannwhitneyu(kd_df.query("kdidx=='low_idx' and cftype=='mul'")[rtype].values,
                                  kd_df.query("kdidx=='low_idx' and cftype=='add'")[rtype].values,
                                  alternative='greater')
        print('p(low mul > low add) = %.4f' %p)
    
        # plt.plot([1.0, 1.0, 1.27, 1.27], [0.73, 0.75, 0.75, 0.73], lw=1, c="k")
        # # plt.plot([0.27, 0.27, 2.27, 2.27], [1.0, 1.05, 1.05, 1.0], lw=1, c="k")
        # plt.text((1.0+1.27)/2, 0.75, convert_pvalue_to_asterisks(p),
        #         ha='center', va='bottom', color="k")
        
        
        ax = plt.gca()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    
        plt.show()
        
        
#%% Analyse one model

def analyse1(load_name, x_tensor, y_tensor, ini_state, test_cd):
    
    this_model = {}
    
    ## load model
    load_path = PATH_SAVE + load_name + '.pth'
    # load_model_name = load_path.split('.')[0]
    model = rRNN(input_size, hidden_size, output_size)
    model.load_state_dict(torch.load(load_path, map_location='cpu'))
    print(model)
        
    ## get model data
    nNeuron = model.hidden_size
    h, r, pred = model(x_tensor, ini_state)
    x_test_np = x_tensor.detach().numpy()
    y_test_np = y_tensor.detach().numpy()

    h_test_np = h.detach().numpy()
    r_test_np = r.detach().numpy()
    p_test_np = pred.detach().numpy()
    
    this_model['data'] = {'h':h_test_np, 'r':r_test_np, 'p':p_test_np}
    
    ## model performance
    mse_list = []
    r2_list = []
    for ic in range(x_test_np.shape[0]):
        mse_c = np.sum((y_test_np[ic, :, :] - p_test_np[ic, :, :])**2, axis=0)/x_test_np.shape[1]
        mse_list.append(mse_c)

        r2 = 1 - (np.sum((p_test_np[ic, :, :] - y_test_np[ic, :, :])**2)/
              np.sum((y_test_np[ic, :, :]
                              -np.tile(np.mean(y_test_np[ic, :, :], axis=0)[np.newaxis, :],
                                       (y_test_np.shape[1], 1)))**2))
        r2_list.append(r2)

    mse_np = np.array(mse_list)
    mse_mean = np.mean(mse_np)
    mse_std = np.std(mse_np)
    # print('%s, MSE = %.4f ± %.4f' % (load_name, mse_mean, mse_std))

    r2_np = np.array(r2_list)
    r2_mean = np.mean(r2_np)
    r2_std = np.std(r2_np)
    # print('%s, r2 = %.4f ± %.4f' % (load_name, r2_mean, r2_std))

    this_model['performance'] = {'mse_mean': mse_mean, 'mse_std': mse_std,
                                 'r2_mean': r2_mean, 'r2_std': r2_std}

    ## fitting PD
    load_path = PATH_SAVE + 'PDfit_' + load_name + '//'
    if operator.eq(os.path.exists(load_path), False):
         os.mkdir(load_path)

    pdfit_contents, pdfit_comparison = fitting_PD(r_test_np, test_cd, 
                                                  load_name, load_path)

    this_model['fitting_PD'] = {'pdfit': pdfit_contents,
                                'pdfit_comparison': pdfit_comparison}

    ## fitting full model
    load_path = PATH_SAVE + 'fit_' + load_name + '//'
    if operator.eq(os.path.exists(load_path), False):
        os.mkdir(load_path)

    full_dr_models = fitting_full_models(r_test_np, pdfit_contents, test_cd, load_name, load_path)
    
    this_model['fitting_Full'] = full_dr_models
    
    ## calculate coefficients mean and std, also similarity
    full_collection_dict = full_dr_models['full_collection']
    coef_in = []
    for iNeuron in range(nNeuron):
        coef_in.append(full_collection_dict[iNeuron]['coef_np'][:, 0:-1])
    coef_mean = np.mean(abs(np.array(coef_in)), axis=0)
    coef_std = np.std(abs(np.array(coef_in)), axis=0)

    fdis = sim.cal_similarity(coef_mean)
    sim.print_all(coef_mean)
    
    this_model['coef'] = {'mean': coef_mean, 'std': coef_std,
                          'fdis': fdis}
    
    ## distance between neural states
    pcs1, evr1 = getting_neural_states_withPCA(
        r_test_np[test_cd['TypeIdx']>0, :, :], GO-600, GO)

    pcs2, evr2 = getting_neural_states_withPCA(
        r_test_np[test_cd['TypeIdx']>0, :, :], MO-100, MO+100)
    
    pcs3, evr3 = getting_neural_states_withPCA(
        r_test_np[test_cd['TypeIdx']>0, :, :], MO2-100, MO2+100)
    
    this_model['pca_results'] = {'pcs': [pcs1, pcs2, pcs3],
                                 'evr': [evr1, evr2, evr3]}

    # group by 1st reach direction
    r1_dir, r1_indices = np.unique(test_cd[test_cd['TypeIdx']>0]['Dir1'],
                                   return_inverse=True)
    #[0, 60, 120, 180, 240, 300]
    r2_dir, r2_indices = np.unique(test_cd.loc[test_cd['TypeIdx']>0, 'CondName'],
                                   return_inverse=True)
    #['180', 'CCW120', 'CCW60', 'CW120', 'CW60']

    this_model['pca_distance'] = {}
    epochlist = ['GO', 'MO', 'MO2']
    for idp, pcs in enumerate([pcs1, pcs2, pcs3]):
        state_points = []
        for icond in range(len(test_cd[test_cd['TypeIdx']>0]['CondName'])):
          state_points.append((pcs[icond, 0], pcs[icond, 1]))
        dist_np = squareform(pdist(state_points, 'euclidean'))
        dist_np = dist_np/dist_np.max() #normalize

        # Grouped by 1st reach, within group
        r2_dis_dict = {}
        for dr2index, dr2name in enumerate(list(r2_dir)):
          templist = []
          for dr1index in range(len(r1_dir)):
            templist.append(dist_np[(r1_indices==dr1index)&(r2_indices==0),
                                    (r1_indices==dr1index)&(r2_indices==dr2index)][0])
          r2_dis_dict[dr2name] = templist
        
        r2_dis_pd = pd.DataFrame(r2_dis_dict)
        d2nameorder = ['CCW60', 'CCW120', '180', 'CW120', 'CW60']
        r2_dis_pd = r2_dis_pd[d2nameorder]
                
        # Grouped by 2nd reach, within group
        r1_dis_dict = {}
        for dr1index, dr1name in enumerate(r1_dir):
          templist = []
          for dr2index in range(len(r2_dir)):
              templist.append(dist_np[(r1_indices==0)&(r2_indices==dr2index),
                                      (r1_indices==dr1index)&(r2_indices==dr2index)][0])
          r1_dis_dict[str(dr1name)] = templist
        
        r1_dis_pd = pd.DataFrame(r1_dis_dict)
    
        this_model['pca_distance'][epochlist[idp]] = {'in1Group': r2_dis_pd, 
                                                      'in21Group': r1_dis_pd}   

    
    # perturbation as knocking down certain nodes
    full_collection_dict = full_dr_models['full_collection']

    dr_idx = test_cd[test_cd['TypeIdx']>0].index.tolist()
    # sr_idx = test_cd[test_cd['TypeIdx']==0].index.tolist()

    kdtype = ['high_idx', 'low_idx', 'random_idx']
    idxtype = ['1st', 'add', 'mul']
    kd_df = pd.DataFrame()
    for iii, cfidx in enumerate(idxtype):
      hi, li, nzi = getting_highlow_indices(full_collection_dict, r_test_np, iii)
      
      for ii, kd_idx_list in enumerate([hi, li, nzi]):
        mseA, mse1, mse2 = getting_inactivation(dr_idx, kd_idx_list, model, 100,
                                                x_tensor, ini_state, y_test_np)
        kdi = pd.DataFrame({'MSE-total': mseA, 'MSE-1': mse1, 'MSE-2': mse2,
                            'kdidx':[kdtype[ii] for i in range(len(mseA))],
                            'cftype': [cfidx for i in range(len(mseA))]})
        kd_df = kd_df.append(kdi, ignore_index=True)
    
    this_model['knock_down_results'] = kd_df
        
    return model, this_model


#%% Plot one model    

def plot1(load_name, model_to_plot, test_cd, y_tensor, plot_neurons=0):
    
    ## See training process
    if os.path.exists(PATH_SAVE + load_name + '.npy'):
        seetraining = np.load(PATH_SAVE + load_name + '.npy', allow_pickle=True).item()
        loss_list = seetraining['losslist']
        plot_training_process(loss_list)
    
    ## Compare output with target
    y_test_np = y_tensor.detach().numpy()
    p_test_np = model_to_plot['data']['p']
    plot_output_compared_with_target(test_cd, p_test_np, y_test_np)
    
    ## Plot PSTH for each node
    r_test_np = model_to_plot['data']['r']
    ## plot_neuronal_PSTH(r_test_np, test_cd, figsavepath=None)
    if plot_neurons==1:
        plot_special_conds(r_test_np)
    
    ## Compare R-square of PD fitting with different directions
    plot_PD_comparison(model_to_plot['fitting_PD']['pdfit_comparison'], load_name)

    ## Mean R-square
    plot_fitting_meanR2(model_to_plot['fitting_Full'], load_name)

    ## Plot time-varying coefficients in full model
    full_collection_dict = model_to_plot['fitting_Full']['full_collection']
    nNeuron = len(full_collection_dict)
    coef_in = []
    for iNeuron in range(nNeuron):
        coef_in.append(full_collection_dict[iNeuron]['coef_np'][:, 0:-1])
    coef_mean = np.mean(abs(np.array(coef_in)), axis=0)
    coef_std = np.std(abs(np.array(coef_in)), axis=0)
    
    plot_fitting_coefficients(coef_mean, coef_std, nNeuron, load_name)
    # plot_fitting_coefficients(model_to_plot['fitting_Full'], load_name)
    
    ## Plot PCA results
    pcs1, pcs2, pcs3 = model_to_plot['pca_results']['pcs']
    evr1, evr2, evr3 = model_to_plot['pca_results']['evr']
    plot_neural_states_together(pcs1, pcs2, pcs3, evr1, evr2, evr3)
        
    ## Plot PCA states' distance
    plot_pca_distance(model_to_plot['pca_distance'])
        
        ## Plot Knock-down result
    plot_kd_performance(model_to_plot['knock_down_results'])

    
#%% 
def main():
    
    global TO, TOFF, GO, MO, TT, MO2, TT2, timepoints, End, nTime, R_circle
    TO, TOFF, GO, MO, TT, MO2, TT2 = 0, 400, 1000, 1200, 1350, 1450, 1600
    timepoints = np.arange(TT2+100+1)
    End = timepoints[-1]
    nTime = len(timepoints)
    R_circle = 0.15 #m
    
    plot_timeline()

    global input_size, hidden_size, output_size
    input_size = 5
    hidden_size = 200
    output_size = 2

    global input_form
    input_form = 'relative'
    
    global factor21, overlap1m, overlap2m
    factor21, overlap1m, overlap2m = 1, 0, 0
 
    plot_example_targets()   
    
    global cmap
    cmap = tuple(map(tuple, np.array([[255, 27, 27, 255], 
                                              [116, 186, 94, 255],
                                              [17, 137, 207, 255], 
                                              [255, 198, 64, 255], 
                                              [25, 210, 210, 255], 
                                              [212, 33, 212, 255]])/255))
    cmap = ListedColormap(cmap, name='color_md')
 
    ## train models
    # train()

    ## generate validation data
    test_targ1 = [0, 60, 120, 180, 240, 300,
                  0, 60, 120, 180, 240, 300,
                  0, 60, 120, 180, 240, 300,
                  0, 60, 120, 180, 240, 300,
                  0, 60, 120, 180, 240, 300,
                  0, 60, 120, 180, 240, 300] # in deg
    test_targ2 = [120, 180, 240, 300, 0,   60,
                  240, 300, 0,   60,  120, 180,
                  180, 240, 300, 0,   60,  120,
                  60,  120, 180, 240, 300, 0,
                  300, 0,   60,  120, 180, 240,
                  -1, -1, -1, -1, -1, -1] # in deg
    
  
    nTest = len(test_targ1) 
    x_tensor, y_tensor, cd = gen_input_target_specific(test_targ1, test_targ2,
                                                       input_form,
                                                       factor=1,
                                                       overlap1=0, overlap2=0)
    random_condid = np.arange(nTest)
    test_cd = cd.loc[random_condid].rename_axis('CondIdx').reset_index()
    x_tensor = x_tensor[random_condid, :, :]
    y_tensor = y_tensor[random_condid, :, :]
    print(test_cd)
    ini_state = torch.zeros((1, hidden_size)).to('cpu')
    
    plot_inputs(test_cd, x_tensor.detach().numpy())
    
    ## analyse one model
    # load trained models
    load_name = 'dr_rnn_1001'

    # analyse
    model, this_model = analyse1(load_name, x_tensor, y_tensor, ini_state, test_cd)
    
    # plot
    plot1(load_name, this_model, test_cd, y_tensor)
    
    
    ## analyse a batch of models
    # batch coding and save
    for i10 in range(10):
        these_models = []
        for i in range(10):
            load_name = 'dr_rnn_%d' % (i10*10+i)
            _, this_model = analyse1(load_name, x_tensor, y_tensor, ini_state, test_cd)
            these_models.append(this_model)
        np.save(PATH_SAVE + 'dr_rnn_analysis_%d.npy' % i10, these_models)
        print(i10*10+i)
        del these_models
        
    # load
    these_models = []
    for i10 in range(10):
        temp = np.load(PATH_SAVE + 'dr_rnn_analysis_%d.npy' % i10, 
                                         allow_pickle=True)
        for i in range(10):
            these_models.append({k:v for k, v in temp[i].items() if k in [
                'performance', 'coef', 'pca_distance', 'knock_down_results']})
    
    ## find best model
    perform_r2_list = []
    perform_mse_list = []
    fdis_list = []
    for imodel in these_models:
        perform_r2_list.append(imodel['performance']['r2_mean'])
        perform_mse_list.append(imodel['performance']['mse_mean'])
        fdis_list.append(imodel['coef']['fdis'])
        
    best_performance = np.argsort(-np.array(perform_r2_list))
        
    best_similarity = np.argsort(np.array(fdis_list))
    
    meanrank = np.mean(np.array([np.argsort(best_performance), 
                                 np.argsort(best_similarity)]), axis=0)
    best_idx = np.argmin(meanrank)
    
    
    ## for best model
    load_name = 'dr_rnn_%d' % best_idx
    
    # analyse
    # model, this_model = analyse1(load_name, x_tensor, y_tensor, ini_state, test_cd)
    temp = np.load(PATH_SAVE + 'dr_rnn_analysis_%d.npy' % (best_idx//10), 
                                     allow_pickle=True)
    best_model = temp[best_idx%10]
    assert best_model['coef']['fdis'] == fdis_list[best_idx]
    
    # plot
    plot1(load_name, best_model, test_cd, y_tensor, plot_neurons=0)
    print(perform_r2_list[best_idx])
    print(fdis_list[best_idx])
    plot1(load_name, best_model, test_cd, y_tensor, plot_neurons=1)
   
    
    ## for all models
    # performance
    all_performance_r2 = [min(perform_r2_list), np.mean(perform_r2_list), 
                       np.std(perform_r2_list), max(perform_r2_list)]
    print(all_performance_r2)
    
    all_performance_mse = [min(perform_mse_list), np.mean(perform_mse_list), 
                       np.std(perform_mse_list), max(perform_mse_list)]
    print(all_performance_mse)
    
    # similarity
    all_fdis = [min(fdis_list), np.mean(fdis_list), max(fdis_list), np.std(fdis_list)]
    print(all_fdis)
    
    # coefficients
    mm_coef = np.array([imodel['coef']['mean'] for imodel in these_models])
    mmean_coef = np.mean(mm_coef, axis=0)
    mstd_coef = np.std(mm_coef, axis=0)
    
    plot_fitting_coefficients(mmean_coef, mstd_coef, 1, 'all')
    
    # PCA distance
    all_pca_distance = these_models[0]['pca_distance'].copy()
    for idx in np.arange(2, len(these_models)):
        for epoch in these_models[0]['pca_distance'].keys():
            for group in these_models[0]['pca_distance']['GO'].keys():
                all_pca_distance[epoch][group] = all_pca_distance[epoch][group].append(
                    these_models[idx]['pca_distance'][epoch][group],
                    ignore_index=True)
    plot_pca_distance(all_pca_distance)        
    
    # Knock-down result
    all_kd = these_models[0]['knock_down_results'].copy()
    for idx in np.arange(1, len(these_models)):
        all_kd = all_kd.append(
                    these_models[idx]['knock_down_results'], ignore_index=True)
    plot_kd_performance(all_kd)
    
    
    ## try 2 triggers
    x_tensor, y_tensor, cd = gen_2tinput_target_specific(test_targ1, test_targ2,
                                                         input_form,
                                                         factor=1,
                                                         overlap1=0, overlap2=0)
    random_condid = np.arange(nTest)
    test_cd = cd.loc[random_condid].rename_axis('CondIdx').reset_index()
    x_tensor = x_tensor[random_condid, :, :]
    y_tensor = y_tensor[random_condid, :, :]
    print(test_cd)
    ini_state = torch.zeros((1, hidden_size)).to('cpu')
    
    plot_inputs(test_cd, x_tensor.detach().numpy())
    
    # train
    # train_2t()
    
    ## analyse one model
    # load trained models
    load_name = 'dr_2trnn_0'

    # analyse
    model, this_model = analyse1(load_name, x_tensor, y_tensor, ini_state, test_cd)
    np.save(PATH_SAVE + 'dr_2trnn_analysis_0.npy', this_model)
    
    # plot
    plot1(load_name, this_model, test_cd, y_tensor)
    

if __name__ == main():
    
    main()