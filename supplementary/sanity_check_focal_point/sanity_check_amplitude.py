# -*- coding: utf-8 -*-
"""Sanity_Check_Amplitude

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1SndzlQwmUGorM8qZGBv7NunVoknWbvDa
"""

import jax.numpy as jnp
import jax
from jax import experimental
from jax.experimental import optimizers
from jax import grad, jit, vmap
from jax import random
from scipy.special import jv
import os
import datetime
import time
import sys

import numpy as np
import pandas as pd
import time
import math
from jax import device_put

def prop_matrix(array_x, array_y, bottom_array_z, top_array_z, geo_arr, P0, k, r0, y):
  bottom_prop = []
  top_prop = []
  for z in range(point_num):
    bottom_dist_map = jnp.sqrt( jnp.power((array_x - geo_arr[y][z*3]), 2) + jnp.power((array_y - geo_arr[y][z*3+1]), 2) + jnp.power((bottom_array_z - geo_arr[y][z*3+2]), 2) )
    top_dist_map = jnp.sqrt( jnp.power((array_x - geo_arr[y][z*3]), 2) + jnp.power((array_y - geo_arr[y][z*3+1]), 2) + jnp.power((top_array_z - geo_arr[y][z*3+2]), 2) )

    bottom_sin_alpha_map = jnp.sqrt( (jnp.power((array_x - geo_arr[y][z*3]), 2) + jnp.power((array_y - geo_arr[y][z*3+1]), 2)) ) / bottom_dist_map
    top_sin_alpha_map = jnp.sqrt( (jnp.power((array_x - geo_arr[y][z*3]), 2) + jnp.power((array_y - geo_arr[y][z*3+1]), 2)) ) / top_dist_map
    
    #bottom_amplitude_map = ( P0 / bottom_dist_map )
    #top_amplitude_map = ( P0 / top_dist_map)

    bottom_amplitude_map = ( 2 * jv(1, k*r0*bottom_sin_alpha_map) * P0 / (k*r0*bottom_sin_alpha_map*bottom_dist_map) )
    top_amplitude_map = ( 2 * jv(1, k*r0*top_sin_alpha_map) * P0 / (k*r0*top_sin_alpha_map*top_dist_map) )
    
    bottom_prop.append( jax.lax.complex( bottom_amplitude_map * jnp.cos(k*bottom_dist_map), bottom_amplitude_map * jnp.sin(k*bottom_dist_map) ) )
    top_prop.append( jax.lax.complex( top_amplitude_map * jnp.cos(k*top_dist_map), top_amplitude_map * jnp.sin(k*top_dist_map) ) )
  
  bottom_prop = device_put(jnp.array(bottom_prop))
  top_prop = device_put(jnp.array(top_prop))

  return bottom_prop, top_prop

@jit
def loss_func(all_phase, top_prop, bottom_prop, target_amp):
    bottom_phase = all_phase[0]
    top_phase = all_phase[1]
    bottom_phase_exp = jax.lax.complex(jnp.cos(bottom_phase), jnp.sin(bottom_phase))
    top_phase_exp = jax.lax.complex(jnp.cos(top_phase), jnp.sin(top_phase))
    loss_array = 0
    for i in range(point_num):
      bottom_point_Re = jnp.real(bottom_phase_exp) * jnp.real(bottom_prop[i]) - jnp.imag(bottom_phase_exp) * jnp.imag(bottom_prop[i])
      bottom_point_Im = jnp.real(bottom_phase_exp) * jnp.imag(bottom_prop[i]) + jnp.imag(bottom_phase_exp) * jnp.real(bottom_prop[i])

      top_point_Re = jnp.real(top_phase_exp) * jnp.real(top_prop[i]) - jnp.imag(top_phase_exp) * jnp.imag(top_prop[i])
      top_point_Im = jnp.real(top_phase_exp) * jnp.imag(top_prop[i]) + jnp.imag(top_phase_exp) * jnp.real(top_prop[i])

      point_Re = jnp.sum( bottom_point_Re + top_point_Re )
      point_Im = jnp.sum( bottom_point_Im + top_point_Im )

      est_amp = jnp.sqrt( point_Re ** 2 + point_Im ** 2 )

      err = est_amp - target_amp[i]

      loss_array += ( (err**2) )
    
    
    return loss_array

@jit
def step(step, opt_state, top_prop, bottom_prop, target_amp):
  value, grads = jax.value_and_grad(loss_func)(get_params(opt_state), top_prop, bottom_prop, target_amp)
  opt_state = opt_update(step, grads, opt_state)
  return value, opt_state

def slow_loss_func_for_verification(all_phase, top_prop, bottom_prop, target_amp):
    bottom_phase = all_phase[0]
    top_phase = all_phase[1]
    bottom_phase_exp = jax.lax.complex(jnp.cos(bottom_phase), jnp.sin(bottom_phase))
    top_phase_exp = jax.lax.complex(jnp.cos(top_phase), jnp.sin(top_phase))
    loss_array = 0
    for i in range(point_num):
      bottom_point_Re = jnp.real(bottom_phase_exp) * jnp.real(bottom_prop[i]) - jnp.imag(bottom_phase_exp) * jnp.imag(bottom_prop[i])
      bottom_point_Im = jnp.real(bottom_phase_exp) * jnp.imag(bottom_prop[i]) + jnp.imag(bottom_phase_exp) * jnp.real(bottom_prop[i])

      top_point_Re = jnp.real(top_phase_exp) * jnp.real(top_prop[i]) - jnp.imag(top_phase_exp) * jnp.imag(top_prop[i])
      top_point_Im = jnp.real(top_phase_exp) * jnp.imag(top_prop[i]) + jnp.imag(top_phase_exp) * jnp.real(top_prop[i])

      point_Re = jnp.sum( bottom_point_Re + top_point_Re )
      point_Im = jnp.sum( bottom_point_Im + top_point_Im )

      est_amp = jnp.sqrt( point_Re ** 2 + point_Im ** 2 )

      err = est_amp - target_amp[i]
      #jax.device_put(x)[idx]
      loss_array += ( (err**2) )
      
    #print(loss_array)
    #loss_avg = jnp.mean(loss_array)
    #print("loss_avg: ", loss_avg)

    return loss_array

amplitudes_perecentage = np.arange(0.8, 1.3, 0.1)
amplitudes_perecentage

N=16
size=[N,N]
pitch=0.0105
r0=0.005
P0=1.98
l_ambda = 346.0 / 40000.0
k = 2.0*jnp.pi/l_ambda

out_phase_path = './iter_phases/'
if not os.path.exists(out_phase_path):
    os.mkdir(out_phase_path)
out_loss_path = './iter_loss/'
if not os.path.exists(out_loss_path):
    os.mkdir(out_loss_path)

array_x = jnp.array([ (i % N - (N//2-0.5))*pitch for i in range(N * N)]).reshape(N, N)
array_y = array_x.T
bottom_array_z = jnp.full((N,N), 0.0)
top_array_z = jnp.full((N,N), 0.23550056)

amplitudes_perecentage = [0.8, 0.9, 1.1, 1.2]

geo_df = pd.read_csv('geometries_001.csv', header=None)
geo_df = geo_df.drop(geo_df.columns[[0]], axis=1)
geo_arr = geo_df.values

amp_df = pd.read_csv('amplitudes_001.csv', header=None)
amp_df = amp_df.drop(amp_df.columns[[0]], axis=1)
amp_arr = amp_df.values

for x in range(len(amplitudes_perecentage)):
  point_num = 1
  fn_num = '{:0=3}'.format(point_num)

  out_phase_p_path = out_phase_path + 'p' + str(point_num) + '/'
  if not os.path.exists(out_phase_p_path):
    os.mkdir(out_phase_p_path)
  out_loss_p_path = out_loss_path + 'p' + str(point_num) + '/'
  if not os.path.exists(out_loss_p_path):
    os.mkdir(out_loss_p_path)

  all_phase_to_csv = []
  loss_list = []
  calc_time_list = []
  for y in range(250):
    start_time = time.time()
    init_bottom_phase = np.random.uniform(0.0, 2.0*np.pi, size)
    init_top_phase = np.random.uniform(0.0, 2.0*np.pi, size)
    all_phase = ([init_bottom_phase, init_top_phase])
    bottom_prop, top_prop = prop_matrix(array_x, array_y, bottom_array_z, top_array_z, geo_arr, P0, k, r0, y)
    target_amp = amp_arr[y]*amplitudes_perecentage[x]

    opt_init, opt_update, get_params = jax.experimental.optimizers.adam(0.1, b1=0.9, b2=0.999, eps=1e-08)
    opt_state = opt_init(all_phase)
    
    iteration_phases = []

    # for ii in range(200):
    for ii in range(250):
      value, opt_state = step(ii, opt_state, top_prop, bottom_prop, target_amp)


    all_phase = get_params(opt_state)   
    calc_time = time.time() - start_time
    #print("--- %s seconds ---" % calc_time)
    #loss_p = loss_func(all_phase, top_prop, bottom_prop, target_amp)
    # loss_c = slow_loss_func_for_verification(all_phase, top_prop, bottom_prop, target_amp)
    #print(loss_p)
    all_phase_1d_arr = np.array( [all_phase[0].T.ravel(), all_phase[1].T.ravel()] ).ravel()
    all_phase_to_csv.append(all_phase_1d_arr)

    calc_time_list.append(calc_time)

    #print(value)
    print("--- Geo num: %d ---" % y)

  all_phase_to_csv = np.array(all_phase_to_csv)
  np.savetxt('phases_' + str(int(amplitudes_perecentage[x]*100)) + '.csv', all_phase_to_csv, delimiter=',')