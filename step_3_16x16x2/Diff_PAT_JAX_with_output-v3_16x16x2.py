## Manuscript: "Acoustic Hologram Optimisation Using Automatic Differentiation"
## Authors: Tatsuki Fushimi*, Kenta Yamamoto*, Yoichi Ochiai
## Corresponding email: tfushimi@slis.tsukuba.ac.jp

# Loading up library
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


# Setting up acoustic transducer arrays
def prop_matrix(array_x, array_y, bottom_array_z, top_array_z, geo_arr, P0, k, r0, y):
  bottom_prop = []
  top_prop = []
  for z in range(point_num):
    bottom_dist_map = jnp.sqrt( jnp.power((array_x - geo_arr[y][z*3]), 2) + jnp.power((array_y - geo_arr[y][z*3+1]), 2) + jnp.power((bottom_array_z - geo_arr[y][z*3+2]), 2) )
    top_dist_map = jnp.sqrt( jnp.power((array_x - geo_arr[y][z*3]), 2) + jnp.power((array_y - geo_arr[y][z*3+1]), 2) + jnp.power((top_array_z - geo_arr[y][z*3+2]), 2) )

    bottom_sin_alpha_map = jnp.sqrt( (jnp.power((array_x - geo_arr[y][z*3]), 2) + jnp.power((array_y - geo_arr[y][z*3+1]), 2)) ) / bottom_dist_map
    top_sin_alpha_map = jnp.sqrt( (jnp.power((array_x - geo_arr[y][z*3]), 2) + jnp.power((array_y - geo_arr[y][z*3+1]), 2)) ) / top_dist_map
    
    bottom_amplitude_map = ( 2 * jv(1, k*r0*bottom_sin_alpha_map) * P0 / (k*r0*bottom_sin_alpha_map*bottom_dist_map) )
    top_amplitude_map = ( 2 * jv(1, k*r0*top_sin_alpha_map) * P0 / (k*r0*top_sin_alpha_map*top_dist_map) )
    bottom_prop.append( jax.lax.complex( bottom_amplitude_map * jnp.cos(k*bottom_dist_map), bottom_amplitude_map * jnp.sin(k*bottom_dist_map) ) )
    top_prop.append( jax.lax.complex( top_amplitude_map * jnp.cos(k*top_dist_map), top_amplitude_map * jnp.sin(k*top_dist_map) ) )
  
  bottom_prop = device_put(jnp.array(bottom_prop))
  top_prop = device_put(jnp.array(top_prop))

  return bottom_prop, top_prop

# Loss function of the optimization function
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

# Optimizer's step function
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

# Setting up the problem
N=16 # Array Size
size=[N,N]
pitch=0.0105 # Pitch Length in m
r0=0.005 # Transducer radius in m
P0=1.98 # Hirayama 1.95 V Pk-PK
l_ambda = 346.0 / 40000.0 # Acoustic Waevlength
k = 2.0*jnp.pi/l_ambda # Wavenumber

# Transducer position
array_x = jnp.array([ (i % N - (N//2-0.5))*pitch for i in range(N * N)]).reshape(N, N)
array_y = array_x.T
bottom_array_z = jnp.full((N,N), 0.0)
top_array_z = jnp.full((N,N), 0.23550056)

for x in range(5):
  point_num = jnp.power(2, x+1)  # number of control points
  fn_num = '{:0=3}'.format(point_num)
  #Load up target amplitudes
  amp_df = pd.read_csv('amplitudes_' + str(fn_num) + '.csv', header=None)
  amp_df = amp_df.drop(amp_df.columns[[0]], axis=1)
  amp_arr = amp_df.values
  #Load up target points
  geo_df = pd.read_csv('geometries_' + str(fn_num) + '.csv', header=None)
  geo_df = geo_df.drop(geo_df.columns[[0]], axis=1)
  geo_arr = geo_df.values

  geo_amount = len(geo_df)

  all_phase_to_csv = []
  loss_list = []
  calc_time_list = []
  for y in range(1000): #for every sample points
    start_time = time.time() # start recording time
    init_bottom_phase = np.random.uniform(0.0, 2.0*np.pi, size) # initial guess set to random
    init_top_phase = np.random.uniform(0.0, 2.0*np.pi, size) # initial guess set to random
    all_phase = ([init_bottom_phase, init_top_phase])
    # Get array propetries
    bottom_prop, top_prop = prop_matrix(array_x, array_y, bottom_array_z, top_array_z, geo_arr, P0, k, r0, y)
    target_amp = amp_arr[y]

    # Setting up optimizer 
    opt_init, opt_update, get_params = jax.experimental.optimizers.adam(0.1, b1=0.9, b2=0.999, eps=1e-08)
    opt_state = opt_init(all_phase)
    
    iteration_phases = []

    # Optimize
    for ii in range(150):
      value, opt_state = step(ii, opt_state, top_prop, bottom_prop, target_amp)

    # Output
    all_phase = get_params(opt_state)   
    calc_time = time.time() - start_time
    print("--- %s seconds ---" % calc_time)
    all_phase_1d_arr = np.array( [all_phase[0].T.ravel(), all_phase[1].T.ravel()] ).ravel()
    all_phase_to_csv.append(all_phase_1d_arr)

    calc_time_list.append(calc_time)
    print("--- Geo num: %d ---" % y)

  all_phase_to_csv = np.array(all_phase_to_csv)
  np.savetxt('phases_' + str(fn_num) + '.csv', all_phase_to_csv, delimiter=',')

  np.savetxt('calc_time_list_' + str(fn_num) + '.csv', calc_time_list, delimiter=',')