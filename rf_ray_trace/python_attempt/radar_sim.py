import numpy as np

import scipy

import cmath as math

#generate array for ionosphere
#electron_density, drift_velocity, power_refl, spectral_width x altitude x distance (1km resolution)
ionosphere_model = np.zeros((4,1000,2000),dtype=float)

#generate array for ray path
#fwd_IQ,rev_IQ,angle_of_prop,altitude,ground_range x sample distance (1km)
radar_ray_current = np.zeros((5,2000),dtype=complex)
radar_ray_future = np.zeros((5,2000),dtype=complex)

#angle radar ray is propagating
angle_of_propagation = float(10.0)

#resolution between samples
range_resolution = float(1.0)

#index of refraction
index_1 = float(1.0)
index_2 = float(1.0)

#variables for plasma freq
plasma_freq_1 = float(0.0)
plasma_freq_2 = float(0.0)

#location of radar for magnetic data
mag_lat_input = float(60.0)
mag_long_input = float(60.0)

#frequency of signal transmitted
input_frequency = float(10000000.0)

#input IQ values
in_IQ = np.zeros((1000),dtype=complex)

#set a test pulse
in_IQ[0:10] = 1+1j

#output IQ values
out_IQ = np.zeros((1000),dtype=complex)

#variables to hold the reflection effects on IQ
mag_reflection = float(0.0)
phase_reflection = float(0.0)

#constants
electron_charge = 1.602e-19
electron_mass = 9.109e-31
perm_free_space = 8.854e-12

#update ionosphere model
for altitude_counter in range(1000):

	#determine the electron density
	ionosphere_model[0,altitude_counter,:] = altitude_counter * 10000000000

#create radar ray path
for radar_range in range(1999):
	#convert the initial angle of propagation to radians
	radar_ray_current[2,0] = (angle_of_propagation*math.pi)/180.0

	#determine the location of the propagating ray
	radar_ray_current[3,radar_range+1] = radar_ray_current[3,radar_range] + range_resolution*math.sin(radar_ray_current[2,radar_range])
	radar_ray_current[4,radar_range+1] = radar_ray_current[4,radar_range] + range_resolution*math.cos(radar_ray_current[2,radar_range])

	#determine the plasma frequency of the present and future path
	plasma_freq_1 = math.sqrt(ionosphere_model[0,int(radar_ray_current[3,radar_range]),int(radar_ray_current[4,radar_range])]) * electron_charge / math.sqrt(perm_free_space*electron_mass)
	plasma_freq_2 = math.sqrt(ionosphere_model[0,int(radar_ray_current[3,radar_range+1]),int(radar_ray_current[4,radar_range+1])]) * electron_charge / math.sqrt(perm_free_space*electron_mass)

	#calculate the change in the index of refraction
	index_1 = 1.0 - ((plasma_freq_1*2*math.pi)**2/((input_frequency*2*math.pi)**2))
	index_2 = 1.0 - ((plasma_freq_2*2*math.pi)**2/((input_frequency*2*math.pi)**2))

	#direction of propagation
	radar_ray_current[2,radar_range+1] = math.asin((index_1/index_2)*math.sin(math.pi/2.0 - radar_ray_current[2,radar_range].real))

	print radar_ray_current[2,radar_range+1]

#loop through time
for sample_time in range(1000):

	#set the first sample to be what is input
	radar_ray_current[0,0] = in_IQ[sample_time]

	#loop through range
	for radar_range in range(1999):

		#determine reflection effects on radar wave
		#phase_reflection = exp(j*2*math.pi*velocity*time/wavelength)
		#mag_reflection = power of irregularities

		#fwd propagation
		radar_ray_future[0,radar_range+1] = (1-mag_reflection)*radar_ray_current[0,radar_range]

		#reverse propagation
		radar_ray_future[1,radar_range] = radar_ray_current[1,radar_range+1] + mag_reflection*radar_ray_current[0,radar_range]

	#record the incoming voltage sample into an array
	out_IQ[sample_time] = radar_ray_current[1,0]
	
	#set the future array to be the current array
	radar_ray_current[0,:] = radar_ray_future[0,:]
	radar_ray_current[1,:] = radar_ray_future[1,:]

print out_IQ
