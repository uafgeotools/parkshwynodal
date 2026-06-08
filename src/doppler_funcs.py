'''Inversion functions for Doppler shift data.'''
import numpy as np
import numpy.linalg as la

class DopplerCalc:
	def __init__(self, v, d0, t0, c):
		self.t0 = t0
		self.v = v
		self.d0 = d0
		self.c = c
		
	def calc_t(self, tprime):
		t0 = self.t0
		v = self.v
		d0 = self.d0
		c = self.c

		diff_t = tprime - t0
		vel_diff = 1 - v**2/c**2
		sqrt_term = np.sqrt((diff_t)**2 - (vel_diff)*((diff_t)**2 - d0**2/c**2))
		t = (diff_t - sqrt_term)/(vel_diff)

		return t

	def calc_fs(self, tprime, ft0p):
		"""
		Calculate the fundamental frequency produced by an aircraft 
		(where the wave is generated) given the model parameters.

		Parameters:
		tprime (float): Time at which a frequency (ft0p) is observed on the 
			station.
		t0 (float):  The time at which the central frequency of the overtones 
			occur.
		ft0p (float): Frequency recorded on the seismometer, picked from the 
			overtone doppler curve.
		v (float): Velocity of the aircraft.
		d0 (float): Distance between the station and the aircraft at the 
			closest approach.
		c (float): Speed of sound.

		Returns:
		fs (float): Source frequency produced by the aircraft. 
			(Frequency at the source.) 
		"""

		v = self.v
		d0 = self.d0
		c = self.c

		t = self.calc_t(tprime)
		fs = ft0p*(1+(v/c)*(v*t)/(np.sqrt(d0**2+(v*t)**2)))
		self.fs = fs
		return fs
	
	def calc_ft(self, times, fs = None):
		"""
		Calculate the frequency at each given time using the model parameters.

		Args:
			times (list): List of time values.
			t0 (float): The time at which the central frequency of the overtones 
				occur, when the aircraft is at the closest approach to the 
				station.
			fs (float): Source frequency produced by the aircraft.
			v (float): Velocity of the aircraft.
			d0 (float): Distance between the station and the aircraft at the 
				closest approach.
			c (float): Speed of sound.

		Returns:
			list: List of calculated frequency values.
		"""
		if fs is None:
			fs = self.fs
		v = self.v
		d0 = self.d0
		c = self.c

		ft = []
		for tprime in times:
			t = self.calc_t(tprime)
			ft0p = fs/(1+(v/c)*(v*t)/(np.sqrt(d0**2+(v*t)**2)))
									
			ft.append(ft0p)
		return ft
	

class DopplerInversion(DopplerCalc):

	def __init__(self, fobs, tobs, mprior, prior_sigma, num_iterations=4, 
			  off_diagonal=False):
		
		self.tobs = tobs
		self.fobs = fobs
		self.mprior = mprior
		self.prior_sigma = prior_sigma
		num_overtones = len(self.mprior[4:])
		self.num_overtones = num_overtones
		self.num_iterations = num_iterations
		self.off_diagonal = off_diagonal

		# internal variable to store predicted frequencies during iterations
		self.sigma = 10  # Default uncertainty in fs measurements
		self.fpred = None
		self.cprior = None
		self.source_closest_approach = None
		self.sound_speed = None
		self.closest_approach_dist = None
		self.source_speed = None
		self.source_frequencies = None
		
		self.num_overtones = len(self.mprior[4:])

	def cprior_setup(self):
		'''
		Setup the prior model covariance matrix.
		Parameters:
		prior_sigma (list): List of standard deviations for the prior model
			parameters. The order should be [v_sigma, d0_sigma,
			t0_sigma, c_sigma, fs_sigma].'''

		v_sigma = self.prior_sigma[0]
		d0_sigma = self.prior_sigma[1]
		t0_sigma = self.prior_sigma[2]
		c_sigma = self.prior_sigma[3]
		fs_sigma = self.prior_sigma[4]

		cprior0 = np.zeros((len(self.mprior), len(self.mprior)))
	
		cprior0[0][0] = v_sigma**2
		cprior0[1][1] = d0_sigma**2
		cprior0[2][2] = t0_sigma**2
		cprior0[3][3] = c_sigma**2

		for row in range(len(cprior0)):
			if row >= 4:
				cprior0[row][row] = fs_sigma**2

		if self.off_diagonal:
			cprior0[4:][2] =  -0.4 * fs_sigma * t0_sigma
			cprior0[0][1] = -0.7 * v_sigma * d0_sigma
			cprior0[0][3] = 0.85 * v_sigma * c_sigma
			cprior0[1][0] = -0.7 * v_sigma * d0_sigma
			cprior0[1][3] = -0.7 * d0_sigma * c_sigma
			cprior0[2][4:] =  -0.4 * fs_sigma * t0_sigma
			cprior0[3][0] = 0.85 * v_sigma * c_sigma
			cprior0[3][1] = -0.7 * d0_sigma * c_sigma

		cprior = cprior0 * (len(self.mprior))

		Cd0 = np.zeros((len(self.fobs), len(self.fobs)), int)
		np.fill_diagonal(Cd0, self.sigma**2)
		Cd = Cd0*(len(self.fobs))

		mnew = np.array(self.mprior)

		return cprior0, cprior, Cd0, Cd, mnew
						
	def df(self, tp):   
		'''
		Calculate the derivatives of f with respect to fs, v, d0, t0 and c.

		Parameters:
		fs (float): Fundamental frequency produced by the aircraft.
		v (float): Velocity of the aircraft.
		d0 (float): Distance of closest approach between the station and 
			the aircraft.
		tp0 (float): Time of that the central frequency of the overtones occur, 
			when the aircraft is at the closest approach to the station.
		c (float): Speed of sound.
		tp (numpy.ndarray): Array of times.

		Returns:
		tuple: A tuple containing the derivatives of f with respect to fs, v, 
			d0, t0 and c.
		'''

		fs = self.source_frequencies
		v = self.source_speed
		d0 = self.closest_approach_dist
		tp0 = self.source_closest_approach
		c = self.sound_speed

		# Pre-compute common subexpressions
		delta_t = tp - tp0
		v_ratio = v / c
		v_ratio_sq = v_ratio ** 2
		
		sqrt_term = np.sqrt((-d0**2 * v**2 + c**2 * (d0**2 + delta_t**2 * v**2)) 
					/ c**4)
		
		denom_term = ((-tp + tp0) * v**2 + c**2 * sqrt_term)
		
		sqrt_l_term = np.sqrt(d0**2 + (c**4 * v**2 * (-tp + tp0 + sqrt_term)**2)
						/ (c**2 - v**2)**2)

		# Derivative with respect to fs
		f_derivefs = (1 / (1 - (c * v_ratio_sq * denom_term) / ((c**2 - v**2) 
					* sqrt_l_term)))

		# Derivative of f with respect to v
		numerator_v = (-fs * v * (-2 * d0**4 * v**4 + d0**2 * delta_t**2 * v**6 
				+ c**6 * delta_t * (2 * d0**2 + delta_t**2 * v**2) * sqrt_term 
				+ c**2 * (4 * d0**4 * v**2 - delta_t**4 * v**6 + d0**2 * delta_t 
			  	* v**4 * (5 * delta_t - 3 * sqrt_term)) - c**4 * (2 * d0**4 
				- 3 * delta_t**3 * v**4 * (-tp + tp0 + sqrt_term) - d0**2 
				* delta_t * v**2 * (-6 * delta_t + sqrt_term))))
		
		denominator_v = (c * (c - v) * (c + v) * sqrt_term * 
			sqrt_l_term * (c * (-tp + tp0) * v**2 + c * v**2 * sqrt_term - 
				c**2 * sqrt_l_term + v**2 * sqrt_l_term)**2)
		
		f_derivev = numerator_v / denominator_v

		# Derivative of f with respect to d0
		numerator_d0 = (fs * d0 * delta_t * (c - v) * v**2 * (c + v) 
				 * denom_term)
		
		denominator_d0 = (c * sqrt_term * sqrt_l_term * (c * (-tp + tp0) * v**2 
				+ c * v**2 * sqrt_term - c**2 * sqrt_l_term + v**2 
				* sqrt_l_term)**2)
		
		f_derived0 = numerator_d0 / denominator_d0

		# Derivative of f with respect to tp0
		numerator_t0 = (fs * d0**2 * (c - v) * v**2 * (c + v) * denom_term)
		
		denominator_t0 = (
			c * sqrt_term * sqrt_l_term * (c * (-tp + tp0) * v**2 + 
				c * v**2 * sqrt_term - c**2 * sqrt_l_term + 
				v**2 * sqrt_l_term)**2)
		
		f_derivet0 = numerator_t0 / denominator_t0

		# Derivative of f with respect to c
		numerator_c = (
			fs * v**2 * (-2 * d0**4 * v**4 + 2 * d0**2 * delta_t**2 * v**6 + 
				c**6 * delta_t * (d0**2 + delta_t**2 * v**2) * sqrt_term + 
				c**2 * (4 * d0**4 * v**2 - delta_t**4 * v**6 + d0**2 * delta_t *
			 	v**4 * (3 * delta_t - 4 * sqrt_term)) - c**4 * 
				(d0**2 + delta_t**2 * v**2) * (2 * d0**2 - 3 * delta_t * v**2 * 
				(-tp + tp0 + sqrt_term))))
		
		denominator_c = (
			c**2 * (c - v) * (c + v) * sqrt_term * sqrt_l_term * (
				c * (-tp + tp0) * v**2 + c * v**2 * sqrt_term - 
				c**2 * sqrt_l_term + v**2 * sqrt_l_term)**2)
		
		f_derivec = numerator_c / denominator_c
		
		return f_derivev, f_derived0, f_derivet0, f_derivec, f_derivefs


	def data_misfit(self):
		"""
		Calculate the data misfit using the predictions and observations.
		MISFIT FUNCTION: least squares, Tarantola (2005), Eq. 6.251

		Args:
			dnew (array): Array of predicted data.
			dobs (array): Array of observed data.
			ndata (int): Number of data points.
			m (array): Posterior model parameters.
			mprior (array): Prior model parameters.
			cprior (array): Covariance matrix for prior model.
			tsigma (float): Uncertainty in fs measurements, Hz

		Returns:
			float: Data misfit value.
		"""
		dnew = self.fpred
		dobs = self.fobs
		ndata = len(self.fobs)
		m = np.array(self.mnew)
		tsigma = self.sigma
		mprior = np.array(self.mprior)
		
		sigma_obs = tsigma * np.ones((ndata))
		cobs0 = np.diag(np.square(sigma_obs))
		Cdfac = ndata
		dnew = np.array(dnew)
		dobs = np.array(dobs)
		cobs = Cdfac * cobs0
		icobs = la.inv(cobs)
		icprior = la.inv(self.cprior)
		print(dnew)
		Sd = 0.5 * (dnew - dobs).T @ icobs @ (dnew - dobs)
		Sm = 0.5 * (m - mprior).T @ icprior @ (m - mprior)
		S = Sd + Sm

		print("Model Misfit:", Sm)
		print("Data Misfit:", Sd)
		print("Total Misfit:", S)
		return Sd

	
	def full_inversion(self, peaks_assos, sigma=3):
		"""
		Performs inversion using all picked overtones. 

		Args:
			fobs (numpy.ndarray): Picked frequency values from individual 
				overtone inversion picks.
			tobs (numpy.ndarray): Picked time values from individual overtone 
				inversion picks.
			peak_assos (list): List of number of peaks associated with each 
				overtone, for indexing the fobs and tobs arrays.
			mprior (numpy.ndarray): Initial guess for the model parameters, 
			mprior[0] = v, mprior[1] = d0, mprior[2] = t0, mprior[3] = c, 
			mprior[4:] = fs_array.
			num_iterations (int): Number of iterations to perform for the 
				inversion.
			sigma (float): Standard deviation for the data picks, default is 3.
			off_diagonal (bool): Whether to include off-diagonal elements 
				in the prior covariance matrix, default is False.

		Returns:
			numpy.ndarray: The inverted parameters for the function f. Velocity 
				of the aircraft, distance of closest approach, time of closest 
				approach, and the fundamental frequency produced by the 
				aircraft.
			numpy.ndarray: The covariance matrix of the inverted parameters.
			numpy.ndarray: The array of the fundamental frequency produced 
				by the aircraft.
		"""
		
		cprior0, cprior, Cd0, Cd, mnew = self.cprior_setup()
		self.cprior = cprior

		self.sigma = sigma

		qv = 0

		while qv < self.num_iterations:
			if np.any(np.isnan(mnew)) and qv == 0:
				return (self.mprior, cprior0, cprior, self.mprior[4:], 
						'Forward Model')
			elif np.any(np.isnan(mnew)):
				mnew = m
				G = G_hold
				Cpost = la.inv(G.T @ la.pinv(Cd) @ G + la.inv(cprior))
				Cpost0 = la.inv(G.T @ la.pinv(Cd0) @ G + la.inv(cprior0))
				return mnew, Cpost0, Cpost, fs_array, self.data_misfit()
			else:
				m = mnew
			v = m[0]
			d0 = m[1]
			t0 = m[2]
			c = m[3]
			fs_array = m[4:]

			self.source_speed = v
			self.closest_approach_dist = d0
			self.source_closest_approach = t0
			self.sound_speed = c

			fpred = []
			G = np.zeros((0, self.num_overtones + 4))
			cum = 0
			for p in range(self.num_overtones):
				new_row = np.zeros(self.num_overtones + 4)
				fs = fs_array[p]
				self.source_frequencies = fs
				
				for j in range(cum, cum + peaks_assos[p]):
					tprime = self.tobs[j]
					t = ((tprime - t0) - np.sqrt((tprime - t0)**2 - 
							(1 - v**2 / c**2) * ((tprime - t0)**2 - 
								d0**2 / c**2))) / (1 - v**2 / c**2)
					ft0p = fs / (1 + (v / c) * (v * t) / (
							np.sqrt(d0**2 + (v * t)**2)))

					f_derivev, f_derived0, f_derivet0, f_derivec, f_derivefs \
						 = self.df(tprime)
					
					new_row[0] = f_derivev
					new_row[1] = f_derived0
					new_row[2] = f_derivet0
					new_row[3] = f_derivec
					new_row[4 + p] = f_derivefs
					G = np.vstack((G, new_row))
							
					fpred.append(ft0p)
			
				cum = cum + peaks_assos[p]

			Gm = G
			self.fpred = fpred
			gamma = (cprior @ Gm.T @ la.inv(Cd) @ (np.array(fpred) - self.fobs) 
					+ (np.array(m) - np.array(self.mprior)))
			H = (np.identity(len(mnew)) + cprior @ Gm.T @ la.inv(Cd) @ Gm)
			dm = -la.inv(H) @ gamma
			mnew = m + dm
			self.mnew = mnew

			unreasonable = (
				[mn for mn in mnew[4:] 
					if mn <= 5 or mn > 375] or
				mnew[0] <= 0 or mnew[0] > 350 or
				mnew[0] >= mnew[3] or
				mnew[1] < 0 or mnew[1] > 1e5 or
				mnew[2] < 10 or mnew[2] > 240 or
				mnew[3] < 200 or mnew[3] > 400
			)

			if unreasonable and qv > 0:
				mnew = m
				G = G_hold
				Cpost = la.inv(G.T @ la.pinv(Cd) @ G + la.inv(cprior))
				Cpost0 = la.inv(G.T @ la.pinv(Cd0) @ G + la.inv(cprior0))

				return mnew, Cpost0, Cpost, fs_array, self.data_misfit()
			
			elif unreasonable and qv == 0:
				return (self.mprior, cprior0, cprior, self.mprior[4:], 
					'Forward Model'
				)
			elif np.nan in mnew:
				return (self.mprior, cprior0, cprior, self.mprior[4:], 
					'Forward Model')
			else:
				G_hold = G.copy()
			fs_array = m[4:]
			qv += 1


		Cpost = la.inv(Gm.T @ la.inv(Cd) @ Gm + la.inv(cprior))
		Cpost0 = la.inv(Gm.T @ la.inv(Cd0) @ Gm + la.inv(cprior0))
		F_m = self.data_misfit()

		return mnew, Cpost0, Cpost, fs_array, F_m
	
	def main(self):
		if self.method == 'full':
			return self.full_inversion(self.peaks_assos, 3)
		else:
			return self.full_inversion(self.peaks_assos, 10)
