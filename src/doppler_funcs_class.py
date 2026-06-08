'''Inversion functions for Doppler shift data.'''
import numpy as np
import numpy.linalg as la
'''Please comment on how well the code adheres to PEP8 standards. Also, please
 comment on places where the code could be made more efficient and have less 
 lines'''
class DopplerInversion:

	def __init__(self, fobs, tobs, mprior, prior_sigma, num_iterations=4,
			  method='full', off_diagonal=False):
		self.fobs = fobs
		self.tobs = tobs
		self.mprior = mprior
		self.prior_sigma = prior_sigma
		num_overtones = len(self.mprior[4:])
		self.num_overtones = num_overtones
		self.num_iterations = num_iterations
		self.method = method
		self.off_diagonal = off_diagonal

		# internal variable to store predicted frequencies during iterations
		self.sigma = 10  # Default uncertainty in f0 measurements
		self.fpred = None
		self.cprior = None
		self.source_closest_approach = None
		self.sound_speed = None
		self.closest_approach_dist = None
		self.source_speed = None
		self.source_frequencies = None
		
		self.num_overtones = len(self.mprior[4:])

	def cprior_setup(self):
		f0_sigma = self.prior_sigma[0]
		v0_sigma = self.prior_sigma[1]
		l_sigma = self.prior_sigma[2]
		t0_sigma = self.prior_sigma[3]
		c_sigma = self.prior_sigma[4]

		cprior0 = np.zeros((len(self.mprior), len(self.mprior)))
	
		cprior0[0][0] = v0_sigma**2
		cprior0[1][1] = l_sigma**2
		cprior0[2][2] = t0_sigma**2
		cprior0[3][3] = c_sigma**2

		for row in range(len(cprior0)):
			if row >= 4:
				cprior0[row][row] = f0_sigma**2

		if self.off_diagonal:
			cprior0[4:][2] =  -0.4 * f0_sigma * t0_sigma
			cprior0[0][1] = -0.7 * v0_sigma * l_sigma
			cprior0[0][3] = 0.85 * v0_sigma * c_sigma
			cprior0[1][0] = -0.7 * v0_sigma * l_sigma
			cprior0[1][3] = -0.7 * l_sigma * c_sigma
			cprior0[2][4:] =  -0.4 * f0_sigma * t0_sigma
			cprior0[3][0] = 0.85 * v0_sigma * c_sigma
			cprior0[3][1] = -0.7 * l_sigma * c_sigma

		cprior = cprior0 * (len(self.mprior))

		Cd0 = np.zeros((len(self.fobs), len(self.fobs)), int)
		np.fill_diagonal(Cd0, self.sigma**2)
		Cd = Cd0*(len(self.fobs))

		mnew = np.array(self.mprior)

		return cprior0, cprior, Cd0, Cd, mnew
						
	def df(self, tp):   
		f0 = self.source_frequencies
		v0 = self.source_speed
		l = self.closest_approach_dist
		tp0 = self.source_closest_approach
		c = self.sound_speed

		# Pre-compute common subexpressions
		delta_t = tp - tp0
		v_ratio = v0 / c
		v_ratio_sq = v_ratio ** 2
		
		sqrt_term = np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + delta_t**2 * v0**2)) 
					/ c**4)
		
		denom_term = ((-tp + tp0) * v0**2 + c**2 * sqrt_term)
		
		sqrt_l_term = np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + sqrt_term)**2)
						/ (c**2 - v0**2)**2)

		# Derivative with respect to f0
		f_derivef0 = (1 / (1 - (c * v_ratio_sq * denom_term) / ((c**2 - v0**2) 
					* sqrt_l_term)))

		# Derivative of f with respect to v0
		numerator_v0 = (-f0 * v0 * (-2 * l**4 * v0**4 + l**2 * delta_t**2 * v0**6 
				+ c**6 * delta_t * (2 * l**2 + delta_t**2 * v0**2) * sqrt_term 
				+ c**2 * (4 * l**4 * v0**2 - delta_t**4 * v0**6 + l**2 * delta_t 
			  	* v0**4 * (5 * delta_t - 3 * sqrt_term)) - c**4 * (2 * l**4 
				- 3 * delta_t**3 * v0**4 * (-tp + tp0 + sqrt_term) - l**2 
				* delta_t * v0**2 * (-6 * delta_t + sqrt_term))))
		
		denominator_v0 = (c * (c - v0) * (c + v0) * sqrt_term * 
			sqrt_l_term * (c * (-tp + tp0) * v0**2 + c * v0**2 * sqrt_term - 
				c**2 * sqrt_l_term + v0**2 * sqrt_l_term)**2)
		
		f_derivev0 = numerator_v0 / denominator_v0

		# Derivative of f with respect to l
		numerator_l = (f0 * l * delta_t * (c - v0) * v0**2 * (c + v0) 
				 * denom_term)
		
		denominator_l = (c * sqrt_term * sqrt_l_term * (c * (-tp + tp0) * v0**2 
				+ c * v0**2 * sqrt_term - c**2 * sqrt_l_term + v0**2 
				* sqrt_l_term)**2)
		
		f_derivel = numerator_l / denominator_l

		# Derivative of f with respect to tp0
		numerator_t0 = (f0 * l**2 * (c - v0) * v0**2 * (c + v0) * denom_term)
		
		denominator_t0 = (
			c * sqrt_term * sqrt_l_term * (c * (-tp + tp0) * v0**2 + 
				c * v0**2 * sqrt_term - c**2 * sqrt_l_term + 
				v0**2 * sqrt_l_term)**2)
		
		f_derivet0 = numerator_t0 / denominator_t0

		# Derivative of f with respect to c
		numerator_c = (
			f0 * v0**2 * (-2 * l**4 * v0**4 + 2 * l**2 * delta_t**2 * v0**6 + 
				c**6 * delta_t * (l**2 + delta_t**2 * v0**2) * sqrt_term + 
				c**2 * (4 * l**4 * v0**2 - delta_t**4 * v0**6 + l**2 * delta_t *
			 	v0**4 * (3 * delta_t - 4 * sqrt_term)) - c**4 * 
				(l**2 + delta_t**2 * v0**2) * (2 * l**2 - 3 * delta_t * v0**2 * 
				(-tp + tp0 + sqrt_term))))
		
		denominator_c = (
			c**2 * (c - v0) * (c + v0) * sqrt_term * sqrt_l_term * (
				c * (-tp + tp0) * v0**2 + c * v0**2 * sqrt_term - 
				c**2 * sqrt_l_term + v0**2 * sqrt_l_term)**2)
		
		f_derivec = numerator_c / denominator_c
		
		return f_derivef0, f_derivev0, f_derivel, f_derivet0, f_derivec


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
			tsigma (float): Uncertainty in f0 measurements, Hz

		Returns:
			float: Data misfit value.
		"""
		dnew = self.fpred
		dobs = self.fobs
		ndata = len(self.fobs)
		m = self.mnew
		tsigma = self.prior_sigma
		mprior = np.array(self.mprior)
		
		sigma_obs = tsigma * np.ones((ndata))
		cobs0 = np.diag(np.square(sigma_obs))
		m = np.array(m)
		mprior = np.array(self.mprior)
		Cdfac = ndata
		dnew = np.array(dnew)
		dobs = np.array(dobs)
		cobs = Cdfac * cobs0
		icobs = la.inv(cobs)
		icprior = la.inv(self.cprior)

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
			mprior[0] = v0, mprior[1] = l, mprior[2] = t0, mprior[3] = c, 
			mprior[4:] = f0_array.
			num_iterations (int): Number of iterations to perform for the 
			inversion.
			sigma (float): Standard deviation for the data picks, default is 3.
			off_diagonal (bool): Whether to include off-diagonal elements in the
			prior covariance matrix, default is False.

		Returns:
			numpy.ndarray: The inverted parameters for the function f. Velocity 
			of the aircraft, distance of closest approach, time of closest 
			approach, and the fundamental frequency produced by the aircraft.
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
				return mnew, Cpost0, Cpost, f0_array, self.data_misfit()
			else:
				m = mnew
			v0 = m[0]
			l = m[1]
			t0 = m[2]
			c = m[3]
			f0_array = m[4:]

			self.source_speed = v0
			self.closest_approach_dist = l
			self.source_closest_approach = t0
			self.sound_speed = c

			fpred = []
			G = np.zeros((0, self.num_overtones + 4))
			cum = 0
			for p in range(self.num_overtones):
				new_row = np.zeros(self.num_overtones + 4)
				f0 = f0_array[p]
				self.source_frequencies = f0
				
				for j in range(cum, cum + peaks_assos[p]):
					tprime = self.tobs[j]
					t = ((tprime - t0) - np.sqrt((tprime - t0)**2 - 
							(1 - v0**2 / c**2) * ((tprime - t0)**2 - 
								l**2 / c**2))) / (1 - v0**2 / c**2)
					ft0p = f0 / (1 + (v0 / c) * (v0 * t) / (
							np.sqrt(l**2 + (v0 * t)**2)))

					f_derivef0, f_derivev0, f_derivel, f_derivet0, f_derivec \
						 = self.df(tprime)
					
					new_row[0] = f_derivev0
					new_row[1] = f_derivel
					new_row[2] = f_derivet0
					new_row[3] = f_derivec
					new_row[4 + p] = f_derivef0
					G = np.vstack((G, new_row))
							
					fpred.append(ft0p)
			
				cum = cum + peaks_assos[p]

			Gm = G
			
			gamma = (cprior @ Gm.T @ la.inv(Cd) @ (np.array(fpred) - self.fobs) 
					+ (np.array(m) - np.array(self.mprior)))
			H = (np.identity(len(mnew)) + cprior @ Gm.T @ la.inv(Cd) @ Gm)
			dm = -la.inv(H) @ gamma
			mnew = m + dm

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

				return mnew, Cpost0, Cpost, f0_array, self.data_misfit()
			
			elif unreasonable and qv == 0:
				return (self.mprior, cprior0, cprior, self.mprior[4:], 
					'Forward Model'
				)
			elif np.nan in mnew:
				return (self.mprior, cprior0, cprior, self.mprior[4:], 
					'Forward Model')
			else:
				G_hold = G.copy()
			f0_array = m[4:]
			qv += 1
			print(mnew)

		Cpost = la.inv(Gm.T @ la.inv(Cd) @ Gm + la.inv(cprior))
		Cpost0 = la.inv(Gm.T @ la.inv(Cd0) @ Gm + la.inv(cprior0))
		F_m = self.data_misfit()

		return mnew, Cpost0, Cpost, mnew[4:], F_m
	
	def main(self):
		if self.method == 'full':
			return self.full_inversion(self.peaks_assos, 3)
		else:
			self.peaks_assos = [1] * self.num_overtones
			return self.full_inversion(self.peaks_assos, 10)
