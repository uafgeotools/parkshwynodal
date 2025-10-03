import numpy as np
import os
import gc
from obspy.clients.nrl import NRL
from scipy.signal import spectrogram, find_peaks
from src.doppler_funcs import make_base_dir, invert_f, full_inversion, get_sta_elevation, load_waveform
from src.main_inv_fig_functions import time_picks, remove_median, plot_spectrogram, plot_spectrum, get_auto_picks_full
import psutil
import numpy.linalg as la
from src.main_inv_fig_functions import remove_median
from src.doppler_funcs import S

def calc_ft(tpr, tprime0, f0, v0, l, c):
    """
    Calculate the frequency at each given time using the model parameters.

    Args:
        tpr (list): List of time values.
        tprime0 (float): The time at which the central frequency of the overtones occur, 
                        when the aircraft is at the closest approach to the station.
        f0 (float): Fundamental frequency produced by the aircraft.
        v0 (float): Velocity of the aircraft.
        l (float): Distance between the station and the aircraft at the closest approach.
        c (float): Speed of sound.

    Returns:
        list: List of calculated frequency values.
    """
    ft = []
    times = calc_t(v0,l,c,tprime0,tpr)
    for t in times:
        f = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
        ft.append(f)
    return ft

def calc_t(v0,l,c,tprime0,tpr):
	t_array = []	
	for tprime in tpr:
		t = ((tprime - tprime0 + l/c)- np.sqrt((tprime - tprime0 + l/c)**2-(1-v0**2/c**2)*((tprime - tprime0 + l/c)**2-l**2/c**2)))/(1-v0**2/c**2)
		t_array.append(t)
	return t_array

def calc_f0(tprime, tprime0, ft0p, v0, l, c):
    """
    Calculate the fundamental frequency produced by an aircraft (where the wave is generated) given the model parameters.

    Parameters:
    tprime (float): Time at which a frequency (ft0p) is observed on the station.
    tprime0 (float): The time at which the central frequency of the overtones occur.
    ft0p (float): Frequency recorded on the seismometer, picked from the overtone doppler curve.
    v0 (float): Velocity of the aircraft.
    l (float): Distance between the station and the aircraft at the closest approach.
    c (float): Speed of sound.

    Returns:
    f0 (float): Fundamental frequency produced by the aircraft. (Frequency at the source.) 
    """
    t = calc_t(v0,l,c,tprime0,[tprime])[0]
    
    # Geometry now makes sense
    f0 = ft0p * (1 + (v0/c) * (v0*t) / np.sqrt(l**2 + (v0*t)**2))

    return f0

def df(f0, v0, l, tprime0, tprime, c):   
    """
    Calculate the derivatives of f with respect to f0, v0, l, tprime0 and c.
    Uses CORRECTED timing that ensures t=0 at closest approach.
    """

    f_derivef0 = 1/(1 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c*(1 - v0**2/c**2)*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)))
    f_derivev0 = f0*(-v0**2*(-v0*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2 + v0**3*((tprime - tprime0 + l/c)**2 - l**2/c**2)*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c**2*(1 - v0**2/c**2)**2*np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2)) - 2*v0**3*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(c**2*(1 - v0**2/c**2)**3))*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c*(1 - v0**2/c**2)*(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)**(3/2)) - 2*v0*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c*(1 - v0**2/c**2)*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)) + v0**3*((tprime - tprime0 + l/c)**2 - l**2/c**2)/(c**3*(1 - v0**2/c**2)*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)*np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2)) - 2*v0**3*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c**3*(1 - v0**2/c**2)**2*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)))/(1 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c*(1 - v0**2/c**2)*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)))**2
    f_derivel = f0*(-v0**2*(-l - v0**2*(-2*(-(1 - v0**2/c**2)*(2*(tprime - tprime0 + l/c)/c - 2*l/c**2)/2 + (tprime - tprime0 + l/c)/c)/np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + 2/c)*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(2*(1 - v0**2/c**2)**2))*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c*(1 - v0**2/c**2)*(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)**(3/2)) - v0**2*(-(-(1 - v0**2/c**2)*(2*(tprime - tprime0 + l/c)/c - 2*l/c**2)/2 + (tprime - tprime0 + l/c)/c)/np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + 1/c)/(c*(1 - v0**2/c**2)*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)))/(1 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c*(1 - v0**2/c**2)*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)))**2
    f_derivetprime0 = f0*(v0**4*(-2 - 2*(-tprime + tprime0 - (1 - v0**2/c**2)*(-2*tprime + 2*tprime0 - 2*l/c)/2 - l/c)/np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2))*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(2*c*(1 - v0**2/c**2)**3*(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)**(3/2)) - v0**2*(-1 - (-tprime + tprime0 - (1 - v0**2/c**2)*(-2*tprime + 2*tprime0 - 2*l/c)/2 - l/c)/np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2))/(c*(1 - v0**2/c**2)*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)))/(1 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c*(1 - v0**2/c**2)*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)))**2
    f_derivec = f0*(-v0**2*(-(-(1 - v0**2/c**2)*(-2*l*(tprime - tprime0 + l/c)/c**2 + 2*l**2/c**3)/2 - l*(tprime - tprime0 + l/c)/c**2 - v0**2*((tprime - tprime0 + l/c)**2 - l**2/c**2)/c**3)/np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) - l/c**2)/(c*(1 - v0**2/c**2)*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)) - v0**2*(-v0**2*(-2*(-(1 - v0**2/c**2)*(-2*l*(tprime - tprime0 + l/c)/c**2 + 2*l**2/c**3)/2 - l*(tprime - tprime0 + l/c)/c**2 - v0**2*((tprime - tprime0 + l/c)**2 - l**2/c**2)/c**3)/np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) - 2*l/c**2)*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(2*(1 - v0**2/c**2)**2) + 2*v0**4*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(c**3*(1 - v0**2/c**2)**3))*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c*(1 - v0**2/c**2)*(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)**(3/2)) + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c**2*(1 - v0**2/c**2)*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)) + 2*v0**4*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c**4*(1 - v0**2/c**2)**2*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)))/(1 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)/(c*(1 - v0**2/c**2)*np.sqrt(l**2 + v0**2*(tprime - tprime0 - np.sqrt(-(1 - v0**2/c**2)*((tprime - tprime0 + l/c)**2 - l**2/c**2) + (tprime - tprime0 + l/c)**2) + l/c)**2/(1 - v0**2/c**2)**2)))**2
    
    return f_derivef0, f_derivev0, f_derivel, f_derivetprime0, f_derivec

def invert_f(mprior, prior_sigma, coords_array, num_iterations = 4, sigma = 10, off_diagonal = False):
    """
    Inverts the function f using the given initial parameters and data array.
    Uses CORRECTED timing that ensures t=0 at closest approach.

    Args:
        mprior (numpy.ndarray): Initial parameters for the function, f[0] = f0, f[1] = v0, f[2] = l, f[3] = tprime0, f[4] = c.
        prior_sigma (list): List of standard deviations for the prior parameters prior_sigma[0] = f0_sigma, prior_sigma[1] = v0_sigma, prior_sigma[2] = l_sigma, prior_sigma[3] = tprime0_sigma, prior_sigma[4] = c_sigma.
        coords_array (numpy.ndarray): Data picks along overtone doppler curve.
        num_iterations (int): Number of iterations to perform.
        sigma (float): Standard deviation for the data picks, default is 10.
        off_diagonal (bool): Whether to include off-diagonal elements in the prior covariance matrix, default is False.

    Returns:
        numpy.ndarray: The inverted parameters for the function f.
        numpy.ndarray: The covariance matrix of the posterior parameters.
        numpy.ndarray: The normalized covariance matrix of the posterior parameters.
        float: The data misfit value.
    """

    dw,_ = coords_array.shape
    fobs = coords_array[:,1]
    tobs = coords_array[:,0]
    n = 0
 
    cprior0 = np.zeros((5,5))
    f0_sigma = prior_sigma[0]
    v0_sigma = prior_sigma[1]
    l_sigma = prior_sigma[2]
    tprime0_sigma = prior_sigma[3]
    c_sigma = prior_sigma[4]

    cprior0[0][0] = f0_sigma**2
    cprior0[1][1] = v0_sigma**2
    cprior0[2][2] = l_sigma**2
    cprior0[3][3] = tprime0_sigma**2
    cprior0[4][4] = c_sigma**2
    if off_diagonal:
        cprior0[0][3] =  -0.4*f0_sigma*tprime0_sigma

        cprior0[1][2] = -0.7*v0_sigma*l_sigma
        cprior0[1][4] = 0.85*v0_sigma*c_sigma

        cprior0[2][1] = -0.7*v0_sigma*l_sigma
        cprior0[2][4] = -0.7*l_sigma*c_sigma

        cprior0[3][0] =  -0.4*f0_sigma*tprime0_sigma

        cprior0[4][1] = 0.85*v0_sigma*c_sigma
        cprior0[4][2] = -0.7*l_sigma*c_sigma

    cprior = cprior0 * (5)
    Cd0 = np.zeros((len(fobs), len(fobs)), int)
    np.fill_diagonal(Cd0, sigma**2)
    Cd = Cd0*(dw)
    mnew = mprior.copy() #mprior is the initial guess for the parameters, mnew is the updated guess

    while n < num_iterations:
        if np.any(np.isnan(mnew)) and n == 0:
            # Handle the case where mnew contains NaN values
            return mprior, cprior0, cprior, 'Forward Model'
        elif np.any(np.isnan(mnew)):
            mnew = m
            G = G_hold
            Cpost = la.inv(G.T@la.pinv(Cd)@G + la.inv(cprior))
            Cpost0 = la.inv(G.T@la.pinv(Cd0)@G + la.inv(cprior0))
            return mnew, Cpost0, Cpost, S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)
        else:
            m = mnew
            
        f0 = m[0]
        v0 = m[1]
        l = m[2]
        tprime0 = m[3]
        c = m[4]

        fpred = []
        G = np.zeros((dw,5)) 

        # CORRECTED: Use new timing and derivatives
        for i in range(0,dw):
            tprime = tobs[i]

            # CORRECTED frequency calculation
            ft0p = calc_ft([tprime], tprime0, f0, v0, l, c)[0]

            # CORRECTED derivatives
            f_derivef0, f_derivev0, f_derivel, f_derivetprime0, f_derivec = df(f0, v0, l, tprime0, tprime, c)
            
            G[i,0:5] = [f_derivef0, f_derivev0, f_derivel, f_derivetprime0, f_derivec]
            fpred.append(ft0p) 
            
        Gm = G
        
        # steepest ascent vector (Eq. 6.307 or 6.312)
        gamma = cprior @ Gm.T @ la.inv(Cd) @ (np.array(fpred) - fobs) + (np.array(m)  - np.array(mprior)) # steepest ascent vector
        #===================================================
        # QUASI-NEWTON ALGORITHM (Eq. 6.319, nu=1)
        # approximate curvature
        H = np.identity(len(mnew)) + cprior @ Gm.T @ la.inv(Cd) @ Gm
        dm = -la.inv(H) @ gamma
        mnew = m + dm

        # Check for unreasonable parameter values
        unreasonable = (
            mnew[0] <= 5 or mnew[0] > 375 or    # f0
            mnew[1] <= 0 or mnew[1] > 350 or     # v0
            mnew[1] >= mnew[4] or  # v0 must be less than c
            mnew[2] < 0 or mnew[2] > 1e5 or      # l
            mnew[3] < 10 or mnew[3] > 240 or      # tprime0
            mnew[4] < 200 or mnew[4] > 400       # c
        )
        
        if unreasonable and n > 0:
            mnew = m
            G = G_hold
            Cpost = la.inv(G.T@la.pinv(Cd)@G + la.inv(cprior))
            Cpost0 = la.inv(G.T@la.pinv(Cd0)@G + la.inv(cprior0))
            return mnew, Cpost0, Cpost, S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)
        elif unreasonable and n == 0:
            return mprior, cprior0, cprior, 'Forward Model'
        elif np.any(np.isnan(mnew)):
            return mprior, cprior0, cprior, 'Forward Model'
        else:
            G_hold = G.copy()
            n += 1
        print(f"Iteration {n}: {mnew}")

    Cpost = la.inv(G.T@la.pinv(Cd)@G + la.inv(cprior))
    Cpost0 = la.inv(G.T@la.pinv(Cd0)@G + la.inv(cprior0))
    F_m = S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)
    del G, G_hold, Gm, H, dm, gamma, fpred, m, n, Cd, Cd0, cprior, cprior0
    return mnew, Cpost0, Cpost, F_m
def full_inversion(fobs, tobs, peaks_assos, mprior, sigma_prior, num_iterations = 4, sigma = 3, off_diagonal = False):
    """
    Performs inversion using all picked overtones. 
    Uses CORRECTED timing that ensures t=0 at closest approach.

    Args:
        fobs (numpy.ndarray): Picked frequency values from individual overtone inversion picks.
        tobs (numpy.ndarray): Picked time values from individual overtone inversion picks.
        peak_assos (list): List of number of peaks associated with each overtone, for indexing the fobs and tobs arrays.
        mprior (numpy.ndarray): Initial guess for the model parameters, mprior[0] = v0, mprior[1] = l, mprior[2] = tprime0, mprior[3] = c, mprior[4:] = f0_array.
        num_iterations (int): Number of iterations to perform for the inversion.
        sigma (float): Standard deviation for the data picks, default is 3.
        off_diagonal (bool): Whether to include off-diagonal elements in the prior covariance matrix, default is False.

    Returns:
        numpy.ndarray: The inverted parameters for the function f. Velocity of the aircraft, distance of closest approach, time of closest approach, and the fundamental frequency produced by the aircraft.
        numpy.ndarray: The covariance matrix of the inverted parameters.
        numpy.ndarray: The array of the fundamental frequency produced by the aircraft.
    """

    w = len(mprior[4:]) #number of overtones

    qv = 0
    cprior0 = np.zeros((len(mprior),len(mprior)))

    f0_sigma = sigma_prior[0]
    v0_sigma = sigma_prior[1]
    l_sigma = sigma_prior[2]
    tprime0_sigma = sigma_prior[3]
    c_sigma = sigma_prior[4]

    if off_diagonal:
        cprior0[4:][2] =  -0.4*f0_sigma*tprime0_sigma

        cprior0[0][1] = -0.7*v0_sigma*l_sigma
        cprior0[0][3] = 0.85*v0_sigma*c_sigma

        cprior0[1][0] = -0.7*v0_sigma*l_sigma
        cprior0[1][3] = -0.7*l_sigma*c_sigma

        cprior0[2][4:] =  -0.4*f0_sigma*tprime0_sigma

        cprior0[3][0] = 0.85*v0_sigma*c_sigma
        cprior0[3][1] = -0.7*l_sigma*c_sigma
    
    for row in range(len(cprior0)):
        if row == 0:
            cprior0[row][row] = v0_sigma**2
        elif row == 1:
            cprior0[row][row] = l_sigma**2
        elif row == 2:
            cprior0[row][row] = tprime0_sigma**2
        elif row == 3:
            cprior0[row][row] = c_sigma**2
        else:
            cprior0[row][row] = f0_sigma**2
    cprior = cprior0 * (len(mprior))

    Cd0 = np.zeros((len(fobs), len(fobs)), float)
    np.fill_diagonal(Cd0, sigma**2)

    Cd = Cd0*(len(fobs))
    mnew = np.array(mprior)

    while qv < num_iterations:
        if np.any(np.isnan(mnew)) and qv == 0:
            # Handle the case where mnew contains NaN values
            return mprior, cprior0, cprior, mprior[4:], 'Forward Model'
        elif np.any(np.isnan(mnew)):
            mnew = m
            G = G_hold
            Cpost = la.inv(G.T@la.pinv(Cd)@G + la.inv(cprior))
            Cpost0 = la.inv(G.T@la.pinv(Cd0)@G + la.inv(cprior0))
            return mnew, Cpost0, Cpost, f0_array, S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)
        else:
            m = mnew
            
        v0 = m[0]
        l = m[1]
        tprime0 = m[2]
        c = m[3]
        f0_array = m[4:]

        fpred = []
        G = np.zeros((0,w+4))
        cum = 0
        
        # CORRECTED: Use new timing and derivatives for each overtone
        for p in range(w):
            new_row = np.zeros(w+4)
            f0 = f0_array[p]
            
            for j in range(cum,cum+peaks_assos[p]):
                tprime = tobs[j]

                # CORRECTED frequency calculation
                ft0p = calc_ft([tprime], tprime0, f0, v0, l, c)[0]

                # CORRECTED derivatives
                f_derivef0, f_derivev0, f_derivel, f_derivetprime0, f_derivec = df(f0, v0, l, tprime0, tprime, c)
                
                new_row[0] = f_derivev0
                new_row[1] = f_derivel
                new_row[2] = f_derivetprime0
                new_row[3] = f_derivec
                new_row[4+p] = f_derivef0
                G = np.vstack((G, new_row))
                        
                fpred.append(ft0p)
        
            cum = cum + peaks_assos[p]

        Gm = G
        
        # steepest ascent vector (Eq. 6.307 or 6.312)
        gamma = cprior @ Gm.T @ la.inv(Cd) @ (np.array(fpred) - fobs) + (np.array(m)  - np.array(mprior)) # steepest ascent vector
        #===================================================
        # QUASI-NEWTON ALGORITHM (Eq. 6.319, nu=1)
        # approximate curvature
        H = np.identity(len(mnew)) + cprior @ Gm.T @ la.inv(Cd) @ Gm
        dm = -la.inv(H) @ gamma
        mnew = m + dm

        # Check for unreasonable parameter values
        f0_unreasonable = any(f0_val <= 5 or f0_val > 375 for f0_val in mnew[4:])
        unreasonable = (
            f0_unreasonable or
            mnew[0] <= 0 or mnew[0] > 350 or     # v0
            mnew[0] >= mnew[3] or  # v0 must be less than c
            mnew[1] < 0 or mnew[1] > 1e5 or      # l
            mnew[2] < 10 or mnew[2] > 240 or      # tprime0
            mnew[3] < 200 or mnew[3] > 400       # c
        )

        if unreasonable and qv > 0:
            mnew = m
            G = G_hold
            Cpost = la.inv(G.T@la.pinv(Cd)@G + la.inv(cprior))
            Cpost0 = la.inv(G.T@la.pinv(Cd0)@G + la.inv(cprior0))
            return mnew, Cpost0, Cpost, f0_array, S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)
        elif unreasonable and qv == 0:
            return mprior, cprior0, cprior, mprior[4:], 'Forward Model'
        elif np.any(np.isnan(mnew)):
            return mprior, cprior0, cprior, mprior[4:], 'Forward Model'
        else:
            # Store the current G matrix for potential rollback
            G_hold = G.copy()
            
        f0_array = m[4:]
        qv += 1
        print(f"Iteration {qv}: {mnew}")

    Cpost = la.inv(Gm.T@la.inv(Cd)@Gm + la.inv(cprior))
    Cpost0 = la.inv(Gm.T@la.inv(Cd0)@Gm + la.inv(cprior0))
    F_m = S(fpred, fobs, len(fobs), mnew, mprior, cprior, sigma)
    del G, G_hold, Gm, H, dm, gamma, fpred, m, Cd, Cd0, cprior, cprior0
    return mnew, Cpost0, Cpost, mnew[4:], F_m

def get_auto_picks_full(peaks, time_peaks, times, frequencies, spec, corridor_width, tprime0, v0, l, c, sigma_prior, vmax):
    """
    Get automatic picks for all overtones.
    Uses CORRECTED timing that ensures t=0 at closest approach.

    Args:
        peaks (list): List of peak frequencies.
        time_peaks (list): List of times corresponding to the peaks.
        times (np.ndarray): Array of time values from fft.
        frequencies (np.ndarray): Array of frequency values from fft.
        spec (np.ndarray): Spectrogram data from fft.
        corridor_width (float): Width of the corridor for picking.
        tprime0 (float): Model parameter for the arrival time.
        v0 (float): Model parameter for the velocity.
        l (float): Model parameter for the distance.
        c (float): Model parameter for the speed of sound.
        sigma_prior (float): Prior uncertainty for the model parameters.
        vmax (float): Maximum amplitude value for peak detection.

    Returns:
        list: List of observed times.
        list: List of observed frequencies.
        list: List of counts of peaks associated with each overtone, for indexing.
        list: List of fundamental frequencies calculated for each peak.
    """

    peaks_assos = []
    fobs = []
    tobs = []
    f0_array = []
  
    for pp in range(len(peaks)):
        tprime = time_peaks[pp]
        ft0p = peaks[pp]
        
        # CORRECTED: Use updated calc_f0 function
        f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
        f0_array.append(f0)

        maxfreq = []
        coord_inv = []
        ttt = []

        # CORRECTED: Use updated calc_ft function
        ft = calc_ft(times, tprime0, f0, v0, l, c)

        for t_f in range(len(times)):
            upper = int(ft[t_f] + corridor_width)
            lower = int(ft[t_f] - corridor_width)
            
            # Find closest index to upper and lower in frequencies array
            lower_index = np.argmin(np.abs(frequencies - lower))
            upper_index = np.argmin(np.abs(frequencies - upper))
            
            if lower < 0:
                lower = 0
            elif lower >= 250:
                continue
            else:
                pass
                
            if upper > 250:
                upper = 250

            tt = spec[lower_index:upper_index, t_f]

            max_amplitude_index,_ = find_peaks(tt, prominence = 15, wlen=10, height=vmax*0.1)
            if len(max_amplitude_index) == 0:
                continue

            maxa = np.argmax(tt[max_amplitude_index])

            # Get the corresponding index into tt
            peak_idx = int(max_amplitude_index[maxa])
            freq_index = peak_idx + int(np.round(lower_index,0))
            # Now map it to frequency
            max_amplitude_frequency = frequencies[freq_index] 

            maxfreq.append(max_amplitude_frequency)
            coord_inv.append((times[t_f], max_amplitude_frequency))
            ttt.append(times[t_f])

        if len(ttt) > 0 and f0 <= 230:
            coord_inv_array = np.array(coord_inv)
            mtest = [f0, v0, l, tprime0, c]
            
            # CORRECTED: Use updated invert_f function
            mtest, _, _, _ = invert_f(mtest, sigma_prior, coord_inv_array, num_iterations=2)
            
            # CORRECTED: Use updated calc_ft with refined parameters
            ft = calc_ft(ttt, mtest[3], mtest[0], mtest[1], mtest[2], mtest[4])
            delf = np.array(ft) - np.array(maxfreq)

            count = 0
            for i in range(len(delf)):
                if np.abs(delf[i]) <= (4):
                    fobs.append(maxfreq[i])
                    tobs.append(ttt[i])
                    count += 1
            peaks_assos.append(count)
            
        elif f0 > 230:
            # For high frequencies, accept all picks without refinement
            for i in range(len(ttt)):
                fobs.append(maxfreq[i])
                tobs.append(ttt[i])
            peaks_assos.append(len(maxfreq))
            
        else:
            peaks_assos.append(0)

    return tobs, fobs, peaks_assos, f0_array

jet = ['B737', 'B738', 'B739', 'B733', 'B763', 'B772', 'B77W', 'B788', 'B789', 'B744', 'B748', 'B77L', 'CRJ2', 'B732', 'A332', 'A359', 'E75S']

nrl = NRL()
window = 120  # seconds before the arrival time to load the waveform
rerun_fig = True #Flag rerun the figures without saving the inversion results = True
mk_picks = False

# Loop through each station in text file that we already know comes within 2km of the nodes
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/node_crossings_db_UTM.txt','r')

for li in file_in.readlines():
    text = li.split(',')
    date = text[0]
    month = int(date[4:6])
    day = date[6:8]
    flight_num = text[1]
    closest_time = float(text[5])
    sta = text[9]
    equip = text[10]
    alt = float(text[6]) 
    speed_gt = float(text[7]) 
    dist_m = float(text[4])   # Distance in meters
    elev = get_sta_elevation(sta)
    height_m = alt - elev
    distance_gt = np.sqrt(dist_m**2 + (height_m)**2) 

    folder_spec = equip + '_spec_c'
    folder_spectrum = equip + '_spectrum_c'
    DIR = '/scratch/irseppi/nodal_data/plane_info/inversion_results_test_check/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    if os.path.exists(DIR):
        continue
    file_name = '/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/' + equip + '_data_picks/inversepicks/2019-0' + str(month) + '-' + str(day) + '/' + str(flight_num) + '/' + str(sta) + '/' + str(closest_time) + '_' + str(flight_num) + '.csv'
    if not os.path.exists(file_name):
        continue

    else:
        coords = []
        with open(file_name, 'r') as file:
            for line in file:
                pick_data = line.split(',')
                coords.append((float(pick_data[0]), float(pick_data[1])))
            if len(pick_data) == 4:
                start_time = float(pick_data[2])
            else:
                file.close() 
                continue

        file.close()  
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage load: {mem:.2f} MB")
    coords_array = np.array(coords)
    if len(coords_array) == 0:
        continue

    elif equip == 'C185':
        start_time = start_time - 120

    c = 320 # Default speed of sound, average of dataset, m/s
    fa = np.max(coords_array[:, 1]) 
    fr = np.min(coords_array[:, 1])
    #insert method to get initial model here
    fm = (fa+fr)/2 

    #find the closest coordinate to f0
    closest_index = np.argmin(np.abs(coords_array[:, 1] - fm))
    f0 = coords_array[closest_index, 1] 
    tprime0 = coords_array[closest_index, 0]  
    t_hold = np.inf
    for i,t in enumerate(coords_array[:, 0]):
        if t != tprime0:
            if (t - tprime0) < t_hold:
                t_hold = abs(t - tprime0)
                second_index = i

    v0 = c*abs(fa-fr) / (2 * f0)
    v0 = -c + c*np.sqrt(f0**2 + (fa-fr)**2) / (fa+fr)
    slope = (coords_array[closest_index,1] - coords_array[second_index,1]) / (coords_array[closest_index,0] - coords_array[second_index,0])
    l = -((f0*v0**2/c)*(1-(v0/c)**2)**(-3/2))/slope 
    m0 = [f0, v0, l, tprime0, c]

    data, fs, torg, title = load_waveform(sta, start_time)
    frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant')
    if len(times) == 0 or len(frequencies) == 0 or len(Sxx) == 0:
        continue
    spec, MDF = remove_median(Sxx)
    print('Initial model:', m0)
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage load: {mem:.2f} MB")

    middle_index =  len(times) // 2
    middle_column = spec[:, middle_index]
    vmin = 0  
    vmax = np.max(middle_column) 

        
    m0 = [f0, v0, l, tprime0, c]
    sigma_prior = [40, 1, 1, 200, 1]
    m,_,_, F_m = invert_f(m0,sigma_prior, coords_array, num_iterations=3)
    m0[0] = m[0]
    m0[3] = m[3]

    tf = np.arange(0, 240, 1)

    sigma_f0 = 150
    sigma_v0 = 100
    sigma_l = 10000
    sigma_tprime0 = 200
    sigma_c = 100


    m0 = [f0, v0, l, tprime0, c]
    sigma_prior = [sigma_f0, sigma_v0, sigma_l, sigma_tprime0, sigma_c]
    m,_,_, F_m = invert_f(m0,[sigma_f0, sigma_v0, sigma_l, sigma_tprime0, sigma_c], coords_array, num_iterations=3)
    v0 = m[1]
    l = m[2]
    tprime0 = m[3]
    c = m[4]
    mprior = []
    mprior.append(v0)
    mprior.append(l)
    mprior.append(tprime0)
    mprior.append(c)

    mprior[2] = tprime0
    mprior[3] = c

    output2 = '/home/irseppi/REPOSITORIES/parkshwynodal/input/Data_Picks/' + equip + '_data_picks/overtonepicks/2019-0' + str(month) + '-' + str(day) + '/' + str(flight_num) + '/' + str(sta) + '/' + str(closest_time) + '_' + str(flight_num) + '.csv'
    if not os.path.exists(output2):
        continue
    else:
        peaks = []
        freqpeak = []
        with open(output2, 'r') as file:
            for line in file:
                pick_data = line.split(',')
                peaks.append(float(pick_data[1]))
                freqpeak.append(float(pick_data[0]))
        file.close()  
    if len(peaks) <= 15:
        corridor_width = 10
    else:
        corridor_width = 5
    try:
        tobs, fobs, peaks_assos, f0_array = get_auto_picks_full(peaks,freqpeak, times, frequencies, spec, corridor_width, tprime0, v0, l, c, sigma_prior, vmax)
    except:
        continue

    if len(fobs) == 0:
        continue

    for o in range(len(f0_array)):
        mprior.append(float(f0_array[o]))

    tobs, fobs, peaks_assos = time_picks(month, day, flight_num, sta, equip, tobs, fobs, closest_time, start_time, spec, times, frequencies, vmin, vmax, len(peaks), peaks_assos, make_picks=mk_picks)
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage: {mem:.2f} MB")
    if abs(slope) < 1:
        sigma_prior = [10, 125, 15000, 30, 100]
    else:
        sigma_prior = [10, 30, 500, 30, 100]
    if equip in jet:
        sigma_prior = [100, 300, 50000, 100, 100]

    m, covm0, covm, f0_array, F_m = full_inversion(fobs, tobs, peaks_assos, mprior, sigma_prior, num_iterations=2, sigma=3, off_diagonal=False)
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage: {mem:.2f} MB")
    v0 = m[0]
    l = m[1]
    tprime0 = m[2]
    c = m[3]

    covm = np.sqrt(np.diag(covm))
    covm0 = np.sqrt(np.diag(covm0))
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage pre: {mem:.2f} MB")
    closest_index = np.argmin(np.abs(tprime0 - times))
    arrive_time = spec[:,closest_index]
    for i in range(len(arrive_time)):
        if arrive_time[i] < 0:
            arrive_time[i] = 0
    print(slope)
    print(speed_gt, distance_gt)
    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/inversion_results_test_check/' + folder_spec + '/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    qnum = plot_spectrogram(data, fs, torg, title, spec, times, frequencies, tprime0, v0, l, c, f0_array, F_m, arrive_time, MDF, covm0, flight_num, middle_index,mprior[2], closest_time, BASE_DIR, plot_show=False, gt = False)
    qnum = "__"
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage spec 1: {mem:.2f} MB")
    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/inversion_results_test_check/' + folder_spectrum + '/20190'+str(month)+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'
    make_base_dir(BASE_DIR)
    plot_spectrum(spec, frequencies, tprime0, v0, l, c, f0_array, arrive_time, fs, closest_index, closest_time, sta, BASE_DIR)
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage spec 2: {mem:.2f} MB")
    if rerun_fig == False:
        output = open('output/inv_results_no_g_truth_test_check/' + equip + '_full_inv_results.csv', 'a')
        output.write(str(date)+','+str(flight_num)+','+str(sta)+','+str(closest_time)+','+str(v0)+','+str(l)+','+str(tprime0)+','+ str(start_time + tprime0) + ','+str(c)+','+str(f0_array)+','+str(covm0)+','+str(qnum)+','+str(c)+','+str(F_m)+',\n') 
        output.close()
    else:
        continue  # Skip saving results if rerun_fig is True
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage post: {mem:.2f} MB")
    # Explicitly delete large variables and collect garbage to free memory
    # Delete all variables and objects that may impact short-term memory
    del data, fs, torg, title
    del frequencies, times, Sxx, spec, MDF
    del coords, coords_array
    del m, covm0, covm, f0_array, F_m, arrive_time, BASE_DIR
    del peaks, freqpeak, tobs, fobs, peaks_assos, mprior
    del date, month, day, flight_num, closest_time, sta, equip
    del alt, speed_gt, dist_m, elev, height_m, distance_gt
    del folder_spec, folder_spectrum, DIR, file_name
    del start_time, c, fa, fr, fm, closest_index, f0, tprime0, t_hold, second_index
    del v0, slope, l, m0, sigma_prior, tf
    del sigma_f0, sigma_v0, sigma_l, sigma_tprime0, sigma_c
    del output2, corridor_width, qnum

    gc.collect()
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2) 
    print(f"Memory usage: {mem:.2f} MB")