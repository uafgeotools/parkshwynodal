import numpy as np
from scipy.signal import find_peaks
import numpy.linalg as la
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from src.main_inv_fig_functions import remove_median
from src.doppler_funcs import S
from matplotlib.ticker import MaxNLocator
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime

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



# Download waveform data from IRIS PH5WS
client = Client("http://service.iris.edu", service_mappings={"dataselect": "http://service.iris.edu/ph5ws/dataselect/1"})
starttime = UTCDateTime("2019-03-04T01:17:22")
endtime = UTCDateTime("2019-03-04T01:21:22")
st = client.get_waveforms("ZE", "1010", "*", "DPZ", starttime, endtime)
tr = st[0]
data = tr.data
torg = tr.times()
fs = int(tr.stats.sampling_rate)
title = f'{tr.stats.network}.{tr.stats.station}.{tr.stats.location}.{tr.stats.channel} − starting {tr.stats["starttime"]}'

# Compute spectrogram
WIN_LEN = 1  # window length, in s
NPER = int(WIN_LEN * fs)
frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=NPER, noverlap=int(NPER * .9), detrend='constant')

spec, MDF = remove_median(Sxx)  # Remove median for better visualization
middle_index = len(times) // 2
middle_column = spec[:, middle_index]
vmin, vmax = 0, np.max(middle_column)
coords = []

coords.append((30.375000000000007, 137.44588744588745))
coords.append((51.094758064516135, 136.0930735930736))
coords.append((78.56048387096774, 136.0930735930736))
coords.append((99.28024193548387, 130.6818181818182))
coords.append((112.77217741935485, 115.12445887445887))
coords.append((121.92741935483872, 109.71320346320346))
coords.append((134.45564516129033, 102.94913419913419))
coords.append((151.32056451612905, 101.59632034632034))
coords.append((169.6310483870968, 100.24350649350649))
coords.append((188.4233870967742, 100.24350649350649))
coords_array = np.array(coords)

peaks = []
freqpeak = []

freqpeak.append(112.77217741935485)
peaks.append(115.80086580086581)
freqpeak.append(113.25403225806451)
peaks.append(136.7694805194805)
freqpeak.append(113.25403225806451)
peaks.append(154.35606060606062)
freqpeak.append(113.25403225806451)
peaks.append(174.6482683982684)
freqpeak.append(114.21774193548387)
peaks.append(190.2056277056277)
freqpeak.append(116.62701612903226)
peaks.append(226.05519480519484)
freqpeak.append(113.25403225806451)
peaks.append(57.62987012987013)
freqpeak.append(112.29032258064518)
peaks.append(96.18506493506493)

# Estimate initial model parameters from picked points
c = 320  # Speed of sound (m/s)
fa, fr = np.max(coords_array[:, 1]), np.min(coords_array[:, 1])  # Max/min frequency
fm = (fa + fr) / 2
closest_index = np.argmin(np.abs(coords_array[:, 1] - fm))
f0, tprime0 = coords_array[closest_index, 1], coords_array[closest_index, 0]
t_hold, second_index = np.inf, None
for i, t in enumerate(coords_array[:, 0]):
    if t != tprime0 and abs(t - tprime0) < t_hold:
        t_hold = abs(t - tprime0)
        second_index = i
v0 = c * abs(fa - fr) / (2 * f0)  # Initial velocity estimate
slope = (coords_array[closest_index, 1] - coords_array[second_index, 1]) / (coords_array[closest_index, 0] - coords_array[second_index, 0])
l = -((f0 * v0 ** 2 / c) * (1 - (v0 / c) ** 2) ** (-3 / 2)) / slope  # Initial length estimate
#l = tprime0 *c
m0 = [f0, v0, l, tprime0, c]
sigma_prior = [40, 1, 1, 200, 1]  # Initial prior uncertainties

# First inversion to refine model - USING CORRECTED TIMING
m, _, _, F_m = invert_f(m0, sigma_prior, coords_array, num_iterations=3)
m0[0], m0[3] = m[0], m[3]

# Second inversion with wider priors
sigma_prior = [150, 100, 10000, 200, 100]
m, _, _, F_m = invert_f(m0, sigma_prior, coords_array, num_iterations=3)
v0, l, tprime0, c = m[1], m[2], m[3], m[4]
mprior = [v0, l, tprime0, c]

# Automatically associate picked peaks with overtone curves - USING CORRECTED TIMING
corridor_width = 10 if len(peaks) <= 15 else 5
tobs, fobs, peaks_assos, f0_array = get_auto_picks_full(peaks, freqpeak, times, frequencies, spec, corridor_width, tprime0, v0, l, c, sigma_prior, vmax)
mprior += [float(f) for f in f0_array]

start_time = 23.14717741935484
end_time = 222.15322580645164

# Filter picks to only those within the selected time window
ftobs, ffobs, peak_ass = [], [], []
cum = 0
for p in range(len(f0_array)):
    count = 0
    for j in range(cum, cum + peaks_assos[p]):
        if start_time <= tobs[j] <= end_time:
            ftobs.append(tobs[j])
            ffobs.append(fobs[j])
            count += 1
    cum += peaks_assos[p]
    peak_ass.append(count)
peaks_assos = peak_ass
tobs, fobs = ftobs, ffobs

# Final inversion using filtered picks - USING CORRECTED TIMING
sigma_prior = [10, 125, 15000, 30, 100] if abs(slope) < 1 else [10, 30, 500, 30, 100]
m, covm0, covm, f0_array, F_m = full_inversion(fobs, tobs, peaks_assos, mprior, sigma_prior, num_iterations=2, sigma=3, off_diagonal=False)
v0, l, tprime0, c = m[0], m[1], m[2], m[3]
Cpost, Cpost0 = np.sqrt(np.diag(covm)), np.sqrt(np.diag(covm0))

# Plot results
closest_index = np.argmin(np.abs(tprime0 - times))
arrive_time = np.clip(spec[:, closest_index], 0, None)
vmin, vmax = np.min(arrive_time), np.max(arrive_time)
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=False, figsize=(8, 10))

# Plot raw waveform
ax1.plot(torg, data, 'k', linewidth=0.5)
ax1.set_title(title)
ax1.margins(x=0)
ax1.set_position([0.125, 0.6, 0.775, 0.3])
ax1.set_ylabel('Counts')

# Plot spectrogram and inversion results
cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
ax2.set_xlabel('Time (s)')
ax2.axvline(x=tprime0, c='#377eb8', ls='--', linewidth=0.7, label=f"t\u2080' = {tprime0:.2f} s")
f0lab = sorted(f0_array)

# CORRECTED: Use updated calc_ft for plotting
for f0 in f0lab:
    ft = calc_ft(times, tprime0, f0, v0, l, c)  # Uses corrected timing
    ax2.plot(times, ft, '#377eb8', ls=(0, (5, 20)), linewidth=0.7)
    
    ax2.scatter(tprime0, f0, color='black', marker='x', s=30)
    
fss = 'x-small'

# Estimate overtone frequency spacing and uncertainty
if len(f0lab) > 1:
    f_range = []
    NTRY = 1000
    for _ in range(NTRY):
        ftry = [f0_array[i - 4] + np.random.uniform(-Cpost0[i], Cpost0[i]) for i in range(4, len(Cpost0))]
        ftry = np.sort(ftry)
        f1 = [ftry[g] - ftry[g - 1] for g in range(1, len(ftry))]
        f_range.append(np.nanmedian(f1))
    med_df = np.nanmedian(f_range)
    mad_df = np.nanmedian(np.abs(f_range - med_df))
else:
    med_df = mad_df = "NaN"

# Format overtone frequencies for display
if len(f0lab) > 10:
    f0lab_lines = [', '.join([f"{f:.2f}" for f in f0lab[i:i + 10]]) for i in range(0, len(f0lab), 10)]
    f0lab_str = '[%s]' % (',\n'.join(f0lab_lines))
else:
    f0lab_str = '[' + ', '.join([f"{f:.2f}" for f in f0lab]) + ']'

# Compose plot title with inversion results and uncertainties
if isinstance(F_m, str):
    misfit_str = f"\n[{F_m}]"
else:
    misfit_str = f"\nMisfit: {F_m:.4f}"
df_str = f", df\u2080 = {med_df:.2f} \u00B1 {mad_df:.2f} Hz" if med_df != "NaN" else ""
ax2.set_title(
    f"t\u2080'= {tprime0:.2f} \u00B1 {Cpost0[2]:.2f} s, v\u2080 = {v0:.2f} \u00B1 {Cpost0[0]:.2f} m/s, "
    f"c = {c:.2f} \u00B1 {Cpost0[3]:.2f} m/s, l = {l:.2f} \u00B1 {Cpost0[1]:.2f} m, \n"
    f"f\u2080 = {f0lab_str} \u00B1 {np.median(Cpost0[3:]):.2f} Hz{df_str}{misfit_str}",
    fontsize=fss
)
ax2.legend(loc='upper right', fontsize='small')
ax2.set_ylabel('Frequency (Hz)')
ax2.margins(x=0)


ax2.set_ylim(0, int(fs / 2))
ax1.tick_params(axis='both', which='major', labelsize=9)
ax2.tick_params(axis='both', which='major', labelsize=9)
ax3.tick_params(axis='both', which='major', labelsize=9)


vmax_freq = np.max(arrive_time)
ax3.grid()

ax3.plot(frequencies, spec[:, closest_index], c='#377eb8')

for pp in range(len(f0_array)):
    f0 = f0_array[pp]
    if fs / 2 < f0:
        continue
    tprime = tprime0
    t = ((tprime - tprime0) - np.sqrt((tprime - tprime0) ** 2 - (1 - v0 ** 2 / c ** 2) * ((tprime - tprime0) ** 2 - l ** 2 / c ** 2))) / (1 - v0 ** 2 / c ** 2)
    ft0p = f0 / (1 + (v0 / c) * (v0 * t) / (np.sqrt(l ** 2 + (v0 * t) ** 2)))
    if np.isnan(ft0p):
        continue
    if ft0p > 250:
        continue

    upper = int(ft0p + 10)
    lower = int(ft0p - 10)
    tt = spec[lower:upper, closest_index]
    if upper > 250 or lower < 0:
        freqp = ft0p
        ampp = np.interp(ft0p, frequencies, arrive_time)
    else:
        ampp = np.max(tt)
        freqp = np.argmax(tt) + lower
    ax3.scatter(freqp, ampp, color='black', marker='x', s=100, zorder=10)
    ax3.text(freqp - 1, ampp + 0.8, f"{freqp:.2f}", fontsize=12, fontweight='bold')

ax3.set_xlim(0, int(fs / 2))
ax3.set_xticks(np.arange(0, int(fs / 2) + 1, 20))

ax3.set_ylim(0, vmax_freq * 1.1)

plt.show()
