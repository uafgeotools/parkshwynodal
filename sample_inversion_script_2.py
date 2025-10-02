import numpy as np
from scipy.signal import find_peaks
import numpy.linalg as la
import obspy
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from src.main_inv_fig_functions import remove_median
from src.doppler_funcs import S
from matplotlib.ticker import MaxNLocator
def calc_ft(times, tprime0, f0, v0, l, c):
    """
    Calculate the frequency at each given time using the model parameters.

    Args:
        times (list): List of time values.
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
    beta = v0 / c
    beta2 = beta * beta
    
    for tprime in times:
        delta = tprime - tprime0
        
        # CORRECTED: t=0 at closest approach
        numerator = delta + l/c - np.sqrt(beta2 * delta**2 + beta2 * (2*l/c) * delta + (l/c)**2)
        denominator = 1 - beta2
        t = numerator / denominator
        
        # Geometry now makes sense: v0*t is displacement from closest approach
        ft0p = f0 / (1 + (v0/c) * (v0*t) / np.sqrt(l**2 + (v0*t)**2))
                                
        ft.append(ft0p)
    return ft
def df(f0, v0, l, tp0, tp, c):   
    """
    Calculate the derivatives of f with respect to f0, v0, l, tprime0 and c.
    Uses CORRECTED timing that ensures t=0 at closest approach.
    """
    
    beta = v0 / c
    beta2 = beta * beta
    delta = tp - tp0
    lc = l / c
    
    # Helper terms
    A = beta2 * delta**2 + beta2 * 2 * lc * delta + lc**2
    sqrtA = np.sqrt(A)
    
    # CORRECTED emitter time
    t = (delta + lc - sqrtA) / (1 - beta2)
    
    # Geometry terms
    v0t = v0 * t
    R = np.sqrt(l**2 + v0t**2)
    cos_theta = v0t / R
    
    # The frequency function
    F = f0 / (1 + beta * cos_theta)
    
    # ===== PARTIAL DERIVATIVES =====
    
    # 1. Derivative with respect to f0
    f_derivef0 = 1 / (1 + beta * cos_theta)
    
    # 2. Derivative with respect to v0
    # Need ∂t/∂v0, ∂cos_theta/∂v0
    dAdbeta = 2 * beta * delta**2 + 2 * beta * 2 * lc * delta
    dbeta_dv0 = 1 / c
    
    dsqrtA_dv0 = (0.5 / sqrtA) * dAdbeta * dbeta_dv0
    dt_dv0 = (-dsqrtA_dv0 - 2 * beta * t * (1 - beta2) / (1 - beta2)**2)
    
    dcos_theta_dv0 = (t + v0 * dt_dv0) / R - (v0t * (v0 * dt_dv0 + v0t * v0 / R**2)) / R**2
    dF_dv0 = -f0 * (dbeta_dv0 * cos_theta + beta * dcos_theta_dv0) / (1 + beta * cos_theta)**2
    
    f_derivev0 = dF_dv0
    
    # 3. Derivative with respect to l
    dA_dl = beta2 * 2 * delta / c + 2 * lc / c
    dsqrtA_dl = (0.5 / sqrtA) * dA_dl
    dt_dl = (1/c - dsqrtA_dl) / (1 - beta2)
    
    dcos_theta_dl = (v0 * dt_dl) / R - (v0t * (l + v0t * v0 * dt_dl / R)) / R**2
    dF_dl = -f0 * beta * dcos_theta_dl / (1 + beta * cos_theta)**2
    
    f_derivel = dF_dl
    
    # 4. Derivative with respect to tp0
    dA_dtp0 = -beta2 * 2 * delta - beta2 * 2 * lc
    dsqrtA_dtp0 = (0.5 / sqrtA) * dA_dtp0
    dt_dtp0 = (-1 - dsqrtA_dtp0) / (1 - beta2)
    
    dcos_theta_dtp0 = (v0 * dt_dtp0) / R - (v0t * (v0t * v0 * dt_dtp0 / R)) / R**2
    dF_dtp0 = -f0 * beta * dcos_theta_dtp0 / (1 + beta * cos_theta)**2
    
    f_derivetprime0 = dF_dtp0
    
    # 5. Derivative with respect to c
    dbeta_dc = -v0 / c**2
    dA_dc = (2 * beta * dbeta_dc * delta**2 + 
             2 * beta * dbeta_dc * 2 * lc * delta +
             beta2 * 2 * delta * (-l/c**2) +
             -2 * lc * l / c**2)
    dsqrtA_dc = (0.5 / sqrtA) * dA_dc
    dt_dc = (-l/c**2 - dsqrtA_dc + 2 * beta * dbeta_dc * t * (1 - beta2)) / (1 - beta2)**2
    
    dcos_theta_dc = (v0 * dt_dc) / R - (v0t * (v0t * v0 * dt_dc / R)) / R**2
    dF_dc = -f0 * (dbeta_dc * cos_theta + beta * dcos_theta_dc) / (1 + beta * cos_theta)**2
    
    f_derivec = dF_dc
    
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
            
            # CORRECTED timing calculation
            beta = v0 / c
            beta2 = beta * beta
            delta = tprime - tprime0
            lc = l / c
            
            A = beta2 * delta**2 + beta2 * 2 * lc * delta + lc**2
            sqrtA = np.sqrt(A)
            t = (delta + lc - sqrtA) / (1 - beta2)
            
            # CORRECTED frequency calculation
            ft0p = f0 / (1 + (v0/c) * (v0*t) / np.sqrt(l**2 + (v0*t)**2))
            
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
                
                # CORRECTED timing calculation
                beta = v0 / c
                beta2 = beta * beta
                delta = tprime - tprime0
                lc = l / c
                
                A = beta2 * delta**2 + beta2 * 2 * lc * delta + lc**2
                sqrtA = np.sqrt(A)
                t = (delta + lc - sqrtA) / (1 - beta2)
                
                # CORRECTED frequency calculation
                ft0p = f0 / (1 + (v0/c) * (v0*t) / np.sqrt(l**2 + (v0*t)**2))
                
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
    beta = v0 / c
    beta2 = beta * beta
    delta = tprime - tprime0
    
    # CORRECTED: t=0 at closest approach
    numerator = delta + l/c - np.sqrt(beta2 * delta**2 + beta2 * (2*l/c) * delta + (l/c)**2)
    denominator = 1 - beta2
    t = numerator / denominator
    
    # Geometry now makes sense
    f0 = ft0p * (1 + (v0/c) * (v0*t) / np.sqrt(l**2 + (v0*t)**2))

    return f0
# Interactive picking of points on spectrogram for overtone curve
def pick_points_on_spectrogram(times, frequencies, spec, vmin, vmax, prompt, axvline=None):
    coords = []
    plt.figure()
    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
    if axvline is not None:
        plt.axvline(x=axvline, c='#377eb8', ls='--')
    def onclick(event):
        if event.xdata is not None and event.ydata is not None:
            coords.append((event.xdata, event.ydata))
            plt.scatter(event.xdata, event.ydata, color='black', marker='x')
            plt.draw()
            print('Clicked:', event.xdata, event.ydata)
    plt.gcf().canvas.mpl_connect('button_press_event', onclick)
    plt.show(block=True)
    return coords

# Interactive picking of single points (overtone peaks)
def pick_single_points(times, frequencies, spec, vmin, vmax, prompt, axvline=None):
    peaks, freqpeak = [], []
    plt.figure()
    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
    if axvline is not None:
        plt.axvline(x=axvline, c='#377eb8', ls='--')
    def onclick(event):
        if event.xdata is not None and event.ydata is not None:
            peaks.append(event.ydata)
            freqpeak.append(event.xdata)
            plt.scatter(event.xdata, event.ydata, color='black', marker='x')
            plt.draw()
            print('Clicked:', event.xdata, event.ydata)
    plt.gcf().canvas.mpl_connect('button_press_event', onclick)
    plt.show(block=True)
    return peaks, freqpeak

# Interactive picking of time window for inversion
def pick_time_window(times, frequencies, spec, vmin, vmax, tobs, fobs):
    set_time = []
    plt.figure()
    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
    plt.scatter(tobs, fobs, color='black', marker='x')
    def onclick(event):
        if event.xdata is not None:
            set_time.append(event.xdata)
            plt.scatter(event.xdata, event.ydata, color='red', marker='x')
            plt.draw()
            print('Clicked:', event.xdata, event.ydata)
    plt.gcf().canvas.mpl_connect('button_press_event', onclick)
    plt.show(block=True)
    return set_time

waveform =  '/home/irseppi/Downloads/IM.IS02.01.CDF.2025.228.20'
tr = obspy.read(waveform)[0]
#tr.trim(tr.stats.starttime + 650, tr.stats.starttime + 900)
tr.trim(tr.stats.starttime + 420, tr.stats.starttime + 540)
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

# User picks overtone curve points
print("Please pick the points on the spectrogram that correspond to the primary overtone of the doppler curves.")
while True:
    coords = pick_points_on_spectrogram(times, frequencies, spec, vmin, vmax, "Pick overtone curve points")
    if input("Do you want to repick your points? (y or n)").lower() != 'y':
        break
coords_array = np.array(coords)

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

# User picks overtone peaks
print("Please pick one point on each overtone, it does not have to be at the center of the doppler.")
while True:
    peaks, freqpeak = pick_single_points(times, frequencies, spec, vmin, vmax, "Pick overtone peaks", axvline=tprime0)
    if input("Do you want to repick your points? (y or n)").lower() != 'y':
        break

# Automatically associate picked peaks with overtone curves - USING CORRECTED TIMING
corridor_width = 10 if len(peaks) <= 15 else 5
tobs, fobs, peaks_assos, f0_array = get_auto_picks_full(peaks, freqpeak, times, frequencies, spec, corridor_width, tprime0, v0, l, c, sigma_prior, vmax)
mprior += [float(f) for f in f0_array]

# User picks time window for inversion
print('Please pick two points on the spectrogram that correspond to the start and end of the time window you want pull data from in the inversion.')
while True:
    set_time = pick_time_window(times, frequencies, spec, vmin, vmax, tobs, fobs)
    if input("Do you want to repick your points? (y or n)").lower() != 'y':
        break
start_time, end_time = set_time[:2]

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
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8, 6))

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
    
    # CORRECTED: Calculate frequency at tprime0 using proper timing
    beta = v0 / c
    beta2 = beta * beta
    delta = tprime0 - tprime0  # = 0
    lc = l / c
    
    A = beta2 * delta**2 + beta2 * 2 * lc * delta + lc**2
    sqrtA = np.sqrt(A)
    t = (delta + lc - sqrtA) / (1 - beta2)  # This should be 0 at closest approach
    
    # CORRECTED: Frequency calculation
    ft0p = f0 / (1 + (v0/c) * (v0*t) / np.sqrt(l**2 + (v0*t)**2))
    ax2.scatter(tprime0, ft0p, color='black', marker='x', s=30)
    
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

# Add colorbar for spectrogram
ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])
cbar = plt.colorbar(mappable=cax, cax=ax3)
cbar.locator = MaxNLocator(integer=True)
cbar.update_ticks()
ax3.set_ylabel('Relative Amplitude (dB)')
ax2.margins(x=0)
#ax2.set_xlim(0, 240)
ax2.set_ylim(0, int(fs / 2))
ax1.tick_params(axis='both', which='major', labelsize=9)
ax2.tick_params(axis='both', which='major', labelsize=9)
ax3.tick_params(axis='both', which='major', labelsize=9)

# Overlay median-detrended frequency (MDF) for reference
spec2 = 10 * np.log10(MDF)
middle_column2 = spec2[:, middle_index]
vmin2, vmax2 = np.min(middle_column2), np.max(middle_column2)
ax4 = fig.add_axes([0.125, 0.11, 0.07, 0.35], sharey=ax2)
ax4.plot(middle_column2, frequencies, c='#ff7f00')
ax4.set_ylim(0, int(fs / 2))
ax4.set_xlim(vmax2 * 1.1, vmin2)
ax4.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
ax4.grid(axis='y')
plt.show()
plt.close()