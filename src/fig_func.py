import gc
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from scipy.signal import find_peaks
from matplotlib.ticker import MaxNLocator
from src.get_save_data import make_base_dir
from src.doppler_funcs import DopplerCalc, DopplerInversion


############################################################################


class SpecPlot:
    def __init__(self, Sxx, sampling_rate, t_wf, data, times, frequencies, 
                 tarrive, spec_window=120):
        
        self.Sxx = Sxx
        self.sampling_rate = sampling_rate
        self.t_wf = t_wf
        self.data = data
        self.times = times
        self.frequencies = frequencies
        self.tarrive = tarrive
        self.start_time = tarrive - spec_window

    ############################################################################

    def remove_median(self):
        """
        Remove the median from the spectrogram.

        Args:
            Sxx (array): The spectrogram data.

        Returns:
            spec: The spectrogram data with the median removed
            MDF: The median removed from the spectrogram
        """

        a, b = self.Sxx.shape

        MDF = np.zeros((a,b))
        for row in range(len(self.Sxx)):
            median = np.median(self.Sxx[row])
            MDF[row, :] = median

        # Avoid log10(0) by replacing zeros with a small positive value
        Sxx_safe = np.where(self.Sxx == 0, 1e-10, self.Sxx)
        MDF_safe = np.where(MDF == 0, 1e-10, MDF)

        spec = 10 * np.log10(Sxx_safe) - (10 * np.log10(MDF_safe))
        self.spec = spec
        self.MDF = MDF
        return spec, MDF

    ############################################################################

    def plot_spectrogram(self, title, t0, v, d0, c, fs_array,
                        F_m, Cpost0, middle_index, file_name=None, 
                        plot_show=True, gt = True):
        """
        Plot and save the waveform, unfiltered, and the spectrogram of the given 
        data. Include the estimated curve using the final model parameters 
        outputs from the inversions.

        Args:
            data (np.ndarray): The waveform data.
            sampling_rate (int): The sampling frequency.
            t_wf (np.ndarray): The time array for the waveform.
            title (str): The title of the plot.
            spec (np.ndarray): The spectrogram data (2D array).
            times (np.ndarray): The time array for the spectrogram.
            frequencies (np.ndarray): The frequency array for the spectrogram.
            t0 (float): The estimated time of aircraft closest approach to the 
                station.
            v (float): The velocity.
            d0 (float): The distance.
            c (float): The speed of sound.
            fs_array (np.ndarray): The array of frequencies.
            F_m (float or str): The data misfit value.
            MDF (np.ndarray): Median removed from spectrogram (2D array).
            Cpost0 (np.ndarray): The normalized posterior covariance matrix.
            middle_index (int): The index of the middle column.
            file_name (str, optional): The name of the file to save the plot. 
                If None, the plot will not be saved. Defaults to None.
            plot_show (bool): If True, show the plot and ask user to provide a 
                quality number. If False, save the plot without showing it. 
            gt (bool): If True, the ground truth is used for the initial model 
                in the inversion.

        """
        print('Plotting spectrogram...')
        t0prime = t0 + d0/c
        if gt:
            type_inv = "[FH/GT]"
        else:
            type_inv = "[FH/NGT]"
        closest_index = np.argmin(np.abs(self.times - t0))
        closest_index = np.argmin(np.abs(t0 - self.times))
        arrive_time = self.spec[:,closest_index]
        for i in range(len(arrive_time)):
            if arrive_time[i] < 0:
                arrive_time[i] = 0
        # Plot settings and calculations
        vmin = np.min(arrive_time) 
        vmax = np.max(arrive_time)

        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     
        ax1.plot(self.t_wf, self.data, 'k', linewidth=0.5)
        ax1.set_title(title)

        ax1.margins(x=0)
        ax1.set_position([0.125, 0.6, 0.775, 0.3]) 
        ax1.set_ylabel('Counts')

        # Plot spectrogram
        cax = ax2.pcolormesh(
            self.times, self.frequencies, self.spec, shading='gouraud', 
            cmap='pink_r', vmin=vmin, vmax=vmax
            )		
        ax2.set_xlabel('Time (s)')

        ax2.axvline(
            x=t0prime, c = '#377eb8', ls = '--', linewidth=0.5, 
            label= "t\u2080' = " + "%.2f" % t0prime +' s')
        ax2.axvline(
            x=t0, c = '#e41a1c', ls = '--', linewidth=0.7,
            label= "t\u2080 = " + "%.2f" % t0 +' s')

        for pp in range(len(fs_array)):
            fs = fs_array[pp]
            doppler_calc = DopplerCalc(v, d0, t0, c)
            ft = doppler_calc.calc_ft(self.times, fs)

            ax2.plot(
                self.times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7
                )
            ax2.scatter(t0prime, fs, color='black', marker='x', s=30, zorder=10)

        text_size = 'x-small'
        fslab = sorted(fs_array)

        if len(fs_array) <= 1:
            med_df = "NaN"
            mad_df = "NaN"
        else:
            #Generate random samples of fs values withing their sigma from the  
            #covariance matrix. Calculate the median of the differences and MAD 
            # to obtain error
            f_range = []
            NTRY = 1000
            for N in range(NTRY):
                ftry = []
                for c_index  in range(4, len(Cpost0)):
                    xmin = fs_array[c_index-4] - Cpost0[c_index]
                    xmax = fs_array[c_index-4] + Cpost0[c_index]
                    xtry = xmin + (xmax-xmin)*np.random.rand()
                    ftry.append(xtry)

                ftry = np.sort(ftry)
                f1 = []
                for g in range(len(ftry)):
                    if g == 0:
                        continue
                    diff = ftry[g] - ftry[g - 1]
                    f1.append(diff)
                med = np.nanmedian(f1)
                f_range.append(med)
            med_df = np.nanmedian(f_range)
            mad_df = np.nanmedian(np.abs(f_range - med_df))

        if len(fslab) > 10:
            # Split fslab into lines of 10 entries each
            fslab_lines = []
            for i in range(0, len(fslab), 10):
                line = ', '.join(["%.2f" % f for f in fslab[i:i+10]])
                fslab_lines.append(line)
            fslab_str = (',\n').join(fslab_lines)
            fslab_str = '[' + fslab_str + ']'
        else:
            fslab_str = '[' + ', '.join(["%.2f" % f for f in fslab]) + ']'

        if isinstance(F_m, str):
            if med_df == "NaN":
                ax2.set_title(
                    "t\u2080 = "+ "%.2f" % t0 + ' \u00B1 ' + "%.2f" % Cpost0[2] 
                    + ' s, v = ' + "%.2f" % v +' \u00B1 ' + "%.2f" % Cpost0[0]
                    + ' m/s, c = ' + "%.2f" % c +' \u00B1 ' + "%.2f" % Cpost0[3] 
                    + ' m/s, d\u2080 = '+ "%.2f" % d0 + ' \u00B1 ' 
                    + "%.2f" % Cpost0[1] + ' m, \n' + 'f\u209B = ' + fslab_str 
                    + ' \u00B1 ' + "%.2f" % np.median(Cpost0[4:]) + ' Hz'
                    + '\n[' + F_m + ']' + ' ' + type_inv, fontsize=text_size)
            else:
                ax2.set_title(
                    "t\u2080 = " + "%.2f" % t0 + ' \u00B1 ' + "%.2f" % Cpost0[2] 
                    + ' s, v = ' + "%.2f" % v +' \u00B1 ' + "%.2f" % Cpost0[0] 
                    + ' m/s, c = ' + "%.2f" % c +' \u00B1 ' + "%.2f" % Cpost0[3] 
                    + ' m/s, d\u2080 = '+ "%.2f" % d0 +' \u00B1 ' 
                    + "%.2f" % Cpost0[1] + ' m, \n' + 'f\u209B = ' 
                    + fslab_str + ' \u00B1 ' + "%.2f" % np.median(Cpost0[4:]) 
                    + ' Hz, df\u209B = ' + "%.2f" % med_df + ' \u00B1 ' 
                    + "%.2f" % mad_df + ' Hz\n[' + F_m + ']' + ' ' 
                    + type_inv, fontsize=text_size)
                
        elif med_df == "NaN":
            ax2.set_title(
                "t\u2080 = "+ "%.2f" % t0 + ' \u00B1 ' + "%.2f" % Cpost0[2] 
                + ' s, v = ' + "%.2f" % v +' \u00B1 ' + "%.2f" % Cpost0[0] 
                + ' m/s, c = ' + "%.2f" % c +' \u00B1 ' + "%.2f" % Cpost0[3] 
                + ' m/s, d\u2080 = '+ "%.2f" % d0 + ' \u00B1 ' 
                + "%.2f" % Cpost0[1] + ' m, \n' + 'f\u209B = ' + fslab_str 
                + ' \u00B1 ' + "%.2f" % np.median(Cpost0[4:]) + ' Hz\nMisfit: ' 
                + "%.4f" % F_m + ' ' + type_inv, fontsize=text_size)
        else:
            ax2.set_title(
                "t\u2080 = " + "%.2f" % t0 + ' \u00B1 ' + "%.2f" % Cpost0[2] 
                + ' s, v = ' + "%.2f" % v +' \u00B1 ' + "%.2f" % Cpost0[0] 
                + ' m/s, c = ' + "%.2f" % c +' \u00B1 ' + "%.2f" % Cpost0[3] 
                + ' m/s, d\u2080 = '+ "%.2f" % d0 +' \u00B1 ' 
                + "%.2f" % Cpost0[1] + ' m, \n' + 'f\u209B = ' 
                + fslab_str + ' \u00B1 ' + "%.2f" % np.median(Cpost0[4:]) 
                + ' Hz, df\u209B = ' + "%.2f" % med_df + ' \u00B1 ' 
                + "%.2f" % mad_df + ' Hz\nMisfit: ' 
                + "%.4f" % F_m + ' ' + type_inv, fontsize=text_size)

        ax2.legend(loc='upper right',fontsize = 'small')
        ax2.set_ylabel('Frequency (Hz)')
        ax1.set_xlim(0,max(self.t_wf))
        ax2.set_xlim(0,max(self.t_wf))
        ax2.margins(x=0)
        ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

        # Set colorbar with integer ticks only
        cbar = plt.colorbar(mappable=cax, cax=ax3)
        cbar.locator = MaxNLocator(integer=True)
        cbar.update_ticks()
        ax3.set_ylabel('Relative Amplitude (dB)')

        ax2.margins(x=0)
        ax2.set_ylim(0, int(self.sampling_rate/2))

        ax1.tick_params(axis='both', which='major', labelsize=9)
        ax2.tick_params(axis='both', which='major', labelsize=9)
        ax3.tick_params(axis='both', which='major', labelsize=9)
        cbar.ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        cbar.update_ticks()

        # Plot overlay
        spec2 = 10 * np.log10(self.MDF)
        middle_column2 = spec2[:, middle_index]
        vmin2 = np.min(middle_column2)
        vmax2 = np.max(middle_column2)

        # Create ax4 and plot on the same y-axis as ax2
        ax4 = fig.add_axes([0.125, 0.11, 0.07, 0.35], sharey=ax2) 
        ax4.plot(middle_column2, self.frequencies, c='#ff7f00')  
        ax4.set_ylim(0, int(self.sampling_rate/2))
        ax4.set_xlim(vmax2*1.1, vmin2) 
        ax4.tick_params(
            left=False, right=False, labelleft=False, labelbottom=False, 
            bottom=False)
        
        ax4.grid(axis='y')

        if plot_show:
            plt.show()     

        if file_name is not None:
            print('Saving figure to:', file_name)
            fig.savefig(file_name, dpi=600)
        plt.close(fig)
        gc.collect()


################################################################################


class GetPicks(SpecPlot):

    def __init__(self, spec_plot, vmin=None, vmax=None, make_picks=True,
                save_picks=False):
        
        self.__dict__.update(spec_plot.__dict__)  # copy all attributes
        self.vmin = vmin
        self.vmax = vmax
        self.make_picks = make_picks
        self.save_picks = save_picks

    ############################################################################

    def doppler_points(self):
        '''Interactive picking of points on spectrogram for overtone curve'''

        # User picks overtone curve points
        print(
            "Please pick the points on the spectrogram that correspond to the "
            "primary overtone of the doppler curves."
        )
        while True:
            coords = []
            plt.figure()
            plt.pcolormesh(
                self.times, self.frequencies, self.spec, shading='gouraud', 
                cmap='pink_r', vmin=self.vmin, vmax=self.vmax
            )
            def onclick(event):
                if event.xdata is not None and event.ydata is not None:
                    coords.append((event.xdata, event.ydata))
                    plt.scatter(
                        event.xdata, event.ydata, color='black', marker='x')
                    plt.draw()
                    print('Clicked:', event.xdata, event.ydata)
            plt.gcf().canvas.mpl_connect('button_press_event', onclick)
            plt.show(block=True)

            if input(
                "Do you want to repick your points? (y or n)"
                ).lower() != 'y':
                break
        return np.array(coords)
    
    ############################################################################

    def overtone_points(self, axvline=None):
        '''Interactive picking of points on spectrogram for overtone peaks'''
        # User picks overtone peaks
        print(
            "Please pick one point on each overtone, it does not have to be "
            "at the center of the doppler."
        )
        while True:
            peaks, time_peaks = [], []
            plt.figure()
            plt.pcolormesh(
                self.times, self.frequencies, self.spec, shading='gouraud', 
                cmap='pink_r', vmin=self.vmin, vmax=self.vmax
            )
            if axvline is not None:
                plt.axvline(x=axvline, c='#377eb8', ls='--')
            def onclick(event):
                if event.xdata is not None and event.ydata is not None:
                    peaks.append(event.ydata)
                    time_peaks.append(event.xdata)
                    plt.scatter(
                        event.xdata, event.ydata, color='black', marker='x'
                        )
                    plt.draw()
                    print('Clicked:', event.xdata, event.ydata)
            plt.gcf().canvas.mpl_connect('button_press_event', onclick)
            plt.show(block=True)
            if input(
                "Do you want to repick your points? (y or n)"
                ).lower() != 'y':
                break
        return peaks, time_peaks

    ############################################################################

    def time_window_points(self, tobs, fobs):
        '''Interactive picking of time window for inversion'''
        # User picks time window for inversion
        print(
            'Please pick two points on the spectrogram that correspond to the '
            'start and end of the time window you want pull data from in the '
            'inversion.'
        )
        while True:
            set_time = []
            plt.figure()
            plt.pcolormesh(
                self.times, self.frequencies, self.spec, shading='gouraud', 
                cmap='pink_r',vmin=self.vmin, vmax=self.vmax
            )
            plt.scatter(tobs, fobs, color='black', marker='x')
            def onclick(event):
                if event.xdata is not None:
                    set_time.append(event.xdata)
                    plt.scatter(
                        event.xdata, event.ydata, color='red', marker='x'
                    )
                    plt.draw()
                    print('Clicked:', event.xdata, event.ydata)
            plt.gcf().canvas.mpl_connect('button_press_event', onclick)
            plt.show(block=True)
            if input(
                "Do you want to repick your points? (y or n)"
                ).lower() != 'y':
                break
        start_time, end_time = set_time[:2]
        return start_time, end_time
    
    ############################################################################

    def single_doppler_data(self, file_name, BASE_DIR):
        """
        Pick the points for the doppler shift. Specific to Seppi 2025 data 
        structure and flightradar25 information needed.

        Args:
            spec (numpy.ndarray): The spectrogram data.
            times (numpy.ndarray): The time array.
            frequencies (numpy.ndarray): The frequency array.
            vmin (float): The minimum amplitude value for the center line of the 
                spectrogram. Used for adjusting colorbar.
            vmax (float): The maximum amplitude value for the center line of the 
                spectrogram. Used for adjusting colorbar.
            month (int): The month of the data.
            day (int): The day of the data.
            flight (int): The flight number.
            sta (int or str): The station identifier.
            equip (str): The equipment identifier.
            closest_time (float): The time of closest approach.
            start_time (float): The start time of the spectrogram, to save for 
                future reference on plotting the spectrogram.
            make_picks (bool): If you come to this function and no picks exist, 
                it will allow you to make new picks.

        Returns:
            list: The list of picks the user picked along the most prominent 
                overtone.
        """

        if Path(f'{BASE_DIR}/{file_name}').exists():
            coords = []
            with open(f'{BASE_DIR}/{file_name}', 'r') as file:
                for line in file:
                    pick_data = line.split(',')
                    coords.append((float(pick_data[0]), float(pick_data[1])))
                if len(pick_data) == 4:
                    start_time = float(pick_data[2])
                else:
                    plt.figure()
                    plt.pcolormesh(
                        self.times, self.frequencies, self.spec, 
                        shading='gouraud', cmap='pink_r', vmin=self.vmin, 
                        vmax=self.vmax
                        )
                    plt.scatter(
                        [coord[0] for coord in coords], 
                        [coord[1] for coord in coords], 
                        color='black', marker='x')
                    plt.show()
                    correct_time = input(
                        "No start time found in file. Do your picks line" \
                        " up with this signal?(y/n): ")
                    if correct_time == 'y': 
                        start_time = self.tarrive - self.spec_window
                        # Rewrite file with start_time as third column 
                        # and move \n to next column
                        with open(f'{BASE_DIR}/{file_name}', 'r') as file:
                            lines = file.readlines()
                        with open(f'{BASE_DIR}/{file_name}', 'w') as file:
                            for line in lines:
                                pick_data = line.strip().split(',')
                                # Only keep first two columns, append start_time 
                                # then move \n to next column
                                if len(pick_data) >= 2:
                                    wline = f'{pick_data[0]},{pick_data[1]}'
                                    wline += f',{start_time},\n'
                                    file.write(wline)
                    else:
                        return [], None
            file.close()  
            return np.array(coords), start_time
        
        elif self.make_picks:
            start_time = self.tarrive - self.spec_window
            coords = self.doppler_points()
            if self.save_picks:
                make_base_dir(BASE_DIR)
                with open(f'{BASE_DIR}/{file_name}', 'w') as file:
                    for coord in coords:
                        file.write(f'{coord[0]},{coord[1]},{start_time},\n')
                file.close()
            return np.array(coords), start_time
        
        else:
            return [], None

    ############################################################################

    def overtone_data(self, t0, file_name, BASE_DIR):
        """
        Pick the points for the overtone shift. Specific to Seppi 2025 data 
        structure and flightradar25 information needed.

        Args:
            spec (numpy.ndarray): The spectrogram data.
            times (numpy.ndarray): The time array.
            frequencies (numpy.ndarray): The frequency array.
            vmin (float): The minimum amplitude value for the center line of the 
                spectrogram. Used for adjusting colorbar.
            vmax (float): The maximum amplitude value for the center line of the 
                spectrogram. Used for adjusting colorbar.
            month (int): The month of the data.
            day (int): The day of the data.
            flight (int): The flight number.
            sta (int or str): The station identifier.
            equip (str): The equipment identifier.
            closest_time (float): The time of closest approach.
            start_time (float): The start time of the spectrogram, to save for 
                future reference on plotting the spectrogram.
            t0 (float): The estimated acoustic wave arrival time.
            tarrive (float): The initial calculated time of acoustic wave 
                arrival.
            make_picks (bool): If you come to this function and no picks exist,     
                it will allow you to make new picks.

        Returns:
            list: List of frequencies picked by user along different overtones.
            list: List of times corresponding to the picked frequencies.
        """

        if Path(f'{BASE_DIR}/{file_name}').exists():

            peaks = []
            time_peaks = []
            with open(f'{BASE_DIR}/{file_name}', 'r') as file:
                for line in file:
                    pick_data = line.split(',')
                    peaks.append(float(pick_data[1]))
                    time_peaks.append(float(pick_data[0]))
            file.close()  
            return peaks, time_peaks

        elif self.make_picks:
            peaks, time_peaks = self.overtone_points(axvline=t0)
            if self.save_picks:
                make_base_dir(BASE_DIR)
                with open(f'{BASE_DIR}/{file_name}', 'w') as file:
                    for p, f in zip(peaks, time_peaks):
                        file.write(f'{f},{p},\n')
                file.close()
            return peaks, time_peaks
        else:
            return [], []
        
    ############################################################################

    def auto_picks_full(self, peaks, time_peaks, corridor_width, t0, v, d0, 
                            c, sigma_prior):
        """
        Get automatic picks for all overtones.

        Args:
            peaks (list): List of peak frequencies.
            time_peaks (list): List of times corresponding to the peaks.
            times (np.ndarray): Array of time values from fft.
            frequencies (np.ndarray): Array of frequency values from fft.
            spec (np.ndarray): Spectrogram data from fft.
            corridor_width (float): Width of the corridor for picking.
            t0 (float): Model parameter for the arrival time.
            v (float): Model parameter for the velocity.
            d0 (float): Model parameter for the distance.
            c (float): Model parameter for the speed of sound.
            sigma_prior (float): Prior uncertainty for the model parameters.
            vmax (float): Maximum amplitude value for peak detection.

        Returns:
            list: List of observed times.
            list: List of observed frequencies.
            list: List of counts of peaks associated with each overtone, 
                for indexing.
            list: List of fundamental frequencies calculated for each peak.
        """

        peaks_assos = []
        fobs = []
        tobs = []
        fs_array = []
    
        for pp in range(len(peaks)):
            tprime = time_peaks[pp]
            ft0p = peaks[pp]

            doppler_calc = DopplerCalc(v, d0, t0, c)
            fs = doppler_calc.calc_fs(tprime, ft0p)
            fs_array.append(fs)

            maxfreq = []
            coord_inv = []
            ttt = []

            ft = doppler_calc.calc_ft(self.times, fs)
 
            for t_f in range(len(self.times)):

                upper = int(ft[t_f] + corridor_width)
                lower = int(ft[t_f] - corridor_width)
                #find closest index to upper and lower in frequencies array
                lower_index = np.argmin(np.abs(self.frequencies - lower))
                upper_index = np.argmin(np.abs(self.frequencies - upper))
                if lower < 0:
                    lower = 0
                elif lower >= 250:
                    continue
                else:
                    pass
                if upper > 250:
                    upper = 250

                tt = self.spec[lower_index:upper_index, t_f]

                max_amplitude_index,_ = find_peaks(
                    tt, prominence = 15, wlen=10, height=self.vmax*0.1)
                
                if len(max_amplitude_index) == 0:
                    continue

                maxa = np.argmax(tt[max_amplitude_index])

                # Get the corresponding index into tt
                peak_idx = int(max_amplitude_index[maxa])
                freq_index = peak_idx + int(np.round(lower_index,0))
                # Now map it to frequency
                max_amplitude_frequency = self.frequencies[freq_index] 

                maxfreq.append(max_amplitude_frequency)
                coord_inv.append((self.times[t_f], max_amplitude_frequency))
                ttt.append(self.times[t_f])

            if len(ttt) > 0 and fs <= 230:
                coord_inv_array = np.array(coord_inv)
                mtest = [v, d0, t0, c, fs]
                invert_f = DopplerInversion(coord_inv_array[:,1], 
                                            coord_inv_array[:,0], mtest, 
                                            sigma_prior, num_iterations=2)
                mtest,_,_,_,_ = invert_f.full_inversion(
                                            [len(coord_inv_array[:,1])]
                                            )
                doppler_calc = DopplerCalc(
                                    mtest[0], mtest[1], mtest[2], mtest[3]
                                    )
                ft = doppler_calc.calc_ft(ttt, mtest[4])
                delf = np.array(ft) - np.array(maxfreq)

                count = 0
                for i in range(len(delf)):
                    if np.abs(delf[i]) <= (4):
                        fobs.append(maxfreq[i])
                        tobs.append(ttt[i])
                        count += 1
                peaks_assos.append(count)
            elif fs > 230:
                for i in range(len(ttt)):
                    fobs.append(maxfreq[i])
                    tobs.append(ttt[i])
                peaks_assos.append(len(maxfreq))
            else:
                peaks_assos.append(0)

        return tobs, fobs, peaks_assos, fs_array
    
    ############################################################################

    def final_data(self, tobs, fobs, fs_array, peaks_assos, file_name=None, 
                   BASE_DIR=None):
        """
        Pick the points for the time window used for the inversion. 
        Specific to Seppi 2025 data structure and flightradar25 information 
        needed.

        Args:
            tobs (list): The time array.
            fobs (list): The frequency array.
            fs_array (list): The sampling frequency array.
            peaks_assos (list or bool): The number of peaks associated with 
                each overtone.
            file_name (str): The file path to save picks.
            BASE_DIR (str): The base directory path.

        Returns:
            list: The time array, including data for all overtones.
            list: The frequency array, including data for all overtones.
            list: The number of data points associated with each overtone,  
                for indexing purposes.
        """

        if file_name is not None and Path(f'{BASE_DIR}/{file_name}').exists():
            set_time = []
            with open(f'{BASE_DIR}/{file_name}', 'r') as file:
                for line in file:
                    pick_data = line.split(',')
                    set_time.append(float(pick_data[0]))
            file.close()  
            if len(set_time) <= 1:
                return tobs, fobs, peaks_assos
            s_time = set_time[0]
            e_time = set_time[1]
            ftobs = []
            ffobs = []
        
            peak_ass = []
            cum = 0
            
            for p in range(len(fs_array)):
                count = 0
                for j in range(cum,cum+peaks_assos[p]):
                    if tobs[j] >= s_time and tobs[j] <= e_time:
                        ftobs.append(tobs[j])
                        ffobs.append(fobs[j])
                        count += 1
                cum = cum + peaks_assos[p]
            
                peak_ass.append(count)
            peaks_assos = peak_ass
            tobs = ftobs
            fobs = ffobs

            return tobs, fobs, peaks_assos

        elif self.make_picks:
            start_time, end_time = self.time_window_points(tobs, fobs)
            if self.save_picks:
                make_base_dir(BASE_DIR)
                with open(f'{BASE_DIR}/{file_name}', 'w') as file:
                    for tt in [start_time, end_time]:
                        file.write(f'{tt},\n')
                file.close()

            # Filter picks to only those within the selected time window
            ftobs, ffobs, peak_ass = [], [], []
            cum = 0
            for p in range(len(fs_array)):
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
            return ftobs, ffobs, peaks_assos

        else:
            return tobs, fobs, peaks_assos