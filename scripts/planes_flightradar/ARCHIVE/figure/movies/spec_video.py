import matplotlib.pyplot as plt
import obspy
import imageio.v2 as imageio	
	
n = "/scratch/naalexeev/NODAL/2019-03-14T23:00:00.000000Z.2019-03-15T00:00:00.000000Z.1132.mseed" 

tr = obspy.read(n)
images = []
for tim in range(1,30,1):
	try:
		tr[2].trim(tr[2].stats.starttime + (42 * 60) + tim, tr[2].stats.starttime + (43 * 60) +tim)
		
		t                  = tr[2].times()

		title    = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'
		
		fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6))     

		ax1.set(xlim=(0, 240), ylim=(-2000, 2000))
		ax1.plot(t, tr, 'k', linewidth=0.5)
		ax1.set_title(title)
		#ax1.axvline(x=tim, c = 'r', ls = '--')

		s,f,tm,im = ax2.specgram(tr[2], scale='dB', vmin = -120, vmax = 50)
		ax2.set_xlabel('Time - Seconds')
		ax2.axvline(x=tim, c = 'k', ls = '--')
		ax2.set_ylabel('Frequency (Hz)')

		ax3 =fig.add_axes([0.88, 0.1, 0.02, 0.37])
		plt.colorbar(mappable=im, cax=ax3, spacing='uniform', label='Relative Amplitude (dB)')

		# Save the plot as a PNG file
		plt.savefig(f'frame_{tim}.png')

		# Add the PNG file to the list of images
		images.append(imageio.imread(f'frame_{tim}.png'))
	except:
		continue
# Create a GIF from the images with 5 frames per second (one frame every 0.2 seconds)
imageio.mimsave('flight.gif', images, fps=4)
