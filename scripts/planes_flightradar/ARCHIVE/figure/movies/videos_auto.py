import pandas as pd
import matplotlib.pyplot as plt
from celluloid import Camera
import matplotlib.animation as animation
from matplotlib import rc
from moviepy.editor import VideoFileClip, CompositeVideoClip
from moviepy.video.io.bindings import mplfig_to_npimage

# Read the data from the text files
# Assuming the text files contain two columns: 'lat' for latitude and 'lon' for longitude
data1 = pd.read_csv('file1.txt', sep="\t")
data2 = pd.read_csv('file2.txt', sep="\t")

fig, ax = plt.subplots()

camera = Camera(fig)
for i in range(len(data1)):
    # Plot data from file1
    ax.scatter(data1['lon'][:i], data1['lat'][:i], color='blue')
    # Plot data from file2 on top of file1 data
    ax.scatter(data2['lon'][:i], data2['lat'][:i], color='red')
    camera.snap()

animation = camera.animate()
animation.save('layered_plot.mp4')

# Read the data from the text files
# Assuming the text files contain two columns: 'lat' for latitude and 'lon' for longitude
data1 = pd.read_csv('file1.txt', sep="\t")
data2 = pd.read_csv('file2.txt', sep="\t")

fig, ax = plt.subplots()

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

def animate(i):
    ax.clear()
    # Plot data from file1
    ax.scatter(data1['lon'][:i], data1['lat'][:i], color='blue')
    # Plot data from file2 on top of file1 data
    ax.scatter(data2['lon'][:i], data2['lat'][:i], color='red')

ani = animation.FuncAnimation(fig, animate, frames=len(data1), repeat=True)

rc('animation', html='html5')
ani.save('layered_plot.mp4', writer=writer)

# Read the data from the text file
# Assuming the text file contains three columns: 'time' for timestamp, 'lat' for latitude and 'lon' for longitude
data = pd.read_csv('file.txt', sep="\t")

# Load the spectrogram video
video = VideoFileClip("spectrogram.mp4")

# Create a figure for the map
fig, ax = plt.subplots()

def make_frame(t):
    ax.clear()
    # Find the data up to time t and plot it on the map
    current_data = data[data['time'] <= t]
    ax.scatter(current_data['lon'], current_data['lat'], color='red')
    return mplfig_to_npimage(fig)

# Create a video clip for the map
map_clip = VideoClip(make_frame, duration=video.duration)

# Combine the spectrogram video and map video
final_clip = CompositeVideoClip([video, map_clip.resize(video.size)])

# Write the result to a file
final_clip.write_videofile("output.mp4")
