import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Define the parameters
rpm = 2500
int_frequency = 3 * rpm / 60  # Convert rpm to Hz
omega = 2 * np.pi * int_frequency  # Angular frequency

# Time array
t = np.linspace(0, 0.1, 1000)  # Zoom in by reducing the time range

# Define the waves with harmonics
def generate_wave(t, omega, phase, harmonics=6):
    wave = np.zeros_like(t)
    for n in range(1, harmonics + 1):
        wave+= (1/n) * np.sin(n * omega * t + phase)
        wave += (1/(n+0.5) * np.sin((n+0.5) * omega * t + phase))
    return wave

wave1 = generate_wave(t, omega, 0)

# Create the figure and axis
fig, ax = plt.subplots()
line1, = ax.plot(t, wave1, label='Wave 1')

# Add the frequency text annotation
frequency_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_ylim(-1.5, 6)
# Update function for animation
def update(frame):
    global wave1, wave2, wave3, wave_sum
    frequency = int_frequency * (1 - 0.001 * frame)
    amplitude = 1 - 0.001 * frame
    omega = 2 * np.pi * frequency
    wave1 = amplitude *generate_wave(t, omega,(1/3)*np.pi)
    line1.set_ydata(wave1)

    # Update the frequency text annotation
    frequency_text.set_text(f'Frequency: {frequency:.2f} Hz')

    return line1, frequency_text

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=1000, interval=200, blit=True)

# Add the frequency text annotation
frequency_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


# Pause the animation
#ani.event_source.stop()

# Display the animation
plt.show()