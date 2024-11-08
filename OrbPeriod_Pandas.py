#!/usr/bin/env python

import sys
import math
import pandas as pd
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

# Load the data from the eclipse file
# User privides file...
if len(sys.argv) == 2:
  infile = sys.argv[1]
# User doesn't provide file...
if len(sys.argv) == 1:
  fileList = glob.glob('/data/ops/opsdata/current/predicts/eclipse/G*')
  latestFile = max(fileList,key=os.path.getctime)
  print(latestFile)
  if query_yes_no("Run on this file?","yes") is True:
    infile = latestFile
  else:
    print("Please provide the name of the eclipse file that you wish to use (with full path):")
    infile = raw_input()
  print("Using default eclipse file:",infile)
# User types nonsense...
if len(sys.argv) > 2:
  print("Script accepts either 1 or 0 arguments! Exiting...")
  exit()

# Function to convert datetime to Modified Julian Date (MJD)
def datetime_to_mjd(dt):
    t = Time(dt, format='datetime')
    return t.mjd

# Function to calculate midEclipse and orbPeriod from a CSV file
def calculate_midEclipse_and_orbPeriod(infile):
    # Load the CSV file
    df = pd.read_csv(infile)

    # Drop the first and last rows
    df = df.drop(index=[0, len(df) - 1])
    
    # Convert the 'Start Time (UTCJFOUR)' and 'Stop Time (UTCJFOUR)' columns to datetime
    df['Start Time (UTCJFOUR)'] = pd.to_datetime(df['Start Time (UTCJFOUR)'], format='%j/%Y %H:%M:%S.%f')
    df['Stop Time (UTCJFOUR)'] = pd.to_datetime(df['Stop Time (UTCJFOUR)'], format='%j/%Y %H:%M:%S.%f')
    
    # Calculate the mid-point between the start and stop times
    midEclipse_dt = (df['Start Time (UTCJFOUR)'] + (df['Stop Time (UTCJFOUR)'] - df['Start Time (UTCJFOUR)']) / 2).tolist()

    # Convert mid-point times to Modified Julian Date (MJD)
    midEclipse = [datetime_to_mjd(dt) for dt in midEclipse_dt]
    
    # Calculate the times between each 'midEclipse' element and the next element in seconds
    orbPeriod = [(midEclipse[i+1] - midEclipse[i]) * 86400 for i in range(len(midEclipse) - 1)]

    return midEclipse, orbPeriod

# Example usage
file_path = infile
midEclipse, orbPeriod = calculate_midEclipse_and_orbPeriod(file_path)

# Print the results
#print("midEclipse:", midEclipse[:5])  # Displaying first 5 elements for brevity
#print("orbPeriod:", max(orbPeriod),min(orbPeriod))    # Displaying first 5 elements for brevity

# Remove the first None value from orbPeriod for plotting
midEclipse_plot = midEclipse[1:-1]
orbPeriod_plot = orbPeriod[1:]

# Determine what mission weeks the eclipse file covers
firstweek = 54622.0 # MJD of start of MW 001
#print(midEclipse[0],midEclipse[-1])
startweek = math.floor((midEclipse[0] - firstweek)/7) + 1
endweek = math.floor((midEclipse[-1] - firstweek)/7) + 1

# Determine the time at which each mission week starts
secweekstart = (startweek)*7 + firstweek
endweekstart = (endweek)*7 + firstweek
i=0
weektimes = []
while i < (endweek - startweek):
    starttime = (startweek + i)*7 + firstweek
    weektimes.append(starttime)
    i+=1

# Load the beta angles file
beta_file_path = '/Users/jeggen/FERMI/FSSCtesting/OrbPeriod/beta_angles.txt'
beta_df = pd.read_csv(beta_file_path, sep=r'\s+', header=None)

# Use the first six columns to form a timestamp and the last column for beta angles
timestamps = beta_df.iloc[:, :6]
beta_angles = beta_df.iloc[:, -1]

# Create datetime objects from the first six columns
timestamps['datetime'] = pd.to_datetime(timestamps.apply(lambda row: f'{row[0]}-{row[1]:02d}-{row[2]:02d} {row[3]:02d}:{row[4]:02d}:{row[5]:02d}', axis=1))

# Convert to Modified Julian Date (MJD)
betaMJD = [datetime_to_mjd(dt) for dt in timestamps['datetime']]
betaAngle = beta_angles.tolist()

# Print the first few entries of betaMJD and betaAngle for verification
#print("Beta arrays:",str(len(betaMJD)),str(len(betaAngle)))
#betaMJD[:5], betaAngle[:5]

# Calculate start times of each mission week
mission_week_start = 54622.0
mission_week_duration = 7  # days
min_mjd = min(midEclipse_plot)
max_mjd = max(midEclipse_plot)

startMW = np.arange(mission_week_start, max_mjd, mission_week_duration)
startMW = startMW[startMW >= min_mjd]
numMW = []
for start in startMW:
    num = int((start - mission_week_start + 1)/mission_week_duration)
    numMW.append(num)

# Plot midEclipse_plot vs orbPeriod_plot and betaMJD vs betaAngle on different axes
fig, ax1 = plt.subplots(figsize=(12, 8))

# Create a secondary y-axis to plot Beta Angle vs. Beta MJD
color = 'tab:red'
ax1.set_ylabel('Beta Angle (degrees)', color=color)
ax1.plot(betaMJD, betaAngle, marker='x', linestyle='-', color=color, label='Beta Angle')
ax1.tick_params(axis='y', labelcolor=color)
ax1.axhline(y=24, color='magenta', linestyle='-.', alpha=0.7)
ax1.axhline(y=-24, color='magenta', linestyle='-.', alpha=0.7)
ax1.axhline(y=14, color='teal', linestyle='-.', alpha=0.7)
ax1.axhline(y=-14, color='teal', linestyle='-.', alpha=0.7)
ax1.legend(loc='lower left')

# Plot Orbital Period vs. Mid Eclipse
ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_xlabel('Modified Julian Date (MJD)')
ax2.set_ylabel('Orbital Period (seconds)', color=color)
ax2.plot(midEclipse_plot, orbPeriod_plot, marker='o', linestyle='-', color=color, label='Orbital Period')
ax2.tick_params(axis='y', labelcolor=color)
ax2.legend(loc='lower right')

# Print a buffer line
print('\n')

# Plot vertical lines for the start times of each mission week and calculate stats
for i in range(len(startMW) - 1):
    start_time = startMW[i]
    end_time = startMW[i + 1]
    ax1.axvline(x=start_time, color='green', linestyle='--', alpha=0.7)

    # If this is the end of the range, plot one more line for the end of the last MW
    if end_time == startMW[-1]:
        ax1.axvline(x=end_time, color='green', linestyle='--', alpha=0.7)
    
    # Calculate statistics within each mission week
    mask = (midEclipse_plot >= start_time) & (midEclipse_plot < end_time)
    mw_orbPeriod = np.array(orbPeriod_plot)[mask]
    mw_betaAngle = np.array(betaAngle)[(np.array(betaMJD) >= start_time) & (np.array(betaMJD) < end_time)]
    
    mean_orbPeriod = round(np.mean(mw_orbPeriod),2)
    min_orbPeriod = round(np.min(mw_orbPeriod),2)
    max_orbPeriod = round(np.max(mw_orbPeriod),2)
    min_betaAngle = round(np.min(mw_betaAngle),2)
    max_betaAngle = round(np.max(mw_betaAngle),2)
    
    # Print the mission week and statistics
    print(f'Mission Week {numMW[i+1]}:')
    print(f'Mean Orbital Period: {mean_orbPeriod}')
    print(f'Start/Stop MJD: {start_time}'+'/'+f'{end_time}')
    print(f'Min/Max Orbital Period: {min_orbPeriod}'+'/'+f'{max_orbPeriod}')
    print(f'Min/Max Beta Angle: {min_betaAngle}'+'/'+f'{max_betaAngle}\n')
    
    # Annotate the mission week number in the plot
    mid_point = (start_time + end_time) / 2
    ax1.text(mid_point, ax1.get_ylim()[1], f'{numMW[i+1]}', ha='center', 
             va='top', fontsize=10, color='black', 
             bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))


# Set limits based on midEclipse_plot
plt.xlim(min(midEclipse_plot), max(midEclipse_plot))

plt.title('Orbital Period and Beta Angle vs. Time')
plt.grid(True)
fig.tight_layout()
plt.show()
