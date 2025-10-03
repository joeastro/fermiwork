#!/usr/bin/env python

# Orbit period (synodic) prediction for Fermi
# The purpose of this script is to use a weekly eclipse file to generate 
# a plot & print values for the expected orbital period for upcoming 
# mission  weeks. This should aid in providing more accurate initial 
# estimates for profile repeat periods prior to running the planning 
# software (i.e. TAKO).
#
# Written by Joe Eggen: Nov. 08, 2024
#
# 10/03/2025 - added feature to print the date/time of a profile change
#              based on beta angle (detect_beta_crossings_in_range).

import os
import sys
import math
import csv
import glob
from datetime import datetime, timedelta
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via input() and return their answer.

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
        choice = input().lower()
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
    infile = input()
  print("Using default eclipse file:",infile)
# User types nonsense...
if len(sys.argv) > 2:
  print("Script accepts either 1 or 0 arguments! Exiting...")
  exit()

# Function to convert datetime to Modified Julian Date (MJD)
def datetime_to_mjd(dt):
    t = Time(dt, format='datetime')
    return t.mjd

def read_csv(infile):
    with open(infile, newline='') as file:
        reader = csv.DictReader(file)
        data = [row for row in reader]
    return data

def calculate_midEclipse_and_orbPeriod(infile):
    # Load the CSV file as a list of dictionaries
    data = read_csv(infile)

    # Drop the first and last rows
    data = data[1:-1]

    # Convert the 'Start Time (UTCJFOUR)' and 'Stop Time (UTCJFOUR)' columns to datetime
    for row in data:
        row['Start Time (UTCJFOUR)'] = datetime.strptime(row['Start Time (UTCJFOUR)'], '%j/%Y %H:%M:%S.%f')
        row['Stop Time (UTCJFOUR)'] = datetime.strptime(row['Stop Time (UTCJFOUR)'], '%j/%Y %H:%M:%S.%f')

    # Calculate the mid-point between the start and stop times
    midEclipse_dt = [(row['Start Time (UTCJFOUR)'] + (row['Stop Time (UTCJFOUR)'] - row['Start Time (UTCJFOUR)']) / 2) for row in data]

    # Convert mid-point times to Modified Julian Date (MJD)
    midEclipse = [datetime_to_mjd(dt) for dt in midEclipse_dt]

    # Calculate the orbital period between each 'midEclipse' element and the next element in seconds
    orbPeriod = [(midEclipse[i + 1] - midEclipse[i]) * 86400.0 for i in range(len(midEclipse) - 1)]

    return midEclipse, orbPeriod

def detect_beta_crossings_in_range(betaMJD, betaAngle, datetimes):
    """
    Detect when beta angle crosses threshold values within a specific time range.
    """
    if len(betaAngle) < 2:
        return
    
    # Define thresholds
    thresholds = [-24, -14, 14, 24]
    
    # Track previous state to detect crossings
    prev_angle = betaAngle[0]
    
    for i in range(1, len(betaAngle)):
        curr_angle = betaAngle[i]
        curr_datetime = datetimes[i]
        day_of_year = curr_datetime.timetuple().tm_yday

        # Check for crossings of each threshold
        for threshold in thresholds:
            # Detect crossing (previous and current angles are on opposite sides of threshold)
            if (prev_angle <= threshold < curr_angle) or (prev_angle >= threshold > curr_angle):
                
                # Format the date/time string
                date_time_str = curr_datetime.strftime("%Y-%m-%d %H:%M:%S")
                
                
                # Determine which range we're entering and print appropriate message
                if threshold == 14 and curr_angle > 14:
                    if curr_angle <= 24:
                        print(f"Switch to +50/-60 deg. asymmetric profile by {date_time_str} (DoY {day_of_year})")
                    else:
                        print(f"Switch to +50 deg. modified sine profile by {date_time_str} (DoY {day_of_year})")
                elif threshold == -14 and curr_angle < -14:
                    if curr_angle >= -24:
                        print(f"Switch to -50/+60 deg. asymmetric profile by {date_time_str} (DoY {day_of_year})")
                    else:
                        print(f"Switch to -50 deg. modified sine profile by {date_time_str} (DoY {day_of_year})")
                elif threshold == 24 and curr_angle > 24:
                    print(f"Switch to +50 deg. modified sine profile by {date_time_str} (DoY {day_of_year})")
                elif threshold == -24 and curr_angle < -24:
                    print(f"Switch to -50 deg. modified sine profile by {date_time_str} (DoY {day_of_year})")
                elif threshold == -24 and curr_angle > -24:
                    print(f"Switch to -50/+60 deg. asymmetric  profile by {date_time_str} (DoY {day_of_year})")
                elif threshold == 14 and curr_angle < 14 and prev_angle > 14:
                    if curr_angle >= -14:
                        print(f"Switch to +-50 deg. symmetric profile by {date_time_str} (DoY {day_of_year})")
                elif threshold == -14 and curr_angle > -14 and prev_angle < -14:
                    if curr_angle <= 14:
                        print(f"Switch to +-50 deg. symmetric profile by {date_time_str} (DoY {day_of_year})")
        prev_angle = curr_angle

# Load the beta angles file
beta_file_path = '/Home/lhea2/gsscops/beta_angles.txt'

# Get unique filename for output plot
figstem = infile.split('/')[-1]
figName = figstem+'.predict.png'

# Example usage
file_path = infile
midEclipse, orbPeriod = calculate_midEclipse_and_orbPeriod(file_path)

# Print the results
#print("midEclipse:", midEclipse[:5])  # Displaying first 5 elements for brevity
#print("orbPeriod:", max(orbPeriod),min(orbPeriod))    # Displaying first 5 elements for brevity

# Remove the first None value from orbPeriod for plotting
midEclipse_plotStr = np.array(midEclipse[1:-1])
midEclipse_plot = midEclipse_plotStr.astype(float)
orbPeriod_plotStr = np.array(orbPeriod[1:])
orbPeriod_plot = orbPeriod_plotStr.astype(float)

# Determine what mission weeks the eclipse file covers
firstweek = 54622.0 # MJD of start of MW 001
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

# Use the first six columns to form a timestamp and the last column for beta angles
# Since we can't use the pandas package the process is somewhat convoluted
def read_first_six_columns(filename):
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile,delimiter=' ')
        data = [row[:6] for row in reader]  # Slice to get the first six columns
    return data
def read_last_column(filename):
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile,delimiter=' ')
        data = [row[-1] for row in reader if row]  # Get the last column of each row, ensure the row is not empty
    return data
timestamps = read_first_six_columns(beta_file_path)
beta_angles = read_last_column(beta_file_path)

# Create datetime objects from the first six columns
date_format = '%Y %m %d %H %M %S'
datetimes = [datetime.strptime(str(row[0])+' '+str(row[1])+' '+ str(row[2])+' '+ str(row[3])+' '+ str(row[4])+' '+ str(row[5]),date_format) for row in timestamps]

# Convert to Modified Julian Date (MJD)
betaMJD = np.array([datetime_to_mjd(dt) for dt in datetimes])
betaAngleStr = np.array(beta_angles)
betaAngle = betaAngleStr.astype(float)

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
ax2.set_ylabel('Synodic Orbital Period (seconds)', color=color)
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
    mw_betaAngle = betaAngle[(np.array(betaMJD) >= start_time) & (np.array(betaMJD) < end_time)]
    mw_betaAngle = np.array(mw_betaAngle)

    # Get beta angle data for this mission week
    beta_mask = (np.array(betaMJD) >= start_time) & (np.array(betaMJD) < end_time)
    mw_betaMJD = np.array(betaMJD)[beta_mask]
    mw_datetimes = np.array(datetimes)[beta_mask]
    
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
    print(f'Min/Max Beta Angle: {min_betaAngle}'+'/'+f'{max_betaAngle}')

    # Check for beta angle threshold crossings within this mission week
    if len(mw_betaAngle) > 1:  # Need at least 2 points to detect crossings
        detect_beta_crossings_in_range(mw_betaMJD, mw_betaAngle, mw_datetimes)
    
    print()  # Add blank line between mission weeks
    
    # Annotate the mission week number in the plot
    mid_point = (start_time + end_time) / 2
    ax1.text(mid_point, ax1.get_ylim()[1], f'{numMW[i+1]}', ha='center', 
             va='top', fontsize=10, color='black', 
             bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))


# Set limits based on midEclipse_plot
plt.xlim(min(midEclipse_plot), max(midEclipse_plot))

plt.title('Synodic Orbital Period and Beta Angle vs. Time')
plt.grid(True)
fig.tight_layout()
#plt.show()
print("Plot of orbit period saved to file: "+figName)
plt.savefig(figName)
