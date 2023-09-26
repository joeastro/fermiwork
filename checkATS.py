#!/usr/bin/env python

####################
# The purpose of this script is to perform basic checks on the real vs. fake
# ATS products for a given week. The script will automatically download the 
# fake poducts but the user will have to manually download the real products 
# as the page they come from requires password access.
#
# Usage: python checkATS.py [MW]
#
# Author: Joe Eggen
# Dec. 12, 2022 - First version
####################

import os
import sys
import urllib.request
from datetime import datetime
import numpy as np

urlBase = 'https://www.slac.stanford.edu/exp/glast/ops/'

# Load the data from the science timeline
if len(sys.argv) == 2:
  missionWeek = sys.argv[1]
else:
  print('Error: please provide the mission week as an arguement')

# Retrieve the fake products from the web
fakePlan = missionWeek + 'plan.fake'
fakeRev = missionWeek + 'rev.fake'

# Keep track of warnings
warningCount = 0

# Check if the fake products are already present, and download them if not
if os.path.exists(fakePlan) is False:
  print('Downloading ' + urlBase + fakePlan)
  urllib.request.urlretrieve(urlBase + fakePlan, fakePlan)
else: 
  print(fakePlan + ' already present. Continuing...')

if os.path.exists(fakeRev) is False:
  print('Downloading ' + urlBase + fakeRev)
  urllib.request.urlretrieve(urlBase + fakeRev, fakeRev)
else: 
  print(fakeRev + ' already present. Continuing...')

# Check for the real products. The user will have to download them maunually
# NOTE: The user can modify the naming convention of the real product files 
#       as they wish, so long as the MW number remains in the filename.
realPlan = 'MW' + missionWeek + '_planned.csv'
realRev = 'MW' + missionWeek + '_review.csv'

# Check if the real products are present
if os.path.exists(realPlan) is False:
  print('File not found:',realPlan)
  print('You must manually download the real ATS products to use this script!')
  sys.exit()
else:
  print(realPlan + ' already present. Continuing...')

if os.path.exists(realRev) is False:
  print('File not found',realRev)
  print('You must manually download the real ATS products to use this script!')
  sys.exit()
else:
  print(realRev + ' already present. Continuing...')

print('')

# The "review" file is the one that specifies LPASTART/LPASTOP command timing, 
# so  retrieve all the lines for those commands
reviewLines = []
with open (realRev, 'rt') as infile:
  for line in infile:
    if 'LPASTART' in line or 'LPASTOP' in line:
      splitline = line.split(',')
      reviewLines.append((splitline[0],splitline[1],splitline[2],splitline[5]))

# The ATS should *not* begin with an LPASTOP command 
if reviewLines[0][3] == 'LPASTART':
  print('ATS begins with LPASTART command... OK')
else:
  print('############################################################')
  print('WARNING!! The ATS does NOT begin with an LPASTART command!')
  print('############################################################')
  warningCount += 1

timingErr = False
n=1
while n + 1 < len(reviewLines):
  #print(n)
  stop = datetime.fromisoformat(reviewLines[n][2])
  start = datetime.fromisoformat(reviewLines[n+1][2])
  #print(stop,start)
  diff = start - stop
  if diff.seconds < 10:
    print('\n')
    print('WARNING!!! Time between LPASTOP and LPASTART commands is '+str(diff.seconds)+' < 10 seconds!')
    print('LPASTOP @ '+str(stop))
    timingErr = True
    warningCount += 1
  n += 2

if timingErr == False:
  print('Timing check... OK.')
  print('  All LPASTOP/LPASTART commands occur at least 10 seconds apart.')
else:
  print('WARNING!!! Timing check FAILED!')
  print('\n')

# Perform the diffs for the real vs fake products
print('Executing diffs...')
print('diff '+realPlan+' '+fakePlan)
print('####################')
os.system('diff '+realPlan+' '+fakePlan)
print('####################\n')

print('diff '+realRev+' '+fakeRev)
print('####################')
os.system('diff '+realRev+' '+fakeRev)
print('####################')

print('ATS check completed')
print('Number of warnings:',str(warningCount))
