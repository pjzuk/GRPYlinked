from sys import argv
import sys

fname = argv[1]
data = []


with open(fname,'r') as dataFile:
  fistline = dataFile.readline()
  for line in dataFile:
    data.append([ float(i) for i in line.split() ])

with open('out_msd.dat','w') as outFile:
  # setup toolbar
  toolbar_width = 10
  sys.stdout.write("0%   DONE ")
  sys.stdout.flush()
  sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line 
 
  outFile.write(str(0) + ' ' + str(0)+ '\n')
  prev = 0
  for steps in range(1,int(len(data)/100.)):
    t = data[steps][1] - data[0][1]
    msd = 0.
    ln = 0.
    for i in range(len(data)-steps):
      ln += 1
      msd += (data[i+steps][3] - data[i][3])**2. +  (data[i+steps][4] - data[i][4])**2. + (data[i+steps][5] - data[i][5])**2.
    msd = msd/ln
    outFile.write(str(t) + ' ' + str(msd) + '\n')
    outFile.flush()
    this = int(float(steps)/float( len( range(1,int(len(data)/100.)) ))*100)
    if ( prev != this):
      prev = this
      sys.stdout.write(str(int(this))+"%")
      sys.stdout.write("\b" * (toolbar_width + 1)) # return to start of line 
      sys.stdout.flush()
  sys.stdout.write("\n") 
