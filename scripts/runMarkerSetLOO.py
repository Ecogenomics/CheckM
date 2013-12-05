import os

for t in ['0.95', '0.96','0.97','0.98','0.99','1.0']:
  os.system('./markerSetLOO.py -u ' + t + ' -s ' + t + ' -r 4')
