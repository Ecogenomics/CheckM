import os

for r in ['Archaea', 'Bacteria']:
  for t in ['1.0', '0.99', '0.98', '0.97', '0.96', '0.95']:
    os.system('./markerSetTest.py -T ' + r + ' -u ' + t + ' -s ' + t)
