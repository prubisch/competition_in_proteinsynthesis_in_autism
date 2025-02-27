import os 
import sys
import numpy as np

MCELL_PATH = os.environ.get('MCELL_PATH','')
sys.path.append(os.path.join(MCELL_PATH, 'lib'))
import mcell as m




print('Import of MCELL was sucessful')
