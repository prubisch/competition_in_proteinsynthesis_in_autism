import os 
import sys
import numpy as np
import mcell as m

MCELL_PATH = os.environ.get('MCELL_PATH','')
sys.path.append(os.path.join(MCELL_PATH, 'lib'))


print('Import of MCELL was sucessful')
