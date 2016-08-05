import numpy as np
from scipy.special import 
from scipy import linalg

# expand around a corner

# interior angle fraction of corner (note, the angle is pi/angle_fraction)
angle_frac = 10

# side lengths
side_A = 1
side_B = 1

# number of boundary samples
N_B = 10
# number of interior samples
N_I = 10
# length of expansion
K = 10

# assemble matrix
A = np.zeros(N_B+N_1,K)

# perform QR decomposition
