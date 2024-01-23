import numpy as np

eps = np.float64(1.0)
for i in range(100):
    eps /= 2
    one = eps+1.0

    # INSERT SOMETHING HERE

DBL_EPSILON = np.finfo(np.float64).eps
assert(DBL_EPSILON == eps)
