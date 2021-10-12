import numpy as np


def mode(x):
	u,c = np.unique(x,return_counts=True)
	return u[c.argmax()]
	
