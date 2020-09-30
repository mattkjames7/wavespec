import numpy as np

def ApplyWindowFunction(t,v,WindowFunction=None,Param=None):
	'''
	Apply a window function to a tim series.
	
	Inputs
	======
	t : float
		Time array 
	v : float 
		Time series data to be windowed
	WindowFunction : None | str
		If None - no window is applied, otherwise the string names the 
		window function to be applied (see below for list of functions)
	Param : float
		Sometimes a window function may be modified by some parameter,
		setting this keyword to None will force the routine to use a 
		default value where needed.
		
	Returns
	=======
	vw : float
		Time series data, v, with the appropriate window function 
		applied to it.
		
	Window Functions
	================
	Function			| Param
	--------------------|-------------
	None				|
	'cosine-bell'		|
	'hamming'			|
	'triangle'			|
	'welch'				|
	'blackman'			|
	'nuttall'			|
	'blackman-nuttall'	|
	'flat-top'			|
	'cosine'			|
	'gaussian'			|
	
	'''


	# get the appropriate window function and parameters
	WF = {	'none':				(_WFNone,0.0),
			'cosine-bell':		(_WFCosineBell,10.0),
			'hamming':			(_WFHamming,50.0),
			'hann':				(_WFHann,50.0),
			'triangle':			(_WFTriangle,50.0),
			'welch':			(_WFWelch,0.0),
			'blackman':			(_WFBlackman,0.0),
			'nuttall':			(_WFNuttall,0.0),
			'blackman-nuttall':	(_WFBlackmanNuttall,0.0),
			'flat-top':			(_WFFlatTop,0.0),
			'cosine':			(_WFCosine,10.0),
			'gaussian':			(_WFGaussian,(0.4,50.0))}
	
	Func,Pdef = WF.get(WindowFunction,(_WFNone,0.0))
	
	#check if anycustom parameters are being used
	if Param is None:
		P = Pdef
	else:
		P = Param
	
	#apply to data
	return Func(t,v,P)
	


def _WFNone(t,v,P):
	'''
	No window function - just return original array.
	'''
	
	return v


def _WFCosineBell(t,v,P=10.0):
	'''
	This will multiply the date by the Split cosine bell function.
	
	Inputs
	======
	t : float
		Time array
	v : float
		Time series to be windowed
	P : float
		Percentage of time series at each end to be part of the cosine.
		The remaining 100 - 2*P % is left unchanged (if you set to 50.0,
		then the whole window has the function applied to it).
		
	Returns
	=======
	out : float
		Windowed version of v
		

	'''

	#get the time range
	t0 = np.nanmin(t)
	t1 = np.nanmax(t)
	tr = t1 - t0
	
	#get a scaled time array
	ts = 100.0*(t - t0)/tr
	
	#work out the indices for each section
	i0 = np.where(ts < P)[0]
	i2 = np.where(ts > (100 - P))[0]

	#calculate the window function
	w = np.ones(t.size,dtype=v.dtype)

	#first section
	w[i0] = 0.5 + 0.5*np.cos((ts[i0]/P + 1.0)*np.pi)
	
	#last section
	w[i2] = 0.5 + 0.5*np.cos(np.pi*(ts[i2] - (100 - P))/P)
	
	
	#multiply by v
	out = v*w
	
	return out


def _WFHamming(t,v,P=50.0):
	'''
	This will multiply the date by the Hamming window function.
	
	Inputs
	======
	t : float
		Time array
	v : float
		Time series to be windowed
	P : float
		Percentage of window to have the function applied at each end
		(if you set to 50.0, then the whole window has the function 
		applied to it).
		
	Returns
	=======
	out : float
		Windowed version of v
	'''
	#get the time range
	t0 = np.nanmin(t)
	t1 = np.nanmax(t)
	tr = t1 - t0
	
	#get a scaled time array
	ts = 100.0*(t - t0)/tr
	
	#work out the indices for each section
	i0 = np.where(ts < P)[0]
	i2 = np.where(ts > (100 - P))[0]
	
	
	#calculate the window function
	w = np.ones(t.size,dtype=v.dtype)

	#first section
	w[i0] = 0.53836 - 0.46164*np.cos(np.pi*ts[i0]/P)
	
	#last section
	w[i2] = 0.53836 - 0.46164*np.cos(np.pi*(1.0 + (ts[i2] - (100 - P))/P))
	
	#multiply by v
	out = v*w

	return out

def _WFHann(t,v,P=50.0):
	'''
	This will multiply the date by the Hann window (sometimes 
	erroneously called the "Hanning" window). This is similar to the 
	Hamming window, but is touches zero at each end.
	
	Inputs
	======
	t : float
		Time array
	v : float
		Time series to be windowed
	P : float
		Percentage of window to have the function applied at each end
		(if you set to 50.0, then the whole window has the function 
		applied to it).
		
	Returns
	=======
	out : float
		Windowed version of v
	
	'''
	#get the time range
	t0 = np.nanmin(t)
	t1 = np.nanmax(t)
	tr = t1 - t0
	
	#get a scaled time array
	ts = 100.0*(t - t0)/tr
	
	#work out the indices for each section
	i0 = np.where(ts < P)[0]
	i2 = np.where(ts > (100 - P))[0]
	
	
	#calculate the window function
	w = np.ones(t.size,dtype=v.dtype)

	#first section
	w[i0] = 0.5 - 0.5*np.cos(np.pi*ts[i0]/P)
	
	#last section
	w[i2] = 0.5 - 0.5*np.cos(np.pi*(1.0 + (ts[i2] - (100 - P))/P))
	
	#multiply by v
	out = v*w

	return out
	
def _WFTriangle(t,v,P=50.0):
	'''
	This will multiply the date by the Triangle window.
	
	Inputs
	======
	t : float
		Time array
	v : float
		Time series to be windowed
	P : float
		Percentage of window to have the function applied at each end
		(if you set to 50.0, then the whole window has the function 
		applied to it).
		
	Returns
	=======
	out : float
		Windowed version of v
	
	'''
	#get the time range
	t0 = np.nanmin(t)
	t1 = np.nanmax(t)
	tr = t1 - t0
	
	#get a scaled time array
	ts = 100.0*(t - t0)/tr
	
	#work out the indices for each section
	i0 = np.where(ts < P)[0]
	i2 = np.where(ts > (100 - P))[0]

	#calculate the window function
	w = np.ones(t.size,dtype=v.dtype)
	
	#first part
	w[i0] = 1.0 - np.abs(ts[i0] - P)/P
	
	#second part
	w[i2] = 1.0 - np.abs(ts[i2] - (100 - P))/P
	
	out = w*v
	
	return out
	
def _WFWelch(t,v,P=5.0):
	'''
	This will multiply the date by the Welch window.
	
	Inputs
	======
	t : float
		Time array
	v : float
		Time series to be windowed
	P : float
		Percentage of window to have the function applied at each end
		(if you set to 50.0, then the whole window has the function 
		applied to it).
		
	Returns
	=======
	out : float
		Windowed version of v
	
	'''
	#get the time range
	t0 = np.nanmin(t)
	t1 = np.nanmax(t)
	tr = t1 - t0
	
	#get a scaled time array
	ts = 100.0*(t - t0)/tr
	
	#work out the indices for each section
	i0 = np.where(ts < P)[0]
	i2 = np.where(ts > (100 - P))[0]

	#calculate the window function
	w = np.ones(t.size,dtype=v.dtype)
	
	#first part
	w[i0] = 1 - ((ts[i0] - P)/P)**2
	
	#second part
	w[i2] = 1 - ((ts[i2] - (100 - P))/P)**2
	
	out = w*v	

	
	return out
	
def _WFBlackman(t,v,P=50.0):
	'''
	This will multiply the date by the Blackman window.
	
	Inputs
	======
	t : float
		Time array
	v : float
		Time series to be windowed
	P : float
		Percentage of window to have the function applied at each end
		(if you set to 50.0, then the whole window has the function 
		applied to it).
		
	Returns
	=======
	out : float
		Windowed version of v
	
	'''
	#get the time range
	t0 = np.nanmin(t)
	t1 = np.nanmax(t)
	tr = t1 - t0
	
	#get a scaled time array
	ts = 100.0*(t - t0)/tr
	
	#work out the indices for each section
	i0 = np.where(ts < P)[0]
	i2 = np.where(ts > (100 - P))[0]

	#calculate the window function
	w = np.ones(t.size,dtype=v.dtype)

	#some constants
	a0 = 7938.0/18608
	a1 = 9240.0/18608
	a2 = 1430.0/18608	

	
	#first part
	w[i0] = a0 - a1*np.cos(np.pi*ts[i0]/P) + a2*np.cos(2*np.pi*ts[i0]/P)
	
	#second part
	w[i2] = a0 - a1*np.cos(np.pi*(ts[i2] - (100 - 2*P))/P) + a2*np.cos(2*np.pi*(ts[i2] - (100 - 2*P))/P)
	
	out = w*v	

	return out
	
def _WFNuttall(t,v,P=50.0):
	'''
	This will multiply the date by the Nuttall window.
	
	Inputs
	======
	t : float
		Time array
	v : float
		Time series to be windowed
	P : float
		Percentage of window to have the function applied at each end
		(if you set to 50.0, then the whole window has the function 
		applied to it).
		
	Returns
	=======
	out : float
		Windowed version of v
	
	'''
	#get the time range
	t0 = np.nanmin(t)
	t1 = np.nanmax(t)
	tr = t1 - t0
	
	#get a scaled time array
	ts = 100.0*(t - t0)/tr
	
	#work out the indices for each section
	i0 = np.where(ts < P)[0]
	i2 = np.where(ts > (100 - P))[0]

	#calculate the window function
	w = np.ones(t.size,dtype=v.dtype)

	#some constants
	a0 = 0.355768
	a1 = 0.487396
	a2 = 0.144232
	a3 = 0.012604

	
	#first part
	w[i0] = a0 - a1*np.cos(np.pi*ts[i0]/P) + a2*np.cos(2*np.pi*ts[i0]/P) - a3*np.cos(3*np.pi*ts[i0]/P)
	
	#second part
	w[i2] = a0 - a1*np.cos(np.pi*(ts[i2] - (100 - 2*P))/P) + a2*np.cos(2*np.pi*(ts[i2] - (100 - 2*P))/P) - a3*np.cos(3*np.pi*(ts[i2] - (100 - 2*P))/P)
	
	out = w*v	

	return out


def _WFBlackmanNuttall(t,v,P=50.0):
	'''
	This will multiply the date by the Blackman-Nuttall window.
	
	Inputs
	======
	t : float
		Time array
	v : float
		Time series to be windowed
	P : float
		Percentage of window to have the function applied at each end
		(if you set to 50.0, then the whole window has the function 
		applied to it).
		
	Returns
	=======
	out : float
		Windowed version of v
	
	'''
	#get the time range
	t0 = np.nanmin(t)
	t1 = np.nanmax(t)
	tr = t1 - t0
	
	#get a scaled time array
	ts = 100.0*(t - t0)/tr
	
	#work out the indices for each section
	i0 = np.where(ts < P)[0]
	i2 = np.where(ts > (100 - P))[0]

	#calculate the window function
	w = np.ones(t.size,dtype=v.dtype)

	#some constants
	a0 = 0.3635819
	a1 = 0.4891775
	a2 = 0.1365995
	a3 = 0.0106411

	
	#first part
	w[i0] = a0 - a1*np.cos(np.pi*ts[i0]/P) + a2*np.cos(2*np.pi*ts[i0]/P) - a3*np.cos(3*np.pi*ts[i0]/P)
	
	#second part
	w[i2] = a0 - a1*np.cos(np.pi*(ts[i2] - (100 - 2*P))/P) + a2*np.cos(2*np.pi*(ts[i2] - (100 - 2*P))/P) - a3*np.cos(3*np.pi*(ts[i2] - (100 - 2*P))/P)
	
	out = w*v	

	return out
	
def _WFFlatTop(t,v,P=0.0):
	'''
	This will multiply the date by the flat top window.
	
	Inputs
	======
	t : float
		Time array
	v : float
		Time series to be windowed
	P : float
		Percentage of window to have the function applied at each end
		(if you set to 50.0, then the whole window has the function 
		applied to it).
		
	Returns
	=======
	out : float
		Windowed version of v
	
	'''
	#get the time range
	t0 = np.nanmin(t)
	t1 = np.nanmax(t)
	tr = t1 - t0
	
	#get a scaled time array
	ts = 100.0*(t - t0)/tr
	
	#work out the indices for each section
	i0 = np.where(ts < P)[0]
	i2 = np.where(ts > (100 - P))[0]

	#calculate the window function
	w = np.ones(t.size,dtype=v.dtype)

	#some constants
	a0 = 0.21557895
	a1 = 0.41663158
	a2 = 0.277263158
	a3 = 0.083578947
	a4 = 0.006947368

	
	#first part
	w[i0] = a0 - a1*np.cos(np.pi*ts[i0]/P) + a2*np.cos(2*np.pi*ts[i0]/P) - a3*np.cos(3*np.pi*ts[i0]/P) + a4*np.cos(4*np.pi*ts[i0]/P)
	
	#second part
	w[i2] = a0 - a1*np.cos(np.pi*(ts[i2] - (100 - 2*P))/P) + a2*np.cos(2*np.pi*(ts[i2] - (100 - 2*P))/P) - a3*np.cos(3*np.pi*(ts[i2] - (100 - 2*P))/P) + a4*np.cos(4*np.pi*(ts[i2] - (100 - 2*P))/P)
	
	out = w*v	

	return out
	
def _WFCosine(t,v,P=50.0):
	'''
	This will multiply the date by the sine window.
	
	Inputs
	======
	t : float
		Time array
	v : float
		Time series to be windowed
	P : float
		Percentage of window to have the function applied at each end
		(if you set to 50.0, then the whole window has the function 
		applied to it).
		
	Returns
	=======
	out : float
		Windowed version of v
	
	'''


	#get the time range
	t0 = np.nanmin(t)
	t1 = np.nanmax(t)
	tr = t1 - t0
	
	#get a scaled time array
	ts = 100.0*(t - t0)/tr
	
	#work out the indices for each section
	i0 = np.where(ts < P)[0]
	i2 = np.where(ts > (100 - P))[0]

	#calculate the window function
	w = np.ones(t.size,dtype=v.dtype)

	#first section
	w[i0] = np.sin(0.5*(ts[i0]/P)*np.pi)
	
	#last section
	w[i2] = np.sin(0.5*np.pi*(ts[i2] - (100 - 2*P))/P)
	
	
	#multiply by v
	out = v*w


	return out
	
def _WFGaussian(t,v,P=(0.5,50.0)):
	'''
	This will multiply the date by the Gaussian window.
	
	Inputs
	======
	t : float
		Time array
	v : float
		Time series to be windowed
	P : tuple
		p[0] : width of the Gaussian
		p[1] : Percentage of window to have the function applied at each 
		end	(if you set to 50.0, then the whole window has the function 
		applied to it).
		
	Returns
	=======
	out : float
		Windowed version of v
	
	'''


	#get the time range
	t0 = np.nanmin(t)
	t1 = np.nanmax(t)
	tr = t1 - t0
	
	#get a scaled time array
	ts = 100.0*(t - t0)/tr
	
	#work out the indices for each section
	i0 = np.where(ts < P[1])[0]
	#i1 = np.where((ts >= P) & (ts <= (100 - P)))[0]
	i2 = np.where(ts > (100 - P[1]))[0]

	#calculate the window function
	w = np.ones(t.size,dtype=v.dtype)

	#first section
	w[i0] = np.exp(-0.5*((ts[i0] - P[1])/(P[0]*P[1]))**2)
	
	#last section
	w[i2] = np.exp(-0.5*(((ts[i2] - (100 - P[1])))/(P[0]*P[1]))**2)
	
	
	#multiply by v
	out = v*w

	
	return out
