import numpy as np

def ApplyWindowFunction(v,WindowFunction=None,Param=None):

	# get the appropriate window function and parameters
	WF = {	'none':				(_WFNone,0.0),
			'cosine-bell':		(_WFCosineBell,10.0),
			'hamming':			(_WFHamming,10.0),
			'triangle':			(_WFTriangle,50.0),
			'welch':			(_WFWelch,0.0),
			'blackman':			(_WFBlackman,0.0),
			'nuttall':			(_WFNuttall,0.0),
			'blackman-nuttall':	(_WFBlackmanNuttall,0.0),
			'flat-top':			(_WFFlatTop,0.0),
			'cosine':			(_WFCosine,10.0),
			'gaussian':			(_WFGaussian,0.5)}
	
	Func,Pdef = WF.get(WindowFunction,(_WFNone,0.0))
	
	#check if anycustom parametersare being used
	if Param is None:
		P = Pdef
	else:
		P = Param
	
	#apply to data
	return Func(v,P)
	


def _WFNone(v,P):
	return v


def _WFCosineBell(v,P=10.0):
	
	l = np.size(v)
	
	t_len = np.int32(P*l/100.0)
	theta = (np.pi/np.float32(t_len))
	
	out = np.array(v)
	i = np.arange(t_len)
	out[i] = v[i]*(0.5+(0.5*np.cos(i*theta + np.pi)))
	
	i = np.arange(l-1-t_len,l)
	out[i] = v[i]*(0.5+(0.5*np.cos((l-1-i)*theta + np.pi)))
 
	return out


def _WFHamming(v,P=10.0):
	
	l = np.size(v)
	
	t_len = np.int32(P*l/100.0)
	theta = (np.pi/np.float32(t_len))
	
	out = np.array(v)
	i = np.arange(0,t_len)
	out[i] = v[i]*(0.46*np.cos(i*theta + np.pi))+0.54
	
	i = np.arange(l-1-t_len,l)
	out[i] = v[i]*(0.46*np.cos((l-1-i)*theta + np.pi))+0.54

	return out
	
def _WFTriangle(v,P=50.0):
	
	l = np.size(v)
	t_len = np.int32(P*l/100.0)
	
	out = np.array(v)
	
	i = np.arange(0,t_len)
	out[i] = (np.float32(i)/t_len)*v[i]
	i = np.arange(l-t_len,l)
	out[i] = ((l-np.float32(i))/t_len)*v[i]
	
	return out
	
def _WFWelch(v,P=0.0):

	l = np.float32(np.size(v))
	out = np.array(v)
	i = np.arange(0,l)
	out[i] = (1-((i-((l-1)/2))/((l-1)/2))**2)*v[i]
	
	return out
	
def _WFBlackman(v,P=0.0):
	
	l = np.size(v)
	out = np.array(v)
	
	a_0 = 7938.0/18608
	a_1 = 9240.0/18608
	a_2 = 1430.0/18608	
	i = np.arange(0,l)
	out[i] = (a_0 - a_1*np.cos((2*np.pi*i)/(l-1)) + a_2*np.cos((4*np.pi*i)/(l-1)))*v[i]

	return out
	
def _WFNuttall(v,P=0.0):
	
	l = np.size(v)
	out = np.array(v)	
	
	a_0 = 0.355768
	a_1 = 0.487396
	a_2 = 0.144232
	a_3 = 0.012604
	i = np.arange(0,l)
	out[i] = (a_0 - a_1*np.cos((2*np.pi*i)/(l-1))+a_2*np.cos((4*np.pi*i)/(l-1))-a_3*np.cos((6*np.pi*i)/(l-1)))*v[i]
	
	return out


def _WFBlackmanNuttall(v,P=0.0):

	l = np.size(v)
	out = np.array(v)	

	a_0 = 0.3635819
	a_1 = 0.4891775
	a_2 = 0.1365995
	a_3 = 0.0106411
	i = np.arange(0,l)
	out[i] = (a_0 - a_1*np.cos((2*np.pi*i)/(l-1))+a_2*np.cos((4*np.pi*i)/(l-1))-a_3*np.cos((6*np.pi*i)/(l-1)))*v[i]	
	
	return out
	
def _WFFlatTop(v,P=0.0):

	l = np.size(v)
	out = np.array(v)	

	a_0 = 1.0
	a_1 = 1.93
	a_2 = 1.29
	a_3 = 0.388
	a_4 = 0.028
	i = np.arange(0,l)
	out[i] = (a_0 - a_1*np.cos((2*np.pi*i)/(l-1))+a_2*np.cos((4*np.pi*i)/(l-1))-a_3*np.cos((6*np.pi*i)/(l-1))+a_4*np.cos((8*np.pi*i)/(l-1)))*v[i]	
	
	return out
	
def _WFCosine(v,P=10.0):
	
	l = np.size(v)
	out = np.array(v)	

	
	t_len = np.int32(P*l/100.0)
	theta = (np.pi/np.float32(t_len))
	
	out = np.array(v)
	i = np.arange(0,t_len)
	out[i] = v[i]*np.cos(i*theta/2 - np.pi/2)
	
	i = np.arange(l-1-t_len,l)
	out[i] = v[i]*np.cos((l-1-i)*theta/2 - np.pi/2)

	return out
	
def _WFGaussian(v,P=0.5):

	l = np.size(v)
	out = np.array(v)	
	i = np.arange(0,l)
	out[i] = v[i]*(np.exp(-0.5*(((i-(l-1)/2.0)/(P*(l-1)/2.0))**2)))
	
	return out
