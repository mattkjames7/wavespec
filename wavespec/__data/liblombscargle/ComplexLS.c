/*This is the C function in which I will try to write my function for the 
complex Lomb-Scargle */
#include "ComplexLS.h"
/* Samuel Wharton's code!!!*/

void ComplexLS(int lent, int lenw, float *time, float *data, float *omega, float *amplitude, float *phase, float *a, float *b)
{
	/*This will calculate the amplitude and phase of the Lomb-Scargle Periodogram
	
	Parameters
	----------
	int lent - length of the time array. 
	int lenw - length of the omega array. 
	float *time - the time array. 
	float *data - the data/amplitude array. 
	float *omega - the angular frequency array. 
	float *amplitude - the LS amplitude array to be filled. 
	float *phase - the LS phase array to be filled. 
	
	*/
	
	int i;
	
	for (i=0; i<lenw; i++){
		
		int t;
		double ss = 0.0, sc = 0.0, costop = 0.0, cosbottom = 0.0;
		double sintop = 0.0, sinbottom = 0.0, tau;
		
		/*Calculate the sum of twosin and twocos*/
		for (t=0; t<lent; t++){
			ss = ss + sin(2.0*omega[i]*time[t]);
			sc = sc + cos(2.0*omega[i]*time[t]); 
			}
		
		tau = (atan2(ss,sc))/(2*omega[i]); 
			
		/*Now calculate the component values*/
		for (t=0; t<lent; t++){
			costop = costop + data[t]*cos(omega[i]*(time[t]-tau));
			cosbottom = cosbottom + pow(cos(omega[i]*(time[t]-tau)),2);
			sintop = sintop + data[t]*sin(omega[i]*(time[t]-tau));
			sinbottom = sinbottom + pow(sin(omega[i]*(time[t]-tau)),2);	
			}

		/*Now calculate a and b*/
		a[i] = (pow(2.0/lent,0.5)*costop)/(pow(cosbottom,0.5));
		b[i] = (pow(2.0/lent,0.5)*sintop)/(pow(sinbottom,0.5));
		
		/*Now work out the phase and amplitude. */
		phase[i] = -atan2(b[i],a[i]) -omega[i]*tau;  
		if (phase[i] > M_PI){phase[i] -= 2*M_PI;}
		if (phase[i] <= -M_PI){phase[i] += 2*M_PI;}
		
		/* I modified this bit so that the amplitude matches FFT */
		//amplitude[i] = (lent/2.0)*pow(pow(a[i],2) + pow(b[i],2),0.5);
		amplitude[i] = pow(pow(a[i],2) + pow(b[i],2),0.5);
		
		}
}
