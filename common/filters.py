"""
A script demonstrating the use of HIGH; LOW; BAND filters with Scipy.
Taken from http://mpastell.com/2010/01/18/fir-with-scipy/

Simple fft-based filtering:
http://www.swharden.com/blog/2009-01-21-signal-filtering-with-python/

Also interesting signal smoothing with windows:
http://www.scipy.org/Cookbook/SignalSmooth
FIR filtering:
http://www.scipy.org/Cookbook/FIRFilter

I got the Savitsky-Golay filter from: 
http://www.scipy.org/Cookbook/SavitzkyGolay

25.7. adding the weiner filter accessor method - from scipy

"""

import scipy.signal as signal
import scipy
import numpy as np
import matplotlib.pyplot as plt

run_examples = False 

def weiner(x,window_len=11):
    """
    Smoothing using the Weiner filter (implemented in SciPy):
    http://en.wikipedia.org/wiki/Wiener_filter
    """
    return signal.weiner(x, window_len)

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]

#Plot frequency and phase response
def mfreqz(b,a=1):
    w,h = signal.freqz(b,a)
    h_dB = 20 * np.log10(abs(h))
    plt.subplot(211)
    plt.plot(w/max(w),h_dB)
    plt.ylim(-150, 5)
    plt.ylabel('Magnitude (db)')
    plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    plt.title(r'Frequency response')
    plt.subplot(212)
    h_Phase = np.unwrap(np.arctan2(np.imag(h),np.real(h)))
    plt.plot(w/max(w),h_Phase)
    plt.ylabel('Phase (radians)')
    plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    plt.title(r'Phase response')
    plt.subplots_adjust(hspace=0.5)

#Plot step and impulse response
def impz(b,a=1):
    l = len(b)
    impulse = np.repeat(0.,l); impulse[0] =1.
    x = np.arange(0,l)
    response = signal.lfilter(b,a,impulse)
    plt.subplot(211)
    plt.stem(x, response)
    plt.ylabel('Amplitude')
    plt.xlabel(r'n (samples)')
    plt.title(r'Impulse response')
    plt.subplot(212)
    step = np.cumsum(response)
    plt.stem(x, step)
    plt.ylabel('Amplitude')
    plt.xlabel(r'n (samples)')
    plt.title(r'Step response')
    plt.subplots_adjust(hspace=0.5)

if run_examples:

    """ LOW PASS FILTER EXAMPLE """
    n = 61
    a = signal.firwin(n, cutoff = 0.3, window = "hamming")
    #Frequency and phase response
    mfreqz(a)
    plt.show()
    #Impulse and step response
    plt.figure(2)
    impz(a)
    plt.show()
    
    """ HIGH PASS FILTER EXAMPLE """
    n = 101
    a = signal.firwin(n, cutoff = 0.3, window = "hanning")
    #Spectral inversion
    a = -a
    a[n/2] = a[n/2] + 1
    plt.figure(3)
    mfreqz(a)
    plt.show()
    
    """ BAND-PASS FILTER EXAMPLE """
    n = 1001
    #Lowpass filter
    a = signal.firwin(n, cutoff = 0.3, window = 'blackmanharris')
    #Highpass filter with spectral inversion
    b = - signal.firwin(n, cutoff = 0.5, window = 'blackmanharris'); b[n/2] = b[n/2] + 1
    #Combine into a bandpass filter
    d = - (a+b); d[n/2] = d[n/2] + 1
    #Frequency response
    plt.figure(4)
    mfreqz(d)
    plt.show()

    # This script demonstrates how to use band-pass (low-pass)  
    # filtering to eliminate electrical noise and static  
    # from signal data!  
    
    ##################  
    ### PROCESSING ###  
    ##################  
    
    xs=np.arange(1,100,.01) #generate Xs (0.00,0.01,0.02,0.03,...,100.0)  
    signal = sin1=np.sin(xs*.3) #(A)  
    sin1=np.sin(xs) # (B) sin1  
    sin2=np.sin(xs*2.33)*.333 # (B) sin2  
    sin3=np.sin(xs*2.77)*.777 # (B) sin3  
    noise=sin1+sin2+sin3 # (C)  
    static = (np.random.random_sample((len(xs)))-.5)*.2 # (D)  
    sigstat=static+signal # (E)  
    rawsignal=sigstat+noise # (F)  
    fft=scipy.fft(rawsignal) # (G) and (H)  
    bp=fft[:]  
    for i in range(len(bp)): # (H-red)  
        if i>=10:bp[i]=0  
    ibp=scipy.ifft(bp) # (I), (J), (K) and (L)  
    
    ################  
    ### GRAPHING ###  
    ################  
    
    h,w=6,2  
    plt.figure(figsize=(12,9))  
    plt.subplots_adjust(hspace=.7)  
    
    plt.subplot(h,w,1);plt.title("(A) Original Signal")  
    plt.plot(xs,signal)  
    
    plt.subplot(h,w,3);plt.title("(B) Electrical Noise Sources (3 Sine Waves)")  
    plt.plot(xs,sin1,label="sin1")  
    plt.plot(xs,sin2,label="sin2")  
    plt.plot(xs,sin3,label="sin3")  
    plt.legend()  
    
    plt.subplot(h,w,5);plt.title("(C) Electrical Noise (3 sine waves added together)")  
    plt.plot(xs,noise)  
    
    plt.subplot(h,w,7);plt.title("(D) Static (random noise)")  
    plt.plot(xs,static)  
    plt.axis([None,None,-1,1])  
    
    plt.subplot(h,w,9);plt.title("(E) Signal + Static")  
    plt.plot(xs,sigstat)  
    
    plt.subplot(h,w,11);plt.title("(F) Recording (Signal + Static + Electrical Noise)")  
    plt.plot(xs,rawsignal)  
    
    plt.subplot(h,w,2);plt.title("(G) FFT of Recording")  
    fft=scipy.fft(rawsignal)  
    plt.plot(abs(fft))  
    plt.text(200,3000,"signals",verticalalignment='top')  
    plt.text(9500,3000,"static",verticalalignment='top',  
               horizontalalignment='right')  
    
    plt.subplot(h,w,4);plt.title("(H) Low-Pass FFT")  
    plt.plot(abs(fft))  
    plt.text(17,3000,"sin1",verticalalignment='top',horizontalalignment='left')  
    plt.text(37,2000,"sin2",verticalalignment='top',horizontalalignment='center')  
    plt.text(45,3000,"sin3",verticalalignment='top',horizontalalignment='left')  
    plt.text(6,3000,"signal",verticalalignment='top',horizontalalignment='left')  
    plt.axvspan(10,10000,fc='r',alpha='.5')  
    plt.axis([0,60,None,None])  
    
    plt.subplot(h,w,6);plt.title("(I) Inverse FFT")  
    plt.plot(ibp)  
    
    plt.subplot(h,w,8);plt.title("(J) Signal vs. iFFT")  
    plt.plot(signal,'k',label="signal",alpha=.5)  
    plt.plot(ibp,'b',label="ifft",alpha=.5)  
    
    plt.subplot(h,w,10);plt.title("(K) Normalized Signal vs. iFFT")  
    plt.plot(signal/max(signal),'k',label="signal",alpha=.5)  
    plt.plot(ibp/max(ibp),'b',label="ifft",alpha=.5)  
    
    plt.subplot(h,w,12);plt.title("(L) Difference / Error")  
    plt.plot(signal/max(signal)-ibp/max(ibp),'k')
    
    
    #savefig("SIG.png",dpi=200)  
    plt.show()  

def savitzky_golay(y, window_size, order, deriv=0):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techhniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m, y, mode='valid')
