import numpy as np
# from pylab import *
import common.sonpy2.son as son
import matplotlib.pyplot as plt
import scipy.signal as sp
# from scipy import *

DEBUG_EXTRACTION = False

def extract_voltage_sample(file_name, total_time=5):
    voltChNr = 2
    voltCh = son.Channel(voltChNr, file_name)
    volt = voltCh.data()
    volt_dt = voltCh.info.dt
    ind_max = len(volt) - 1
    
    if ind_max*volt_dt < total_time:
        assert False, "Voltage length is too small for the desired sample size."
    elif ind_max*volt_dt/2. < total_time:
        volt = volt[int(-total_time/volt_dt)]
    else:
        volt = volt[int(ind_max/2.):(int(ind_max/2.)+int(total_time/volt_dt))]
    
    return volt, volt_dt

def extract_steps(filename, setup, pre_ms, post_ms, step_ms=0, case="box"):

    markNr = 32

    if setup == 'DMD':
        currChNr = 1
        voltChNr = 2
    elif setup == 'SEC':
        currChNr = 6
        voltChNr = 7
    elif setup == 'current_injection':
        currChNr = 1
        voltChNr = 2
            
    else:
        print('please, specify setup correctly! ...aborting')
    
        
    currCh = son.Channel(currChNr, filename)
    curr = currCh.data()
    voltCh = son.Channel(voltChNr, filename)
    volt = voltCh.data()
    # get the sampling intervals (dt) for the current and voltage channels:
    curr_dt = currCh.info.dt
    volt_dt = voltCh.info.dt
    
    try:
        mark = son.Channel(markNr, filename)
        currMarker = mark.data(timeunits = 'ticks')

        #print zip(*currMarker)
        # since _marker.py zips data into a list, one has to use the unzip 
        # operator '*' to recover the initial data types:


        [MarkerTimings, hexMark] = zip(*currMarker) 
        #print hexMark

        # MarkerTimings come in tics. Current data are sampled each nth tic, where
        # n is the product of timePerADC and channel.divide. Division of MarkerTimings
        # by this number gives proper indices for data extraction:
        
        # The old version of software was storing the divide in ch.info.divide, the new one in ch.info.lChanDvd

        if case in ["combo", "exp", "combo2"]:
            MarkerTimings = MarkerTimings[::2]
            hexMark = hexMark[::2]
            
            hexMark=[int(x.encode("hex")) for x in hexMark]
            # Here I filter out all the nonstandard stimulus types (long pulses, bridge) 
            tmp = np.array(zip(hexMark,MarkerTimings))
            tmp = tmp[tmp[:,0]<6]
            MarkerTimings = tmp[:,1]
            
        try:
            MarkTCurr = MarkerTimings/(currCh.info.divide*currCh.fhead.timePerADC)
            # same for voltage data:
            MarkTVolt = MarkerTimings/(voltCh.info.divide*voltCh.fhead.timePerADC)

        except:
            MarkTCurr = MarkerTimings/(currCh.info.lChanDvd*currCh.fhead.timePerADC)
            # same for voltage data:
            MarkTVolt = MarkerTimings/(voltCh.info.lChanDvd*voltCh.fhead.timePerADC)

    except:  
        #applies if marker channel was not recorded
        raise
        print 'no marker channel found... trying to extract step onsets from current channel'
        if curr_dt == volt_dt:
            [b,a] = sp.butter(3,0.05) # define low-pass filter parameters
            filtCurr = sp.lfilter(b,a,curr)[200:] # ommit filter artifact at beginning of data
            CurrAmpl = max(filtCurr)-min(filtCurr)
            stepInd = np.nonzero(filtCurr<max(filtCurr)-(CurrAmpl/2))
            temp = stepInd[0][np.nonzero(np.diff(stepInd)>1)[1]+1]
            MarkTCurr = np.insert(temp,0,stepInd[0][0])
            MarkTVolt = MarkTCurr
            MarkerTimings = MarkTCurr*curr_dt
 
            # extract step length:
            temp2 = stepInd[0][np.nonzero(np.diff(stepInd)>1)[1]]
            stepEnd = np.append(temp2,stepInd[0][-1])
            steps = stepEnd - MarkTCurr
            step = round(np.mean(steps))
            step_ms = int(round(step*1000*curr_dt))
            print 'extracted {0} steps of {1}ms length'.format(len(steps),step_ms)
            
        else:
            print 'current and voltage smapling rates do not match! Aborting...'
            return

    # check if sampling interval is extracted correctly:
    currLength = round(currCh.info.dt*len(curr), 3)
    fileLength = round(currCh.fhead.maxFTime*currCh.fhead.usPerTime*currCh.fhead.dTimeBase,3)

    if currLength == fileLength:
        print('sampling rate correct')
    else:
        print('sampling rate does not match with max Time!')
        print'\n', 'data length = ', currLength, ', file length = ', fileLength
    
    pre = round((float(pre_ms)/1000)/curr_dt)
    post = round((float(post_ms)/1000)/curr_dt)

    #
    # if step length is not specified (step_ms set to zero), the following
    # code tries to extract the duration of the current step
    #

    if step_ms == 0:
        print 'step length not specified... trying to extract from current traces'
        search = round((500.0/1000)/curr_dt) # 500 ms interval for step search
        prelCurrents = np.empty(((pre  + post + search), len(MarkTCurr)))
        print len(curr), len(prelCurrents[:,0])
        for n, m in enumerate(MarkTCurr):
            prelCurrents[:,n] = curr[m - pre : m  + post + search]

        meanPrelCurr = np.mean(prelCurrents,1)

        currAmp = np.mean(meanPrelCurr[pre + round(0.001/curr_dt) : pre + round(0.004/curr_dt)])
        bline = np.mean(meanPrelCurr[0 : pre])
        halfAmp = currAmp - bline

        stepInd = np.nonzero(meanPrelCurr > int(bline+halfAmp))
        step = stepInd[0][-1] - stepInd[0][0]
        print 'extracted ', round((step*curr_dt)*1000), 'ms step length', '\n'

    else:
        step = round((float(step_ms)/1000)/curr_dt)
    
    Currents = np.empty(((pre + step + post), len(MarkTCurr)))
    n = 0
    for m in MarkTCurr:
        Currents[:,n] = curr[m - pre : m + step + post]
        n += 1

    Voltages = np.empty(((pre + step + post), len(MarkTVolt)))
    n = 0
    for m in MarkTVolt:
        Voltages[:,n] = volt[m - pre : m + step + post]
        n += 1

    if DEBUG_EXTRACTION:
        plt.plot(curr/10.)
        plt.plot(volt)
        plt.plot(MarkTVolt,np.ones(len(MarkTVolt))*(-40), "*")
        plt.show()
        
    # exract current and voltage averages:
    
    meanVolt = np.mean(Voltages,1)
    meanCurrent = np.mean(Currents,1)
       
    return MarkerTimings, Currents, Voltages, meanCurrent, meanVolt, curr_dt, volt_dt, step

def extractSpikesForRasters(filename, setup):

    # extracts spike times and timings of trigger pulses which mark the onset
    # of current injection periods as well as onset of current steps for extraction
    # of passive cell parameters, and the trigger pulses which mark the timings of
    # injected pulse packets

    # returns:
    # 

    # expects:
    # the filename of an .smr-file containing the raw data:
    # channel 1 with injected noise currents
    # channel 2 with the voltage responses
    # channel 3 with analog data with trigger pulses

    # trigger pulses should have different amplitudes:
    # 2 V for onset of stimulation periods
    # 4 V for onset of parameter test pulses
    # 5 V for pulse packets within the noise

    # respective stimulation files can be generated with the matlab-file
    # StationaryJMIP_spikeDyn.m

    # C. Boucsein 15.01.2013


    spikeThr = 0 # set spike detection threshold in mV
    trigThr = 1 # set trigger detection threshold in V
    

    if setup == 'DMD':
        voltChNr = 2
        trigChNr = 3
    elif setup == 'SEC':
        voltChNr = 7
        trigChNr = 8
    else:
        print('please, specify setup correctly! ...aborting')
    
    voltCh = son.Channel(voltChNr, filename)
    volt = voltCh.data()
    trigCh = son.Channel(trigChNr, filename)
    trig = trigCh.data()

    # get the sampling intervals (dt) for the current and voltage channels:
    volt_dt = voltCh.info.dt
    trig_dt = trigCh.info.dt

    #
    # detect spikes in voltage traces
    #

    temp = np.nonzero(volt>spikeThr) # find voltage values above threshold
    idx = np.insert(np.nonzero(np.diff(temp)[0]>1)[0]+1,0,0) # indices of first threshold crossings rel.
    # to indices above threhold
    spikeTtics = temp[0][idx] # extract indices of spikes in tics
    spikeTsec = spikeTtics * volt_dt
    
    #
    # detect trigger pulses in trigger channel
    #

    aboveT = np.nonzero(trig>trigThr) # find voltage values above threshold
    trigidx = np.insert(np.nonzero(np.diff(aboveT)[0]>1)[0]+1,0,0) # indices of first threshold crossings rel. to indices above threhold
    trigTtics = aboveT[0][trigidx] # extract indices of trigger events in tics
    trigTsec = trigTtics * trig_dt
    trigAmp = trig[trigTtics + 5]

    offsetTtics = trigTtics[np.nonzero(trigAmp < 2.5)]
    offsetTsec = offsetTtics * trig_dt

    temp5 = np.nonzero(trigAmp < 3.5)
    onsetTtics = trigTtics[temp5[0][np.nonzero(trigAmp[temp5]>2.5)]]

    print 'onsetTtics: ', onsetTtics
    print 'diff(onsetTtics): ', np.diff(onsetTtics)
    
    temp2 = np.nonzero(trigAmp < 4.5)
    paratrigTtics = trigTtics[temp2[0][np.nonzero(trigAmp[temp2]>3.5)]]

    temp3 = np.nonzero(trigAmp < 5.5)
    PPTtics = trigTtics[temp3[0][np.nonzero(trigAmp[temp3]>4.5)]]

    stimduration = (offsetTtics[0] - onsetTtics[0])*trig_dt

    data = {}
        
    data['spikeT'] = {}
    data['metaData'] = {}
    data['metaData']['onsetTtics'] = onsetTtics
    data['metaData']['offsetTtics'] = offsetTtics
    data['metaData']['paratrigTtics'] = paratrigTtics
    data['metaData']['PPTtics'] = PPTtics
    data['metaData']['volt_dt'] = volt_dt
    data['metaData']['trig_dt'] = trig_dt
    data['metaData']['stimduration'] = stimduration


    # for l in metaBlockLines:
        # if len(lines[l]) > 2:
            # data['run{0}'.format(k)]['metaData'][lines[l][1][0:-1]] = [lines[l][2]]

        # now fill spike time data

    trial = 0

    for m in onsetTtics:
        print 'trial nr. ', trial
        firstParaTic = paratrigTtics[np.nonzero(paratrigTtics > m)[0]][0]
        temp4 = np.nonzero(spikeTtics > m)
        spikeTimes = (spikeTtics[temp4][np.nonzero(spikeTtics[temp4]<firstParaTic)] - m) * volt_dt
        data['spikeT']['trial{0}'.format(trial)] = spikeTimes
        trial += 1

        if trial == 0:
            temp5 = np.nonzero(PPTtics > m)
            PPTimes = (PPTtics[temp5][np.nonzero(PPTtics[temp5]<firstParaTic)] - m) * trig_dt
            data['PPTimes'] = PPTimes


    print 'extracted {0} trials from data file {1}'.format(trial, filename)

    # conduct = float(data['run{0}'.format(k)]['metaData']['g'][0][0:-2])
    # ampl = float(data['run{0}'.format(k)]['metaData']['amplitude'][0][0:-2])
    # print 'conductance: {0}, amplitude: {1}'.format(conduct, ampl)



    return data

#
#
# plotting routines:
#
#


def plot_avCurrStepResp(filename, condition, meanVolt, fitData, Rin, tauFit_ms, DC, Cap, pre, step, volt_dt, ax0, ay0, dx, dy, xwidth, xswidth, ywidth, yswidth, plotPos):

    #
    # plot voltage traces and fit
    # 

    n,m=plotPos

    xAxis = np.linspace(1., len(meanVolt), len(meanVolt))
    axesPos = [ax0+(m*(dx+xwidth)), 1-ay0-n*(dy+ywidth), xwidth, ywidth]
    ax = plt.axes(axesPos)
    plt.plot(meanVolt, color = 'k')
    plt.plot(xAxis[pre:pre+fitData.size], fitData)

    ax.set_autoscale_on(False)
    # ax.set_ylim(np.min(values)*1.3, np.max(values)*1.3)
    # ax.set_ylim((wcCurr[pre + (0.1*pre):pre + (1.8*pre),:].min())*1.3, (wcCurr[pre + (0.1*pre):pre + (1.8*pre),:].max()))
    # ax.set_ylim((wcCurr[2*pre:3*pre,:].min())*1.5, (wcCurr[2*pre:3*pre,:].max())*1.3)
    # ax.set_xlim(pre*0.5, wcCurr.shape[0])
    ax.set_frame_on(False)
    ax.set_xticks([])
    ax.set_yticks([])

    #
    # plot scale bars
    #
    plt.hlines(ax.get_ylim()[0], pre*0.52, pre*0.52 + 0.1/volt_dt, colors = 'k', linestyles = 'solid', linewidth = 3)
    plt.text(pre*0.5, ax.get_ylim()[0] - (ax.get_ylim()[1]-ax.get_ylim()[0])*0.08,'100 ms', fontsize = 8,)
    plt.vlines(pre*0.52, ax.get_ylim()[0], ax.get_ylim()[0] + 1, colors = 'k', linestyles = 'solid', linewidth = 3)
    minYax = ax.get_ylim()[0]
    maxYax = ax.get_ylim()[1]
    plt.text(pre*0.1, minYax + 0.3,'1 mV', fontsize = 8, rotation = 'vertical')
    plt.text(pre+1.2*step, maxYax-(0.3*(maxYax-minYax)), 'Rin = {0}  MOhm'.format(Rin), fontsize = 12)
    plt.text(pre+1.2*step, maxYax-(0.4*(maxYax-minYax)), 'Tau = {0} ms'.format(tauFit_ms), fontsize = 12)
    plt.text(pre+1.2*step, maxYax-(0.5*(maxYax-minYax)), 'Cap = {0} pF'.format(Cap), fontsize = 12)
    plt.text(pre+1.2*step, maxYax-(0.6*(maxYax-minYax)), 'DC = {0} pA'.format(DC), fontsize = 12)
    plt.text(pre+step*0.5, ax.get_ylim()[1] + (ax.get_ylim()[1]-ax.get_ylim()[0])*0.03, 'file {0},{1} baseline mp: {2:.1f}mV'.format(filename, '\n', condition), fontsize = 10)
    # plt.text(pre*0.5, ax.get_ylim()[1] + (ax.get_ylim()[1]-ax.get_ylim()[0])*0.03, ('passive responses ' + conditions[n]), fontsize = 10)

    return

def plot_singleCurrStepResp(filename, condition, Voltages, meanVolt, fitData, Rin, tauFit_ms, plateau, DC, Cap, pre, step, volt_dt, ax0, ay0, dx, dy, xwidth, xswidth, ywidth, yswidth, plotPos):

    #
    # plot voltage traces and fit
    # 

    n,m=plotPos

    xAxis = np.linspace(1., len(meanVolt), len(meanVolt))
    axesPos = [ax0+(m*(dx+xwidth)), 1-ay0-n*(dy+ywidth), xwidth, ywidth]   
    ax = plt.axes(axesPos)

    for k in np.arange(np.shape(Voltages)[1]):
        plt.plot(xAxis,Voltages[:,k], color = '0.8')

    plt.plot(meanVolt, color = 'k', linewidth = 3)
    plt.plot(xAxis[pre+step:pre+step+fitData.size], fitData, color = 'g', linewidth = 3)

    ax.set_autoscale_on(False)
    plotSpace = (np.max(meanVolt) - np.min(meanVolt))*0.4
    ax.set_ylim(np.min(meanVolt)-plotSpace, np.max(meanVolt)+plotSpace)
    # ax.set_ylim((wcCurr[pre + (0.1*pre):pre + (1.8*pre),:].min())*1.3, (wcCurr[pre + (0.1*pre):pre + (1.8*pre),:].max()))
    # ax.set_ylim((wcCurr[2*pre:3*pre,:].min())*1.5, (wcCurr[2*pre:3*pre,:].max())*1.3)
    # ax.set_xlim(pre*0.5, wcCurr.shape[0])
    ax.set_frame_on(False)
    ax.set_xticks([])
    ax.set_yticks([])

    #
    # plot scale bars
    #

    # specify size of voltage and time scalebar:
    voltSclb = 5 # in mV
    timeSclb = 0.1 # in seconds

    xRange = ax.get_xlim()[1]
    yRange = ax.get_ylim()[1] - ax.get_ylim()[0]

    # scalebar coordinates (Lx: left x, Rx: right x, Ly: lower y, Uy: upper y):
    sclbLx = ax.get_xlim()[1]-(pre)
    sclbRx = sclbLx + timeSclb/volt_dt
    sclbLy = ax.get_ylim()[1]-yRange*0.5
    sclbUy = sclbLy + voltSclb

    # plt.hlines(ax.get_ylim()[0], pre*0.52, pre*0.52 + 0.1/volt_dt, colors = 'k', linestyles = 'solid', linewidth = 3)
    plt.hlines(sclbLy, sclbLx, sclbRx, colors = 'k', linestyles = 'solid', linewidth = 3)
    plt.text(sclbLx + 0.02/volt_dt, sclbLy - 0.1*yRange,'{0} ms'.format(timeSclb*1000), fontsize = 8,)
    plt.vlines(sclbLx, sclbLy, sclbUy, colors = 'k', linestyles = 'solid', linewidth = 3)
    # minYax = ax.get_ylim()[0]
    # maxYax = ax.get_ylim()[1]
    plt.text(sclbLx - 0.06*xRange, sclbLy + 0.3,'{0} mV'.format(voltSclb), fontsize = 8, rotation = 'vertical')
    # plt.text(pre+1.2*step, maxYax-(0.3*(maxYax-minYax)), 'Rin = {0}  MOhm'.format(Rin), fontsize = 12)
    # plt.text(pre+1.2*step, maxYax-(0.4*(maxYax-minYax)), 'Tau = {0} ms'.format(tauFit_ms), fontsize = 12)
    # plt.text(pre+1.2*step, maxYax-(0.5*(maxYax-minYax)), 'Cap = {0} pF'.format(Cap), fontsize = 12)
    # plt.text(pre+1.2*step, maxYax-(0.6*(maxYax-minYax)), 'DC = {0} pA'.format(DC), fontsize = 12)
    plt.text(pre+step*0.5, ax.get_ylim()[1] - (ax.get_ylim()[1]-ax.get_ylim()[0])*0.0, 'file {0},{1} baseline mp: {2:.1f}mV {3} plateuiness: {4}'.format(filename, '\n', condition, '\n', plateau), fontsize = 10)
    # plt.text(pre*0.5, ax.get_ylim()[1] + (ax.get_ylim()[1]-ax.get_ylim()[0])*0.03, ('passive responses ' + conditions[n]), fontsize = 10)

    return

def downSampleWithPeak(yValues, sampleFactor):

    #
    # routine for down sampling of time series data containing one short transient event
    # which otherwise does not make it to the display because of aliasing.
    #
    # expects the time series data as numpy array, where first dimension has length of time
    # series and second dimension is number of trials, and the desired factor for down sampling
    #
    # returns the down-sampled time series, and x-axis vectors for each of them
    # to enable proper plotting
    #

    indVector = np.arange(0, yValues.shape[0]-sampleFactor, sampleFactor)

    # generate empty arrays for y and x values:

    downsYValues = np.zeros((len(indVector),yValues.shape[1]))
    xAxis = np.zeros((len(indVector),yValues.shape[1]), int)

    # loop over traces and fill empty arrays for downsampled y and x values

    for m in np.arange(yValues.shape[1]):
        peakInd = np.nonzero(yValues[:,m] == np.max(yValues[:,m]))[0][0]
        shift = np.mod(peakInd, sampleFactor)
        xAxis[:,m] = indVector + shift
        downsYValues[:,m] = yValues[:,m][xAxis[:,m]]

    return downsYValues, xAxis


#
#
# fitting routines:
#
#

def fit_cell_parameters(meanVolt, meanCurrent, volt_dt, pre_tics, step_tics, fitInt_tics):

    pre_tics = pre_tics + 180
    bline = np.mean(meanVolt[0:pre_tics])
    Vampl = min(meanVolt) - bline

    #
    # optimize fit interval: find minimum and take 75% of the time to reach mininmum
    # as the fit interval (due to Ih-activation, it is misleading to take the time
    # of the minimum itself)
    #

    minTime=np.nonzero(meanVolt==meanVolt.min(0))
    fitInt_tics=(int(minTime[0][0])-pre_tics)*0.75
    fitInt_ms=int(round(fitInt_tics*volt_dt*1000))
    print('fit interval: {0} ms'.format(fitInt_ms))
    xfitInt = np.linspace(1., fitInt_tics, fitInt_tics)
    stepInt = np.linspace(1., step_tics, step_tics)
    tau_ms = 20
    tau_tics = round((float(tau_ms)/1000)/volt_dt)

    # Ufit = bline - (-Vampl*(1-np.exp(-fitInt/tau_tics)))
    
    #
    # Fit mono-exponential function to voltage response:
    #

    fitfunc = lambda p, x: bline - (-p[0]*(1-np.exp(-x/p[1]))) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    p0 = [Vampl, tau_tics]  # initial values for fit
    p1, success = sp.optimize.leastsq(errfunc, p0[:], args=(xfitInt, meanVolt[pre_tics:pre_tics+fitInt_tics]))
    fitData = fitfunc(p1,stepInt)

    voltAmp = (p1[0]*(1-np.exp(-step_tics/p1[1])))

    #
    # calculate current amplitude and Rin
    #

    currbline = np.mean(meanCurrent[0:pre_tics])
    currMax = np.mean(meanCurrent[pre_tics:pre_tics+step_tics])
    currAmp = (currbline-currMax)/1000
    # Rin = int(-round(voltAmp/currAmp))
    Rin = -round(voltAmp/currAmp)
    tauFit_tics = p1[1]
    # tauFit_ms = int(round(tauFit_tics*volt_dt*1000))
    tauFit_ms = tauFit_tics*volt_dt*1000
    DC = int(round(currbline))
    Cap = int(round((tauFit_ms/Rin)*1000))
    Rin = int(Rin)
    tauFit_ms = int(tauFit_ms)
    # print(Cap)

    return fitData, Rin, tauFit_ms, DC, Cap


def fitExpRise(xValues, yValues, startGuess):

    # 
    # startGuess needs to be array of the four parameters of the fit function
    # A*exp(k*x+n)+m in that order
    # where
    # A: y-axis intercept
    # k: kinkyness
    # n: x-axis offset
    # m: y-axis offset
    # 

    #
    # first sort xValues in ascending order and arrange yValues accordingly
    #

    xValuesSort = np.sort(xValues)
    ind = xValues.argsort()
    yValuesSort = yValues[ind]

    #
    # if startGuess is empty, start parameters are estimated from data:
    #

    print 'yValuesSort: ', yValuesSort

    if bool(startGuess) == False:
        print 'trying to guess start parameters...'
        A = np.max(yValues)
        k = 0.05
        # n = xValues[(yValues>(b+0.4*(a/2-b))) & (yValues<(b+0.6*(a/2-b)))][0]
        n = xValues[-1]
        m = np.min(yValues)
        startGuess = [A,k,n,m]
        print 'start guess: ', startGuess

    p0 = startGuess  # initial values for fit

    #
    # Fit mono-exponential rise function to yValues:
    #

    fitfunc = lambda p, x: p[0]*np.exp(p[1]*x-p[2])+p[3] # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function

    p1, success = sp.optimize.leastsq(errfunc, p0[:], args=(xValuesSort, yValuesSort))

    # 
    # calculate a graph of the fitted function of length xValues but with 
    # small sample interval for proper plotting
    #

    fitSteps = np.linspace(xValuesSort[0], xValuesSort[-1], 10000)

    fitData = fitfunc(p1,fitSteps)

    return fitData, p1


def fitBoltzRise(xValues, yValues, startGuess):

    # 
    # startGuess needs to be array of the four parameters of the fit function 
    # (here Boltzmann-function) in the order [a, b, k, x0]
    #  f(x) = b + (a-b)/(1+exp(k*(x-x0))) where:
    # a: upper limit of function values
    # b: lower limit of function values
    # k: steepness of sigmoidal rise (if k > 0 -> falling sigmoidal,
    #    if k < 0 -> rising sigmoidal
    # x0: x-axis location where half-maximal value of Boltzmann-function
    # is reached
    #

    #
    # first sort xValues in ascending order and arrange yValues accordingly
    #

    xValuesSort = np.sort(xValues)
    ind = xValues.argsort()
    yValuesSort = yValues[ind]

    #
    # if startGuess is empty, start parameters are estimated from data:
    #

    if bool(startGuess) == False:
        print 'trying to guess start parameters...'
        a = 2*np.max(yValues)
        b = np.min(yValues)
        k = -1
        try:
            x0 = xValues[(yValues>(b+0.4*(a/2-b))) & (yValues<(b+0.6*(a/2-b)))][0]
        except:
            x0 = xValues[-1]
        startGuess = [a,b,k,x0]
        print 'start guess: ', startGuess

    p0 = startGuess  # initial values for fit

    # startGuess = [1000,0,-1.05,-60]

    #
    # Fit mono-exponential rise function to yValues:
    #

    # fitfunc = lambda p, x: p[0]*np.exp(p[1]*x+p[2])+p[3] # Target function
    # fitfunc = lambda p, x: p[0]/(np.exp(-p[2]*x+p[3]))+p[1] # Target function
    # fitfunc = lambda p, x: p[0]/(np.exp(p[1]*x+p[2])) # Target function
    fitfunc = lambda p, x: p[1]+(p[0]-p[1])/(1+np.exp(p[2]*(x-p[3]))) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function

    p1, success = sp.optimize.leastsq(errfunc, p0[:], args=(xValuesSort, yValuesSort))

    # 
    # calculate a graph of the fitted function of length xValues but with 
    # small sample interval for proper plotting
    #

    fitSteps = np.linspace(xValuesSort[0], xValuesSort[-1], 10000)

    fitData = fitfunc(p1,fitSteps)

    return fitData, p1

