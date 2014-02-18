import numpy, pylab, os, shelve, zshelve
# import PyVision.src.h5shelve as h5shelve
# import PyVision.components.sonpy.son as son
import sonpy.son as son
from numpy import nan
import Tkinter, tkFileDialog
from threading import Thread
import tempfile
import scipy.io as io
from NeuroTools.parameters import ParameterRange
from NeuroTools.parameters import ParameterSpace
from NeuroTools.sandbox import make_name
import scipy.signal as signal
from scipy.signal import lfilter,buttord,butter
import shutil, mdp
import PyVision.src.analysis.receptive_field_mapping_core as rf
from mpl_toolkits.mplot3d import Axes3D
import gc
import scipy.ndimage as nd
import figures

def sort_during_go(onsets,offsets,times,data):
    during_go = []
    
    take_go = numpy.zeros(len(data),dtype=int)
    take_no_go = numpy.ones(len(data),dtype=int)
    
    for index in range(len(onsets)):
        take = (times>=onsets[index])&(times<offsets[index])
        take2 = numpy.where(take==True)
        take_go[take2] = 1
        take_no_go[take2] = 0
        
    return take_go.astype(bool),take_no_go.astype(bool)
        
def zscore(data):
    return (data-data.mean())/data.std()

def blip_finder(self,threshold):
    corner_freq = 4.
    stop_freq = 6.
    data = self.filter(corner_freq,stop_freq)
    data = zscore(data)
    data2 = numpy.where(data>=threshold,1.,0.)
    data_diff = numpy.diff(data2)
    
    up = numpy.where(data_diff==1)[0]
    down = numpy.where(data_diff==-1)[0]
    
    up_time = (up/self.sampling_rate)
    down_time = (down/self.sampling_rate)
    
    width = down-up
    width =  (width/self.sampling_rate)*1000.
    amplitude = numpy.array([data[up[i]:down[i]].max() for i in range(len(up)-1)])
    
    return width[0:-1], amplitude,up_time[0:-1],down_time[0:-1]


def blip_finder_go_no_go(file,self,threshold):
    w,a,up,down = blip_finder(self,threshold)
    take_go, take_no_go = sort_during_go(file.onsets,file.offsets,up,w)
    return w[take_go],w[take_no_go],a[take_go],a[take_no_go]

def spike_sort_dialog(spike_sort_function,channel,threshold,window,max_clusters,frac):
    # are we happy?
    input = raw_input('Are you happy with the sorting? default is yes, if not write: no ')
    if input == 'no':
        pylab.close('all')
        input = raw_input('Please specify the cut out window, old was %f,%f '%(window[0],window[1]))
        if not len(input) == 0:
            input = input.split(',')
            window = numpy.array(input).astype(float)
            
        input = raw_input('Please give threshold, old was %f '%threshold)
        if not len(input) == 0:
            threshold = float(input)
        
        input = raw_input('Max clusters, old was %d '%max_clusters)
        if not len(input) == 0:
            max_clusters = int(input)

        input = raw_input('PCA frac, old was %d '%frac)
        if not len(input) == 0:
            frac = float(input)

        #input = raw_input('PCA components, old was %d '%components)
        #if not len(input) == 0:
        #    components = int(input)
        # we restart the spike sorting
        spike_sort_function(channels=[channel],threshold=threshold,window=window,pca=True,threshold_detection=True,plot=True,max_clusters=max_clusters,frac=frac)
    else:
       input = raw_input('You are happy, great. Please provide a list of unit id which you considere to be valid neurons ')
       if not len(input) == 0:
           input = input.split(',')
           valid_units = numpy.array(input).astype(int)
       else:
           valid_units = []
       return valid_units
        
        
def savetxt(filename, data, format, delimiter):
    """
    Due to the lack of savetxt in older versions of numpy
    we provide a cut-down version of that function.
    """
    f = open(filename, 'w')
    f.write('%d\n'%data.shape[1])
    for row in data:
        f.write(delimiter.join([format%val for val in row]) + '\n')
    f.close()

def pca_spike_waves(spike_waves,components=2,frac=0.):
    Pcomponents,  Trans, fracVar = pylab.prepca(spike_waves,frac=frac)
    spike_waves_pca = Pcomponents#[0:components,:]
    return spike_waves_pca



def savebinary(filename,spike_waves):
    from scipy.io.numpyio import fwrite, fread

    fd = open(filename, 'wb')
    fwrite(fd, spike_waves.size, spike_waves)
    fd.close()
    

def Klusters(spike_waves,timestamps,channel=1,tempdir='',max_clusters=3,frac=1.0):
    # pca
    spike_waves_pca = mdp.pca(spike_waves.T,**{'output_dim':frac}).T
    
    dims_new_data = numpy.array(spike_waves_pca.shape)
    dims_new_data[0] += 1
    new_data = numpy.zeros(dims_new_data)
    new_data[0:-1,:] = spike_waves
    new_data[-1,:] = timestamps
    
    # write the data to the KlustaKwik
    klustakwik_filename = tempdir+'klustakwik_tmp'
    savetxt(klustakwik_filename+'.fet.%d'%channel,spike_waves_pca.T,'%s',' ')
    savebinary(klustakwik_filename+'.spk.%d'%channel,spike_waves.T)
    print "Spike files are in : ",tempdir
    raw_input("Please press any button when sorting is done.")
    
def KlustaKwik(spike_waves,timestamps=None,channel=1,pca=True,tempdir='',max_clusters=3,frac=1.0):
    if not spike_waves.shape[1] >1:
        return
    # print('KlustaKwik of channel %d'%self.info.chan)
    #if pca:
    #    spike_waves_pca  = pca_spike_waves(spike_waves,components=components)
    #    data = spike_waves_pca
    #else:
    data = spike_waves
    
    # peak trough height ratio
    #data_norm = data.copy()
    #data_norm = data_norm/numpy.abs(data_norm).max()
    #trough = data_norm.min(axis=0)
    #height = data_norm.max(axis=0)
    #pth_ratio = trough/height
    
    print('PCA')
    spike_waves_pca = mdp.pca(data.T,**{'output_dim':frac}).T
    if spike_waves_pca.shape[1] == 1:
        print "We only have once dim,we need at least 2"
    data = spike_waves_pca
    if timestamps is not None:
        dims_new_data = numpy.array(data.shape)
        dims_new_data[0] += 1
        new_data = numpy.zeros(dims_new_data)
        new_data[0:-1,:] = data
        new_data[-1,:] = timestamps

        print "With timestamps"
        print "Dim of data: ",data.shape
        print "Dim of data with timestamps: ",new_data.shape
    # 
    #from scipy.stats import shapiro
    
    #features = data.shape[0]
    #w_values = numpy.zeros(features)
    #for feature in range(features):
    #    w_values[feature] = shapiro(data[feature,:])[0]
    #    i_select = numpy.argsort(w_values)[0:components]

    #data = data[i_select,:]
    #print data.shape
    # data = mdp.fastica(data.T,dtype='float32').T
    
    #    print data.shape
    #    spike_waves_pca = data
    
    # combine data
    #data = numpy.append(pth_ratio.reshape((1,pth_ratio.shape[0])),data,axis=0)
    
    # write the data to the KlustaKwik
    klustakwik_filename = tempdir+'klustakwik_tmp'
    savetxt(klustakwik_filename+'.fet.%d'%channel,data.T,'%s',' ')
    # run
    line = 'KlustaKwik %s %d -MaxPossibleClusters %d -UseFeatures all -Screen 0'%(klustakwik_filename,channel,max_clusters)
    print(line)
    os.system(line)
    # load clusters
    clusters = numpy.loadtxt(klustakwik_filename+'.clu.%d'%channel,skiprows=1)
    print "Done loading clusters"
    return clusters,spike_waves_pca

class MKK(Thread):
    def __init__(self,klustakwik_filename,channel,max_clusters):
        Thread.__init__(self)
        self.klustakwik_filename = klustakwik_filename
        self.channel = channel
        self.max_clusters = max_clusters
        self.clusters = None
        
    def run(self):
        # run
        line = 'KlustaKwik %s %d -MaxPossibleClusters %d -UseFeatures all -Screen 0'%(self.klustakwik_filename,self.channel,self.max_clusters)
        print(line)
        os.system(line)
        # load clusters
        self.clusters = numpy.loadtxt(self.klustakwik_filename+'.clu.%d'%self.channel,skiprows=1)
        
def KlustaKwik_MultiParameter(spike_waves,channel=1,pca=True,tempdir='',max_clusters=[3,6,9],components=2):

    
    if not spike_waves.shape[1] >1:
        return
    # print('KlustaKwik of channel %d'%self.info.chan)
    print('a')
    #if pca:
    #    spike_waves_pca  = pca_spike_waves(spike_waves,components=components)
    #    if 1: return [],[]
    #    data = spike_waves_pca
    #else:
    data = spike_waves
    #    spike_waves_pca = data
    print('b')
    # write the data to the KlustaKwik
    klustakwik_filename = 'klustakwik_tmp'
    klustakwik_filename_fet = klustakwik_filename+'.fet.%d'%channel
    savetxt(tempdir+klustakwik_filename_fet,data.T,'%s',' ')
    print('c')
    threads = []
    for max_cluster in max_clusters:
        klustakwik_filename_tmp = 'max_cluster_%d_%s'%(max_cluster,klustakwik_filename)
        shutil.copyfile(tempdir+klustakwik_filename_fet,tempdir+klustakwik_filename_tmp+'.fet.%d'%channel)
        thread = MKK(tempdir+klustakwik_filename_tmp,channel,max_cluster)
        thread.start()
        threads.append(thread)

    clusters = {}
    for thread in threads:
        thread.join()
        clusters[thread.max_clusters] = thread.clusters
        
    return clusters

class ThresholdDetectionThread(Thread):
    def __init__(self,file,channel,threshold,window):
        Thread.__init__(self)
        self.file = file
        self.channel = channel
        self.threshodl= threshold
        self.window = window
        
    def run(self):
        self.file.threshold_detection(channels=[self.channel],threshold=self.threshold,window=self.window)

class FileGroup(object):
    def __init__(self,filenames=None,stimulus_format='VE'):
        self.files = {}
        self.stimulus_format = stimulus_format
        self.tempdir = tempfile.mkdtemp()+'/'
        if filenames == None:
            print('Please select the smr files')
            filenames = tkFileDialog.askopenfilenames()
        self.filenames = filenames
        for index, filename in enumerate(filenames):
            self.files[index] = File(filename,stimulus_format=stimulus_format)
    
    def set_valid_units(self,channel,units=[]):
        for file in self.files.values():
            file.set_valid_units(channel,units=units)
        
    def plot_preview(self,xlim=[1,5.]):
        for file in self.files.values():
            print(file)
            file.plot_preview(xlim=xlim)
    
    def plot_psth_across_experiments(self,channels,binsize,units=None,onsets=None,valid_units=False):
        for file in self.files.values():
            file.plot_psth_across_experiments(channels,binsize,units=units,onsets=onsets,valid_units=valid_units)
        
    def chan_list(self):
        for file in self.files.values():
            print(file)
            file.chan_list()
            
    def save_sorted_data(self):
        for file in self.files.values():
            print('Save file: '+file.name)
            file.save_sorted_data()
            
    def load_sorted_data(self,load_raw=True):
        for file in self.files.values():
            print('Load file: '+file.name)
            file.load_sorted_data(load_raw=load_raw)
        
    
    def plot_spike_clusters(self,channels=None):
        for file in self.files.values():
            file.plot_spike_clusters(channels=channels)
    
    def merge_clusters(self,channel,clusters,new_cluster_id=None):
        for file in self.files.values():
            file.merge_clusters(channel,clusters,new_cluster_id)
            
    def spike_sorting(self,channels=None,threshold=5.,window=[0.5,1.5],pca=True,threshold_detection=True,plot=True,max_clusters=3,frac=1.0,minimal_isi=0.0,with_timestamps=False):
        
        gc.collect()
        for channel in channels:
            spike_waves_dict = {}
            spike_times_dict = {}
            timestamps_dict = {}
            keys = []
            number_of_spikes = {}
            for key, file in self.files.items():
                keys.append(key)
                print('File: '+file.name)
                file.check_if_channels_are_loaded([channel])

            # threshold detection
            if threshold_detection:
                thread_list = []
                for key, file in self.files.items():
                    # first we need to do a threshold detection
                    # thread = ThresholdDetectionThread(file,channel,threshold,window)
                    thread = Thread(target=file.threshold_detection(channels=[channel],threshold=threshold,window=window,minimal_isi=minimal_isi))
                    thread.start()
                    thread_list.append(thread)
                for thread in thread_list:
                    thread.join()
                
            for key, file in self.files.items():
                spike_waves_dict[key] = file.channels[channel].spike_waves
                spike_times_dict[key] = file.channels[channel].spike_times
                number_of_spikes[key] = file.channels[channel].spike_waves.shape[1]
                dims = int(numpy.sum(file.channels[channel].window_tick))

                if with_timestamps:
                    timestamps_dict[key] = file.channels[channel].spike_ticks
            
            spike_waves_all_files = numpy.empty((dims,0)) 
            spike_times_all_files = numpy.empty((0))
            timestamps_all_files = numpy.empty((0))
            for key in keys:
                spike_waves_all_files = numpy.append(spike_waves_all_files,spike_waves_dict[key],axis=1)
                spike_times_all_files = numpy.append(spike_times_all_files,spike_times_dict[key])
                if with_timestamps:
                    timestamps_all_files = numpy.append(timestamps_all_files,timestamps_dict[key])
                    
            if with_timestamps:
                timestamps = timestamps_all_files
                
            else:
                timestamps = None
            print('KlustaKwik')
            clusters,spike_waves_pca = KlustaKwik(spike_waves_all_files,channel=channel,pca=pca,tempdir=self.tempdir,max_clusters=max_clusters,frac=frac,timestamps=timestamps)
            
            
            if plot:
                unique_clusters = numpy.unique(clusters)
                # isi
                pylab.figure()
                pylab.get_current_fig_manager().window.wm_geometry("+0+0")
                
                x = numpy.ceil(numpy.sqrt(len(unique_clusters)))
                y = numpy.ceil(numpy.sqrt(len(unique_clusters)))
                for index, cluster in enumerate(unique_clusters):
                    pylab.subplot(x,y,index+1)
                    take = clusters == cluster
                    spikes = spike_times_all_files[take]*1000.# to get it into ms
                    isi = numpy.diff(spikes)
                    isi_dist = numpy.histogram(isi,bins=numpy.arange(0,100.,1.))
                    pylab.plot(isi_dist[1][0:-1],isi_dist[0])
                    pylab.title('C: %d u: %d #s %d'%(channel,cluster,spikes.shape[0]))
                # spike waves
                pylab.figure()
                pylab.get_current_fig_manager().window.wm_geometry("+650+0")
                ylim = (spike_waves_all_files.min(),spike_waves_all_files.max())
                    
                for index, cluster in enumerate(unique_clusters):
                    pylab.subplot(x,y,index+1)
                    take = clusters == cluster
                    data = spike_waves_all_files[:,take]
                    pylab.plot(self.files[0].channels[channel].get_spike_waves_times(data=data),data)
                    pylab.title('C: %d u: %d #s %d'%(channel,cluster,data.shape[1]))
                    pylab.ylim(ylim)
                    
                # clusters
                fig = pylab.figure()
                pylab.get_current_fig_manager().window.wm_geometry("+1300+0")
                #if components >2:
                #    ax = Axes3D(fig)
                for cluster in numpy.unique(clusters):
                    take = clusters == cluster
                    #if components >2:
                    #    ax.scatter(spike_waves_pca[0,take],spike_waves_pca[1,take],spike_waves_pca[2,take])
                    #else:
                    pylab.plot(spike_waves_pca[0,take],spike_waves_pca[1,take],'.',label='unit %d, #spikes %d'%(cluster,numpy.sum(take)))
                pylab.legend(numpoints=1)
                pylab.title('Channel: %d'%(channel))

            #valid_units =spike_sort_dialog(self.spike_sorting,channel,threshold,window,max_clusters,components)
            #if valid_units is not None:
            #    for key, file in self.files.items():
            #        file.channels[channel].valid_units = valid_units
                
            # are we happy?
            input = raw_input('Are you happy with the sorting? default is no, write yes to leave ')
            if not input == 'yes':
                threshold_detection = False
                #pylab.close('all')
                input = raw_input('Please specify the cut out window, old was %f,%f '%(window[0],window[1]))
                if not len(input) == 0:
                    input = input.split(',')
                    window = numpy.array(input).astype(float)
                    threshold_detection = True
                input = raw_input('Please give threshold, old was %f '%threshold)
                if not len(input) == 0:
                    threshold = float(input)
                    threshold_detection = True
                input = raw_input('Max clusters, old was %d '%max_clusters)
                if not len(input) == 0:
                    max_clusters = int(input)
                
                input = raw_input('PCA frac, old was %f '%frac)
                if not len(input) == 0:
                    frac = float(input)
                pylab.close('all')
                self.spike_sorting(channels=[channel],threshold=threshold,window=window,pca=pca,threshold_detection=threshold_detection,plot=True,max_clusters=max_clusters,frac=frac)
                return None
            elif input == 'yes':
                input = raw_input("Great sorting is done, Please give the units which you think are valid spikes: format: 1,2,4:  ")
                if not len(input) == 0:
                    valid_units = numpy.array(input.split(",")).astype('int')
                else:
                    valid_units = []
                
                #input = raw_input('Great sorting down. Do you want to merge clusters manually? '))pylab.close('all')
                #if input == 'yes':
                #    input = raw_input(''))
                
                    
            count = 0
            # it is very important that we iterate over the "keys" as this will ensure that the order is kept
            for key in keys:
                file = self.files[key]
                clusters_tmp = clusters[count:count+number_of_spikes[key]]
                file.channels[channel].spike_waves_pca = spike_waves_pca[:,count:count+number_of_spikes[key]]
                file.channels[channel].valid_units = valid_units
                count += number_of_spikes[key]
                file.channels[channel].clusters = clusters_tmp
                for cluster in numpy.unique(clusters_tmp):
                    take = clusters_tmp==cluster
                    file.channels[channel].spike_waves_clusters[cluster] = file.channels[channel].spike_waves[:,take]
                    file.channels[channel].spike_times_clusters[cluster] = file.channels[channel].spike_times[take]
                    if not numpy.sum(take) == 0:
                        file.channels[channel].spike_waves_pca_clusters[cluster] = file.channels[channel].spike_waves_pca[:,take]
                    else:
                        file.channels[channel].spike_waves_pca_clusters[cluster] = numpy.empty((file.channels[channel].spike_waves_pca.shape[0],0))

    def spike_sorting_klusters(self,channels=None,threshold=5.,window=[0.5,1.5],pca=True,threshold_detection=True,plot=True,max_clusters=3,frac=1.0,minimal_isi=0.0,with_timestamps=False):
        
        gc.collect()
        for channel in channels:
            spike_waves_dict = {}
            spike_times_dict = {}
            timestamps_dict = {}
            keys = []
            number_of_spikes = {}
            for key, file in self.files.items():
                keys.append(key)
                print('File: '+file.name)
                file.check_if_channels_are_loaded([channel])

            # threshold detection
            if threshold_detection:
                thread_list = []
                for key, file in self.files.items():
                    # first we need to do a threshold detection
                    # thread = ThresholdDetectionThread(file,channel,threshold,window)
                    thread = Thread(target=file.threshold_detection(channels=[channel],threshold=threshold,window=window,minimal_isi=minimal_isi))
                    thread.start()
                    thread_list.append(thread)
                for thread in thread_list:
                    thread.join()
                
            for key, file in self.files.items():
                spike_waves_dict[key] = file.channels[channel].spike_waves
                spike_times_dict[key] = file.channels[channel].spike_times
                number_of_spikes[key] = file.channels[channel].spike_waves.shape[1]
                dims = int(numpy.sum(file.channels[channel].window_tick))

                if with_timestamps:
                    timestamps_dict[key] = file.channels[channel].spike_ticks
            
            spike_waves_all_files = numpy.empty((dims,0)) 
            spike_times_all_files = numpy.empty((0))
            timestamps_all_files = numpy.empty((0))
            for key in keys:
                spike_waves_all_files = numpy.append(spike_waves_all_files,spike_waves_dict[key],axis=1)
                spike_times_all_files = numpy.append(spike_times_all_files,spike_times_dict[key])
                if with_timestamps:
                    timestamps_all_files = numpy.append(timestamps_all_files,timestamps_dict[key])
                    
            if with_timestamps:
                timestamps = timestamps_all_files
                
            else:
                timestamps = None
            print('Klusters')
            Klusters(spike_waves_all_files,timestamps,channel=channel,tempdir=self.tempdir)
            
def jkljlk():
            if plot:
                unique_clusters = numpy.unique(clusters)
                # isi
                pylab.figure()
                pylab.get_current_fig_manager().window.wm_geometry("+0+0")
                
                x = numpy.ceil(numpy.sqrt(len(unique_clusters)))
                y = numpy.ceil(numpy.sqrt(len(unique_clusters)))
                for index, cluster in enumerate(unique_clusters):
                    pylab.subplot(x,y,index+1)
                    take = clusters == cluster
                    spikes = spike_times_all_files[take]*1000.# to get it into ms
                    isi = numpy.diff(spikes)
                    isi_dist = numpy.histogram(isi,bins=numpy.arange(0,100.,1.))
                    pylab.plot(isi_dist[1][0:-1],isi_dist[0])
                    pylab.title('C: %d u: %d #s %d'%(channel,cluster,spikes.shape[0]))
                # spike waves
                pylab.figure()
                pylab.get_current_fig_manager().window.wm_geometry("+650+0")
                ylim = (spike_waves_all_files.min(),spike_waves_all_files.max())
                    
                for index, cluster in enumerate(unique_clusters):
                    pylab.subplot(x,y,index+1)
                    take = clusters == cluster
                    data = spike_waves_all_files[:,take]
                    pylab.plot(self.files[0].channels[channel].get_spike_waves_times(data=data),data)
                    pylab.title('C: %d u: %d #s %d'%(channel,cluster,data.shape[1]))
                    pylab.ylim(ylim)
                    
                # clusters
                fig = pylab.figure()
                pylab.get_current_fig_manager().window.wm_geometry("+1300+0")
                #if components >2:
                #    ax = Axes3D(fig)
                for cluster in numpy.unique(clusters):
                    take = clusters == cluster
                    #if components >2:
                    #    ax.scatter(spike_waves_pca[0,take],spike_waves_pca[1,take],spike_waves_pca[2,take])
                    #else:
                    pylab.plot(spike_waves_pca[0,take],spike_waves_pca[1,take],'.',label='unit %d, #spikes %d'%(cluster,numpy.sum(take)))
                pylab.legend(numpoints=1)
                pylab.title('Channel: %d'%(channel))

            #valid_units =spike_sort_dialog(self.spike_sorting,channel,threshold,window,max_clusters,components)
            #if valid_units is not None:
            #    for key, file in self.files.items():
            #        file.channels[channel].valid_units = valid_units
                
            # are we happy?
            input = raw_input('Are you happy with the sorting? default is no, write yes to leave ')
            if not input == 'yes':
                threshold_detection = False
                #pylab.close('all')
                input = raw_input('Please specify the cut out window, old was %f,%f '%(window[0],window[1]))
                if not len(input) == 0:
                    input = input.split(',')
                    window = numpy.array(input).astype(float)
                    threshold_detection = True
                input = raw_input('Please give threshold, old was %f '%threshold)
                if not len(input) == 0:
                    threshold = float(input)
                    threshold_detection = True
                input = raw_input('Max clusters, old was %d '%max_clusters)
                if not len(input) == 0:
                    max_clusters = int(input)
                
                input = raw_input('PCA frac, old was %f '%frac)
                if not len(input) == 0:
                    frac = float(input)
                pylab.close('all')
                self.spike_sorting(channels=[channel],threshold=threshold,window=window,pca=pca,threshold_detection=threshold_detection,plot=True,max_clusters=max_clusters,frac=frac)
                return None
            elif input == 'yes':
                input = raw_input("Great sorting is done, Please give the units which you think are valid spikes: format: 1,2,4:  ")
                if not len(input) == 0:
                    valid_units = numpy.array(input.split(",")).astype('int')
                else:
                    valid_units = []
                
                #input = raw_input('Great sorting down. Do you want to merge clusters manually? '))pylab.close('all')
                #if input == 'yes':
                #    input = raw_input(''))
                
                    
            count = 0
            # it is very important that we iterate over the "keys" as this will ensure that the order is kept
            for key in keys:
                file = self.files[key]
                clusters_tmp = clusters[count:count+number_of_spikes[key]]
                file.channels[channel].spike_waves_pca = spike_waves_pca[:,count:count+number_of_spikes[key]]
                file.channels[channel].valid_units = valid_units
                count += number_of_spikes[key]
                file.channels[channel].clusters = clusters_tmp
                for cluster in numpy.unique(clusters_tmp):
                    take = clusters_tmp==cluster
                    file.channels[channel].spike_waves_clusters[cluster] = file.channels[channel].spike_waves[:,take]
                    file.channels[channel].spike_times_clusters[cluster] = file.channels[channel].spike_times[take]
                    if not numpy.sum(take) == 0:
                        file.channels[channel].spike_waves_pca_clusters[cluster] = file.channels[channel].spike_waves_pca[:,take]
                    else:
                        file.channels[channel].spike_waves_pca_clusters[cluster] = numpy.empty((file.channels[channel].spike_waves_pca.shape[0],0))

class File(son.File):
    
    
    def __init__(self,filename=None,cpus=1,stimulus_format='VE'):
        if filename == None:
            print('Please select the smr file')
            filename = tkFileDialog.askopenfilename()
        son.File.__init__(self,filename)
        self.channels = {}
        self.cpus = cpus
        self.stimulus_format = stimulus_format
        self.tempdir = tempfile.mkdtemp()+'/'
        self._have_stimulus_parameters = False
        # self.onsets
        self.file_attributes = ['parameter_space','parameter_ranges','parameter_space_indexes','_have_stimulus_parameters','parameters_filename','presented_frames','frame_times','frame_duration','pre_post_stimulus_time','stimulus_time','stimulus_refreshrate','onsets','offsets','p']
        
        # self.path = self.name[0:self.name.rfind('/')+1]
        self.path = os.path.dirname(self.name)+'/'
        # self.exp_name = self.name[self.name.rfind('/')+1:-4]
        self.exp_name = os.path.basename(self.name)[0:-4]
        self.figure_path = self.path+'figures/'
        if not os.path.exists(self.figure_path): os.mkdir(self.figure_path)
            
    def set_valid_units(self,channel,units=[]):
        self.channels[channel].set_valid_units(units) 
    
    def merge_clusters(self,channel,clusters,new_cluster_id=None):
        self.channels[channel].merge_clusters(clusters,new_cluster_id=new_cluster_id)
    
    def assign_sorted_WaveMarks_to_channels(self,WaveMark_Channels=[]):
        
        for combination in WaveMark_Channels:
            wavemark_channel_number = combination[0]
            adc_channel_number = combination[1]
            self.get_channels([wavemark_channel_number,adc_channel_number])
            
            wm_channel = self.channels[wavemark_channel_number]
            adc_channel = self.channels[adc_channel_number]
            
            adc_channel.spike_times =  wm_channel.spike_times
            adc_channel.clusters =  wm_channel.clusters
            adc_channel.spike_waves =  wm_channel.spike_waves
            adc_channel.spike_waves_time =  wm_channel.spike_waves_time
            adc_channel.sort_clustered_data()
            
    def calculate_stimulus_onset_offset(self):

        #if self.stimulus_format == 'VE':
        stimulus_onset_frame = int(self.pre_post_stimulus_time*self.stimulus_refreshrate)
        stimulus_offset_frame = int((self.pre_post_stimulus_time+self.stimulus_time)*self.stimulus_refreshrate)
        stimulus_onset_times = self.trigger_on_stimulus_frame(stimulus_onset_frame)
        stimulus_offset_times = self.trigger_on_stimulus_frame(stimulus_offset_frame)
        #else:
            
        return stimulus_onset_times,stimulus_offset_times
        
    def trigger_on_stimulus_frame(self,frame=0):

        frame_onset_times = []
        for index in range(self.onsets.shape[0]):
            onset = self.onsets[index]
            offset = self.offsets[index]
            frame_times_tmp = self.frame_times[(self.frame_times>=onset)&(self.frame_times<offset)]
            frame_onset_times.append(frame_times_tmp[frame])

        return numpy.array(frame_onset_times)
            
    def plot_preview(self,xlim=[1,5.]):

        adc_channels = self.header.return_adc_chan_list()
        print adc_channels
        self.get_channels(channels=adc_channels,type=None)
        self.plot_channels(channels=adc_channels,xlim=xlim)
        
    def get_frame_duration(self):
        self.frame_times
        self.frame_duration = numpy.diff(self.frame_times).mean()*1000.
    frame_duration = property(get_frame_duration)
    
    def select_stimulus_parameters(self):
        print('Please select the stimulus parameter file ')
        filename = tkFileDialog.askopenfilename()
        self.parameters_filename = filename
        if filename.endswith('shelve') and self.stimulus_format=='VE':
            s = shelve.open(filename)
            self.parameter_space = s['space']
            self.parameter_ranges = s['ranges']
            self.parameter_space_indexes = {}
            experiments = range(len(s['data'].keys()))
            for experiment_number in experiments:
                index = s['data'][str(experiment_number)]['index']
                self.parameter_space_indexes[experiment_number] = numpy.array(index)
                # self.parameter_space_indexes.append(index)
            self._have_stimulus_parameters = True
            if s['data']['0'].has_key('presented_frames'):
                self.presented_frames = s['data']['0']['presented_frames']
            self.pre_post_stimulus_time = s['data']['0']['p']['stimuli']['pre_post_stimulus_time']/1000. # to get it to s
            self.stimulus_time = s['data']['0']['p']['stimuli']['stimulus_time']/1000. #
            self.stimulus_refreshrate = s['data']['0']['p']['stimuli']['refreshrate']
            self.p = s['data']['0']['p']
            s.close()
        elif filename.endswith('mat') and self.stimulus_format == 'PT':
            d = io.loadmat(filename)
            monitor_refreshrate = float(d['frameRate'][0][0])
            if d.has_key('TRframes'):
                self.frame_duration = (1./monitor_refreshrate)*d['TRframes'][0][0]
                self.stimulus_refreshrate = monitor_refreshrate/d['TRframes'][0][0]
                number_of_stimulus_frames = d['moviedurationframe'][0][0]/d['TRframes'][0][0]
                self.frame_times = (numpy.arange(0,number_of_stimulus_frames,1)*(self.frame_duration))+d['movieduration1'][0][0]
                
            if d.has_key('TRframe'):
                self.frame_duration = (1./monitor_refreshrate)*d['TRframe'][0][0]
                self.stimulus_refreshrate = monitor_refreshrate/d['TRframe'][0][0]
                # frame times should be in sec
                number_of_stimulus_frames = d['movieDurationFrames2'][0][0]/d['TRframe'][0][0]
                self.frame_times = (numpy.arange(0,number_of_stimulus_frames,1)*(self.frame_duration))+d['movieDurationSecs1'][0][0]

            if d.has_key('Rnoiseimg'):
                self.presented_frames = d['Rnoiseimg']
            if d.has_key('movieduration1'):
                self.pre_post_stimulus_time = float(d['movieduration1'][0][0])
            elif d.has_key('movieDurationSecs1'):
                self.pre_post_stimulus_time = float(d['movieDurationSecs1'][0][0])
            if d.has_key('movieduration'):
                self.stimulus_time = float(d['movieduration'][0][0])
            elif d.has_key('movieDurationSecs2'):
                self.stimulus_time = float(d['movieDurationSecs2'][0][0])
            
            
            self._have_stimulus_parameters = True
            
    def get_stimulus_onsets_offsets(self,channel=None):
        '''
        '''
        if self.stimulus_format == 'VE':
            if channel == None:
                channel = 34
            self.get_channels(channels=[channel])
            self.channels[channel].value
            self.onsets = self.channels[channel].value[0::2]
            self.offsets = self.channels[channel].value[1::2]
        elif self.stimulus_format == 'PT':
            if channel == None:
                channel = 4
            self.get_channels(channels=[channel],type='StimulusChannelPT')
            self.onsets = self.channels[channel].onsets
            self.offsets = self.channels[channel].offsets
    onsets = property(get_stimulus_onsets_offsets)
        
    def get_frame_times(self,channel=35):
        if self.stimulus_format == 'VE':
            self.get_channels(channels=[channel])
            self.channels[channel].value
            self.frame_times = self.channels[channel].value
        elif self.stimulus_format == 'PT':
            print('Not possible, because frame times are not stored with PT')
    frame_times = property(get_frame_times)
    
    def save_sorted_data(self):
        filename = self.name[0:-4]+'_sorted_spikes'
        #s = shelve.open(filename+'.shelve')
        data_file = zshelve.open(filename+'.zshelve')
        self.save_name = filename+'.zshelve'
        #data.writeback = True
        data = {}
        if not data.has_key('channels'):
            data['channels'] = {}
            data['file'] = {}
#        if self._have_stimulus_parameters:
            # to_save = ['parameter_space','parameter_ranges','parameter_space_indexes','_have_stimulus_parameters','parameters_filename','presented_frames','frame_times','frame_duration']
        for attribute in self.file_attributes:
            if hasattr(self,attribute):
                data['file'][attribute] = eval('self.%s'%attribute)
            # data['parameter_space'] = self.parameter_space
            # data['parameter_ranges'] = self.parameter_ranges
            # data['parameter_space_indexes'] = self.parameter_space_indexes
            # data['_have_stimulus_parameters'] = True
            #data['parameters_filename'] = self.parameters_filename
            #if hasattr(self,'presented_frames'):
            #    data['presented_frames'] = self.presented_frames
                
        for channel_key, channel in self.channels.items():
            print('save channel ',channel_key)
            data['channels'][channel_key] = channel.as_dict()
        #if format == 'shelve':
            #s = shelve.open(filename+'.shelve')
            #s['file'].update(data['file'])
            #s['channels'].update(data['channels'])
        data_file.update(data)
        data_file.sync()
        data_file.close()
        #else:
        #    print('Sorry so far only shelve is supported')
        #io.savemat(filename+'.mat',data,do_compression=True)
        #return data
    
    def load_sorted_data(self,load_raw=True):
        filename = self.name[0:-4]+'_sorted_spikes'
        s = zshelve.open(filename+'.zshelve')
        
        for attribute in self.file_attributes:
            if s['file'].has_key(attribute):
                print('Loading '+attribute)
                print 'load_raw ',load_raw
                exec("self.%s = s['file'][attribute]"%(attribute))
        
        # channels
        channels = numpy.array(s['channels'].keys()).astype(int)
        self.check_if_channels_are_loaded(channels,load_raw=load_raw)
        for channel_key in channels:
            channel = self.channels[channel_key]
            channel.set_from_dict(s['channels'][channel_key])
            if channel._thresholds_detected:
                channel.sort_clustered_data()
        s.close()

    def check_if_channels_are_loaded(self,channel_numbers,load_raw=True):
        for channel_number in channel_numbers:
            if not channel_number in self.channels.keys():
                self.get_channels(channels=[channel_number],load_raw=load_raw)
                
    def select_channels(self):
        self.chan_list()
        input = raw_input('Please give channel numbers which you would like to load (use , as seperator). ')
        input = input.split(',')
        channel_numbers = numpy.array(input).astype(int)
        self.check_if_channels_are_loaded(channel_numbers)
                    
        return channel_numbers
                      
    def get_channels(self,channels=None,type=None,load_raw=True):
        if channels==None:
            channels = self.select_channels()
        
        for channel in channels:
            print('Loading channel %d'%channel)
            if type == None:
                type_to_load = son._ChannelInfo(self.header,channel).type()
                if type_to_load == 'AdcMark':
                    chan = AdcMarkChannel(channel,self.name)
                elif type_to_load == 'EventRise':
                    chan = EventRiseChannel(channel,self.name)
                else:
                    chan = Channel(channel,self.name,load_raw=load_raw)
            else:
                chan = eval('%s(channel,self.name)'%type)
            chan.tempdir = self.tempdir
            self.channels[channel] = chan
    
    def plot_channels(self,channels,xlim=None):
        if channels==None:
            channels = select_channels()
        self.check_if_channels_are_loaded(channels)
        for channel_key in channels:
            self.channels[channel_key].plot(xlim=xlim)
            pylab.plot(self.onsets,numpy.ones(len(self.onsets)),'og')
            pylab.plot(self.offsets,numpy.ones(len(self.offsets)),'or')
    
    def plot_spike_clusters(self,channels):
        for channel_key in channels:
            if self.channels.has_key(channel_key):
                channel = self.channels[channel_key]
                channel.plot_spike_clusters()

    def sort_spikes_to_experiments(self,spike_times,onsets,offsets):
        if not self._have_stimulus_parameters:
            self.select_stimulus_parameters()
            
        trial_index = self.parameter_space[1].index('trials')
        dims = self.parameter_space[0]

        spike_times_experiments = numpy.empty(dims,dtype=object)
        #spike_times_experiments_pre_stimulus = numpy.empty(dims,dtype=object)
        #spike_times_experiments_stimulus = numpy.empty(dims,dtype=object)
        #spike_times_experiments_post_stimulus = numpy.empty(dims,dtype=object)
        
        experiments = len(self.parameter_space_indexes.keys())

        #stimulus_onsets, stimulus_offsets = self.calculate_stimulus_onset_offset()
        
        for experiment_number in range(experiments):
            # select the spike times
            take = (spike_times>=onsets[experiment_number])&(spike_times<offsets[experiment_number])
            #take_pre_stimulus = (spike_times>=onsets[experiment_number])&(spike_times<stimulus_onsets[experiment_number])
            #take_stimulus = (spike_times>=stimulus_onsets[experiment_number])&(spike_times<stimulus_offsets[experiment_number])
            #take_post_stimulus = (spike_times>=stimulus_offsets[experiment_number])&(spike_times<offsets[experiment_number])
            
            
            index = list(self.parameter_space_indexes[experiment_number])
            #spike_times_experiments_pre_stimulus[tuple(index)] = spike_times[take_pre_stimulus]-onsets[experiment_number]
            #spike_times_experiments_stimulus[tuple(index)] = spike_times[take_stimulus]-onsets[experiment_number]
            #spike_times_experiments_post_stimulus[tuple(index)] = spike_times[take_post_stimulus]-onsets[experiment_number]
            
            spike_times_experiments[tuple(index)] = spike_times[take]-onsets[experiment_number]
            
            #spike_times = spike_times[non_take]
        return spike_times_experiments#_pre_stimulus,spike_times_experiments_stimulus,spike_times_experiments_post_stimulus
    
    def remove_trials_dimension(self,data):
        data_tmp = numpy.array([])
        for i in range(data.shape[0]):
            data_tmp = numpy.append(data_tmp,data[i])
        return data_tmp
        
    def plot_psth_across_experiments(self,channels,binsize,units=None,onsets=None,valid_units=False,plot_spikes=True,sigma=0.,return_psths=False):
        return figures.plot_psth_across_experiments(self,channels,binsize,units=units,onsets=onsets,valid_units=valid_units,plot_spikes=plot_spikes,sigma=sigma,return_psths=return_psths)
     
    def plot_1d_tuning(self,channels,binsize,units=None,onsets=None,valid_units=False,plot_spikes=True):
        return figures.plot_1d_tuning(self,channels=channels,binsize=binsize,units=units,onsets=onsets,valid_units=valid_units,plot_spikes=plot_spikes)


   
    def _old(self):
        '''
        Binsize in seconds
        so far only 2d
        '''
        
        if onsets == None:
            self.onsets
            onsets = self.onsets
        
        trial_index = self.parameter_space[1].index('trials')
        duration = (self.offsets-onsets).mean()
        
        bins = numpy.arange(0.,duration,binsize) 
        for channel_key in channels:
            channel = self.channels[channel_key]
            if units == None:
                if len(channel.valid_units) <1:
                    units_tmp = channel.spike_times_clusters.keys()
                else:
                    units_tmp = channel.valid_units
            else:
                units_tmp = units
            for unit in units_tmp:
                if not channel.spike_times_clusters.has_key(unit):
                    continue
                spike_times = channel.spike_times_clusters[unit]
                spike_times_sorted = self.sort_spikes_to_experiments(spike_times,onsets,self.offsets)
                # spike_times_pre_stimulus,spike_times_stimulus,spike_times_post_stimulus = self.sort_spikes_to_experiments(spike_times,onsets,self.offsets)
                
                experiments = {}
                for key,value in  self.parameter_ranges.items():
                    experiments[key.replace('.','_')] = ParameterRange(value)
                experiments.pop('trials')
                experiments = ParameterSpace(experiments)
                psths = {}
                dims,labels = experiments.parameter_space_dimension_labels()
                ranges = experiments.get_ranges_values()
                #fr_pre_stimulus = numpy.empty(dims)
                #fr_stimulus = numpy.empty(dims)
                #fr_post_stimulus = numpy.empty(dims)
                for experiment in experiments.iter_inner():
                    index = experiments.parameter_space_index(experiment)
                    name = make_name(experiment,experiments.range_keys())
                    pylab.figure()
                    all_spike_times = self.remove_trials_dimension(spike_times_sorted[index])
                    #all_spike_times_pre_stimulus = self.remove_trials_dimension(spike_times_pre_stimulus[index])
                    #all_spike_times_stimulus = self.remove_trials_dimension(spike_times_stimulus[index])
                    #all_spike_times_post_stimulus = self.remove_trials_dimension(spike_times_post_stimulus[index])

                    psth = numpy.histogram(all_spike_times,bins=bins)
                    psths[name] = psth
                    #psth_pre_stimulus = numpy.histogram(all_spike_times_pre_stimulus,bins=bins)
                    #psth_stimulus = numpy.histogram(all_spike_times_stimulus,bins=bins)
                    #psth_post_stimulus = numpy.histogram(all_spike_times_post_stimulus,bins=bins)
                    # psths[name] = psth
                    #take_stimulus = (all_spike_times>=self.pre_post_stimulus_time)&(all_spike_times<(self.pre_post_stimulus_time+self.stimulus_time))
                    #take_baseline_pre = (all_spike_times<self.stimulus_time)
                    #take_baseline_post = (all_spike_times>=(self.pre_post_stimulus_time+self.stimulus_time))
                    
                    #fr[index] = len(all_spike_times[take_stimulus])
                    #fr_baseline_pre[index] = len(all_spike_times[take_baseline_pre]) 
                    #fr_baseline_post[index] = len(all_spike_times[take_baseline_post])
                    #fr_pre_stimulus[index] = len(all_spike_times_pre_stimulus)
                    #fr_stimulus[index] = len(all_spike_times_stimulus)
                    #fr_post_stimulus[index] = len(all_spike_times_post_stimulus)
                    #psth_min = psth[0].min()
                    #psth_max = psth[0].max()
                    psth_y = psth[0]
                    pylab.fill_between(psth[1][0:-1:1],psth_y,edgecolor='gray',color='gray')
                    #pylab.fill_between(psth_pre_stimulus[1][0:-1:1],psth_pre_stimulus[0],edgecolor='gray',color='gray')
                    #pylab.fill_between(psth_stimulus[1][0:-1:1],psth_stimulus[0],edgecolor='red',color='red')
                    #pylab.fill_between(psth_post_stimulus[1][0:-1:1],psth_post_stimulus[0],edgecolor='gray',color='gray')
                    
                    pylab.ylabel('spike count in bin (%g s)'%binsize,color='red')
                    pylab.xlabel('Time (s)')
                    pylab.yticks(color='red')
                    pylab.xlim(-0.1,duration+0.1) # +- 100 ms
                    #pylab.ylim(psth_min,psth_max+0.1*psth_max)
                    pylab.twinx()
                    
                    trials = spike_times_sorted[index].shape[0]
                    for trial in range(spike_times_sorted[index].shape[0]):
                        spike_times_tmp = spike_times_sorted[index][trial]
                        pylab.plot(spike_times_tmp,numpy.ones(len(spike_times_tmp))*trial,'.',color='green')
                        # pre
                        #spike_times_tmp = spike_times_pre_stimulus[index][trial]
                        #pylab.plot(spike_times_tmp,numpy.ones(len(spike_times_tmp))*trial,'.',color='green')
                        # stimulus
                        #spike_times_tmp = spike_times_stimulus[index][trial]
                        #pylab.plot(spike_times_tmp,numpy.ones(len(spike_times_tmp))*trial,'.',color='black')
                        # post stimulus
                        #spike_times_tmp = spike_times_post_stimulus[index][trial]
                        #pylab.plot(spike_times_tmp,numpy.ones(len(spike_times_tmp))*trial,'.',color='blue')
                        
                    if trial == 0:
                        pylab.ylabel('Trial')
                    pylab.yticks(color='black')
                    pylab.ylim(0+-0.1*trials,trials+0.1*trials)
                    pylab.xlim(-0.1,duration+0.1) # +- 100 ms
                    pylab.title('Channel: %d unit: %d'%(channel.info.chan,unit))
                    pylab.figtext(0.5,0.02,name,fontsize=10.,ha='center')
                # all psths
                pylab.figure()
                psth_max = 0
                
                for name, psth in psths.items():
                    if psth[0].max() > psth_max:
                        psth_max = psth[0].max()
                    pylab.plot(psth[1][0:-1:1],psth[0],label=name)
                pylab.xlim(-0.1,duration+0.1) # +- 100 ms
                pylab.ylim(-5,psth_max+0.1*psth_max)
                pylab.legend()
                
                # firing rate
                #pylab.figure()
                #pylab.plot(fr_stimulus.squeeze(),'r',label='stimulus')
                #pylab.plot(fr_pre_stimulus.squeeze(),color='0',label='baseline pre')
                #pylab.plot(fr_post_stimulus.squeeze(),color='0.5',label='baseline post')
                #index = numpy.where(numpy.array(dims)>1)[0][0]
                #label = labels[index]
                #xtickslabels = ranges[label]
                #xticks = pylab.xticks()[0]
                #pylab.xticks(xticks,xtickslabels)
                #pylab.xlabel(label)
                #pylab.ylabel('spike count')
                #pylab.legend()
 
    def sort_spikes_to_trials(self,spike_times,onsets,offsets):
        #self.onsets
        trials = onsets.shape[0]
        spike_times_trials = {}
        for trial in range(trials):
            take = (spike_times>=onsets[trial])&(spike_times<offsets[trial])
            non_take = take == False
            spike_times_trials[trial] = spike_times[take]-onsets[trial]
            spike_times = spike_times[non_take]
        return spike_times_trials
    
    def flatten_dict(self,data):
        return_data = numpy.array([])
        for value in data.values():
            return_data = numpy.append(return_data,value)
        return return_data.flatten()
    
    def plot_psth_across_trials(self,channels,binsize,onsets=None,units=None):
        '''
        Binsize in seconds
        '''
        if onsets == None:
            onsets = self.onsets
            
        onsets = numpy.array(onsets)-0.5 # 0.5 seconds before
        offsets = numpy.array(self.offsets) + 0.5
        #offsets = self.offsets
        trials = onsets.shape[0]
        duration = (offsets-onsets).mean()
        print duration
        
        bins = numpy.arange(0.,duration,binsize) 
        for channel_key in channels:
            channel = self.channels[channel_key]
            if units == None:
                units_tmp = channel.spike_times_clusters.keys()
            else:
                units_tmp = units
            for unit in units_tmp:
                if not channel.spike_times_clusters.has_key(unit):
                    continue
                spike_times = channel.spike_times_clusters[unit]
                
                spike_times_experiments = self.sort_spikes_to_trials(spike_times,onsets,offsets)
                
                pylab.figure()
                all_spike_times = self.flatten_dict(spike_times_experiments)
                all_spike_times.sort()
                psth = numpy.histogram(all_spike_times,bins=bins)
                psth_min = psth[0].min()
                psth_max = psth[0].max()
                psth_x = psth[1][0:-1:1].repeat(2)
                psth_x[1::2] += binsize
                psth_y = psth[0].repeat(2)
                pylab.fill_between(psth_x,psth_y,edgecolor='red',color='red')
                #pylab.step(psth[1][0:-1:1],psth[0],color='black',where='post')
                pylab.ylabel('spike count in bin (%g s)'%binsize,color='red')
                pylab.xlabel('Time (s)')
                pylab.yticks(color='red')
                pylab.xlim(-0.1,duration+0.1) # +- 100 ms
                pylab.ylim(psth_min,psth_max+0.1*psth_max)
                pylab.twinx()
                
                for trial in range(trials):
                    spike_times_tmp = spike_times_experiments[trial]
                    pylab.plot(spike_times_tmp,numpy.ones(len(spike_times_tmp))*trial,'.',color='black')
                    if trial == 0:
                        pylab.ylabel('Trial')
                pylab.yticks(color='black')
                pylab.ylim(0+-0.1*trials,trials+0.1*trials)
                pylab.xlim(-0.1,duration+0.1) # +- 100 ms
                pylab.title('Channel: %d unit: %d #spikes %d'%(channel.info.chan,unit,all_spike_times.shape[0]))
                
            
    def threshold_detection(self,channels=None,threshold=5.,window=[0.5,1.5],minimal_isi=0.0):
        channel_keys = channels
        if channel_keys is None:
            channel_keys = self.channels.keys()
        else:
            self.check_if_channels_are_loaded(channels)
        # if no channels were selected
        if len(channel_keys) == 0:
            self.get_channels()
            channel_keys = self.channels.keys()
        
        for channel_key in channel_keys:
            channel = self.channels[channel_key]
            channel.spike_waves_clusters = {}
            channel.spike_times_clusters = {}
            channel.spike_waves_pca_clusters = {}
            channel.threshold_detection(threshold=threshold,window=window,minimal_isi=minimal_isi)

    def calculate_receptive_field_properties(self,channel=None,unit=None,lags=numpy.arange(-4.,-1,1,dtype='int'),normalize_spike_times_to_stimulus_onset=False,sigma=(1.,1.,1.),f_values=[0.,1.],mean_value=0.5,save_fig=False,prefix=''):
        return figures.calculate_receptive_field_properties(self,channel=channel,unit=unit,lags=lags,normalize_spike_times_to_stimulus_onset=normalize_spike_times_to_stimulus_onset,sigma=sigma,f_values=f_values,mean_value=mean_value,save_fig=save_fig,prefix=prefix)
    
    def show_reverse_correlation(self,channel=None,unit=None,lags=numpy.arange(-4.,-1,1,dtype='int'),lim_data=True,normalize_spike_times_to_stimulus_onset=False,sigma=(1.,1.,0.),df=True,f_values=[0.,1.],mean_value=0.5):
        if not hasattr(self,'presented_frames'):
            print('The file does not have presented frames')
            return 
        #
        
        return figures.show_reverse_correlation(self,channel=channel,unit=unit,lags=lags,lim_data=lim_data,normalize_spike_times_to_stimulus_onset=normalize_spike_times_to_stimulus_onset,sigma=sigma,df=df,f_values=f_values,mean_value=mean_value)
        """
        self.frame_times
        self.onsets
        self.frame_duration
        if normalize_spike_times_to_stimulus_onset:
            spike_times_tmp = self.flatten_dict(self.sort_spikes_to_trials(self.channels[channel].spike_times_clusters[unit],self.onsets))
        else:
            spike_times_tmp = self.channels[channel].spike_times_clusters[unit]
        if not df:
            rc,count = rf.reverse_correlation_stimulus_times(spike_times_tmp,self.presented_frames,self.frame_times,lags=lags)
            rf.show_receptive_field(rc,count,self.frame_duration,lags,lim_data=lim_data,sigma=sigma,p=self.p)
        else:
            rc,count = rf.reverse_correlation_stimulus_times_different_filters(spike_times_tmp,self.presented_frames,self.frame_times,lags=lags,values=f_values,mean_value=mean_value)
            # noise sd
            index = range(self.presented_frames.shape[2])
            index_p = numpy.random.permutation(index)
            rc_noise,count_noise = rf.reverse_correlation_stimulus_times_different_filters(spike_times_tmp,self.presented_frames[:,:,index_p],self.frame_times,lags=lags,values=f_values,mean_value=mean_value)
            rf.show_receptive_field_different_filter(rc,count,self.frame_duration,lags,lim_data=True,sigma=sigma,p=self.p)
        pylab.figtext(0.5,0.98,'File: %s Channel: %d unit: %d'%(self.name.split('/')[-1],channel,unit))
        return self.name.split('/')[-1][:-4]+'_Ch%d_unit%d'%(channel,unit),rc,rc_nosie


        """
    def show_receptive_fields_of_all_channels_and_units(self,sigma=(1.,1.,0.)):
        #pylab.ioff()
        for channelnumber,channel in self.channels.items():
            for unit in channel.valid_units:
                figures.show_reverse_correlation_max_response(self,channel=channelnumber,unit=unit,sigma=sigma,df=True,savefig=True)
                pylab.close('all')

        
        #pylab.ion()
  

    def show_reverse_correlation_max_response(self,channel=None,unit=None,lags=numpy.arange(-4.,-1,1,dtype='int'),lim_data=True,normalize_spike_times_to_stimulus_onset=False,sigma=(1.,1.,0.),df=True,f_values=[0.,1.],mean_value=0.5):
        if not hasattr(self,'presented_frames'):
            print('The file does not have presented frames')
            return 
        #
        
        return figures.show_reverse_correlation_max_response(self,channel=channel,unit=unit,lags=lags,lim_data=lim_data,normalize_spike_times_to_stimulus_onset=normalize_spike_times_to_stimulus_onset,sigma=sigma,df=df,f_values=f_values,mean_value=mean_value)





      
class SpikeSortThread(Thread):
    def __init__(self,channel,threshold,window,pca,threshold_detection,minimal_isi):
        Thread.__init__(self)
        self.channel = channel
        self.threshold = threshold
        self.window = window
        self.pca = pca
        self.threshold_detection = threshold_detection
        self.minimal_isi = minimal_isi
    def run(self):
        if self.threshold_detection or not self.channel._thresholds_detected:
            self.channel.threshold_detection(threshold=self.threshold,window=self.window,minimal_isi=self.minimal_isi)
        self.channel.KlustaKwik(pca=self.pca)

class Channel(son.Channel):
    
    def __init__(self,channel, filename,load_raw=True):
        son.Channel.__init__(self,channel, filename)
        if load_raw:
            self.value = self.data()
        self.spike_times = numpy.empty((0))
        self.spike_ticks = numpy.empty((0))
        self.tempdir = ''
        self._thresholds_detected = False
        self.spike_waves_clusters = {}
        self.spike_times_clusters = {}
        self.valid_units = []
    
    def assign_sorted_WaveMarks_to_channels(self,WaveMark_Channels=[]):
        
        for combination in WaveMark_Channels:
            wavemark_channel_number = combintation[0]
            adc_channel_number = combintation[1]
            self.get_channels([wavemark_channel_number,adc_channel_number])
            
            wm_channel = self.channels[wavemark_channel_number]
            adc_channel = self.channels[adc_channel_number]
            
            adc_channel.spike_times =  wm_channel.spike_times
            adc_channel.clusters =  wm_channel.clusters
            adc_channel.spike_waves =  wm_channel.spike_waves
            adc_channel.spike_waves_time =  wm_channel.spike_waves_time
           
    def set_valid_units(self,units=[]):
        self.valid_units = units

    def renumber_clusters(self,start_id=1):

        cluster_ids = numpy.unique(self.clusters)
        takes = {}
        for index, cluster_id in enumerate(cluster_ids):
            takes[index] = self.clusters == cluster_id

        for index, take in takes.items():
            self.clusters[take] = index
        self.sort_clustered_data()    
        
    def merge_clusters(self,clusters,new_cluster_id=None):
        if new_cluster_id == None:
            clusters = numpy.array(clusters)
            new_cluster_id = int(''.join(clusters.astype('str')))

        for cluster in clusters:
            self.clusters[self.clusters==cluster] = new_cluster_id

        self.sort_clustered_data()
        
        
    def select_valid_units(self,plot=True):
        if plot:
            self.plot_spike_clusters()
        
        input = raw_input('Please provide a list of unit ids which you considere to be valid neurons ')
        if not len(input) == 0:
            input = input.split(',')
            valid_units = numpy.array(input).astype(int)
        else:
            valid_units = []
        self.valid_units = valid_units
        

    def _filtfilt(self,b, a, input_vector):
        """input_vector has shape (n,1)"""
        forward = lfilter(b, a, input_vector, axis=0)
        return numpy.flipud(lfilter(b, a, numpy.flipud(forward), axis = 0))
    
    def filter(self,corner_freq,stop_freq):
        '''
        '''
        nyq_freq = self.sampling_rate/2.0
        wp = corner_freq/nyq_freq
	ws = stop_freq/nyq_freq
	gpass = 0.1
	gstop = 1.0
	order, Wn = buttord(wp, ws, gpass, gstop)
	if corner_freq > stop_freq:
            # lowpass
            b, a = butter(order, Wn, btype='high')
        elif corner_freq < stop_freq:
            # highpass
            b, a = butter(order, Wn, btype='low')
            
        return self._filtfilt(b,a,self.value)
    
    def get_value(self):
        self.value = self.data()
    value = property(get_value)
    
    def as_dict(self):
        data = {}
        #if self._thresholds_detected and self.type() == 'Adc':
        to_save = ['spike_waves','spike_times','spike_ticks','window','spike_waves_time','median','mean','clusters','threshold','effective_threshold','_thresholds_detected','spike_waves_pca','valid_units']
        for attribute in to_save:
            if hasattr(self,attribute):
                data[attribute] = eval('self.%s'%attribute)
        return data
    
    def set_from_dict(self,data):
        for key,value in data.items():
            exec('self.%s = value'%key)
            
    
    def get_time(self):
        return numpy.arange(0.,self.value.shape[0])/self.sampling_rate
    time = property(get_time)
    
    def get_spike_waves_times(self,data=None):
        if data == None:
            return self.spike_waves_time.repeat(self.spike_waves.shape[1]).reshape((self.spike_waves_time.shape[0],self.spike_waves.shape[1]))
        else:
            return self.spike_waves_time.repeat(data.shape[1]).reshape((self.spike_waves_time.shape[0],data.shape[1]))
    spike_waves_times = property(get_spike_waves_times)

    def get_dt(self):
        # in s
        return 1./self.sampling_rate 
    dt = property(get_dt)
    
    def get_sampling_rate(self):
        '''get_sampling_rate(filename,chan_num)\n
        This function extracts the sampling rate from the file and channel header information.\n
        expects: name of spike2 file filename, channel number chan_num
        returns: sampling rate sample_rate.'''
        chan_info = self.info
        head = self.fhead
        sample_rate = 0.0
        allowed_chan_kinds = [1,6,7,9]
        if chan_info.kind[0] in allowed_chan_kinds:
	    if head.systemID[0] < 6: # older spike2 versions
                sample_rate = 1./(chan_info.divide[0]*head.usPerTime[0]*head.timePerADC[0]*1e-6)
	    else: # newer spike2 versions
                sample_rate = 1./(chan_info.lChanDvd[0]*head.usPerTime[0]*head.dTimeBase[0])
        else:
	    print 'sampling rate not determinable - invalid channel type!' # insert exception handle later...
	    sample_rate = nan
        return sample_rate
    sampling_rate = property(get_sampling_rate)
    
    def _savetxt(self,filename, data, format, delimiter):
        """
        Due to the lack of savetxt in older versions of numpy
        we provide a cut-down version of that function.
        """
        f = open(filename, 'w')
        f.write('%d\n'%data.shape[1])
        for row in data:
            f.write(delimiter.join([format%val for val in row]) + '\n')
        f.close()
        
    def pca_spike_waves(self,components=2):
        Pcomponents,  Trans, fracVar = pylab.prepca(self.spike_waves)
        self.spike_waves_pca = Pcomponents[0:components,:]
        return self.spike_waves_pca
        
    def KlustaKwik(self,pca=True,components=2,max_spike_clusters=3):
        
        if not self.spike_times.shape[0] >1:
            return
        print('KlustaKwik of channel %d'%self.info.chan)
        if pca:
            data = self.pca_spike_waves(components=components)
        else:
            self.spike_waves_pca = None
            
        # write the data to the KlustaKwik
        self._klustakwik_filename = self.tempdir+'klustakwik_tmp'
        self._savetxt(self._klustakwik_filename+'.fet.%d'%self.info.chan,data.T,'%s',' ')
        # run
        os.system('KlustaKwik %s %d -MaxClusters %d -UseFeatures all -Screen 0'%(self._klustakwik_filename,self.info.chan,max_spike_clusters))
        # load clusters
        self.clusters = numpy.loadtxt(self._klustakwik_filename+'.clu.%d'%self.info.chan,skiprows=1)
        self.sort_clustered_data()
        
    def sort_clustered_data(self):
        self.spike_waves_clusters = {}
        self.spike_times_clusters = {}
        self.spike_waves_pca_clusters = {}
        if hasattr(self,'clusters'):
            for cluster in numpy.unique(self.clusters):
                take = self.clusters==cluster
                if numpy.sum(take) == 0:
                    # no spike in this cluster
                    pass
                else:
                    self.spike_waves_clusters[cluster] = self.spike_waves[:,take]
                    self.spike_times_clusters[cluster] = self.spike_times[take]
                    if hasattr(self,'spike_waves_pca'):
                        self.spike_waves_pca_clusters[cluster] = self.spike_waves_pca[:,take]
        else:
            print('Channel %d has no cluster'%self.info.chan)
        gc.collect()
    def threshold_detection(self,threshold=5.,window=[0.5,1.5],minimal_isi=0.):
        print('Threshold detection of channel %d'%self.info.chan)
        self._thresholds_detected = True
        # window in ms
        self.threshold = threshold
        self.window = window
        self.window_tick = numpy.round(numpy.array(window)/(self.dt*1000.))
        self.spike_waves_time = numpy.arange(-self.window_tick[0],self.window_tick[1],1)*self.dt*1000.
        self.median = numpy.median(numpy.abs(self.value))
        self.mean = numpy.mean(self.value)
        self.spike_waves = numpy.empty((int(sum(self.window_tick)),0))
        # thresholding
        self.effective_threshold = self.mean-(threshold*self.median)


        # we filter the channel with 1kHz to avoid double spike detection
        self.corner_freq = 1000.
        self.stop_freq = 3500.
        value_tmp = self.filter(self.corner_freq,self.stop_freq)
        # spike detection
        self._below_threshold = numpy.where(value_tmp < self.effective_threshold)[0]
        
        if len(self._below_threshold) <= 0:
            pass
        else:
            take = numpy.where(numpy.diff(self._below_threshold)>1)[0]+1
            take = numpy.append(0,take)
            spike_times_tmp = self.time[self._below_threshold][take]
            spike_ticks_tmp = self._below_threshold[take]
            
            spike_waves = numpy.empty((self.spike_waves_time.shape[0],spike_times_tmp.shape[0]))
            spike_times = numpy.empty(spike_times_tmp.shape[0])
            spike_ticks = numpy.empty(spike_times_tmp.shape[0])
            count = 0
            shifts = numpy.empty(spike_times_tmp.shape[0])
            
            for spike_tick in spike_ticks_tmp:
                if spike_tick-self.window_tick[0]<0:
                    continue
                tmp = self.value[spike_tick-self.window_tick[0]:spike_tick+self.window_tick[1]]
                shift_min = numpy.where(tmp==tmp.min())[0][0]
                shift_max = numpy.where(tmp==tmp[self.window_tick[0]::].max())[0][0]
                
                if shift_min < shift_max:
                    shift = shift_min
                else:
                    shift=0
                    
                cutout_shift = shift-(self.window_tick[0])
                tmp = self.value[cutout_shift+spike_tick-self.window_tick[0]:cutout_shift+spike_tick+self.window_tick[1]]
                if tmp.shape[0] == spike_waves.shape[0]:
                    spike_waves[:,count] = tmp
                    spike_times[count] = spike_times_tmp[count]
                    spike_ticks[count] = spike_tick
                    shifts[count] = shift
                    count += 1
            
            self.spike_waves = spike_waves[:,0:count]
            self.spike_times = spike_times[0:count]
            self.spike_ticks = spike_ticks[0:count]
            self.shifts = shifts[0:count]
            #if minimal_isi >0.:
                # remove spike times which follow a spike below minimal_isi
            #    lager_than_minimal_isi = numpy.diff(self.spike_times)>=minimal_isi
                # we have to add one True because diff reduces the data set
            #    lager_than_minimal_isi = numpy.append(True,lager_than_minimal_isi)
            #    self.spike_waves = self.spike_waves[:,lager_than_minimal_isi]
            #    self.spike_times = self.spike_times[lager_than_minimal_isi]
             #   self.spike_ticks = self.spike_ticks[lager_than_minimal_isi]
            
            print('Found %d spikes'%self.spike_times.shape[0])

    def plot(self,xlim=None):
        '''
        xlim = [1.,4.] for plotting from 1 to 4th second 
        '''
        pylab.figure()
        if xlim is not None:
            take = (self.time>=xlim[0])&(self.time<xlim[1])
            pylab.plot(self.time[take],self.value[take])
        else:
            pylab.plot(self.time,self.value)
        
        if self._thresholds_detected:
            for unit, spike_times in self.spike_times_clusters.items():
                pylab.plot(spike_times,numpy.ones(len(spike_times))*self.effective_threshold,'o',label='unit %d'%unit)
            if xlim is not None:
                pylab.xlim(xlim)
            pylab.legend(numpoints=1)
        pylab.title('Channel %d'%self.info.chan)
                
    def plot_spike_clusters(self):
        unit_keys = self.spike_waves_clusters.keys()
        # isi
        pylab.figure()
        x = numpy.ceil(numpy.sqrt(len(unit_keys)))
        y = numpy.ceil(numpy.sqrt(len(unit_keys)))
        for index, unit in enumerate(unit_keys):
            pylab.subplot(x,y,index+1)
            spikes = self.spike_times_clusters[unit]*1000.# to get it into ms
            isi = numpy.diff(spikes)
            isi_dist = numpy.histogram(isi,bins=numpy.arange(0,100.,1.))
            pylab.plot(isi_dist[1][0:-1],isi_dist[0])
            pylab.title('Channel: %d unit: %d #spikes %d'%(self.info.chan,unit,spikes.shape[0]))
        # spike waves
        pylab.figure()
        ylim = (self.spike_waves.min(),self.spike_waves.max())
        
        for index, unit in enumerate(unit_keys):
            pylab.subplot(x,y,index+1)
            data = self.spike_waves_clusters[unit]
            pylab.plot(self.get_spike_waves_times(data=data),data)
            pylab.title('Channel: %d unit: %d #spikes %d'%(self.info.chan,unit,data.shape[1]))
            pylab.ylim(ylim)
            pylab.figtext(0.5,0.01,self.fhead.name,ha='center')
        if len(self.spike_waves_pca_clusters.keys())>0:
            # cluster
            pylab.figure()
            for index, unit in enumerate(unit_keys):
                data = self.spike_waves_pca_clusters[unit]
                if data.shape[0] > 0:
                    pylab.plot(data[0,:],data[1,:],'.',label='unit %d, #spikes %d '%(unit,data.shape[1]))
            pylab.legend(numpoints=1)
            pylab.title('Channel: %d'%(self.info.chan))

class AdcMarkChannel(Channel):
    def __init__(self,channel, filename,load_raw=True):
        Channel.__init__(self,channel, filename,load_raw=load_raw)
        self.spike_times_from_markers()
        self.sort_clustered_data()
    
    def spike_times_from_markers(self):
        data = self.value
        
        number_spike_times = len(data)
        window = data[0][2].shape[0]
        
        self.spike_times = numpy.empty(number_spike_times)
        self.clusters = numpy.empty(number_spike_times)
        self.spike_waves = numpy.empty((window,number_spike_times))
        self.spike_waves_time = numpy.arange(0,window,1)*self.dt*1000.
        
        for index, spike in enumerate(data):
            self.spike_times[index] = spike[0]
            self.clusters[index] = spike[1][0]
            self.spike_waves[:,index] = spike[2]
            
   
class EventRiseChannel(Channel):
    pass
        

class StimulusChannelPT(Channel):
    
    def get_offsets(self):
        return self.onsets+numpy.diff(self.onsets).mean()
    offsets = property(get_offsets)
    
    def get_onsets(self):
        # we compute where in the data we have large changes
        thr = 200
        # since diff works in the interval we loose one data point
        changes = numpy.diff(self.value)
        # check where the diff is lager than the threshold
        changes_more_thr = numpy.where(changes>thr)[0]
        # 
        interval_changes_more_thr = numpy.diff(changes_more_thr)
        # where are these intervals lager than 230
        thr = 100
        onsets_intervals = numpy.where(interval_changes_more_thr>thr)[0]
        # +1 because diff gives value for n-n+1, but we need n+1
        onsets = changes_more_thr[0]+1
        onsets = numpy.append(onsets,changes_more_thr[onsets_intervals+1])
        
        # onsets are calculated, we override the property(get_onsets)
        self.onsets = numpy.array(onsets)/self.sampling_rate
        return numpy.array(onsets)/self.sampling_rate
    onsets = property(get_onsets)
