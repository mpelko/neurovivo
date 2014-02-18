from __init__ import __doc__
import numpy

def _get_block_headers(fhead, chanInfo):
    '''
    Returns a matrix containing the Son data block headers for a channel.
    'fhead' is a FileHeader instance, and 'chanInfo' is a ChannelInfo instance.

    The returned header in memory contains, for each disk block,
    a column with rows 0-4 representing:
         Offset to start of block in file
         Start time in clock ticks
         End time in clock ticks
         Chan number
         Items
    '''
    from scipy import io, zeros, arange, take
    if chanInfo.firstblock == -1:
        raise ValueError, 'No data on channel %i' % chanInfo.chan
    succBlock = 1
    # Pre-allocate memory for header data:
    header = zeros([6, chanInfo.blocks], int)
    # Get first data block:
    fhead.fid.seek(chanInfo.firstblock)
    ####################################MAX

    # Last and next block pointers, start and end times in clock ticks
    header[0:4, 0] = numpy.fromfile(fhead.fid, numpy.int32, count=4)
    # Channel number and number of items in block:
    header[4:6, 0] = numpy.fromfile(fhead.fid, numpy.int16, count=2)

    ####################################MAX
    """
    # Last and next block pointers, start and end times in clock ticks
    header[0:4, 0] = io.fread(fhead.fid, 4, 'l')
    print("header[0:4, 0]")
    print(header[0:4, 0])
    # Channel number and number of items in block:
    header[4:6, 0] = io.fread(fhead.fid, 2, 'h')
    """
    ####################################MAX
    ####################################MAX
    #print "header[succBlock, 0]", header[succBlock, 0]
    ####################################MAX
    # If only one block:
    if header[succBlock, 0] == -1: header[0, 0] = int(chanInfo.firstblock)
    # Loop if more blocks:
    else:
        fhead.fid.seek(header[succBlock, 0])
        for i in arange(1, chanInfo.blocks):
            ####################################MAX

            header[0:4, i] = numpy.fromfile(fhead.fid, numpy.int32, count=4)
            header[4:6, i] = numpy.fromfile(fhead.fid, numpy.int16, count=2)

            ####################################MAX
            """
            header[0:4, i] = io.fread(fhead.fid, 4, 'l')
            header[4:6, i] = io.fread(fhead.fid, 2, 'h')
            """
            ####################################MAX
            if header[succBlock, i] > 0:
                fhead.fid.seek(header[succBlock, i])
            header[0, i-1] = header[0, i]
        # Replace predBlock for previous column:
        header[0, -1] = header[1, -2]
    # Delete succBlock data:
    header = take(header, (0,2,3,4,5),axis=0)
    return header

class FileHeader:
    '''
    Reads the file header for a Son file.
    >>> fhead = FileHeader(filename)
    '''
    def __init__(self, filename):
        from scipy import io, zeros
        self.name = filename
        self.fid = open(filename, 'rb')
        self.fid.seek(0)
        ####################################MAX

        self.systemID   = numpy.fromfile(self.fid, numpy.int16, count=1)
        self.copyright  = self.fid.read(10)
        self.creator    = self.fid.read(8)
        self.usPerTime  = numpy.fromfile(self.fid, numpy.int16, count=1)
        self.timePerADC = numpy.fromfile(self.fid, numpy.int16, count=1)
        self.filestate  = numpy.fromfile(self.fid, numpy.int16, count=1)
        self.firstdata  = numpy.fromfile(self.fid, numpy.int32, count=1)
        self.channels   = numpy.fromfile(self.fid, numpy.int16, count=1)
        ####################################MAX
        # print "self.channels", self.channels
        ####################################MAX
        self.chansize   = numpy.fromfile(self.fid, numpy.int16, count=1)
        self.extraData  = numpy.fromfile(self.fid, numpy.int16, count=1)
        self.buffersize = numpy.fromfile(self.fid, numpy.int16, count=1)
        self.osFormat   = numpy.fromfile(self.fid, numpy.int16, count=1)
        self.maxFTime   = numpy.fromfile(self.fid, numpy.int32, count=1)
        self.dTimeBase  = numpy.fromfile(self.fid, numpy.float64, count=1)
        if self.systemID < 6: self.dTimeBase = 1e-6
        self.timeDate = {
                'Detail' : numpy.fromfile(self.fid, numpy.int8, count=6),
                'Year' :   numpy.fromfile(self.fid, numpy.int16, count=1)}
        ####################################MAX
        """
        self.systemID   = io.fread(self.fid, 1, 'h')
        self.copyright  = self.fid.read(10)
        self.creator    = self.fid.read(8)
        self.usPerTime  = io.fread(self.fid, 1, 'h')
        self.timePerADC = io.fread(self.fid, 1, 'h')
        self.filestate  = io.fread(self.fid, 1, 'h')
        self.firstdata  = io.fread(self.fid, 1, 'l')
        self.channels   = io.fread(self.fid, 1, 'h')
        self.chansize   = io.fread(self.fid, 1, 'h')
        self.extraData  = io.fread(self.fid, 1, 'h')
        self.buffersize = io.fread(self.fid, 1, 'h')
        self.osFormat   = io.fread(self.fid, 1, 'h')
        self.maxFTime   = io.fread(self.fid, 1, 'l')
        self.dTimeBase  = io.fread(self.fid, 1, 'd')
        print("self.dTimeBase")
        print(self.dTimeBase)

        if self.systemID < 6: self.dTimeBase = 1e-6
        self.timeDate = {
		        'Detail' : io.fread(self.fid, 6, 'B'),
		        'Year' :   io.fread(self.fid, 1, 'h')}
        print("self.timeDate")
        print(self.timeDate)
        """

        ####################################MAX
        if self.systemID < 6:
            self.timeDate['Detail'] = zeros(6)
            self.timeDate['Year'] = 0
        pad = self.fid.read(52)
        self.fileComment = {}
        pointer = self.fid.tell()
        for i in range(1, 6):
            ####################################MAX

            bytes = numpy.fromfile(self.fid, numpy.int8, count=1)

            ####################################MAX
            """		
            bytes = io.fread(self.fid, 1, 'B')
            """
            ####################################MAX
            self.fileComment[i] = self.fid.read(bytes)
            pointer = pointer + 80
            self.fid.seek(pointer)
	
    def __del__(self):
        self.fid.close()

    def getSamplingRateInfo(self, IDChannel):
        info = _ChannelInfo(self, IDChannel)
        return info.idealRate

    def chan_list(self):
        '''
        Prints a table with details of active channels in a Son file.
        '''
        print '-' * 38
        print '%-5s %-10s %-10s %s' %('Chan', 'Type', 'Title', 'Comment')
        print '-' * 38
        for chan in range(1, self.channels+1):
            info = _ChannelInfo(self, chan)
            if info.kind > 0:
                print '%-5i %-10s %-10s %s'	%(
                        chan, info.type(), info.title, info.comment)
        print '-' * 38
	
    def return_chan_list(self):
        '''
        Returns a list with details of active channels in a Son file.
        '''
        channels = []
        for chan in range(1, self.channels+1):
            info = _ChannelInfo(self, chan)
            if info.kind > 0:
                channels.append(chan)
			
        return channels
		
    def _ticks_to_seconds(self, timestamp, timeunits):
        '''
        Scales a timestamp vector.

        'timeunits' is one of
                'ticks': returns the time in base clock ticks.
                'microseconds', 'milliseconds' or 'seconds': scales the
                output to the appropriate unit.
        '''
        if timeunits is 'ticks':
            pass
        elif timeunits is 'microseconds':
            timestamp = timestamp * (self.usPerTime * self.dTimeBase * 1e6)
        elif timeunits is 'milliseconds':
            timestamp = timestamp * (self.usPerTime * self.dTimeBase * 1e3)
        elif timeunits is 'seconds': # default, time in seconds
            timestamp = timestamp * (self.usPerTime * self.dTimeBase)
        else: raise ValueError
        return timestamp #, timeunits

    def __repr__(self):
        return 'File header from Son file %s' % self.name

class _ChannelInfo:
    '''
    Reads the Son file channel header for a channel.
    info = ChannelInfo(fhead, channel)
    where
	    'fhead' is an instance containing the file header, and
	    'channel' is the channel number.

    The resulting instance will contain values which follow the CED
    disk header structure.
    '''
    def __init__(self, fhead, channel):
        from scipy import io
        self.chan = channel
        # Offset due to file header and preceding channel headers:
        base = 512 + (140*(channel-1))
        fhead.fid.seek(base)

        ####################################MAX      
        self.delSize      = numpy.fromfile(fhead.fid, numpy.int16, count=1)
        self.nextDelBlock = numpy.fromfile(fhead.fid, numpy.int32, count=1)
        self.firstblock   = numpy.fromfile(fhead.fid, numpy.int32, count=1)
        ####################################MAX
        # print "self.firstblock", self.firstblock
        ####################################MAX
        self.lastblock    = numpy.fromfile(fhead.fid, numpy.int32, count=1)
        self.blocks       = numpy.fromfile(fhead.fid, numpy.int16, count=1)
        ####################################MAX
        # print "self.blocks", self.blocks
        ####################################MAX
        self.nExtra       = numpy.fromfile(fhead.fid, numpy.int16, count=1)
        self.preTrig      = numpy.fromfile(fhead.fid, numpy.int16, count=1)
        self.free0        = numpy.fromfile(fhead.fid, numpy.int16, count=1)
        self.phySz        = numpy.fromfile(fhead.fid, numpy.int16, count=1)
        self.maxData      = numpy.fromfile(fhead.fid, numpy.int16, count=1)
        bytes             = numpy.fromfile(fhead.fid, numpy.int8, count=1)
        pointer           = fhead.fid.tell()
        self.comment      = fhead.fid.read(bytes)
        fhead.fid.seek(pointer+71)
        self.maxChanTime  = numpy.fromfile(fhead.fid, numpy.int32, count=1)
        self.lChanDvd     = numpy.fromfile(fhead.fid, numpy.int32, count=1)
        self.phyChan      = numpy.fromfile(fhead.fid, numpy.int16, count=1)
        bytes             = numpy.fromfile(fhead.fid, numpy.int8, count=1)
        pointer           = fhead.fid.tell()
        self.title        = fhead.fid.read(bytes)
        fhead.fid.seek(pointer+9)
        self.idealRate    = numpy.fromfile(fhead.fid, numpy.float32, count=1)
        ####################################MAX
        self.kind         = numpy.fromfile(fhead.fid, numpy.uint8, count=1)
        ####################################MAX
        #self.kind         = numpy.fromfile(fhead.fid, numpy.int8, count=1)
        ####################################MAX
        pad               = numpy.fromfile(fhead.fid, numpy.bool_, count=1)

        if self.kind in [1, 6]:
            self.scale    = numpy.fromfile(fhead.fid, numpy.float32, count=1)
            self.offset   = numpy.fromfile(fhead.fid, numpy.float32, count=1)
            bytes         = numpy.fromfile(fhead.fid, numpy.int8, count=1)
            pointer       = fhead.fid.tell()
            self.units    = fhead.fid.read(bytes)
            fhead.fid.seek(pointer+5)
            if fhead.systemID < 6: 
                ####################################MAX
                self.divide = numpy.fromfile(fhead.fid, numpy.int16, count=1)
                self.dt = self.divide*fhead.usPerTime*fhead.timePerADC*1e-6
                ####################################MAX
                #self.divide = numpy.fromfile(fhead.fid, numpy.int32, count=1)
                ####################################MAX
            else: 
                ####################################MAX
                self.interleave = numpy.fromfile(fhead.fid, numpy.int16, count=1)
                self.dt = self.lChanDvd*fhead.usPerTime*fhead.dTimeBase
                ####################################MAX
                #self.interleave = numpy.fromfile(fhead.fid, numpy.int32, count=1)
                ####################################MAX
        elif self.kind in [7, 9]:
            self.min       = numpy.fromfile(fhead.fid, numpy.float32, count=1)
            self.max       = numpy.fromfile(fhead.fid, numpy.float32, count=1)
            bytes          = numpy.fromfile(fhead.fid, numpy.int8, count=1)
            pointer        = fhead.fid.tell()
            self.units     = fhead.fid.read(bytes)
            fhead.fid.seek(pointer+5)
            if fhead.systemID < 6: 
                ####################################MAX
                self.divide = numpy.fromfile(fhead.fid, numpy.int16, count=1)
                self.dt = self.divide*fhead.usPerTime*fhead.timePerADC*1e-6
                ####################################MAX
                #self.divide = numpy.fromfile(fhead.fid, numpy.int32, count=1)
                ####################################MAX
            else: 
                ####################################MAX
                self.interleave = numpy.fromfile(fhead.fid, numpy.int16, count=1)
                self.dt = self.lChanDvd*fhead.usPerTime*fhead.dTimeBase
                ####################################MAX
                #self.interleave = numpy.fromfile(fhead.fid, numpy.int32, count=1)
                ####################################MAX

        elif self.kind in [4]:
            self.initLow   = numpy.fromfile(fhead.fid, numpy.int8, count=1)
            self.nextLow   = numpy.fromfile(fhead.fid, numpy.int8, count=1)
		
        ##################################
        """
        self.delSize      = io.fread(fhead.fid, 1, 'h')
        self.nextDelBlock = io.fread(fhead.fid, 1, 'l')
        self.firstblock   = io.fread(fhead.fid, 1, 'l')
        self.lastblock    = io.fread(fhead.fid, 1, 'l')
        self.blocks       = io.fread(fhead.fid, 1, 'h')
        self.nExtra       = io.fread(fhead.fid, 1, 'h')
        self.preTrig      = io.fread(fhead.fid, 1, 'h')
        self.free0        = io.fread(fhead.fid, 1, 'h')
        self.phySz        = io.fread(fhead.fid, 1, 'h')
        self.maxData      = io.fread(fhead.fid, 1, 'h')
        bytes             = io.fread(fhead.fid, 1, 'B')
        pointer           = fhead.fid.tell()
        self.comment      = fhead.fid.read(bytes)
        fhead.fid.seek(pointer+71)
        self.maxChanTime  = io.fread(fhead.fid, 1, 'l')
        self.lChanDvd     = io.fread(fhead.fid, 1, 'l')
        self.phyChan      = io.fread(fhead.fid, 1, 'h')
        bytes             = io.fread(fhead.fid, 1, 'B')
        pointer           = fhead.fid.tell()
        self.title        = fhead.fid.read(bytes)
        fhead.fid.seek(pointer+9)
        self.idealRate    = io.fread(fhead.fid, 1, 'f')
        self.kind         = io.fread(fhead.fid, 1, 'B')
        pad               = io.fread(fhead.fid, 1, 'b')

        if self.kind in [1, 6]:
            self.scale    = io.fread(fhead.fid, 1, 'f')
            self.offset   = io.fread(fhead.fid, 1, 'f')
            bytes         = io.fread(fhead.fid, 1, 'B')
            pointer       = fhead.fid.tell()
            self.units    = fhead.fid.read(bytes)
            fhead.fid.seek(pointer+5)
            if fhead.systemID < 6: self.divide = io.fread(fhead.fid, 1, 'l')
            else: self.interleave = io.fread(fhead.fid, 1, 'l')
        elif self.kind in [7, 9]:
            self.min       = io.fread(fhead.fid, 1, 'f')
            self.max       = io.fread(fhead.fid, 1, 'f')
            bytes          = io.fread(fhead.fid, 1, 'B')
            pointer        = fhead.fid.tell()
            self.units     = fhead.fid.read(bytes)
            fhead.fid.seek(pointer+5)
            if fhead.systemID < 6: self.divide = io.fread(fhead.fid, 1, 'l')
            else: self.interleave = io.fread(fhead.fid, 1, 'l')
        elif self.kind in [4]:
            self.initLow   = io.fread(fhead.fid, 1, 'B')
            self.nextLow   = io.fread(fhead.fid, 1, 'B')
        """	
        ##################################
        try:
            self.units        = self.units.strip()
        except AttributeError:
            pass

    def type(self):
        '''
        Adc        16-bit integer waveform data
        EventFall  event data, times taken on low edge of pulse
        EventRise  event data, times taken on high edge of pulse
        EventBoth  event data, times taken on both edges of pulse
        Marker     an event time plus four identifying bytes
        AdcMark    16-bit integer waveform transient shapes
        RealMark   array of real numbers attached to a marker
        TextMark   string of text attached to a marker
        RealWave   32-bit floating point waveforms
        '''
        if   self.kind == 1: type = 'Adc'
        elif self.kind == 2: type = 'EventFall'
        elif self.kind == 3: type = 'EventRise'
        elif self.kind == 4: type = 'EventBoth'
        elif self.kind == 5: type = 'Marker'
        elif self.kind == 6: type = 'AdcMark'
        elif self.kind == 7: type = 'RealMark'
        elif self.kind == 8: type = 'TextMark'
        elif self.kind == 9: type = 'RealWave'
        return type

    def __repr__(self):
        return 'Information from channel %i' % self.chan
	
	

class Channel:
    '''
    Creates an instance whose attributes are those of a channel from a
    Son file.

    >>> chan = Channel(chanNumber, filename)
    chan.info will contain the channel information, and
    chan.fhead will contain the SON file header

    Useful calls:
    >>> chan.type()
    >>> chan.data()
    '''
    def __init__(self, channel, filename):
        self.fhead = FileHeader(filename)
        self.info = _ChannelInfo(self.fhead, channel)
        self.blockheader = _get_block_headers(self.fhead, self.info)

    def __repr__(self):
        return  'Channel %i from %s' % (self.info.chan, self.fhead.name)

    def data(self, startTime=-1, stopTime=-1, timeunits='seconds', as_float=True):
        '''
        Returns channel data from a Son file.

        If only a fragment of data is required, limits (in blocks) can be
        defined in 'start' and 'stop'.

        If only 'start' is set, data from that block will be returned.

        'timeunits' can be 'ticks', 'microseconds', 'milliseconds', or
        'seconds', and cause times to be scaled to the appropriate unit.

        'as_float' only makes sense for waveform channels (ie channels of
        type AdcMark, RealMark, Adc, or RealWave). If True, returns data as
        floating point values (scaling and applying offset for Adc/AdcMark
        channel data). Else, data will be of type int16.
        '''
        startBlock, stopBlock =self._time_to_blocks(startTime, stopTime)
        #		print("startTime = ", startTime, "stopTime = ", stopTime)
        #		print("startBlock = ", startBlock, "stopBlock = ", stopBlock)
        if self.info.kind in [1, 9]:
            from _waveform import data
            return data(self, startBlock, stopBlock, timeunits, as_float)
        elif self.info.kind in [2, 3, 4, 5, 6, 7, 8]:
            from _marker import data
            return data(self, startBlock, stopBlock, timeunits)

    def type(self):
        '''
        Adc        16-bit integer waveform data
        EventFall  event data, times taken on low edge of pulse
        EventRise  event data, times taken on high edge of pulse
        EventBoth  event data, times taken on both edges of pulse
        Marker     an event time plus four identifying bytes
        AdcMark    16-bit integer waveform transient shapes
        RealMark   array of real numbers attached to a marker
        TextMark   string of text attached to a marker
        RealWave   32-bit floating point waveforms
        '''
        return self.info.type()

    def __del__(self):
        self.fhead.fid.close()
	
    def _time_to_blocks(self,startTime, stopTime):
        ####################################MAX
        startBlock = -1
        stopBlock = -1
        ####################################MAX
        if startTime != -1:
            startTics = int((float(startTime)/(self.fhead.usPerTime*self.fhead.timePerADC))/self.fhead.dTimeBase)
            #print startTics
            startBlock = numpy.searchsorted(self.blockheader[1,:],startTics)
            #print start
            if startBlock == self.blockheader.shape[1]:
                print 'start time exceeds file length... aborting'
                return 0
        if startBlock == 0:
            startBlock =1
        if stopTime != -1:
            stopTics = int((float(stopTime)/(self.fhead.usPerTime*self.fhead.timePerADC))/self.fhead.dTimeBase)
            #print stopTics
            stopBlock = numpy.searchsorted(self.blockheader[2,:],stopTics)+1
            #print stop
            if stopBlock == self.blockheader.shape[1]:
                print 'stop time exceeds file length... aborting'
                return 0
        return startBlock,stopBlock

class File:
    '''
    Creates an instance whose attributes are those of a Son file.

    Useful calls:
    >>> f = File(filename)       # then:
    >>> f.chanlist()             # will print a list of available channels, and
    >>> f.getchannel(chanNumber) # will create a channel instance named after the
	                             # channel number provided (see class Channel).
    '''
    def __init__(self, filename):
        self.name = filename
        self.header = FileHeader(self.name)

    def __repr__(self):
        return  'Son file %s' % self.name

    def get_channel(self, channel):
        '''
        Set an attribute of current Son file to the channel provided.
        '''
        attr = 'ch%02i' % channel
        #		print(attr)
        self.__dict__[attr] = Channel(channel, self.name)

    def chan_list(self):
        '''
        Prints a table with details of active channels in a Son file.
        '''
        self.header.chan_list()

    def getSamplingRateInfo(self, IDChannel):
        return self.header.usPerTime*1000000

    def getNbChannelsInfo(self):
        nbChannels =0
        for chan in range(1, self.header.channels+1):
            info = _ChannelInfo(self.header, chan)
            if info.kind > 0:
                nbChannels = nbChannels+1
        return nbChannels




