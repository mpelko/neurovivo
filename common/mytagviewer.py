from Cheetah.Template import Template
from morphforge.core import SettingsMgr
from morphforge.core import Flatten
from morphforge.core import unit


class MyTagViewer(object):
    def __init__(self, results, tags=(["Voltage"],  ["Current"],
    ["Conductance"],["State"], ["State (Tau)"], ["State (SteddyState)"],
    ["Tau"]), figsize=(12, 8), timeranges=(None,), filename=None,
    annotation=None, jitter=0.0, tracecolour=None, legendlabel="${result.name} - ${trace.comment}", 
    haslegend=True,sorttraces = True, colors=None):
        #from morphforge.core.pylab import mpl
        #from morphforge.core.pylab import mpl
        
        self.sorttraces = sorttraces
        self.haslegend = haslegend
        # Filter out tag sets which actually have relevant traces:
        allTraceObjs = Flatten([ result.traces for result in results ]) 
        
        self.tags = []
        for tgs in tags:
            validTraces = allTraceObjs[:] 
            for tag in tgs:
                validTraces = [ trc for trc in validTraces if tag in trc.tags ]
            if validTraces != []: self.tags.append(tgs)
        
        self.legendlabel = legendlabel
        self.results = results
        
        self.figsize = figsize
        
        self.ignoreInitial = unit("5:ms")
        self.axisPadding = 0.3
        self.axisPadding = 0.1
        
        self.jitter = float(jitter)
        self.timeranges = timeranges #if timeranges else [None] 
        
        self.tracecolour = tracecolour
        self.colors = colors
        
        self.annotation = annotation
        self.filename = filename
        
        self.fig = None
        self.subfigs = []
        self.Render()
    
    def getPlotMinMax(self, trace, timeTrace, plotCurrentMin, plotCurrentMax):    
        timeUnits = timeTrace.units
        dataUnits = trace.units

        trMinRaw = min(trace.rescale(dataUnits).magnitude[ timeTrace.magnitude > self.ignoreInitial.rescale(timeUnits).magnitude ])
        trMaxRaw = max(trace.rescale(dataUnits).magnitude[ timeTrace.magnitude > self.ignoreInitial.rescale(timeUnits).magnitude ])
        trMin = trMinRaw * dataUnits
        trMax = trMaxRaw * dataUnits

        plotCurrentMin = trMin if plotCurrentMin == None or trMin < plotCurrentMin else plotCurrentMin  
        plotCurrentMax = trMax if plotCurrentMax == None or trMax > plotCurrentMax else plotCurrentMax
        return plotCurrentMin, plotCurrentMax
        
    def Render(self):
        from morphforge.core import pylab_wrapper
        self.fig = pylab_wrapper.mpl.figure(figsize=self.figsize)
        
        print "Time Ranges:", self.timeranges
        
        nGraphs = len(self.tags)
        nTimeRanges = len(self.timeranges)
        print self.colors 
        
        print self.tags
        
        for i, tags in enumerate(self.tags):
            for iT, timeRange in enumerate(self.timeranges):
                subfigIndex = (i) * nTimeRanges + iT + 1
                print nGraphs, nTimeRanges, subfigIndex, i
                s = self.fig.add_subplot(nGraphs, nTimeRanges, subfigIndex)
                self.subfigs.append((s))
                i = 0
                # Find Traces with the relevant tags - start with complete list then slowly remove:
                spMax = None
                spMin = None
                for result in self.results:
                    validTraces = result.traces
                    for tag in tags:
                        validTraces = [ trc for trc in validTraces if tag in trc.tags ]
                    
                    traceList =sorted(validTraces, key=lambda t:t.name) if self.sorttraces else validTraces 
                    print tags, validTraces
                    currentUnit = traceList[0]._data.units
                    for trc in traceList:
                        spMin, spMax = self.getPlotMinMax(trc._data, trc._time, spMin, spMax)
                        jitterFactor = float(i) * self.jitter
                        i = i + 1.0
                        print jitterFactor
                        
                        label = Template(self.legendlabel, [{"trace":trc, "result":result}]).respond()

                        if self.colors != None:
                            self.tracecolour = self.colors[ int(i)-1 ]
                        
                        if self.tracecolour != None:
                            pylab_wrapper.mpl.plot(trc._time, trc._data.rescale(currentUnit).magnitude + jitterFactor, color=self.tracecolour, label=label)
                        else:
                            pylab_wrapper.mpl.plot(trc._time, trc._data.rescale(currentUnit).magnitude,  label=label)
                        
                
                pylab_wrapper.mpl.xlabel("Time (ms)")
                pylab_wrapper.mpl.ylabel("-".join(tags) + " (%s)"%currentUnit.dimensionality.string )
                pylab_wrapper.mpl.grid("on")
                
                rg = spMax - spMin 
                pylab_wrapper.mpl.ylim((spMin - rg * self.axisPadding, spMax + rg * self.axisPadding))

                if timeRange:
                    pylab_wrapper.mpl.xlim(timeRange)

                if self.haslegend:
                  pylab_wrapper.mpl.legend()
       
        if self.filename:
            self.fig.savefig(self.filename)
            if SettingsMgr.showAllPlots(): pylab_wrapper.mpl.show()
        else:
            pass
            #mpl.show()

__all__ = ["MyTagViewer"]


