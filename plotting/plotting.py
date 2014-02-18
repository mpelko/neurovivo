import matplotlib.pylab as pyl
import numpy as np

colors = ['tomato', 'RoyalBlue', 'SeaGreen', 'Sienna', "MediumVioletRed"]
# -----------------------------------------------------------------------------
def plot_with_std(x, y, std, color="b", alpha=0.4, label=None, axes=None, rectify=False, *args, **kwargs):
# -----------------------------------------------------------------------------
    """
    Plots the 2D function with the std area around.
    """
    import matplotlib.pylab as pyl

    std = np.array(std); x = np.array(x); y = np.array(y)
    x_f = np.concatenate((x,x[::-1]))
    y_low = y - std
    y_high = y + std
    if not rectify == False:
        y[y<=rectify]=rectify
        y_low[y_low<=rectify]=rectify
        y_high[y_high<=rectify]=rectify
    y_f = np.concatenate((y_low, y_high[::-1]))
    if not axes == None:
        pyl=axes
    pyl.fill(x_f, y_f, color=color, alpha=alpha, *args, **kwargs)
    pyl.plot(x, y, color=color, linewidth=2, label=label, *args, **kwargs)
    #pyl.plot(x, y_high, color=color, linewidth=1, *args, **kwargs)
    #pyl.plot(x, y_low, color=color, linewidth=1, *args, **kwargs)


# -----------------------------------------------------------------------------
def plot_spike_trains(trains, binary=True, full_x=True, time=None):
# -----------------------------------------------------------------------------
    """
    WARNING: This is an obsolete function, only kept for legacy code. Use 
    plot_spike_population instead.
    
    Visulizes the spike trains. trains are given as a list of spike trains in 
    binary code if binary is True or in sparce code otherwise.
    """
    if binary:
        from correlations.correlated_input import binary_to_sparse_train
    
    if binary:
        time = len(trains[0])
    else:
        assert time, "Specify the total time of spike trains!"
        
    y_width = 1
    for i, train in enumerate(trains):
        if binary:
            train = np.array(binary_to_sparse_train(train)[0])
        y_position = np.zeros(len(train)) + i * y_width
        pyl.plot(train, y_position,".", color = "blue")
    #pyl.xlabel("t [ms]")
    #pyl.ylabel("N (neuron spike trains)")
    pyl.ylim((-y_width,(i + 1) * y_width))
    if full_x:
        pyl.xlim(0, time)


# -----------------------------------------------------------------------------
def plot_train_correlation(train1, train2, spike_width = 1):
# -----------------------------------------------------------------------------
    from correlations.correlation_analysis import train_correlation

    correlation = train_correlation(train1, train2, spike_width)
    d_length = (len(correlation)/2)
    pyl.plot(range(-d_length, d_length+1), correlation)
    pyl.xlabel("t")
    pyl.ylabel("corr")


# -----------------------------------------------------------------------------
def plot_trace(trace, *args, **kwargs):
# -----------------------------------------------------------------------------
    #from common.trace import Trace
    assert type(trace).__name__ == 'Trace', "Input must be of the type Trace."

    if trace.label:
        pyl.plot(trace._time, trace._data, *args, label=trace.label, **kwargs)
    else:
        pyl.plot(trace._time, trace._data, *args, **kwargs)
    pyl.xlabel("time [{0}]".format(str(trace._time.units).split()[-1]))
    if trace.name:
        pyl.ylabel("{0} [{1}]".format(trace.name,\
                                      str(trace._data.units).split()[-1]))

# -----------------------------------------------------------------------------
def plot_simpletrace(simpletrace, *args, **kwargs):
# -----------------------------------------------------------------------------
    #from common.trace import Trace
    assert type(simpletrace).__name__ == 'SimpleTrace', "Input must be of the type SimpleTrace."

    time = np.arange(len(simpletrace._data))*simpletrace._dt
    pyl.plot(time, simpletrace._data, *args, **kwargs)
    pyl.xlabel("time [{0}]".format("ms"))
    pyl.ylabel("Voltage [mV]")
    #if simpletrace._name:
    #    pyl.ylabel("{0} [{1}]".format(simpletrace._name, "mV"))

# -----------------------------------------------------------------------------
def plot_histogram(data, bin_borders, errors=None, *args, **kwargs):
# -----------------------------------------------------------------------------
    """
    Plot the function y(h) in the histogram shape.
    """
    assert len(data) + 1 == len(bin_borders), "Bad input."
    width = bin_borders[1] - bin_borders[0]
    pyl.bar(bin_borders[:-1], data, width, *args, **kwargs)
    if not errors == None:
        assert len(data) == len(errors)
        for i, error in enumerate(errors):
            x = bin_borders[i] + width/2.
            pyl.plot([x,x],[data[i]+error, data[i]-error], "r")

# -----------------------------------------------------------------------------
def plot_spike_population(spike_population, full_x=True):
# -----------------------------------------------------------------------------
    """
    Visulizes the spike trains of the population.
    """
    time = spike_population._total_time
    y_width = 1
    for i, train in enumerate(spike_population.spike_trains):
        y_position = np.zeros(len(train._spike_times)) + i * y_width
        pyl.plot(train._spike_times, y_position,".", color = "blue")
    #pyl.xlabel("t [ms]")
    #pyl.ylabel("N (neuron spike trains)")
    pyl.ylim((-y_width,(i + 1) * y_width))
    if full_x:
        pyl.xlim(0, time)

# -----------------------------------------------------------------------------
def plot_meshgrid(sweep_parameters_list, results, interpolation="nearest", bar_label=None, title=None, **kwargs):
# -----------------------------------------------------------------------------
    """
    sweep_parameters_list: a list of dictionaries {sweep_parameter_name:values}
    results: a list of results. Requires homogeneous spacing of sweep parameter values.
    """
    
    "TODO: Check that parameter values are homogeneously spaced out." 
    
    x = sweep_parameters_list[0][sweep_parameters_list[0].keys()[0]]
    y = sweep_parameters_list[1][sweep_parameters_list[1].keys()[0]]
    dims = (len(x),len(y))
    Z = np.array(results).reshape(dims)
    
    #X,Y = pyl.meshgrid(x,y)
    
    ax = pyl.subplot(111)
    y_d = (y[1]-y[0])/2.
    x_d = (x[1]-x[0])/2.
    im = pyl.imshow(Z, cmap=pyl.cm.jet, origin="lower", extent=(min(y)-y_d,max(y)+y_d,min(x)-x_d,max(x)+x_d), 
                    aspect="auto", **kwargs)
    im.set_interpolation(interpolation)
    cbar = pyl.colorbar()
    if bar_label:
        cbar.set_label(bar_label)
    pyl.xlabel(sweep_parameters_list[1].keys()[0])
    pyl.ylabel(sweep_parameters_list[0].keys()[0])
    if title:
        pyl.title(title)

