import matplotlib as mpl
mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble='\usepackage{color}')


font = {'family' : 'serif',
        'serif':['Palatino'],
        'weight' : 'bold',
        'size'   : 16}

colors = ['tomato', 'RoyalBlue', 'SeaGreen', 'Sienna', "MediumVioletRed", "yellow"]
mpl.rc('font', **font)
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['axes.color_cycle'] = colors
