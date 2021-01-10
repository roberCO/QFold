import random
import re
import numpy as np

from bokeh.plotting import figure, show, output_file
from bokeh.palettes import Turbo256
from bokeh.models import Legend, LegendItem, ColumnDataSource, LabelSet, Label
from collections import OrderedDict
from bokeh.palettes import Dark2_5 as palette
import itertools
from bokeh.models import SingleIntervalTicker, LinearAxis
import numpy as np

def plot_q_vs_c(data):
    
    output_file("TTS comparison quantum vs random.html")

    plot_tts = figure(
        title='TTS comparison quantum vs random', 
        x_axis_type="log", 
        y_axis_type="log", 
        y_range=(1/6, 6),
        x_range=(1, 10**6)
        )

    # fix xaxis to cross y axis in the point (1,1)
    plot_tts.xaxis.fixed_location = 1

    for protein_key in data:

        min_tts_q = data[protein_key]['min_tts_q']
        min_tts_c = data[protein_key]['min_tts_c']

        relation = min_tts_c / min_tts_q
        min_tts = min(min_tts_q, min_tts_c)

        if 'minifold' in protein_key:
            color_point = 'red'
        elif 'random' in protein_key:
            color_point = 'blue'
        elif 'original' in protein_key:
            color_point = 'green'
        else:
            color_point = 'black'

        size_point = int(data[protein_key]['number_bits'])**1.75

        # plot dipeptides
        if data[protein_key]['number_aas'] == 2:

            # paint the point in the plot (only the point, the label is plotted out of the loop)
            plot_tts.circle(min_tts, relation, size=size_point, fill_color=color_point, fill_alpha=0.6, line_color=color_point)

        elif data[protein_key]['number_aas'] == 3:

            plot_tts.triangle(min_tts, relation, size=size_point, fill_color=color_point, fill_alpha=0.6, line_color=color_point)

        elif data[protein_key]['number_aas'] == 4:

            plot_tts.square(min_tts, relation, size=size_point, fill_color=color_point, fill_alpha=0.6, line_color=color_point)


    # plot the plane for classical and quantum region
    plot_tts.quad(top=[10**10], bottom=[1], left=[1],right=[10**10], color="green", fill_alpha=0.1)
    plot_tts.quad(top=[1], bottom=[10**-10], left=[1], right=[10**10], color="red", fill_alpha=0.1)
    
    # create the label of each axis
    plot_tts.yaxis.axis_label='relation classical/quantum'
    # the x_label is created appart from the axis label due a bokeh bug
    x_label = Label(x=50, y=270, x_units='screen', y_units='screen', text=' TTS ', render_mode='css')

    # add a fake multiline figure to the plot. It is necessary to add a legend of the classical\quantum area
    r = plot_tts.multi_line([[0, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, 0]], color=["green", "red"], line_alpha=0.3, line_width=20)
    
    # add the legend of the fake multiline for the c\q regions
    legend = Legend(items=[
        LegendItem(label="quantum", renderers=[r], index=0),
        LegendItem(label="classical", renderers=[r], index=1),
    ])

    # add all legends to the plot
    plot_tts.add_layout(x_label, 'left')
    plot_tts.add_layout(legend, 'right')

    # position all legends
    plot_tts.legend.location = "center_right"
    plot_tts.title.align = 'center'

    show(plot_tts)

def plot_q_vs_c_slope(data):

    output_file("TTS slope quantum vs random.html")

    width = 700
    height = 450

    plot_q_c_slop = figure(
        #title='Evolution of tts with different steps', # Usually graphs do not have title
        x_axis_label='Classical min(TTS)', 
        y_axis_label='Quantum min(TTS)',
        x_axis_type="log",
        y_axis_type="log",
        x_range=(2*10**1, 3*10**4), 
        y_range=(2*10**1, 3*10**4),
        plot_height=height,
        plot_width=width)

    al = []
    bl = []

    for init in ['minifold', 'random']:

        x_point = []
        y_point = []
        line_color = []
        marker = []
        legend = []
        size = []


        for protein_key in data:

            if init in protein_key:

                x_point.append(data[protein_key]['min_tts_c'])
                y_point.append(data[protein_key]['min_tts_q'])
                line_color.append('red' if 'minifold' in protein_key else 'blue')

                if re.search("[A-Z]{4}", protein_key):
                    marker.append('square')
                elif re.search("[A-Z]{3}", protein_key):
                    marker.append('triangle')
                elif re.search("[A-Z]{2}", protein_key):
                    marker.append('circle')
                else:
                    raise ValueError('The string {} does not fit the dipeptide, tripeptide or tetrapeptide description.'.format(protein_key))

                # TODO change this
                legend.append('minifold' if 'minifold initialization' in protein_key else 'random initialization')
                size.append(str(int(re.findall('_[0-9]_', protein_key)[0][1])*2))

        logcx = np.log(x_point)
        logqy = np.log(y_point)

        model = np.polynomial.polynomial.polyfit(logcx, logqy, 1)
        logb, a = model
        b = np.exp(logb)

        al.append(a)
        bl.append(b)

        # 100 linearly spaced numbers
        x_fit = np.linspace(2e1, 3e4)

        # the function, which is y = x^2 here
        y_fit = b*x_fit**a
        print('a,b',a,b)

        source = ColumnDataSource(dict(x = x_point, y = y_point, line_color=line_color, marker=marker, legend=legend, size = size))

        #plot_q_c_slop.triangle(min_tts_c, min_tts_q, size=10, line_color='red', color='transparent')      
        #plot_q_c_slop.circle(min_tts_c, min_tts_q, size=10, line_color='blue', color='transparent')

        line_color = 'red' if init == 'minifold' else 'blue'    
        x_fit = list(x_fit)
        #fit_source = ColumnDataSource(dict(x = x_fit, y = y_fit, line_color='green', legend='y='+str(b)+'*x**'+str(a)))
        plot_q_c_slop.line(x_fit, y_fit, line_color=line_color) #
        plot_q_c_slop.scatter(x="x", y="y", line_color="line_color", fill_alpha=0, marker="marker", source=source, size = "size") #legend_group='legend',
    
    x_diag = [1, 10**6]
    y_diag = [1, 10**6]
    plot_q_c_slop.line(x_diag, y_diag, line_width=2, line_color='gray', line_dash="dashed")

    plot_q_c_slop.yaxis.major_label_orientation = "vertical"
    plot_q_c_slop.xgrid.grid_line_color = None
    plot_q_c_slop.ygrid.grid_line_color = None

    x = -1
    y = -1
    rdb = plot_q_c_slop.diamond(x, y, color = 'blue')
    rlb = plot_q_c_slop.line(x, y, color = 'blue')

    rdr = plot_q_c_slop.diamond(x, y, color = 'red')
    rlr = plot_q_c_slop.line(x, y, color = 'red')

    rcg = plot_q_c_slop.circle(x, y, line_color="black")
    rtg = plot_q_c_slop.triangle(x, y, line_color="black")
    rsg = plot_q_c_slop.square(x, y, line_color="black")

    rdg2 = plot_q_c_slop.diamond(x, y, line_color="black", size = 4)
    rdg3 = plot_q_c_slop.diamond(x, y, line_color="black", size = 6)
    rdg4 = plot_q_c_slop.diamond(x, y, line_color="black", size = 8)
    rdg5 = plot_q_c_slop.diamond(x, y, line_color="black", size = 10)

    legend = Legend(items=[
        ("random initialization", [rdb, rlb]),
        (str(np.round(bl[1],3))+'*x**'+str(np.round(al[1],3)), [rlb]),
        ("minifold initialization", [rdr, rlr]),
        (str(np.round(bl[0],3))+'*x**'+str(np.round(al[0],3)), [rlr]),

        ("dipeptides", [rcg]),
        ("tripeptides", [rtg]),
        ("tetrapeptides", [rsg]),

        ("2 bits", [rdg2]),
        ("3 bits", [rdg3]),
        ("4 bits", [rdg4]),
        ("5 bits", [rdg5])
    ], location=(.35*width, .28*height), background_fill_alpha = 0, border_line_alpha = 0)

    plot_q_c_slop.add_layout(legend, 'left')

    citation = Label(x=1/3*width, y=1/5*height, x_units='screen', y_units='screen',
                    text='Quantum advantage', render_mode='canvas',
                    border_line_color='black', border_line_alpha=0.0,
                    background_fill_color='white', background_fill_alpha=0.0)

    plot_q_c_slop.add_layout(citation)

    citation = Label(x=1/4*width, y=2/3*height, x_units='screen', y_units='screen',
                    text='Classical advantage', render_mode='canvas',
                    border_line_color='black', border_line_alpha=0.0,
                    background_fill_color='white', background_fill_alpha=0.0)

    plot_q_c_slop.add_layout(citation)



    show(plot_q_c_slop)

def plot_q_opt_step(data):

    output_file("Quantum optimal step.html")
    quantum_optimal_step = {}

    for protein_key in data:

        protein_name = protein_key.split('_')[0]
        bits = protein_key.split('_')[1]
        init_method = protein_key.split('_')[2]
        step = data[protein_key]['min_tts_q_step']

        if len(protein_name) == 2:

            if protein_name+'-'+init_method in quantum_optimal_step.keys():
                quantum_optimal_step[protein_name+'-'+init_method][bits] = step
            else:
                quantum_optimal_step[protein_name+'-'+init_method] = {bits: step}

    p = figure(plot_width=400, plot_height=400, x_range=(3, 10), x_axis_type=None)
    ticker = SingleIntervalTicker(interval=1, num_minor_ticks=10)
    xaxis = LinearAxis(ticker=ticker)
    p.add_layout(xaxis, 'below')
    # create a color iterator
    colors = itertools.cycle(palette)    
    
    for protein, color in zip(quantum_optimal_step, colors):

        ordered_quantum_steps = OrderedDict(sorted(quantum_optimal_step[protein].items()))


        p.line(list(ordered_quantum_steps.keys()), list(ordered_quantum_steps.values()), line_width=2, color=color, legend=protein)

        #p.circle([1, 2], [3, 4], size=20, color="navy", alpha=0.5)


    show(p)