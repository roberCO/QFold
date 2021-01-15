import random
import re
import numpy as np
import itertools
from collections import OrderedDict

from bokeh.plotting import figure, show, output_file, gridplot
from bokeh.io import export_png, export_svgs
from bokeh.palettes import Turbo256
from bokeh.models import Legend, LegendItem, ColumnDataSource, LabelSet, Label, Marker, renderers
from bokeh.palettes import Dark2_5 as palette
from bokeh.models import SingleIntervalTicker, LinearAxis, Slope
from bokeh.models.tickers import FixedTicker


def plot_q_vs_c(data):
    
    output_file("TTS comparison quantum vs random.html")

    width = 800
    height = 450

    x_range = (2**5, 2**11)
    x_min, x_max = x_range
    plot_tts = figure(
        x_axis_type="log", 
        y_axis_type="log", 
        y_range=(2*10**1, 3*10**4),
        x_range=(x_min, x_max),
        plot_height=height,
        plot_width=width
        )

    plot_tts.xaxis.ticker = [2**6, 2**8, 2**10]

    # fix xaxis to cross y axis in the point (1,1)
    #plot_tts.xaxis.fixed_location = 1

    x_point = []
    c_point = []
    q_point = []
    minirandom = []

    # filter by fixed beta
    data = dict(filter(lambda item: item[1]['schedule'] == 'fixed', data.items()))

    for protein_key in data:

        min_tts_q = data[protein_key]['min_tts_q']
        min_tts_c = data[protein_key]['min_tts_c']

        c_point.append(min_tts_c)
        q_point.append(min_tts_q)

        #relation = min_tts_c / min_tts_q
        #min_tts = min(min_tts_q, min_tts_c)

        if 'minifold' in protein_key:
            color_point = 'red'
            minirandom.append('minifold')
        elif 'random' in protein_key:
            color_point = 'blue'
            minirandom.append('random')
        elif 'original' in protein_key:
            color_point = 'green'
        else:
            color_point = 'black'

        size_point = int(data[protein_key]['number_bits'])*3
        space_size = 2**(int(data[protein_key]['number_bits']) * (2*int(data[protein_key]['number_aas']) -2))
        x_point.append(space_size)

        # plot dipeptides
        if data[protein_key]['number_aas'] == 2:

            # paint the point in the plot (only the point, the label is plotted out of the loop)
            #plot_tts.circle(min_tts, relation, size=size_point, fill_color=color_point, fill_alpha=0.6, line_color=color_point)
            plot_tts.circle(space_size, min_tts_q, size=size_point, fill_color=color_point, fill_alpha=1, line_color=color_point)
            plot_tts.circle(space_size, min_tts_c, size=size_point, line_color=color_point, fill_alpha=0, fill_color='white')

        elif data[protein_key]['number_aas'] == 3:

            #plot_tts.triangle(min_tts, relation, size=size_point, fill_color=color_point, fill_alpha=0.6, line_color=color_point)
            plot_tts.triangle(space_size, min_tts_q, size=size_point, fill_color=color_point, fill_alpha=1, line_color=color_point)
            plot_tts.triangle(space_size, min_tts_c, size=size_point, line_color=color_point, fill_alpha=0, fill_color='white')

        elif data[protein_key]['number_aas'] == 4:

            #plot_tts.square(min_tts, relation, size=size_point, fill_color=color_point, fill_alpha=0.6, line_color=color_point)
            plot_tts.square(space_size, min_tts_q, size=size_point, fill_color=color_point, fill_alpha=1, line_color=color_point)
            plot_tts.square(space_size, min_tts_c, size=size_point, line_color=color_point, fill_alpha=0, fill_color='white')

    # plot the plane for classical and quantum region
    #plot_tts.quad(top=[10**10], bottom=[1], left=[1],right=[10**10], color="green", fill_alpha=0.1)
    #plot_tts.quad(top=[1], bottom=[10**-10], left=[1], right=[10**10], color="red", fill_alpha=0.1)
    #'''
    x_fit_rand = []
    q_fit_rand = []
    c_fit_rand = []
    x_fit_mini = []
    q_fit_mini = []
    c_fit_mini = []

    assert(len(x_point) == len(q_point))
    for i in range(len(x_point)):
        if minirandom[i] == 'random':
            x_fit_rand.append(x_point[i])
            q_fit_rand.append(q_point[i])
            c_fit_rand.append(c_point[i])
        elif minirandom[i] == 'minifold':
            x_fit_mini.append(x_point[i])
            q_fit_mini.append(q_point[i])
            c_fit_mini.append(c_point[i])
        
    logx_fit_rand = np.log(x_fit_rand)
    logq_fit_rand = np.log(q_fit_rand)
    logc_fit_rand = np.log(c_fit_rand)
    logx_fit_mini = np.log(x_fit_mini)
    logq_fit_mini = np.log(q_fit_mini)
    logc_fit_mini = np.log(c_fit_mini)
    

    modelqr = np.polynomial.polynomial.polyfit(logx_fit_rand, logq_fit_rand, 1)
    modelqm = np.polynomial.polynomial.polyfit(logx_fit_mini, logq_fit_mini, 1)
    modelcr = np.polynomial.polynomial.polyfit(logx_fit_rand, logc_fit_rand, 1)
    modelcm = np.polynomial.polynomial.polyfit(logx_fit_mini, logc_fit_mini, 1)
    models = [modelcm, modelcr, modelqm, modelqr]

    model_renders = {}
        
    for model, mod, ini in zip(models, ['c', 'c', 'q', 'q'], ['m', 'r', 'm', 'r']):
        logb, a = model
        b = np.exp(logb)
        # 100 linearly spaced numbers
        x_min, x_max = x_range
        x_fit = np.linspace(x_min, x_max)

        # the function, which is y = x^2 here
        y_fit = b*x_fit**a
        print('a,b',a,b)
        line_color = 'red' if 'm' in ini else 'blue'
        line_dash = 'solid' if 'q' in mod else 'dashed'  
        x_fit = list(x_fit)

        plot_tts.line(x_fit, y_fit, line_color=line_color, line_dash = line_dash)

        x = -1
        y = -1
        model_renders[(mod,ini)] = plot_tts.line(x, y, line_color=line_color, line_dash = line_dash)

        y0 = 0
        if 'c' in mod:
            y0 += 280
        if 'r' in ini:
            y0 -= 25

        background_fill_alpha = 0.1 if line_dash == 'solid' else 0

        citation = Label(x=270, y=80 + y0, x_units='screen', y_units='screen',
            text='TTS = '+str(np.round(b, 3))+'* size**'+str(np.round(a, 3)), text_font_size='12px', render_mode='canvas',
            border_line_color=line_color, border_line_alpha=0.0, border_line_dash = line_dash,
            background_fill_color=line_color, background_fill_alpha=background_fill_alpha, text_color = line_color)

        plot_tts.add_layout(citation)

    red_circle = plot_tts.circle(x, y, color = 'red')
    blue_circle = plot_tts.circle(x, y, color = 'blue')

    red_circumphere = plot_tts.circle(x, y, line_color = 'red', color = 'white', fill_alpha = 0)
    blue_circumphere = plot_tts.circle(x, y, line_color = 'blue', color = 'white', fill_alpha= 0)
    #'''
    
    # create the label of each axis
    plot_tts.yaxis.axis_label='TTS'
    plot_tts.xaxis.axis_label='Size of the configuration space'

    legend = Legend(items=[
        ("minifold initialization, quantum", [model_renders[('q','m')], red_circle]),
        ("random initialization, quantum", [model_renders[('q','r')], blue_circle]),
        ("minifold initialization, classical", [model_renders[('c','m')], red_circumphere]),
        ("random initialization, classical", [model_renders[('c','r')], blue_circumphere])
    ], background_fill_alpha = 0, border_line_alpha = 0, location=(.37*width, .6*height))

    # add all legends to the plot
    plot_tts.add_layout(legend, 'left')

    # position all legends
    #plot_tts.legend.location = "top_left"
    #plot_tts.title.align = 'center'

    plot_tts.output_backend = "svg"

    export_svgs(plot_tts, filename="size_vs_TTS_slope.svg")

    show(plot_tts)

def plot_q_vs_c_tts_ratio(data):
    
    output_file("TTS comparison quantum vs random.html")

    data = dict(filter(lambda item: item[1]['schedule'] == 'fixed', data.items()))

    width = 800
    height = 450

    plot_tts = figure( 
        x_axis_type="log", 
        y_axis_type="log", 
        y_range=(1/3, 50),
        #x_range=(20, 1e4),
        plot_height=height,
        plot_width=width,
        )

    # fix xaxis to cross y axis in the point (1,1)
    plot_tts.xaxis.fixed_location = 1

    for protein_key in data:

        min_tts_q = data[protein_key]['min_tts_q']
        min_tts_c = data[protein_key]['min_tts_c']
        schedule = data[protein_key]['schedule']
        initializer = data[protein_key]['initializer']
        space_size = 2**(int(data[protein_key]['number_bits']) * (2*int(data[protein_key]['number_aas']) -2))

        if schedule == 'fixed':
            color_point = 'blue'
        elif schedule == 'Boltzmann' or schedule == 'logarithmic':
            color_point = 'green'
        elif schedule == 'Cauchy' or schedule == 'linear':
            color_point = 'red'
        elif schedule == 'exponential':
            color_point = 'orange'
        elif schedule == 'geometric':
            color_point = 'black'
        else:
            raise ValueError(schedule)

        if initializer == 'minifold':
            fill_alpha = .5
        elif initializer == 'random':
            fill_alpha = 0
        else:
            raise ValueError(initializer)

        relation = min_tts_c / min_tts_q
        min_tts = min(min_tts_q, min_tts_c)

        size_point = int(data[protein_key]['number_bits'])*3

        # plot dipeptides
        if data[protein_key]['number_aas'] == 2:

            plot_tts.circle(min_tts_q, relation, size=size_point, fill_color=color_point, fill_alpha=fill_alpha, line_color=color_point)

        elif data[protein_key]['number_aas'] == 3:

            plot_tts.triangle(min_tts_q, relation, size=size_point, fill_color=color_point, fill_alpha=fill_alpha, line_color=color_point)

        elif data[protein_key]['number_aas'] == 4:

            plot_tts.square(min_tts_q, relation, size=size_point, fill_color=color_point, fill_alpha=fill_alpha, line_color=color_point)

    # plot the plane for classical and quantum region
    plot_tts.quad(top=[10**10], bottom=[1], left=[1],right=[10**10], color="#003366", fill_alpha=0.05)
    plot_tts.quad(top=[1], bottom=[10**-10], left=[1], right=[10**10], color="#996633", fill_alpha=0.05)
    
    # create the label of each axis
    plot_tts.yaxis.axis_label='ratio classical/quantum min TTS'
    plot_tts.xaxis.axis_label='min(classical min TTS, quantum min TTS)'
    # the x_label is created appart from the axis label due a bokeh bug
    x_label = Label(x=100, y=70, x_units='screen', y_units='screen', text='quantum TTS', text_font_style = 'italic', text_font_size = "9pt", render_mode='canvas')

    # add a fake multiline figure to the plot. It is necessary to add a legend of the classical\quantum area
    r = plot_tts.multi_line([[0, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, 0]], color=["#003366", "#996633"], line_alpha=0.3, line_width=20)
    
    r_dot = plot_tts.dot(0, 0, line_color="black")
    r_circle = plot_tts.circle(0, 0, line_color="black", fill_alpha = 0)
    r_triangle = plot_tts.triangle(0, 0, line_color="black", fill_alpha = 0)
    r_square = plot_tts.square(0, 0, line_color="black", fill_alpha = 0)
    r_diamond = plot_tts.diamond(0, 0, line_color="black", fill_alpha = 0)

    r_hex_color = plot_tts.hex(0, 0, line_color=['black', 'green', 'red', 'orange', 'blue'], fill_alpha = 0)

    r_hex_init = plot_tts.hex(0, 0, color='blue', fill_alpha = [0,1])

    r_bits = plot_tts.hex(0, 0, line_color="white", size = [3, 6, 9, 12, 15], fill_alpha = 0)

    # add the legend of the fake multiline for the c\q regions
    legend = Legend(items=[
        LegendItem(label="qTTS < cTTS", renderers=[r], index=0),
        LegendItem(label="cTTS < qTTS", renderers=[r], index=1),
        #LegendItem(label="fixed schedule", renderers=[r_hex_color], index=0),
        #LegendItem(label="Boltzmann schedule", renderers=[r_hex_color], index=1),
        #LegendItem(label="Cauchy schedule", renderers=[r_hex_color], index=2),
        #LegendItem(label="exponential schedule", renderers=[r_hex_color], index=3),
        #LegendItem(label="geometric schedule", renderers=[r_hex_color], index=4),
        LegendItem(label="random initialization", renderers=[r_hex_init], index=0),
        LegendItem(label="minifold initialization", renderers=[r_hex_init], index=1),
        LegendItem(label="dipeptides", renderers=[r_circle]),
        LegendItem(label="tripeptides", renderers=[r_triangle]),
        LegendItem(label="tetrapeptides", renderers=[r_square]),
        LegendItem(label="size = # rotation bits", renderers = [r_bits])
    ])

    # add all legends to the plot
    plot_tts.add_layout(x_label, 'left')
    plot_tts.add_layout(legend, 'right')

    # position all legends
    plot_tts.legend.location = "center_right"
    plot_tts.title.align = 'center'

    plot_tts.output_backend = "svg"
    export_svgs(plot_tts, filename="TTS_ratio.svg")

    show(plot_tts)


def TTSplotter(data, schedule, width = 800, height = 450, title = None):
    # filter by fixed beta
    data = dict(filter(lambda item: item[1]['schedule'] == schedule, data.items()))

    x_range = (2*10**1, 3*10**4)

    plot_q_c_slop = figure(
        #title='Evolution of tts with different steps', # Usually graphs do not have title
        x_axis_type="log",
        y_axis_type="log",
        x_range= x_range, 
        y_range=(2*10**1, 3*10**4),
        plot_height=height,
        plot_width=width,
        title = title)

    plot_q_c_slop.background_fill_alpha = 0

    if schedule == 'Boltzmann' or schedule == 'logarithmic':
        plot_q_c_slop.yaxis.axis_label = 'Quantum min(TTS)'
        plot_q_c_slop.yaxis.axis_label_text_font_size = "15pt"
    elif schedule == 'fixed':
        plot_q_c_slop.yaxis.axis_label = 'Quantum min(TTS)'
    else:
        plot_q_c_slop.yaxis.axis_label_text_font_size = "0pt"
        plot_q_c_slop.yaxis.major_label_text_font_size = "0pt"
        #plot_q_c_slop.yaxis.ticker = FixedTicker(ticks=[])

    plot_q_c_slop.xaxis.axis_label = 'Classical min(TTS)'

    if schedule != 'fixed':
        plot_q_c_slop.xaxis.axis_label_text_font_size = "15pt"
        plot_q_c_slop.xaxis.major_label_text_font_size = "15pt"
        plot_q_c_slop.title.text_font_size = '15pt'

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

                legend.append('minifold' if 'minifold initialization' in protein_key else 'random initialization')
                size.append(str(int(re.findall('_[0-9]_', protein_key)[0][1])*4)) if schedule != 'fixed' else size.append(str(int(re.findall('_[0-9]_', protein_key)[0][1])*3))

        logcx = np.log(x_point)
        logqy = np.log(y_point)

        model = np.polynomial.polynomial.polyfit(logcx, logqy, 1)
        logb, a = model
        b = np.exp(logb)

        al.append(a)
        bl.append(b)

        # 100 linearly spaced numbers
        x_fit = np.linspace(x_range[0]/1e2, x_range[1]*1e2)

        # the function, which is y = x^2 here
        y_fit = b*x_fit**a
        print('a,b',a,b)

        source = ColumnDataSource(dict(x = x_point, y = y_point, line_color=line_color, marker=marker, legend=legend, size = size))

        #plot_q_c_slop.triangle(min_tts_c, min_tts_q, size=10, line_color='red', color='transparent')      
        #plot_q_c_slop.circle(min_tts_c, min_tts_q, size=10, line_color='blue', color='transparent')

        line_color = 'red' if init == 'minifold' else 'blue'    
        x_fit = list(x_fit)
        #fit_source = ColumnDataSource(dict(x = x_fit, y = y_fit, line_color='green', legend='y='+str(b)+'*x**'+str(a)))
        plot_q_c_slop.line(x_fit, y_fit, line_color=line_color) 
        # Failed attempt to picture an slope: we do not want a line, but a line when we are in log scale
        # slope = Slope(gradient=a, y_intercept=b, line_color=line_color, line_dash='dashed', line_width=3.5)
        # plot_q_c_slop.add_layout(slope)

        text_font_size = '16px' if schedule != 'fixed' else '12px'

        y0 = 25 if init == 'minifold' else 0
        x0 = 300 if (schedule == 'Cauchy' or schedule == 'linear' or schedule == 'exponential' or schedule == 'geometric') else 0
        citation = Label(x=width-600 + x0, y=20 + y0, x_units='screen', y_units='screen',
            text='qTTS = '+str(np.round(b, 3))+'*cTTS^'+str(np.round(a, 3)), text_font_size= text_font_size, render_mode='canvas',
            border_line_color=line_color, border_line_alpha=0.0, border_line_dash = 'solid',
            background_fill_color=line_color, background_fill_alpha=0, text_color = line_color)

        plot_q_c_slop.add_layout(citation)

        plot_q_c_slop.scatter(x="x", y="y", line_color="line_color", fill_alpha=0, marker="marker", source=source, size = "size") #legend_group='legend',

    x_diag = [1, 10**6]
    y_diag = [1, 10**6]
    plot_q_c_slop.line(x_diag, y_diag, line_width=2, line_color='gray', line_dash="dashed")

    plot_q_c_slop.yaxis.major_label_orientation = "vertical"
    plot_q_c_slop.xgrid.grid_line_color = None
    plot_q_c_slop.ygrid.grid_line_color = None

    return plot_q_c_slop, al, bl

def generate_legend(plot, x0, y0, position = True, schedule = 'fixed'):
    x = -1
    y = -1
    rdb = plot.diamond(x, y, color = 'blue')
    rlb = plot.line(x, y, color = 'blue')

    rdr = plot.diamond(x, y, color = 'red')
    rlr = plot.line(x, y, color = 'red')

    rcg = plot.circle(x, y, line_color="black", fill_alpha = 0)
    rtg = plot.triangle(x, y, line_color="black", fill_alpha = 0)
    rsg = plot.square(x, y, line_color="black", fill_alpha = 0)

    rdg2 = plot.diamond(x, y, line_color="black", size = 4, fill_alpha = 0)
    rdg3 = plot.diamond(x, y, line_color="black", size = 6, fill_alpha = 0)
    rdg4 = plot.diamond(x, y, line_color="black", size = 8, fill_alpha = 0)
    rdg5 = plot.diamond(x, y, line_color="black", size = 10, fill_alpha = 0)

    if position == True:
        location = (x0, y0)
    else:
        location = 'top_left'

    legend = Legend(items=[
        ("random initialization", [rdb, rlb]),
        #(r'qTTS ='+str(np.round(bl[1],3))+r'*cTTS^'+exp1, [rlb]),
        ("minifold initialization", [rdr, rlr]),
        #(r'qTTS ='+str(np.round(bl[0],3))+r'*cTTS^'+exp0, [rlr]),
        ("size = # rotation bits", []),

        ("dipeptides", [rcg]),
        ("tripeptides", [rtg]),
        ("tetrapeptides", [rsg])
    ], location=location, background_fill_alpha = 0, border_line_alpha = 0)

    plot.add_layout(legend, 'left')
    if schedule != 'fixed':
        plot.legend.label_text_font_size = '15pt'

    return plot


def plot_q_vs_c_slope(data):

    output_file("TTS slope quantum vs random.html")
    width = 800
    height = 450

    plot_q_c_slop, al, bl = TTSplotter(data, schedule = 'fixed', width = width, height = height)

    plot_q_c_slop = generate_legend(plot_q_c_slop, .3*width, .4*height)

    citation = Label(x=1/3*width, y=1/7*height, x_units='screen', y_units='screen',
                    text='Quantum advantage', render_mode='canvas', text_font_size = '12px',
                    border_line_color='black', border_line_alpha=0.0,
                    background_fill_color='white', background_fill_alpha=0.0)

    plot_q_c_slop.add_layout(citation)

    citation = Label(x=1/4*width, y=5/6*height, x_units='screen', y_units='screen',
                    text='Classical advantage', render_mode='canvas', text_font_size = '12px',
                    border_line_color='black', border_line_alpha=0.0,
                    background_fill_color='white', background_fill_alpha=0.0)

    plot_q_c_slop.add_layout(citation)

    show(plot_q_c_slop)

    plot_q_c_slop.output_backend = "svg"
    export_svgs(plot_q_c_slop, filename="fixed_beta_TTS_slope.svg")

def plot_q_vs_c_slope_var(data):

    output_file("schedules.html")

    width = 750
    height = 500
    width2 = 450

    # create a new plot
    s1, al1, bl1 = TTSplotter(data, schedule = 'logarithmic', width = width, height = height, title = 'Boltzmann/logarithmic')
    s1 = generate_legend(s1, .4*width, .4*height, position = True, schedule = 'Boltzmann')
    s1.output_backend = "svg"
    export_svgs(s1, filename="Boltzmann_beta_TTS_slope.svg")

    # create another one
    s2, al2, bl2 = TTSplotter(data, schedule = 'linear', width = width2, height = height, title = 'Cauchy/linear')
    s2.output_backend = "svg"
    export_svgs(s2, filename="Cauchy_beta_TTS_slope.svg")

    # create another
    s3, al3, bl3 = TTSplotter(data, schedule = 'geometric', width = width2, height = height, title = 'Geometric')
    s3.output_backend = "svg"
    export_svgs(s3, filename="geometric_beta_TTS_slope.svg")

    # create final one
    s4, al4, bl4 = TTSplotter(data, schedule = 'exponential', width = width2, height = height, title = 'Exponential')
    s4.output_backend = "svg"
    export_svgs(s4, filename="exponential_beta_TTS_slope.svg")

    # put all the plots in a grid layout
    p = gridplot([[s1, s2, s3, s4]])

    #p.output_backend = "svg"
    export_png(p, filename="var_beta_TTS_slope.png")

    # show the results
    #show(p)

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