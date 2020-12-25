import random
import re

from bokeh.plotting import figure, show, output_file
from bokeh.palettes import Turbo256
from bokeh.models import Legend, LegendItem, ColumnDataSource, LabelSet, Label

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
            plot_tts.square(min_tts, relation, size=size_point, fill_color=color_point, fill_alpha=0.6, line_color=color_point)

        elif data[protein_key]['number_aas'] == 3:

            plot_tts.circle(min_tts, relation, size=size_point, fill_color=color_point, fill_alpha=0.6, line_color=color_point)


        elif data[protein_key]['number_aas'] == 4:

            plot_tts.triangle(min_tts, relation, size=size_point, fill_color=color_point, fill_alpha=0.6, line_color=color_point)


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

    plot_q_c_slop = figure(
        title='Evolution of tts with different steps', 
        x_axis_label='Classical min(TTS)', 
        y_axis_label='Quantum min(TTS)',
        x_axis_type="log",
        y_axis_type="log",
        x_range=(10**2, 10**4), 
        y_range=(10**2, 10**4))

    x_point = []
    y_point = []
    line_color = []
    marker = []
    legend = []

    for protein_key in data:

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
    
        legend.append('minifold' if 'minifold' in protein_key else 'random')
        size.append(re.findall('_[0-9]_', protein_key)[0][1])

    source = ColumnDataSource(dict(x = x_point, y = y_point, line_color=line_color, marker=marker, legend=legend, size=size))
        #plot_q_c_slop.triangle(min_tts_c, min_tts_q, size=10, line_color='red', color='transparent')
            
        #plot_q_c_slop.circle(min_tts_c, min_tts_q, size=10, line_color='blue', color='transparent')
        

    plot_q_c_slop.scatter(x="x", y="y", size=10, line_color="line_color", fill_alpha=0, marker="marker", legend_group='legend', source=source)
    
    x_line = [1, 10**6]
    y_line = [1, 10**6]
    plot_q_c_slop.line(x_line, y_line, line_width=2, line_color='gray', line_dash="dashed")

    plot_q_c_slop.yaxis.major_label_orientation = "vertical"
    plot_q_c_slop.xgrid.grid_line_color = None
    plot_q_c_slop.ygrid.grid_line_color = None

    show(plot_q_c_slop)
