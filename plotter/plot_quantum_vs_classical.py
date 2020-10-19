import random

from bokeh.plotting import figure, show, output_file
from bokeh.palettes import Turbo256
from bokeh.models import Legend, LegendItem, ColumnDataSource, LabelSet, Label

def plot_q_vs_c(data):
    
    output_file("TTS comparison quantum vs random.html")
    plot_tts = figure(
        title='TTS comparison quantum vs random', 
        x_axis_type="log", 
        y_axis_type="log", 
        y_range=(1/5, 5),
        x_range=(1, 10**6)
        )

    # fix xaxis to cross y axis in the point (1,1)
    plot_tts.xaxis.fixed_location = 1

    # it stores data to create the labels
    x = []
    y = []
    point_names = []

    for protein_key in data:

        min_tts_q = data[protein_key]['min_tts_q']
        min_tts_c = data[protein_key]['min_tts_c']

        relation = min_tts_c / min_tts_q
        min_tts = min(min_tts_q, min_tts_c)

        protein_name = protein_key.split('_')[0]+protein_key.split('_')[1]+protein_key.split('_')[2][0]+'_'+protein_key.split('_')[3]+'_'+protein_key.split('_')[4]

        # save data to create the label of each point
        x.append(min_tts)
        y.append(relation)
        point_names.append(protein_name)

        # paint the point in the plot (only the point, the label is plotted out of the loop)
        plot_tts.square(min_tts, relation, line_width=1, line_color=random.sample(Turbo256,15))

    # add label of each point
    source = ColumnDataSource(dict(x=x,y=y,point_names=point_names))
    labels = LabelSet(x='x', y='y', text='point_names', level='glyph', x_offset=5, y_offset=5, text_font_size="6pt", source=source, render_mode='canvas')
    plot_tts.add_layout(labels)

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