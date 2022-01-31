from bokeh.io import output_file, show
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
from bokeh.models import FactorRange
from bokeh.layouts import row

def plot_m_vs_r(data):

    protein_keys = data.keys()

    precision_minifold = []
    precision_random = [] 
    tts_minifold = []
    tts_random =  []

    color_precision = []
    color_tts = []

    proteins = []

    for protein_key in protein_keys:
        
        protein = protein_key.split('_')[0]+'_'+protein_key.split('_')[1]+'_'+protein_key.split('_')[3]+'_'+protein_key.split('_')[4]
        if not protein in proteins:
            proteins.append(protein)

        if 'minifold' in protein_key:
            precision_minifold.append(data[protein_key]['precision'])
            tts_minifold.append((-1)*data[protein_key]['min_tts'])

        elif 'random' in protein_key:
            precision_random.append(data[protein_key]['precision'])
            tts_random.append((-1)*data[protein_key]['min_tts'])


    output_file("random_vs_minifold.html")

    metrics = ["random", "minifold", "random", "minifold"]

    y_precision = [ (protein, metric) for protein in proteins for metric in metrics[:2]]
    y_tts = [ (protein, metric) for protein in proteins for metric in metrics[2:]]

    counts_precision = sum(zip(precision_random, precision_minifold), ())
    counts_tts = sum(zip(tts_random, tts_minifold), ())

    for _ in proteins:

        # append color for minifold precision
        color_precision.append('#219c2b')
        # append color for random precision
        color_precision.append('#42f551')

        # append color for minifold tts
        color_tts.append('#ad1717')
        # append color for random tts
        color_tts.append('#fc2b2b')
        
    source_precision = ColumnDataSource(data=dict(y=y_precision, counts=counts_precision, color=color_precision))
    source_tts = ColumnDataSource(data=dict(y=y_tts, counts=counts_tts, color=color_tts))

    plot_precision = figure(x_axis_label="Precision", x_range=(0, 1), y_range=FactorRange(*y_precision), plot_height=900, title="Precision comparison  minifold vs random",
                toolbar_location=None)
    plot_tts = figure(x_axis_label="TTS", x_range=(-1200, 0), y_range=FactorRange(*y_tts), plot_height=900, title="TTS comparison minifold vs random",
                toolbar_location=None)

    #p.hbar_stack(metrics, y='proteins', height=0.9, color='#008000', source=ColumnDataSource(precision), legend_label= precision_names)
    plot_precision.hbar(y='y', right='counts', height=0.9, fill_color='color', source=source_precision)

    #p.hbar_stack(metrics, y='proteins', height=0.9, color='#ff0000', source=ColumnDataSource(tts), legend_label=tts_names)
    plot_tts.hbar(y='y', right='counts', height=0.9, fill_color='color', source=source_tts)


    plot_precision.y_range.range_padding = 0.1
    plot_precision.ygrid.grid_line_color = None
    plot_precision.axis.minor_tick_line_color = None
    plot_precision.outline_line_color = None

    plot_tts.y_range.range_padding = 0.1
    plot_tts.ygrid.grid_line_color = None
    plot_tts.axis.minor_tick_line_color = None
    plot_tts.outline_line_color = None
    plot_tts.yaxis.visible = False

    plot_precision.title.align = 'center'
    plot_tts.title.align = 'center'



    show(row(plot_tts, plot_precision))