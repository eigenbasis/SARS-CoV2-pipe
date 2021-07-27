from bokeh.models import HoverTool, ColumnDataSource
from bokeh.plotting import figure
from bokeh.io import save, output_file
import sys
from operator import methodcaller

path_to_depth_report_file, output_file_path = sys.argv[1], sys.argv[2]


def bar_plot(depth_dict):
    '''The function is used to generate coverage depth plot from samtools depth report file.'''
    source = ColumnDataSource(data=depth_dict)

    p=figure(
        x_axis_label = "Genome position", 
        y_axis_label = 'Coverage depth', 
        plot_width=1750, 
        plot_height=500
    )
   
    hover = HoverTool(
    tooltips=[("Genome position", "@{Genome position}"), ('Coverage depth', "@{Coverage depth}")],
    attachment = 'right'
    )
    
    p.xaxis.major_label_orientation = "vertical"
    p.add_tools(hover)
    p.vbar(x='Genome position', top='Coverage depth',  width=0.5, bottom=0, fill_color=('green'), source=source, name = 'vbar')
    output_file(output_file_path)
    save(p, output_file_path)

with open(path_to_depth_report_file, "r+") as depth_file:
    data_dict = dict(list(map(lambda x:list(map(int,x[1:])), list(map(methodcaller("split", "\t"), depth_file.read().split('\n')))))[:-1])
    depth_dict = {'Genome position':list(data_dict.keys()), 'Coverage depth':list(data_dict.values())}

bar_plot(depth_dict)