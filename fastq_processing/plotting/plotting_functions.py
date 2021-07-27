from bokeh.models import HoverTool, ColumnDataSource, Legend
from bokeh.plotting import figure
from bokeh.palettes import magma

x_label = 'Parauga ņemšanas datums'
y_label = 'Paraugu skaits'
plot_width = 1750
plot_height = 1000

def bar_plot(data_dict, unq_dict, x_range_name):
    '''The function is used to generate bar plot of pre-defined format
    from data_dict where number of cases are mapped to sampling date.'''
    
    source = ColumnDataSource(data=data_dict)
    p=figure(
        x_range=unq_dict[x_range_name],
        x_axis_label = x_label,
        y_axis_label = y_label,
        plot_width=plot_width,
        plot_height=plot_height 
    )
    
    hover = HoverTool(
    tooltips=[(x_label, "@sampling_date"), (y_label, "@number_of_cases")],
    attachment = 'right'
    )
    
    p.xaxis.major_label_orientation = "vertical"
    p.add_tools(hover) 
    p.vbar(x='sampling_date', top='number_of_cases',  width=0.5, bottom=0, fill_color=('green'), source=source)

    return p

def stacked_bar_plot(data_dict, unq_dict, strata='lineage', title='', units='abs'):
    '''The function is used to generate stacked bar plot from 
    data_dict where data are stratified by lineage and by date.'''
    
    source = ColumnDataSource(data=data_dict)
    legend = sorted(list(data_dict.keys()))
    legend.remove('sampling_date')
    legend.remove('total_cases')

    if units == 'rel': 
        y_label = 'Paraugu daļa (%)'
        total_label = 'Paraugu skaits'
    else: 
        y_label = total_label = 'Paraugu skaits'
 
    palette = magma(len(legend))    

    p=figure(title= title, x_range=unq_dict['sampling_date'], x_axis_label = x_label, y_axis_label = y_label, 
        plot_width=plot_width, plot_height=plot_height, toolbar_location = 'below', name='figure')
    
    p.add_layout(Legend(name = 'legend'), 'right')
    p.vbar_stack(legend, x='sampling_date', width=0.9, source=source, name = 'background', color = "#e9ecf0", line_color = '#acaeb0')
    
    renderers = p.vbar_stack(legend, x='sampling_date', width=0.9, source=source, legend_label=legend, name = legend, color = palette)
   
    for r in renderers:
        name_1 = r.name
        hover = HoverTool(tooltips=[
                (x_label,'@sampling_date'),
                (strata, f'{name_1}'),
                (y_label, '@$name'),
                (f'{total_label} dienā', '@total_cases')], 
                renderers=[r])
        hover.toggleable = False
        p.add_tools(hover)
    
    p.legend.click_policy="hide"
    p.xaxis.major_label_orientation = "vertical" 

    return p