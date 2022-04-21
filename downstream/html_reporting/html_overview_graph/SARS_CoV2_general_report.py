#import libraries

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import os
import time

#take files from the same directory
__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

#get current time for file names
timestr = time.strftime("%Y%m%d-%H%M%S")
datestr = time.strftime("%Y%m%d")
datestrread = time.strftime("%Y-%m-%d")

#importing and prepping the first data
file = os.path.join(__location__, 'filtered_mutations.csv')

#path = f"C:\Users\elza\Documents\rakus\filtered_mutations.csv"
df = pd.read_csv(file)

#transform the date so that plotly can comprehend it and plot
df = df.drop_duplicates()
df = df.dropna()

def df_to_plotly(df):
    return {'z': df["Fill"], 
            'x': df["sampling_date"], 
            'y': df["Mutation"] }

fig = go.Figure(data=go.Heatmap(df_to_plotly(df)))
fig.update_layout(title="Mutāciju apkopojums",
                  yaxis={"title": 'Mutācija'},
                  xaxis={"title": 'Paraugu ņemšanas datums',"tickangle": 45},
                  yaxis_nticks=len(df["Mutation"]))
fig.update_yaxes(showticklabels=False)

#importing and prepping the second data

file_area = os.path.join(__location__, 'summary_file_02_02_2022.csv')
df_area = pd.read_csv(file_area)
dates_to_exclude = ['2021-03-29', '0', '', 'NaT']
df_area_clean = df_area.loc[~df_area['sampling_date'].isin(dates_to_exclude)]
df_area_clean['sampling_date'] =  pd.to_datetime(df_area_clean['sampling_date'], format='%Y-%m-%d')
df_area_clean = df_area_clean.loc[(df_area_clean['sampling_date'] >= '2021-01-01')]
df_area_clean = df_area_clean.loc[df_area_clean['lineage'] != 'None']

#group the data for plotting

df_area_clean['week_start'] = df_area_clean['sampling_date'].dt.to_period('W').apply(lambda r: r.start_time)
df_area_summary = df_area_clean.groupby( [ "week_start", "lineage"] ).size().to_frame(name = 'count').reset_index()
df_area_summary = df_area_summary.sort_values(['lineage', 'week_start'], ascending=[True, True])

fig2 = px.area(df_area_summary, x="week_start", y="count", color="lineage", groupnorm='fraction', labels={ # replaces default labels by column name
                "week_start": "Nedēļa",  "count": "Paraugu daļa", "lineage": "Paveidi"
            },)

#building the bar chart
df_bar = df_area_summary
df_bar['weekly_relative'] = df_bar.groupby('week_start').transform(lambda x: x/x.sum())
df_bar['max_relative'] = df_bar.groupby('lineage')['weekly_relative'].transform('max')
df_bar['palette'] = np.where(df_bar['max_relative']>= 0.5, True, False)

fig3 = px.bar(df_bar, x='week_start', y='weekly_relative', color='lineage', barmode='relative', labels={ # replaces default labels by column name
                "week_start": "Nedēļa",  "count": "Paraugu daļa", "lineage": "Paveidi"
            },)

#put the bar chart behind the area chart
fig4 = px.area(df_bar, x='week_start', y='weekly_relative'
            ,color = 'lineage'
            #, line_group='lineage', color='max_relative'
            , color_discrete_sequence=px.colors.qualitative.Alphabet
            , labels={
                "week_start": "Nedēļa",  "weekly_relative": "Paraugu daļa", "lineage": "Paveidi"
            })
#fig4.for_each_trace(lambda trace: trace.update(fillcolor = trace.line.color))
#fig4.update_traces(name=df_bar['lineage'], showlegend = True)
fig4.add_bar(x=df_bar['week_start'], y=df_bar['weekly_relative'], opacity=0, alignmentgroup='week_start', 
                    customdata = df_bar["lineage"],
                    hovertemplate='Paraugu daļa: %{y:.2f}'+'<br>Nedēļa: %{x}'+'<br>Paveids: %{customdata: .3f}')

#defining a function that will combine the plots

def figures_to_html(figs, filename="dashboard_" + datestr + ".html"):
    dashboard = open(filename, 'w', encoding="utf-8")
    dashboard.write('<html><head></head><body>' + "\n")
    dashboard.write("<div><h1>SARS-CoV-2 mutāciju un paveidu apkopojums</h1></div>")
    for fig in figs:
        inner_html = fig.to_html().split('<body>')[1].split('</body>')[0]
        dashboard.write(inner_html)
    dashboard.write("<div><p>Ziņojums sagatavots: " + datestrread + "</p></div>")
    dashboard.write("</body></html>" + "\n")

#creating the joint report
figures_to_html([fig, fig4])

#export just one plot

#fig4.write_html(str("sample_bar_inv_" + timestr + ".html"))

#fig4.show()
