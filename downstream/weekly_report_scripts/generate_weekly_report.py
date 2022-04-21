from docxtpl import DocxTemplate, InlineImage
from docx.shared import Cm
import time, pandas as pd, sys, re

#Import template document + paths to files with data_1 & plots
template = DocxTemplate('/home/jevgenijs/Desktop/cov_flow/downstream/weekly_report_template.docx')
found_path = sys.argv[1] #path to the process_samples.py found_samples report (to get sampling date & lineages for each sample)
image_path= "/".join(found_path.split('/')[:-1])

#Generate data_1 lineage data_1 table
data_2 = pd.read_csv(found_path)
data_2.fillna("0",inplace=True)

voc_list = []
data_3 = data_2.groupby(data_2['lineage']).count()
exclude_dict = {}
for lin in data_3.index:
    if "AY." in lin and lin not in voc_list:
        exclude_dict[lin] = "B.1.617.2"
    elif "BA.1" in lin and lin not in voc_list:
        exclude_dict[lin] = 'BA.1'
    elif "BA.2" in lin and lin not in voc_list:
        exclude_dict[lin] = 'BA.2'


#Combine data_1 to report
sample_count= len(data_2)
by_lab_count = data_2.groupby(data_2['seq_institution']).count()
seq_lab_list = [f'{by_lab_count.loc[id]["receiving_lab_sample_id"]} no Eurofins' if 'Eurofins' in id else f'{by_lab_count.loc[id]["receiving_lab_sample_id"]} no NMRL' for id in list(by_lab_count.index)]
seq_lab_string = ', '.join(seq_lab_list)
lineage_data_table=[{'lineage':exclude_dict[lin] if lin in exclude_dict.keys() else lin,'sublin':lin if lin in exclude_dict.keys() else 'n/a','count':count} for lin, count in zip(data_3.index, data_3.sampling_date)]
total_data_table = [
    {'sample_number':number+1,
    'sample_id':sample,
    'sampling_date':date if date != '0' else 'n/a',
    'testing_lab':lab,
    'lineage':lin, 
    'sequenced':sequenced,
    'district':district if district != '0' else 'n/a'} for number,sample,date,lab,lin,sequenced,district in zip(data_2.index,data_2.receiving_lab_sample_id, data_2.sampling_date, data_2.testing_lab, data_2.lineage, data_2.seq_institution, data_2.district)]
total_df = pd.DataFrame(total_data_table)
data_2 = data_2[(data_2.sampling_date != '0') & (data_2.sampling_date != '44284')]
seq_date_range= f"no {min(data_2.sampling_date)} līdz {max(data_2.sampling_date)}"
batch_quality_1= InlineImage(template,f'{image_path}/plot1.png', width = Cm(16.35), height=Cm(5.7))
batch_quality_2=InlineImage(template,f'{image_path}/plot2.png', width = Cm(16.35), height=Cm(5.7))
batch_quality_3= InlineImage(template,f'{image_path}/plot3.png', width = Cm(16.35), height=Cm(5.7))
batch_quality_4=InlineImage(template,f'{image_path}/plot4.png', width = Cm(16.35), height=Cm(5.7))
batch_lineages= InlineImage(template,f'{image_path}/plot5.png', width = Cm(16.48), height=Cm(13.58))
batch_S_mut= InlineImage(template,f'{image_path}/plot6.png', width = Cm(17.19), height=Cm(16.49))
seq_sum_stats= InlineImage(template,f'{image_path}/plot7.png', width = Cm(17.19), height=Cm(12))
district_plot=InlineImage(template,f'{image_path}/plot8.png', width = Cm(17.19), height=Cm(17.31))


#Declare template variables
context = {
    'seq_date_range': seq_date_range,
    'seq_lab_string': seq_lab_string,
    'sample_count': sample_count,
    'total_data_table': total_data_table,
    'lineage_data_table': lineage_data_table,
    'plot_1': batch_quality_1,
    'plot_2': batch_quality_2,
    'plot_3': batch_quality_3,
    'plot_4': batch_quality_4,
    'plot_5': batch_lineages,
    'plot_6': batch_S_mut,
    'plot_7': seq_sum_stats,
    'plot_8': district_plot
    }

#Render automated report
template.render(context)
template.save(f'{image_path}/SARS-CoV2_sekv_pārskats_{time.strftime("%Y%m%d")}.docx')
total_df.to_csv(f'{image_path}/SARS-CoV2_sekv_pārskats_dati_{time.strftime("%Y%m%d")}.csv', header=True, index=False)
