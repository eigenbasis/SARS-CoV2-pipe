from docxtpl import DocxTemplate, InlineImage
import datetime as dt
import csv

# create a document object from the template
doc = DocxTemplate("test_tempplate_C_seq.docx")

# read in the data to be filled
def csv_to_dic(filepath):
  content = []
  with open(filepath) as csvfile:
    csv_reader = csv.reader(csvfile)
    headers = next(csv_reader)

    for row in csv_reader:
      row_data = {key: value for key, value in zip(headers, row)}
      content.append(row_data)

  return content

content = csv_to_dic("filler_for_testing.csv")

for i in content:
    print("Reporting sample: ", i["sample_nr"])
    context = {
    "report_date": dt.datetime.now().strftime("%d-%b-%Y"),
    "sample_nr": i["sample_nr"],
    "sample_date": i["sample_date"],
    "test_material": i["test_material"],
    "seq_date": i["seq_date"],
    "seq_inst": i["seq_inst"],
    "lineage": i["lineage"],
    "gen_len": i["gen_len"],
    "gen_N_perc": i["gen_N_perc"],
    "gen_GC": i["gen_GC"],
    "avg_cov": i["avg_cov"],
    }

    doc.render(context)
    filename = ("filled_temp_" + i["sample_nr"] + ".docx")
    doc.save(filename)
