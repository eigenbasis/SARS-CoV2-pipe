# This script generates an .html report about sequencing results

import csv
import datetime as dt

# 1. Set up multiple variables to store the titles, text within the report
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
reprot_date = dt.datetime.now().strftime("%d/%m/%Y")

#read in the style sheet
with open('style.txt', 'r') as file:
    style = file.read().replace('\n', '')

#read in the base64 encoded logo picture
with open('image.txt', 'r') as file:
    logo = file.read().replace('\n', '')

#generate the reports
for i in content:
    print("Reporting sample: ", i["sample_nr"])

    # 2. Combine them together using a long f-string
    html = f'''
        <html>
            <head>
                <title>{i["sample_nr"]}</title>
                <style>{style}</style>
            </head>
            <body>
                <div class="flex-container">
                    <div>
                    <img src="{logo}" >
                    </div>
                    <div class="signature">
                        <p><strong>Laboratorijas dienests <br>
                            SIA Rīgas Austrumu klīniskā universitātes slimnīca</strong> <br>
                            Linezera iela 3, k-4, Rīga, LV-1006, Latvija <br>
                            Paraugu pieņemšanas daļa (tālr.: 67014573; 29114493) <br>
                            Sekvencēšanas daļa (tālr.: 67014523)</p>
                    </div>
                </div>
                <h1>SARS-CoV-2 RNS pilna genoma sekvencēšana</h1>
                <h2>Parauga Nr.  {i["sample_nr"]}</h2>
                <div class="stats">
                    <p>Atskaites izveidošanas datums: {reprot_date}<br>
                    Parauga ņemšanas datums: {i["sample_date"]}<br>
                    Sekvencēšanas datums: {i["seq_date"]}<br>
                    Testējamais materiāls: <strong>{i["test_material"]}</strong><br>
                    Sekvencēšanas institūcija: <strong>{i["seq_inst"]}</strong></p>
                </div>
                <div class="conclusion">
                    <h2>Slēdziens</h2>
                    <p>Testētajā materiālā noteiktā SARS-CoV-2 vīrusa RNS atbilst <strong>{i["lineage"]}</strong> (Omikron, 21K) līnijai 
                    <a href="https://covariants.org/variants/21K.Omicron" target="_blank">https://covariants.org/variants/21K.Omicron</a>
                    </p>
                </div>
                <div class="results">
                    <h2>Rezultāti</h2>
                    <p>Paveids: <strong>{i["lineage"]}</strong><br>
                    Genoma garums: {i["gen_len"]}<br>
                    Genoma N %: {i["gen_N_perc"]}<br>
                    Genoma GC saturs: {i["gen_GC"]}<br>
                    Vidējais pārklājums: <strong>{i["avg_cov"]}</strong></p>
                </div>
                <div class="method">
                    <h2>Metodes apraksts un sekvencēšanas kvalitātes rādītāji</h2>
                    <p>Pilna genoma sekvence iegūta izmantojot CovidSeq bibliotēku sagatavošanas komplektu un sekvencēšana
                        izmantojot MidOutput 300 ciklu reaģentu kārtridžu un standarta plūsmas šūnu uz NextSeq 550 Dx platformas.
                        Sekvenču dati tika apstrādāti, izmantojot SIA RAKUS Laboratorijas dienesta izveidotu darbplūsmu (<a href="https://github.com/NMRL/SARS-CoV2_assembly" target="_blank">https://github.com/NMRL/SARS-CoV2_assembly</a>). Iegūtie sekvenču dati tika filtrēti izmantojot cutadapt 2.31 un fastp 0.20.1, lai atbrīvotos no zemas kvalitātes nolasījumiem. Filtrētās sekvences tika pielīdzinātas Wuhan-Hu-1 SARS-CoV-2 izolāta references genomam (Accession number: MN908947.3) ar bwa 0.7.17-r1198-dirty mem un samtools 1.12 [4]. ARTIC v3 amplifikācijas oligonukleotīdi tika noņemti, izmantojot iVar
                        1.3.1 un bedtools v.2.30.00. Pēc kvalitātes kontroles nolasījumi tika atkārtoti pielīdzināti Wuhan-Hu-1 references genomam, izmantojot abra 0.97. Variantu identificēšana un consensus sekvences izveidošana tika veikta ar freebayes v0.9.21, vcflib 1.0.2/vcffilter un
                        Consensusfixer 0.4. Variantu anotācijai tika izmantots snpEff 5.0e un PANGO celmu piešķiršana tika veikta, izmantojot Pangolin.
                    </p>
                </div>
                <br>
                <div class="persons">
                    <p>Sekvencēšanu veica: Jūlija Čevere, Kristīne Liepiņa, Anastasija Žuravļova<br>
                    Bioinformātiķi: Reinis Vangravs, Jevgēnijs Bodrenko Atbildīgā persona: Ģirts Šķenders
                    </p>
                </div>
                <br>
                <div class="footer">
                    <p>Sekvencēšanas pārskats ir sastādīts elektroniski un tam ir informatīvs raksturs.</p>
                </div>
            </body>
        </html>
        '''
    # 3. write out the generated html string to file and encode it in order to maintain special characters 
    filename = ("filled_temp_" + i["sample_nr"] + ".html")
    with open(filename, 'w', encoding="utf-16") as f:
        f.write(html)
