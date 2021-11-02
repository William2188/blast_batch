import requests
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import re
import lxml.html as lh
from bs4 import BeautifulSoup
import pandas as pd
def gender_genie(text):
    url = 'https://www.tcdb.org/progs/blast.php'
    caption = ''

    form_data = {
        'BLAST': 'blastp',
        'EXPECT': '1',
        'DESCRIPTIONS': '50',
        'ALIGNMENTS': '50',
        'MATRIX': 'BLOSUM62',
        'OGAP': '11',
        'EGAP': '1',
        'FILTER': 'no',
        'SEQUENCE': text,
        'submit': 'Click here to BLAST',
    }

    response = requests.post(url, data=form_data)

    #tree = lh.document_fromstring(response.content)
    soup = BeautifulSoup(response.text, 'html.parser')
    #return tree.xpath("")

    return soup.find_all('pre')[1]

def match_record(lines, rx_name):
    rx_name = rx_name.lower()[0:rx_name.find('transport')]
    rx_name = rx_name[0:rx_name.find('(')]
    rx_name = rx_name.replace('l-','')
    rx_name = rx_name.replace('l ', '')
    rx_name = rx_name.replace('reversible', '')
    rx_name = rx_name.replace('probable', '')
    rx_name = rx_name.replace('export', '')
    rx_name = rx_name.replace('via', '')
    rx_name = rx_name.replace('atpase', '')
    rx_name = rx_name.replace('atpa', '')
    rx_name = rx_name.replace('d ', '')
    rx_name = rx_name.lstrip()
    rx_name = rx_name.rstrip()
    

    state = 0
    output = []
    record = ''
    for line in lines:
        if line == '':
            pass
        elif ">gnl" in line and state == 0:
            state = 1
            record = line
        elif "Query" in line and state == 1:
            state = 0
            output.append(record)
        elif state == 1:
            record = record + ' ' + line

        else:
            pass
    Partial_Ratio = 0
    closest_str = ''
    expect_val = 1.0
    for item in output:
        expect=(re.findall('Expect = (\d+.\d+e-\d+)', item))
        expect_val1 = 1.0
        for x in expect:
            expect_val1 = float(x)
        name_str = item.split('|')[5].split('\t')[0]


        name_str = " ".join(name_str.split()[1:-1])
        name_str = name_str[0:name_str.find('Length')].strip()

        if expect_val1 < 1.0e-05:
            #print (expect_val1, expect_val1 < 1.0e-05)
            #print(item)
            #print (name_str + ', expect = ' + str(expect_val))
            #rx_name = '3 aminobutyrate transport in via proton symport'

            #rx_name = rx_name.split('transport')[0].strip()
            #Token_Sort_Ratio = fuzz.token_sort_ratio(rx_name, name_str)
            Partial_Ratio1 = fuzz.partial_ratio(rx_name.lower(), name_str.lower())
            #if Partial_Ratio1 > 80:
            #    print('string to search: ' + rx_name.lower() + '\t found cloest string: ' + name_str.lower() + '\t score: ' + str(Partial_Ratio1))
            if Partial_Ratio1 >= Partial_Ratio:
                Partial_Ratio = Partial_Ratio1
                closest_str = name_str.lower()
                expect_val = expect_val1
    #print('string to search: ' + rx_name.lower() + '\t found cloest string: ' + name_str.lower() + '\t score: ' + str(Partial_Ratio))
    return rx_name.lower(), closest_str.lower(), Partial_Ratio, expect_val

def process_datafile(model_fname):
    file_lines = []
    gr_list = []
    sequence_list = []
    with open(model_fname) as f:
        for line in f:
            file_lines.append(line)
    seq = ''
    for line in file_lines:
        if line[0] == '>':
            gr_list.append(line)

            if seq != '':
                sequence_list.append(seq)

                seq = ''
        else:
            seq = seq + line.strip('\n')
    if seq != '':
        sequence_list.append(seq)
    return gr_list, sequence_list


if __name__ == '__main__':
    #results = gender_genie('MEPSPSVTLEPQPAGPPSAESGLRRSMGPRHLVMIAMGGVIGSGLFLSSGYTISQAGPLG\
    #AVIAYLIGSFVVYLVMACLGELAIAYPVSGAFHIYAARSIGPATGFATAWLYWLCWAVAI\
    #GSEFTAAGLLMQRWFPGIDVWVWCVVFAAILFASNAVSSRFFGESEFWFSIVKVGTIIVL\
    #IVLGGAALVGFHPLAASGNHPFLFENFNTPGGLFPNGFTGVLVTALAVFYAFSGSELIGV\
    #AAGETADPATSIPKAMRTTVFRLLIFFVGAIAVIAATLPFDQVGVDESPFVSVLSSIGIP\
    #FAADIMNFVIITALLSAGNSGLYSCARMLYSLSEEGHAPRAFRRLTRRGIPMIALSVSML\
    #GGLASLISSVVAPETVYLVLVSVAGFAVVGVWMSITASHFFHRRSFVKNGGNVAALSYRA\
    #PLYPLVPILAFSLCFISLIGIAFDPNQVAALYFGIPFVGACYAFFYFKYGRGARTAELA')
    #lines = results.text.split('\n')

    df = pd.read_excel('Transport collab.xlsx')
    df['gr_rule'].fillna('', inplace=True)
    last_file_name = ''
    df_output = pd.DataFrame()
    for index, row in df.iterrows():
        row['gr_rule_corrected2'] = ''
        rx_name = str(row['rxn_name'])
        model = row['model']
        model_fname = 'data/'+str(model)+'.faa'
        print(model_fname)
        gr_list = []
        sequence_list = []
        if(model_fname != last_file_name):
            gr_list, sequence_list = process_datafile(model_fname)
        gr_rules_str = str(row['gr_rule']).replace('(','').replace(')','').replace('AND', 'or').replace('and','or')
        print("row" + " " + str(index) + "\n" + gr_rules_str)
        if(gr_rules_str != ''):
            gr_rules = gr_rules_str.split('or')
            for gr in gr_rules:
                gr=gr.strip()
                #print(gr)
                score = 0
                for idx in range(len(gr_list)):
                    if gr in gr_list[idx]:
                        results = gender_genie(sequence_list[idx])
                        lines = results.text.split('\n')
                        print(gr, rx_name, match_record(lines, rx_name))
                        temp1, temp2, score, expect_val = match_record(lines, rx_name)

                        if score >= 70:
                            if row['gr_rule_corrected2'] == '':
                                row['gr_rule_corrected2'] = gr
                            else:
                                row['gr_rule_corrected2'] = row['gr_rule_corrected2'] + ' or ' + gr
                        break
            print(row['gr_rule_corrected2'])

        pd2=(pd.DataFrame([row]))
        if df_output.empty:
            df_output=pd2
        else:
            df_output=df_output.append(pd2)

        df_output.to_excel('Transport collab updated.xlsx', index=False)












