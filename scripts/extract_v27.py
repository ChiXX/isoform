import pandas as pd
import numpy as np

data = pd.read_csv('sp_list.csv')
samples = data.SAMPLE.to_list()
table = ['v27/transcripts_extract.csv']

def get_info(l, *info):
    try:
        info_dic =  {i.strip().split(' ')[0]:i.strip().split(' ')[1].strip('"') for i in l.split('\t')[-1][:-3].split(';')}
    except:
        info_dic = {i.strip().split(' ')[0]:i.strip().split(' ')[1].strip('"') for i in l.split('\t')[-1][:-2].split(';')}
    return ','.join([info_dic[i] for i in info])

def extract_t(sp):
    wdir = data[data.SAMPLE == sp].DataFilesFolder.to_list()[0]
    tmp = pd.read_csv('v27'+wdir+'/transcript.gtf', sep='\t', comment='#', header=None)
    tmp = tmp[tmp[2] == 'transcript'][[8]]
    tmp[['gene_id', 'transcript_id', sp]] = pd.DataFrame(tmp[8].apply(lambda x: get_info(x, 'gene_id', 'transcript_id', 'FPKM')))[8].str.split(',', expand=True)
    tmp['ID'] = tmp.apply(lambda x: x['gene_id']+'_'+x['transcript_id'], axis=1)
    tmp[sp] = tmp[sp].astype(float)
    tmp.drop(['gene_id', 'transcript_id', 8], axis=1, inplace=True)
    tmp.drop_duplicates(subset='ID', keep='first', inplace=True)
    return tmp

def extract_g(sp):
    wdir = data[data.SAMPLE == sp].DataFilesFolder.to_list()[0] 
    tmp = pd.read_csv('v27'+wdir+'/gene.tsv', sep='\t')[['Gene ID', 'Gene Name', 'FPKM']]
    tmp.columns = ['ID', 'GENE', sp]
    tmp['ID'] = tmp.apply(lambda x: x['GENE']+'_'+x['ID'], axis=1)
    tmp.drop('GENE', axis=1, inplace=True)
    tmp.drop_duplicates(subset='ID', keep='first', inplace=True)
    return tmp


def extract_transcript(t):
    data = extract_t(samples[0])
    for sp in samples[1:]:
        try:
            data = pd.merge(data, extract_t(sp), on='ID', how='left')
        except:
            data[sp] = np.nan
    data.to_csv(t, index=False)

        
def extract_gene(t):
    data = extract_g(samples[0])
    for sp in samples[1:]:
        try:
            data = pd.merge(data, extract_g(sp), on='ID', how='left')
        except:
            data[sp] = np.nan
    data.to_csv(t, index=False)

for t in table:
    if 'transcripts' in t:
        extract_transcript(t)
    if 'genes' in t:
        extract_gene(t)
        
       
