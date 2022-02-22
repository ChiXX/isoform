import pandas as pd

samples = list(set([i.split('/')[-1].split('.')[0] for i in snakemake.input]))
table = snakemake.output

def get_info(l, *info):
    try:
        info_dic =  {i.strip().split(' ')[0]:i.strip().split(' ')[1].strip('"') for i in l.split('\t')[-1][:-3].split(';')}
    except:
        info_dic = {i.strip().split(' ')[0]:i.strip().split(' ')[1].strip('"') for i in l.split('\t')[-1][:-2].split(';')}
    return ','.join([info_dic[i] for i in info])

def extract_t(sp, wdir):
    tmp = pd.read_csv(wdir+'/'+sp+'/'+sp+'.gtf', sep='\t', comment='#', header=None)
    tmp = tmp[tmp[2] == 'transcript'][[8]]
    tmp[['gene_id', 'transcript_id', sp]] = pd.DataFrame(tmp[8].apply(lambda x: get_info(x, 'gene_id', 'transcript_id', 'FPKM')))[8].str.split(',', expand=True)
    tmp['ID'] = tmp.apply(lambda x: x['gene_id']+'_'+x['transcript_id'], axis=1)
    tmp[sp] = tmp[sp].astype(float)
    tmp.drop(['gene_id', 'transcript_id', 8], axis=1, inplace=True)
    tmp.drop_duplicates(subset='ID', keep='first', inplace=True)
    return tmp

def extract_g(sp, wdir):
    tmp = pd.read_csv(wdir+'/'+sp+'/'+sp+'.tsv', sep='\t')[['Gene ID', 'Gene Name', 'FPKM']]
    tmp.columns = ['ID', 'GENE', sp]
    tmp['ID'] = tmp.apply(lambda x: x['GENE']+'_'+x['ID'], axis=1)
    tmp.drop('GENE', axis=1, inplace=True)
    tmp.drop_duplicates(subset='ID', keep='first', inplace=True)
    return tmp


def extract_transcript(t):
    wdir = t.split('/')[0]
    data = extract_t(samples[0], wdir)
    for sp in samples[1:]:
        data = pd.merge(data, extract_t(sp, wdir), on='ID', how='left')
    data.to_csv(t, index=False)

        
def extract_gene(t):
    wdir = t.split('/')[0]
    data = extract_g(samples[0], wdir)
    for sp in samples[1:]:
        data = pd.merge(data, extract_g(sp, wdir), on='ID', how='left')
    data.to_csv(t, index=False)

for t in table:
    if 'transcripts' in t:
        extract_transcript(t)
    if 'genes' in t:
        extract_gene(t)
        
       
