import os
import numpy as np
import pandas as pd

configfile: "config.yaml"

def if_exist(path):
    p = path.split('/')
    if os.path.exists('~/scanb_data/'+p[2]+'/'+p[3]+'/'+p[4]+'/'+'alignment.bam'):
        exist.append(path)

def sample_path(wc):
    return '~/scanb_data/'+'/'.join(str(gene[gene['SAMPLE'] == wc.sp]['DataFilesFolder']).split('/')[2:-1])+'/'+'alignment.bam'

def which_gtf(wc):
    if wc.which == 'normal':
        return config["gene_coordinate"]
    elif wc.which == 'denovo':
        return config["denovo_coordinate"]


exist = []
gene = pd.read_csv(config["cohort_table"])[['SAMPLE', 'DataFilesFolder']]
gene['SAMPLE'] = gene['SAMPLE'].apply(lambda x: x.split('.')[0])
gene = gene[gene.SAMPLE.apply(lambda x: True if x in config["samples"] else False)]
gene.DataFilesFolder.apply(lambda x: if_exist(x))
gene = gene[gene.DataFilesFolder.apply(lambda x: True if x in exist else False)]



EVENT = 'Ki67_fake_event'
DAYS = 'Ki67_fake_days'

rule all:
    input:
        'result/'+EVENT[:-6].lower()+'.jpg',
        #'result/denovo-ki67_censored_ts.jpg',
        #'result/denovo-ki67_censored_gs.jpg',
        #'result/normal-ki67_censored_ts.jpg',
        #'result/normal-ki67_censored_gs.jpg',
        #'result/v27-rfi_age_trm_gs.jpg',
        #'denovo/transcripts/'+EVENT.split("_")[0]+'-best.feature',
        #'denovo/genes/'+EVENT.split("_")[0]+'-best.feature',
        #'normal/transcripts/'+EVENT.split("_")[0]+'-best.feature',
        #'normal/genes/'+EVENT.split("_")[0]+'-best.feature',




rule pre_assembly:
    input:
        bam = sample_path,
        gtf = config["gene_coordinate"]
    output:
        gtf = "pre_assembly/{sp}.gtf",
    threads: 2
    resources:
        runtime = "01:00:00"
    shell:
        "stringtie {input.bam} -l {wildcards.sp} -p 2 -G {input.gtf} -o {output.gtf}"

rule merge:
    input:
#        bam = expand("pre_assembly/{sp}.gtf", sp=config["samples"]),
        gtf = config["gene_coordinate"]
    output:
        gtf = "pre_assembly/merged.gtf"
    resources:
        runtime = "01:00:00"
    shell:
        "ls pre_assembly/*.gtf | head -10 > pre_assembly/mergelist.txt; "
        "stringtie --merge -G {input.gtf} -c 1 -o {output.gtf} pre_assembly/mergelist.txt"

rule assembly:
    input:
        bam = sample_path,
        gtf = which_gtf
    output:
        gtf = "{which}/{sp}/{sp}.gtf",
        csv = "{which}/{sp}/{sp}.tsv",
        cov = "{which}/{sp}/{sp}_cov.gtf"
    threads: 2
    resources:
        runtime = "02:00:00"
    shell:
        "stringtie {input.bam} -l {wildcards.sp} -p 2 -e -G {input.gtf} -o {output.gtf} -A {output.csv} -C {output.cov}"

rule extract:
    input:
        expand('{which}/{sp}/{sp}.tsv', sp=config["samples"], which=config["mode"])
    output:
        '{which}/{set}.csv'
    resources:
        runtime = "10:00:00"
    script:
        'scripts/extract.py'


def which_feature(wc):
    if wc.n == '1':
        return  config["sample_info"]
    else:
        return '{which}/{set}/'+EVENT.split("_")[0]+'-features.Round'+str(int(wc.n)-1)


rule pre_selection:
    input:
        data = '{which}/{set}.csv',
        info = config["sample_info"],
        feature = which_feature
    output:
        '{which}/{set}/{batch}.R{n,[0-9]+}'
    threads: 16
    resources:
        runtime = "2:00:00"
    params:
        event = EVENT,
        days = DAYS
    script:
        "scripts/pre_selection.py"

def how_many_batches(wc):
    if wc.n == '1':
        with open(wc.which+'/'+wc.set+'.csv') as f:
            i = 0
            for l in f:
                i += 1
        # splits = i//(113*1.6)+1
        splits = i//(336*1.6)+1
        batches = [str(b) for b in range(1, int(np.log(1-0.95)/np.log(1-(1/splits))))]
        return expand('{which}/{set}/{batch}.R{n}', which=wc.which, set=wc.set, batch=batches, n = wc.n)
    else:
        i = len(pd.read_csv(wc.which+'/'+wc.set+'/'+EVENT.split("_")[0]+'-features.Round'+str(int(wc.n)-1)).Feature)
        splits = i//(113*1.6)+1
        if splits == 1:
            splits = 4
        batches = [str(b) for b in range(1, int(np.log(1-0.95)/np.log(1-(1/splits))))]
        return expand('{which}/{set}/{batch}.R{n}', which=wc.which, set=wc.set, batch=batches, n = wc.n)



checkpoint multivariant_analysis:
    input:
        how_many_batches
    output:
        '{which}/{set}/'+EVENT.split("_")[0]+'-features.Round{n}',
        '{which}/{set}/'+EVENT.split("_")[0]+'-batches.Round{n}',
        '{which}/{set}/'+EVENT.split("_")[0]+'-describe.Round{n}'
    threads: 1
    resources:
        runtime = "1:00:00"
    script:
        "scripts/filter.py"

def check_result(wc):
    if 'n' not in wc:
        wc = {'which':wc.which, 'set':wc.set, 'n':1}
    with checkpoints.multivariant_analysis.get(**wc).output[2].open() as f:
        mean_ci = float(f.readline().split(',')[0][6:])
        if wc['n'] != 1:
            with open(os.path.join(wc['which'], wc['set'], EVENT.split("_")[0]+'-describe.Round'+str(wc['n']-1))) as p:
                mean_ci_pre = float(p.readline().split(',')[0][6:])
        else:
            mean_ci_pre = 0
        # for test
        if mean_ci > 0 and wc['n'] < 9:
            wc['n'] += 1
            return check_result(wc)
        # for test

        #if mean_ci > mean_ci_pre and wc['n'] < 20:
        #    wc['n'] += 1
        #    return check_result(wc)
        else:
            return os.path.join(wc['which'], wc['set'], EVENT.split("_")[0]+'-features.Round'+str(wc['n']-1))




rule how_many_rounds:
    input:
        check_result
    output:
        '{which}/{set}/'+EVENT.split("_")[0]+'-best.feature',
        '{which}/{set}/'+EVENT.split("_")[0]+'-conclusion.csv'
    script:
         'scripts/conclude1.py'

rule conclusion:
    input:
        'denovo/transcripts/'+EVENT.split("_")[0]+'-conclusion.csv',
        'denovo/genes/'+EVENT.split("_")[0]+'-conclusion.csv',
        'normal/transcripts/'+EVENT.split("_")[0]+'-conclusion.csv',
        'normal/genes/'+EVENT.split("_")[0]+'-conclusion.csv',
    output:
        'result/'+EVENT.split("_")[0]+'_conculsion.csv'
    script:
        'scripts/conclude2.py'



rule onerun_selection:
    input:
        data = '{which}/{set}.csv',
        info = config["sample_info"]
    output:
        '{which}/{set}/'+EVENT.split("_")[0]+'_1run_{which}_{set}_grids_result.csv',
        '{which}/{set}/'+EVENT.split("_")[0]+'_1run_{which}_{set}_best_coefs.csv'
    threads: 16
    resources:
        runtime = "48:00:00"
    params:
        event = EVENT,
        days = DAYS
    script:
        "scripts/one_run.py"




rule cox_uni:
    input:
        data = '{which}/{set}.csv',
        info = config["sample_info"]
    output:
        '{which}/{set}/'+EVENT.split("_")[0]+'_uni.csv'
    threads: 16
    resources:
        runtime = "48:00:00"
    params:
        event = EVENT,
        days = DAYS
    script:
        "scripts/cox_unitest.py"




rule cox_fdr:
    input:
        '{which}/{set}/'+EVENT.split("_")[0]+'_uni.csv'
    output:
        '{which}/{set}/'+EVENT.split("_")[0]+'_uni_fdr.csv'
    threads: 1
    resources:
        runtime = "01:00:00"
    script:
        "scripts/FDR.R"


rule rename:
    input:
        '{which}/{set}/'+EVENT.split("_")[0]+'_uni.csv',
        '{which}/{set}/'+EVENT.split("_")[0]+'_1run_{which}_{set}_best_coefs.csv',
        '{which}/{set}/'+EVENT.split("_")[0]+'-best.feature'
    output:
        ''
    shell:
        ''

rule cox_analysis:
    input:
        data = lambda wc: '{which}/genes.csv' if 'gs' in wc.path else '{which}/transcripts.csv',
        info = config["sample_info"],
    output:
        '{which}/{path}/cox_'+EVENT.split("_")[0]+'_{method}_analysis.txt'
    threads: 16
    resources:
        runtime = "01:00:00"
    params:
        event = EVENT,
        days = DAYS
    script:
        'scripts/cox_analysis.py'


rule cox_plot:
    input:
        '{which}/{path}/cox_'+EVENT.split("_")[0]+'_uni_analysis.txt',
        '{which}/{path}/cox_'+EVENT.split("_")[0]+'_1run_analysis.txt',
        '{which}/{path}/cox_'+EVENT.split("_")[0]+'_ens_analysis.txt',
        config["sample_info"]
    output:
        'result/{which}-{path}.jpg'
    params:
        event = EVENT,
        days = DAYS
    script:
        'scripts/cox_plot.py'

rule cox_plot2:
    input:
        'result/denovo-'+EVENT.split("_")[0].lower()+'_age_trm_ts.jpg',
        'result/denovo-'+EVENT.split("_")[0].lower()+'_age_trm_gs.jpg',
        'result/normal-'+EVENT.split("_")[0].lower()+'_age_trm_ts.jpg',
        'result/normal-'+EVENT.split("_")[0].lower()+'_age_trm_gs.jpg',
        'result/v27-'+EVENT.split("_")[0].lower()+'_age_trm_gs.jpg',
        config["sample_info"]
    output:
        'result/'+EVENT[:-6]+'.jpg'
    params:
        event = EVENT,
        days = DAYS
    script:
        'scripts/cox_plot2.py'



rule cox_plot3:
    input:
        'result/denovo-'+EVENT.split("_")[0].lower()+'_censored_ts.jpg',
        'result/denovo-'+EVENT.split("_")[0].lower()+'_censored_gs.jpg',
        'result/normal-'+EVENT.split("_")[0].lower()+'_censored_ts.jpg',
        'result/normal-'+EVENT.split("_")[0].lower()+'_censored_gs.jpg',
        config["sample_info"]
    output:
        'result/'+EVENT[:-6].lower()+'.jpg'
    params:
        event = EVENT,
        days = DAYS
    script:
        'scripts/cox_plot3.py'
