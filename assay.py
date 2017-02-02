import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import urllib.request
import zipfile
import argparse
import re
import os
import sys

COLUMNS = ['PUBCHEM_AID', 'PUBCHEM_RESULT_TAG', 'PUBCHEM_SID', 'PUBCHEM_CID',\
        'PUBCHEM_ACTIVITY_OUTCOME', 'PUBCHEM_ACTIVITY_SCORE']
PANEL_COLUMNS = ['PUBCHEM_AID', 'RESULT_PANEL_ID', 'PUBCHEM_RESULT_TAG', 'PUBCHEM_SID', 'PUBCHEM_CID',\
        'PUBCHEM_ACTIVITY_OUTCOME', 'PUBCHEM_ACTIVITY_SCORE']
EXTRAS = ['PUBCHEM_ACTIVITY_URL', 'PUBCHEM_ASSAYDATA_COMMENT']
OUTCOME = {'Active':1, 'Inactive':0, 'Unspecified':3, 'Inconclusive':4, 'Probe':5}

filename2aid = lambda x: int(os.path.split(x)[-1].replace('.csv.gz', ''))

def download(outdir='data/assay'):
    annotation = 'ftp://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/Aid2Annotation.gz'
    annotation = pd.read_csv(annotation, sep='\t', compression='gzip')
    annotation = pd.pivot_table(annotation, index='AID', columns='Title', values='Annotation', aggfunc=lambda x: ', '.join(x))
    annotation = annotation[annotation.notnull().sum().sort_values(ascending=False).index]

    depositor = 'http://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/Aid2DepositorName.gz'
    depositor = pd.read_csv(depositor, sep='\t', compression='gzip')
    annotation = pd.merge(depositor, annotation.reset_index(), on='AID')

    annotation = annotation.sort_values('AID')
    annotation.to_pickle('data/pkl/annotation.pkl')

    base = 'http://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Data/%s'

    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    for line in urllib.request.urlopen(base % '').readlines():
        filename = re.search('\d+_\d+.zip', str(line))
        if not filename:
            continue

        filename = filename.group()

        if os.path.exists(os.path.join(outdir, filename)):
            continue

        print(base % filename)
        try:
            urllib.request.urlretrieve(base % filename, os.path.join(outdir, filename))
        except Exception as e:
            print(e)

def panel(aid, data):
    data = data.set_index('PUBCHEM_RESULT_TAG')
    outcomes = data.ix['RESULT_PANEL_ID', :].str.endswith('_OUTCOME')
    outcomes = outcomes[outcomes.values == True].index
    scores = data.ix['RESULT_PANEL_ID', :].str.endswith('_AC')
    scores = scores[scores.values == True].index

    for outcome, score in zip(outcomes, scores):
        df = data.reset_index()[PANEL_COLUMNS[:-2] + [outcome,score]]
        panel_id = data.ix['RESULT_PANEL_ID', outcome].replace('_OUTCOME', '')
        df['RESULT_PANEL_ID'] = panel_id
        df.columns = PANEL_COLUMNS
        yield panel_id, df

def build(aids, name, indir='data/assay', outdir='data/pkl', other=False):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    result = []

    for filename in sorted(os.listdir(indir)):
        archive = zipfile.ZipFile(os.path.join(indir, filename))

        for csv in sorted(archive.namelist(), key=filename2aid):
            aid = filename2aid(csv)
            if not other and aid not in aids: continue
            if other and aid in aids: continue

            data = archive.open(csv)
            data = pd.read_csv(data, sep=',', compression='gzip', low_memory=False)
            data['PUBCHEM_AID'] = aid
            data['RESULT_PANEL_ID'] = None

            if data['PUBCHEM_RESULT_TAG'].str.contains('PANEL').sum() > 0:
                for panel_id, df in panel(aid, data):
                    cond = df['PUBCHEM_RESULT_TAG'].astype(str).str.startswith('RESULT')
                    df = df[-cond].copy()
                    df['PUBCHEM_RESULT_TAG'] = df['PUBCHEM_RESULT_TAG'].astype(int)
                    df['PUBCHEM_ACTIVITY_OUTCOME'] = df['PUBCHEM_ACTIVITY_OUTCOME'].map(lambda x: OUTCOME.get(x, x))

                    result.append(df)
                    print(name, aid, panel_id, df.shape)
            else:
                cond = data['PUBCHEM_RESULT_TAG'].astype(str).str.startswith('RESULT')
                data = data.ix[-cond, COLUMNS]
                data['PUBCHEM_RESULT_TAG'] = data['PUBCHEM_RESULT_TAG'].astype(int)
                data['PUBCHEM_ACTIVITY_OUTCOME'] = data['PUBCHEM_ACTIVITY_OUTCOME'].map(lambda x: OUTCOME.get(x, x))

                result.append(data)
                print(name, aid, data.shape)

    if len(result) > 0:
        result = pd.concat(result)
        result.to_pickle(os.path.join(outdir, '%s.pkl' % name.replace('/','-')))

def build_others(indir='data/assay', outdir='data/pkl'):
    annotation = pd.read_pickle(os.path.join(outdir, 'annotation.pkl'))
    aids = annotation['AID']

    filename = 'others.pkl'
    if not os.path.exists(os.path.join(outdir, filename)):
        build(aids, 'others', indir, outdir, other=True)

    if os.path.exists(os.path.join(outdir, filename)):
        data = pd.read_pickle(os.path.join(outdir, filename))        
        print('others', data.shape)

def build_core(indir, outdir, column):
    outdir = os.path.join(outdir, column)
    annotation = pd.read_pickle('data/pkl/annotation.pkl')
    df = annotation.groupby(column).size().sort_values(ascending=False)
    df = df.to_frame().rename(columns={0:'count'}).reset_index()

    for i,row in df.iterrows():
        cond = annotation[column] == row[column]
        aids = annotation.ix[cond, 'AID']
        filename = '%s.pkl' % row[column].replace('/','-')
        if not os.path.exists(os.path.join(outdir, filename)):
            build(aids, row[column], indir, outdir)

        if os.path.exists(os.path.join(outdir, filename)):
            data = pd.read_pickle(os.path.join(outdir, filename))        
            print(row[column], data.shape)

def depositor(indir='data/assay', outdir='data/pkl'):
    build_core(indir, outdir, 'DepositorName')

def assay_type(indir='data/assay', outdir='data/pkl'):
    build_core(indir, outdir, 'Assay Type')

def assay_format(indir='data/assay', outdir='data/pkl'):
    build_core(indir, outdir, 'Assay Format')

def assay_method(indir='data/assay', outdir='data/pkl'):
    build_core(indir, outdir, 'Assay Detection Method')

def cell_type(indir='data/assay', outdir='data/pkl'):
    build_core(indir, outdir, 'Assay Cell Type')

def summary():
    if not os.path.exists('figure'):
        os.makedirs('figure')

    annotation = pd.read_pickle('data/pkl/annotation.pkl')

    indir = 'data/pkl'
    writer = pd.ExcelWriter(os.path.join('data', 'summary.xlsx'))
    for category in ['DepositorName','Assay Type','Assay Format','Assay Detection Method','Assay Cell Type']:
        mat = []
        for basename in annotation[category].dropna().unique():

            filename = os.path.join(indir, category, '%s.pkl' % basename.replace('/','-'))
            if os.path.exists(filename):
                df = pd.read_pickle(filename)
                assays = df['PUBCHEM_AID'].nunique()
                outcomes = df.groupby('PUBCHEM_ACTIVITY_OUTCOME').size().to_dict()
                outcome_labels = sorted(OUTCOME.items(), key=lambda x: x[0])
                outcomes = [outcomes.get(v,0) for k,v in outcome_labels]
                outcome_labels = [x[0] for x in outcome_labels]

                if 'RESULT_PANEL_ID' in df.columns:
                    cond = df['RESULT_PANEL_ID'].notnull()
                    panels = df.ix[cond, 'PUBCHEM_AID'].nunique()
                else:
                    panels = 0

                # abbreviation
                basename = basename.replace(',', ',\n')

                mat.append([category, basename, assays, panels, df.shape[0]] + outcomes)
        df = pd.DataFrame(mat)
        df.columns = ['Category',category,'Assays','Panels','Data Points'] + outcome_labels
        df = df.sort_values('Data Points')
        df.to_excel(writer, sheet_name=category, index=False)

        top = min(20, df.shape[0])
        df.iloc[-top:].plot(kind='barh', x=category, y='Data Points', figsize=(10,8), 
                fontsize=10, title='Data points by %s (Top %d)' % (category, top))
        plt.legend(loc='lower right', fontsize=10)
        plt.tight_layout()
        plt.savefig('figure/%s.png' % category)

    writer.close()

def pivot(indir='data/pkl', outdir='data/pkl'):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    df = pd.read_pickle(os.path.join(indir, 'result.pkl'))
    df = pd.pivot_table(df, index='PUBCHEM_AID', columns='PUBCHEM_ACTIVITY_OUTCOME', values='PUBCHEM_RESULT_TAG', aggfunc='count')
    outcome = {v:k for k,v in OUTCOME.items()}
    df.columns = df.columns.map(lambda x: outcome[x])

    df = df.fillna(0).astype(int)
    df.to_pickle(os.path.join(outdir, 'pivot.pkl'))

def subset():
    df = pd.read_pickle('data/pkl/result.pkl')
    df = df[['PUBCHEM_AID', 'PUBCHEM_CID', 'PUBCHEM_SID', 'PUBCHEM_ACTIVITY_OUTCOME']]
    print('Data potins', df.shape, 'Assays', df['PUBCHEM_AID'].nunique(), 'Unique compounds', df['PUBCHEM_CID'].nunique(), 'Unique substances', df['PUBCHEM_SID'].nunique())

    df = df[df['PUBCHEM_CID'].notnull()]
    df['PUBCHEM_CID'] = df['PUBCHEM_CID'].astype(int)
    print('Compound ID is not null')
    print('Data potins', df.shape, 'Assays', df['PUBCHEM_AID'].nunique(), 'Unique compounds', df['PUBCHEM_CID'].nunique())

    cond = np.logical_or(df['PUBCHEM_CID'].notnull(), df['PUBCHEM_SID'].notnull())
    df = df[cond]
    print('Compound ID is not null OR Substance ID is not null')
    print('Data potins', df.shape, 'Assays', df['PUBCHEM_AID'].nunique(), 'Unique compounds', df['PUBCHEM_CID'].nunique(), 'Unique substances', df['PUBCHEM_SID'].nunique())

    df = df[df['PUBCHEM_ACTIVITY_OUTCOME'].isin([OUTCOME['Active'],OUTCOME['Inactive']])]
    print('OUTCOME is Active or Inactive')
    print('Data potins', df.shape, 'Assays', df['PUBCHEM_AID'].nunique(), 'Unique compounds', df['PUBCHEM_CID'].nunique(), 'Unique substances', df['PUBCHEM_SID'].nunique())

    df.to_pickle('data/pkl/subset.pkl')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--download', help='download PCBA zip', action='store_true')
    parser.add_argument('--build', help='build PCAB descriptions and results pickle', action='store_true')
    parser.add_argument('--pivot', help='write pivot table to excel file', action='store_true')
    parser.add_argument('--subset', help='build subset for analysis', action='store_true')
    parser.add_argument('--summary', help='show data summary', action='store_true')

    parser.add_argument('--depositor', help='build depositor dataset', action='store_true')
    parser.add_argument('--assay_type', help='build assay type dataset', action='store_true')
    parser.add_argument('--assay_format', help='build assay format dataset', action='store_true')
    parser.add_argument('--assay_method', help='build assay method dataset', action='store_true')
    parser.add_argument('--cell_type', help='build assay cell type dataset', action='store_true')
    parser.add_argument('--others', help='build other assay dataset', action='store_true')
    args = parser.parse_args()

    if args.download:
        download()
    elif args.build:
        build() 
    elif args.pivot:
        pivot()
    elif args.subset:
        subset()
    elif args.summary:
        summary()
    elif args.depositor:
        depositor()
    elif args.assay_type:
        assay_type()
    elif args.assay_format:
        assay_format()
    elif args.assay_method:
        assay_method()
    elif args.cell_type:
        cell_type()
    elif args.others:
        build_others()
    else:
        parser.print_help()
