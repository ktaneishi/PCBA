import pandas as pd
import numpy as np
import urllib.request
import zipfile
import argparse
import gzip
import re
import os
import sys

columns = ['PUBCHEM_AID', 'PUBCHEM_RESULT_TAG', 'PUBCHEM_SID', 'PUBCHEM_CID',\
        'PUBCHEM_ACTIVITY_OUTCOME', 'PUBCHEM_ACTIVITY_SCORE']
source = ['PUBCHEM_ACTIVITY_URL']
comment = ['PUBCHEM_ASSAYDATA_COMMENT']
OUTCOME = {'Active':1, 'Inactive':0, 'Unspecified':3, 'Inconclusive':4, 'Probe':5}

def download(outdir='data/assay'):
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

def build(indir='data/assay', outdir='pkl'):
    if os.path.exists(os.path.join(outdir, 'description.pkl')) and \
            os.path.exists(os.path.join(outdir, 'result.pkl')):
        return

    result, desc = [], []

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for filename in sorted(os.listdir(indir)):
        archive = zipfile.ZipFile(os.path.join(indir, filename))

        for csv in sorted(archive.namelist(), key=lambda x: int(os.path.split(x)[-1].replace('.csv.gz', ''))):
            aid = int(os.path.split(csv)[-1].replace('.csv.gz', ''))
            data = archive.open(csv)
            data = gzip.open(data)
            data = pd.read_csv(data, sep=',', low_memory=True)

            data['PUBCHEM_AID'] = aid

            # description
            cond = data['PUBCHEM_RESULT_TAG'].astype(str).str.startswith('RESULT')
            desc.append(data.ix[cond, columns])

            # result
            data = data.ix[-cond, columns]
            data['PUBCHEM_RESULT_TAG'] = data['PUBCHEM_RESULT_TAG'].astype(int)
            data['PUBCHEM_ACTIVITY_OUTCOME'] = data['PUBCHEM_ACTIVITY_OUTCOME'].map(lambda x: OUTCOME.get(x, x))
            result.append(data[columns])

            print(aid, data.shape)

    desc = pd.concat(desc)
    desc.to_pickle(os.path.join(outdir, 'description.pkl'))
    del desc
    result = pd.concat(result)
    result.to_pickle(os.path.join(outdir, 'result.pkl'))

def pivot(indir='pkl', outdir='pkl'):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    df = pd.read_pickle(os.path.join(indir, 'result.pkl'))
    df = pd.pivot_table(df, index='PUBCHEM_AID', columns='PUBCHEM_ACTIVITY_OUTCOME', values='PUBCHEM_RESULT_TAG', aggfunc='count')
    outcome = {v:k for k,v in OUTCOME.items()}
    df.columns = df.columns.map(lambda x: outcome[x])

    df = df.fillna(0).astype(int)
    df.to_pickle(os.path.join(outdir, 'pivot.pkl'))

def subset(threshold=100):
    df = pd.read_pickle('pkl/result.pkl')
    df = df[['PUBCHEM_AID', 'PUBCHEM_CID', 'PUBCHEM_SID', 'PUBCHEM_ACTIVITY_OUTCOME']]
    print('Data potins', df.shape, 'Assays', df['PUBCHEM_AID'].nunique(), 'Unique compounds', df['PUBCHEM_CID'].nunique(), 'Unique substances', df['PUBCHEM_SID'].nunique())

    '''
    df = df[df['PUBCHEM_CID'].notnull()]
    df['PUBCHEM_CID'] = df['PUBCHEM_CID'].astype(int)
    print('Compound ID is not null')
    print('Data potins', df.shape, 'Assays', df['PUBCHEM_AID'].nunique(), 'Unique compounds', df['PUBCHEM_CID'].nunique())
    '''

    cond = np.logical_or(df['PUBCHEM_CID'].notnull(), df['PUBCHEM_SID'].notnull())
    df = df[cond]
    print('Compound ID is not null OR Substance ID is not null')
    print('Data potins', df.shape, 'Assays', df['PUBCHEM_AID'].nunique(), 'Unique compounds', df['PUBCHEM_CID'].nunique(), 'Unique substances', df['PUBCHEM_SID'].nunique())

    df = df[df['PUBCHEM_ACTIVITY_OUTCOME'].isin([OUTCOME['Active'],OUTCOME['Inactive']])]
    print('OUTCOME is Active or Inactive')
    print('Data potins', df.shape, 'Assays', df['PUBCHEM_AID'].nunique(), 'Unique compounds', df['PUBCHEM_CID'].nunique(), 'Unique substances', df['PUBCHEM_SID'].nunique())

    pivot = pd.read_pickle('pkl/pivot.pkl')
    cond = np.logical_and(pivot['Active'] >= threshold, pivot['Inactive'] >= threshold)
    df = df[df['PUBCHEM_AID'].isin(pivot[cond].index)]
    print('Including greater than equal %d results' % threshold)
    print('Data potins', df.shape, 'Assays', df['PUBCHEM_AID'].nunique(), 'Unique compounds', df['PUBCHEM_CID'].nunique(), 'Unique substances', df['PUBCHEM_SID'].nunique())

    df.to_pickle('pkl/subset.pkl')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--download', help='download PCBA zip', action='store_true')
    parser.add_argument('--build', help='build PCAB descriptions and results pickle', action='store_true')
    parser.add_argument('--pivot', help='write pivot table to excel file', action='store_true')
    parser.add_argument('--subset', help='build subset for analysis', action='store_true')
    args = parser.parse_args()

    if args.download:
        download()
    elif args.build:
        build() 
    elif args.pivot:
        pivot()
    elif args.subset:
        subset()
    else:
        parser.print_help()
