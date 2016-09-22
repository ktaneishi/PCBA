import pandas as pd
import urllib.request
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import gzip
import re
import os
import sys
import time

def download(outdir='data/compound'):
    base = 'http://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/%s'

    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    for line in urllib.request.urlopen(base % '').readlines():
        filename = re.search('Compound_\d+_\d+.sdf.gz', str(line))
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
        time.sleep(10)

def ecfp(indir='data/compound', outdir='data/ecfp', radius=2, nBits=1024):
    outdir = '%s_%d_%d' % (outdir, radius*2, nBits)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for filename in sorted(os.listdir(indir)):
        print(filename)
        if os.path.exists(os.path.join(outdir, filename.replace('sdf.gz','pkl'))):
            continue

        start, end = map(int, filename.replace('.sdf.gz','').split('_')[1:3])

        if start == 102125001:
            continue

        cid, fp = [], []
        f = gzip.open(os.path.join(indir, filename))
        for mol in Chem.ForwardSDMolSupplier(f):
            if not mol:
                continue
            cid.append(int(mol.GetProp('PUBCHEM_COMPOUND_CID')))
            fp.append(AllChem.GetHashedMorganFingerprint(
                    mol, radius=radius, nBits=nBits, useChirality=False,
                    useBondTypes=True, useFeatures=False))
        f.close()
        fp = pd.DataFrame.from_dict(fp)
        fp.index = cid
        fp.to_pickle(os.path.join(outdir, filename.replace('sdf.gz','pkl')))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--download', help='download compounds gzipped SDF', action='store_true')
    parser.add_argument('--ecfp', help='generate ECFP from compounds', action='store_true')
    parser.add_argument('--radius', help='ECFP radius', action='store', type=int, default=2)
    parser.add_argument('--nBits', help='ECFP number of bits', action='store', type=int, default=1024)
    args = parser.parse_args()

    if args.download:
        download()
    if args.ecfp:
        ecfp(radius=args.radius, nBits=args.nBits)
    else:
        parser.print_help()
