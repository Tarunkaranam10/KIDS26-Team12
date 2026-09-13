#!/usr/bin/env python3
"""Prepare bounded beta matrices. Technical probe allowlist must be supplied for real modeling."""
import argparse,csv,pathlib,math,json
from acquire_tcga import digest
ROOT=pathlib.Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--master',type=pathlib.Path,default=ROOT/'data/processed/master_samples.tsv');p.add_argument('--gdc-dir',type=pathlib.Path,default=ROOT/'data/raw/gdc/beta');p.add_argument('--published-matrix',type=pathlib.Path);p.add_argument('--probes',type=pathlib.Path);p.add_argument('--smoke-probes',type=int,default=0);p.add_argument('--out',type=pathlib.Path,default=ROOT/'data/processed/beta.tsv');a=p.parse_args()
 if not a.probes and a.smoke_probes<=0:p.error('Supply --probes frozen technical allowlist, or --smoke-probes for engineering only')
 if a.probes and a.smoke_probes:p.error('Choose one probe-selection mode')
 allow=set(a.probes.read_text().split()) if a.probes else None
 rows=list(csv.DictReader(a.master.open(encoding='utf-8-sig'),delimiter='\t'));samples=[x['sample_id'] for x in rows]
 if len(set(samples))!=len(samples):raise ValueError('Duplicate samples')
 def valid(x):
  if x.lower() in ['na','nan','null','']:return 'NA'
  v=float(x)
  if not math.isfinite(v) or not 0<=v<=1:raise ValueError('Invalid beta')
  return x
 columns={};probes=[]
 if a.published_matrix:
  with a.published_matrix.open() as f:
   reader=csv.reader(f,delimiter='\t');header=next(reader);idx=[]
   for s in samples:
    candidates=[i for i,h in enumerate(header) if h[:16]==s[:16]]
    if len(candidates)!=1:raise ValueError('Missing/ambiguous matrix specimen: '+s)
    idx.append(candidates[0])
   a.out.parent.mkdir(parents=True,exist_ok=True)
   with a.out.open('w',newline='') as o:
    w=csv.writer(o,delimiter='\t');w.writerow(['probe_id']+samples)
    for row in reader:
     if not row[0].startswith('cg') or (allow is not None and row[0] not in allow):continue
     if row[0] in columns:raise ValueError('Duplicate probe')
     columns[row[0]]=True;w.writerow([row[0]]+[valid(row[i]) for i in idx]);probes.append(row[0])
     if a.smoke_probes and len(probes)>=a.smoke_probes:break
 else:
  for m in rows:
   f=a.gdc_dir/m['id']/m['filename']
   if digest(f)!=m['md5']:raise ValueError('Input MD5 mismatch: '+str(f))
   vals={}
   with f.open() as stream:
    for row in csv.reader(stream,delimiter='\t'):
     if not row or not row[0].startswith('cg'):continue
     if allow is not None and row[0] not in allow:continue
     if row[0] in vals:raise ValueError('Duplicate probe')
     vals[row[0]]=valid(row[1])
     if a.smoke_probes and len(vals)>=a.smoke_probes:break
   columns[m['sample_id']]=vals
  probes=sorted(set.union(*(set(x) for x in columns.values())))
  a.out.parent.mkdir(parents=True,exist_ok=True)
  with a.out.open('w',newline='') as f:
   w=csv.writer(f,delimiter='\t');w.writerow(['probe_id']+samples)
   for cg in probes:w.writerow([cg]+[columns[s].get(cg,'NA') for s in samples])
 if not probes:raise ValueError('No probes')
 a.out.with_suffix('.provenance.json').write_text(json.dumps({'engineering_only':bool(a.smoke_probes),'n_samples':len(samples),'n_probes':len(probes),'probe_allowlist_sha256':digest(a.probes,'sha256') if a.probes else None,'track':'historical-publication' if a.published_matrix else 'current-GDC-SeSAMe','sha256':digest(a.out,'sha256')},indent=2))
 print(len(samples),'samples',len(probes),'probes')
if __name__=='__main__':main()
