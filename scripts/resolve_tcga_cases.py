#!/usr/bin/env python3
"""Resolve TCGA participant barcodes to authoritative GDC project IDs."""
import argparse,csv,json,pathlib
from acquire_tcga import API,request,write_tsv
ROOT=pathlib.Path(__file__).resolve().parents[1]

def participants_from_matrix(matrix):
 with matrix.open(encoding='utf-8-sig') as f:header=next(csv.reader(f,delimiter='\t'))
 ids=sorted({x[:12] for x in header if x.startswith('TCGA-')})
 if not ids:raise ValueError('No TCGA barcodes found in matrix header')
 return ids

def resolve(ids,chunk_size=400):
 rows=[]
 for start in range(0,len(ids),chunk_size):
  chunk=ids[start:start+chunk_size]
  payload={'filters':{'op':'in','content':{'field':'submitter_id','value':chunk}},
           'format':'JSON','size':len(chunk),'fields':'submitter_id,project.project_id'}
  hits=json.loads(request(API+'/cases',payload))['data']['hits']
  for h in hits:
   project=h.get('project',{}).get('project_id','')
   if project.startswith('TCGA-'):
    rows.append({'patient_id':h['submitter_id'],'cancer_type':project.removeprefix('TCGA-'),'project_id':project})
 return sorted(rows,key=lambda x:x['patient_id'])

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--matrix',type=pathlib.Path,default=ROOT/'data/raw/pancanatlas/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv')
 p.add_argument('--out',type=pathlib.Path,default=ROOT/'config/published_case_projects.tsv')
 a=p.parse_args();ids=participants_from_matrix(a.matrix);rows=resolve(ids)
 write_tsv(a.out,rows,['patient_id','cancer_type','project_id'])
 missing=sorted(set(ids)-{x['patient_id'] for x in rows})
 print(json.dumps({'matrix_participants':len(ids),'resolved':len(rows),'unresolved':len(missing)}))
 if missing:print('First unresolved:',','.join(missing[:10]))
if __name__=='__main__':main()
