#!/usr/bin/env python3
"""Build metadata from the downloaded historical 450K matrix header without loading the matrix."""
import argparse,csv,pathlib,collections
from acquire_tcga import write_tsv
from build_master import barcode,read
ROOT=pathlib.Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--matrix',type=pathlib.Path,default=ROOT/'data/raw/pancanatlas/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv');p.add_argument('--published',type=pathlib.Path,default=ROOT/'data/raw/pancanatlas');p.add_argument('--case-map',type=pathlib.Path,default=ROOT/'config/published_case_projects.tsv');p.add_argument('--out',type=pathlib.Path,default=ROOT/'config/published_beta_metadata.tsv');a=p.parse_args()
 with a.matrix.open() as f:header=next(csv.reader(f,delimiter='\t'))
 sub=collections.defaultdict(set)
 aliases={'OVCA':'OV','AML':'LAML'}
 for row in read(a.published/'TCGASubtype.20170308.tsv'):sub[row['pan.samplesID'][:12]].add(aliases.get(row['cancer.type'],row['cancer.type']))
 case_map={}
 if a.case_map.exists():
  for row in read(a.case_map):
   if row['patient_id'] in case_map and case_map[row['patient_id']]!=row['cancer_type']:raise ValueError('Conflicting GDC project mapping')
   case_map[row['patient_id']]=row['cancer_type']
 manifest=read(ROOT/'config/panimmune_source_manifest.tsv');source=next(r for r in manifest if r['filename']==a.matrix.name)
 rows=[]
 for sample in header:
  if not sample.startswith('TCGA-'):continue
  patient,_,kind=barcode(sample)
  if len(sample)<16:raise ValueError('Matrix lacks vial resolution; adjudicate before joining')
  cancer=case_map.get(patient,next(iter(sub[patient])) if len(sub[patient])==1 else 'UNRESOLVED')
  cancer_source='GDC_case_project' if patient in case_map else ('PanCanAtlas_subtype' if cancer!='UNRESOLVED' else 'unresolved')
  rows.append({'id':source['id'],'filename':source['filename'],'md5':source['md5'],'size':source['size'],'patient_id':patient,'sample_id':sample,'cancer_type':cancer,'cancer_type_source':cancer_source,'platform':'Illumina Human Methylation 450','workflow':'PanCanAtlas_historical_published_beta','access':'open'})
 if not rows:raise ValueError('No TCGA barcodes found in matrix header')
 write_tsv(a.out,rows,list(rows[0]));print(len(rows),'matrix samples indexed; use build_master --metadata',a.out)
if __name__=='__main__':main()
