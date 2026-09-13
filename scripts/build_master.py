#!/usr/bin/env python3
"""Audited TCGA specimen-level join. Ambiguous replicates are excluded, never averaged."""
import argparse,csv,collections,pathlib,re,math,json
from acquire_tcga import write_tsv
ROOT=pathlib.Path(__file__).resolve().parents[1]
def read(p):
 with pathlib.Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def barcode(s):
 if not re.fullmatch(r'TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}(?:-[0-9]{2}[A-Z]?(?:-[A-Z0-9]+)*)?',s):raise ValueError('Invalid TCGA barcode: '+s)
 return s[:12],s[:15],s[13:15] if len(s)>=15 else ''
def num(s):
 try:
  v=float(s);return v if math.isfinite(v) else None
 except (ValueError,TypeError):return None
def build(metadata,labels,subtypes,quality,purity):
 targets=collections.defaultdict(list);types=collections.defaultdict(set);qc=collections.defaultdict(list);pur=collections.defaultdict(list)
 for x in labels:
  patient,key,kind=barcode(x['sampleID']);targets[key].append(x)
 aliases={'OVCA':'OV','AML':'LAML'}
 for x in subtypes:types[x['pan.samplesID'][:12]].add(aliases.get(x['cancer.type'],x['cancer.type']))
 for x in quality:
  if 'methylation450' in x['platform'].lower():qc[x['aliquot_barcode'][:16]].append(x)
 for x in purity:pur[x['sample'][:15]].append(x)
 # A normal, recurrent, or metastatic array must not make an otherwise unique
 # primary tumor look duplicated. Count primary-tumor candidates only; retain
 # the stricter rule that multiple primary arrays/specimens need adjudication.
 primary_counts=collections.Counter(
  x['patient_id'] for x in metadata if barcode(x['sample_id'])[2]=='01')
 out=[];audit=[]
 for m in metadata:
  patient,key,kind=barcode(m['sample_id']);reason='';target=targets[key];q=qc[m['sample_id'][:16]]
  if patient!=m['patient_id']:reason='patient_metadata_mismatch'
  elif kind!='01':reason='not_primary_tumor'
  elif primary_counts[patient]!=1:reason='ambiguous_multiple_primary_files_or_specimens_per_patient'
  elif len(target)!=1:reason='missing_or_ambiguous_HRD_label'
  elif not m.get('cancer_type') or m['cancer_type']=='UNRESOLVED':reason='missing_tumor_type'
  elif types[patient] and (len(types[patient])!=1 or m['cancer_type'] not in types[patient]):reason='conflicting_tumor_type'
  elif not q:reason='missing_published_methylation_quality_annotation'
  elif any(x['Do_not_use'].lower()=='true' or x['AWG_excluded_because_of_pathology']=='1.0' for x in q):reason='quality_excluded'
  if not reason:
   vals=[num(target[0][x]) for x in ['hrd-loh','lst1','ai1','HRD']]
   if any(v is None or v<0 for v in vals):reason='invalid_HRD_label'
   elif abs(sum(vals[:3])-vals[3])>1e-6:reason='component_sum_mismatch'
  audit.append({**m,'status':'excluded' if reason else 'included','reason':reason or 'unique_primary_specimen_match;label_resolves_sample_type_not_vial'})
  if reason:continue
  pp=pur[key];pp=pp[0] if len(pp)==1 else {}
  out.append({**m,'HRD_LOH':vals[0],'LST':vals[1],'TAI':vals[2],'HRDsum':vals[3],'purity':pp.get('purity',''),'ploidy':pp.get('ploidy',''),'label_source':'PanImmune:66dd07d7-6366-4774-83c3-5ad1e22b177e','match_resolution':'sample_type_15char_unique_methylation_specimen','quality_annotation':'published_450K_no_exclusion','partition':'locked_CNS' if m['cancer_type'] in ['GBM','LGG'] else 'development'})
 return out,audit
def main():
 p=argparse.ArgumentParser();p.add_argument('--metadata',type=pathlib.Path,default=ROOT/'config/gdc_beta_dev.tsv');p.add_argument('--published',type=pathlib.Path,default=ROOT/'data/raw/pancanatlas');p.add_argument('--out',type=pathlib.Path,default=ROOT/'data/processed');a=p.parse_args()
 rows,audit=build(read(a.metadata),read(a.published/'TCGA.HRD_withSampleID.txt'),read(a.published/'TCGASubtype.20170308.tsv'),read(a.published/'merged_sample_quality_annotations.tsv'),read(a.published/'TCGA_mastercalls.abs_tables_JSedit.fixed.txt'))
 write_tsv(a.out/'matching_audit.tsv',audit,list(audit[0]) if audit else ['status','reason'])
 if not rows:raise SystemExit('No unambiguous eligible matches. Inspect matching_audit.tsv; do not weaken joins silently.')
 write_tsv(a.out/'master_samples.tsv',rows,list(rows[0]));print(json.dumps({'included':len(rows),'excluded':len(audit)-len(rows),'cancers':dict(collections.Counter(r['cancer_type'] for r in rows))}))
if __name__=='__main__':main()
