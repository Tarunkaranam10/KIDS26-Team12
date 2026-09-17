#!/usr/bin/env python3
"""Build source-linked evidence matrix; curated JSON supplies reviewed facts."""
import csv,json,pathlib,re,argparse
FIELDS=['paper_id','Paper','Year','Category','Link','local_file','pages_sections','review_status','delta_category','Cancer_cohort','N','Assay','Input_datatype','Ground_truth_HRD_definition','HRD_calculation_method','LOH_definition','LST_definition','TAI_definition','Upstream_allele_specific_CN_method','Purity_ploidy_correction','CNV_method','Methylation_preprocessing','Raw_IDAT_MU_requirement','RNA_preprocessing','Model_class','Feature_selection','Training_strategy','Validation_strategy','Cross_validation_design','External_validation','Performance_metrics','Calibration_uncertainty','Single_sample_compatible','Cross_cancer_compatible','Code_software','Synthetic_data_relevance','Major_strengths','Major_limitations','Direct_relevance_KIDS26','Implementation_ideas','documented_fact','inference','recommendation']
NEW_METHOD_FIELDS=['LOH_definition','LST_definition','TAI_definition','Upstream_allele_specific_CN_method','Purity_ploidy_correction','Raw_IDAT_MU_requirement','Calibration_uncertainty']
def main():
 p=argparse.ArgumentParser();p.add_argument('--source',type=pathlib.Path,default=pathlib.Path('..'));a=p.parse_args();root=pathlib.Path(__file__).resolve().parents[1]
 notes={}
 for f in sorted((root/'config').glob('literature_notes*.json')):notes.update(json.loads(f.read_text(encoding='utf-8-sig')))
 files=list((a.source/'papers').glob('*.pdf'))
 def norm(s):return re.sub('[^a-z0-9]','',s.lower())
 rows=[]
 for n,src in enumerate(csv.DictReader((a.source/'HRD_methylation_classifier_literature_map.csv').open(encoding='utf-8-sig')),1):
  row={k:'Not established in reviewed material' for k in FIELDS};row.update(paper_id=f'L{n:02}',Paper=src['Title'],Year=src['Year'],Category=src['Category'],Link=src['Link'],local_file='',pages_sections='',review_status='Unavailable locally; metadata only; no methods inferred')
  matches=[f for f in files if f.stem[:4]==src['Year'] and norm(f.stem[5:])[:60]==norm(src['Title'])[:60]]
  if matches:row.update(local_file='../papers/'+matches[0].name,review_status='Local abstract and relevant methods/results reviewed; supplementary gaps retained')
  if str(n) in notes:row.update(notes[str(n)])
  else:row.update(documented_fact='Listed in supplied map; full text unavailable locally.',inference='No methodological inference from title.',recommendation='Verify relevant methods before relying on this reference.',Direct_relevance_KIDS26=src['Category'])
  if 'Newly available' in row['review_status']:
   for key in NEW_METHOD_FIELDS:
    if row[key]=='Not established in reviewed material':row[key]='Absent from reviewed paper or not reported'
  rows.append(row)
 extras_path=root/'config/literature_extras.json'
 if extras_path.exists():
  for extra in json.loads(extras_path.read_text(encoding='utf8')):
   row={k:'Not established in reviewed material' for k in FIELDS};row.update(extra)
   for key in NEW_METHOD_FIELDS:
    if row[key]=='Not established in reviewed material':row[key]='Absent from reviewed paper or not reported'
   rows.append(row)
 with (root/'docs/literature_evidence.csv').open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=FIELDS,lineterminator='\n');w.writeheader();w.writerows(rows)
 (root/'docs/literature_evidence.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False),encoding='utf8')
 lines=['# Literature evidence: paper-by-paper digest','Evidence labels distinguish source claims from KIDS26 judgments. Unknown is not a negative finding. Page numbers are physical PDF pages including covers. Original CSV titles/years are retained, with discrepancies noted. This is a methods-focused first pass, not a systematic review or a complete supplementary-methods audit.']
 for row in rows:lines += [f"## {row['paper_id']} - {row['Paper']}",f"Source: {row['local_file'] or row['Link']}; {row['pages_sections']}".rstrip(),f"Status: {row['review_status']}",f"Delta category: {row['delta_category']}",f"**Documented fact:** {row['documented_fact']}",f"**Inference:** {row['inference']}",f"**Recommendation:** {row['recommendation']}"]
 (root/'docs/01_LITERATURE_EVIDENCE.md').write_text('\n\n'.join(lines),encoding='utf8');print(len(rows),'rows;',sum(bool(r['local_file']) for r in rows),'local PDFs mapped')
if __name__=='__main__':main()
