#!/usr/bin/env python3
"""Render an auditable cohort flow and per-cancer HRDsum table."""
import argparse,csv,collections,pathlib,statistics
ROOT=pathlib.Path(__file__).resolve().parents[1]
def read(path):
 with path.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))
def quantile(values,p):
 x=sorted(values);pos=(len(x)-1)*p;lo=int(pos);hi=min(lo+1,len(x)-1);return x[lo]+(x[hi]-x[lo])*(pos-lo)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--audit',type=pathlib.Path,required=True);p.add_argument('--master',type=pathlib.Path,required=True);p.add_argument('--report',type=pathlib.Path,default=ROOT/'docs/16_COHORT_QC.md');p.add_argument('--by-cancer',type=pathlib.Path,default=ROOT/'docs/cohort_qc_by_cancer.csv');a=p.parse_args()
 audit,master=read(a.audit),read(a.master);reasons=collections.Counter(x['reason'] for x in audit if x['status']=='excluded');types=collections.Counter(x['sample_id'][13:15] for x in audit);groups=collections.defaultdict(list)
 for x in master:groups[x['cancer_type']].append(float(x['HRDsum']))
 rows=[]
 for cancer,vals in sorted(groups.items()):rows.append({'cancer_type':cancer,'N':len(vals),'HRDsum_min':min(vals),'HRDsum_Q1':quantile(vals,.25),'HRDsum_median':statistics.median(vals),'HRDsum_Q3':quantile(vals,.75),'HRDsum_max':max(vals),'HRDsum_ge_42':sum(v>=42 for v in vals)})
 a.by_cancer.parent.mkdir(parents=True,exist_ok=True)
 with a.by_cancer.open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 miss={k:sum(not x.get(k,'').strip() for x in master) for k in ['HRDsum','HRD_LOH','LST','TAI','purity','ploidy','cancer_type']}
 lines=['# Historical TCGA cohort QC','','This report is generated from the bounded real matrix header plus downloaded public labels/QC. It is a metadata-level audit; beta values have not been downloaded.','','## Construction flow','',f'- Matrix columns starting with TCGA: **{len(audit):,}**',f'- Unique participants: **{len({x["patient_id"] for x in audit}):,}**',f'- Primary-tumor columns (sample type 01): **{types["01"]:,}**',f'- Included after precedence-ordered rules: **{len(master):,}**',f'- Locked CNS (GBM/LGG): **{sum(x["partition"]=="locked_CNS" for x in master):,}**','','Exclusions are precedence ordered, so each column appears once:','']
 lines += [f'- `{k}`: {v:,}' for k,v in reasons.most_common()]
 lines += ['','Multiple normal/recurrent/metastatic arrays do not disqualify a unique primary. Multiple primary arrays/specimens for a participant remain excluded pending deterministic QC adjudication. Published HRD labels resolve to 15-character sample type, while the methylation matrix retains full aliquots; this is recorded as a coarser match.','','## Included-cohort missingness','']+[f'- {k}: {v:,} / {len(master):,}' for k,v in miss.items()]
 lines += ['','## Per-cancer HRDsum distribution','','| Cancer | N | Min | Q1 | Median | Q3 | Max | N >= 42 |','|---|---:|---:|---:|---:|---:|---:|---:|']
 for x in rows:lines.append(f'| {x["cancer_type"]} | {x["N"]} | {x["HRDsum_min"]:.1f} | {x["HRDsum_Q1"]:.1f} | {x["HRDsum_median"]:.1f} | {x["HRDsum_Q3"]:.1f} | {x["HRDsum_max"]:.1f} | {x["HRDsum_ge_42"]} |')
 lines += ['','The value 42 is descriptive here and is not assumed to be a validated pediatric or pan-cancer cutoff. Final cohort release still requires beta-row extraction, platform/probe QC, and review of exclusions.']
 a.report.write_text('\n'.join(lines)+'\n',encoding='utf-8');print(f'{len(master)} included; {len(rows)} cancers')
if __name__=='__main__':main()
