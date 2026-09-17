#!/usr/bin/env python3
"""Render an auditable cohort flow and per-cancer HRDsum table.

=============================================================================
scripts/audit_cohort.py - Turn the matching audit + master table into a
                          human-readable cohort QC report.
=============================================================================

USAGE
    python scripts/audit_cohort.py \
        --audit  data/processed/matching_audit.tsv \
        --master data/processed/master_samples.tsv \
        [--report docs/16_COHORT_QC.md] \
        [--by-cancer docs/cohort_qc_by_cancer.csv]

PURPOSE
-------
scripts/build_master.py emits two files: matching_audit.tsv (one row for EVERY
candidate matrix column, with its disposition) and master_samples.tsv (only the
rows that survived). This script reads both and renders the flow between them
as Markdown, so a reader can see exactly how 9,664 candidate columns became
7,707 analysable specimens, and why each excluded column was dropped.

This is a reporting script. It computes no model inputs and makes no cohort
decisions - it only describes decisions already made upstream. Re-run it
whenever the master table is rebuilt.

WHY THE EXCLUSION TABLE MATTERS
--------------------------------
build_master.py applies its exclusion rules in a strict precedence order, so
each column appears exactly once in the audit with exactly one reason. That
property is what makes the counts below sum correctly and makes the cohort
reproducible from the inputs. If you change the rule order upstream, these
numbers change and this report must be regenerated.

KNOWN ISSUES (audited 2026-09-16)
----------------------------------
  * The narrative line "beta values have not been downloaded" is hardcoded in
    the report text below and is NOW FALSE - data/processed/beta.tsv exists
    (26 GB, 336,480 probes x 7,707 samples). Update that sentence.
  * quantile() here is a hand-rolled linear-interpolation implementation with
    no unit test anywhere, and it feeds the published per-cancer table.
  * Line building `types` slices sample_id[13:15] positionally without going
    through build_master.barcode(), so a malformed ID is silently bucketed as
    non-primary rather than raising.
  * rows[0] is dereferenced without an emptiness guard; an empty master table
    would raise IndexError rather than a useful message.
  * There is no check that --audit and --master came from the SAME
    build_master.py run. Mismatched inputs would produce a plausible-looking
    but incoherent report.
=============================================================================
"""
import argparse,csv,collections,pathlib,statistics
ROOT=pathlib.Path(__file__).resolve().parents[1]

def read(path):
 """Read a TSV into a list of dicts. utf-8-sig strips a BOM if present."""
 with path.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))

def quantile(values,p):
 """Linear-interpolated quantile, p in [0,1].

 Matches numpy's default 'linear' method. Implemented by hand because this
 project is deliberately standard-library-only. NOT covered by any test - see
 the KNOWN ISSUES note above.
 """
 x=sorted(values);pos=(len(x)-1)*p;lo=int(pos);hi=min(lo+1,len(x)-1);return x[lo]+(x[hi]-x[lo])*(pos-lo)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--audit',type=pathlib.Path,required=True);p.add_argument('--master',type=pathlib.Path,required=True);p.add_argument('--report',type=pathlib.Path,default=ROOT/'docs/16_COHORT_QC.md');p.add_argument('--by-cancer',type=pathlib.Path,default=ROOT/'docs/cohort_qc_by_cancer.csv');a=p.parse_args()

 # audit  = every candidate column and its fate; master = the survivors.
 # reasons = histogram of exclusion causes, over excluded rows only.
 # types   = histogram of TCGA sample-type codes taken from characters 13:15 of
 #           the barcode ('01' = primary tumour, '11' = normal, and so on).
 audit,master=read(a.audit),read(a.master);reasons=collections.Counter(x['reason'] for x in audit if x['status']=='excluded');types=collections.Counter(x['sample_id'][13:15] for x in audit);groups=collections.defaultdict(list)

 # Collect HRDsum values per cancer type from the INCLUDED cohort only.
 for x in master:groups[x['cancer_type']].append(float(x['HRDsum']))

 # --- Per-cancer distribution table ---------------------------------------
 # Five-number summary plus a count above 42. Reporting the full spread rather
 # than a mean matters here: HRDsum is strongly right-skewed in most tissues,
 # so a mean alone would be misleading.
 rows=[]
 for cancer,vals in sorted(groups.items()):rows.append({'cancer_type':cancer,'N':len(vals),'HRDsum_min':min(vals),'HRDsum_Q1':quantile(vals,.25),'HRDsum_median':statistics.median(vals),'HRDsum_Q3':quantile(vals,.75),'HRDsum_max':max(vals),'HRDsum_ge_42':sum(v>=42 for v in vals)})
 a.by_cancer.parent.mkdir(parents=True,exist_ok=True)
 with a.by_cancer.open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

 # --- Missingness in the included cohort ----------------------------------
 # Counts empty/whitespace-only cells per field. purity and ploidy are the ones
 # that actually vary (178 rows lack ABSOLUTE calls); the label fields should be
 # complete by construction, so a non-zero count there indicates an upstream bug.
 # NOTE: an absent COLUMN and an empty CELL are counted identically here.
 miss={k:sum(not x.get(k,'').strip() for x in master) for k in ['HRDsum','HRD_LOH','LST','TAI','purity','ploidy','cancer_type']}

 # --- Build the Markdown report -------------------------------------------
 # The construction flow reads top to bottom as a funnel: all columns -> unique
 # participants -> primary tumours -> included -> of which locked.
 # FIXME: the "beta values have not been downloaded" clause below is stale.
 lines=['# Historical TCGA cohort QC','','This report is generated from the bounded real matrix header plus downloaded public labels/QC. It is a metadata-level audit; beta values have not been downloaded.','','## Construction flow','',f'- Matrix columns starting with TCGA: **{len(audit):,}**',f'- Unique participants: **{len({x["patient_id"] for x in audit}):,}**',f'- Primary-tumor columns (sample type 01): **{types["01"]:,}**',f'- Included after precedence-ordered rules: **{len(master):,}**',f'- Locked CNS (GBM/LGG): **{sum(x["partition"]=="locked_CNS" for x in master):,}**','','Exclusions are precedence ordered, so each column appears once:','']

 # Most common exclusion reason first, so the dominant attrition cause is
 # immediately visible.
 lines += [f'- `{k}`: {v:,}' for k,v in reasons.most_common()]

 # This paragraph documents two genuinely important subtleties: (1) a patient
 # with extra normal/recurrent arrays is NOT disqualified - only ambiguity among
 # PRIMARY tumours excludes them; (2) the HRD labels key on a 15-character
 # sample-type ID while the methylation matrix carries full 16+ character
 # aliquots, so every join is a coarser sample-level match, not proof of
 # identical physical material.
 lines += ['','Multiple normal/recurrent/metastatic arrays do not disqualify a unique primary. Multiple primary arrays/specimens for a participant remain excluded pending deterministic QC adjudication. Published HRD labels resolve to 15-character sample type, while the methylation matrix retains full aliquots; this is recorded as a coarser match.','','## Included-cohort missingness','']+[f'- {k}: {v:,} / {len(master):,}' for k,v in miss.items()]

 lines += ['','## Per-cancer HRDsum distribution','','| Cancer | N | Min | Q1 | Median | Q3 | Max | N >= 42 |','|---|---:|---:|---:|---:|---:|---:|---:|']
 for x in rows:lines.append(f'| {x["cancer_type"]} | {x["N"]} | {x["HRDsum_min"]:.1f} | {x["HRDsum_Q1"]:.1f} | {x["HRDsum_median"]:.1f} | {x["HRDsum_Q3"]:.1f} | {x["HRDsum_max"]:.1f} | {x["HRDsum_ge_42"]} |')

 # Explicit disclaimer on the 42 cutoff. It appears in the literature as an
 # HRD-high threshold for ADULT cancers under specific assay conditions; it is
 # used here purely as a descriptive column and must not be read as a validated
 # pan-cancer or pediatric decision boundary.
 lines += ['','The value 42 is descriptive here and is not assumed to be a validated pediatric or pan-cancer cutoff. Final cohort release still requires beta-row extraction, platform/probe QC, and review of exclusions.']
 a.report.write_text('\n'.join(lines)+'\n',encoding='utf-8');print(f'{len(master)} included; {len(rows)} cancers')
if __name__=='__main__':main()
