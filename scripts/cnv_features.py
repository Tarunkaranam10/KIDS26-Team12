#!/usr/bin/env python3
"""Explicitly total-CN proxies, NOT canonical scarHRD. Coordinates: 0-based half-open.

=============================================================================
scripts/cnv_features.py - Total-copy-number summary features from segments.
=============================================================================

USAGE
    python scripts/cnv_features.py <segments.tsv> <output.tsv> [--threshold 0.2]

INPUT SCHEMA (tab-separated, one row per segment)
    sample_id, chromosome, arm, start, end, log2r
        chromosome  autosomes only, '1' through '22'
        arm         'p' or 'q' - segments MUST already be split at the
                    centromere; this script will not do it for you
        start, end  0-based half-open genomic coordinates
        log2r       log2 ratio of observed to expected total copy number

WHY THE DOCSTRING SHOUTS "NOT canonical scarHRD"
-------------------------------------------------
This is the single most important thing to understand about this file, and the
reason the gated CNV branch exists at all.

Canonical HRDsum = HRD_LOH + LST + TAI. All three components require
ALLELE-AWARE data - you must be able to tell the maternal copy from the
paternal copy:

    HRD_LOH  loss of heterozygosity. A copy-neutral LOH event (2 copies of one
             allele, 0 of the other) has total copy number 2, exactly like a
             normal diploid region. Total CN alone literally cannot see it.
    TAI      telomeric allelic imbalance. "Imbalance" is a statement about the
             ratio between the two alleles. Total CN has no such concept.
    LST      large-scale state transitions, defined on allele-specific states.

This script computes features from TOTAL copy number only. The companion file
scripts/simulate_scars.py contains a machine-checked proof of the resulting
non-identifiability: signal(1,1,purity) and signal(2,0,purity) return the
IDENTICAL log2 ratio for every purity, so a balanced diploid region and a
copy-neutral LOH region are indistinguishable here.

Therefore the outputs below are PROXIES and covariates, never substitutes for
the reference labels. The most dangerous field is `methyl_LST_like_conservative`
- the "_like" and "_conservative" suffixes are load-bearing. It counts
large-segment transitions in total CN, which is a related but strictly weaker
signal than true allelic LST. Do not rename it to 'LST'. Do not compare it to
published LST values as though they measured the same thing.

ALSO CRITICAL: do not feed these features to a model that is predicting HRDsum
from methylation. The reference HRDsum labels are themselves derived from
allele-aware SNP6 segmentation - using CN-derived features as predictors of a
CN-derived label is circular and will produce impressive, meaningless numbers.
These belong in a clearly separated, explicitly gated comparison arm.

DESIGN NOTE: every fraction is normalised by `covered_bp` (the total length of
supplied segments) rather than by genome length. Two samples with different
assay coverage therefore remain comparable, but the denominator is coverage,
not the genome - which is why the field is named covered_bp and the fractions
carry the _covered suffix.
=============================================================================
"""
import argparse,csv,collections,math,pathlib,statistics

def features(segments,threshold=0.2,min_long=10_000_000,max_gap=3_000_000):
 """Compute total-CN summary features for ONE sample's segments.

 Parameters
     segments   iterable of dict-like rows for a single sample_id
     threshold  |log2r| at or above which a segment counts as altered. 0.2 is
                roughly a 15% copy-number change - conventional for array data,
                but it IS a tunable and results are sensitive to it.
     min_long   minimum segment length (10 Mb) for the LST-like transition
                count. Mirrors the "large-scale" in LST and suppresses the
                noisy short-segment calls that dominate raw array output.
     max_gap    maximum distance (3 Mb) between two segments for a transition
                between them to count. Prevents calling a "transition" across a
                large uncovered region where we simply have no data.

 Raises on any structural problem rather than silently skipping rows - a
 malformed segment file should stop the run, not quietly shrink the output.
 """
 if threshold<=0:raise ValueError('Positive amplitude threshold required')
 rows=[]
 for s in segments:
  # Copy and coerce. dict(s) avoids mutating the caller's row.
  r=dict(s);r['start']=int(r['start']);r['end']=int(r['end']);r['log2r']=float(r['log2r'])
  # Half-open interval means end must be strictly greater than start; a
  # zero-length segment is meaningless and would skew length-weighted sums.
  if r['start']<0 or r['end']<=r['start'] or not math.isfinite(r['log2r']):raise ValueError('Invalid segment')
  # Arms must already be split at the centromere. A segment spanning the
  # centromere would produce a spurious transition count, so this is a hard
  # requirement rather than something inferred here.
  if r['arm'] not in ['p','q']:raise ValueError('Provide centromere-split p/q arms')
  # Autosomes only: sex chromosomes have sample-dependent expected copy number
  # and would need separate handling.
  if str(r['chromosome']) not in [str(c) for c in range(1,23)]:raise ValueError('Autosomes 1-22 required')
  rows.append(r)
 if not rows:raise ValueError('No segments')

 # Deterministic genomic order; the adjacency scans below depend on it.
 rows.sort(key=lambda r:(int(r['chromosome']),r['start']))

 # Overlapping segments would double-count territory in every length-weighted
 # statistic, so reject rather than attempt a merge.
 for a,b in zip(rows,rows[1:]):
  if a['chromosome']==b['chromosome'] and a['end']>b['start']:raise ValueError('Overlapping segments')

 length=lambda r:r['end']-r['start']
 covered=sum(map(length,rows));altered=sum(length(r) for r in rows if abs(r['log2r'])>=threshold)

 # --- LST-like transition count ------------------------------------------
 # Count boundaries where TWO LONG neighbouring segments on the SAME ARM, close
 # enough together, differ in amplitude by at least the threshold.
 #
 # Conservative by construction: same-arm only (no centromere-spanning calls),
 # both sides long, gap bounded. This is a total-CN shadow of allelic LST and
 # will systematically miss copy-neutral events, which are invisible here.
 transitions=0
 for a,b in zip(rows,rows[1:]):
  if a['chromosome']==b['chromosome'] and a['arm']==b['arm'] and length(a)>=min_long and length(b)>=min_long and b['start']-a['end']<=max_gap and abs(a['log2r']-b['log2r'])>=threshold:transitions+=1

 return {'covered_bp':covered,                         # denominator for all fractions
 'segment_count':len(rows),                            # fragmentation proxy
 'FGA_covered':altered/covered,                        # fraction of genome altered
 'gain_fraction_covered':sum(length(r) for r in rows if r['log2r']>=threshold)/covered,
 'loss_fraction_covered':sum(length(r) for r in rows if r['log2r']<=-threshold)/covered,
 # Length-weighted mean |log2r|: unlike FGA this responds to the MAGNITUDE of
 # change, not just whether a threshold was crossed.
 'amplitude_burden':sum(length(r)*abs(r['log2r']) for r in rows)/covered,
 'median_segment_bp':statistics.median(map(length,rows)),
 'large_altered_segments':sum(length(r)>=min_long and abs(r['log2r'])>=threshold for r in rows),
 # See the docstring: this is NOT LST. The name is deliberately awkward.
 'methyl_LST_like_conservative':transitions}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('segments',type=pathlib.Path);p.add_argument('output',type=pathlib.Path);p.add_argument('--threshold',type=float,default=.2);a=p.parse_args()
 # Group rows by sample. utf-8-sig strips a UTF-8 BOM if the file came from
 # Excel or a Windows tool, which would otherwise corrupt the first column name.
 groups=collections.defaultdict(list)
 with a.segments.open(encoding='utf-8-sig') as f:
  for r in csv.DictReader(f,delimiter='\t'):groups[r['sample_id']].append(r)
 # sorted() makes output row order deterministic across runs.
 out=[{'sample_id':s,**features(v,a.threshold)} for s,v in sorted(groups.items())]
 if not out:raise ValueError('Empty input')
 a.output.parent.mkdir(parents=True,exist_ok=True)
 with a.output.open('w',newline='') as f:
  # lineterminator='\n' pins LF. csv's default is '\r\n', which previously
  # leaked CRLF into generated tables in this repository.
  w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
if __name__=='__main__':main()
