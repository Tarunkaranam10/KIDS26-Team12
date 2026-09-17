#!/usr/bin/env python3
"""Analytic identifiability fixtures plus noisy probe-grid observations. No clinical data.

=============================================================================
scripts/simulate_scars.py - Synthetic allele-specific fixtures that demonstrate
                            what total copy number can and cannot recover.
=============================================================================

USAGE
    python scripts/simulate_scars.py [--out results/simulation] [--seed 26]
                                     [--noise 0.05] [--step-mb 1]
                                     [--breakpoint-shift-mb 0]

OUTPUTS (three files in --out)
    allelic_truth.tsv            ground-truth segments with major/minor CN
    probe_observations.tsv       simulated noisy log2r + minor-allele fraction
                                 on a regular probe grid, at four purities
    expected_and_limitations.json hand-reasoned expected scar counts, plus a
                                 machine-checked non-identifiability result

WHY THIS FILE EXISTS
--------------------
It is a self-contained, zero-dependency demonstration of the central
methodological constraint of the whole project: methylation-derived TOTAL copy
number cannot recover the allele-aware quantities that HRDsum is built from.

The proof is the `nonidentifiability` field written at the end. It asserts that
signal(1,1,purity)[0] == signal(2,0,purity)[0] for every purity tested - a
normal balanced diploid region and a copy-neutral LOH region produce the
IDENTICAL log2 ratio. No amount of downstream cleverness can separate them
using total CN alone, because the information is not present in the signal.

Use this when someone proposes deriving LOH or TAI from beta-value-inferred
copy number. It takes seconds to run and settles the argument.

NO REAL DATA IS INVOLVED. Every number here is generated from the toy model
below. The fixtures are hand-reasoned, NOT outputs of a pinned scarHRD run -
see the `note` field in the JSON, which says so explicitly.
=============================================================================
"""
import argparse,csv,json,math,pathlib,random

def signal(a,b,purity):
 """Forward model: allele-specific copy number -> observed array signal.

 Parameters
     a, b    major and minor allele copy number (non-negative integers)
     purity  tumour fraction in [0,1]; the rest is normal diploid contamination

 Returns
     (log2_ratio, minor_allele_fraction)

 The mixture: observed total copy number is the purity-weighted tumour
 contribution purity*(a+b) plus the normal contribution 2*(1-purity). Dividing
 by the diploid expectation of 2 and taking log2 gives the standard log2 ratio.

 THE KEY PROPERTY: the first return value depends on a and b ONLY through their
 SUM. So (a=1,b=1) and (a=2,b=0) - balanced diploid versus copy-neutral LOH -
 are mathematically indistinguishable in log2r at every purity. The second
 return value, the minor-allele fraction, DOES separate them, which is exactly
 why allele-aware assays are required for LOH and TAI.
 """
 if a<0 or b<0 or not 0<=purity<=1:raise ValueError('Invalid allele count or purity')
 total=purity*(a+b)+2*(1-purity)
 if total<=0:return None,None
 return math.log2(total/2),(purity*b+(1-purity))/total

# Chromosome 1 toy length 100 Mb; centromere [45,55) Mb. A/B are allele-specific CN.
#
# Each fixture is a list of (start_Mb, end_Mb, major_CN, minor_CN) and is chosen
# to isolate one specific discrimination problem:
#
#   balanced               negative control: no scars of any kind.
#   copy_neutral_LOH       THE critical case. total CN = 2 everywhere, identical
#                          to 'balanced' in log2r, but minor CN = 0 on the first
#                          segment. Invisible to total-CN methods.
#   balanced_gain          total CN changes (2 -> 4) with alleles still balanced:
#                          a real LST-like boundary with no LOH and no AI.
#   telomeric_AI           allelic imbalance (2:1) touching the p-telomere -
#                          the defining TAI configuration.
#   whole_chromosome_LOH   LOH across the ENTIRE chromosome. Counts as 0 for
#                          scar purposes because canonical definitions exclude
#                          whole-chromosome events (they arise from missegregation,
#                          not from the double-strand-break repair failure that
#                          HRD measures). A naive implementation would wrongly
#                          count this - hence the fixture.
#   centromere_crossing_AI imbalance on both arms. Tests that transitions are
#                          NOT counted across the centromere.
FIXTURES={
 'balanced':[(0,45,1,1),(55,100,1,1)],
 'copy_neutral_LOH':[(0,20,2,0),(20,45,1,1),(55,100,1,1)],
 'balanced_gain':[(0,20,2,2),(20,45,1,1),(55,100,1,1)],
 'telomeric_AI':[(0,20,2,1),(20,45,1,1),(55,100,1,1)],
 'whole_chromosome_LOH':[(0,45,2,0),(55,100,2,0)],
 'centromere_crossing_AI':[(0,30,1,1),(30,45,2,1),(55,70,2,1),(70,100,1,1)]}

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=pathlib.Path,default=pathlib.Path('results/simulation'));p.add_argument('--seed',type=int,default=26);p.add_argument('--noise',type=float,default=.05);p.add_argument('--step-mb',type=float,default=1);p.add_argument('--breakpoint-shift-mb',type=float,default=0);a=p.parse_args()
 if a.step_mb<=0 or a.noise<0 or abs(a.breakpoint_shift_mb)>5:raise ValueError('Invalid grid/noise/shift')
 # Seeded RNG: identical inputs always produce identical outputs.
 a.out.mkdir(parents=True,exist_ok=True);rng=random.Random(a.seed);segs=[];obs=[]

 for name,intervals in FIXTURES.items():
  # Ground truth: exact segment boundaries with allele-specific copy numbers.
  for start,end,major,minor in intervals:
   segs.append(dict(sample_id=name,chromosome='1',start=int(start*1e6),end=int(end*1e6),total_cn=major+minor,major_cn=major,minor_cn=minor))

  # Simulated observations at four purities. Low purity compresses log2r toward
  # zero (normal contamination dilutes the tumour signal), so 0.2 shows how
  # quickly real scars become undetectable in impure samples.
  for purity in [.2,.5,.8,1.]:
   x=a.step_mb/2   # sample at bin centres, not edges
   while x<100:
    # Probe-grid degradation and shifted observation relative to true breakpoints.
    # --breakpoint-shift-mb deliberately misaligns the observation grid from the
    # true boundaries, simulating the real situation where array probes do not
    # sit exactly on breakpoints. Useful for testing boundary robustness.
    lookup=x-a.breakpoint_shift_mb
    hit=next((v for v in intervals if v[0]<=lookup<v[1]),None)
    if hit:
     ratio,baf=signal(hit[2],hit[3],purity)
     # Gaussian noise on log2r mimics array measurement error. Note that
     # ideal_minor_fraction is stored WITHOUT noise: it is the theoretical
     # allele fraction, recorded to show what an allele-aware assay would see.
     obs.append(dict(sample_id=name,purity=purity,position=int(x*1e6),log2r=ratio+rng.gauss(0,a.noise),ideal_minor_fraction=baf))
    x+=a.step_mb

 for name,rows in [('allelic_truth.tsv',segs),('probe_observations.tsv',obs)]:
  with (a.out/name).open('w',newline='') as f:
   # lineterminator='\n' pins LF; csv's default '\r\n' previously leaked CRLF
   # into generated tables in this repository.
   w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

 # --- Expected values and their caveats ------------------------------------
 # These counts are HAND-REASONED from the canonical scar definitions, not
 # produced by running scarHRD. The `note` field says so, and it should stay
 # there: treating these as verified ground truth would be exactly the kind of
 # unearned confidence this project is trying to avoid.
 expected={'coordinate_system':'0-based half-open; toy chr1 length 100Mb; centromere 45-55Mb',
 'balanced':{'LOH':0,'TAI':0,'LST':0},'copy_neutral_LOH':{'LOH':1,'TAI':1,'LST_allelic_boundary':1},
 'balanced_gain':{'LOH':0,'TAI':0,'LST':1},'telomeric_AI':{'LOH':0,'TAI':1,'LST':1},
 # Whole-chromosome events are excluded by canonical scar definitions.
 'whole_chromosome_LOH':{'LOH':0,'TAI':0,'LST':0},
 # Two arm-internal transitions; the centromere gap must not add a third.
 'centromere_crossing_AI':{'LOH':0,'TAI':0,'LST':2},
 'note':'Hand-reasoned toy expectations; NOT outputs from pinned scarHRD. Canonical implementations may differ at boundaries/gaps. Validate separately before asserting exact equivalence.',
 # THE HEADLINE RESULT. Machine-checked at every run: balanced diploid (1,1)
 # and copy-neutral LOH (2,0) give identical log2 ratios at all purities.
 # If this is ever False, the forward model has been broken.
 'nonidentifiability':all(signal(1,1,p)[0]==signal(2,0,p)[0] for p in [.2,.5,.8,1.])}
 (a.out/'expected_and_limitations.json').write_text(json.dumps(expected,indent=2));print(len(segs),'segments;',len(obs),'probe observations; total-CN nonidentifiability verified')
if __name__=='__main__':main()
