#!/usr/bin/env python3
"""Analytic identifiability fixtures plus noisy probe-grid observations. No clinical data."""
import argparse,csv,json,math,pathlib,random

def signal(a,b,purity):
 if a<0 or b<0 or not 0<=purity<=1:raise ValueError('Invalid allele count or purity')
 total=purity*(a+b)+2*(1-purity)
 if total<=0:return None,None
 return math.log2(total/2),(purity*b+(1-purity))/total
# Chromosome 1 toy length 100 Mb; centromere [45,55) Mb. A/B are allele-specific CN.
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
 a.out.mkdir(parents=True,exist_ok=True);rng=random.Random(a.seed);segs=[];obs=[]
 for name,intervals in FIXTURES.items():
  for start,end,major,minor in intervals:
   segs.append(dict(sample_id=name,chromosome='1',start=int(start*1e6),end=int(end*1e6),total_cn=major+minor,major_cn=major,minor_cn=minor))
  for purity in [.2,.5,.8,1.]:
   x=a.step_mb/2
   while x<100:
    # Probe-grid degradation and shifted observation relative to true breakpoints.
    lookup=x-a.breakpoint_shift_mb
    hit=next((v for v in intervals if v[0]<=lookup<v[1]),None)
    if hit:
     ratio,baf=signal(hit[2],hit[3],purity)
     obs.append(dict(sample_id=name,purity=purity,position=int(x*1e6),log2r=ratio+rng.gauss(0,a.noise),ideal_minor_fraction=baf))
    x+=a.step_mb
 for name,rows in [('allelic_truth.tsv',segs),('probe_observations.tsv',obs)]:
  with (a.out/name).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
 expected={'coordinate_system':'0-based half-open; toy chr1 length 100Mb; centromere 45-55Mb',
 'balanced':{'LOH':0,'TAI':0,'LST':0},'copy_neutral_LOH':{'LOH':1,'TAI':1,'LST_allelic_boundary':1},
 'balanced_gain':{'LOH':0,'TAI':0,'LST':1},'telomeric_AI':{'LOH':0,'TAI':1,'LST':1},
 'whole_chromosome_LOH':{'LOH':0,'TAI':0,'LST':0},
 'centromere_crossing_AI':{'LOH':0,'TAI':0,'LST':2},
 'note':'Hand-reasoned toy expectations; NOT outputs from pinned scarHRD. Canonical implementations may differ at boundaries/gaps. Validate separately before asserting exact equivalence.',
 'nonidentifiability':all(signal(1,1,p)[0]==signal(2,0,p)[0] for p in [.2,.5,.8,1.])}
 (a.out/'expected_and_limitations.json').write_text(json.dumps(expected,indent=2));print(len(segs),'segments;',len(obs),'probe observations; total-CN nonidentifiability verified')
if __name__=='__main__':main()
