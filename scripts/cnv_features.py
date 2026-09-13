#!/usr/bin/env python3
"""Explicitly total-CN proxies, NOT canonical scarHRD. Coordinates: 0-based half-open."""
import argparse,csv,collections,math,pathlib,statistics

def features(segments,threshold=0.2,min_long=10_000_000,max_gap=3_000_000):
 if threshold<=0:raise ValueError('Positive amplitude threshold required')
 rows=[]
 for s in segments:
  r=dict(s);r['start']=int(r['start']);r['end']=int(r['end']);r['log2r']=float(r['log2r'])
  if r['start']<0 or r['end']<=r['start'] or not math.isfinite(r['log2r']):raise ValueError('Invalid segment')
  if r['arm'] not in ['p','q']:raise ValueError('Provide centromere-split p/q arms')
  if str(r['chromosome']) not in [str(c) for c in range(1,23)]:raise ValueError('Autosomes 1-22 required')
  rows.append(r)
 if not rows:raise ValueError('No segments')
 rows.sort(key=lambda r:(int(r['chromosome']),r['start']))
 for a,b in zip(rows,rows[1:]):
  if a['chromosome']==b['chromosome'] and a['end']>b['start']:raise ValueError('Overlapping segments')
 length=lambda r:r['end']-r['start']
 covered=sum(map(length,rows));altered=sum(length(r) for r in rows if abs(r['log2r'])>=threshold)
 transitions=0
 for a,b in zip(rows,rows[1:]):
  if a['chromosome']==b['chromosome'] and a['arm']==b['arm'] and length(a)>=min_long and length(b)>=min_long and b['start']-a['end']<=max_gap and abs(a['log2r']-b['log2r'])>=threshold:transitions+=1
 return {'covered_bp':covered,'segment_count':len(rows),'FGA_covered':altered/covered,
 'gain_fraction_covered':sum(length(r) for r in rows if r['log2r']>=threshold)/covered,
 'loss_fraction_covered':sum(length(r) for r in rows if r['log2r']<=-threshold)/covered,
 'amplitude_burden':sum(length(r)*abs(r['log2r']) for r in rows)/covered,
 'median_segment_bp':statistics.median(map(length,rows)),
 'large_altered_segments':sum(length(r)>=min_long and abs(r['log2r'])>=threshold for r in rows),
 'methyl_LST_like_conservative':transitions}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('segments',type=pathlib.Path);p.add_argument('output',type=pathlib.Path);p.add_argument('--threshold',type=float,default=.2);a=p.parse_args()
 groups=collections.defaultdict(list)
 with a.segments.open(encoding='utf-8-sig') as f:
  for r in csv.DictReader(f,delimiter='\t'):groups[r['sample_id']].append(r)
 out=[{'sample_id':s,**features(v,a.threshold)} for s,v in sorted(groups.items())]
 if not out:raise ValueError('Empty input')
 a.output.parent.mkdir(parents=True,exist_ok=True)
 with a.output.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t');w.writeheader();w.writerows(out)
if __name__=='__main__':main()
