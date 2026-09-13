#!/usr/bin/env python3
"""Build a versioned, autosomal, technically masked HM450/EPIC-v1 probe bridge."""
import argparse,csv,gzip,hashlib,json,pathlib,urllib.request
ROOT=pathlib.Path(__file__).resolve().parents[1]
SOURCES={
 'HM450':('HM450.hg19.manifest.202209.tsv.gz','https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotationV1/main/Anno/HM450/archive/202209/HM450.hg19.manifest.tsv.gz'),
 'EPIC':('EPIC.hg19.manifest.202209.tsv.gz','https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotationV1/main/Anno/EPIC/archive/202209/EPIC.hg19.manifest.tsv.gz')}

def sha256(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()

def fetch(path,url):
 path.parent.mkdir(parents=True,exist_ok=True);part=path.with_suffix(path.suffix+'.part')
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'KIDS26-research/0.1'}),timeout=60) as r,part.open('wb') as f:
  for block in iter(lambda:r.read(1024*1024),b''):f.write(block)
 part.replace(path)

def load(path):
 with gzip.open(path,'rt',encoding='utf-8') as f:
  rows={r['probeID']:r for r in csv.DictReader(f,delimiter='\t') if r['probeID'].startswith('cg')}
 required={'probeID','CpG_chrm','CpG_beg','CpG_end','MASK_general'}
 if not rows or not required.issubset(next(iter(rows.values()))):raise ValueError(f'Unexpected annotation schema: {path}')
 return rows

def truthy(value):return value.strip().upper() in {'TRUE','T','1','YES'}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--download',action='store_true')
 p.add_argument('--annotation-dir',type=pathlib.Path,default=ROOT/'data/raw/annotation')
 p.add_argument('--out',type=pathlib.Path,default=ROOT/'config/shared_autosomal_probes.txt');a=p.parse_args()
 paths={k:a.annotation_dir/name for k,(name,url) in SOURCES.items()}
 for k,path in paths.items():
  if not path.exists():
   if not a.download:raise FileNotFoundError(f'{path} missing; rerun with --download')
   fetch(path,SOURCES[k][1])
 hm,ep=load(paths['HM450']),load(paths['EPIC']);common=set(hm)&set(ep);counts={'hm450_cg':len(hm),'epic_cg':len(ep),'common_id':len(common)};keep=[]
 autosomes={str(i) for i in range(1,23)}
 for probe in sorted(common):
  x,y=hm[probe],ep[probe];chrom=x['CpG_chrm'].removeprefix('chr')
  if chrom not in autosomes:continue
  if (x['CpG_chrm'],x['CpG_beg'],x['CpG_end'])!=(y['CpG_chrm'],y['CpG_beg'],y['CpG_end']):continue
  if truthy(x['MASK_general']) or truthy(y['MASK_general']):continue
  keep.append(probe)
 if not keep:raise ValueError('No probes survived bridge filters')
 a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text('\n'.join(keep)+'\n',encoding='ascii')
 meta={'bridge':'HM450_to_EPIC_v1_hg19','annotation_release':'Zhou InfiniumAnnotation archive/202209','filters':['probe ID present on both arrays','cg probes only','autosomes chr1-chr22','identical hg19 CpG coordinate','MASK_general false on both arrays'],'counts':{**counts,'retained':len(keep)},'output_sha256':sha256(a.out),'sources':{k:{'url':SOURCES[k][1],'path':str(path.relative_to(ROOT)),'bytes':path.stat().st_size,'sha256':sha256(path)} for k,path in paths.items()},'warning':'Valid for EPIC v1 only. Confirm PBTP array generation; EPIC v2 requires a separately versioned bridge.'}
 a.out.with_suffix('.metadata.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8');print(json.dumps(meta,indent=2))
if __name__=='__main__':main()
