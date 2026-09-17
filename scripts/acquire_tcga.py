#!/usr/bin/env python3
"""Public TCGA acquisition. Python >=3.10, standard library only. Run from repo root."""
import argparse,csv,hashlib,json,pathlib,re,time,urllib.request,urllib.parse,urllib.error
ROOT=pathlib.Path(__file__).resolve().parents[1]
API='https://api.gdc.cancer.gov'
TIERS={
 'labels':['TCGA.HRD_withSampleID.txt','TCGASubtype.20170308.tsv','merged_sample_quality_annotations.tsv','TCGA_mastercalls.abs_tables_JSedit.fixed.txt'],
 'beta':['jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv'],
 'cnv':['TCGA_mastercalls.abs_segtabs.fixed.txt','broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg'],
 'rna':['EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv']}
def digest(p,kind='md5'):
 h=hashlib.new(kind)
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def request(url,payload=None):
 req=urllib.request.Request(url,data=None if payload is None else json.dumps(payload).encode(),headers={'User-Agent':'KIDS26-research/0.1','Content-Type':'application/json'})
 for attempt in range(4):
  try:
   with urllib.request.urlopen(req,timeout=60) as r:return r.read()
  except (urllib.error.URLError,TimeoutError):
   if attempt==3:raise
   time.sleep(2**attempt)
def download(row,out):
 out.mkdir(parents=True,exist_ok=True)
 name=row['filename']
 if pathlib.Path(name).name!=name or '/' in name or '\\' in name:raise ValueError('Unsafe filename')
 dest=out/name;part=out/(name+'.part');size=int(row['size'])
 if dest.exists():
  if dest.stat().st_size==size and digest(dest)==row['md5']:return 'verified-existing'
  raise ValueError(f'Existing file failed verification: {dest}; move it aside after inspection')
 for attempt in range(4):
  try:
   offset=part.stat().st_size if part.exists() else 0
   if offset>size:raise ValueError('Partial file larger than manifest size')
   if offset<size:
    req=urllib.request.Request(API+'/data/'+row['id'],headers={**({'Range':f'bytes={offset}-'} if offset else {}),'Accept-Encoding':'identity','User-Agent':'KIDS26-research/0.1'})
    with urllib.request.urlopen(req,timeout=60) as res:
     append=offset>0 and res.status==206
     if res.status==206:
      cr=re.fullmatch(r'(?:bytes )?(\d+)-(\d+)/(\d+)',res.headers.get('Content-Range',''))
      if not cr or int(cr[1])!=offset or int(cr[3])!=size or int(cr[2])<offset:raise ValueError('Invalid range response')
     if 'text/html' in res.headers.get('Content-Type',''):raise ValueError('HTML instead of data')
     with part.open('ab' if append else 'wb') as f:
      for block in iter(lambda:res.read(1024*1024),b''):f.write(block)
   if part.stat().st_size!=size or digest(part)!=row['md5']:raise ValueError(f'Checksum/size failure: {part}. Inspect and remove only this .part before retrying.')
   part.replace(dest)
   # newline='\n' pins LF: write_text defaults to os.linesep translation, which emits CRLF on
   # Windows and breaks `sha256sum -c` by putting a stray \r inside the filename field.
   dest.with_suffix(dest.suffix+'.sha256').write_text(digest(dest,'sha256')+'  '+dest.name+'\n',newline='\n')
   return 'downloaded-verified'
  except (urllib.error.URLError,TimeoutError,OSError):
   if attempt==3:raise
   time.sleep(2**attempt)
def write_tsv(p,rows,fields):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',newline='',encoding='utf8') as f:
  # csv's default lineterminator is '\r\n'; pin LF so generated manifests are not CRLF.
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n');w.writeheader();w.writerows(rows)
def clean_row(r):
 """Strip stray whitespace/CR from manifest fields. Upstream source manifests are CRLF, and a
 trailing \\r in the last column corrupts any value taken from it (notably size, and filename
 whenever column order differs)."""
 return {k:(v or '').strip() for k,v in r.items()}
def published(tier,run,out):
 rows={}
 for name in ['panimmune','celloforigin']:
  with (ROOT/'config'/f'{name}_source_manifest.tsv').open(encoding='utf-8-sig') as f:
   for r in csv.DictReader(f,delimiter='\t'):r=clean_row(r);rows[r['filename']]=r
 selected=[rows[n] for n in TIERS[tier]]
 write_tsv(ROOT/'config'/f'tcga_{tier}.manifest.tsv',selected,['id','filename','md5','size'])
 print(json.dumps({'tier':tier,'files':len(selected),'bytes':sum(int(r['size']) for r in selected),'download':run},indent=2))
 for row in selected:
  print(row['filename'], download(row,out) if run else row['size'],flush=True)
def select_development(rows,limit,per_project=False):
 if per_project:
  groups={}
  for row in rows:groups.setdefault(row['cancer_type'],set()).add(row['patient_id'])
  chosen=set().union(*(set(sorted(ids)[:limit]) for ids in groups.values()))
 else:
  chosen=set();order=[]
  for row in rows:
   if row['patient_id'] not in chosen:chosen.add(row['patient_id']);order.append(row['patient_id'])
  chosen=set(order[:limit])
 return [row for row in rows if row['patient_id'] in chosen]

def discover(kind,projects,limit,run,out,per_project=False):
 # GDC methylation arrays are filed ONLY as 'Masked Intensities' (idat) or 'Methylation Beta Value'.
 # There is no 'Raw Intensities' methylation data_type: that value belongs to Affymetrix SNP6 and
 # GeneChip expression arrays, so querying it for a methylation platform always returns zero.
 dtype={'beta':'Methylation Beta Value','idat':'Masked Intensities','masked-idat':'Masked Intensities'}[kind]
 filters={'op':'and','content':[{'op':'in','content':{'field':'cases.project.project_id','value':projects}},{'op':'in','content':{'field':'data_type','value':[dtype]}},{'op':'in','content':{'field':'platform','value':['Illumina Human Methylation 450']}},{'op':'in','content':{'field':'access','value':['open']}}]}
 payload={'filters':filters,'format':'JSON','size':1000,'from':0,'sort':'file_id:asc','fields':'file_id,file_name,md5sum,file_size,access,data_type,data_format,platform,analysis.workflow_type,cases.submitter_id,cases.project.project_id,cases.samples.submitter_id,cases.samples.sample_type'}
 hits=[]
 while True:
  obj=json.loads(request(API+'/files',payload));hits+=obj['data']['hits']
  if len(hits)>=obj['data']['pagination']['total']:break
  payload['from']=len(hits)
 raw=out/'discovery';raw.mkdir(parents=True,exist_ok=True)
 (raw/f'{kind}_query.json').write_text(json.dumps(payload,indent=2))
 (raw/f'{kind}_response.json').write_text(json.dumps(hits,indent=2))
 rows=[]
 for h in hits:
  cases=h.get('cases',[])
  if len(cases)!=1:continue
  c=cases[0];samples=c.get('samples',[])
  if len(samples)!=1 or samples[0].get('sample_type')!='Primary Tumor':continue
  s=samples[0];patient=c['submitter_id']
  row={'id':h['file_id'],'filename':h['file_name'],'md5':h['md5sum'],'size':h['file_size'],'patient_id':patient,'sample_id':s['submitter_id'],'cancer_type':c['project']['project_id'].removeprefix('TCGA-'),'platform':h.get('platform'),'workflow':h.get('analysis',{}).get('workflow_type','unknown'),'access':h['access']}
  rows.append(row)
 # Keep all files for selected patients to preserve IDAT color pairs; beta duplicates are resolved by audited matching later.
 selected=select_development(rows,limit,per_project)
 fields=['id','filename','md5','size','patient_id','sample_id','cancer_type','platform','workflow','access']
 write_tsv(ROOT/'config'/f'gdc_{kind}_all.tsv',rows,fields)
 write_tsv(ROOT/'config'/f'gdc_{kind}_dev.tsv',selected,fields)
 write_tsv(ROOT/'config'/f'gdc_{kind}_dev.manifest.tsv',selected,['id','filename','md5','size'])
 print(f'{len(rows)} eligible files; {len(selected)} development files; {sum(int(r["size"]) for r in selected)} bytes')
 if not rows:print('No eligible open files. Check GDC facets/portal; do not interpret zero as proof raw assays do not exist.')
 if run:
  for r in selected:print(r['filename'],download(r,out/kind/r['id']),flush=True)
def main():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
 a=sub.add_parser('published');a.add_argument('--tier',choices=TIERS,default='labels');a.add_argument('--download',action='store_true');a.add_argument('--out',type=pathlib.Path,default=ROOT/'data/raw/pancanatlas')
 a=sub.add_parser('discover');a.add_argument('--kind',choices=['beta','idat','masked-idat'],default='beta');a.add_argument('--projects',nargs='+',default=['TCGA-BRCA','TCGA-OV']);a.add_argument('--limit',type=int,default=12,help='Total patients, or patients per project with --limit-per-project');a.add_argument('--limit-per-project',action='store_true');a.add_argument('--download',action='store_true');a.add_argument('--out',type=pathlib.Path,default=ROOT/'data/raw/gdc')
 a=p.parse_args()
 if a.command=='published':published(a.tier,a.download,a.out)
 else:
  if a.limit<1:p.error('--limit must be positive')
  discover(a.kind,a.projects,a.limit,a.download,a.out,a.limit_per_project)
if __name__=='__main__':main()
