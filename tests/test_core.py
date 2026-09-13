import pathlib,sys,unittest,tempfile,hashlib
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'scripts'))
from build_master import barcode,build
from cnv_features import features
from simulate_scars import signal
from acquire_tcga import download,select_development
class IdentityTests(unittest.TestCase):
 def inputs(self):
  m=[dict(patient_id='TCGA-AB-1234',sample_id='TCGA-AB-1234-01A',cancer_type='BRCA',id='file',filename='x',md5='')]
  h=[{'sampleID':'TCGA-AB-1234-01','ai1':'2','lst1':'3','hrd-loh':'4','HRD':'9'}]
  t=[{'pan.samplesID':'TCGA-AB-1234','cancer.type':'BRCA'}]
  q=[{'platform':'HumanMethylation450','aliquot_barcode':'TCGA-AB-1234-01A-01D-0000-01','Do_not_use':'False','AWG_excluded_because_of_pathology':'0.0'}]
  return m,h,t,q,[]
 def test_good_join_and_component_order(self):
  out,audit=build(*self.inputs());self.assertEqual((out[0]['HRD_LOH'],out[0]['LST'],out[0]['TAI']),(4,3,2))
 def test_no_patient_only_join(self):
  x=self.inputs();x[1][0]['sampleID']='TCGA-AB-1234-02';self.assertFalse(build(*x)[0])
 def test_replicate_excluded(self):
  x=self.inputs();x[0].append(dict(x[0][0]));self.assertFalse(build(*x)[0])
 def test_nonprimary_array_does_not_exclude_unique_primary(self):
  x=self.inputs();x[0].append({**x[0][0],'sample_id':'TCGA-AB-1234-11A'})
  out,audit=build(*x);self.assertEqual(len(out),1)
  self.assertEqual([r['reason'] for r in audit if r['sample_id'].endswith('11A')],['not_primary_tumor'])
 def test_bad_sum_excluded(self):
  x=self.inputs();x[1][0]['HRD']='12';self.assertFalse(build(*x)[0])
 def test_authoritative_metadata_type_allowed_when_subtype_row_absent(self):
  x=self.inputs();x[2].clear();self.assertEqual(len(build(*x)[0]),1)
 def test_conflicting_tumor_type_excluded(self):
  x=self.inputs();x[0][0]['cancer_type']='OV';self.assertFalse(build(*x)[0])
 def test_missing_qc_excluded(self):
  x=self.inputs();x[3].clear();self.assertFalse(build(*x)[0])
 def test_quality_exclusion(self):
  x=self.inputs();x[3][0]['Do_not_use']='True';self.assertFalse(build(*x)[0])
 def test_invalid_barcode(self):
  with self.assertRaises(ValueError):barcode('not-a-patient')
class SignalTests(unittest.TestCase):
 def test_total_cn_cannot_identify_loh(self):
  for p in [.2,.5,1]:
   self.assertEqual(signal(1,1,p)[0],signal(2,0,p)[0]);self.assertNotEqual(signal(1,1,p)[1],signal(2,0,p)[1])
 def test_length_weighted_features(self):
  s=[dict(chromosome='1',arm='p',start=0,end=10_000_000,log2r=.4),dict(chromosome='1',arm='p',start=10_000_000,end=40_000_000,log2r=0)]
  x=features(s);self.assertEqual(x['FGA_covered'],.25);self.assertEqual(x['methyl_LST_like_conservative'],1)
 def test_no_centromere_transition(self):
  s=[dict(chromosome='1',arm='p',start=0,end=20_000_000,log2r=.4),dict(chromosome='1',arm='q',start=22_000_000,end=42_000_000,log2r=0)]
  self.assertEqual(features(s)['methyl_LST_like_conservative'],0)
 def test_overlap_rejected(self):
  s=[dict(chromosome='1',arm='p',start=0,end=20,log2r=0),dict(chromosome='1',arm='p',start=10,end=30,log2r=.4)]
  with self.assertRaises(ValueError):features(s)
class DownloadTests(unittest.TestCase):
 def test_balanced_development_selection_preserves_patient_files(self):
  rows=[{'patient_id':p,'cancer_type':c,'id':str(i)} for i,(p,c) in enumerate([('a','A'),('a','A'),('b','A'),('c','B'),('c','B'),('d','B')])]
  out=select_development(rows,1,True);self.assertEqual([(x['patient_id'],x['id']) for x in out],[('a','0'),('a','1'),('c','3'),('c','4')])
 def test_verified_existing_and_corruption(self):
  with tempfile.TemporaryDirectory() as d:
   p=pathlib.Path(d);(p/'x').write_bytes(b'abc');row=dict(filename='x',size='3',md5=hashlib.md5(b'abc').hexdigest(),id='unused')
   self.assertEqual(download(row,p),'verified-existing');(p/'x').write_bytes(b'bad')
   with self.assertRaises(ValueError):download(row,p)
 def test_path_traversal_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   with self.assertRaises(ValueError):download(dict(filename='../bad',size='0'),pathlib.Path(d))
if __name__=='__main__':unittest.main()
