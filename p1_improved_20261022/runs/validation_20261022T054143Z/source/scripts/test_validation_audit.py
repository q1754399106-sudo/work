"""Small synthetic controls for audit pairing, duplicates and uncertainty."""
import itertools,unittest
from analyze_validation import coverage,crossed_bootstrap,bad_solver

class AuditTests(unittest.TestCase):
 def setUp(self):
  self.jobs=[dict(model=dict(name=m),calibration_seeds=[1,2],test_seeds=[11,12],suites=[dict(name='Main',scenarios=['s1','s2'])]) for m in ['m1','m2']]
  self.plan=dict(jobs=self.jobs,methods=['Frozen','CombinedUncertaintyV2'],expected_rows=64)
  self.rows=[]
  for mi,m in enumerate(['m1','m2']):
   for c,o,s,method,route in itertools.product([1,2],[11,12],['s1','s2'],self.plan['methods'],['Direct','Split']):
    base=1+mi+.1*c+.01*o
    self.rows.append(dict(Model=m,Suite='Main',Scenario=s,CalibrationSeed=str(c),Seed=str(o),Method=method,Route=route,RelativeError=str(base*(.8 if method=='CombinedUncertaintyV2' else 1))))
 def test_exact_coverage_and_duplicate_replacement(self):
  self.assertTrue(coverage(self.rows,self.plan)['complete'])
  bad=self.rows[:-1]+[self.rows[0]];c=coverage(bad,self.plan)
  self.assertFalse(c['complete']);self.assertEqual(len(c['missing']),1);self.assertEqual(len(c['duplicates']),1)
 def test_paired_constant_gain_is_exact_under_every_axis_draw(self):
  result=crossed_bootstrap(self.rows,self.jobs,'Frozen',draws=80)
  for x in [result['relative_improvement'],*result['crossed_percentile_95_interval']]:self.assertAlmostEqual(x,.2)
 def test_all_scale_failure_or_nonfinite_cannot_pass(self):
  row=dict(FinalFlag='0',NormalResidual='1e-8',AllSelectionSolverFlags='0',AllSelectionMaxNormalResidual='1e-8')
  self.assertFalse(bad_solver(row,2e-7))
  for key,value in [('AllSelectionSolverFlags','1'),('AllSelectionMaxNormalResidual','3e-7'),('NormalResidual','nan'),('NormalResidual','-1')]:
   self.assertTrue(bad_solver(dict(row,**{key:value}),2e-7))
if __name__=='__main__':unittest.main()
