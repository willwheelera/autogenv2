''' 
Runs a set of tests to check new updates are not breaking any features.
This is also a good demo of advanced features.
'''

from crystal import CrystalWriter
from crystalmanager import CrystalManager
from qwalkmanager import QWalkManager
from variance import VarianceWriter,VarianceReader
from linear import LinearWriter,LinearReader
from dmc import DMCWriter,DMCReader
from trialfunc import SlaterJastrow
from autorunner import RunnerPBS
import numpy as np

###################################################################################################################
# Individual tests.
def h2_crystal_equil_test():
  jobs=[]

  cwriter=CrystalWriter({
      'xml_name':'../../BFD_Library.xml',
      'cutoff':0.2,
      'spin_polarized':False
    })
  cwriter.set_struct_fromxyz(open('structures/h2_equil.xyz','r').read())

  cman=CrystalManager(
      name='crys',
      path='test_h2eq_crys',
      writer=cwriter,
      runner=RunnerPBS(
          queue='secondary',
          nn=1,np=16,
          walltime='0:01:00'
        )
    )
  jobs.append(cman)

  var=QWalkManager(
      name='var',
      path=cman.path,
      writer=VarianceWriter(),
      reader=VarianceReader(),
      runner=RunnerPBS(
          nn=1,np=16,queue='secondary',walltime='0:05:00'
        ),
      trialfunc=SlaterJastrow(cman,kpoint=0)
    )

  jobs.append(var)

  lin=QWalkManager(
      name='linear',
      path=cman.path,
      writer=LinearWriter(),
      reader=LinearReader(),
      runner=RunnerPBS(
          nn=1,np=16,queue='secondary',walltime='0:10:00'
        ),
      trialfunc=SlaterJastrow(slatman=cman,jastman=var,kpoint=0)
    )

  jobs.append(lin)

  for kidx in cman.qwfiles['kpoints']:
    kpt=cman.qwfiles['kpoints'][kidx]
    dmc=QWalkManager(
        name='dmc_%d'%kidx,
        path=cman.path,
        writer=DMCWriter(),
        reader=DMCReader(),
        runner=RunnerPBS(
            nn=1,np=16,queue='secondary',walltime='0:10:00',
          ),
        trialfunc=SlaterJastrow(slatman=cman,jastman=lin,kpoint=kidx)
      )
    jobs.append(dmc)

  return jobs

def h2_crystal_stretch_test():
  jobs=[]

  cwriter=CrystalWriter({
      'xml_name':'../../BFD_Library.xml',
      'cutoff':0.2,
      'spin_polarized':True,
      'initial_spins':[1,-1]
    })
  cwriter.set_struct_fromxyz(open('structures/h2_stretch.xyz','r').read())

  cman=CrystalManager(
      name='crys',
      path='test_h2st_crys',
      writer=cwriter,
      runner=RunnerPBS(
          queue='secondary',
          nn=1,np=16,
          walltime='0:01:00'
        )
    )
  jobs.append(cman)

  var=QWalkManager(
      name='var',
      path=cman.path,
      writer=VarianceWriter(),
      reader=VarianceReader(),
      runner=RunnerPBS(
          nn=1,np=16,queue='secondary',walltime='0:05:00'
        ),
      trialfunc=SlaterJastrow(cman,kpoint=0)
    )

  jobs.append(var)

  lin=QWalkManager(
      name='linear',
      path=cman.path,
      writer=LinearWriter(),
      reader=LinearReader(),
      runner=RunnerPBS(
          nn=1,np=16,queue='secondary',walltime='0:10:00'
        ),
      trialfunc=SlaterJastrow(slatman=cman,jastman=var,kpoint=0)
    )

  jobs.append(lin)

  for kidx in cman.qwfiles['kpoints']:
    kpt=cman.qwfiles['kpoints'][kidx]
    dmc=QWalkManager(
        name='dmc_%d'%kidx,
        path=cman.path,
        writer=DMCWriter(),
        reader=DMCReader(),
        runner=RunnerPBS(
            nn=1,np=16,queue='secondary',walltime='0:10:00',
          ),
        trialfunc=SlaterJastrow(slatman=cman,jastman=lin,kpoint=kidx)
      )
    jobs.append(dmc)

  return jobs

def si_crystal_test():
  ''' Simple tests that check PBC is working Crystal, and that QMC can be performed on the result.'''
  jobs=[]

  cwriter=CrystalWriter({
      'xml_name':'../../BFD_Library.xml',
      'cutoff':0.2,
      'kmesh':(3,3,3),
      'spin_polarized':False
    })
  cwriter.set_struct_fromcif(open('structures/si.cif','r').read(),primitive=True)

  cman=CrystalManager(
      name='crys',
      path='test_si_crys',
      writer=cwriter,
      runner=RunnerPBS(
          queue='secondary',
          nn=1,np=16,
          walltime='0:10:00'
        )
    )
  jobs.append(cman)

  var=QWalkManager(
      name='var',
      path=cman.path,
      writer=VarianceWriter(),
      reader=VarianceReader(),
      runner=RunnerPBS(
          nn=1,np=16,queue='secondary',walltime='0:10:00'
        ),
      trialfunc=SlaterJastrow(cman,kpoint=0)
    )

  jobs.append(var)

  lin=QWalkManager(
      name='linear',
      path=cman.path,
      writer=LinearWriter(),
      reader=LinearReader(),
      runner=RunnerPBS(
          nn=1,np=16,queue='secondary',walltime='1:00:00'
        ),
      trialfunc=SlaterJastrow(slatman=cman,jastman=var,kpoint=0)
    )

  jobs.append(lin)

  for kidx in cman.qwfiles['kpoints']:
    kpt=cman.qwfiles['kpoints'][kidx]
    dmc=QWalkManager(
        name='dmc_%d'%kidx,
        path=cman.path,
        writer=DMCWriter(),
        reader=DMCReader(),
        runner=RunnerPBS(
            nn=1,np=16,queue='secondary',walltime='1:00:00',
          ),
        trialfunc=SlaterJastrow(slatman=cman,jastman=lin,kpoint=kidx)
      )
    jobs.append(dmc)

  return jobs

def si_pyscf_test():
  ''' Simple tests that check PBC is working Crystal, and that QMC can be performed on the result.'''
  jobs=[]

  pwriter=PySCFPBCWriter({
      'gmesh':[16,16,16],
      'cif':open('si.cif','r').read()
    })
  pman=PySCFManager(
      name='scf',
      path='test_si_pyscf',
      writer=pwriter,
      runner=PySCFRunnerPBS(
          queue='secondary',
          nn=1,
          np=16,
          walltime='4:00:00'
        )
    )
  jobs.append(pman)

  var=QWalkManager(
      name='var',
      path=pman.path,
      writer=VarianceWriter(),
      reader=VarianceReader(),
      runner=RunnerPBS(
          nn=1,queue='secondary',walltime='0:10:00'
        ),
      trialfunc=SlaterJastrow(pman,kpoint=0)
    )
  jobs.append(var)

  lin=QWalkManager(
      name='lin',
      path=pman.path,
      writer=LinearWriter(),
      reader=LinearReader(),
      runner=RunnerPBS(
          nn=1,queue='secondary',walltime='0:30:00'
        ),
      trialfunc=SlaterJastrow(slatman=pman,jastman=var,kpoint=0)
    )
  jobs.append(lin)

  return jobs

def mno_test():
  ''' Test spinful calculations in PBC and real kpoint sampling. Also test trialfunc dependency. '''
  cwriter=CrystalWriter({
      'xml_name':'../../BFD_Library.xml',
      'cutoff':0.2,
      'kmesh':(3,3,3),
      'initial_spins':[1,0],
      'total_spin':5
    })
  cwriter.set_struct_fromcif(open('structures/mno.cif','r').read(),primitive=True)
  cwriter.set_options({})

  cman=CrystalManager(
      name='crys',
      path='test_mno_crys',
      writer=cwriter,
      runner=RunnerPBS(
          queue='secondary',
          nn=1,np=20,
          walltime='2:00:00'
        )
    )

  var=QWalkManager(
      name='var',
      path=cman.path,
      writer=VarianceWriter(),
      reader=VarianceReader(),
      runner=RunnerPBS(
          nn=1,np=20,queue='secondary',walltime='0:10:00'
        ),
      trialfunc=SlaterJastrow(cman,kpoint=0)
    )

  lin=QWalkManager(
      name='linear',
      path=cman.path,
      writer=LinearWriter(),
      reader=LinearReader(),
      runner=RunnerPBS(
          nn=1,np=20,queue='secondary',walltime='1:00:00'
        ),
      trialfunc=SlaterJastrow(slatman=cman,jastman=var,kpoint=0)
    )

  jobs=[cman,var,lin]
  for kidx in cman.qwfiles['kpoints']:
    kpt=cman.qwfiles['kpoints'][kidx]
    dmc=QWalkManager(
        name='dmc_%d'%kidx,
        path=cman.path,
        writer=DMCWriter(),
        reader=DMCReader(),
        runner=RunnerPBS(
            nn=1,np=16,queue='secondary',walltime='1:00:00',
          ),
        trialfunc=SlaterJastrow(slatman=cman,jastman=lin,kpoint=kidx)
      )
    jobs.append(dmc)
  return jobs

###################################################################################################################
# Test suites
def quick_crystal():
  ''' Fast-running tests of main crystal features. '''
  jobs=[]

  # Spinless molecule
  jobs+=h2_crystal_equil_test()

  # Spinless crystal
  jobs+=si_crystal_test()

  # Spinful molecule
  jobs+=h2_crystal_stretch_test()

  # Spinful crystal
  jobs+=mno_test()

  return jobs

def quick_pyscf():
  ''' Fast-running tests of main pyscf features. '''
  # Spinless molecule

  # Spinless crystal

  # Spinful molecule

  # Spinful crystal
  return [] #TODO

def thorough_crystal():
  ''' Longer tests that look over additional crystal features, compared to the quick tests.
  You should probably run the quick tests before running these.'''
  return [] # TODO

def thorough_pyscf():
  ''' Longer tests that look over additional pyscf features, compared to the quick tests.
  You should probably run the quick tests before running these.'''
  return [] # TODO

def run_tests():
  jobs=[]

  # Tests that are run.
  jobs+=quick_crystal()

  for job in jobs:
    job.nextstep()

if __name__=='__main__':
  run_tests()
