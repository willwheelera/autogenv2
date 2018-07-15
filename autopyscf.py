from __future__ import print_function
import os
import shutil as sh
from pymatgen.io.cif import CifParser
import pyscf
from pyscf.scf.uhf import UHF,mulliken_meta
from copy import deepcopy
from autopaths import paths

##########################################################
# Tools for handling basis set input.
from xml.etree.ElementTree import ElementTree
def edit_xml_basis(xml_name,symbols,min_exp=0.2,naug=2,alpha=3,
                       cutoff=0.2,basis_name='vtz',
                       nangular={"s":1,"p":1,"d":1,"f":1,"g":0}
                       ):
  ''' Read a basis from an xml and edit according to the cutoff and augment recipe.

  Args:
    xml (str): path to xml to modify (not in place).
    symbols (list): elements to include in the basis.
    min_exp (float): smallest exponent for augmented basis. 
    naug (int): number of basis elements to add.
    alpha (float): each successive augment has exponent multiplied by this. 
    nangular (dict): number of BFD orbitals to keep. 
  Returns
    dict: Basis coefficients organized into a dict.
  '''

  allbasis={}

  for symbol in symbols:
    allbasis[symbol]=[]
    transition_metals=["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"]
    if symbol in transition_metals:
      nangular['s']=max(nangular['s'],2)
    tree = ElementTree()
    tree.parse(xml_name)
    element = tree.find('./Pseudopotential[@symbol="{}"]'.format(symbol))
    basis_path = './Basis-set[@name="{}"]/Contraction'.format(basis_name)

    # add in the first nangular basis functions.
    found_orbitals = []  
    for contraction in element.findall(basis_path):
      angular = contraction.get('Angular_momentum')
      if found_orbitals.count(angular) >= nangular[angular]:
        continue
      nterms = 0
      basis_sec={'angular':angular,'primitives':[]}
      for basis_term in contraction.findall('./Basis-term'):
        exp = basis_term.get('Exp')
        coeff = basis_term.get('Coeff')
        if float(exp) > cutoff:
          basis_sec['primitives'].append((exp, coeff))
          nterms+=1
      if nterms > 0:
        allbasis[symbol].append(basis_sec)
        found_orbitals.append(angular)

    angular_uncontracted=['s','p']
    if symbol in transition_metals:
      angular_uncontracted.append('d')

    for angular in angular_uncontracted:
      for i in range(0,naug):
        exp=min_exp*alpha**i
        basis_sec={'angular':angular,'primitives':[(exp,1.0)]}
        allbasis[symbol].append(basis_sec)
  return allbasis

####################################################
def format_basis(bdict):
  outlines=['{']
  for atom in sorted(bdict):
    outlines.append("    '%s':pyscf.gto.basis.parse('''"%atom)
    for belement in bdict[atom]:
      outlines.append("      %s %s"%(atom,belement['angular']))
      for exp,coef in belement['primitives']:
        outlines.append("        %s %s"%(exp,coef))
    outlines.append("    '''),")
  outlines.append('  }')
  return '\n'.join(outlines)

####################################################
class PySCFWriter:
  def __init__(self,options={}):
    self.basis='bfd_vtz'
    self.charge=0
    self.completed=False
    self.xc="" #Any valid input for PySCF. This gets put into the 'xc' variable
    self.diis_start_cycle=1
    self.ecp="bfd"
    self.level_shift=0.0
    self.conv_tol=1e-10
    self.max_cycle=50
    self.method='ROHF' 
    self.postHF=False   
    self.direct_scf_tol=1e-10
    self.spin=0
    self.xyz=""
    
    # ncore: %d     -- Number of core states.
    # ncas: %d      -- Number of states in active space. 
    # nelec: (%d,%d)-- Number of up (x) and down (y) electrons in active space.
    # tol: %g       -- tolerance on coefficient for det to be included in QWalk calculations.
    # method: %s    -- CASSCF or CASCI.
    self.cas={}

    # Used to name functions from conversion.

    # Default chosen by method at runtime.
    self.dm_generator=None

    self.set_options(options)
    
  #-----------------------------------------------
    
  def set_options(self, d):
    saved=deepcopy(self.__dict__)
    selfdict=self.__dict__

    # Check important keys are set. 
    for k in d.keys():
      if not k in selfdict.keys():
        raise AssertionError("Invalid option. %s not a keyword for PySCFWriter"%k)
      selfdict[k]=d[k]

    # If charge and spin should have same parity.
    if selfdict['charge']%2!=selfdict['spin']%2:
      reason='Spin and charge should both be even or both be odd. Charge=%d, spin=%d.'\
          %(selfdict['charge'],selfdict['spin'])
      selfdict['charge']=saved['charge']
      selfdict['spin']=saved['spin']
      raise AssertionError("Invalid option. "+reason)

    # If postHF got set, new options are required input.
    if self.postHF==True:
      for key in ['ncore','nelec','ncas','tol','method']:
        assert key in self.cas.keys(),"%s missing from 'cas' settings! "%key+\
            "Make sure all of 'ncore','nelec','ncas','tol','method' are set."

  #-----------------------------------------------
  def is_consistent(self,other):
    # dm_generator currently gets printed differently because of the dictionaries.
    skipkeys = ['completed','chkfile','dm_generator']
    for otherkey in other.__dict__.keys():
      if otherkey not in self.__dict__.keys():
        print(self.__class__.__name__,': other is missing a key.')
        return False
    for selfkey in self.__dict__.keys():
      if selfkey not in other.__dict__.keys():
        print(self.__class__.__name__,': self is missing a key.')
        return False
    for key in self.__dict__.keys():
      if self.__dict__[key]!=other.__dict__[key] and key not in skipkeys:
        print(self.__class__.__name__,": Different keys [{}] = \n{}\n or \n {}"\
            .format(key,self.__dict__[key],other.__dict__[key]))
        return False
    return True
    
  #-----------------------------------------------
  def pyscf_input(self,fname,chkfile):
    f=open(fname,'w')
    add_paths=[]

    # Figure out correct default initial guess (if not set).
    if self.dm_generator is None:
      if self.method in ['RKS','RHF','ROKS','ROHF']:
        self.dm_generator=dm_from_rhf_minao()
      elif self.method in ['UKS','UHF']:
        self.dm_generator=dm_from_uhf_minao()
      else:
        print(self.__class__.__name__,": Warning--default guess not set for method=%s.\n Trying UHF."%self.method)
        self.dm_generator=dm_from_uhf_minao()

    # Format the basis
    if type(self.basis)!=str:
      basisstr=format_basis(self.basis)
    else:
      basisstr="'%s'"%self.basis

    outlines=[
        "import sys",
        "sys.path.append('%s')"%paths['pyscf'],
        "import pyscf",
        "from pyscf import gto,scf,mcscf,fci,lib",
        "from pyscf.scf import RHF, ROHF, UHF",
        "from pyscf.dft.rks import RKS",
        "from pyscf.dft.roks import ROKS",
        "from pyscf.dft.uks import UKS",
        "basis=%s"%basis,
        "mol=gto.Mole(verbose=4)",
        "mol.build(atom='''"+self.xyz+"''',",
        "basis=basis,"%basisstr,
        "ecp='%s')"%self.ecp,
        "mol.charge=%i"%self.charge,
        "mol.spin=%i"%self.spin,
        "m=%s(mol)"%self.method,
        "m.max_cycle=%d"%self.max_cycle,
        "m.direct_scf_tol=%g"%self.direct_scf_tol,
        "m.chkfile='%s'"%chkfile,
        "m.conv_tol=%g"%self.conv_tol,
        "m.diis_start_cycle=%d"%self.diis_start_cycle
      ] + self.dm_generator

    if self.level_shift>0.0:
      outlines+=["m.level_shift=%g"%self.level_shift]
    
    if self.xc!="":
      outlines+=['m.xc="%s"'%self.xc]

    outlines+=["print('E(HF) =',m.kernel(init_dm))"]
    outlines+=['print ("HF_done")']        

    if self.postHF :
      outlines += ["mc=mcscf.%s(m, ncas=%i, nelecas=(%i, %i),ncore= %i)"%( 
                   self.cas['method'], self.cas['ncas'], self.cas['nelec'][0], 
                   self.cas['nelec'][1], self.cas['ncore']), 
                   "mc.direct_scf_tol=%f"%self.direct_scf_tol,
                   "mc.kernel()",
                   'print ("PostHF_done")']
    outlines += ['print ("All_done")']

    f.write('\n'.join(outlines))

    self.completed=True
     
####################################################
  
####################################################
class PySCFPBCWriter:
  def __init__(self,options={}):
    self.charge=0
    self.cif=''
    self.supercell=None
    self.completed=False
    self.xc="pbe,pbe" #Any valid input for PySCF. This gets put into the 'xc' variable
    self.diis_start_cycle=1
    self.cell_precision=1e-8 
    self.ecp="bfd"
    self.level_shift=0.0
    self.conv_tol=1e-7
    self.max_cycle=50
    self.method='RKS' 
    self.direct_scf_tol=1e-7
    self.spin=0
    self.gmesh=None
    self.xyz=""
    self.latticevec=""
    self.kpts=[2,2,2]
    self.basis='bfd_vtz'
    self.remove_linear_dep=False # float for setting tolerance.
    
    # Default chosen by method at runtime.
    self.dm_generator=None

    self.set_options(options)
  #-----------------------------------------------

  def from_cif(self,cifstring,primitive=True,supercell=[[1,0,0],[0,1,0],[0,0,1]]):
    struct=CifParser.from_string(cifstring).get_structures(primitive=primitive)[0]
    struct.make_supercell(supercell)
    struct=struct.as_dict()
    self.latticevec=" "
    for a in struct['lattice']['matrix']:
      for b in a:
        self.latticevec+= str(b)+ " "

    #print(struct['sites'])
    self.xyz=""
    self.elements=set()
    for s in struct['sites']:
      if len(s['species']) > 1:
        print("More than one species per site.. Taking the first.")
      self.xyz+=s['species'][0]['element']+" " + " ".join(map(str,s['xyz'])) + "\n"
      self.elements.add(s['species'][0]['element'])
    
  #-----------------------------------------------
    
  def set_options(self, d):
    selfdict=self.__dict__

    # Check important keys are set. 
    for k in d.keys():
      if not k in selfdict.keys():
        raise AssertionError("Error:",k,"not a keyword for PySCFWriter")
      selfdict[k]=d[k]

    # Must be done after bdf_library is set.
    if 'cif' in d.keys():
      if self.supercell is None:
        self.from_cif(d['cif'])
      else:
        self.from_cif(d['cif'],supercell=self.supercell)


  #-----------------------------------------------
  # Obsolete with update options?
  def is_consistent(self,other):
    issame,diff=self._check_keys(other,
        skipkeys = ['completed','chkfile','dm_generator'])
    if not issame:
      print("Inconsistency: {} from one doesn't match {} from other."\
          .format(diff['self'],diff['other']))
    return issame

  #-----------------------------------------------
      
  def pyscf_input(self,fname,chkfile):
    f=open(fname,'w')
    add_paths=[]

    if type(self.basis)!=str:
      basisstr=format_basis(self.basis)
    else:
      basisstr="%s"%self.basis

    # Figure out correct default initial guess (if not set).
    if self.dm_generator is None:
      if self.method in ['RKS','RHF','ROHF']:
        self.dm_generator=dm_from_rhf_minao()
      elif self.method in ['UKS','UHF']:
        self.dm_generator=dm_from_uhf_minao()
      else:
        print("Warning: default guess not set for method=%s.\n Trying UHF."%self.method)
        self.dm_generator=dm_from_uhf_minao()

    outlines=[
        "import sys",
        "sys.path.append('%s')"%paths['pyscf'],
        "import pyscf",
        "import numpy",
        "from pyscf.scf.addons import remove_linear_dep_",
        "from pyscf.pbc import gto,scf",
        "from pyscf.pbc.scf import KRHF as RHF",
        "from pyscf.pbc.scf import KUHF as UHF",
        "from pyscf.pbc.dft import KRKS as RKS",
        "from pyscf.pbc.dft import KUKS as UKS"
      ]

    if self.gmesh is not None:
      gmesh=["  gs="+str(self.gmesh)+","]
    else:
      gmesh=[]

    #The basis
    outlines+=["basis=%s"%basisstr]
    # The cell/molecule
    outlines+=[
        "mol=gto.M(verbose=4,",
        "  atom='''"+self.xyz+"''',",
        ]+gmesh+[
        "  a='''"+str(self.latticevec) +"''',",
        "  basis=basis,",
        "  precision=%s,"%self.cell_precision,
        "  spin=%i,"%self.spin,
        "  ecp='%s')"%self.ecp,
        "mol.charge=%i"%self.charge
        ]
    #Set up k-points
    outlines+=['kpts=mol.make_kpts('+str(self.kpts) + ')']
    
    #Mean field
    outlines+=[
        "m=%s(mol,kpts)"%self.method,
        "m.max_cycle=%d"%self.max_cycle,
        "m.direct_scf_tol=%g"%self.direct_scf_tol,
        "m.chkfile='%s'"%chkfile,
        "m.conv_tol=%g"%self.conv_tol,
        "m.diis_start_cycle=%d"%self.diis_start_cycle
      ] 
      
    outlines+=self.dm_generator
    if self.method in ['UKS','UHF']:
      outlines+=['dm_kpts= numpy.array([[init_dm[0] for k in range(len(kpts))],' +\
                          '[init_dm[1] for k in range(len(kpts))]])'] 
    else: 
      outlines+=['dm_kpts= [init_dm for k in range(len(kpts))]']

    if self.level_shift>0.0:
      outlines+=["m.level_shift=%g"%self.level_shift]
    
    if self.xc!="":
      outlines+=['m.xc="%s"'%self.xc]

    if self.remove_linear_dep:
      lincutoff=1e-8 if(self.remove_linear_dep is True) else self.remove_linear_dep
      outlines+=['m=remove_linear_dep_(m,lindep=%f)'%lincutoff]

    outlines+=["print('E(HF) =',m.kernel(numpy.array(dm_kpts)))"]
    outlines += ['print ("All_done")']
    f.write('\n'.join(outlines))

    self.completed=True
    
####################################################
class PySCFReader:
  def __init__(self):
    self.output={}
    self.completed=False

  #------------------------------------------------
  def read_chkfile(self,chkfile):
    ''' Read all data from the chkfile.'''
    ret={}
    mol=pyscf.lib.chkfile.load_mol(chkfile)

    # TODO density matrix for mcscf parts.
    # I don't think those results are saved in the chkfile.
    uhf=UHF(mol)
    dm=uhf.from_chk(chkfile)
    ret['basis_labels']=mol.spheric_labels(fmt=False)
    ret['density_matrix']=dm

    for key in ('scf','mcscf'):
      ret[key]=pyscf.lib.chkfile.load(chkfile,key)
    return ret
          
  ##------------------------------------------------
  # This restart check only works for MCSCF. I don't need that.
  #def check_restart(self, outfile):
  #  lines = open(outfile,'r').read().split('\n')
  #  if ('HF_done' in lines) and  ('All_done' not in lines):
  #    return True
  #  if 'SCF not converged.' in lines:
  #    return True
  #  return False

  #------------------------------------------------
  def check_restart(self,outfile):
    ''' Check if a restart is needed to complete. '''
    # Note: this only checks the restart an SCF run.
    restart=True
    for line in open(outfile,'r'):
      if "converged SCF energy" in line: restart=False
    return restart

  #------------------------------------------------
  def collect(self,outfile,chkfile):
    self.output={}
    self.output['file']=outfile

    converged=False
    lines = open(outfile,'r').read().split('\n')
    for line in reversed(lines):
      if "converged SCF energy" in line: 
        converged=True
        self.output.update(self.read_chkfile(chkfile))
        self.output['chkfile']=chkfile
        self.output['conversion']=[]
        break
    if converged:
      self.completed=True
      return 'done'
    else:
      return 'killed'

  #------------------------------------------------
  def write_summary(self):
    print("#### Variance optimization")
    for f,out in self.output.items():
      nruns=len(out)
      print(f,"Number of runs",nruns)
      for run in out:
        print("dispersion",run['sigma'])

# You should be also able to use the DM read into autogenv2, but I haven't
# thought about how precisely this should work.
####################################################
def dm_from_rhf_minao():
  return ["init_dm=pyscf.scf.rhf.init_guess_by_minao(mol)"]
      
####################################################
def dm_from_uhf_minao():
  return ["init_dm=pyscf.scf.uhf.init_guess_by_minao(mol)"]

####################################################
def dm_set_spins(atomspins,double_occ={},startdm=None):
  ''' startdm should be the location of a chkfile if not None. '''
  if startdm is None:
    dmstarter='pyscf.scf.uhf.init_guess_by_minao(mol)'
  else:
    dmstarter=startdm

  return [
    "atomspins=%r"%atomspins,
    "double_occ=%r"%double_occ,
    "init_dm=%s"%dmstarter,
    "print(init_dm[0].diagonal())",
    "for atmid, (shl0,shl1,ao0,ao1) in enumerate(mol.offset_nr_by_atom()):",
    "  if atmid < len(atomspins) and atomspins[atmid]!=0:",
    "    opp=int((atomspins[atmid]+1)/2)",
    "    s=(opp+1)%2",
    "    sym=mol.atom_pure_symbol(atmid)",
    "    print(sym,atmid,s,opp)",
    "    docc=[]",
    "    if sym in double_occ:",
    "      docc=double_occ[sym]",
    "    for ii,i in enumerate(range(ao0,ao1)):",
    "      if ii not in docc:",
    "        init_dm[opp][i,i]=0.0",
  ]

####################################################
def dm_from_chkfile(chkfile):
  """ Read a dm from a chkfile produced by a PySCF calculation. 

  It's preferrable to use absolute file paths, because the working directory
  will change to the folder where the driver is executed when this line is read."""
  return ["init_dm=m.from_chk('%s')"%chkfile]
