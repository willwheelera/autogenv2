''' 
Routines for extracting crystal parameters and results and writing them into qwalk input files.
A casual user would be interested in convert_crystal. 
'''

from __future__ import division,print_function
import numpy as np
import sys

def error(message,errortype):
  print(message)
  exit(errortype)

periodic_table = [
  "h","he","li","be","b","c","n","o","f","ne","na","mg","al","si","p","s","cl","ar",
  "k","ca","sc","ti","v","cr","mn","fe","co","ni","cu","zn","ga","ge","as","se","br",
  "kr","rb","sr","y","zr","nb","mo","tc","ru","rh","pd","ag","cd","in","sn","sb","te",
  "i","xe","cs","ba","la","ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm",
  "yb","lu","hf","ta","w","re","os","ir","pt","au","hg","tl","pb","bi","po","at","rn",
  "fr","ra","ac","th","pa","u","np","pu","am","cm","bk","cf","es","fm","md","no","lr",
  "rf","db","sg","bh","hs","mt","ds","rg","cp","uut","uuq","uup","uuh","uus","uuo"
]

###############################################################################
# Convenince method
def convert_crystal(
    base="qwalk",
    propoutfn="prop.in.o",
    realonly=False,
    nvirtual=50):
  """
  Uses rountines in this library to convert crystal files into qwalk files in one call.
  Files are named by [base]_[kindex].sys etc.
  Args:
    base (str): see naming.
    propoutfn (str): name of either crystal or properties output file.
    realonly (bool): whether to only the real kpoints.
    nvirtual (int): number of virtual orbtials to include in orbitals section.
  Returns:
    dict: files produced by this call.
  """
  # kfmt='coord' is probably a bad thing because it doesn't always work and can 
  # lead to unexpected changes in file name conventions.

  # keeps track of the files that get produced.
  files={}

  info, lat_parm, ions, basis, pseudo = read_gred()
  eigsys = read_kred(info,basis)

  if eigsys['nspin'] > 1:
    eigsys['totspin'] = read_outputfile(propoutfn)
  else:
    eigsys['totspin'] = 0

  # Useful quantities.
  basis['ntot'] = int(round(sum(basis['charges'])))
  basis['nmo']  = sum(basis['nao_shell']) # = nao
  eigsys['nup'] = int(round(0.5 * (basis['ntot'] + eigsys['totspin'])))
  eigsys['ndn'] = int(round(0.5 * (basis['ntot'] - eigsys['totspin'])))

  maxmo_spin=min(max(eigsys['nup'],eigsys['ndn'])+nvirtual,basis['nmo'])

  #  All the files that will get produced.
  files={
      'kpoints':{},
      'basis':base+".basis",
      'jastrow2':base+".jast2",
      'orbplot':{},
      'orb':{},
      'sys':{},
      'slater':{}
    }
  write_basis(basis,ions,files['basis'])
  write_jast2(lat_parm,ions,files['jastrow2'])
 
  for kpt in eigsys['kpt_coords']:
    if eigsys['ikpt_iscmpx'][kpt] and realonly: continue
    kidx=eigsys['kpt_index'][kpt]
    files['kpoints'][kidx]=kpt
    files['orbplot'][kidx]="%s_%d.plot"%(base,kidx)
    files['slater'][kidx]="%s_%d.slater"%(base,kidx)
    files['orb'][kidx]="%s_%d.orb"%(base,kidx)
    files['sys'][kidx]="%s_%d.sys"%(base,kidx)
    write_slater(basis,eigsys,kpt,
        outfn=files['slater'][kidx],
        orbfn=files['orb'][kidx],
        basisfn=files['basis'],
        maxmo_spin=maxmo_spin)
    write_orbplot(basis,eigsys,kpt,
        outfn=files['orbplot'][kidx],
        orbfn=files['orb'][kidx],
        basisfn=files['basis'],
        sysfn=files['sys'][kidx],
        maxmo_spin=maxmo_spin)
    write_orb(eigsys,basis,ions,kpt,files['orb'][kidx],maxmo_spin)
    write_sys(lat_parm,basis,eigsys,pseudo,ions,kpt,files['sys'][kidx])

  return files

###############################################################################
def read_gred(gred="GRED.DAT"):
  ''' Read the structure, basis, and pseudopotential from the GRED.DAT file.
  Args:
    gred (str): path to GRED.DAT file.
  Returns:
    tuple: (info,lat_parm,ions,basis,pseudo). Each are dictionaries with the following infomation:
      info (information useful for KRED.DAT),
      lat_parm (lattice parameters), 
      ions (ion positions and charges), 
      basis (basis set definition), 
      pseudo (pseudopotential defintion).
  '''
  lat_parm = {}
  ions = {}
  basis = {}
  pseudo = {}

  gred = open(gred,'r').read()

  # Fix numbers with no space between them.
  gred = gred.replace("-"," -")
  gred = gred.replace("E -","E-") 

  gred_words = gred.split()
  nparms = [int(w) for w in gred_words[1:4]]
  cursor = 4

  # These follow naming of cryapi_inp (but "inf" -> "info").
  info = [int(w) for w in gred_words[cursor          :cursor+nparms[0]]]
  itol = [int(w) for w in gred_words[cursor+nparms[0]:cursor+nparms[1]]]
  par  = [int(w) for w in gred_words[cursor+nparms[1]:cursor+nparms[2]]]
  cursor += sum(nparms)

  lat_parm['struct_dim'] = int(info[9])

  # Lattice parameters.
  lat_parm['latvecs'] = \
      np.array(gred_words[cursor:cursor+9],dtype=float).reshape(3,3).T.round(15)
  if (lat_parm['latvecs'] > 100).any():
    print("Lattice parameter larger than 100 A! Reducing to 100.")
    print("If this is a dimension < 3 system, there is no cause for alarm.")
    print("Otherwise if this is a problem for you, please generalize crystal2qmc.")
    lat_parm['latvecs'][lat_parm['latvecs']>100] = 100.
  cursor += 9
  prim_trans= np.array(gred_words[cursor:cursor+9],dtype=float).reshape(3,3)
  cursor += 9
  lat_parm['conv_cell'] = prim_trans.dot(lat_parm['latvecs'])
  cursor += info[1] + 48*48 + 9*info[1] + 3*info[1] # Skip symmetry part.

  # Lattice "stars" (?) skipped.
  cursor += info[4]+1 + info[78]*3 + info[4]+1 + info[4]+1 + info[78] + info[78]*3

  # Some of ion information.
  natoms = info[23]
  ions['charges'] = [float(w) for w in gred_words[cursor:cursor+natoms]]
  cursor += natoms
  # Atom positions.
  atom_poss = np.array(gred_words[cursor:cursor+3*natoms],dtype=float)
  ions['positions'] = atom_poss.reshape(natoms,3)
  cursor += 3*natoms

  # Basis information (some ion information mixed in).
  nshells = info[19]
  nprim   = info[74]
  # Formal charge of shell.
  basis['charges'] = np.array(gred_words[cursor:cursor+nshells],dtype=float)
  cursor += nshells
  # "Adjoined gaussian" of shells.
  basis['adj_gaus'] = np.array(gred_words[cursor:cursor+nshells],dtype=float)
  cursor += nshells
  # Position of shell.
  shell_poss = np.array(gred_words[cursor:cursor+3*nshells],dtype=float)
  basis['positions'] = shell_poss.reshape(nshells,3)
  cursor += 3*nshells
  # Primitive gaussian exponents.
  basis['prim_gaus'] = np.array(gred_words[cursor:cursor+nprim],dtype=float)
  cursor += nprim
  # Coefficients of s, p, d, and (?).
  basis['coef_s'] = np.array(gred_words[cursor:cursor+nprim],dtype=float)
  cursor += nprim
  basis['coef_p'] = np.array(gred_words[cursor:cursor+nprim],dtype=float)
  cursor += nprim
  basis['coef_dfg'] = np.array(gred_words[cursor:cursor+nprim],dtype=float)
  cursor += nprim
  basis['coef_max'] = np.array(gred_words[cursor:cursor+nprim],dtype=float)
  cursor += nprim
  # Skip "old normalization"
  cursor += 2*nprim
  # Atomic numbers.
  ions['atom_nums'] = np.array(gred_words[cursor:cursor+natoms],dtype=int)
  cursor += natoms
  # First shell of each atom (skip extra number after).
  basis['first_shell'] = np.array(gred_words[cursor:cursor+natoms],dtype=int)
  cursor += natoms + 1
  # First primitive of each shell (skips an extra number after).
  basis['first_prim'] = np.array(gred_words[cursor:cursor+nshells],dtype=int)
  cursor += nshells + 1
  # Number of prims per shell.
  basis['prim_shell'] = np.array(gred_words[cursor:cursor+nshells],dtype=int)
  cursor += nshells
  # Type of shell: 0=s,1=sp,2=p,3=d,4=f.
  basis['shell_type'] = np.array(gred_words[cursor:cursor+nshells],dtype=int)
  cursor += nshells
  # Number of atomic orbtials per shell.
  basis['nao_shell'] = np.array(gred_words[cursor:cursor+nshells],dtype=int)
  cursor += nshells
  # First atomic orbtial per shell (skip extra number after).
  basis['first_ao'] = np.array(gred_words[cursor:cursor+nshells],dtype=int)
  cursor += nshells + 1
  # Atom to which each shell belongs.
  basis['atom_shell'] = np.array(gred_words[cursor:cursor+nshells],dtype=int)
  cursor += nshells

  # Pseudopotential information.
  # Pseudopotential for each element.
  pseudo_atom = np.array(gred_words[cursor:cursor+natoms],dtype=int)
  cursor += natoms
  cursor += 1 # skip INFPOT
  ngauss = int(gred_words[cursor])
  cursor += 1
  headlen = int(gred_words[cursor])
  cursor += 1
  # Number of pseudopotentials.
  numpseudo = int(gred_words[cursor])
  cursor += 1
  # Exponents of r^l prefactor.
  r_exps = -1*np.array(gred_words[cursor:cursor+ngauss],dtype=int)
  cursor += ngauss
  # Number of Gaussians for angular momenutum j
  n_per_j = np.array(gred_words[cursor:cursor+headlen],dtype=int)
  cursor += headlen
  # index of first n_per_j for each pseudo.
  pseudo_start = np.array(gred_words[cursor:cursor+numpseudo],dtype=int)
  cursor += numpseudo + 1
  # Actual floats of pseudopotential.
  exponents = np.array(gred_words[cursor:cursor+ngauss],dtype=float)
  cursor += ngauss
  prefactors = np.array(gred_words[cursor:cursor+ngauss],dtype=float)
  cursor += ngauss
  # Store information nicely.
  npjlen = int(headlen / len(pseudo_start))
  for aidx,atom in enumerate(ions['atom_nums']):
    psidx = pseudo_atom[aidx]-1
    start = pseudo_start[psidx]
    if psidx+1 >= len(pseudo_start): end = ngauss
    else                           : end = pseudo_start[psidx+1]
    if atom not in pseudo.keys():
      pseudo[atom] = {}
      pseudo[atom]['prefactors'] = prefactors[start:end]
      pseudo[atom]['r_exps'] = r_exps[start:end]
      pseudo[atom]['n_per_j'] = n_per_j[npjlen*psidx:npjlen*(psidx+1)]
      pseudo[atom]['exponents'] = exponents[start:end]

  ## Density matrix information.
  # This is impossible to figure out.  See `cryapi_inp.f`.
  #atomic_charges = np.array(gred_words[cursor:cursor+natoms],dtype=float)
  #cursor += natoms
  #mvlaf = info[55] #???
  ## Skip symmetry information.
  #cursor += mvlaf*4 + info[19]*info[1] + 
  #print("atomic_charges",atomic_charges)

  return info, lat_parm, ions, basis, pseudo

###############################################################################
# Reads in kpoints and eigen{values,vectors} from KRED.DAT.
def read_kred(info,basis,kred="KRED.DAT"):
  ''' Read the KRED and provide information about the CRYSTAL solutions. 
  Args:
    info (dict): should be produced by read_gred.
    basis (dict): also from read_gred, the basis set info.
    kred (str): path to KRED.DAT.
  Returns:
    eigsys (dict): orbitals from the SCF calculation. 
  '''

  charcount=0
  eigsys = {
      'nkpts_dir':None,
      'recip_vecs':None,
      'kpt_index':None,
      'ikpt_iscmpx':None,
      'kpt_weights':None,
      'nspin':None,
      'eigvals':None,
      'eig_weights':None,
      'kpt_file_start':{},
      'kred':kred
    }

  kred = open(kred)
  kred_words = []
  for lin in kred:
    charcount+=len(lin)
    # Stop at eigenvectors. The second condition is to avoid stopping early for gamma-only calculations.
    if lin=='          0          0          0\n' and len(kred_words)>13:
      break # We'll do the eigenvectors one at a time to save memory.
    kred_words += lin.split()
  eigsys['kpt_file_start'][(0,0,0)]=[charcount]

  cursor = 0

  # Number of k-points in each direction.
  eigsys['nkpts_dir'] = np.array([int(w) for w in kred_words[cursor:cursor+3]])
  cursor += 3
  # Total number of inequivilent k-points.
  nikpts = int(kred_words[cursor])
  cursor += 1
  # Reciprocal basis.
  recip_vecs = np.array(kred_words[cursor:cursor+9],dtype=float)
  eigsys['recip_vecs'] = recip_vecs.reshape(3,3)
  cursor += 9
  # Inequivilent k-point coord in reciprocal basis.
  ikpt_coords = np.array(kred_words[cursor:cursor+3*nikpts],int)
  ikpt_coords = list(map(tuple,ikpt_coords.reshape(nikpts,3)))
  # Useful to compare to old output format.
  eigsys['kpt_index'] = dict(zip(ikpt_coords,range(len(ikpt_coords))))
  cursor += 3*nikpts
  # is complex (0) or not (1), converted to True (if complex) or False
  ikpt_iscmpx = \
    np.array([int(w) for w in kred_words[cursor:cursor+nikpts]]) == 0
  eigsys['ikpt_iscmpx'] = dict(zip(ikpt_coords,ikpt_iscmpx))
  cursor += nikpts
  # Skip symmetry information.
  cursor += 9*48
  # Geometric weight of kpoints.
  eigsys['kpt_weights'] = np.array(kred_words[cursor:cursor+nikpts],dtype=float)
  cursor += nikpts
  # Eigenvalues: (how many) = (spin) * (number of basis) * (number of kpoints)
  eigsys['nspin'] = info[63]+1
  nevals = eigsys['nspin']*info[6]*nikpts
  eigsys['eigvals'] = np.array(kred_words[cursor:cursor+nevals],dtype=float)
  cursor += nevals
  # Weights of eigenvales--incorperating Fermi energy cutoff.
  nbands = int(round(nevals / nikpts / eigsys['nspin']))
  eigsys['eig_weights'] = np.array(kred_words[cursor:cursor+nevals],dtype=float)\
      .reshape(nikpts,eigsys['nspin'],nbands)
  cursor += nevals

  # Information about eigenvectors.
  #nkpts  = np.prod(eigsys['nkpts_dir'])
  eigsys['nao'] = sum(basis['nao_shell'])
  eigsys['nbands'] = int(round(nevals / nikpts / eigsys['nspin']))
  
  # Here we simply mark where the eigenvectors are for later lookup.
  for line in kred:
    llen=len(line)
    charcount+=llen
    if llen==34:
      kpt=tuple([int(i) for i in line.split()])
      if kpt in eigsys['kpt_file_start']:
        eigsys['kpt_file_start'][kpt].append(charcount)
      else:
        eigsys['kpt_file_start'][kpt] = [charcount]

  ## It's probably true that kpt_coords == ikpt_coords, with repitition for spin
  ## up and spin down, because we only read in inequivilent kpoints. However,
  ## ordering might be different, and the ordering is correct for kpt_coords.
  ## If there are bugs, this might be a source.
  eigsys['kpt_coords'] = ikpt_coords # kpt_coords
  #eigsys['eigvecs'] = eigvecs

  return eigsys

###############################################################################
# Look up an eigenvector from KRED.DAT.
# TODO: Further reduction in memory usage can be had by specifying nvirtual here.
def eigvec_lookup(kpt,eigsys,spin=0):
  ''' Look up eigenvector at kpt from KRED.DAT using information from eigsys about where the eigenvectors start and end.
  Args:
    kpt (tuple of int): Kpoint coordinates.
    eigsys (dict): data from read_kred. 
    iscomplex (bool): Is the kpoint complex.
    spin (int): desired spin component. 
  Returns:
    array: eigenstate indexed by [band, ao]
  '''
  ncpnts = int(eigsys['nbands']* eigsys['nao'])
  if eigsys['ikpt_iscmpx'][kpt]:
    ncpnts *= 2
  linesperkpt = ncpnts//4 + int(ncpnts%4>0)

  kredf = open(eigsys['kred'],'r')
  kredf.seek(eigsys['kpt_file_start'][kpt][spin])
  eigvec=[]
  for li,line in enumerate(kredf):
    if li==linesperkpt: break
    eigvec+=line.split()

  eigvec = np.array(eigvec,dtype=float)

  if eigsys['ikpt_iscmpx'][kpt]:
    eigvec=eigvec.reshape(ncpnts//2,2)
    eigvec=eigvec[:,0] + eigvec[:,1]*1j # complexify.

  return eigvec.reshape(eigsys['nbands'],eigsys['nao'])

###############################################################################
def read_outputfile(fname = "prop.in.o"):
  ''' Reads total spin from output file. 
  Args: 
    fname (str): either crystal or properties output.
  Returns:
    int: spin of the system.
  '''
  fin = open(fname,'r')
  for line in fin:
    if "SUMMED SPIN DENSITY" in line:
      spin = float(line.split()[-1])
  if abs(round(spin) - spin) > 1e-8:
    print("Warning: spin %f is not close to integer!"%spin)
    print("  I'm rounding this to %d."%int(round(spin)))
  spin = int(round(spin))
  return spin

###############################################################################
def find_basis_cutoff(lat_parm):
  if lat_parm['struct_dim'] > 0:
    latvecs = lat_parm['latvecs']
    cutoff_divider = 2.000001
    cross01 = np.cross(latvecs[0], latvecs[1])
    cross12 = np.cross(latvecs[1], latvecs[2])
    cross02 = np.cross(latvecs[0], latvecs[2])

    heights = [0,0,0]
    heights[0]=abs(np.dot(latvecs[0], cross12)/np.dot(cross12,cross12)**.5)
    heights[1]=abs(np.dot(latvecs[1], cross02)/np.dot(cross02,cross02)**.5)
    heights[2]=abs(np.dot(latvecs[2], cross01)/np.dot(cross01,cross01)**.5)
    return min(heights)/cutoff_divider
  else:
    return 7.5

###############################################################################
def write_slater(basis,eigsys,kpt,outfn,orbfn,basisfn,maxmo_spin=-1):
  ntot = basis['ntot']
  nmo  = basis['nmo']
  nup  = eigsys['nup']
  ndn  = eigsys['ndn']
  if maxmo_spin < 0:
    maxmo_spin=nmo
  uporbs = np.arange(nup)+1
  dnorbs = np.arange(ndn)+1
  if eigsys['nspin'] > 1:
    dnorbs += maxmo_spin
  if eigsys['ikpt_iscmpx'][kpt]: orbstr = "corbitals"
  else:                          orbstr = "orbitals"
  uporblines = ["{:5d}".format(orb) for orb in uporbs]
  width = 10
  for i in reversed(range(width,len(uporblines),width)):
    uporblines.insert(i,"\n ")
  dnorblines = ["{:5d}".format(orb) for orb in dnorbs]
  for i in reversed(range(width,len(dnorblines),width)):
    dnorblines.insert(i,"\n ")
  outlines = [
      "slater",
      "{0} {{".format(orbstr),
      "  magnify 1",
      "  nmo {0}".format(dnorbs[-1]),
      "  orbfile {0}".format(orbfn),
      "  include {0}".format(basisfn),
      "  centers { useglobal }",
      "}",
      "detwt { 1.0 }",
      "states {",
      "  # Spin up orbitals.", 
      "  " + " ".join(uporblines),
      "  # Spin down orbitals.",
      "  " + " ".join(dnorblines),
      "}"
    ]
  with open(outfn,'w') as outf:
    outf.write("\n".join(outlines))

###############################################################################
def write_orbplot(basis,eigsys,kpt,outfn,orbfn,basisfn,sysfn,maxmo_spin=-1):
  ntot = basis['ntot']
  nmo  = basis['nmo']
  nup  = eigsys['nup']
  ndn  = eigsys['ndn']
  uporbs = np.arange(nup)+1
  dnorbs = np.arange(ndn)+1
  if maxmo_spin <0 and eigsys['nspin'] > 1:
    dnorbs += nmo
  else:
    dnorbs += maxmo_spin
  if eigsys['ikpt_iscmpx'][kpt]: orbstr = "corbitals"
  else:                          orbstr = "orbitals"
  uporblines = ["{:5d}".format(orb) for orb in uporbs]
  width = 10
  for i in reversed(range(width,len(uporblines),width)):
    uporblines.insert(i,"\n ")
  dnorblines = ["{:5d}".format(orb) for orb in dnorbs]
  for i in reversed(range(width,len(dnorblines),width)):
    dnorblines.insert(i,"\n ")
  outlines_prefix = [
      "method { ",
      "plot",
      "{0} {{".format(orbstr),
      "  magnify 1",
      "  nmo {0}".format(dnorbs[-1]),
      "  orbfile {0}".format(orbfn),
      "  include {0}".format(basisfn),
      "  centers { useglobal }",
      "}",
      "plotorbitals {",
    ]
  outlines_postfix= ["}","}","include "+sysfn]
  with open(outfn+".up.plot",'w') as outf:
    outf.write("\n".join(outlines_prefix+ [" ".join(uporblines)] + outlines_postfix))
  with open(outfn+".dn.plot",'w') as outf:
    outf.write("\n".join(outlines_prefix+ [" ".join(dnorblines)] + outlines_postfix))

###############################################################################
# f orbital normalizations are from 
# <http://winter.group.shef.ac.uk/orbitron/AOs/4f/equations.html>
def normalize_eigvec(eigvec,basis):
  snorm = 1./(4.*np.pi)**0.5
  pnorm = snorm*(3.)**.5
  dnorms = [
      .5*(5./(4*np.pi))**.5,
      (15./(4*np.pi))**.5,
      (15./(4*np.pi))**.5,
      .5*(15./(4.*np.pi))**.5,
      (15./(4*np.pi))**.5
    ]
  fnorms = [
      ( 7./(16.*np.pi))**.5,
      (21./(32.*np.pi))**.5,
      (21./(32.*np.pi))**.5,
      (105./(16.*np.pi))**.5, # xyz
      (105./(4.*np.pi))**.5,
      (35./(32.*np.pi))**.5,
      (35./(32.*np.pi))**.5
    ]

  # Duplicate coefficients for complex, and if multiple basis elements are d.
  # This is to align properly with the d-components of eigvecs.
  tmp = [[f for f in dnorms] for i in range(sum(basis['shell_type']==3))]
  dnorms = []
  for l in tmp: dnorms += l
  dnorms = np.array(dnorms)
  # Likewise for f.
  tmp = [[f for f in fnorms] for i in range(sum(basis['shell_type']==4))]
  fnorms = []
  for l in tmp: fnorms += l
  fnorms = np.array(fnorms)

  ao_type = []
  for sidx in range(len(basis['shell_type'])):
    ao_type += \
      [basis['shell_type'][sidx] for ao in range(basis['nao_shell'][sidx])]
  ao_type = np.array(ao_type)

  if any(ao_type==1):
    error("sp orbtials not implemented in normalize_eigvec(...)","Not implemented")

  eigvec[:,ao_type==0] *= snorm
  eigvec[:,ao_type==2] *= pnorm
  eigvec[:,ao_type==3] *= dnorms
  eigvec[:,ao_type==4] *= fnorms

  return eigvec
      
###############################################################################
def write_orb(eigsys,basis,ions,kpt,outfn,maxmo_spin=-1):
  outf=open(outfn,'w')
  if maxmo_spin < 0:
    maxmo_spin=basis['nmo']

  eigvecs=[normalize_eigvec(eigvec_lookup(kpt,eigsys,spin),basis) for spin in range(eigsys['nspin'])]
  atidxs = np.unique(basis['atom_shell'])-1
  nao_atom = np.zeros(atidxs.size,dtype=int)
  for shidx in range(len(basis['nao_shell'])):
    nao_atom[basis['atom_shell'][shidx]-1] += basis['nao_shell'][shidx]
  #nao_atom = int(round(sum(basis['nao_shell']) / len(ions['positions'])))
  coef_cnt = 1
  totnmo = maxmo_spin*eigsys['nspin'] #basis['nmo'] * eigsys['nspin']
  for moidx in np.arange(totnmo)+1:
    for atidx in atidxs+1:
      for aoidx in np.arange(nao_atom[atidx-1])+1:
        outf.write(" {:5d} {:5d} {:5d} {:5d}\n"\
            .format(moidx,aoidx,atidx,coef_cnt))
        coef_cnt += 1
  eigvec_flat = [e[0:maxmo_spin].flatten() for s in range(eigsys['nspin']) for e in eigvecs[s]]
  print_cnt = 0
  outf.write("COEFFICIENTS\n")
  if eigsys['ikpt_iscmpx'][kpt]: #complex coefficients
    for eigv in eigvec_flat: #zip(eigreal_flat,eigimag_flat):
      for r,i in zip(eigv.real,eigv.imag):
        outf.write("({:<.12e},{:<.12e}) "\
            .format(r,i))
        print_cnt+=1
        if print_cnt%5==0: outf.write("\n")
  else: #Real coefficients
    for eigr in eigvec_flat:
      for r in eigr:
        outf.write("{:< 15.12e} ".format(r))
        print_cnt+=1
        if print_cnt%5==0: outf.write("\n")

###############################################################################
# TODO Generalize to no pseudopotential.
def write_sys(lat_parm,basis,eigsys,pseudo,ions,kpt,outfn):
  outlines = []
  min_exp = min(basis['prim_gaus'])
  cutoff_length = (-np.log(1e-8)/min_exp)**.5
  basis_cutoff = find_basis_cutoff(lat_parm)
  cutoff_divider = basis_cutoff*2.0 / cutoff_length

  if lat_parm['struct_dim'] != 0:
    outlines += [
        "system { periodic",
        "  nspin {{ {} {} }}".format(eigsys['nup'],eigsys['ndn']),
        "  latticevec {",
      ]
    for i in range(3):
      outlines.append("    {:< 15} {:< 15} {:< 15}".format(*lat_parm['latvecs'][i]))
    outlines += [
        "  }",
        "  origin { 0 0 0 }",
        "  cutoff_divider {0}".format(cutoff_divider),
        "  kpoint {{ {:4}   {:4}   {:4} }}".format(
            *(np.array(kpt)/eigsys['nkpts_dir']*2.)
          )
      ]
  else: # is molecule.
    outlines += [
        "system { molecule",
        "  nspin {{ {} {} }}".format(eigsys['nup'],eigsys['ndn']),
      ]
  for aidx in range(len(ions['positions'])):
    if ions['atom_nums'][aidx]-200-1 < 0:
      error("All-electron calculations not implemented yet.","Not implemented")
    outlines.append(
      "  atom {{ {0} {1} coor {2} }}".format(
        periodic_table[ions['atom_nums'][aidx]-200-1], # Assumes ECP.
        ions['charges'][aidx],
        "{:< 15} {:< 15} {:< 15}".format(*ions['positions'][aidx])
      )
    )
  outlines.append("}")
  done = []
  for elem in pseudo.keys():
    atom_name = periodic_table[elem-200-1]
    n_per_j = pseudo[elem]['n_per_j']
    numL = sum(n_per_j>0)

    for i in range(1,len(n_per_j)):
      if (n_per_j[i-1]==0)and(n_per_j[i]!=0):
        error("ERROR: Weird pseudopotential, please generalize write_sys(...).",
              "Not implemented.")

    n_per_j = n_per_j[n_per_j>0]
    order = list(np.arange(n_per_j[0],sum(n_per_j))) + \
            list(np.arange(n_per_j[0])) 
    exponents   = pseudo[elem]['exponents'][order]
    prefactors  = pseudo[elem]['prefactors'][order]
    r_exps      = pseudo[elem]['r_exps'][order]
    if numL > 2: aip = 12
    else:        aip =  6
    npjline = n_per_j[1:].tolist()+[n_per_j[0]]
    outlines += [
        "pseudo {",
        "  {}".format(atom_name),
        "  aip {:d}".format(aip),
        "  basis {{ {}".format(atom_name),
        "    rgaussian",
        "    oldqmc {",
        "      0.0 {:d}".format(numL),
        "      "+' '.join(["{}" for i in range(numL)]).format(*npjline)
      ]
    cnt = 0
    for eidx in range(len(exponents)):
      outlines.append("      {:d}   {:<12} {:< 12}".format(
        r_exps[cnt]+2,
        float(exponents[cnt]),
        float(prefactors[cnt])
      ))
      cnt += 1
    outlines += ["    }","  }","}"]
  with open(outfn,'w') as outf:
    outf.write("\n".join(outlines))

###############################################################################
def write_jast2(lat_parm,ions,outfn):
  basis_cutoff = find_basis_cutoff(lat_parm)
  atom_types = [periodic_table[eidx-200-1] for eidx in ions['atom_nums']]
  atom_types=set(atom_types)
  outlines = [
      "jastrow2",
      "group {",
      "  optimizebasis",
      "  eebasis {",
      "    ee",
      "    cutoff_cusp",
      "    gamma 24.0",
      "    cusp 1.0",
      "    cutoff {0}".format(basis_cutoff),
      "  }",
      "  eebasis {",
      "    ee",
      "    cutoff_cusp",
      "    gamma 24.0",
      "    cusp 1.0",
      "    cutoff {0}".format(basis_cutoff),
      "  }",
      "  twobody_spin {",
      "    freeze",
      "    like_coefficients { 0.25 0.0 }",
      "    unlike_coefficients { 0.0 0.5 }",
      "  }",
      "}",
      "group {",
      "  optimize_basis",
    ]
  for atom_type in atom_types:
    outlines += [
      "  eibasis {",
      "    {0}".format(atom_type),
      "    polypade",
      "    beta0 0.2",
      "    nfunc 3",
      "    rcut {0}".format(basis_cutoff),
      "  }"
    ]
  outlines += [
      "  onebody {",
    ]
  for atom_type in atom_types:
    outlines += [
      "    coefficients {{ {0} 0.0 0.0 0.0}}".format(atom_type),
    ]
  outlines += [
      "  }",
      "  eebasis {",
      "    ee",
      "    polypade",
      "    beta0 0.5",
      "    nfunc 3",
      "    rcut {0}".format(basis_cutoff),
      "  }",
      "  twobody {",
      "    coefficients { 0.0 0.0 0.0 }",
      "  }",
      "}"
    ]
  with open(outfn,'w') as outf:
    outf.write("\n".join(outlines))

###############################################################################
def write_basis(basis,ions,outfn):
  hybridized_check = 0.0
  hybridized_check += sum(abs(basis['coef_s'] * basis['coef_p']))
  hybridized_check += sum(abs(basis['coef_p'] * basis['coef_dfg']))
  hybridized_check += sum(abs(basis['coef_s'] * basis['coef_dfg']))
  if hybridized_check > 1e-10:
    error("Hybridized AOs (like sp) not implmemented in write_basis(...)",
          "Not implemented.")

  # If there's no hybridization, at most one of coef_s, coef_p, and coef_dfg is
  # nonzero. Just add them, so we have one array.
  coefs = basis['coef_s'] + basis['coef_p'] + basis['coef_dfg']

  shell_type = np.tile("Unknown...",basis['shell_type'].shape)
  typemap = ["S","SP","P","5D","7F_crystal","G","H"]
  for i in range(5): shell_type[basis['shell_type']==i] = typemap[i]

  cnt = 0
  aidx = 0
  atom_type = ions['atom_nums'][aidx]
  outlines = [
      "basis {",
      "  {0}".format(periodic_table[atom_type-200-1]),
      "  aospline",
      "  normtype CRYSTAL",
      "  gamess {"
    ]
  done_atoms = set([atom_type])
  for sidx in range(len(shell_type)):
    nprim = basis['prim_shell'][sidx]
    new_aidx = basis['atom_shell'][sidx]-1

    new_atom_type = ions['atom_nums'][new_aidx]
    if aidx != new_aidx:
      if new_atom_type in done_atoms:
        cnt+=nprim
        continue
      else:
        outlines += ["  }","}"]
        atom_type = new_atom_type
        aidx = new_aidx
        outlines += [
            "basis {",
            "  {0}".format(periodic_table[atom_type-200-1]),
            "  aospline",
            "  normtype CRYSTAL",
            "  gamess {"
          ]
        done_atoms.add(atom_type)

    outlines.append("    {0} {1}".format(shell_type[sidx],nprim))
    for pidx in range(nprim):
      outlines.append("      {0} {1} {2}".format(
        pidx+1,
        basis['prim_gaus'][cnt],
        coefs[cnt]
      ))
      cnt += 1
  outlines += ["  }","}"]
  with open(outfn,'w') as outf:
    outf.write("\n".join(outlines))

###############################################################################
def write_moanalysis():
  raise NotImplementedError()
  return None

if __name__ == "__main__":
  from argparse import ArgumentParser
  parser=ArgumentParser('Convert a crystal file (defaults in brackets).')
  parser.add_argument('-b','--base',type=str,default='qwalk',
      help="[='qw'] First part of string for file names.")
  parser.add_argument('-p','--propout',type=str,default='prop.in.o',
      help="[='prop.in.o'] Name of file containing net spin. Either crystal or properties stdout.")
  parser.add_argument('-r','--real',type=bool,default=False,
      help="[=False] Convert only real kpoints.")
  parser.add_argument('-v','--nvirtual',type=int,default=50,
      help="[=50] Number of unoccupied or virtual orbitals to allow access to.")
  args=parser.parse_args()

  convert_crystal(args.base,args.propout,args.kset,args.nvirtual)

