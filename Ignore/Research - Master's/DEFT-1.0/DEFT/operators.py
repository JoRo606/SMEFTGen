import collections
import itertools
import sympy
from sympy.matrices.sparse import SparseMatrix

import DEFT
from DEFT import frac,Field,Term,Invariant,Index,indexConsistency,mash

class Terms:
  def __getitem__(self,index):
    return (self.terms[index],self.weights[index])
  def __add__(self,otherterms):
    return Terms(self.terms+otherterms.terms,self.weights+otherterms.weights)
  def __mul__(self,scalar):
    return Terms(self.terms,[scalar*w for w in self.weights])
  def __hash__(self):
    hsh=0
    for t in self.terms: hsh |= hash(t)
    for w in self.weights: hsh |= hash(w)
    return hsh
  def collectTerms(self):
    uniqueTerms=[]
    uniqueWeights=[]
    for t,w in self:
      try:
        ind=uniqueTerms.index(t)
      except ValueError:
        uniqueTerms.append(t)
        uniqueWeights.append(w)
        continue
      uniqueWeights[ind]+=uniqueTerms[ind].compareSign(t)*w
    return Terms(uniqueTerms,uniqueWeights)
  def filterByWeight(self,weightFunc):
    filttermsnweights=[(t,w) for t,w in self if weightFunc(w)]
    return Terms(*zip(*filttermsnweights))
  def filterByTerm(self,termFunc):
    filttermsnweights=[(t,w) for t,w in self if termFunc(t)]
    return Terms(*zip(*filttermsnweights))
  def shortstr(self):
    return '  '.join([str(w)+' '+t.shortstr() for w,t in zip(self.weights,self.terms)])
  def textRep(self,suppress=[]):
    return '\n'.join(str(w)+'\n'+t.textRep(suppress) for w,t in zip(self.weights,self.terms))
  def texRep(self,suppress=[]):
    #return ' ~ '.join([str(w)+' '+t.texRep(suppress) for w,t in zip(self.weights,self.terms)])
    return ' +~ '.join([str(w)+' '+t.texRep(suppress) for w,t in zip(self.weights,self.terms)])
  def __init__(self,tms,wgts):
    self.terms=list(tms)
    self.weights=list(wgts)

class Relation:
  def __eq__(self,otherrel):
    return self.id==otherrel.id #and self.type==otherrel.type and set(zip([mash(t) for t in self.terms],self.weights))==set(zip([mash(t) for t in otherrel.terms],otherrel.weights))
  def __str__(self):
    return self.type+': '+'   '.join(t.shortstr() for t in self.terms)+('\n'+'  '.join(str(w) for w in self.weights) if self.weights is not None else '')
  def texRep(self):
    return self.type+': '+' ~ '.join(str(w)+' '+t.texRep() for t,w in zip(self.terms,self.weights) )
  ##maybe dodgy to conjugate weights?
  #def conjugate(self,flipWeights=False):
  #  newrel=Relation(self.type,[t.conjugate() for t in self.terms],None if self.weights is None else self.weights[:])
  #  if flipWeights and self.weights is not None:
  #    for i in range(len(self.terms)):
  #      newrel.weights[i]=newrel.weights[i]*newrel.terms[i].isign()*self.terms[i].isign()
  #  return newrel
  def removeDups(self):
    uniqueTerms=[]
    uniqueWeights=[]
    for i,t in enumerate(self.terms):
      try:
        ind=uniqueTerms.index(t)
      except ValueError:
        uniqueTerms.append(t)
        uniqueWeights.append(self.weights[i])
        continue
      uniqueWeights[ind]+=self.weights[i]
    self.terms=uniqueTerms
    self.weights=uniqueWeights
    self.removeZeroWeights()
    return
  def removeZeroWeights(self):
    nonZeroTs=[]
    nonZeroWs=[]
    for i,w in enumerate(self.weights): self.weights[i] = sympy.simplify(w)
    for i,t in enumerate(self.terms):
      if self.weights[i]==0: continue
      nonZeroTs.append(t)
      nonZeroWs.append(self.weights[i])
    self.terms=nonZeroTs
    self.weights=nonZeroWs
    return
  def removeUnsyms(self):
    nonZeroTs=[]
    nonZeroWs=[]
    for i,t in enumerate(self.terms):
      if not t.checkSyms(): continue
      nonZeroTs.append(t)
      nonZeroWs.append(self.weights[i])
    self.terms=nonZeroTs
    self.weights=nonZeroWs
    return
  def addbetter(self,deadrel,deadterm):
    if id(self)==id(deadrel): return
    selfindex=next((i for i,t in enumerate(self.terms) if mash(deadterm)==mash(t)),None)
    deadindex=next((i for i,t in enumerate(deadrel.terms) if mash(deadterm)==mash(t)),None)
    selfwght=self.weights[ selfindex ]
    deadwght=deadrel.weights[ deadindex ]
    #remove deadterm
    self.terms.pop(selfindex)
    self.weights.pop(selfindex)
    for deadindex2,tm in enumerate(deadrel.terms):
      if deadindex==deadindex2: continue
      selfindex=next((i for i,t in enumerate(self.terms) if mash(tm)==mash(t)),None)
      if selfindex is None:
        self.terms.append(tm)
        self.weights.append( -selfwght/deadwght * deadrel.weights[deadindex2])
      else:
        self.weights[selfindex]=sympy.simplify(self.weights[selfindex]-selfwght/deadwght * deadrel.weights[deadindex2])
    self.type+=deadrel.type
    return
  def __init__(self,typ,terms,weights):
    if len(terms)!=len(weights): raise Exception('There isn\'t one weight per term')
    self.type=typ
    self.terms=terms
    self.weights=weights
    self.id=None


def IBPRels(term):
  rels=[]
  for fieldi,f in enumerate(term):
    if len(f.Dindices)==0: continue
    rel=[]
    for fieldj in range(len(term._fields)):
      t=term._DEFTcopy()
      dis=t[fieldi].Dindices.pop()
      t[fieldi].remove(dis[0])
      t[fieldi].remove(dis[1])
      t[fieldj].Dindices.insert(0,dis)
      #fudge to keep signs consistent
      if len(t[fieldj].Dindices)>1:
        posinindexlist=t[fieldj].indices.index(t[fieldj].Dindices[1][0])
        t[fieldj].indices.insert(posinindexlist,dis[1])
        t[fieldj].indices.insert(posinindexlist,dis[0])
      else:
        t[fieldj].indices.append(dis[0])
        t[fieldj].indices.append(dis[1])
      rel.append(t)
    rels.append(rel)
  return rels

def removeDups(lst):
  newLst=[]
  for i in lst:
    i._clear()
    if i not in newLst: newLst.append(i)
  return newLst

def removeDupsWMash(lst):
  newLst=[]
  for t in lst: t._clear()
  #for t in lst: t._mash=None
  for t in lst:
    if mash(t) not in [mash(tm) for tm in newLst]: newLst.append(t)
  return newLst


def generateAllFreeTerms(listOfFields,intermediateTermCheck,secretSpeedUps=False,finalTermCheck=None):
  def makeTermFromClasses(lst):
    t=Term()
    for f in lst:
      t.add(f(t))
    return t
  def makeTermFromFields(lst):
    t=Term()
    for f in lst:
      t.add(f)
    return t
  def conjugateAll(fieldlst):
    fieldLists=[[f._DEFTcopy() for f in fieldlst] for i in range(2**len(fieldlst) - 1)]
    fieldLists.append(fieldlst)
    for i in range(2**len(fieldlst)):
      for j,f in enumerate(fieldLists[i]):
        if (i/2**j)% 2==1: f.conjugate()
      #terms[i].refreshUOnes()
    return fieldLists
  def diffAll(term):
    terms=[term._DEFTcopy() for i in range(len(term._fields))]
    for i in range(len(terms)):
      #terms[i][i+1].differentiate()
      terms[i][i].differentiate()
      terms[i].dimension+=frac(1,1)
    return terms
  n=1
  listOfTerms=[]
  while True:
    newTerms=[]
    combs=list(itertools.combinations_with_replacement(listOfFields,n))
    ##filter combs by intermediate term check
    filteredcombs=[c for c in combs if intermediateTermCheck(makeTermFromClasses(c))]
    if len(list(combs))>0 and len(filteredcombs)==0: break
    for c in filteredcombs:
      t=makeTermFromClasses(c)
      if not intermediateTermCheck(t):
        continue
      #print t.shortstr()
      fieldList=[f() for f in c]
      if secretSpeedUps: newTerms+=filter(zeroUOnes,removeDups([makeTermFromFields(fl) for fl in conjugateAll(fieldList)]))
      else: newTerms+=removeDups([makeTermFromFields(fl) for fl in conjugateAll(fieldList)])
    #if secretSpeedUps:
    #  if len(newTerms)==0 and n>4: break
    #else:
    #  if len(newTerms)==0: break
    TBDTerms=newTerms
    while True:
      derivativeTerms=[]
      for t in TBDTerms:
        derivativeTerms+=diffAll(t)
      derivativeTerms=filter(intermediateTermCheck,derivativeTerms)
      if len(derivativeTerms)==0: break
      derivativeTerms=removeDups(derivativeTerms)
      newTerms+=derivativeTerms
      TBDTerms=derivativeTerms
    n+=1
    listOfTerms+=newTerms
  if finalTermCheck is not None: listOfTerms=filter(finalTermCheck,listOfTerms)
  return listOfTerms

def zeroUOnes(term,ignore=['B','L']):
  return all(frac.frac[0]==0 for key,frac in term.uones.iteritems() if key not in ignore)

def combComp(inds,ngroup):
  LCR=[([],[],inds)]
  for n in range(ngroup):
    newLCR=[]
    for lcr in LCR:
      for i in range(len(lcr[2])-1):
        newLCR.append( (lcr[0]+lcr[2][:i],lcr[1]+[lcr[2][i]],lcr[2][i+1:]) )
      if len(lcr[2])>0:
        newLCR.append( (lcr[0]+lcr[2][:-1],lcr[1]+[lcr[2][-1]],[]) )
    LCR=newLCR
  newLCR=[]
  for lcr in LCR:
    newLCR.append( (lcr[1],lcr[0]+lcr[2]) )
  return newLCR
def splits(inds,ngroup,allowPossibleDups=False):
  splits=[((),inds)]
  totalsplits=splits
  for i in range(len(inds)//ngroup):
    newsplits=[]
    for s in splits:
      newbitsandrest=combComp(s[1],ngroup)
      for newbit,rest in newbitsandrest:
        newsplits.append( ( s[0]+(newbit,),rest) )
    splits=newsplits
    totalsplits+=newsplits
  if allowPossibleDups: newsplits=[]
  else: newsplits=set()
  for s in totalsplits:
    if allowPossibleDups: newsplits.append( frozenset( frozenset(g) for g in s[0] ) )
    else: newsplits.add( frozenset( frozenset(g) for g in s[0] ) )
  return newsplits

def fastContractOneGroup(groupinds,allowPossibleDups=False):
  ninds=DEFT.groupsDimFund[groupinds[0].group]
  upperinds=[i for i in groupinds if i.upper]
  lowerinds=[i for i in groupinds if not i.upper]
  def pe(gi):
    if len(gi)==0: return []
    #if gi[0].group=='su3' or gi[0].group=='sug': ninds=3
    #else: ninds=2
    return list(itertools.combinations(gi,ninds))
  def mutuallyExclusive(prod):
    inds=list(itertools.chain.from_iterable(list(p) for p in prod))
    return len(set(inds))==len(list(inds))
  def pegproper(inds):
    if len(inds)==0: return [()]
    #if inds[0].group=='su3' or inds[0].group=='sug': ninds=3
    #else: ninds=2
    return list(splits(inds,ninds,allowPossibleDups))
  def peg(pe,inds):
    possibleGroups=[()]
    for n in range(1,len(inds) // 2 + 1):
      possibleGroups+=[tuple(p) for p in itertools.combinations(pe,n) if mutuallyExclusive(p)]
    return possibleGroups
  possibleEUGroups=pegproper(upperinds)
  possibleELGroups=pegproper(lowerinds)
  possibleContractions=[]
  for pegu in possibleEUGroups:
    upperflat=list(itertools.chain.from_iterable(pegu))
    uppercomp=[i for i in upperinds if i not in upperflat]
    for pegl in possibleELGroups:
      lowerflat=list(itertools.chain.from_iterable(pegl))
      lowercomp=[i for i in lowerinds if i not in lowerflat]
      if len(uppercomp)!=len(lowercomp): continue
      if len(uppercomp)==0: possibleContractions.append( (pegu,pegl,()) )
      else:
        possibleLowerPerms=itertools.permutations(lowercomp)
        for plp in possibleLowerPerms:
           pdg=tuple(zip(uppercomp,plp))
           possibleContractions.append( (pegu,pegl,pdg) )
  return possibleContractions

def mashConverter(itr):
  if hasattr(itr, '__iter__'):
    if isinstance(itr,frozenset):
      return sorted([mashConverter(i) for i in itr],key=str)
    else:
      return tuple(mashConverter(i) for i in itr)
  else:
    return itr
def mashConverter2(itr):
  if hasattr(itr, '__iter__'):
    if isinstance(itr,frozenset):
      return sorted([mashConverter(i) for i in itr])
    else:
      return tuple(sorted([mashConverter(i) for i in itr]))
  else:
    return itr
def mashConverter3(itr):
  return [(mashConverter2(i[0]),mashConverter2(i[1]),i[2]) for i in itr]

def fastContract(term,allowPossibleDups=False,intermediateSimplify=False):
  freeIndices=term.freeIndices()
  #print 'freeinds',freeIndices
  freeIndices.sort(key=lambda x: x.group)
  indicesByGroup=itertools.groupby(freeIndices,key=lambda x : x.group)
  indicesByGroup2=itertools.groupby(freeIndices,key=lambda x : x.group)
  #print 'groupinds',[list(groupinds) for k,groupinds in indicesByGroup2]
  possibleContractions=[fastContractOneGroup(list(groupinds),allowPossibleDups) for k,groupinds in indicesByGroup]
  contractedTerms=[]
  terms=[term._DEFTcopy()]
  for pcg in possibleContractions:
    newterms=[]
    for pc in pcg:
      for t in terms:
        newterm=t._DEFTcopy()
        ntfis=newterm.freeIndices()
        for inds in pc[0]:
          actualinds=[next(i for i in ntfis if i==iv) for iv in inds]
          newterm.addInvariant( Invariant('epsu','\epsilon',actualinds) )
        for inds in pc[1]:
          actualinds=[next(i for i in ntfis if i==iv) for iv in inds]
          newterm.addInvariant( Invariant('epsl','\epsilon',actualinds) )
        for inds in pc[2]:
          actualinds=[next(i for i in ntfis if i==iv) for iv in inds]
          newterm.addInvariant( Invariant('del','\delta',actualinds) )
        if intermediateSimplify:
          if newterm.checkSymsSimp(): newterms.append( newterm)
          else: continue
        else: newterms.append( newterm)
    terms=removeDupsWMash(newterms)
  return terms


#BCTI compliant 2016-02-26
def fierzRelsEUEU(term,group):
  def findIndex(dummy,t):
    return next(i for i in t.dummyIndices() if i==dummy)
  epsu = [x for x in term._invariants if x.name=='epsu' and x.indices[0].group==group] 
  rels=[]
  for eupair in itertools.combinations(epsu,2):
    newterm1=term._DEFTcopy()
    newterm2=term._DEFTcopy()
    for eu in eupair:
        for newterm in [newterm1,newterm2]:
          newterm._clear()
          newterm._invariants.remove(eu)
    newterm1.addInvariant(Invariant('epsu','\epsilon',[findIndex(eupair[0].indices[0],newterm1),findIndex(eupair[1].indices[0],newterm1)]))
    newterm1.addInvariant(Invariant('epsu','\epsilon',[findIndex(eupair[0].indices[1],newterm1),findIndex(eupair[1].indices[1],newterm1)]))
    newterm2.addInvariant(Invariant('epsu','\epsilon',[findIndex(eupair[0].indices[0],newterm2),findIndex(eupair[1].indices[1],newterm2)]))
    newterm2.addInvariant(Invariant('epsu','\epsilon',[findIndex(eupair[0].indices[1],newterm2),findIndex(eupair[1].indices[0],newterm2)]))
    rels.append( ([term,newterm1,newterm2],[-1*term.isign(),1*newterm1.isign(),-1*newterm2.isign()]) )
  return rels

#BCTI compliant 2016-02-26
def fierzRelsELEL(term,group):
  def findIndex(dummy,t):
    return next(i for i in t.dummyIndices() if i==dummy)
  epsl = [x for x in term._invariants if x.name=='epsl' and x.indices[0].group==group] 
  rels=[]
  for elpair in itertools.combinations(epsl,2):
    newterm1=term._DEFTcopy()
    newterm2=term._DEFTcopy()
    for el in elpair:
        for newterm in [newterm1,newterm2]:
          newterm._clear()
          newterm._invariants.remove(el)
    newterm1.addInvariant(Invariant('epsl','\epsilon',[findIndex(elpair[0].indices[0],newterm1),findIndex(elpair[1].indices[0],newterm1)]))
    newterm1.addInvariant(Invariant('epsl','\epsilon',[findIndex(elpair[0].indices[1],newterm1),findIndex(elpair[1].indices[1],newterm1)]))
    newterm2.addInvariant(Invariant('epsl','\epsilon',[findIndex(elpair[0].indices[0],newterm2),findIndex(elpair[1].indices[1],newterm2)]))
    newterm2.addInvariant(Invariant('epsl','\epsilon',[findIndex(elpair[0].indices[1],newterm2),findIndex(elpair[1].indices[0],newterm2)]))
    rels.append( ([term,newterm1,newterm2],[-1*term.isign(),1*newterm1.isign(),-1*newterm2.isign()]) )
  return rels

#BCTI compliant 2016-02-26
def fierzRelsEUEL(term,group):
  def findIndex(dummy,t):
    return next(i for i in t.dummyIndices() if i==dummy)
  epsu = [x for x in term._invariants if x.name=='epsu' and x.indices[0].group==group] 
  epsl = [x for x in term._invariants if x.name=='epsl' and x.indices[0].group==group]
  rels=[]
  for eu in epsu:
    for el in epsl:
      newTerms=[]
      newWeights=[]
      for elindexperm in itertools.permutations(el.indices):
        indexpairs=zip(eu.indices,elindexperm)
        newterm=term._DEFTcopy()
        newterm._clear()
        newterm._invariants.remove(eu)
        newterm._invariants.remove(el)
        for ip in indexpairs:
          newterm.addInvariant(Invariant('del','\delta',[findIndex(ip[0],newterm),findIndex(ip[1],newterm)]))
        newTerms.append(newterm)
        newWeights.append( newterm.isign() * (1 if newterm._cmpParity(el.indices,elindexperm) else -1) )
      newTerms.append(term)
      newWeights.append( term.isign() )
      rels.append( (newTerms,newWeights) )
  return rels

#BCTI compliant 2016-02-26
def fierzRelsDelEL(term,group):
  def findIndex(dummy,t):
    return next(i for i in t.dummyIndices() if i==dummy)
  dls = [x for x in term._invariants if x.name=='del' and x.indices[0].group==group] 
  epsl = [x for x in term._invariants if x.name=='epsl' and x.indices[0].group==group] 
  rels=[]
  for dl in dls:
    for el in epsl:
      newterm1=term._DEFTcopy()
      newterm2=term._DEFTcopy()
      for newterm in [newterm1,newterm2]:
        newterm._clear()
        newterm._invariants.remove(dl)
        newterm._invariants.remove(el)
      newterm1.addInvariant(Invariant('del','\delta',[findIndex(dl.indices[0],newterm1),findIndex(el.indices[0],newterm1)]))
      newterm1.addInvariant(Invariant('epsl','\epsilon',[findIndex(dl.indices[1],newterm1),findIndex(el.indices[1],newterm1)]))
      newterm2.addInvariant(Invariant('del','\delta',[findIndex(dl.indices[0],newterm2),findIndex(el.indices[1],newterm2)]))
      newterm2.addInvariant(Invariant('epsl','\epsilon',[findIndex(dl.indices[1],newterm2),findIndex(el.indices[0],newterm2)]))
      #rels.append( ([term,newterm1,newterm2],[-1,1,-1]) )
      rels.append( ([term,newterm1,newterm2],[-1*term.isign(),1*newterm1.isign(),-1*newterm2.isign()]) )
  return rels

#BCTI compliant 2016-02-26
def fierzRelsDelEU(term,group):
  def findIndex(dummy,t):
    return next(i for i in t.dummyIndices() if i==dummy)
  dls = [x for x in term._invariants if x.name=='del' and x.indices[0].group==group] 
  epsu = [x for x in term._invariants if x.name=='epsu' and x.indices[0].group==group] 
  rels=[]
  for dl in dls:
    for eu in epsu:
      newterm1=term._DEFTcopy()
      newterm2=term._DEFTcopy()
      for newterm in [newterm1,newterm2]:
        newterm._clear()
        newterm._invariants.remove(dl)
        newterm._invariants.remove(eu)
      newterm1.addInvariant(Invariant('del','\delta',[findIndex(eu.indices[1],newterm1),findIndex(dl.indices[1],newterm1)]))
      newterm1.addInvariant(Invariant('epsu','\epsilon',[findIndex(eu.indices[0],newterm1),findIndex(dl.indices[0],newterm1)]))
      newterm2.addInvariant(Invariant('del','\delta',[findIndex(eu.indices[0],newterm2),findIndex(dl.indices[1],newterm2)]))
      newterm2.addInvariant(Invariant('epsu','\epsilon',[findIndex(eu.indices[1],newterm2),findIndex(dl.indices[0],newterm2)]))
      rels.append( ([term,newterm1,newterm2],[-1*term.isign(),1*newterm1.isign(),-1*newterm2.isign()]) )
  return rels

#BCTI compliant 2016-02-18
def switchDerivatives(term):
  def findIndex(dummy,t):
    return next(i for i in t.dummyIndices() if i==dummy)
  rels=[]
  for fieldi,f in enumerate(term):
    if len(f.Dindices)<2: continue
    for i in range(len(f.Dindices)-1):
      rel=[]
      relweights=[]
      a=f.Dindices[i][0]
      adot=f.Dindices[i][1]
      b=f.Dindices[i+1][0]
      bdot=f.Dindices[i+1][1]
      template=term._DEFTcopy()
      template._fields[fieldi].Dindices.pop(i)
      template._fields[fieldi].Dindices.pop(i)
      #template._fields[fieldi].dimension-=frac(2,1)
      #for Bfs,Y,gs in [('B','Y',sympy.symbols('gp'))]:
      for Y,gs in DEFT.U1GaugeGroupsCoupling.iteritems():
        Bfs=DEFT.gaugeGroupsFields[Y]
      #for Bfs,Y,gs in [('F','Z',sympy.symbols('gz'))]:
      #for Bfs,Y,gs in []:#[('B','Y',sympy.symbols('gp')),('F','Z',sympy.symbols('gz'))]:
        try:
          if f.U1Dict[Y]==frac(0,1): continue
        except KeyError:
          continue
        #F term
        t=template._DEFTcopy()
        fs=Bfs()
        fs.indices=[findIndex(a,t),findIndex(b,t)]
        fs.symmetries[0].indices=fs.indices
        t._fields.insert(fieldi,fs)
        fs._term=t
        t._invariants.append(Invariant('epsu','\epsilon',[findIndex(adot,t),findIndex(bdot,t)]))
        for ind in [a,adot,b,bdot]:
          t._fields[fieldi+1].remove(findIndex(ind,t))
        tmpsgn=t.contractInvariantIndices()
        t.fixTwoIndexInvs()
        rel.append(t)
        relweights.append( tmpsgn * sympy.I * gs * f.U1Dict[Y].frac[0] / f.U1Dict[Y].frac[1] / 2)
        #F bar term
        t=template._DEFTcopy()
        fs=Bfs()
        fs.conjugate()
        fs.indices=[findIndex(adot,t),findIndex(bdot,t)]
        fs.symmetries[0].indices=fs.indices
        t._fields.insert(fieldi,fs)
        fs._term=t
        t._invariants.append(Invariant('epsu','\epsilon',[findIndex(a,t),findIndex(b,t)]))
        for ind in [a,adot,b,bdot]:
          t._fields[fieldi+1].remove(findIndex(ind,t))
        tmpsgn=t.contractInvariantIndices()
        t.fixTwoIndexInvs()
        rel.append(t)
        relweights.append( tmpsgn * sympy.I * gs * f.U1Dict[Y].frac[0] / f.U1Dict[Y].frac[1]/  2)
      for grp,gs in DEFT.NAGaugeGroupsCoupling.iteritems():
        Wfs=DEFT.gaugeGroupsFields[grp]
      #for Wfs,grp,gs in [('W','su2',sympy.symbols('g')),('G','su3',sympy.symbols('gs'))]:
        for grpind in [gi for gi in f.indices if gi.group==grp]:
          #F term
          t=template._DEFTcopy()
          digup=Index(grp,term.smallestNewLabel(grp),True,True)
          diglow=Index(grp,term.smallestNewLabel(grp),False,True)
          fs=Wfs()
          fs.indices=([findIndex(grpind,t),diglow] if grpind.upper else [digup,findIndex(grpind,t)])+[findIndex(a,t),findIndex(b,t)]
          fs.symmetries[0].indices=fs.indices[-2:]
          fs.symmetries[1].indices=[(findIndex(grpind,t) if grpind.upper else digup),(findIndex(grpind,t) if not grpind.upper else diglow)]
          t._fields.insert(fieldi,fs)
          fs._term=t
          if grpind.upper:
            t[fieldi+1].replaceIndexInField(findIndex(grpind,t),digup)
            t[fieldi+1].replaceIndexInSym(findIndex(grpind,t),digup)
          else:
            t[fieldi+1].replaceIndexInField(findIndex(grpind,t),diglow)
            t[fieldi+1].replaceIndexInSym(findIndex(grpind,t),diglow)
          t._invariants.append(Invariant('del','\delta',[digup,diglow]))
          t._invariants.append(Invariant('epsu','\epsilon',[findIndex(adot,t),findIndex(bdot,t)]))
          for ind in [a,adot,b,bdot]:
            t._fields[fieldi+1].remove(findIndex(ind,t))
          tmpsgn=t.contractInvariantIndices()
          t.fixTwoIndexInvs()
          rel.append(t)
          relweights.append( (1 if grpind.upper else -1) * tmpsgn * sympy.I * gs /2)
          #F bar term
          t=template._DEFTcopy()
          digup=Index(grp,term.smallestNewLabel(grp),True,True)
          diglow=Index(grp,term.smallestNewLabel(grp),False,True)
          fs=Wfs()
          fs.conjugate()
          fs.indices=([diglow,findIndex(grpind,t)] if grpind.upper else [findIndex(grpind,t),digup])+[findIndex(adot,t),findIndex(bdot,t)]
          fs.symmetries[0].indices=fs.indices[-2:]
          fs.symmetries[1].indices=[(findIndex(grpind,t) if grpind.upper else digup),(findIndex(grpind,t) if not grpind.upper else diglow)]
          t._fields.insert(fieldi,fs)
          fs._term=t
          if grpind.upper:
            t[fieldi+1].replaceIndexInField(findIndex(grpind,t),digup)
            t[fieldi+1].replaceIndexInSym(findIndex(grpind,t),digup)
          else:
            t[fieldi+1].replaceIndexInField(findIndex(grpind,t),diglow)
            t[fieldi+1].replaceIndexInSym(findIndex(grpind,t),diglow)
          t._invariants.append(Invariant('del','\delta',[digup,diglow]))
          t._invariants.append(Invariant('epsu','\epsilon',[findIndex(a,t),findIndex(b,t)]))
          for ind in [a,adot,b,bdot]:
            t._fields[fieldi+1].remove(findIndex(ind,t))
          tmpsgn=t.contractInvariantIndices()
          t.fixTwoIndexInvs()
          rel.append(t)
          relweights.append( (1 if grpind.upper else -1) *  tmpsgn * sympy.I * gs /2)
      ###higher derivative terms
      if i>0:
        newrel=[]
        newrelweights=[]
        for k,rterm in enumerate(rel):
          if len(rterm[fieldi+1].Dindices)>i: posinindexlist=rterm[fieldi+1].indices.index(rterm[fieldi+1].Dindices[i][0])
          else: posinindexlist=0
          prods=itertools.product([fieldi,fieldi+1], repeat=i)
          for p in prods:
            t=rterm._DEFTcopy()
            dIndices=t[fieldi+1].Dindices[:i]
            t[fieldi+1].Dindices=t[fieldi+1].Dindices[i:]
            for di in dIndices:
              t[fieldi+1].remove(di[1])
              t[fieldi+1].remove(di[0])
            fsdinds=[di for j,di in enumerate(dIndices) if p[j]==fieldi]
            fielddinds=[di for j,di in enumerate(dIndices) if p[j]==fieldi+1]
            t[fieldi].Dindices=fsdinds+t[fieldi].Dindices
            t[fieldi+1].Dindices=fielddinds+t[fieldi+1].Dindices
            t[fieldi].indices+=list(itertools.chain.from_iterable(fsdinds))
            t[fieldi+1].indices=t[fieldi+1].indices[:posinindexlist]+list(itertools.chain.from_iterable(fielddinds))+t[fieldi+1].indices[posinindexlist:]
            newrel.append(t)
            newrelweights.append(relweights[k])
        rel=newrel
        relweights=newrelweights
      ###switched term
      t=term._DEFTcopy()
      t[fieldi].Dindices[i],t[fieldi].Dindices[i+1]= t[fieldi].Dindices[i+1],t[fieldi].Dindices[i]
      #crap to flip indices around also in f.indices, so that it doesn't cock up signs
      dotpos1=t[fieldi].indices.index(t[fieldi].Dindices[i][1])
      dotpos2=t[fieldi].indices.index(t[fieldi].Dindices[i+1][1])
      t[fieldi].indices.pop(dotpos1)
      t[fieldi].indices.insert(dotpos1,t[fieldi].Dindices[i+1][1])
      t[fieldi].indices.pop(dotpos2)
      t[fieldi].indices.insert(dotpos2,t[fieldi].Dindices[i][1])
      pos1=t[fieldi].indices.index(t[fieldi].Dindices[i][0])
      pos2=t[fieldi].indices.index(t[fieldi].Dindices[i+1][0])
      t[fieldi].indices.pop(pos1)
      t[fieldi].indices.insert(pos1,t[fieldi].Dindices[i+1][0])
      t[fieldi].indices.pop(pos2)
      t[fieldi].indices.insert(pos2,t[fieldi].Dindices[i][0])
      rel.append(t)
      relweights.append( -1 )
      #unaltered term
      rel.append(term)
      relweights.append( 1 )
      rels.append( (rel,relweights) )
  return rels

from DEFT.fields import *
listOfSMFields=[H,B,ER,LL,QL,UR,DR,W,G]
listOfSMFields3F=[H,B,ER1,ER2,ER3,LL1,LL2,LL3,QL1,QL2,QL3,UR1,UR2,UR3,DR1,DR2,DR3,W,G]
listOfSMFields3G=[H,B,ERi,LLi,QLi,URi,DRi,W,G]
listOfSpurionFields3G=[Gs3,Gs6,Gs8,Gs10,Gs15p,Gs15,Gs24,Gs27]


def evenNumbersOfIndsBetter(term):
  allIndices=[]
  for f in term: allIndices+=f.indices
  for group,size in DEFT.groupsDimFund.iteritems():
    upperCount=sum(1 for i in allIndices if i.group==group and i.upper)
    lowerCount=sum(1 for i in allIndices if i.group==group and not i.upper)
    if (upperCount % size)!=(lowerCount % size): return False
  return True

def generateAllTerms(listOfFields,intermediateTermCheck,intermediateSimplify=False,finalTermCheck=None):
  #def actualTermCheck(t):
  #  return termCheck(t) and zeroUOnes(t) and evenNumbersOfIndsBetter(t)
  print "INFO: Generating free terms"
  terms= generateAllFreeTerms(listOfFields,intermediateTermCheck,True,finalTermCheck)
  #print "INFO: Filtering by U(1)s"
  terms=filter(zeroUOnes,terms)
  #print "INFO: Filtering by numbers of indices"
  terms=filter(evenNumbersOfIndsBetter,terms)
  #print 'uncontermstr',[t.shortstr() for t in terms]
  print "INFO: Contracting"
  contractedTerms=[]
  for iterm,term in enumerate(terms):
    printProgressBar(iterm,len(terms))
    contractedTerms+=fastContract(term,False,True)
  print '\r'+'*'*25+'contracted'+'*'*26
  contractedTerms=[t for t in contractedTerms if t.checkSyms()]
  #print 'contermstr',[t.shortstr() for t in contractedTerms]
  #print 'INFO: Generated',len(contractedTerms),'terms'
  return contractedTerms

def dimLEQ(n):
  def fun(term):
    return term.dimension <= frac(n,1)
  return fun


def makePDFDoc(latex,filename='DEFT'):
  import os
  import subprocess
  starttext=r'''
\documentclass[a4paper,12pt]{article}
\pdfoutput=1
\usepackage{fullpage}
\begin{document}   

'''
  endtext=r'''

\end{document}
'''

  texname=filename+'.tex'
  outfile=open(texname,'w')
  outfile.write(starttext+latex+endtext)
  outfile.close()
  subprocess.call('pdflatex -interaction=batchmode '+texname,shell=True)
  for ext in ['.log','.aux']:
    subprocess.call('rm '+filename+ext,shell=True)
  if os.path.isfile(filename+'.pdf'):
    subprocess.call('rm '+filename+'.tex',shell=True)
  return

def funcDiff(terms,field,coeffs=None):
  eom=[]
  weights=[]
  strweights=[]
  for i,fullt in enumerate(terms):
    newbit=fullt.funcDiff(field)
    newterms=[]
    newweights=[]
    for w,t in newbit:
      if not indexConsistency(t):
        print fullt.textRep()
        print fullt.longstr()
        print t.longstr()
        indexConsistency(t,True)
        raise Exception('Indices inconsistent in funcDiff term')
      t.fixDerivativeIndices()
      prevt=next((j for j,tm in enumerate(newterms) if t==tm),None)
      if prevt is None:
        newweights.append(w)
        newterms.append(t)
      else:
        newweights[prevt]=newweights[prevt]+w*newterms[prevt].compareSign(t)
    eom+=newterms
    weights+=newweights
    if coeffs is not None:
      strweights+=[coeffs[i]]*len(newterms)
  if coeffs is None:
    return eom,weights
  else:
    return eom,weights,strweights

def expandFieldStrength(term,fieldi):
  fieldDict={'B':AB}
  try: potentialClass=fieldDict[term[fieldi]._basename]
  except KeyError: return
  ######first term
  newterm1=term._DEFTcopy()
  potential=potentialClass(newterm1)
  fs=newterm1[fieldi]
  newterm1._fields[fieldi]=potential
  potential.Dindices=fs.Dindices
  potential.indices+=itertools.chain.from_iterable( [d[0],d[1]] for d in fs.Dindices )
  if fs.isConjugate:
    newlorL=Index('lorL',newterm1.smallestNewLabel('lorL'),False,True)
    potential.Dindices.append( (newlorL,fs.indices[0]) )
    potential.indices+=[newlorL,fs.indices[0]]
    potential.indices[1]=fs.indices[1]
    potential.indices[0].dummy=True
    newterm1._invariants.append( Invariant('epsl','\epsilon',[newlorL,potential.indices[0]]) )
  else:
    newlorR=Index('lorR',newterm1.smallestNewLabel('lorR'),False,True)
    potential.Dindices.append( (fs.indices[0],newlorR) )
    potential.indices+=[fs.indices[0],newlorR]
    potential.indices[0]=fs.indices[1]
    potential.indices[1].dummy=True
    newterm1._invariants.append( Invariant('epsl','\epsilon',[newlorR,potential.indices[1]]) )
  ######second term
  newterm2=term._DEFTcopy()
  potential=potentialClass(newterm2)
  fs=newterm2[fieldi]
  newterm2._fields[fieldi]=potential
  potential.Dindices=fs.Dindices
  potential.indices+=itertools.chain.from_iterable( [d[0],d[1]] for d in fs.Dindices )
  if fs.isConjugate:
    newlorL=Index('lorL',newterm2.smallestNewLabel('lorL'),False,True)
    potential.Dindices.append( (newlorL,fs.indices[1]) )
    potential.indices+=[newlorL,fs.indices[1]]
    potential.indices[1]=fs.indices[0]
    potential.indices[0].dummy=True
    newterm2._invariants.append( Invariant('epsl','\epsilon',[newlorL,potential.indices[0]]) )
  else:
    newlorR=Index('lorR',newterm2.smallestNewLabel('lorR'),False,True)
    potential.Dindices.append( (fs.indices[1],newlorR) )
    potential.indices+=[fs.indices[1],newlorR]
    potential.indices[0]=fs.indices[0]
    potential.indices[1].dummy=True
    newterm2._invariants.append( Invariant('epsl','\epsilon',[newlorR,potential.indices[1]]) )
  return [newterm1,newterm2]


#@profile
def makeSwitchDerivativeRelations(terms):
  mashDict=dict((mash(t),t) for t in terms)
  def makeSDRel(cts,trms,wghts):
    properterms=[]
    properweights=[]
    for i,c in enumerate(cts):
      if not indexConsistency(c):
        print c.shortstr()
        print c.textRep()
        indexConsistency(c,True)
        for f in c:
          print '************'
          for i in f.indices: print i,id(i)
          print '***'
          for s in f.symmetries: print s
        print c.longstr()
        raise Exception('ERROR: indices inconsistent from switchDer')
      if not c.checkSyms(): continue
      c.fixDerivativeIndices()
      try: newt=mashDict[mash(c)]
      except KeyError:
        #print 'WARNING: missing term in makeSDRel'
        continue
      properterms.append(newt)
      properweights.append(wghts[i] * newt.compareSign(c) )
    return Relation('SwitchDerivative',properterms,properweights)
  conjrels=[]
  for t in terms:
    #if any(mash(t)==mash(conjr.terms[-2]) or mash(t)==mash(conjr.terms[-1]) for conjr in conjrels): continue
    conjlists=switchDerivatives(t)
    for conjterms,conjweights in conjlists:
      conjrels.append( makeSDRel(conjterms,terms,conjweights) )
  return removeDupsFromRelList(conjrels)

def makeIBPRelations(terms):
  mashDict=dict((mash(t),t) for t in terms)
  def makeIBPRel(ibpts,trms):
    properterms=[]
    properweights=[]
    for i in ibpts:
      if not indexConsistency(i): raise Exception('ERROR: trying to make IBP rel out of term with inconsistent inds')
      #if not i.checkSyms(): continue
      try: newt=mashDict[mash(i)]
      except KeyError:
        #print 'WARNING: IBP rel generated absent term. Ignoring'
        continue
      properterms.append(newt)
      i.fixDerivativeIndices()
      properweights.append( newt.compareSign(i) )
    return (Relation('IBP',properterms,properweights),'')
  ibprels=[]
  for t in terms:
    if not t.checkSyms(): continue
    ibptermsList=IBPRels(t)
    for ibpterms in ibptermsList:
      #if len(ibpterms)==1 and len(t._fields)>1: continue
      ir,tx= makeIBPRel(ibpterms,terms) 
      ibprels.append( ir )
  return removeDupsFromRelList(ibprels)

def makeFuncDiffandEOMRelations(fdFields,dimFourTerms,dimFourCoeffs,dimNTerms):
  eomrels=[]
  for f in fdFields:
    eom,weights,strweights=funcDiff(dimFourTerms,f,dimFourCoeffs)
    #print '...done funcdiff',f.getName()
    if len(eom)==0: continue
    sweights=sympy.symbols(' '.join(strweights))
    for j in range(len(eom)): weights[j]*=sweights[j]
    ##generate bianchi idents for vectors
    if f.isVector():
      othereom=[next(eombit for eombit in eom if len(eombit._fields)==1 and not eombit[0].isConjugate),next(eombit for eombit in eom if len(eombit._fields)==1 and eombit[0].isConjugate)]
      otherweights=[weights[eom.index(eombit)] for eombit in eom]
      otherweights[1]=otherweights[1]*(-1)
    eomrels+=makeEOMRelations(f.__class__.__name__,eom,dimNTerms,weights)
    if f.isVector():
      eomrels+=makeEOMRelations(f.__class__.__name__+'p',othereom,dimNTerms,otherweights)
  for r in eomrels:
    if not all(t in dimNTerms for t in r.terms):
      raise Exception('ERROR: Generated term in EOM relation which is not in original list')
  return removeDupsFromRelList(eomrels)

def makeFDFields(fields):
  classAndInstDict=dict((f,f()) for f in fields)
  fdFields=[]
  for k,v in classAndInstDict.iteritems():
    if v.isFieldStrength():
      fdFields.append(v._vectorPotential())
    elif v.isSelfConjugate:
      fdFields.append(v)
    else:
      fdFields.append(v)
      conjField=k()
      conjField.conjugate()
      fdFields.append(conjField)
  return fdFields

'''
def makeEOMRelationsBetter(name,eom,terms,weights):
  words=False
  for t in terms: t._clear()#=None
  for t in eom: t._clear()#=None
  allInvariantss=[eombit._invariants for eombit in eom]
  for eombit in eom:
    if words: print eombit.shortstr()
    eombit._invariants=[inv for inv in eombit._invariants if all(i.dummy for i in inv.indices)]
  relevantTerms=[ [t for t in terms if t.contains(eombit)] for eombit in eom]
  for i,eombit in enumerate(eom):
    eombit._invariants=allInvariantss[i]
  if words:
    print '**************aargharaaghereghdgh'
    print name,[len(rt) for rt in relevantTerms]
  relevantTermsAndQuotients=[ [(t,t.quotientPlaceholderBetter(eombit)) for t in relevantTerms[i]] for i,eombit in enumerate(eom)]
  relevantTermsAndQuotients=[ [(t,q) for t,q in rti if len(q)>0] for rti in relevantTermsAndQuotients ]
  if words: print name,[len(rt) for rt in relevantTermsAndQuotients]
  eomrels=[]
  def compareQuotients(q1,q2):
    for w1,q in q1:
      if q==q2: return q.compareSign(q2)*w1
    return 0
  ##make list of all possible quotients
  allPossibleQuotients=removeDupsWMash(list(itertools.chain.from_iterable( itertools.chain.from_iterable([q for w,q in qu] for t,qu in rtqi) for rtqi in relevantTermsAndQuotients) ))
  if words: print [q.shortstr() for q in allPossibleQuotients]
  for quotient in allPossibleQuotients: quotient.fixDerivativeIndices()
  for quotient in allPossibleQuotients:
    product=relevantTermsAndQuotients[:]
    if words: print 'quot',quotient.shortstr()# for w,q in quotient]
    for i in range(len(eom)):
      newprodi=[]
      for t,qu in product[i]:
        for thew,theq in qu: theq.fixDerivativeIndices()
        wght=compareQuotients(qu,quotient)
        if wght !=0: newprodi.append((wght*weights[i],t))
      product[i]=newprodi
      if words: print [str(w)+' '+t.shortstr() for w,t in product[i]]
    product=[p for p in product if len(p)>0]
    for p in itertools.product(*product):
      wts,tms=zip(*p)
      eomrels.append(Relation('EOM'+name,list(tms),list(wts)))
  return eomrels
'''

def makeEOMRelations(name,eom,terms,weights):
  words=False
  for t in terms: t._clear()#=None
  for t in eom: t._clear()#=None
  allInvariantss=[eombit._invariants for eombit in eom]
  for eombit in eom:
    if words: print eombit.shortstr()
    eombit._invariants=[inv for inv in eombit._invariants if all(i.dummy for i in inv.indices)]
  relevantTerms=[ [t for t in terms if t.contains(eombit)] for eombit in eom]
  for i,eombit in enumerate(eom):
    eombit._invariants=allInvariantss[i]
  if words:
    print '**************aargharaaghereghdgh'
    print name,[len(rt) for rt in relevantTerms]
  relevantTermsAndQuotients=[ [(t,t.quotientPlaceholderBetter(eombit)) for t in relevantTerms[i]] for i,eombit in enumerate(eom)]
  relevantTermsAndQuotients=[ [(t,q) for t,q in rti if len(q)>0] for rti in relevantTermsAndQuotients ]
  if words: print name,[len(rt) for rt in relevantTermsAndQuotients]
  eomrels=[]
  def compareQuotients(q1,q2):
    for w1,q in q1:
      if q==q2: return q.compareSign(q2)*w1
    return 0
  ##make list of all possible quotients
  allPossibleQuotients=removeDupsWMash(list(itertools.chain.from_iterable( itertools.chain.from_iterable([q for w,q in qu] for t,qu in rtqi) for rtqi in relevantTermsAndQuotients) ))
  if words: print [q.shortstr() for q in allPossibleQuotients]
  for quotient in allPossibleQuotients: quotient.fixDerivativeIndices()
  for quotient in allPossibleQuotients:
    product=relevantTermsAndQuotients[:]
    if words: print 'quot',quotient.shortstr()# for w,q in quotient]
    for i in range(len(eom)):
      newprodi=[]
      for t,qu in product[i]:
        for thew,theq in qu: theq.fixDerivativeIndices()
        wght=compareQuotients(qu,quotient)
        if wght !=0: newprodi.append((wght*weights[i],t))
      product[i]=newprodi
      if words: print [str(w)+' '+t.shortstr() for w,t in product[i]]
    product=[p for p in product if len(p)>0]
    for p in itertools.product(*product):
      wts,tms=zip(*p)
      eomrels.append(Relation('EOM'+name,list(tms),list(wts)))
  ####ugly hack
  def badRelation2(rel):
    kintermfield=next(e for e in eom if len(e._fields)==1)[0]._basename
    gaugefieldgroup=next((group for group in DEFT.NAGaugeGroupsCoupling.keys() if DEFT.gaugeGroupsFields[group]._basename==kintermfield),None)
    if gaugefieldgroup is None: return False
    relterms=sorted(rel.terms,key=lambda t: 100*sum(len(f.Dindices) for f in t)+sum(1 if f._basename==kintermfield else 0 for f in t))
    return (sum(1 if 'epsu' in inv.name and inv.group()==gaugefieldgroup else 0 for inv in relterms[-1]._invariants)>0 and sum(1 if 'epsl' in inv.name and inv.group()==gaugefieldgroup else 0 for inv in relterms[-1]._invariants)>0)
  eomrels=[r for r in eomrels if not badRelation2(r)]
  return eomrels

def removeRedundantRels(rs):
  for r in rs:
    r._mash=set( mash(t) for t in r.terms )
  newrs=[]
  for r1 in rs:
    doAdd=True
    for r2 in newrs:
      if r1._mash==r2._mash:
        r1indextor2index=dict( (i,next(j for j in range(len(r2.terms)) if mash(r1.terms[i])==mash(r2.terms[j]))) for i in range(len(r1.terms)) )
        #print r1.type,r2.type
        if r1.weights is None and r2.weights is None:
          doAdd=False
          break
        weightratios=set( sympy.simplify( r1.weights[i]/r2.weights[j]*r1.terms[i].compareSign(r2.terms[j]) ) for i,j in r1indextor2index.iteritems() )
        if len(weightratios)==1: doAdd=False
        else:
          r2.terms=[r2.terms[r1indextor2index[i]] for i in range(len(r1.terms))]
          r2.weights=[r2.weights[r1indextor2index[i]] for i in range(len(r1.weights))]
          print 'r1',r1,'\nr2',r2
      if not doAdd: break
    if doAdd: newrs.append(r1)
  return newrs

def extraFierz(terms):
  def fieldContent(term):
    return set(f.getName() for f in term)
  def findIndex(dummy,t):
    return next(i for i in t.dummyIndices() if i==dummy)
  termGroups=[list(g) for k,g in itertools.groupby(terms,fieldContent)]
  rels=[]
  for grp,dim in DEFT.groupsDimFund.iteritems():
    if dim <=2: continue
    for tg in termGroups:
      #print tg[0].shortstr()
      fieldIndices=list(itertools.chain.from_iterable(f.indices for f in tg[0]))
      sunupper=[i for i in fieldIndices if i.group==grp and i.upper]
      sunlower=[i for i in fieldIndices if i.group==grp and not i.upper]
      #print sunupper,sunlower
      if len(sunupper)>dim and len(sunupper)==len(sunlower):
        newTerms=[]
        newWeights=[]
        for elindexperm in itertools.permutations(sunlower):
          indexpairs=zip(sunupper,elindexperm)
          newterm=tg[0]._DEFTcopy()
          newterm._clear()
          newterm._invariants=[inv for inv in newterm._invariants if inv.group()!=grp]
          for ip in indexpairs:
            newterm.addInvariant(Invariant('del','\delta',[findIndex(ip[0],newterm),findIndex(ip[1],newterm)]))
          newTerms.append(newterm)
          newWeights.append( newterm.isign() * (1 if newterm._cmpParity(sunlower,elindexperm) else -1) )
        rels.append( (newTerms,newWeights) )
  return rels

#@profile
def makeFierzRelations(terms):
  mashDict=dict((mash(t),t) for t in terms)
  #@profile
  def makeFierzRel(fierztrmsnwghts,trms,type):
    properterms=[]
    properweights=[]
    for ti,wi in zip(fierztrmsnwghts[0],fierztrmsnwghts[1]):
      if not indexConsistency(ti): raise Exception('making fierz rel with term of inconsistent inds')
      if not ti.checkSyms(): continue
      ti._clear()
      try: newt=mashDict[mash(ti)]
      except KeyError:
        #print 'WARNING: term generated by fierzing is absent from term list'
        continue
      properterms.append(newt)
      ti.fixDerivativeIndices()
      properweights.append(wi * ti.isign() * newt.compareSign(ti))
    return Relation(type,properterms,properweights)
  fierzrels=[]
  typesOfTwoIndexFierz=[fierzRelsEUEU,fierzRelsELEL,fierzRelsDelEU,fierzRelsDelEL,fierzRelsEUEL]
  typesOfAnyIndexFierz=[fierzRelsEUEL]
  n=0
  for t in terms:
    #if not t.checkSyms(): continue
    for g,dim in sorted(DEFT.groupsDimFund.iteritems(),key=lambda x:x[0]):
    #for g in ['lorL','lorR','su2']:
      if dim==2:
        for ft in typesOfTwoIndexFierz:
          fierzrels+=[makeFierzRel( rellist , terms,ft.__name__) for rellist in ft(t,g)]
    #for g in ['su3','sug']:
      else:
        for ft in typesOfAnyIndexFierz:
          fierzrels+=[makeFierzRel( rellist , terms,ft.__name__) for rellist in ft(t,g)]
    n+=1
    #if n>100: break
  fierzrels+=[makeFierzRel(rellist,terms,'fierzSUNExtra') for rellist in extraFierz(terms)]
  return removeDupsFromRelList(fierzrels)

def metric(term):
  def conjugateDerivatives(t):
    cd=0
    for f in t:
      if f.isConjugate: cd+=len(f.Dindices)
    return cd
  things=[
('sum(len(f.Dindices) for f in term)',5),
('len([i for i in term._invariants if \'eps\' in i.name])',4),
('sum(len(f.indices) for f in term)',4),
('(0 if term.conjugate()==term else 1)',1)]#,
#('len([f for f in term if f.isConjugate])',3)
#]
  met=0
  for thing in things:
    met = met << thing[1]
    met = met ^ (eval(thing[0]) % 2**thing[1] )
  return met

def removeDupsFierzSym(terms):
  contractedFierzTerms=removeDups(list(itertools.chain.from_iterable([fierzTrans(t) for t in terms])))
  return filter(lambda x:x.checkSyms(),contractedFierzTerms)

def fierzTrans(term):
  def findIndex(dummy,t):
    return next(i for i in t.dummyIndices() if i==dummy)
  epsu = [x for x in term._invariants if x.name=='epsu'] 
  epsl = [x for x in term._invariants if x.name=='epsl']
  terms=[term]
  for eu in epsu:
    newTerms=[]
    for el in epsl:
      if el not in terms[0]._invariants: continue
      if eu.indices[0].group==el.indices[0].group:
        for elindexperm in itertools.permutations(el.indices):
          indexpairs=zip(eu.indices,elindexperm)
          for t in terms:
            newterm=t._DEFTcopy()
            for ip in indexpairs:
              newterm._invariants.append( Invariant('del','\delta',[findIndex(ip[0],newterm),findIndex(ip[1],newterm)]) )
            newterm._clear()
            newterm._invariants.remove(eu)
            newterm._invariants.remove(el)
            newTerms.append(newterm)
        terms=newTerms
        break
  return terms

def makeCoeffs(terms):
  coeffs=[]
  for t in terms:
    name='c'+t.shortstr().replace(' ','')
    if name in coeffs:
      name=name+'p'
      while name in coeffs:
        name=name+'p'
    coeffs.append(name)
  return coeffs

def makeMatrix(tms,rls,mashes,subs=[]):
  rows=[]
  for r in rls:
    newrow=[0]*len(tms)
    for w,t in zip(r.weights,r.terms):
      newrow[mashes[mash(t)]]=w
    rows.append(newrow)
  if len(subs)>0:
    for i in range(len(rows)):
      for j in range(len(rows[i])):
        if isinstance(rows[i][j],int): continue
        rows[i][j]=rows[i][j].subs(subs)
  ##filter
  rows=[row for row in rows if len(set(row))>1 or set(row).pop()!=0]
  ###sort
  rows.sort(key = lambda r: next(i for i,e in enumerate(r) if e!=0))
  return rows

#@profile
def betterReducer(subsetTms,rls,termsToIgnore):
  if len(subsetTms)==0: return [],rls
  mashesToIgnore=set(mash(t) for t in termsToIgnore)
  relevantTerms=[t for t in subsetTms if mash(t) not in mashesToIgnore]
  relevantMashes=dict((mash(t),i) for i,t in enumerate(relevantTerms))
  workingRels=[r for r in rls if all(mash(t) in relevantMashes.keys() for t in r.terms)]
  if len(workingRels)==0: return [],rls
  notWorkingRels=[r for r in rls if any(mash(t) not in relevantMashes.keys() for t in r.terms)]
  rows=makeMatrix(relevantTerms,workingRels,relevantMashes)
  print 'INFO: matrix size',len(relevantTerms),len(rows)
  print 'INFO: row echelon reduction ...',
  symmat=sympy.Matrix(rows)
  rowechelon=symmat.rref(simplify=True)[0].tolist()
  print 'done'
  newrowechelon=[r for r in rowechelon if len(set(r))>1 or r[0]!=0]
  newrelations=[Relation('DEFT',[t for i,t in enumerate(relevantTerms) if row[i]!=0],[w for w in row if w!=0]) for row in newrowechelon]
  for r in notWorkingRels: r._mashes=set(mash(t) for t in r.terms)
  goneTerms=[]
  print 'INFO: eliminating terms in other relations ...',
  for rr in newrelations:
    termToGo=rr.terms[0]
    for r in notWorkingRels:
      if mash(termToGo) not in r._mashes: continue
      r.addbetter(rr,termToGo)
    goneTerms.append(termToGo)
  print 'done'
  remainingRels=newrelations+[r for r in notWorkingRels if len(r.terms)>0]
  return goneTerms,remainingRels


###definitions of splitters
def DsnamesEpssSplitter(term):
  fieldNames=['D'+str(len(f.Dindices))+f._name for f in term]
  fieldNames.sort()
  return ''.join(fieldNames)+' epsu'+str(sum(1 if 'epsu' in inv.name else 0 for inv in term._invariants))+' epsl'+str(sum(1 if 'epsl' in inv.name else 0 for inv in term._invariants))
def namesDsEpssSplitter(term):
  fieldNames=[f._name for f in term]
  fieldNames.sort()
  return ''.join(fieldNames)+' D'+str(sum(len(f.Dindices) for f in term))+' epsu'+str(sum(1 if 'epsu' in inv.name else 0 for inv in term._invariants))+' epsl'+str(sum(1 if 'epsl' in inv.name else 0 for inv in term._invariants))
def namesDsSplitter(term):
  fieldNames=[f._name for f in term]
  fieldNames.sort()
  return ''.join(fieldNames)+' D'+str(sum(len(f.Dindices) for f in term))
def matterSplitter(term):
  fieldNames=[f._name for f in term if f.isFermion() or f._basename=='H']
  fieldNames.sort()
  return ''.join(fieldNames)
def unitSplitter(term):
  return '1'
###end splitter defs


def divideByRels(subsetTms,rls,nSubsets):
  mashes=dict((mash(t),i) for i,t in enumerate(subsetTms))
  workingRels=[r for r in rls if all(mash(t) in mashes.keys() for t in r.terms)]
  relsByInds=[[mashes[mash(t)] for t in r.terms] for r in workingRels]
  subsetIndss=[]
  for n in range(nSubsets-1):
    subsetIndss.append(set(relsByInds[0]))
    while len(subsetIndss[-1])<len(subsetTms)/nSubsets and len(relsByInds)>0:
      #print 'running while loop'
      for j,ri in enumerate(relsByInds):
        #print 'running for loop',j
        if not any(i in subsetIndss[-1] for i in ri): continue
        subsetIndss[-1].update(ri)
        relsByInds.pop(j)
        break
      else: break
  subsetIndsLists=[[i for i in range(len(subsetTms)) if i in subsetInds] for subsetInds in subsetIndss]
  subsetIndsLists.append([i for i in range(len(subsetTms)) if all(i not in subsetInds for subsetInds in subsetIndss)]) 
  return [[subsetTms[i] for i in subsetIndsList] for subsetIndsList in subsetIndsLists]

#@profile
def betterTermEliminator(tms,rls):
  maxSubsetSize=20
  for t in tms: t._clear()
  if len(set(mash(t) for t in tms))<len(tms):
    raise Exception('ERROR: terms\' mashes are not unique')
  tms.sort(key=metric,reverse=True)
  fullMashes=dict((mash(t),i) for i,t in enumerate(tms))
  goneTms=[]
  for splitterFunc in [DsnamesEpssSplitter,namesDsEpssSplitter,namesDsSplitter,matterSplitter,unitSplitter]:
    print 'INFO: using splitterFunc',splitterFunc.__name__
    differentSplitterStrings=set(splitterFunc(t) for t in tms)
    for splitStr in differentSplitterStrings:
      subsetTms=[t for t in tms if splitterFunc(t)==splitStr]
      if len(subsetTms)>maxSubsetSize:
        subsets=divideByRels(subsetTms,rls,len(subsetTms)/maxSubsetSize)
        for sub in subsets:
          newGoneTms,rls=betterReducer(sub,rls,goneTms)
          goneTms=goneTms+newGoneTms
      newGoneTms,rls=betterReducer(subsetTms,rls,goneTms)
      goneTms=goneTms+newGoneTms
  return tms,rls

def categorise(terms):
  def termstr(term):
    fnames=list(collections.Counter(f._name for f in term).iteritems())
    fnames.sort(key=lambda x: x[0])
    tstr=' '.join(str(k)+('^'+str(v) if v>1 else '') for k,v in fnames)
    nD=sum(len(f.Dindices) for f in term)
    if nD>0: tstr+=(' D^'+str(nD) if nD>1 else ' D')
    return tstr
  termcount=list(collections.Counter(termstr(t) for t in terms).iteritems())
  termcount.sort(key=lambda x: x[0])
  for k,v in termcount:
    print v,k
  return
 

def analyse(tms,rls,DEFTFieldStr):
  fullMashes=dict((mash(t),i) for i,t in enumerate(tms))
  mask=set( min(fullMashes[mash(t)] for t in r.terms) for r in rls )
  makeCatalogue=True
  if makeCatalogue:
    print 'INFO: making mini catalogue'
    tex=''
    for i,t in enumerate(tms):
      tex+=r'\noindent'+'\n'+r'\bf{'+str(i)+'}~~~'+t.texRep()+r'  \texttt{'+t.__repr__()+'} \n\n'
    makePDFDoc(tex,DEFTFieldStr+'minicatalogue')
    tex=''
    for r in rls:
      tex+=r.texRep()+'\n\n'
      tex+='[ '+' '.join(str(fullMashes[mash(t)]) for t in r.terms)+' ]\n\n'
    makePDFDoc(tex,DEFTFieldStr+'minirelcatalogue')
  if len(mask)!=len(rls):
    raise Exception('ERROR: Number of reduced row echelon rows does not equal number of terms eliminated')
  print 'INFO: A viable basis follows'
  tex=''
  print 'INFO: Not in basis'
  for i,t in enumerate(tms):
    if i in mask:
      print t.shortstr(),i
  print '****'
  print 'INFO: In basis'
  for i,t in enumerate(tms):
    if i not in mask:
      print t.shortstr(),i
      tex+=t.texRep()+'\n\n'
  makeCatalogue=True
  if makeCatalogue: makePDFDoc(tex,DEFTFieldStr+'remTerms')
  kvs=list(collections.Counter([sum(len(f.Dindices) for f in t) for i,t in enumerate(tms) if i not in mask]).iteritems())
  kvs.sort(key=lambda x:x[0])
  for k,v in kvs:
    categorise([t for i,t in enumerate(tms) if i not in mask and sum(len(f.Dindices) for f in t)==k])
  print 'INFO: #derivatives, #terms'
  print '\n'.join('  '+str(k)+'\t'+str(v) for k,v in kvs)
  print '-------'
  print 'Tot\t'+str(sum(v for k,v in kvs))
  countBandL=False
  if countBandL:
    bvterms=[]
    lvterms=[]
    bminuslvterms=[]
    for i,t in enumerate(tms):
      if i not in mask:
        bn=t.uones['B']
        ln=t.uones['L']
        if ln!=frac(0,1): lvterms.append(t)
        if bn!=frac(0,1): bvterms.append(t)
        if (bn-ln)!=frac(0,1): bminuslvterms.append(t)
    print 'INFO: Total,B viol,L viol,B-L viol'
    print '     ',len(tms)-len(mask),len(bvterms),len(lvterms),len(bminuslvterms)
  return

def makeBasisListFromMonomials(tms,rls,bfileName):
  basisfile=open(bfileName,'w')
  fullMashes=dict((mash(t),i) for i,t in enumerate(tms))
  mask=[fullMashes[mash(r.terms[0])] for r in rls]
  for i in range(len(tms)):
    if i not in mask:
      basisfile.write(r'[(1,'+str(i)+r')]'+'\n')
  basisfile.close()

def makeBasisListFromBDict(bdict,bfileName):
  basisfile=open(bfileName,'w')
  for name,tms in bdict.iteritems():
    tmsstr=','.join("('{0}','{1}')".format(str(w),t.__repr__()) for w,t in zip(tms.weights,tms.terms))
    basisfile.write('[{0}]:{1}\n'.format(tmsstr,name))
  basisfile.close()

def readBasisList(tms,bfileName,fields=None):
  basisfile=open(bfileName,'r')
  basislines=basisfile.readlines()
  basisfile.close()
  basislines=[bl for bl in basislines if bl[0]!='#']
  basisDict=collections.OrderedDict()
  for bl in basislines:
    blsplit=bl.strip().split(':')
    wghts,inds=zip(*eval(blsplit[0]))
    wghts=[sympy.sympify(w) for w in wghts]
    tlist=[]
    wlist=[]
    for loc,i in enumerate(inds):
      if isinstance(i,int):
        tlist.append(tms[i])
        wlist.append(wghts[loc])
      else:
        pt=Term.unrepr(i,fields)
        tmatch=next(t for t in tms if pt==t)
        tlist.append(tmatch)
        wlist.append(wghts[loc]*tmatch.compareSign(pt))
    terms=Terms(tlist,wlist)
    if len(blsplit)==2:
      termname=blsplit[1]
    else:
      termname=terms.shortstr()
    while termname in basisDict.keys():
      termname+='p'
    basisDict[termname]=terms
  if len(basisDict)!=len(basislines):
    raise Exception('ERROR: more than one element of the basis list has the same name')
  return basisDict

def basisifyBetter(tms,oldrowechelon,basisDict):
  if len(tms)-len(oldrowechelon)-len(basisDict)!=0:
    raise Exception('ERROR: Wrong number of terms in basis list')
  mashes=dict((mash(t),i) for i,t in enumerate(tms))
  rows=[]
  for ts in basisDict.values():
    newrow=[0]*len(tms)
    for w,t in zip(ts.weights,ts.terms):
      newrow[mashes[mash(t)]]=w
    rows.append(newrow)
  basisMatrix=sympy.Matrix(rows)
  unbasisVectors=basisMatrix.nullspace(simplify=True)
  complementBasisMatrix=sympy.Matrix( [uv.T for uv in unbasisVectors] )
  ##each row of newExpressedInOld gives a new term expressed as a linear combination of the old ones (the tms)
  ##inverting newExpressedInOld as a sparsematrix object can result in nans and zoos
  #newExpressedInOld=SparseMatrix( (complementBasisMatrix,basisMatrix) )
  newExpressedInOld=sympy.Matrix( (complementBasisMatrix,basisMatrix) )
  #make unbasisDict
  unbasisDict=collections.OrderedDict()
  for ubv in unbasisVectors:
    ubvl=ubv.T.tolist()[0]
    #print 'ubvl',ubvl
    wghts,inds=[w for w in ubvl if w!=0],[i for i,w in enumerate(ubvl) if w!=0]
    terms=Terms([tms[i] for i in inds],wghts)
    termname=terms.shortstr()
    while termname in unbasisDict.keys():
      termname+='p'
    unbasisDict[termname]=terms
  #from the above comes the unbasisDict
  #oldrowechelonMat=SparseMatrix(makeMatrix(tms,oldrowechelon,mashes))
  oldrowechelonMat=sympy.Matrix(makeMatrix(tms,oldrowechelon,mashes))
  print 'INFO: inverting conversion matrix'
  newExpressedInOldInv=newExpressedInOld.inv()
  newExpressedInOldInv=sympy.Matrix(newExpressedInOldInv.tolist())
  oldrowechelonMat.tolist()
  oldrowechelonMat=sympy.Matrix(oldrowechelonMat.tolist())
  newrowechelonMat=oldrowechelonMat*newExpressedInOldInv
  print 'INFO: putting in RREF form'
  newrowechelonMat=sympy.Matrix(newrowechelonMat)
  newrowechelonMat=DEFTRREF(newrowechelonMat)
  if not sympy.eye(len(unbasisDict)).equals(newrowechelonMat[:,:len(unbasisDict)]):
    mask=[]
    for i in range(newrowechelonMat.shape[0]):
      mask.append( min(j for j in range(newrowechelonMat.shape[1]) if newrowechelonMat[i,j]!=0) )
    for j in range(newrowechelonMat.shape[1]):
      if j not in mask:
        print '********',j,'***********'
        if j<len(unbasisDict):
          print unbasisDict.values()[j].textRep()
        else:
          print basisDict.values()[j-len(unbasisDict)].textRep()
    #sympy.pprint(newrowechelonMat[:len(unbasisDict),:len(unbasisDict)])
    #sympy.pprint(newrowechelonMat[-len(basisDict):,-len(basisDict):])
    raise Exception('ERROR: The specified terms do not form a basis')
  return list(unbasisDict.values())+list(basisDict.values()),newrowechelonMat.tolist()

def changeBasis(tms,oldrowechelon,startBasisDict,endBasisDict):
  if len(tms)-len(oldrowechelon)-len(startBasisDict)!=0 or len(startBasisDict)-len(endBasisDict)!=0:
    raise Exception('ERROR: Wrong number of terms in basis list')
  mashes=dict((mash(t),i) for i,t in enumerate(tms))
  def basisToMon(basisDict):
    rows=[]
    for ts in basisDict.values():
      newrow=[0]*len(tms)
      for w,t in zip(ts.weights,ts.terms):
        newrow[mashes[mash(t)]]=w
      rows.append(newrow)
    return sympy.Matrix(rows)
  startBasisToMon=basisToMon(startBasisDict)
  endBasisToMon=basisToMon(endBasisDict)
  unbasisVectors=endBasisToMon.nullspace(simplify=True)
  complementBasisMatrix=sympy.Matrix( [uv.T for uv in unbasisVectors] )
  ##each row of newExpressedInOld gives a new term expressed as a linear combination of the old ones (the tms)
  endBasisAndUnbasisToMon=sympy.Matrix( (complementBasisMatrix,endBasisToMon) )
  print 'INFO: inverting endBasisDict',
  monToEndBasisAndUnbasis=endBasisAndUnbasisToMon.inv()
  print '...done'
  #make unbasisDict
  unbasisDict=collections.OrderedDict()
  for ubv in unbasisVectors:
    ubvl=ubv.T.tolist()[0]
    #print 'ubvl',ubvl
    wghts,inds=[w for w in ubvl if w!=0],[i for i,w in enumerate(ubvl) if w!=0]
    terms=Terms([tms[i] for i in inds],wghts)
    termname=terms.shortstr()
    while termname in unbasisDict.keys():
      termname+='p'
    unbasisDict[termname]=terms
  #from the above comes the unbasisDict
  #print unbasisDict
  oldrowechelonMat=sympy.Matrix(makeMatrix(tms,oldrowechelon,mashes))
  #print oldrowechelonMat.shape
  #print len(unbasisDict),len(tms)
  #divideAndConquer(oldrowechelonMat)
  newrowechelonMat=oldrowechelonMat*monToEndBasisAndUnbasis
  newrowechelonMat=DEFTRREF(newrowechelonMat)
  if not sympy.eye(len(unbasisDict)).equals(newrowechelonMat[:,:len(unbasisDict)]):
    raise Exception('ERROR: The specified terms do not form a basis')
  startBasisToEndBasis=startBasisToMon*monToEndBasisAndUnbasis[:,len(unbasisDict):]-startBasisToMon*monToEndBasisAndUnbasis[:,:len(unbasisDict)]*newrowechelonMat[:,len(unbasisDict):]
  #sympy.pprint(startBasisToEndBasis)
  #for i in range(startBasisToEndBasis.shape[0]):
  #  for j in range(startBasisToEndBasis.shape[1]):
  #    print i,j
  #    print startBasisToEndBasis[i,j]
  #    print startBasisToEndBasis[i,j].expand().simplify()
  return list(unbasisDict.values())+list(endBasisDict.values()),startBasisToEndBasis.tolist()

def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 50, fill = '*'):
  """
  Call in a loop to create terminal progress bar
  @params:
      iteration   - Required  : current iteration (Int)
      total       - Required  : total iterations (Int)
      prefix      - Optional  : prefix string (Str)
      suffix      - Optional  : suffix string (Str)
      decimals    - Optional  : positive number of decimals in percent complete (Int)
      length      - Optional  : character length of bar (Int)
      fill        - Optional  : bar fill character (Str)
  """
  percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
  filledLength = int(length * iteration // total)
  bar = fill * filledLength + '-' * (length - filledLength)
  print '\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix)+'\r',
  ## Print New Line on Complete
  #if iteration == total: print 'sfasdlfqa'

def getMask(mat):
  mask=[]
  for i in range(mat.shape[0]):
    mask.append( set(j for j in range(mat.shape[1]) if mat[i,j]!=0) )
  return mask

#@profile
def divideAndConquer(mat,heuristicColSets=[]):
  def superSimplify(mat):
    for i in range(mat.shape[0]):
      for j in range(mat.shape[1]):
        mat[i,j]=mat[i,j].expand().simplify()
    return
  def runSmallRREF(mat,mask,rows,cols):
    smallMat=sympy.Matrix([[mat[i,j] for j in cols] for i in rows])
    #sympy.pprint(smallMat)
    #print 'INFO: start smallRREF',smallMat.shape,#rows,cols,
    smallMat=smallMat.rref(simplify=True)[0]
    superSimplify(smallMat)
    smallMask=getMask(smallMat)
    #sympy.pprint(smallMat)
    #print 'tzsmallmat',testZero(smallMat)
    #print '... done'
    for si,i in enumerate(rows):
      for sj,j in enumerate(cols):
        mat[i,j]=smallMat[si,sj]
      mask[i]=set(cols[sj] for sj in smallMask[si])
    return
  def listSplit(lst,ln):
    nSublists=(len(lst)-1)//ln + 1
    subLists=[lst[ln*i:ln*(i+1)] for i in range(nSublists-1)]
    subLists.append( lst[ln*(nSublists-1):] )
    return subLists
  def runSmallRREFWrapper(mat,mask,rows,cols):
    rowMax=10
    if len(rows)>rowMax:
      subrows=listSplit(rows,rowMax)
      for sr in subrows:
        sc=sorted(list(set(itertools.chain.from_iterable(mask[i] for i in sr))))
        runSmallRREF(mat,mask,sr,sc)
    runSmallRREF(mat,mask,rows,cols)
  def trivialElims(mat):#,doneRows):
    mask=getMask(mat)
    doneRows=[i for i,m in enumerate(mask) if len(m)<=1]
    zeroCols=set(m[0] for m in mask if len(m)==1)
    for i in [i for i in range(mat.shape[0]) if i not in doneRows]:
      for j in zeroCols:
        mat[i,j]=0
    return doneRows
  def normaliseRows(mat,mask):
    for i,m in enumerate(mask):
      if len(m)==0: continue
      norm=mat[i,min(m)]
      for j in m: mat[i,j]=mat[i,j]/norm
    return
  def testZero(mat):
    for i in range(mat.shape[0]):
      print 'tz',i
      for j in range(mat.shape[1]):
        if (mat[i,j]!=0)!=(sympy.simplify(mat[i,j])!=0): return False
    return True
  #doneRows=trivialElims(mat)
  mask=getMask(mat)
  for hcscount,hcs in enumerate(heuristicColSets):
    printProgressBar(hcscount,len(heuristicColSets))
    #print 'hcs',hcs
    runSmallRREF(mat,mask,[i for i,m in enumerate(mask) if len(m)>1 and m.issubset(hcs)],sorted(list(set(itertools.chain.from_iterable(m for m in mask if len(m)>1 and m.issubset(hcs))))))
  if len(heuristicColSets)>0: print '\r'+'*'*19+'run heuristic reduction'+'*'*19
  count=collections.Counter(min(m) for m in mask if len(m)>1)
  commonMinCols=[k for k,v in count.iteritems() if v>1]
  while len(commonMinCols)>0:
    mc=min(commonMinCols)
    #print mc
    printProgressBar(mc,mat.shape[1])
    runSmallRREF(mat,mask,[i for i,m in enumerate(mask) if len(m)>1 and min(m)==mc],sorted(list(set(itertools.chain.from_iterable(m for m in mask if len(m)>1 and min(m)==mc)))))
    count=collections.Counter(min(m) for m in mask if len(m)>1)
    commonMinCols=[k for k,v in count.iteritems() if v>1]
  print '\r'+'*'*25+'in REF form'+'*'*25
  mincols=[min(m) for m in mask if len(m)>0]
  overlappingLeadingCols=[mc for mc in mincols if sum(1 if mc in m else 0 for m in mask)>1]
  while len(overlappingLeadingCols)>0:
    mc=min(overlappingLeadingCols)
    printProgressBar(mc,mat.shape[1])
    runSmallRREF(mat,mask,[i for i,m in enumerate(mask) if mc in m],sorted(list(set(itertools.chain.from_iterable(m for m in mask if mc in m)))))
    multiplicities=collections.Counter(itertools.chain.from_iterable(m for m in mask))
    overlappingLeadingCols=[min(m) for m in mask if len(m)>0 and multiplicities[min(m)]>1]
  print '\r'+'*'*25+'in RREF form'+'*'*24
  #newmask=getMask(mat)
  #print (mask==newmask)
  normaliseRows(mat,mask)
  return mask

def DEFTRREF(mat):
  def makeOrderedMat(mat,mask):
    newMat=sympy.zeros(*mat.shape) 
    rowPerm=range(mat.shape[0])
    rowPerm.sort(key=lambda x: min(mask[x]) if len(mask[x])>0 else mat.shape[1])
    for i,pi in enumerate(rowPerm):
      for j in mask[pi]:
        newMat[i,j]=mat[pi,j]
    return newMat
  mask=divideAndConquer(mat)
  return makeOrderedMat(mat,mask)
 
def matrixToRels(tms,mat,mask):
  rels=[]
  for i,m in enumerate(mask):
    if len(m)==0: continue
    rels.append( Relation('DEFT',[tms[j] for j in sorted(list(m))],[mat[i,j] for j in sorted(list(m))]) )
  return rels

def DEFTRREFRels(tms,rs):
  for t in tms: t._clear()
  if len(set(mash(t) for t in tms))<len(tms):
    raise Exception('ERROR: terms\' mashes are not unique')
  tms.sort(key=metric,reverse=True)
  mashes=dict((mash(t),i) for i,t in enumerate(tms))
  print 'INFO: making matrix'
  mat=sympy.Matrix(makeMatrix(tms,rs,mashes))
  print 'INFO: making heuristic term sets'
  heuristicColSets=[]
  for splitterFunc in [DsnamesEpssSplitter,namesDsEpssSplitter,namesDsSplitter,matterSplitter]:
    differentSplitterStrings=sorted(list(set(splitterFunc(t) for t in tms)))
    for splitStr in differentSplitterStrings[:-1]:
      cols=set(j for j,t in enumerate(tms) if splitterFunc(t)==splitStr)
      heuristicColSets.append(cols)
  print 'INFO: putting in RREF form'
  mask=divideAndConquer(mat,heuristicColSets)
  return tms,matrixToRels(tms,mat,mask)

def removeDupsFromRelList(relations):
  for r in relations: r.removeDups()
  for r in relations: r.removeUnsyms()
  return [r for r in relations if len(r.terms)>0]

