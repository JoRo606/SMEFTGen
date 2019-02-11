#!/usr/bin/env python

import sys
sys.path.insert(0,'..')

import itertools
import os
import sympy
import random
import DEFT.operators as op
import DEFT
from DEFT import frac,indexConsistency,mash

import unittest

class RankTestCase(unittest.TestCase):

  dimCutoff=6

  subsstr=[['cBB', '-1/32'], ['cBcBc', '1/32'], ['cFF', '-1/32'], ['cFcFc', '1/32'], ['cWW', '-1/16'], ['cWcWc', '1/16'], ['cGG', '-1/16'], ['cGcGc', '1/16'], ['cDeRceR', '0'], ['ceRcDeR', 'I'], ['cDlLclL', '0'], ['clLcDlL', 'I'], ['cDqLcqL', '0'], ['cqLcDqL', 'I'], ['cDuRcuR', '0'], ['cuRcDuR', 'I'], ['cDdRcdR', '0'], ['cdRcDdR', 'I'], ['cDDHcH', '0'], ['cDHcDH', '-1/2'], ['cHcDDH', '0'], ['cDDScS', '0'], ['cDScDS', '-1/2'], ['cScDDS', '0'], ['cDDSS', '0'], ['cDSDS', '-1/4'], ['cSDDS', '0'],['cDqL4cqL4', '0'], ['cqL4cDqL4', 'I'], ['cDuR4cuR4', '0'], ['cuR4cDuR4', 'I'], ['cDdR4cdR4', '0'], ['cdR4cDdR4', 'I'], ['cG4G4', '-1/32'], ['cG4cG4c', '1/32']]
  subs=[(sympy.sympify(s[0]),sympy.sympify(s[1])) for s in subsstr]
  reprsubs=[('e[ab] e[cd] B[da] B[cb]','-1/32'),('e[ab] e[cd] B![db] B![ca]','1/32'),('e[ab] e[cd] d[ef] d[gh] W[daeh] W[cbgf]','-1/16'),('e[ab] e[cd] d[ef] d[gh] W![dbhe] W![cafg]','1/16'),('e[ab] e[cd] d[ef] d[gh] G[daeh] G[cbgf]','-1/16'),('e[ab] e[cd] d[ef] d[gh] G![dbhe] G![cafg]','1/16'),('d[ab] d[cd] D[bd]eR![a] eR[c]','0'),('d[ab] d[cd] eR![a] D[bd]eR[c]','I'),('d[ab] d[cd] d[ef] D[bd]lL![cf] lL[ae]','0'),('d[ab] d[cd] d[ef] lL![cf] D[bd]lL[ae]','I'),('d[ab] d[cd] d[ef] d[gh] D[bd]qL![cfh] qL[aeg]','0'),('d[ab] d[cd] d[ef] d[gh] qL![cfh] D[bd]qL[aeg]','I'),('d[ab] d[cd] d[ef] D[bd]uR![af] uR[ce]','0'),('d[ab] d[cd] d[ef] uR![af] D[bd]uR[ce]','I'),('d[ab] d[cd] d[ef] D[bd]dR![af] dR[ce]','0'),('d[ab] d[cd] d[ef] dR![af] D[bd]dR[ce]','I'),('e[ab] e[cd] d[ef] D[bd]D[ac]H![f] H[e]','0'),('e[ab] e[cd] d[ef] D[bd]H![f] D[ac]H[e]','-1/2'),('e[ab] e[cd] d[ef] H![f] D[bd]D[ac]H[e]','0')]
  reprsubs=[('e[ab] e[cd] B[da] B[cb]','-1/32'),('e[ab] e[cd] B![db] B![ca]','1/32'),('e[ab] e[cd] d[ef] d[gh] W[daeh] W[cbgf]','-1/16'),('e[ab] e[cd] d[ef] d[gh] W![dbhe] W![cafg]','1/16'),('e[ab] e[cd] d[ef] d[gh] G[daeh] G[cbgf]','-1/16'),('e[ab] e[cd] d[ef] d[gh] G![dbhe] G![cafg]','1/16'),('d[ab] d[cd] D[bd]eR![a] eR[c]','0'),('d[ab] d[cd] eR![a] D[bd]eR[c]','I'),('d[ab] d[cd] d[ef] D[bd]lL![cf] lL[ae]','0'),('d[ab] d[cd] d[ef] lL![cf] D[bd]lL[ae]','I'),('d[ab] d[cd] d[ef] d[gh] D[bd]qL![cfh] qL[aeg]','0'),('d[ab] d[cd] d[ef] d[gh] qL![cfh] D[bd]qL[aeg]','I'),('d[ab] d[cd] d[ef] D[bd]uR![af] uR[ce]','0'),('d[ab] d[cd] d[ef] uR![af] D[bd]uR[ce]','I'),('d[ab] d[cd] d[ef] D[bd]dR![af] dR[ce]','0'),('d[ab] d[cd] d[ef] dR![af] D[bd]dR[ce]','I'),('e[ab] e[cd] d[ef] D[bd]D[ac]H![f] H[e]','0'),('e[ab] e[cd] d[ef] D[bd]H![f] D[ac]H[e]','1/2'),('e[ab] e[cd] d[ef] H![f] D[bd]D[ac]H[e]','0'),('e[ab] e[cd] d[ef] d[gh] G4[daeh] G4[cbgf]','-1/32'),('e[ab] e[cd] d[ef] d[gh] G4![dbhe] G4![cafg]','1/32'),('d[ab] d[cd] d[ef] d[gh] D[bd]qL4![cfh] qL4[aeg]','0'),('d[ab] d[cd] d[ef] d[gh] qL4![cfh] D[bd]qL4[aeg]','I'),('d[ab] d[cd] d[ef] D[bd]uR4![af] uR4[ce]','0'),('d[ab] d[cd] d[ef] uR4![af] D[bd]uR4[ce]','I'),('d[ab] d[cd] d[ef] D[bd]dR4![af] dR4[ce]','0'),('d[ab] d[cd] d[ef] dR4![af] D[bd]dR4[ce]','I'),('e[ab] e[cd] D[ac]S[] D[bd]S[]','1/4'),('e[ab] e[cd] S[] D[ac]D[bd]S[]','0')]
  listOfFields=[op.H,op.B,op.W,op.LL,op.ER,op.QL4,op.UR4,op.DR4,op.G4,op.QL,op.UR,op.DR,op.G,op.phi]
  termSubs=[(DEFT.Term.unrepr(r,listOfFields),sympy.sympify(w)) for r,w in reprsubs]
  ###,('e[ab] d[cd] H![d] eR![b] lL[ac]','cHceRclL'),('e[ab] d[cd] H[c] eR[b] lL![ad]','cHeRlLc'),('e[ab] e[cd] d[ef] H![d] qL![bcf] uR[ae]','cHcqLcuR'),('e[ab] e[cd] d[ef] H[d] qL[bce] uR![af]','cHqLuRc'),('e[ab] d[cd] d[ef] H[c] qL![bdf] dR[ae]','cHqLcdR'),('e[ab] d[cd] d[ef] H![d] qL[bce] dR![af]','cHcqLdRc'),('d[ab] d[cd] H![b] H![d] H[a] H[c]','cHcHcHH')]
  subPrime=False
  subPrime=True
  subZero=False

  def relrepr(self,rel):
    #return frozenset( (str(w),mash(t)) for w,t in zip(rel.weights,rel.terms)  )
    return [ (str(w),repr(t)) for w,t in zip(rel.weights,rel.terms)  ]
  def relmash(self,rel):
    #return frozenset( (str(w),mash(t)) for w,t in zip(rel.weights,rel.terms)  )
    return [ (str(w),mash(t)) for w,t in zip(rel.weights,rel.terms)  ]
  def makeTmsRelsRREF(self,DEFTFields,n):
    #generate terms
    dimFourTerms=op.generateAllTerms(DEFTFields,op.dimLEQ(4))
    self.assertTrue(all(indexConsistency(t) for t in dimFourTerms),'indices of dim4 terms inconsistent')
    dimFourTerms=[t for t in dimFourTerms if t.dimension==frac(4,1)]
    dimNTerms=op.generateAllTerms(DEFTFields,op.dimLEQ(n))
    self.assertTrue(all(indexConsistency(t) for t in dimNTerms),'indices of dimN terms inconsistent')
    dimNTerms=[t for t in dimNTerms if t.dimension==frac(n,1)]

    self.assertFalse(any(t.fixDerivativeIndices() for t in dimFourTerms),'indices of dim4 terms needed rearranging')
    self.assertFalse(any(t.fixDerivativeIndices() for t in dimNTerms),'indices of dimN terms needed rearranging')
    
    #make relations
    EOMCoeffSubs=[]
    if n>4:
      dimFourTermsSimple=op.removeDupsFierzSym(dimFourTerms)
      dimFourCoeffsSimple=op.makeCoeffs(dimFourTermsSimple)
      fdFields=op.makeFDFields(DEFTFields)
      eomrels=op.makeFuncDiffandEOMRelations(fdFields,dimFourTermsSimple,dimFourCoeffsSimple,dimNTerms)
      for t,w in zip(dimFourTermsSimple,dimFourCoeffsSimple):
        matcht,sub=next(((t2,w2) for t2,w2 in self.termSubs if t==t2),(None,None))
        if matcht is not None: EOMCoeffSubs.append( (sympy.sympify(w),sub*matcht.compareSign(t)) )
    else: eomrels=[]
    frs=op.makeFierzRelations(dimNTerms)
    sds=op.makeSwitchDerivativeRelations(dimNTerms)
    ibps=op.makeIBPRelations(dimNTerms)

    relations=frs+ibps+sds+eomrels
   
    ##check sign of terms in relations 
    for r in relations:
      self.assertTrue(all( next(tm for tm in dimNTerms if mash(tm)==mash(t)).compareSign(t)==1 for t in r.terms ),'sign of terms in relations differ from sign of terms in terms list')
    ##check weights of relations are sympy objects
    self.assertTrue( all(all(isinstance(w,sympy.Basic) for w in r.weights) for r in relations) )

    #substitute remaining coeffs with prime numbers or zeroes to speed up row reduction
    if self.subPrime or self.subZero:
      allCoeffs=set(itertools.chain.from_iterable(itertools.chain.from_iterable(w.free_symbols for w in r.weights) for r in relations))
      subbedCoeffs=[s[0] for s in self.subs]
      remainingCoeffs=list(allCoeffs-set(subbedCoeffs))
      if self.subPrime: listOfInts=[23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163]
      else: listOfInts=[0]*30
      if len(remainingCoeffs)>len(listOfInts): raise Exception('I don\'t know enough primes!')
      random.shuffle(listOfInts)
      extraSubs=[(co,sympy.Integer(pr)) for co,pr in zip(remainingCoeffs,listOfInts[:len(remainingCoeffs)])]
    else: extraSubs=[]

    #perform subs
    for r in relations:
      for i,w in enumerate(r.weights): r.weights[i]=sympy.simplify(w).subs(EOMCoeffSubs+extraSubs)
    #perform row reduction
    tms,rrefrls=op.DEFTRREFRels(dimNTerms,relations)
    return tms,relations,rrefrls

  def test_HLL_numbers(self):
    fields=[op.H,op.LL]
    answerDict={1:0,2:1,3:0,4:3,5:2,6:6,7:6,8:18}
    #answerDict={6:6,7:6,8:18}
    for dim,answer in answerDict.iteritems():
      if dim>self.dimCutoff: continue
      tms,rels,rrefrels=self.makeTmsRelsRREF(fields,dim)
      self.assertEqual(len(tms)-len(rrefrels),answer,'There should be '+str(answer)+' operators at dim '+str(dim)+', not '+str(len(tms)-len(rrefrels))+' operators')
 
  def test_phiHLL_numbers(self):
    fields=[op.phi,op.H,op.LL]
    answerDict={1:1,2:2,3:2,4:6,5:5,6:12,7:21,8:48}
    #answerDict={6:12,7:21,8:48}
    for dim,answer in answerDict.iteritems():
      if dim>self.dimCutoff: continue
      tms,rels,rrefrels=self.makeTmsRelsRREF(fields,dim)
      self.assertEqual(len(tms)-len(rrefrels),answer,'There should be '+str(answer)+' operators at dim '+str(dim)+', not '+str(len(tms)-len(rrefrels))+' operators')
  def test_BER_numbers(self):
    fields=[op.B,op.ER]
    answerDict={1:0,2:0,3:0,4:3,5:0,6:1,7:0,8:5}
    #answerDict={6:1,7:0,8:5}
    for dim,answer in answerDict.iteritems():
      if dim>self.dimCutoff: continue
      tms,rels,rrefrels=self.makeTmsRelsRREF(fields,dim)
      self.assertEqual(len(tms)-len(rrefrels),answer,'There should be '+str(answer)+' operators at dim '+str(dim)+', not '+str(len(tms)-len(rrefrels))+' operators')
  def test_HBWLLER_numbers(self):
    fields=[op.H,op.B,op.W,op.LL,op.ER]
    answerDict={1:0,2:1,3:0,4:10,5:2,6:23,7:12,8:179}
    for dim,answer in answerDict.iteritems():
      if dim>self.dimCutoff: continue
      tms,rels,rrefrels=self.makeTmsRelsRREF(fields,dim)
      self.assertEqual(len(tms)-len(rrefrels),answer,'There should be '+str(answer)+' operators at dim '+str(dim)+', not '+str(len(tms)-len(rrefrels))+' operators')
  def test_SM_numbers(self):
    fields=op.listOfSMFields
    answerDict={1:0,2:1,3:0,4:19,5:2,6:84,7:30,8:993}
    #answerDict={6:84,7:30,8:993}
    for dim,answer in answerDict.iteritems():
      if dim>self.dimCutoff: continue
      tms,rels,rrefrels=self.makeTmsRelsRREF(fields,dim)
      self.assertEqual(len(tms)-len(rrefrels),answer,'There should be '+str(answer)+' operators at dim '+str(dim)+', not '+str(len(tms)-len(rrefrels))+' operators')
 
  def test_SU4_numbers(self):
    fields=[op.H,op.B,op.W,op.LL,op.ER,op.QL4,op.UR4,op.DR4]
    answerDict={1:0,2:1,3:0,4:17,5:2,6:68}
    for dim,answer in answerDict.iteritems():
      if dim>self.dimCutoff: continue
      tms,rels,rrefrels=self.makeTmsRelsRREF(fields,dim)
      self.assertEqual(len(tms)-len(rrefrels),answer,'There should be '+str(answer)+' operators at dim '+str(dim)+', not '+str(len(tms)-len(rrefrels))+' operators')
  
  def test_SU4Gauge_numbers(self):
    DEFT.NAGaugeGroupsCoupling['su4']=sympy.symbols('g4')
    DEFT.gaugeGroupsFields['su4']=op.G4
    fields=[op.H,op.B,op.W,op.LL,op.ER,op.QL4,op.UR4,op.DR4,op.G4]
    answerDict={1:0,2:1,3:0,4:19,5:2,6:76}
    for dim,answer in answerDict.iteritems():
      if dim>self.dimCutoff: continue
      tms,rels,rrefrels=self.makeTmsRelsRREF(fields,dim)
      self.assertEqual(len(tms)-len(rrefrels),answer,'There should be '+str(answer)+' operators at dim '+str(dim)+', not '+str(len(tms)-len(rrefrels))+' operators')

if __name__ == '__main__':
  unittest.main()
