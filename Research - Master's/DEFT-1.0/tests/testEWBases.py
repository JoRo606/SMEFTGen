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
  subPrime=False
  #subPrime=True
  subZero=False

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
  def test_changeBasisEW(self):
    SILHBasisFileText='''[('-1/4','e[ab] e[cd] d[ef] d[gh] H![f] H![h] D[ad]H[e] D[bc]H[g]'),('-1/4','e[ab] e[cd] d[ef] d[gh] D[ad]H![f] D[bc]H![h] H[e] H[g]'),('-1/2','e[ab] e[cd] d[ef] d[gh] D[ad]H![h] H![f] D[bc]H[e] H[g]')]:H
[('-1/4','e[ab] e[cd] d[ef] d[gh] H![f] H![h] D[ad]H[e] D[bc]H[g]'),('-1/4','e[ab] e[cd] d[ef] d[gh] D[ad]H![f] D[bc]H![h] H[e] H[g]'),('1/2','e[ab] e[cd] d[ef] d[gh] D[ad]H![h] H![f] D[bc]H[e] H[g]')]:T
[('-1/(4*g**2)','e[ab] e[cd] e[ef] e[gh] d[ij] d[kl] D[fh]W[ebil] D[cg]W[dakj]'),('-1/(4*g**2)','e[ab] e[cd] e[ef] e[gh] d[ij] d[kl] D[ag]W![hfli] D[bc]W![dejk]'),('-1/(2*g**2)','e[ab] e[cd] e[ef] e[gh] d[ij] d[kl] D[de]W![hfli] D[bg]W[ackj]')]:2W
[('-I/(2*g)','e[ab] e[cd] e[ef] d[gh] d[ij] H![j] D[df]H[g] D[be]W[acih]'),('I/(2*g)','e[ab] e[cd] e[ef] d[gh] d[ij] D[df]H![j] H[g] D[be]W[acih]'),('-I/(2*g)','e[ab] e[cd] e[ef] d[gh] d[ij] H![j] D[af]H[g] D[bc]W![dehi]'),('I/(2*g)','e[ab] e[cd] e[ef] d[gh] d[ij] D[af]H![j] H[g] D[bc]W![dehi]')]:W
[('-I*gp/(4*g**2)','e[ab] e[cd] e[ef] d[gh] H![h] D[df]H[g] D[be]B[ac]'),('I*gp/(4*g**2)','e[ab] e[cd] e[ef] d[gh] D[df]H![h] H[g] D[be]B[ac]'),('-I*gp/(4*g**2)','e[ab] e[cd] e[ef] d[gh] H![h] D[af]H[g] D[bc]B![de]'),('I*gp/(4*g**2)','e[ab] e[cd] e[ef] d[gh] D[af]H![h] H[g] D[bc]B![de]')]:B
[('I/g','e[ab] e[cd] e[ef] d[gh] d[ij] D[df]H![j] D[ae]H[g] W[bcih]'),('I/g','e[ab] e[cd] e[ef] d[gh] d[ij] D[af]H![j] D[bd]H[g] W![cehi]')]:HW
[('-1/g','e[ab] e[cd] e[ef] d[gh] d[ij] D[df]H![j] D[ae]H[g] W[bcih]'),('1/g','e[ab] e[cd] e[ef] d[gh] d[ij] D[af]H![j] D[bd]H[g] W![cehi]')]:HW~
[('I*gp/(2*g**2)','e[ab] e[cd] e[ef] d[gh] D[df]H![h] D[ae]H[g] B[bc]'),('I*gp/(2*g**2)','e[ab] e[cd] e[ef] d[gh] D[af]H![h] D[bd]H[g] B![ce]')]:HB
[('-gp/(2*g**2)','e[ab] e[cd] e[ef] d[gh] D[df]H![h] D[ae]H[g] B[bc]'),('gp/(2*g**2)','e[ab] e[cd] e[ef] d[gh] D[af]H![h] D[bd]H[g] B![ce]')]:HB~
[('-I*g/2','e[ab] e[cd] e[ef] d[gh] d[ij] d[kl] W[ebgl] W[daih] W[fckj]'),('I*g/2','e[ab] e[cd] e[ef] d[gh] d[ij] d[kl] W![fdlg] W![bchi] W![eajk]')]:3W
[('g/2','e[ab] e[cd] e[ef] d[gh] d[ij] d[kl] W[ebgl] W[daih] W[fckj]'),('g/2','e[ab] e[cd] e[ef] d[gh] d[ij] d[kl] W![fdlg] W![bchi] W![eajk]')]:3W~
[('gp**2/(2*g**2)','e[ab] e[cd] d[ef] H![f] H[e] B[cb] B[ad]'),('-gp**2/(2*g**2)','e[ab] e[cd] d[ef] H![f] H[e] B![db] B![ac]')]:gamma
[('I*gp**2/(2*g**2)','e[ab] e[cd] d[ef] H![f] H[e] B[cb] B[ad]'),('I*gp**2/(2*g**2)','e[ab] e[cd] d[ef] H![f] H[e] B![db] B![ac]')]:gamma~
[('cHcHcHH','d[ab] d[cd] d[ef] H![b] H![d] H![f] H[a] H[c] H[e]')]:6'''
    WarsawBasisFileText='''[('-1/2','e[ab] e[cd] d[ef] d[gh] H![f] H![h] D[ad]D[bc]H[e] H[g]'),('-1/2','e[ab] e[cd] d[ef] d[gh] D[ad]D[bc]H![f] H![h] H[e] H[g]'),('-1','e[ab] e[cd] d[ef] d[gh] D[ad]H![f] H![h] D[bc]H[e] H[g]')]:phibox
[('-1/2','e[ab] e[cd] d[ef] d[gh] D[ad]H![h] H![f] D[bc]H[e] H[g]')]:phiD
[('I','d[ab] d[cd] d[ef] d[gh] H![f] D[bd]H[e] lL![ch] lL[ag]'),('-I','d[ab] d[cd] d[ef] d[gh] D[bd]H![f] H[e] lL![ch] lL[ag]')]:phil1
[('2*I','d[ab] d[cd] d[ef] d[gh] H![h] D[bd]H[e] lL![cf] lL[ag]'),('-2*I','d[ab] d[cd] d[ef] d[gh] D[bd]H![h] H[e] lL![cf] lL[ag]'),('-I','d[ab] d[cd] d[ef] d[gh] H![f] D[bd]H[e] lL![ch] lL[ag]'),('I','d[ab] d[cd] d[ef] d[gh] D[bd]H![f] H[e] lL![ch] lL[ag]')]:phil3
[('-I/8','e[ab] e[cd] e[ef] d[gh] d[ij] d[kl] W[ebgl] W[daih] W[fckj]'),('I/8','e[ab] e[cd] e[ef] d[gh] d[ij] d[kl] W![fdlg] W![bchi] W![eajk]')]:W
[('1/8','e[ab] e[cd] e[ef] d[gh] d[ij] d[kl] W[ebgl] W[daih] W[fckj]'),('1/8','e[ab] e[cd] e[ef] d[gh] d[ij] d[kl] W![fdlg] W![bchi] W![eajk]')]:W~
[('1/4','e[ab] e[cd] d[ef] d[gh] d[ij] H![f] H[e] W[cbgj] W[adih]'),('-1/4','e[ab] e[cd] d[ef] d[gh] d[ij] H![f] H[e] W![dbjg] W![achi]')]:phiW
[('I/4','e[ab] e[cd] d[ef] d[gh] d[ij] H![f] H[e] W[cbgj] W[adih]'),('I/4','e[ab] e[cd] d[ef] d[gh] d[ij] H![f] H[e] W![dbjg] W![achi]')]:phiW~
[('1/4','e[ab] e[cd] d[ef] d[gh] H![h] H[e] B[cb] W[adgf]'),('-1/4','e[ab] e[cd] d[ef] d[gh] H![h] H[e] B![db] W![acfg]')]:phiWB
[('I/4','e[ab] e[cd] d[ef] d[gh] H![h] H[e] B[cb] W[adgf]'),('I/4','e[ab] e[cd] d[ef] d[gh] H![h] H[e] B![db] W![acfg]')]:phiWB~
[('2','e[ab] e[cd] d[ef] d[gh] lL![df] lL![ch] lL[ae] lL[bg]')]:ll
[('1/8','e[ab] e[cd] d[ef] H![f] H[e] B[cb] B[ad]'),('-1/8','e[ab] e[cd] d[ef] H![f] H[e] B![db] B![ac]')]:phiB
[('I/8','e[ab] e[cd] d[ef] H![f] H[e] B[cb] B[ad]'),('I/8','e[ab] e[cd] d[ef] H![f] H[e] B![db] B![ac]')]:phiB~
[('1','d[ab] d[cd] d[ef] H![b] H![d] H![f] H[a] H[c] H[e]')]:phi'''
    targetChangeBasisMatInMathematicaFormat='''{{-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1/2, -1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -2, 0, 0, -g^2/gp^2, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {6, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8}, {0, 0, 0, 0, 0, 0, 0, 0, 0, (1/4)/g, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1/4)/g, 0, 0, 0}, {0, 0, 0, 1, -1, -1, 0, 1, 0, 0, 0, 1/4, 0, 0}, {0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 1/4, 0}, {0, 0, 0, 0, g/gp, 0, 0, -g/gp, 0, 0, 0, -1/4*g/gp, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, -g/gp, 0, 0, 0, -1/4*g/gp, 0}, {-6, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1/4)*g^2/gp^2, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1/4)*g^2/gp^2, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, cHcHcHH^(-1)}}'''
    targetChangeBasisMatInSympyFormat='''[[-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1/2, -1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -2, 0, 0, -g**2/gp**2, 0, 0, 0, 0, 0, 0, 0, 0, 0], [6, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(4*g), 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(4*g), 0, 0, 0], [0, 0, 0, 1, -1, -1, 0, 1, 0, 0, 0, 1/4, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 1/4, 0], [0, 0, 0, 0, g/gp, 0, 0, -g/gp, 0, 0, 0, -g/(4*gp), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -g/gp, 0, 0, 0, -g/(4*gp), 0], [-6, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g**2/(4*gp**2), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g**2/(4*gp**2), 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/cHcHcHH]]'''
    fields=[op.H,op.B,op.W,op.LL]
    warsawfilename='testDEFTWarsawBasis.txt'
    silhfilename='testDEFTSILHBasis.txt'
    with open(warsawfilename,'w') as bfile:
      bfile.write(WarsawBasisFileText)
    with open(silhfilename,'w') as bfile:
      bfile.write(SILHBasisFileText)
    targetChangeBasisMat=sympy.sympify(targetChangeBasisMatInSympyFormat)
    ts,rs,rrefrls=self.makeTmsRelsRREF(fields,6)
    startBasisDict=op.readBasisList(ts,warsawfilename,fields)
    endBasisDict=op.readBasisList(ts,silhfilename,fields)
    endtms,startToEndMat=op.changeBasis(ts,rrefrls,startBasisDict,endBasisDict)
    sympy.pprint(startToEndMat)
    os.remove(warsawfilename)
    os.remove(silhfilename)
    ####pretty print basis conversion
    startTerms=startBasisDict.keys()
    endTerms=endBasisDict.keys()
    for j,et in enumerate(endTerms):
      print et,'=',' + '.join([str(startToEndMat[i][j])+'*'+st for i,st in enumerate(startTerms) if startToEndMat[i][j]!=0])
    ####end pretty print
    self.assertEquals(startToEndMat,targetChangeBasisMat)
  
if __name__ == '__main__':
  unittest.main()
