#!/usr/bin/env python

import sys
sys.path.insert(0,'..')

import os
import itertools
import sympy
import random
import DEFT.operators as op
from DEFT import indexConsistency,Term,frac

import unittest

class TermReprTestCase(unittest.TestCase):

  basisfiletext='''[('-1/4','e[ab] e[cd] d[ef] d[gh] H![f] H![h] D[ad]H[e] D[bc]H[g]'),('-1/4','e[ab] e[cd] d[ef] d[gh] D[ad]H![f] D[bc]H![h] H[e] H[g]'),('-1/2','e[ab] e[cd] d[ef] d[gh] D[ad]H![h] H![f] D[bc]H[e] H[g]')]:H
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
[('cHcHcHH','d[ab] d[cd] d[ef] H![b] H![d] H![f] H[a] H[c] H[e]')]:6
'''

  def test_SM_Dim6_Repr(self):
    fields=op.listOfSMFields
    dim6Terms=op.generateAllTerms(fields,op.dimLEQ(6))
    self.assertTrue(all(indexConsistency(t) for t in dim6Terms),'indices of terms inconsistent')
    dim6TermsREPR=[t.__repr__() for t in dim6Terms]
    dim6TermsRecons=[Term.unrepr(tREPR,fields) for tREPR in dim6TermsREPR]
    for t,trecons in zip(dim6Terms,dim6TermsRecons):
      self.assertEqual(t,trecons)
      self.assertEqual(t.compareSign(trecons),1)

  def test_basis_file_IO(self):
    fields=[op.H,op.B,op.W,op.LL]
    bfilename='testReprBasis.txt'
    bfilename2='testReprBasis2.txt'
    with open(bfilename,'w') as bfile:
      bfile.write(self.basisfiletext)
    dim6Terms=op.generateAllTerms(fields,op.dimLEQ(6))
    dim6Terms=[t for t in dim6Terms if t.dimension==frac(6,1)]
    bDict=op.readBasisList(dim6Terms,bfilename,fields)
    op.makeBasisListFromBDict(bDict,bfilename2)
    with open(bfilename2,'r') as bfile2:
      basisfiletextRecons=bfile2.read()
    os.remove(bfilename)
    os.remove(bfilename2)
    self.assertEqual(self.basisfiletext,basisfiletextRecons)
    
    
if __name__ == '__main__':
  unittest.main()
