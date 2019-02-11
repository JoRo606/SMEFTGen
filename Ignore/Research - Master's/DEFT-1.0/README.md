## DEFT

The code requires Python2.7+ with the sympy package installed. Example scripts are in the tests directory.

The code is not the most user-friendly, nor documented at all, and if you want to do something not covered by the example scripts, it might be easier to send me an email (dwsuth at ucsb dot edu). I'm happy to help if possible.

## Implementing a new basis

Suppose one has created a list of Term instances (operators), called `terms`, in DEFT. Suppose also that one knows algebraic expressions for a set of basis operators in terms of the former:

b1= c1\*t1 + ... + cN\*tN

and so on. Here b1 is a basis operator, t1-tN are all in `terms`, and c1-cN are complex coefficients. 
This can be encoded in a text file as follows. Each line corresponds to a basis operator and is of the form

`[(sympyStr1,index1/reprStr1),...,(sympyStrN,indexN/reprStrN)]:nameStr`

where

* nameStr is a string that names the basis operator b1
* each sympyStr is a string literal which will be passed to sympy.sympify to generate a sympy expression for the corresponding coefficient c1-cN
* each reprStr is the `__repr__` of the corresponding Term t1-tN in `terms`. Alternatively, an integer index gives the position of the appropriate Term in `terms`

See tests/gSMDim6WarsawRecons.txt and tests/gSMDim6SILMRecons.txt for examples of this syntax.

Having constructed the text file, it can be read in with the following syntax:

`basisDict=DEFT.operators.readBasisList(terms,"/path/to/file.txt",fields)`

where fields is a list of the field classes whose instances appear in `terms`. See e.g. the method `test_changeBasisSM()` in tests/testSMBases.py. The returned `basisDict` has the form

`collections.OrderedDict([(nameStr1,termsInst1),(nameStr2,termsInst2),...])`

where each nameStr is the name of a basis operator, and each termsInst is an instance of the Terms class, which encodes the linear expression for the corresponding basis operator in its attributes

`termsInst.terms=[t1,...,tN]`
`termsInst.weights=[c1,...,cN]`

Then, to change between two bases encoded by startBasisDict and endBasisDict, one calls

`endtms,startToEndMat=DEFT.operators.changeBasis(terms,rrefrls,startBasisDict,endBasisDict)`

with `rrefrls` encoding the reduced row echelon form of the redundancy relations for the list `terms` (see tests/testSMBases.py for how to produce this object). `startToEndMat` expresses one operator in the `startBasisDict` as a sum of the operators in `endBasisDict` via

`startBasisDict.keys()[i] = str(startToEndMat[i][j])+endBasisDict.keys()[j]`

where one sums over the index j for the full expression.
