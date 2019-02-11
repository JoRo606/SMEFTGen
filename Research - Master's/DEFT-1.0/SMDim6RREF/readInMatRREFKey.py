import sys
sys.path.insert(0,'..')

import DEFT
import DEFT.operators

##read key file
with open('matRREFkey','r') as f:
  lines=f.readlines()

##create term objects from repr strings
terms=[DEFT.Term.unrepr(line.split(':')[1].strip(),DEFT.operators.listOfSMFields) for line in lines]

##print out as tex
tex='\n\n'.join([str(i)+' ~ : ~ '+t.texRep() for i,t in enumerate(terms)])
DEFT.operators.makePDFDoc(tex,'listOfSMDim6Monomials')
