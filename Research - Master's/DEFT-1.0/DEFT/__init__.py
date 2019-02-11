import collections
import itertools
import sympy
import math

groupsDimFund={'lorL':2,'lorR':2,'su2':2,'su3':3,'su4':4}
U1GaugeGroupsCoupling={'Y':sympy.symbols('gp')}
NAGaugeGroupsCoupling={'su2':sympy.symbols('g'),'su3':sympy.symbols('gs')}

def mash(term):
  return term.__mash__()

class Index:
  def __str__(self):
    return self.group+' '+str(self.label)+(' U' if self.upper else ' L')+(' D' if self.dummy else ' F')
  def __repr__(self):
    return self.__str__()
  def __eq__(self,i):
    return self.group==i.group and self.label==i.label and self.upper==i.upper and self.dummy==i.dummy and self.value==i.value
  def __hash__(self):
    return hash(self.group) | hash(self.label) | hash(self.upper) | hash(self.dummy) | hash(self.value)
  def _DEFTcopy(self):
    return Index(self.group,self.label,self.upper,self.dummy,self.value)
  def __init__(self,group,label,upper=True,dummy=False,value=None):
    self.group=group
    self.label=label
    self.upper=upper
    self.dummy=dummy
    self.value=value
class frac:
  def gcd(self,a, b):
    if b>a: a,b=b,a
    while b:      
      a, b = b, a % b
    return a
  def lcm(self,a, b):
    return a * b / self.gcd(a, b)
  def _simp(self):
    gcd=self.gcd(self.frac[0],self.frac[1])
    gcd=(gcd if gcd>0 else -gcd)
    self.frac=(self.frac[0]/gcd,self.frac[1]/gcd)
    return
  def __str__(self):
    return str(self.frac[0]) if self.frac[1]==1 else str(self.frac[0])+'/'+str(self.frac[1])
  def __add__(self,frac2):
    lcm=self.lcm(self.frac[1],frac2.frac[1])
    res=frac(self.frac[0]*lcm/self.frac[1] + frac2.frac[0]*lcm/frac2.frac[1],lcm)
    res._simp()
    return res
  def __sub__(self,frac2):
    lcm=self.lcm(self.frac[1],frac2.frac[1])
    res=frac(self.frac[0]*lcm/self.frac[1] - frac2.frac[0]*lcm/frac2.frac[1],lcm)
    res._simp()
    return res
  def __mul__(self,frac2):
    gcd1=self.gcd(self.frac[0],frac2.frac[1])
    gcd2=self.gcd(self.frac[1],frac2.frac[0])
    totgcd=gcd1*gcd2
    return frac(self.frac[0]*frac2.frac[0]/totgcd,self.frac[1]*frac2.frac[1]/totgcd)
  def __div__(self,frac2):
    gcd1=self.gcd(self.frac[0],frac2.frac[0])
    gcd2=self.gcd(self.frac[1],frac2.frac[1])
    totgcd=gcd1*gcd2
    return frac(self.frac[0]*frac2.frac[1]/totgcd,self.frac[1]*frac2.frac[0]/totgcd)
  def __cmp__(self,frac2):
    lcm=self.lcm(self.frac[1],frac2.frac[1])
    return self.frac[0]*lcm/self.frac[1] - frac2.frac[0]*lcm/frac2.frac[1]
  def texRep(self):
    return str(self.frac[0]) if self.frac[1]==1 else r'\frac{'+str(self.frac[0])+r'}{'+str(self.frac[1])+'}'
  def _DEFTcopy(self):
    return frac(self.frac[0],self.frac[1])
  def __init__(self,num,den):
    if den<0:
      self.frac=(-1*num,-1*den)
    else:
      self.frac=(num,den)

class Symmetry:
  def __init__(self,field,inds):
    self.indices=inds
    self.field=field

class Unitary(Symmetry):
  def sym(self):
    if self.field.getTerm() is None: return True
    if len(self.field.Dindices)>0: return True
    holyIndex=self.indices[0]
    term=self.field.getTerm()
    invs=term._invariants
    invmatch=next((inv for inv in invs if holyIndex in inv.indices and inv.name=='del'), None)
    if invmatch is None: return True
    otherind=next(i for i in invmatch.indices if i!=holyIndex)
    otherfield=next(f for f in term if otherind in f.indices)
    if not (otherfield._basename==self.field._basename and otherfield.isConjugate!=self.field.isConjugate and len(otherfield.Dindices)==0): return True
    term.removeField(otherfield)
    term.removeField(self.field)
    term.removeInvariant(invmatch)
    term.refreshUOnes()
    return True
  def _DEFTcopy(self,newfield,newindices):
    return Unitary(newfield,[next(i for i in newindices if i==ind) for ind in self.indices])
  def __str__(self):
    return 'Unitary:\n'+'\n'.join('  '+str(i)+' '+str(id(i)) for i in self.indices)  
  def __init__(self,field,inds):
    if len(inds)!=1: raise Exception('Unitary symmetry needs one index only')
    Symmetry.__init__(self,field,inds)
    #self.indices=inds
    #self.indices=[i]

class Traceless(Symmetry):
  def sym(self):
    if self.field.getTerm() is None: return True
    for inv in self.field.getTerm()._invariants:
      if inv.name=='del' and all(ind in inv.indices for ind in self.indices): return False
    return True
  def _DEFTcopy(self,newfield,newindices):
    return Traceless(newfield,[next(i for i in newindices if i==ind) for ind in self.indices])
  def __str__(self):
    return 'Traceless:\n'+'\n'.join('  '+str(i)+' '+str(id(i)) for i in self.indices)  
  def __init__(self,field,inds):
    if len(inds)!=2: raise Exception('Traceless symmetry needs exactly two indices')
    Symmetry.__init__(self,field,inds)
    #self.indices=inds

class Symmetric(Symmetry):
  def sym(self):
    if self.field.getTerm() is None: return True
    for inv in self.field.getTerm()._invariants:
      if 'eps' in inv.name and len([i for i in self.indices if i in inv.indices])>1: return False
    return True
  def _DEFTcopy(self,newfield,newindices):
    return Symmetric(newfield,[next(i for i in newindices if i==oldi) for oldi in self.indices])
  def __str__(self):
    return 'Symmetric:\n'+'\n'.join('  '+str(i)+' '+str(id(i)) for i in self.indices)  
  def __init__(self,field,inds):
    Symmetry.__init__(self,field,inds)
    #self.indices=inds

groupdefs = {
'su4':{'1':(0,0,[]),'4':(1,0,[]),'4bar':(0,1,[]),'15':(1,1,[(Traceless,[0,1])])},
'su3':{'1':(0,0,[]),'3':(1,0,[]),'3bar':(0,1,[]),'8':(1,1,[(Traceless,[0,1])])},
'sug':{'1':(0,0,[]),'3':(1,0,[]),'3bar':(0,1,[]),'6':(2,0,[(Symmetric,[0,1])]),'8':(1,1,[(Traceless,[0,1])]),'10':(3,0,[(Symmetric,range(3))]),'15p':(4,0,[(Symmetric,range(4))]),'15':(2,1,[(Symmetric,range(2))]),'24':(1,3,[(Symmetric,range(1,4))]),'27':(2,2,[(Symmetric,range(2)),(Symmetric,range(2,4))])},
'su2':{'1':(0,0,[]),'2':(1,0,[]),'3':(1,1,[(Traceless,[0,1])]),'7':(6,0,[(Symmetric,range(6))])},
'lorL':{'1':(0,0,[]),'2':(1,0,[]),'2l':(0,1,[]),'3':(0,2,[(Symmetric,range(2))])},
'lorR':{'1':(0,0,[]),'2':(1,0,[]),'2l':(0,1,[]),'3':(0,2,[(Symmetric,range(2))])}
}

class Field:
  def makeIndices(self):
    groups=self.NADict.keys()
    groups.sort()
    for group in groups:
      recentIs=[]
      for i in range(groupdefs[group][self.NADict[group]][0]):
        no=(i+1 if self._term is None else self._term.smallestNewLabel(group)+i)
        self.indices.append( Index(group,no,True) )
        recentIs.append(self.indices[-1])
      for i in range(groupdefs[group][self.NADict[group]][1]):
        no=(i+1 if self._term is None else self._term.smallestNewLabel(group)+i)+groupdefs[group][self.NADict[group]][0]
        self.indices.append( Index(group,no,False) )
        recentIs.append(self.indices[-1])
      for sym,inds in groupdefs[group][self.NADict[group]][2]:
        self.symmetries.append(  sym(self,[recentIs[j] for j in inds])  )
  def addGaugeFuncDicts(self):
    NAgroups=NAGaugeGroupsCoupling.keys()
    indgroups=set(i.group for i in self.indices)
    for indgroup in indgroups:
      if indgroup in NAgroups:
        if isinstance(self,gaugeGroupsFields[indgroup]): self.funcDiffDict[gaugeGroupsFields[indgroup]._vectorPotential._basename]='self._FDNAFS(field,\''+indgroup+'\')'
        else: self.funcDiffDict[gaugeGroupsFields[indgroup]._vectorPotential._basename]='self._FDNA(field,\''+indgroup+'\')'
    U1groups=U1GaugeGroupsCoupling.keys()
    for u1g,charge in self.U1Dict.iteritems():
      if u1g in U1groups and charge!=frac(0,1):
        self.funcDiffDict[gaugeGroupsFields[u1g]._vectorPotential._basename]='self._FDU1(field,\''+u1g+'\')'
    return      
  def conjugate(self):
    for i in self.indices:
      if i.group=='lorL': i.group='lorR'
      elif i.group=='lorR': i.group='lorL'
      else: i.upper = not i.upper
    self.Dindices=[(y,x) for x,y in self.Dindices]
    if self.isSelfConjugate: return
    if self.isConjugate: self._name=self._name[:-1]
    else: self._name+='c'
    #self.tex='{'+self.tex+'}^c'
    for k,v in self.U1Dict.iteritems():
      self.U1Dict[k]=frac(-v.frac[0],v.frac[1])
    if self.isConjugate: self.funcDiffDict[self._basename]=self.funcDiffDict.pop(self._conjugateName(self._basename))
    else: self.funcDiffDict[self._conjugateName(self._basename)]=self.funcDiffDict.pop(self._basename)
    self.isConjugate = not self.isConjugate
  def _conjugateName(self,n):
    if n[-1]=='c': return n[:-1]
    else: return n+'c'
  def replaceIndexInField(self,old,new):
    i=self.indices.index(old)
    self.indices.pop(i)
    self.indices.insert(i,new)
    self.Dindices=[(new if i==old else i,new if j==old else j) for i,j in self.Dindices]
    return
  def replaceIndexInSym(self,old,new):
    for sym in self.symmetries:
      try:
        i=sym.indices.index(old)
        sym.indices.pop(i)
        sym.indices.insert(i,new)
      except ValueError:
        pass
    return
  def smallestNewLabel(self,group):
    if self._term is not None: return self._term.smallestNewLabel(group)
    labels=[i.label for i in self.indices if i.group==group]
    if len(labels)>0: return max(labels)+1
    else: return 1
  def differentiate(self):
    if self._term is None:
      self.indices.append( Index('lorL',self.smallestNewLabel('lorL'),False) )
      self.indices.append( Index('lorR',self.smallestNewLabel('lorR'),False) )
      self.Dindices.append( (self.indices[-2],self.indices[-1]) )
    else:
      self.indices.append( Index('lorL',self._term.smallestNewLabel('lorL'),False) )
      self.indices.append( Index('lorR',self._term.smallestNewLabel('lorR'),False) )
      self.Dindices.append( (self.indices[-2],self.indices[-1]) )
      self._term._clear()
    self.dimension+=frac(1,1)
  #BCTI compliant 2016-02-19
  def _FDSelf(self,field):
    #if len(self.Dindices)==0: return [Term([one()])] 
    #print 'fdself',field
    newf=one()
    for di in self.Dindices:
      newf.indices.append(di[0])
      newf.indices.append(di[1])
      newf.Dindices.append( (di[0],di[1]) )
      newf.dimension+=frac(1,1)
    rt=Term()
    rt.add(newf,True)
    invpairs=zip(field.indices,self.indices[:len(field.indices)])
    for invp in invpairs:
      flippedinv=invp[0]._DEFTcopy()
      flippedinv.upper=not flippedinv.upper
      rt._invariants.append(Invariant('del','\delta',[flippedinv,invp[1]]))
    if self.isFermion() and self._term is not None:
      if ([f.isFermion() for f in self._term._fields[:self._term._fields.index(self)]].count(True) % 2 == 0): return [(1,rt)]
      else: return [(-1,rt)]
    else:
      return [(1,rt)]
  #BCTI compliant 2016-02-19
  def _FDFS(self,field):
    #upperi=(self.indices[0] if self.indices[0].upper else self.indices[1])
    #loweri=(self.indices[1] if self.indices[0].upper else self.indices[0])
    alpha=self.indices[0]
    beta=self.indices[1]
    fieldlorL=next(i for i in field.indices if 'lorL'==i.group)
    fieldlorR=next(i for i in field.indices if 'lorR'==i.group)
    if len(self.indices)>2 and 'su' in self.indices[2].group:
      gaugeupperi=(self.indices[2] if self.indices[2].upper else self.indices[3])
      gaugeloweri=(self.indices[3] if self.indices[2].upper else self.indices[2])
      fieldgaugeupperi=next(i for i in field.indices if 'su' in i.group and i.upper)
      fieldgaugeloweri=next(i for i in field.indices if 'su' in i.group and not i.upper)
    lorR= (alpha.group=='lorR')
    newf1=one()
    newf1._term=self._term
    newf1.differentiate()
    newf1._term=None
    ####account for extra derivs acting on field strength
    tempDinds=newf1.Dindices
    newf1.Dindices=[]
    for dL,dR in self.Dindices:
      newdL=dL._DEFTcopy()
      newdR=dR._DEFTcopy()
      newf1.Dindices.append( (newdL,newdR) )
      newf1.indices+=[newdL,newdR]
      newf1.dimension+=frac(1,1)
    newf1.Dindices+=tempDinds
    ###end accounting
    newf2=newf1._DEFTcopy()
    ##first term
    dlorLi=newf1.indices[0]
    dlorRi=newf1.indices[1]
    newt1=Term()
    newt1.add(newf1,True)
    dlorLi.dummy=True
    dlorRi.dummy=True
    lorLind=fieldlorL._DEFTcopy()
    lorRind=fieldlorR._DEFTcopy()
    lorLind.upper=not lorLind.upper
    lorRind.upper=not lorRind.upper
    if len(self.indices)>2 and 'su' in self.indices[2].group:
      upperind=fieldgaugeloweri._DEFTcopy()
      upperind.upper=not upperind.upper
      lowerind=fieldgaugeupperi._DEFTcopy()
      lowerind.upper=not lowerind.upper
    if lorR:
      newf1.indices[1]=alpha._DEFTcopy()
      newf1.Dindices[-1]=(dlorLi,newf1.indices[1])
      newt1._invariants.append( Invariant('epsl',r'\epsilon',[dlorLi,lorLind]) )
      #newt1._invariants.append( Invariant('del',r'\delta',[upperi._DEFTcopy(),dlorRi]) )
      newt1._invariants.append( Invariant('del',r'\delta',[lorRind,beta._DEFTcopy()]) )
    else:
      #newt1._invariants.append( Invariant('del',r'\delta',[upperi._DEFTcopy(),dlorLi]) )
      newf1.indices[0]=alpha._DEFTcopy()
      newf1.Dindices[-1]=(newf1.indices[0],dlorRi)
      newt1._invariants.append( Invariant('epsl',r'\epsilon',[dlorRi,lorRind]) )
      newt1._invariants.append( Invariant('del',r'\delta',[lorLind,beta._DEFTcopy()]) )
    if len(self.indices)>2 and 'su' in self.indices[2].group:
      newt1._invariants.append( Invariant('del',r'\delta',[gaugeupperi._DEFTcopy(),lowerind]) )
      newt1._invariants.append( Invariant('del',r'\delta',[gaugeloweri._DEFTcopy(),upperind]) )
    ##second term
    dlorLi=newf2.indices[0]
    dlorRi=newf2.indices[1]
    newt2=Term()
    newt2.add(newf2,True)
    dlorLi.dummy=True
    dlorRi.dummy=True
    lorLind=fieldlorL._DEFTcopy()
    lorRind=fieldlorR._DEFTcopy()
    lorLind.upper=not lorLind.upper
    lorRind.upper=not lorRind.upper
    if len(self.indices)>2 and 'su' in self.indices[2].group:
      upperind=fieldgaugeloweri._DEFTcopy()
      upperind.upper=not upperind.upper
      lowerind=fieldgaugeupperi._DEFTcopy()
      lowerind.upper=not lowerind.upper
    if lorR:
      newf2.indices[1]=beta._DEFTcopy()
      newf2.Dindices[-1]=(dlorLi,newf2.indices[1])
      newt2._invariants.append( Invariant('epsl',r'\epsilon',[lorLind,dlorLi]) )
      #newt2._invariants.append( Invariant('del',r'\delta',[loweri._DEFTcopy(),dlorRi]) )
      newt2._invariants.append( Invariant('del',r'\delta',[lorRind,alpha._DEFTcopy()]) )
    else:
      newf2.indices[0]=beta._DEFTcopy()
      newf2.Dindices[-1]=(newf2.indices[0],dlorRi)
      #newt2._invariants.append( Invariant('del',r'\delta',[loweri._DEFTcopy(),dlorLi]) )
      newt2._invariants.append( Invariant('epsl',r'\epsilon',[lorRind,dlorRi]) )
      newt2._invariants.append( Invariant('del',r'\delta',[lorLind,alpha._DEFTcopy()]) )
    if len(self.indices)>2 and 'su' in self.indices[2].group:
      newt2._invariants.append( Invariant('del',r'\delta',[gaugeupperi._DEFTcopy(),lowerind]) )
      newt2._invariants.append( Invariant('del',r'\delta',[gaugeloweri._DEFTcopy(),upperind]) )
    newf1._term=newt1
    newf2._term=newt2
    return [(1,newt1),(-1,newt2)]
  #BCTI compliant 2016-02-19
  def _FDU1(self,field,group):
    newts=[]
    #symboldict={'Y':sympy.symbols('gp'),'Z':sympy.symbols('gz')}
    for i in range(len(self.Dindices)):
      newf=self._DEFTcopy()
      di=newf.Dindices.pop(i)
      newf.remove(di[0])
      newf.remove(di[1])
      newt=Term()
      newt.add(newf,True)
      lorLind=field.indices[0]._DEFTcopy()
      lorRind=field.indices[1]._DEFTcopy()
      lorLind.upper=not lorLind.upper
      lorRind.upper=not lorRind.upper
      newt._invariants.append( Invariant('del',r'\delta',[di[0],lorLind]) )
      newt._invariants.append( Invariant('del',r'\delta',[di[1],lorRind]) )
      #print 'fieldFDU1,1',newt.longstr()
      newts.append(( sympy.I * U1GaugeGroupsCoupling[group] * self.U1Dict[group].frac[0] / self.U1Dict[group].frac[1] ,newt))
    return newts
  #BCTI compliant 2016-02-19
  def _FDNA(self,field,group):
    newts=[]
    #symboldict={'su2':sympy.symbols('g'),'su3':sympy.symbols('gs')}
    #Ndict={'su2':2,'su3':3}
    fieldlorL=next(i for i in field.indices if 'lorL'==i.group)
    fieldlorR=next(i for i in field.indices if 'lorR'==i.group)
    fieldgaugeupperi=next(i for i in field.indices if 'su' in i.group and i.upper)
    fieldgaugeloweri=next(i for i in field.indices if 'su' in i.group and not i.upper)
    for i in range(len(self.Dindices)):
      for naind in [p for p in self.indices if p.group==group]:
        #term with uncontracted SU(N) indices
        newf=self._DEFTcopy()
        di=newf.Dindices.pop(i)
        newf.remove(di[0])
        newf.remove(di[1])
        thisindexlocation=newf.indices.index(naind)
        thisindex=newf.indices.pop(thisindexlocation)
        newlabel=self.smallestNewLabel(group)
        newf.indices.insert(thisindexlocation, Index(group,newlabel,naind.upper,True) )
        #####below line is new
        newf.replaceIndexInSym(thisindex,newf.indices[thisindexlocation])
        newt=Term()
        newt.add(newf,True)
        lorLind=fieldlorL._DEFTcopy()
        lorRind=fieldlorR._DEFTcopy()
        lorLind.upper=not lorLind.upper
        lorRind.upper=not lorRind.upper
        newt._invariants.append( Invariant('del',r'\delta',[di[0],lorLind]) )
        newt._invariants.append( Invariant('del',r'\delta',[di[1],lorRind]) )
        upperind=fieldgaugeloweri._DEFTcopy()
        upperind.upper=not upperind.upper
        lowerind=fieldgaugeupperi._DEFTcopy()
        lowerind.upper=not lowerind.upper
        if naind.upper:
          #newt._invariants.append( Invariant('del',r'\delta',[thisindex,fieldgaugeupperi._DEFTcopy()]) )
          #newt._invariants.append( Invariant('del',r'\delta',[newf.indices[thisindexlocation],fieldgaugeloweri._DEFTcopy()]) )
          newt._invariants.append( Invariant('del',r'\delta',[thisindex,lowerind]) )
          newt._invariants.append( Invariant('del',r'\delta',[newf.indices[thisindexlocation],upperind]) )
        else:
          #newt._invariants.append( Invariant('del',r'\delta',[thisindex,fieldgaugeloweri._DEFTcopy()]) )
          #newt._invariants.append( Invariant('del',r'\delta',[newf.indices[thisindexlocation],fieldgaugeupperi._DEFTcopy()]) )
          newt._invariants.append( Invariant('del',r'\delta',[thisindex,upperind]) )
          newt._invariants.append( Invariant('del',r'\delta',[newf.indices[thisindexlocation],lowerind]) )
        newts.append((  (1 if naind.upper else -1)*sympy.I*NAGaugeGroupsCoupling[group] ,newt))
        #term with contracted SU(N) indices to make traceless
        newf=self._DEFTcopy()
        di=newf.Dindices.pop(i)
        newf.remove(di[0])
        newf.remove(di[1])
        newt=Term()
        newt.add(newf,True)
        lorLind=fieldlorL._DEFTcopy()
        lorRind=fieldlorR._DEFTcopy()
        lorLind.upper=not lorLind.upper
        lorRind.upper=not lorRind.upper
        newt._invariants.append( Invariant('del',r'\delta',[di[0],lorLind]) )
        newt._invariants.append( Invariant('del',r'\delta',[di[1],lorRind]) )
        upperind=fieldgaugeloweri._DEFTcopy()
        upperind.upper=not upperind.upper
        lowerind=fieldgaugeupperi._DEFTcopy()
        lowerind.upper=not lowerind.upper
        newt._invariants.append( Invariant('del',r'\delta',[upperind,lowerind]) )
        newts.append((  (-1 if naind.upper else 1)*sympy.I*NAGaugeGroupsCoupling[group]/groupsDimFund[group] ,newt))
    return newts
  def _FDNAFS(self,field,group):
    return self._FDNA(field,group)+self._FDFS(field)
  def remove(self,index):
    self.indices.remove(index)
  def funcDiff(self,field):
    return eval(self.funcDiffDict[field.getName()])
  def getTerm(self):
    return self._term
  def getName(self):
    return 'D'*len(self.Dindices)+self._name
  def __str__(self):
    return self.getName()+'\n'+'\n'.join([str(i) for i in self.indices])+'\nU(1) charges\n'+'\n'.join([k+' '+str(self.U1Dict[k]) for k in self.U1Dict.keys()])
  def _DEFTcopy(self):
    newu1dict=dict((k,v._DEFTcopy()) for k,v in self.U1Dict.iteritems())
    newnadict=dict((k,v) for k,v in self.NADict.iteritems())
    newf=Field(self.dimension._DEFTcopy(),newnadict,newu1dict,self._name,self.tex,self._term)
    newf._basename=self._basename
    newf.isConjugate=self.isConjugate
    newf.isSelfConjugate=self.isSelfConjugate
    newf.funcDiffDict=dict(self.funcDiffDict)
    newf.indices=[i._DEFTcopy() for i in self.indices]
    newf.Dindices=[(next(i for i in newf.indices if i==li),next(i for i in newf.indices if i==ri)) for li,ri in self.Dindices]
    newf.pDindices=[(next(i for i in newf.indices if i==li),next(i for i in newf.indices if i==ri)) for li,ri in self.pDindices]
    newf.symmetries=[s._DEFTcopy(newf,newf.indices) for s in self.symmetries]
    return newf
  def isFermion(self):
    return (sum(1 for i in self.indices if 'lor' in i.group) % 2)==1
  def isVector(self):
    return sum(1 for i in self.indices if 'lorL'==i.group)==1 and sum(1 for i in self.indices if 'lorR'==i.group)==1
  def isFieldStrength(self):
    return sum(1 for i in self.indices if i.group in ['lorL','lorR'])==2 and not self.isVector()
  def __init__(self,dimension,nonabel,uones,name,tex,term=None):
    self.dimension=dimension
    self.NADict = nonabel
    self.U1Dict = uones
    self._name = name
    #self._basename = name
    self.tex = tex
    self.indices = []
    self.Dindices = []
    self.pDindices = []
    self.symmetries = []
    self._term = term
    self.isConjugate=False
    self.isSelfConjugate=False
    self.funcDiffDict={name:'self._FDSelf(field)'}
    self.makeIndices()
    self.addGaugeFuncDicts()

#from fields import Field
class one(Field):
  _basename='1'
  def __init__(self,term=None):
    Field.__init__(self,frac(0,1),{},{},'1','1',term)
    self.isSelfConjugate=True

class Plc(Field):
  def add(self,field):
    if self._term is not None: self._term._clear()
    self.indices+=field.indices
    self.dimension+=field.dimension
    for k,v in field.U1Dict.items():
      try: self.U1Dict[k]+=v
      except KeyError: self.U1Dict[k]=v
    return
  def removeDummies(self,invariants):
    for inv in invariants:
      if any(ind not in self.indices for ind in inv.indices): continue
      for ind in inv.indices: self.indices.remove(ind)
    return
    #newIs=[]
    #for i in self.indices:
    #  if i not in newIs: newIs.append(i)
    #self.indices=newIs
  def explode(self):
    if self._term is None: return
    newplcs=[]
    for i,ind in enumerate(self.indices):
      newplc=Field(frac(0,1),{},{},'plc'+str(i),'plc'+str(i),self._term)
      newplc._basename='plc'
      newplc.indices.append(ind)
      newplcs.append(newplc)
    self._term._fields.remove(self)
    self._term._fields+=newplcs
    return
  def _DEFTcopy(self):
    newf=Plc(self._term)
    newf.indices=[i._DEFTcopy() for i in self.indices]
    return newf
  _basename='Plc'
  def __init__(self,term=None):
    Field.__init__(self,frac(0,1),{},{},'Plc','Plc',term)

class Invariant:
  def __str__(self):
    return self.name+'\n'+'\n'.join([str(i) for i in self.indices])
  def __eq__(self,otherinv):
    return self.name==otherinv.name and self.tex==otherinv.tex and self.indices==otherinv.indices
  def __hash__(self):
    hsh=hash(self.name)
    for i in self.indices: hsh |= hash(i)
    return hsh
  def group(self):
    return self.indices[0].group
  def _DEFTcopy(self,newindices):
    invinds=[]
    extrainds=[]
    for n in range(len(self.indices)):
      invind=next((i for i in newindices if i==self.indices[n]),None)
      if invind is None:
        invind=self.indices[n]._DEFTcopy()
        extrainds.append(invind)
      invinds.append(invind)
    return (Invariant(self.name,self.tex,invinds),extrainds)
  def __init__(self,name,tex,indices=[]):
    self.name=name
    self.tex=tex
    self.indices=indices

class Term:
  def __getitem__(self,index):
    return self._fields[index]
  def _permuteDicts(self):
    if self._pds is not None: return self._pds
    #hash(self)
    dicts=[{'del':'del','epsu':'epsu','epsl':'epsl'}]
    cnt = self._getCounter()
    for f in cnt.keys():
      if cnt[f]==1:
        for d in dicts: d[f+'_0']=f+'_0'
      else:
        names=[f+'_'+str(i) for i in range(cnt[f])]
        permutenames=list(itertools.permutations(names))#[1:]
        newDicts=[]
        for pn in permutenames:
          newBit = dict((names[i],pn[i]) for i in range(cnt[f]))
          for d in dicts:
            cd = d.copy()
            cd.update(newBit)
            newDicts.append(cd)
        dicts=newDicts
    for d in dicts:
      extraDstuff={}
      for k,v in d.items():
        for i in range(k.count('D')):
          extraDstuff[str(i)+'_'+k]=str(i)+'_'+v
      d.update(extraDstuff)
    self._pds=dicts
    return self._pds
  def _permuteNL(self,nl,pd):
    newnl=collections.Counter()
    for tup in nl.keys():
      newtup=( (frozenset(pd[fstr] for fstr in tup[0][0]),frozenset(pd[fstr] for fstr in tup[0][1])) ,tup[1],tup[2])
      newnl[newtup]=nl[tup]
    return newnl
  def _permNL(self):
    if self._pnls is not None: return self._pnls
    self._pnls=[]
    nl=self._getNumberedLinks(self._invariants)
    for pd in self._permuteDicts():
      newnl=collections.Counter()
      for tup in nl.keys():
        newtup=( (frozenset(pd[fstr] for fstr in tup[0][0]),frozenset(pd[fstr] for fstr in tup[0][1])) ,tup[1],tup[2])
        newnl[newtup]=nl[tup]
      self._pnls.append(newnl)
    return self._pnls
  def _cmpParity(self,l1,l2):
    isPos=True
    for pair in itertools.combinations(l1,2):
      if l2.index(pair[0]) > l2.index(pair[1]): isPos = not isPos
    return isPos
  def sign(self):#,newfields=None):
    #newfields=sorted(self._fields,key= lambda x : x.getName())
    def findField(ind):
      return next(f for f in self if ind in f.indices)
    def findFields(inds):
      return [next(f for f in self if i in f.indices) for i in inds]
    rankDict={}
    for f in self:
      rankDict[id(f)]={'epsu':0,'epsl':0,'del':0}
    for inv in self._invariants:
      for i in inv.indices: rankDict[id(findField(i))][inv.name]+=1
    ranking=dict((k,100*v['epsu']+10*v['epsl']+v['del']) for k,v in rankDict.iteritems())
    newfields=sorted(self._fields,key=lambda x: x.getName()+str(ranking[id(x)]))
    degeneracies=collections.Counter(ranking.values())
    for r,degen in degeneracies.iteritems():
      if degen==1: continue
      degenerateFields=[f for f in newfields if ranking[id(f)]==r]
      connections={}
      for df in degenerateFields:
        cons=[]
        for i in df.indices:
          inv=next(iv for iv in self._invariants if i in iv.indices)
          cons.append( inv.name+inv.indices[0].group+''.join(sorted([f.getName() for f in findFields(inv.indices)])) )
        cons.sort()
        connections[id(df)]= ''.join(cons)
      newfields=[f for f in newfields if ranking[id(f)] < r]+sorted(degenerateFields,key=lambda x: connections[id(x)])+[f for f in newfields if ranking[id(f)] > r] 
    newindsids=list(itertools.chain.from_iterable( [id(ind) for ind in f.indices] for f in newfields ))
    isPos=self._cmpParity([f for f in newfields if f.isFermion()],self._fields)
    for inv in self._invariants:
      if 'eps' in inv.name and not self._cmpParity([id(ind) for ind in inv.indices],newindsids): isPos = not isPos
    return isPos
  def isign(self):
    return 1 if self.sign() else -1
  def fixDerivativeIndices(self):
    self._clear()
    changed=False
    for f in self:
      if len(f.Dindices)==0: continue
      dindflat=list(itertools.chain.from_iterable(f.Dindices))
      #print dindflat
      oldindices=f.indices[:]
      for i in dindflat: f.indices.remove(i)
      f.indices+=dindflat
      if oldindices!=f.indices: changed=True
    return changed
  def _getSign(self,fieldList,trm):
    newindsids=list(itertools.chain.from_iterable( [id(ind) for ind in f.indices] for f in fieldList ))
    ###add free indices
    freeinds=[i for i in list(itertools.chain.from_iterable( inv.indices for inv in trm._invariants )) if not i.dummy]
    freeinds.sort(key=str)
    newindsids=newindsids+[id(i) for i in freeinds]
    ###
    isPos=trm._cmpParity([f for f in fieldList if f.isFermion()],trm._fields)
    for inv in trm._invariants:
      if 'eps' in inv.name and not trm._cmpParity([id(ind) for ind in inv.indices],newindsids): isPos = not isPos
    return isPos
  def compareSign(self,othert,force=False):
    if not othert._getCounter()==self._getCounter(): raise Exception('comparing sign of unequal operators')
    setOfSigns=set()
    for m,pd in enumerate(self._permuteDicts()):
      if not othert._permNL()[0]==self._permNL()[m]: continue
      selffd=self._getFieldDict()
      othertfd=othert._getFieldDict()
      #order self fields according to equality with othert fields
      newselffields=[ next(fd for fd in self._fields if pd[selffd[id(fd)]]==othertfd[id(f)]) for f in othert._fields ]
      setOfSigns.add(-1 if (self._getSign(othert._fields,othert) ^ self._getSign(newselffields,self)) else 1)
    if len(setOfSigns)>1: raise Exception('operator zero by sym')
    elif len(setOfSigns)==0: raise Exception('comparing sign of unequal operators')
    else: return setOfSigns.pop()
  def zeroBySym(self):
    self.fixDerivativeIndices()
    try: self.compareSign(self)
    except Exception: return True
    else: return False
  def __eq__(self,othert):
    if not self._getCounter()==othert._getCounter(): return False
    if not self._getLinks(self._invariants)==othert._getLinks(othert._invariants): return False
    nl1=self._getNumberedLinks(self._invariants)
    nl2=othert._getNumberedLinks(othert._invariants)
    if nl1==nl2: return True
    for pnl in self._permNL():
      if nl2==pnl: return True
    return False
  def __hash__(self):
    if self._hash is not None: return self._hash
    pds=self._permuteDicts()
    nl=self._getNumberedLinks(self._invariants)
    hashes=[hash(tuple(sorted(self._permuteNL(nl,pd).items()))) for pd in pds]
    self._hash=hash(tuple(sorted(self._getCounter().items())))
    for h in hashes:
      self._hash |= h
    return self._hash
  def __mash__(self):
    if self._mash is not None: return self._mash
    #hash(self)
    pds=self._permuteDicts()
    nl=self._getNumberedLinks(self._invariants)
    mashes=[frozenset(self._permuteNL(nl,pd).items()) for pd in pds]
    self._mash=set(frozenset(self._getCounter().items()))
    for l in mashes:
      self._mash.add(l)
    self._mash=frozenset(self._mash)
    return self._mash
  def contains(self,othert):
    if not self._compareCounters(othert._getCounter(),self._getCounter()): return False
    nl2=othert._getNumberedLinks(othert._invariants)
    for pnl in self._permNL():
      if self._compareCounters(nl2,pnl): return True
    return False
#################
#################
#####
##if you're about to read the following method, I'm really really sorry. DS
#####
  def quotientPlaceholderBetter(self,othert):
    #othert is a smaller term, possibly with free indices, whic we hope to find in this term and replace with a placeholder with the same 'free' indices
    def getFields(inds,trm):
      return [next(f for f in trm if i in f.indices) for i in inds]
    #if not indexConsistency(othert): raise Exception('Indices of term given as argument are inconsistent')
    #if not indexConsistency(self): raise Exception('Indices of this term are inconsistent')
    self.fixDerivativeIndices()
    othert.fixDerivativeIndices()
    for grp,dim in groupsDimFund.iteritems():
      if dim<=2: continue
      if sum(1 if 'epsu' in inv.name and inv.group()==grp else 0 for inv in self._invariants)>0 and sum(1 if 'epsl' in inv.name and inv.group()==grp else 0 for inv in self._invariants)>0: return []
    #check field content is a subset
    if not self._compareCounters(othert._getCounter(),self._getCounter()):
      return []
    #will ultimately return a version of the following list
    quotients=[]
    #divide othert indices into classes
    AllIn=[inv for inv in othert._invariants if all(i.dummy for i in inv.indices)]
    PartIn=[inv for inv in othert._invariants if not all(i.dummy for i in inv.indices) and any(i.dummy for i in inv.indices)]
    NotIn=[inv for inv in othert._invariants if not any(i.dummy for i in inv.indices)]
    nl2=othert._getNumberedLinks(AllIn)
    words=False
    if words: print 'in quotientPlaceholder',self.shortstr(),othert.shortstr()
    for m,pd in enumerate(self._permuteDicts()):
      if not self._compareCounters(nl2,self._permNL()[m]): continue
      #OK, we have a ball game. Make a copy of self and put in canonical form
      quot=self._DEFTcopy()
      ipd=dict((v, k) for k, v in pd.iteritems())
      ifd=dict((v, k) for k, v in quot._getFieldDict().iteritems())
      if words: print 'ifd',ifd
      if words: print 'ipd',ipd
      othertFieldDict=othert._getFieldDict()
      quotFieldIDs=[ ifd[ipd[othertFieldDict[id(f)]]] for f in othert ]
      if words: print 'quotFieldIDs',quotFieldIDs,'quotNFs',[ipd[othertFieldDict[id(f)]] for f in othert]
      QAllIn=[inv for inv in quot._invariants if all(id(f) in quotFieldIDs for f in getFields(inv.indices,quot))]
      QPartIn=[inv for inv in quot._invariants if not all(id(f) in quotFieldIDs for f in getFields(inv.indices,quot)) and any(id(f) in quotFieldIDs for f in getFields(inv.indices,quot))]
      QNotIn=[inv for inv in quot._invariants if not any(id(f) in quotFieldIDs for f in getFields(inv.indices,quot))]
      if words:
        print 'QAllIn'
        for inv in QAllIn: print inv
        print 'QPartIn'
        for inv in QPartIn: print inv
        print 'QNotIn'
        for inv in QNotIn: print inv
      #check for equality of numberedLinks from AllIn and QAllIn with this permuteDict
      if nl2!=quot._permuteNL( quot._getNumberedLinks(QAllIn) , pd ): continue
      if len(QAllIn)!=len(AllIn): raise Exception('Mismatch between number of invariants in quotient placeholder')
      ###stuff for sgn
      #put quot in canonical form with othert bits to the right
      quot._fields=[f for f in quot if id(f) not in quotFieldIDs]+[next(f for f in quot._fields if id(f)==idf) for idf in quotFieldIDs]
      for inv in quot._invariants:
        if 'eps' not in inv.name: continue
        inv.indices=[i for i in inv.indices if id(getFields([i],quot)[0]) not in quotFieldIDs]+[i for i in inv.indices if id(getFields([i],quot)[0]) in quotFieldIDs]
      isPos=(True if self.compareSign(quot)==1 else False)
      othertIndsIds=list(itertools.chain.from_iterable( [id(ind) for ind in f.indices] for f in othert ))
      quotIndsIds=list(itertools.chain.from_iterable( [id(ind) for ind in f.indices] for f in quot ))
      for inv in QAllIn:
        if 'eps' in inv.name and not quot._cmpParity([id(ind) for ind in inv.indices],quotIndsIds): isPos = not isPos
      for inv in AllIn:
        if 'eps' in inv.name and not othert._cmpParity([id(ind) for ind in inv.indices],othertIndsIds): isPos = not isPos
      ####done sign stuff for now
      othertFieldIndices=itertools.chain.from_iterable(f.indices for f in othert)
      othertInvIndices=itertools.chain.from_iterable(inv.indices for inv in othert._invariants)
      othertSortedFreeFieldIndices=sorted([i for i in othertFieldIndices if not i.dummy],key=str)
      othertSortedFreeInvIndices=sorted([i for i in othertInvIndices if not i.dummy],key=str)
      othertSortedFreeIndices=sorted(othertSortedFreeFieldIndices+othertSortedFreeInvIndices,key=str)
      #fullyInsideOthertInds=list(itertools.chain.from_iterable(inv.indices for inv in fullyInsideOthertInvs))
      if words:
        print '**********\nstarting point othert\n***********',othert.longstr()
        print '**********\nstarting point quot\n***********',quot.longstr()
        print '*******************'
      ###find the possible coresponding quot indices for free indices belong to a field of othert
      othertToQuotIndexDict={}
      for i in othertSortedFreeFieldIndices:
        othertfield=next(f for f in othert if i in f.indices)
        quotfield=next(f for f in quot if id(f)==ifd[ipd[othertFieldDict[id(othertfield)]]])
        othertToQuotIndexDict[i]=[ind for ind in quotfield.indices if ind.group==i.group and ind.upper==i.upper and not any(ind in iv.indices for iv in QAllIn)]
        if words:
          print 'othertfield',othertfield
          print 'quotfield',quotfield
          print 'othert free index',i
          print 'quot index',othertToQuotIndexDict[i]
      ###for each othert invariant with free indices, identify corresponding quotient invariants
      othertToQuotPartInDict={}
      for inv in PartIn:
        othertToQuotPartInDict[inv]=[]
        for iv in QPartIn:
          if len(iv.indices)!=len(inv.indices) or iv.group()!=inv.group(): continue
          isRightInv=True
          for f in othert:
            quotField=next( fd for fd in quot if id(fd)==ifd[ipd[othertFieldDict[id(f)]]] )
            othertSet=collections.Counter( (i.group,i.upper) for i in inv.indices if i in f.indices )
            quotSet=collections.Counter( (i.group,i.upper) for i in iv.indices if i in quotField.indices )
            if othertSet!=quotSet:
              isRightInv=False
              break
            else:
              if words: print othertSet,quotSet
          if isRightInv:
            othertToQuotPartInDict[inv].append(iv)
        ###if none of the quotient invariants match, give up
        if len(othertToQuotPartInDict[inv])==0: break
        if words:
          print 'correspondingQuotInv',othertToQuotPartInDict[inv]
          #print 'namecheck',correspondingQuotInv.name,inv.name
      if any(len(v)==0 for v in othertToQuotPartInDict.values()): continue
      if len(NotIn)>1: print 'WARNING: This code is untested on more than one NotIn inv'
      othertToQuotNotInDict={}
      for inv in NotIn:
        othertToQuotNotInDict[inv]=[(True,iv) for iv in QNotIn if iv.group()==inv.group() and len(iv.indices)==len(inv.indices) and not any(id(f) in quotFieldIDs for f in getFields(iv.indices,quot))]
        if len(othertToQuotNotInDict[inv])==0: break
        if inv.group()=='su2': othertToQuotNotInDict[inv]=othertToQuotNotInDict[inv]+[(False,iv) for dummy,iv in othertToQuotNotInDict[inv]]
      if any(len(v)==0 for v in othertToQuotNotInDict.values()): continue
      #####new bit to make lots of quots
      possibleFreeIndexArrangements=itertools.product(*tuple(othertToQuotIndexDict[i] for i in othertSortedFreeFieldIndices))
      possibleFreeIndexArrangements=[x for x in possibleFreeIndexArrangements if len(x)==len(set(x))]
      if words:
        print 'pfia', possibleFreeIndexArrangements
        print '*****'
        for x in possibleFreeIndexArrangements: print x
      possiblePartInArrangements=itertools.product(*tuple(othertToQuotPartInDict[inv] for inv in PartIn))
      possiblePartInArrangements=[x for x in possiblePartInArrangements if len(x)==len(set(x))]
      if words:
        print 'ppia',possiblePartInArrangements
        print '*****'
        for x in possiblePartInArrangements: print x
      possibleNotInArrangements=itertools.product(*tuple(othertToQuotNotInDict[inv] for inv in NotIn))
      possibleNotInArrangements=[x for x in possibleNotInArrangements if len(x)>0 and len(zip(*x)[1])==len(set(zip(*x)[1]))]
      if len(possibleNotInArrangements)==0: possibleNotInArrangements=[()]
      if words:
        print 'pnia',possibleNotInArrangements
        print '*****'
        for x in possibleNotInArrangements: print x
      def getIndex(ind,term):
        for f in term:
          for i in f.indices:
            if i==ind: return i
        raise Exception('ERROR: can\'t find index')
      def getInv(inv,term):
        for iv in term._invariants:
          if iv==inv: return iv
        raise Exception('ERROR: can\'t find invariant')
      for pfia in possibleFreeIndexArrangements:
        for ppia in possiblePartInArrangements:
          for pnia in possibleNotInArrangements:
            newquot=quot._DEFTcopy()
            newisPos=isPos
            newifd=dict((k,next(id(newquot._fields[i]) for i,f in enumerate(quot._fields) if id(f)==v)) for k, v in ifd.iteritems())
            if words: print 'newifd',newifd
            if words: print 'ipd',ipd
            newquot._clear()
            placeholder=Plc(newquot)
            newquot._fields.append(placeholder)
            placeholder.indices=[ getIndex(i,newquot) for i in pfia ]
            othertToNewQuotIndexDict=dict((i,j) for i,j in zip(othertSortedFreeFieldIndices,placeholder.indices))
            for inv,Qinv in zip(PartIn,ppia):
              if words: print 'othert part inv',inv
              partialInvFields=[]
              newInds=[]
              smallestNewLabel=newquot.smallestNewLabel(inv.group())
              for i in inv.indices:
                if i.dummy:
                  othertfield=next(f for f in othert if i in f.indices)
                  quotfield=next(f for f in newquot if id(f)==newifd[ipd[othertFieldDict[id(othertfield)]]])
                  partialInvFields.append(quotfield)
                else:
                  newi=Index(i.group,smallestNewLabel+len(newInds),i.upper,True)
                  newInds.append(newi)
                  othertToNewQuotIndexDict[i]=newi
              correspondingQuotInv=getInv(Qinv,newquot)
              ####sign stuff
              quotIndsIds=list(itertools.chain.from_iterable( [id(ind) for ind in f.indices] for f in newquot ))
              if 'eps' in correspondingQuotInv.name and not newquot._cmpParity([id(ind) for ind in correspondingQuotInv.indices if any(i in f.indices for f in partialInvFields)],quotIndsIds): newisPos = not newisPos
              if 'eps' in inv.name and not othert._cmpParity([id(ind) for ind in inv.indices if ind.dummy],othertIndsIds): newisPos = not newisPos
              ###end sign stuff
              newinv=Invariant('eps','\epsilon',[i for i in correspondingQuotInv.indices if all( i not in f.indices for f in partialInvFields )]+newInds )
              if words: print 'newinv',len(newInds),len([i for i in correspondingQuotInv.indices if all( i not in f.indices for f in partialInvFields )])
              if newinv.indices[0].upper: newinv.name='epsu'
              else: newinv.name='epsl'
              newquot._invariants.append(newinv)
              placeholder.indices+=newInds
              newquot._invariants.remove(correspondingQuotInv)
            #numberedFields=othertFieldDict.values()
            if words: print 't0',newquot.longstr()
            for nf in othertFieldDict.values():
              for f in newquot._fields[:]:
                if id(f)==newifd[ipd[nf]]:
                  if words: print 'removing',ipd[nf],id(f)
                  newquot._fields.remove(f)
            for inv in QAllIn:
              if words: print 'removing inv',inv
              newquot._invariants.remove(getInv(inv,newquot))
            #for inv in QPartIn:
            #  newquot._invariants.remove(getInv(inv,newquot))
            if words: print 't1',indexConsistency(newquot)
            if words: print 't1',newquot.longstr()
            if words: print 'pniat1',pnia
            for inv,flipQinv in zip(NotIn,pnia):
              flip,Qinv=flipQinv
              newInds=[]
              correspondingQuotInv=getInv(Qinv,newquot)
              smallestNewLabel=newquot.smallestNewLabel(inv.group())
              #if not flip: inv.indices.reverse()
              for i in inv.indices:
                newi=Index(i.group,smallestNewLabel+len(newInds),i.upper,True)
                newInds.append(newi)
                othertToNewQuotIndexDict[i]=newi
              #if not flip: inv.indices.reverse()
              if not flip: newInds.reverse()
              if len(newInds)!=len(correspondingQuotInv.indices): raise Exception('ERROR: found wrong invariant in quot and don\'t know what\'s wrong')
              if correspondingQuotInv.name=='del' and inv.name=='del':
                correspondingQuotInv.indices.reverse()
                if not flip:
                  #newInds.reverse()
                  newisPos= not newisPos
              else:
                print 'WARNING: this code is untested on notin epsilons in quotient'
              for i0,i1 in zip(correspondingQuotInv.indices,newInds):
                newinv=Invariant('eps','\epsilon',[i0,i1] )
                if words: print 'newinv'#,len(newInds),len([i for i in correspondingQuotInv.indices if all( i not in f.indices for f in partialInvFields )])
                if newinv.indices[0].upper: newinv.name='epsu'
                else: newinv.name='epsl'
                newquot._invariants.append(newinv)
              if 'eps' in correspondingQuotInv.name and inv.name=='del':
                #find newinv which contains the first index of the original epsilon
                firstindexinv=newquot._invariants[-2]
                #if this inv is an epsilon, don't change sign, else do change sign
                if len(set(i.upper for i in firstindexinv.indices))!=1: newisPos = not newisPos
              placeholder.indices+=newInds
              newquot._invariants.remove(correspondingQuotInv)
            newquot.fixTwoIndexInvs()
            if words: print 't2a',indexConsistency(newquot)
            if words: print 't2a',newquot.longstr()
            ###sign stuff
            placeholder.indices=[ othertToNewQuotIndexDict[i] for i in othertSortedFreeIndices ]
            ###end sign stuff
            placeholder.explode()
            isGood=indexConsistency(newquot)
            if not isGood:
              print 't3',indexConsistency(newquot,True)
              print newquot.longstr()
              print self.longstr()
              print othert.longstr()
              raise Exception('ERROR: I made a quotient with inconsistent indices!')
            ####clear cache 
            newquot._clear()
            quotients.append((1 if newisPos else -1,newquot))
            ##if othert has no free indices, all other permutations should be indistinguishable
            if len(othertSortedFreeIndices)==0: break
            #if len(othertSortedFreeIndices)==0 or len(quotients)>5: break
    if len(quotients)<2: return quotients
    else:
      uniqueQs=[]
      for w,q in quotients:
        for wu,qu in uniqueQs:
          if qu==q:
            if qu.compareSign(q)*w==wu: break
            else: raise Exception('ERROR: quotients differ by sign. Awooga')
        else:
          uniqueQs.append( (w,q) )
      return uniqueQs
  def _compareCounters(self,c1,c2):
    return all(c1[k]<=c2[k] for k in c1)
  def _getCounter(self):
    return collections.Counter([f.getName() for f in self])
  def _getLinks(self,invariants):
    def getFields(inds):
      fields=[]
      for i in inds:
        for f in self:
          for no,dp in enumerate(f.Dindices):
            if i in dp:
              fields.append(str(no)+'_'+f.getName())
              break
          else:
            if i in f.indices: fields.append(f.getName())
      return fields
    links=[]
    for inv in invariants:
      links.append( (frozenset(getFields(inv.indices)),inv.name,inv.group() ) )
    return collections.Counter(links)
  def _getFieldDict(self):
    fieldDict={}
    counter=collections.Counter()
    #print 'start'
    for f in self:
      #print 'tick'
      fieldDict[id(f)]=f.getName()+'_'+str(counter[f.getName()])
      counter[f.getName()]+=1
    return fieldDict
  def _getNumberedLinks(self,invariants):
    def getFieldSets(inds,fDict):
      upperfields=[]
      lowerfields=[]
      for i in inds:
        f=next((fd for fd in self._fields if i in fd.indices),None)
        if f is None: continue
        for no,dp in enumerate(f.Dindices):
          if i in dp:
            if i.upper: upperfields.append(str(no)+'_'+fDict[id(f)])
            else: lowerfields.append(str(no)+'_'+fDict[id(f)])
            break
        else:
          if i.upper: upperfields.append(fDict[id(f)])
          else: lowerfields.append(fDict[id(f)])
      return frozenset(upperfields),frozenset(lowerfields)
    links=[]
    for inv in invariants:
      links.append( ( getFieldSets(inv.indices,self._getFieldDict()) ,inv.name,inv.group()) )
    return collections.Counter(links)
  def smallestNewLabel(self,group,dummy=False):
    indices=[]
    for f in self: indices+=f.indices
    for inv in self._invariants: indices+=inv.indices
    #indices=(self.freeIndices() if not dummy else self.dummyIndices())
    labels=[i.label for i in indices if i.group==group]
    if len(labels)>0: return max(labels)+1
    else: return 1
  def freeIndices(self):
    allIndices=[]
    for f in self: allIndices+=f.indices
    return [i for i in allIndices if not i.dummy]
  def dummyIndices(self):
    allIndices=[]
    for f in self: allIndices+=f.indices
    return [i for i in allIndices if i.dummy]
  def reindex(self,eombit):
    eomindices=[]
    for f in eombit: eomindices+=f.indices
    for inv in eombit._invariants: eomindices+=inv.indices
    eomindices=set(eomindices)
    selfindices=[]
    for f in self: selfindices+=f.indices
    for inv in self._invariants: selfindices+=inv.indices
    existingIndicesDict={group:[i.label for i in selfindices if i.group==group] for group in groupsDimFund.keys()}
    conversionDict={}
    for i in eomindices:
      #print 'ri i in eomindices',i
      if i.label in existingIndicesDict[i.group]:
        conversionDict[i]=next(j for j in range(1,100) if j not in existingIndicesDict[i.group])
        #print 'ri eiDict',existingIndicesDict
      else: conversionDict[i]=i.label
      #print 'ri conDict',conversionDict
      existingIndicesDict[i.group].append(conversionDict[i])
    return conversionDict
  def applyIndexConversion(self,conDict):
    self._clear()
    #selfindices=[]
    #for f in self: selfindices+=f.indices
    #for inv in self._invariants: selfindices+=inv.indices
    #selfindices=set(selfindices)
    #conversionTable
    #for i in selfindices: i.label=conDict[i]
    conversionTable=conDict.iteritems()
    for i,lbl in conversionTable: i.label=lbl
    return
  def add(self,f,force=False):
    self._clear()
    for no,i in enumerate(f.indices):
      if i in self.freeIndices()+f.indices[:no] and not i.dummy:
        i.label=max([ind.label for ind in self.freeIndices()+f.indices[:no] if ind.group==i.group])+1
      elif i.dummy:
        if force: pass
        else: raise Exception('help help dummy index')
    self.dimension+=f.dimension
    for u1 in f.U1Dict.keys():
      if u1 in self.uones.keys():
        self.uones[u1]+=f.U1Dict[u1]
      else:
        self.uones[u1]=f.U1Dict[u1]
    self._fields.append(f)
    f._term=self
    return
  def removeField(self,field):
    self._clear()
    i=self._fields.index(field)
    self.remove(i)
    return
  def remove(self,i):
    self._clear()
    self.dimension-=self._fields[i].dimension
    for u1 in self._fields[i].U1Dict.keys():
      self.uones[u1]-=self._fields[i].U1Dict[u1]
    self._fields[i]._term=None
    self._fields.pop(i)
    return
  def addTerm(self,t,force=False):
    self._clear()
    for f in t: self.add(f,force)
    for inv in t._invariants: self.addInvariant(inv,force)
    return
  def checkSyms(self):
    for f in self:
      for s in f.symmetries:
        if not s.sym(): return False
    if self.zeroBySym(): return False
    else: return True
  def checkSymsSimp(self):
    for f in self:
      for s in f.symmetries:
        if not s.sym(): return False
    return True
  def refreshUOnes(self):
    self._clear()
    self.uones={}
    for f in self:
      for u1 in f.U1Dict.keys():
        if u1 in self.uones.keys():
          self.uones[u1]+=f.U1Dict[u1]
        else:
          self.uones[u1]=f.U1Dict[u1]
  def fixTwoIndexInvs(self):
    for inv in self._invariants:
      if len(inv.indices)!=2: continue
      i0=inv.indices[0]
      i1=inv.indices[1]
      if i0.upper and i1.upper:
        inv.name='epsu'
        inv.tex=r'\epsilon'
        continue
      if i0.upper and not i1.upper:
        inv.name='del'
        inv.tex=r'\delta'
        continue
      if not i0.upper and i1.upper:
        inv.name='del'
        inv.tex=r'\delta'
        inv.indices.reverse()
        continue
      if not i0.upper and not i1.upper:
        inv.name='epsl'
        inv.tex=r'\epsilon'
        continue
    return
  def contractInvariantIndices(self):
    def removeSame():
      def eqInv(i1,i2):
        if id(i1)==id(i2): return False
        for indseq in itertools.permutations(i2.indices):
          if i1.indices==list(indseq): return True
        return False
      sgn=1
      for inv in self._invariants[:]:
        if inv not in self._invariants: continue
        match = next( (invtest for invtest in self._invariants if eqInv(inv,invtest) ),None)
        if match is not None:
          if 'eps' not in inv.name or 'eps' not in match.name:
            print self.longstr()
            raise Exception('ERROR: found two invariants with same indices but not of eps type')
          self._invariants.remove(inv)
          self._invariants.remove(match)
          sgn *= (-1 if self._cmpParity(inv.indices,match.indices) else 1) * math.factorial(len(inv.indices))
      return sgn
    def getOtherOne(invar,index):
      return next(i for i in invar.indices if not i==index)
    def makeNewTwoInv(inv1,inv2):
      commonIndices=set(inv1.indices).intersection(set(inv2.indices))
      if len(commonIndices)!=1:
        print self.longstr()
        raise Exception('ERROR: contracting two invariants without one common index')
      commonIndex=commonIndices.pop()
      inv1freeI=getOtherOne(inv1,commonIndex)
      inv2freeI=getOtherOne(inv2,commonIndex)
      if inv1.name=='epsl':
        if inv2.name=='epsl': raise Exception('ERROR: Beaver convention broken')
        elif inv2.name=='del': newname='epsl'
        elif inv2.name=='epsu': newname='del'
      if inv1.name=='epsu':
        if inv2.name=='epsu': raise Exception('ERROR: Beaver convention broken')
        elif inv2.name=='del': newname='epsu'
        elif inv2.name=='epsl': newname='del'
      if inv1.name=='del':
        if inv2.name=='del': newname='del'
        elif inv2.name=='epsl': newname='epsl'
        elif inv2.name=='epsu': newname='epsu'
      newsign=(-1 if 'eps' in inv1.name and commonIndex==inv1.indices[0] else 1)*(-1 if 'eps' in inv2.name and commonIndex==inv2.indices[1] else 1)
      return newsign,Invariant(newname,r'\epsilon' if 'eps' in newname else r'\delta',[inv1freeI,inv2freeI])
    def makeNewMultiInv(inv1,inv2):
      raise Exception('ERROR: I have no idea why this was called')
    self._clear()
    fieldIndices=list(itertools.chain.from_iterable(f.indices for f in self))
    sgn=1
    sgn*=removeSame()
    while not all(all(i in fieldIndices for i in inv.indices) for inv in self._invariants if len(inv.indices)==2):
      ##find a pair inv1 and inv2 which share an index
      for inv1 in [iv for iv in self._invariants if len(iv.indices)==2]:
        inv2=next((iv for iv in self._invariants if id(iv)!=id(inv1) and any(i in inv1.indices for i in iv.indices)),None)
        if inv2 is None: continue
        else: break
      if inv2 is None: break
      ##make a new one by contracting indices
      if len(inv2.indices)==2: newsgn,newinv=makeNewTwoInv(inv1,inv2)
      else: newsgn,newinv=makeNewMultiInv(inv1,inv2)
      sgn*=newsgn
      self._invariants.remove(inv1)
      self._invariants.remove(inv2)
      self._invariants.append(newinv)
      sgn*=removeSame()
    return sgn
    unchanged=False
    sgn=1
    while not unchanged:
      unchanged=True
      sgn*= removeSame()
      for inv in [i for i in self._invariants if len(i.indices)==2]:
        if not unchanged: break
        if inv not in self._invariants: continue
        match0=next((iv for iv in self._invariants if inv.indices[0] in iv.indices and id(iv)!=id(inv) and len(iv.indices)==2 ),None)
        if match0 is not None:
          newi1=getOtherOne(match0,inv.indices[0])
          newi2=getOtherOne(inv,inv.indices[0])
          self._invariants.remove(match0)
          self._invariants.remove(inv)
          self._invariants.append( makeNewTwoIndexInvar(newi1,newi2,inv.name,match0.name) )
          if 'eps' in match0.name and 'eps' in inv.name: sgn *= (1 if match0.indices[1]==inv.indices[0] else -1)
          unchanged=False
          continue
        match1=next((iv for iv in self._invariants if inv.indices[1] in iv.indices and id(iv)!=id(inv) and len(iv.indices)==2 ),None)
        if match1 is not None:
          newi1=getOtherOne(match1,inv.indices[1])
          newi2=getOtherOne(inv,inv.indices[1])
          self._invariants.remove(match1)
          self._invariants.remove(inv)
          self._invariants.append( makeNewTwoIndexInvar(newi2,newi1,inv.name,match1.name) )
          if 'eps' in match1.name and 'eps' in inv.name: sgn *= (1 if match1.indices[0]==inv.indices[1] else -1)
          unchanged=False
          continue
      for inv in [i for i in self._invariants if len(i.indices)!=2]:
        if inv not in self._invariants: continue
        for ind in inv.indices:
          match=next((iv for iv in self._invariants if id(iv)!=id(inv) and ind in iv.indices and iv.name=='del'),None)
          if match is not None:
            indloc=inv.indices.index(ind)
            inv.indices.pop(indloc)
            inv.indices.insert(indloc, getOtherOne(match,ind) )
            unchanged=False
            continue
    return sgn
  #BCTI compliant 2016-02-24
  ###add reverse fields 2016-04-22
  def conjugate(self):
    newt=self._DEFTcopy()
    for f in newt:
      f.conjugate()
    for inv in newt._invariants:
      if inv.group() in ['lorL','lorR']: continue
      if inv.name=='del': inv.indices.reverse()
      else:
        if inv.name=='epsu': inv.name='epsl'
        else: inv.name='epsu'
        ##for sign convention eps^{abcd...} = -eps_{abcd...}
        x,y=inv.indices[0],inv.indices[1]
        inv.indices[0],inv.indices[1]=y,x
    newt._fields.reverse()
    newt._clear()
    return newt
  def CPconjugate(self):
    newt=self.conjugate()
    newt._fields.reverse()
    for f in newt:
      for j,i in enumerate(f.indices[:]):
        if 'lor' in i.group:
          newlabel=newt.smallestNewLabel(i.group)
          newindex=Index(i.group,newlabel,i.upper,True)
          newinvariant= Invariant('epsu' if i.upper else 'epsl',r'\epsilon',[i,newindex])
          f.replaceIndexInField(i,newindex)
          f.replaceIndexInSym(i,newindex)
          newt._invariants.append(newinvariant)
    sgn=newt.contractInvariantIndices()
    return sgn,newt
  def removeOnes(self,term):
    newfields=[f for f in term if f.getName()!='1']
    term._clear()
    term._fields=newfields
    return term
  def _patchUpInvariantIndices(self):
    for f in self:
      for i in f.indices:
        inv=next((iv for iv in self._invariants if i in iv.indices),None)
        if inv is None: continue
        else:
          indpos=inv.indices.index(i)
          inv.indices.pop(indpos)
          inv.indices.insert(indpos,i)
    for inv in self._invariants:
      for i in inv.indices:
        inv2=next((iv for iv in self._invariants if i in iv.indices and id(iv)!=id(inv)),None)
        if inv2 is None: continue
        else:
          indpos=inv2.indices.index(i)
          inv2.indices.pop(indpos)
          inv2.indices.insert(indpos,i)
    return
  def _contractDeltasGen(self):
    self._clear()
    fieldIndices=[]
    for f in self: fieldIndices+=f.indices
    for inv in self._invariants[:]:
      if len(inv.indices)!=2 or inv.indices[0].upper!=inv.indices[1].upper or inv.name!='del' or all(i in fieldIndices for i in inv.indices) or not any(i in fieldIndices for i in inv.indices): continue
      #if len(inv.indices)!=2 or inv.indices[0].upper!=inv.indices[1].upper or inv.indices[0].dummy==inv.indices[1].dummy: continue
      freeOne=next(i for i in inv.indices if i not in fieldIndices)
      dummyOne=next(i for i in inv.indices if i in fieldIndices)
      field=next(f for f in self if dummyOne in f.indices)
      field.replaceIndexInField(dummyOne,freeOne)
      field.replaceIndexInSym(dummyOne,freeOne)
      self._invariants.remove(inv)
    return
  def _contractDeltas(self):
    self._clear()
    for inv in self._invariants[:]:
      if len(inv.indices)!=2 or inv.indices[0].upper!=inv.indices[1].upper or inv.indices[0].dummy==inv.indices[1].dummy or inv.name!='del': continue
      #if len(inv.indices)!=2 or inv.indices[0].upper!=inv.indices[1].upper or inv.indices[0].dummy==inv.indices[1].dummy: continue
      freeOne=next(i for i in inv.indices if not i.dummy)
      dummyOne=next(i for i in inv.indices if i.dummy)
      field=next(f for f in self if dummyOne in f.indices)
      field.replaceIndexInField(dummyOne,freeOne)
      field.replaceIndexInSym(dummyOne,freeOne)
      self._invariants.remove(inv)
    return
  def _fixDummyLabels(self):
    self._clear()
    fieldIndices=[]
    for f in self: fieldIndices+=f.indices
    for inv in self._invariants:
      for i in inv.indices:
        if (i in fieldIndices and not i.dummy) or (i not in fieldIndices and i.dummy):
          i.dummy=not i.dummy
          i.label=self.smallestNewLabel(i.group)
    return
  def _IBP(self,term):
    nD=sum(len(f.Dindices) for f in term)
    prods=itertools.product(range(len(term._fields)), repeat=nD)
    newTerms=[]
    for p in prods:
      t=term._DEFTcopy()
      dIndices=[]
      for f in t:
        dIndices+=f.Dindices
        for di in f.Dindices:
          f.remove(di[1])
          f.remove(di[0])
        f.Dindices=[]
      for i in p:
        dis=dIndices.pop()
        t[i].indices.append(dis[0])
        t[i].indices.append(dis[1])
        t[i].Dindices.append(dis)
      newTerms.append(t)
    return newTerms
  def funcDiff(self,field):
    listOfTerms=[]
    balls=False
    #balls=True
    if balls: print 'termfuncdiff0aaaaa',self.shortstr()
    if balls: print 'termfuncdiff0aaaaa',self.longstr()
    for i,f in enumerate(self):
      try: newsubts=f.funcDiff(field)
      except KeyError: pass
      else:
        for wgt,subt in newsubts:
          newt=self._DEFTcopy()
          newt.remove(i)
          if balls: print 'termfuncdiff0a',subt.longstr()
          #if len(subt._fields)==0: pass
          #else: newt.addTerm(subt,True)
          newt.addTerm(subt,True)
          if balls: print 'termfuncdiff0b',subt.longstr()
          #rearrange
          for n in range(len(subt._fields)):
            newt._fields.insert(i, newt._fields.pop() )
          if balls: print 'termfuncdiff1', newt.longstr()
          sgn=newt.contractInvariantIndices()
          if balls: print 'termfuncdiff2', newt.longstr()
          if balls: print 'termfuncdiff3', newt.longstr()
          newt._patchUpInvariantIndices()
          newt._contractDeltas()
          newt.refreshUOnes()
          if balls: print 'termfuncdiff4',newt.longstr()
          if not indexConsistency(newt): raise Exception('New func diff term has inconsistent indices')
          if not any(f._basename=='1' and len(f.Dindices)!=0 for f in newt):
            listOfTerms.append((wgt*sgn,newt))
          else:
            assert sum(1 if f._basename=='1' else 0 for f in newt)==1
            evenDeriv=(len(next(f.Dindices for f in newt if f._basename=='1')) % 2 ==0)
            ibpterms=self._IBP(newt)
            if balls: print [t.shortstr() for t in ibpterms]
            listOfTerms+=[(wgt*sgn*(1 if evenDeriv else -1),t) for t in ibpterms if not any(f._basename=='1' and len(f.Dindices)!=0 for f in t)]
    eom=[(w,self.removeOnes(t)) for w,t in listOfTerms]
    canonicaleom=[]
    for w,t in eom:
      isPos=True
      for inv in t._invariants:
        if inv.name=='del' or all(i.dummy for i in inv.indices): continue
        oldindices=inv.indices
        inv.indices=sorted([i for i in oldindices if not i.dummy],key=str)+[i for i in oldindices if i.dummy]
        if not self._cmpParity(inv.indices,oldindices): isPos = not isPos
      canonicaleom.append( (w*(1 if isPos else -1),t) )
    return canonicaleom
  def addInvariant(self,inv,force=False):
    self._clear()
    for i in inv.indices:
      if not force: i.dummy=True
      if self.dummyIndices().count(i)>1: i.label=self.smallestNewLabel(i.group,True)
    self._invariants.append(inv)
  def removeInvariant(self,inv):
    self._clear()
    self._invariants.remove(inv)
    for i in inv.indices:
      i.dummy=False
  def __str__(self):
    return '*****\n'+' '.join([f.getName() for f in self])+'\n*****\nfree field indices\n'+'\n'.join([str(i) for i in self.freeIndices()])+'\ndummy field indices\n'+'\n'.join([str(i) for i in self.dummyIndices()])+'\nU(1) charges\n'+'\n'.join([k+' '+str(self.uones[k]) for k in self.uones.keys()])
  def __repr__(self):
    #def nextIntGen():
    #  n=1
    #  while True:
    #    yield n
    #    n+=1
    #intGen=nextIntGen()
    iStrList=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    iStrList.reverse()
    indexKey={}
    invStrs=[]
    for inv in self._invariants:
      invStr=inv.name[0]+'['
      for i in inv.indices:
        indexKey[i]=iStrList.pop()
        invStr+=indexKey[i]
      invStr+=']'
      invStrs.append(invStr)
    fStrs=[]
    for f in self:
      for i in f.indices:
        try: indexKey[i]
        except ValueError: indexKey[i]=iStrList.pop()
      diflat=[]
      fStr=''
      for di in f.Dindices:
        fStr+='D['+indexKey[di[0]]+indexKey[di[1]]+']'
        diflat+=[di[0],di[1]]
      fStr+=f._basename+('!' if f.isConjugate else '')+'['
      for i in f.indices:
        if i in diflat: continue
        fStr+=indexKey[i]
      fStr+=']'
      fStrs.append(fStr)
    return ' '.join(invStrs)+(' ' if len(invStrs)>0 else '')+' '.join(fStrs)
  @staticmethod
  def unrepr(reprstr,listOfFields):
    strs=reprstr.split(' ')
    invstrs=[s for s in strs if s.split('[')[0] in 'de']
    fieldstrs=[s for s in strs if s.split('[')[0] not in 'de']
    term=Term()
    indexKey={}
    for fs in fieldstrs:
      bits=fs.replace(']','[').split('[')[:-1]
      isConjugate=('!'==bits[-2][-1])
      basename=bits[-2].strip('!')
      fclass=next(f for f in listOfFields if f._basename==basename)
      f=fclass()
      if isConjugate: f.conjugate()
      refIndices=f.indices
      if len(refIndices)!=len(bits[-1]): raise Exception('Field {0} has incompatible number of indices'.format(basename))
      f.indices=[]
      term.add(f)
      for i,lett in enumerate(bits[-1]):
        indexKey[lett]=Index(refIndices[i].group,term.smallestNewLabel(refIndices[i].group),refIndices[i].upper,True)
        f.indices.append(indexKey[lett])
        f.replaceIndexInSym(refIndices[i],indexKey[lett])
      for lettpair in [bit for i,bit in enumerate(bits[:-2]) if i%2==1]:
        Llett,Rlett=lettpair[0],lettpair[1]
        indexKey[Llett]=Index('lorL',term.smallestNewLabel('lorL'),False,True)
        indexKey[Rlett]=Index('lorR',term.smallestNewLabel('lorR'),False,True)
        f.dimension+=frac(1,1)
        f.Dindices.append( (indexKey[Llett],indexKey[Rlett]) )
        f.indices+=[indexKey[Llett],indexKey[Rlett]]
    for ivs in invstrs:
      iname,iletts=ivs.strip(']').split('[')
      term._invariants.append( Invariant(('epsu' if indexKey[iletts[0]].upper else 'epsl') if iname=='e' else 'del','\epsilon' if iname=='e' else '\delta',[indexKey[l] for l in iletts]) )
    if not indexConsistency(term): raise Exception('indices of constructed term are inconsistent')
    return term 
  def shortstr(self):
    return ' '.join([f.getName() for f in self])
  def longstr(self):
    return self.__str__()+'\nfields\n'+'\n'.join(f.getName()+' '+str(id(f))+'\n  '+'\n  '.join(str(i)+' '+str(id(i)) for i in f.indices) for f in self)+'\ninvariants\n'+'\n'.join(inv.name+'\n  '+'\n  '.join(str(i)+' '+str(id(i)) for i in inv.indices) for inv in self._invariants)
  def _fieldTexRep(self,field,indexKey):
    derivativeTerm=''
    for di in field.Dindices:
      try:
        indexBit=indexKey[di[0]]+' '+indexKey[di[1]]
      except KeyError:
        derivativeTerm+='D '
        continue
      derivativeTerm+='D'+('^{' if di[0].upper else '_{')
      derivativeTerm+=indexBit
      derivativeTerm+='} '
      
    fieldTerm=derivativeTerm+(r'{\bar{'+field.tex+'}}' if field.isConjugate else '{'+field.tex+'}')
    up=''
    down=''
    dindicesflat=[di[0] for di in field.Dindices]+[di[1] for di in field.Dindices]
    for i in field.indices:
      if i in dindicesflat: continue
      try:
        newbit=indexKey[i]+' '
      except KeyError: continue
      if i.upper: up+=newbit
      else: down+=newbit
    return fieldTerm+'^{'+up+'}_{'+down+'}'
  def _invTexRep(self,inv,indexKey):
    up=''
    down=''
    for i in inv.indices:
      try:
        newbit=indexKey[i]+' '
      except KeyError:
        newbit=''
      if i.upper: down+=newbit
      else: up+=newbit
    return inv.tex+'^{'+up+'}_{'+down+'}'
  def textRep(self,suppress=[]):
    from operator import add as opadd
    indexDict={'su3':['A','B','C','D','E','F'],
'su4':['A','B','C','D','E','F','H','I','J','K'],
'su2':['a','b','d','e','f','g','h','i','j','k','l','m','n','o','p'],
'sug':['p','q','r','s','t','u','v'],
'lorL':['1','2','3','4','5','6','7','8'],
'lorR':['!','"','#','$','%','&','*','+']}
    indexKey={}
    for inv in self._invariants:
      if inv.indices[0].group in suppress: continue
      if inv.name=='del':
        indexKey[inv.indices[0]]=(str(inv.indices[0].value) if inv.indices[0].value is not None else indexDict[inv.group()].pop(0))
        indexKey[inv.indices[1]]=(str(inv.indices[1].value) if inv.indices[1].value is not None else indexKey[inv.indices[0]] )
      else:
        for i in inv.indices:
          indexKey[i]=(str(i.value) if i.value is not None else indexDict[inv.group()].pop(0) )
    for f in self:
      for i in f.indices:
        if not i.dummy: indexKey[i]=indexDict[i.group].pop(0)
    topLine=''
    middleLine=''
    bottomLine=''
    lines=[topLine,middleLine,bottomLine]
    for inv in self._invariants:
      if inv.name=='del': continue
      lines=map(opadd,lines,[' ','e',' '])
      for i in inv.indices:
        lines=map(opadd,lines,([' ',' ',indexKey[i]] if i.upper else [indexKey[i],' ',' ']))
    #if len(lines[0])>0: lines=map(opadd,lines,[' ',' ',' '])
    for f in self:
      for di in f.Dindices:
        lines=map(opadd,lines,[' ','D',' '])
        for i in di:
          lines=map(opadd,lines,([' ',' ',indexKey[i]] if not i.upper else [indexKey[i],' ',' ']))
      nm=f._name
      dindicesflat=[di[0] for di in f.Dindices]+[di[1] for di in f.Dindices]
      lines=map(opadd,lines,[' '*len(nm),nm,' '*len(nm)])
      for i in f.indices:
        if i in dindicesflat: continue
        lines=map(opadd,lines,([' ',' ',indexKey[i]] if not i.upper else [indexKey[i],' ',' ']))
    return lines[0]+'\n'+lines[1]+'\n'+lines[2]+'\n'
    #return '$'+' '.join([self._invTexRep(inv,indexKey) for inv in self._invariants if inv.name!='del' and inv.indices[0].group not in suppress])+' '+' '.join([self._fieldTexRep(f,indexKey) for f in self])+'$'

  def texRep(self,suppress=[]):
    indexDict={'su3':['A','B','C','D','E','F','H','I','J','K'],
'su4':['A','B','C','D','E','F','H','I','J','K'],
'su2':['a','b','d','e','f','g','h','i','j','k','l','m','n','o','p'],
'sug':['p','q','r','s','t','u','v'],
'lorL':[r'\alpha',r'\beta',r'\gamma',r'\delta',r'\mu',r'\nu',r'\rho',r'\sigma'],
'lorR':[r'\dot\alpha',r'\dot\beta',r'\dot\gamma',r'\dot\delta',r'\dot\mu',r'\dot\nu',r'\dot\rho',r'\dot\sigma']}
    indexKey={}
    for inv in self._invariants:
      if inv.indices[0].group in suppress: continue
      if inv.name=='del':
        indexKey[inv.indices[0]]=(str(inv.indices[0].value) if inv.indices[0].value is not None else indexDict[inv.group()].pop(0) )
        indexKey[inv.indices[1]]=(str(inv.indices[1].value) if inv.indices[1].value is not None else indexKey[inv.indices[0]] )
      else:
        for i in inv.indices:
          indexKey[i]=(str(i.value) if i.value is not None else indexDict[inv.group()].pop(0) )
    for f in self:
      for i in f.indices:
        if not i.dummy: indexKey[i]=indexDict[i.group].pop(0)
    return '$'+' '.join([self._invTexRep(inv,indexKey) for inv in self._invariants if inv.name!='del' and inv.indices[0].group not in suppress])+' '+' '.join([self._fieldTexRep(f,indexKey) for f in self])+'$'
  def _DEFTcopy(self):
    newt=Term()
    newt._fields=[f._DEFTcopy() for f in self]
    for f in newt._fields: f._term=newt
    newt.dimension=self.dimension._DEFTcopy()
    newt.uones=dict((k,v._DEFTcopy()) for k,v in self.uones.iteritems())
    newt._hash=self._hash
    newindices=[]
    for f in newt: newindices+=f.indices
    newt._invariants=[]
    for inv in self._invariants:
      newinv,extrainds=inv._DEFTcopy(newindices)
      newt._invariants.append(newinv)
      newindices+=extrainds
    return newt
  def _clear(self):
    self._hash=None
    self._mash=None
    self._fieldDict=None
    self._pds=None
    self._pnls=None
    return
  def __init__(self,listOfFields=[]):
    self._fields=[]
    self.dimension=frac(0,1)
    self.uones=dict()
    self._invariants=[]
    self._hash=None
    self._mash=None
    self._fieldDict=None
    self._pds=None
    self._pnls=None
    for f in listOfFields:
      self.add(f)


def indexConsistency(trm,pr=False):
  def DindicesInIndices(term):
    for f in term:
      for pr in f.Dindices:
        if pr[0] not in f.indices: return False
        if pr[1] not in f.indices: return False
    return True
  def fieldIndicesUnique(term):
    fis=[]
    for f in term: fis+=f.indices
    return len(fis)==len(set(fis))
  def InvIndicesUnique(term):
    iis=[]
    for inv in term._invariants: iis+=inv.indices
    return len(iis)==len(set(iis))
  def InvIndicesIdsUnique(term):
    iis=[]
    for inv in term._invariants: iis+=inv.indices
    for i in iis:
      if any(j==i and id(j)!=id(i) for j in iis): return False
    return True
  def InvIndicesSameGroup(term):
    sets=[set(i.group for i in inv.indices) for inv in term._invariants]
    return all(len(s)==1 for s in sets)
  def FisFandDisD(term):
    fis=[]
    for f in term: fis+=f.indices
    iis=[]
    for inv in term._invariants: iis+=inv.indices
    for fi in fis:
      if fi in iis and not fi.dummy: return False
      if fi not in iis and fi.dummy: return False
    for ii in iis:
      if ii in fis and not ii.dummy: return False
      if ii not in fis and ii.dummy: return False
    return True
  def InvIndexIds(term):
    fieldInds=list(itertools.chain.from_iterable(f.indices for f in term))
    fieldInds=[i for i in fieldInds if i.dummy]
    fieldIds=set(id(i) for i in fieldInds)
    invInds=list(itertools.chain.from_iterable(inv.indices for inv in term._invariants))
    invIds=set(id(i) for i in invInds if i.dummy)
    return fieldIds.issubset(invIds)
  def SymIndexIds(term):
    for f in term:
      fids=[id(ind) for ind in f.indices]
      #print 'sii',id(f)
      #print 'sii',fids
      for sym in f.symmetries:
        if any(id(ind) not in fids for ind in sym.indices): return False
    return True
  def delOrder(term):
    return all(inv.indices[0].upper and not inv.indices[1].upper for inv in term._invariants if inv.name=='del')
  tests=[DindicesInIndices,fieldIndicesUnique,InvIndicesUnique,InvIndicesIdsUnique,InvIndicesSameGroup,FisFandDisD,InvIndexIds,SymIndexIds,delOrder]
  results=[]
  for test in tests:
    result=test(trm)
    if pr: print test.__name__+'\t'+('True' if result else 'False')
    results.append(result)
  return all( results )



import fields# import B,W,G
gaugeGroupsFields={'Y':fields.B,'su2':fields.W,'su3':fields.G}

deftstr=r'''
 ______   _______  _______ _________
(  __  \ (  ____ \(  ____ \\__   __/
| (  \  )| (    \/| (    \/   ) (   
| |   ) || (__    | (__       | |   
| |   | ||  __)   |  __)      | |   
| |   ) || (      | (         | |   
| (__/  )| (____/\| )         | |   
(______/ (_______/|/          )_(   
  does effective field theory(ish)  '''
