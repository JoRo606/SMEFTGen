#from DEFT.operators import Field
from DEFT import Field,frac,Unitary

class Gs3(Field):
  _basename='gs3'
  def __init__(self,term=None):
    Field.__init__(self,frac(0,1),{'sug':'3','su3':'1','su2':'1','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'gs3','gs3',term)
    self.irrepDim=3
    #self.makeIndices()

class Gs6(Field):
  _basename='gs6'
  def __init__(self,term=None):
    Field.__init__(self,frac(0,1),{'sug':'6','su3':'1','su2':'1','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'gs6','gs6',term)
    self.irrepDim=6
    #self.makeIndices()

class Gs8(Field):
  _basename='gs8'
  def __init__(self,term=None):
    Field.__init__(self,frac(0,1),{'sug':'8','su3':'1','su2':'1','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'gs8','gs8',term)
    self.isSelfConjugate=True
    self.irrepDim=8
    #self.makeIndices()

class Gs10(Field):
  _basename='gs10'
  def __init__(self,term=None):
    Field.__init__(self,frac(0,1),{'sug':'10','su3':'1','su2':'1','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'gs10','gs10',term)
    self.irrepDim=10
    #self.makeIndices()

class Gs15p(Field):
  _basename='gs15p'
  def __init__(self,term=None):
    Field.__init__(self,frac(0,1),{'sug':'15p','su3':'1','su2':'1','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'gs15p','gs15p',term)
    self.irrepDim=15
    #self.makeIndices()

class Gs15(Field):
  _basename='gs15'
  def __init__(self,term=None):
    Field.__init__(self,frac(0,1),{'sug':'15','su3':'1','su2':'1','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'gs15','gs15',term)
    self.irrepDim=15
    #self.makeIndices()

class Gs24(Field):
  _basename='gs24'
  def __init__(self,term=None):
    Field.__init__(self,frac(0,1),{'sug':'24','su3':'1','su2':'1','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'gs24','gs24',term)
    self.irrepDim=24
    #self.makeIndices()

class Gs27(Field):
  _basename='gs27'
  def __init__(self,term=None):
    Field.__init__(self,frac(0,1),{'sug':'27','su3':'1','su2':'1','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'gs27','gs27',term)
    self.isSelfConjugate=True
    self.irrepDim=27
    #self.makeIndices()

class ER(Field):
  _basename='eR'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,1),'L':frac(1,1)},'eR','e_R',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'

class ER1(Field):
  _basename='eR1'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,1)},'eR1','e_R^1',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'

class ER2(Field):
  _basename='eR2'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,1)},'eR2','e_R^2',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'

class ER3(Field):
  _basename='eR3'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,1)},'eR3','e_R^3',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'

class ERi(Field):
  _basename='eRi'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'sug':'3','su3':'1','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,1)},'eRi','e_R',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'

class LLS(Field):
  _basename='lLS'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'1','lorL':'2','lorR':'1'},{'Y':frac(-1,2)},'lLS','LS_L',term)
    ##self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'


class LLY(Field):
  _basename='lLY'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(0,1),'L':frac(1,1)},'lLY','L_L',term)
    #self.makeIndices()
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'

class LL(Field):
  _basename='lL'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(-1,2),'L':frac(1,1)},'lL','L_L',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'

class LL1(Field):
  _basename='lL1'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(-1,2)},'lL1','L_L^1',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'

class LL2(Field):
  _basename='lL2'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(-1,2)},'lL2','L_L^2',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'

class LL3(Field):
  _basename='lL3'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(-1,2)},'lL3','L_L^3',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'

class LLi(Field):
  _basename='lLi'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'sug':'3','su3':'1','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(-1,2)},'lLi','L_L',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'

class QB(Field):
  _basename='qB'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(1,6),'B':frac(1,3)},'qB','Q_B',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'

class QL(Field):
  _basename='qL'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(1,6),'B':frac(1,3)},'qL','Q_L',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class QL4(Field):
  _basename='qL4'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su4':'4','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(1,6),'B':frac(1,3)},'qL4','Q_L',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'

class QL1(Field):
  _basename='qL1'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(1,6)},'qL1','Q_L^1',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class QL2(Field):
  _basename='qL2'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(1,6)},'qL2','Q_L^2',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class QL3(Field):
  _basename='qL3'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(1,6)},'qL3','Q_L^3',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class QLi(Field):
  _basename='qLi'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'sug':'3','su3':'3','su2':'2','lorL':'2','lorR':'1'},{'Y':frac(1,6)},'qLi','Q_L',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class URY(Field):
  _basename='uRY'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(0,1),'B':frac(1,3)},'uRY','u_R',term)
    #self.makeIndices()
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class UB(Field):
  _basename='uB'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(2,3),'B':frac(1,3)},'uB','u_B',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'

class UR(Field):
  _basename='uR'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(2,3),'B':frac(1,3)},'uR','u_R',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class UR4(Field):
  _basename='uR4'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su4':'4','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(2,3),'B':frac(1,3)},'uR4','u_R',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'

class UR1(Field):
  _basename='uR1'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(2,3)},'uR1','u_R^1',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class UR2(Field):
  _basename='uR2'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(2,3)},'uR2','u_R^2',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class UR3(Field):
  _basename='uR3'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(2,3)},'uR3','u_R^3',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class URi(Field):
  _basename='uRi'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'sug':'3','su3':'3','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(2,3)},'uRi','u_R',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class DB(Field):
  _basename='dB'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,3),'B':frac(1,3)},'dB','d_B',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'

class DR(Field):
  _basename='dR'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,3),'B':frac(1,3)},'dR','d_R',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class DR4(Field):
  _basename='dR4'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su4':'4','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,3),'B':frac(1,3)},'dR4','d_R',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'

class DR1(Field):
  _basename='dR1'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,3)},'dR1','d_R^1',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class DR2(Field):
  _basename='dR2'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,3)},'dR2','d_R^2',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class DR3(Field):
  _basename='dR3'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'3','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,3)},'dR3','d_R^3',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class DRi(Field):
  _basename='dRi'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'sug':'3','su3':'3','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,3)},'dRi','d_R',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AG']='self._FDNA(field,\'su3\')'

class H(Field):
  _basename='H'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'2','lorL':'1','lorR':'1'},{'Y':frac(1,2)},'H','H',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'

class h(Field):
  _basename='h'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'1','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'h','h',term)
    self.isSelfConjugate=True
    #self.makeIndices()

class U(Field):
  _basename='U'
  def __init__(self,term=None):
    Field.__init__(self,frac(0,1),{'su3':'1','su2':'2','lorL':'1','lorR':'1'},{'Y':frac(1,2)},'U','U',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AW']='self._FDNA(field,\'su2\')'
    self.symmetries.append(Unitary(self,[self.indices[0]]))

class S(Field):
  _basename='S'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'1','lorL':'1','lorR':'1'},{'Y':frac(1,2),'Z':frac(0,1)},'S','S',term)
    #self.makeIndices()
#    self.funcDiffDict['AB']='self._FDU1(field,\'Y\')'
#    self.funcDiffDict['AF']='self._FDU1(field,\'Z\')'

class phi(Field):
  _basename='S'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'1','lorL':'1','lorR':'1'},{},'S','\phi',term)
    self.isSelfConjugate=True
    #self.makeIndices()

class phiY(Field):
  _basename='S'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'1','lorL':'1','lorR':'1'},{'Y':frac(1,2)},'S','\phi',term)
    #self.makeIndices()

class psi(Field):
  _basename='F'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'1','lorL':'2','lorR':'1'},{},'F','\psi',term)
    #self.makeIndices()

class psiY(Field):
  _basename='F'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'1','lorL':'2','lorR':'1'},{'Y':frac(-1,2)},'F','\psi',term)
    #self.makeIndices()

class psiYR(Field):
  _basename='E'
  def __init__(self,term=None):
    Field.__init__(self,frac(3,2),{'su3':'1','su2':'1','lorL':'1','lorR':'2'},{'Y':frac(-1,1)},'E','\psi_R',term)
    #self.makeIndices()

class AB(Field):
  _basename='AB'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'1','lorL':'2l','lorR':'2l'},{'Y':frac(0,1)},'AB','AB',term)
    self.isSelfConjugate=True
    #self.makeIndices()

class B(Field):
  _basename='B'
  _vectorPotential=AB
  def __init__(self,term=None):
    Field.__init__(self,frac(2,1),{'su3':'1','su2':'1','lorL':'3','lorR':'1'},{'Y':frac(0,1),'Z':frac(0,1)},'B','B',term)
    #self.isSelfConjugate=True
    #self.makeIndices()
    self.funcDiffDict['AB']='self._FDFS(field)'

class AF(Field):
  _basename='AF'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'1','lorL':'2l','lorR':'2l'},{'Z':frac(0,1)},'AF','AF',term)
    self.isSelfConjugate=True
    #self.makeIndices()

class F(Field):
  _basename='F'
  _vectorPotential=AF
  def __init__(self,term=None):
    Field.__init__(self,frac(2,1),{'su3':'1','su2':'1','lorL':'3','lorR':'1'},{'Y':frac(0,1),'Z':frac(0,1)},'F','F',term)
    #self.isSelfConjugate=True
    #self.makeIndices()
    self.funcDiffDict['AF']='self._FDFS(field)'

class AW(Field):
  _basename='AW'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'3','lorL':'2l','lorR':'2l'},{'Y':frac(0,1)},'AW','AW',term)
    self.isSelfConjugate=True
    #self.makeIndices()

class W(Field):
  _basename='W'
  _vectorPotential=AW
  def __init__(self,term=None):
    Field.__init__(self,frac(2,1),{'su3':'1','su2':'3','lorL':'3','lorR':'1'},{'Y':frac(0,1)},'W','W',term)
    #self.isSelfConjugate=True
    #self.makeIndices()
#    self.funcDiffDict['AW']='self._FDNAFS(field,\'su2\')'

class Phi(Field):
  _basename='P'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'7','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'P','\Phi',term)
    #self.isSelfConjugate=True
    #self.makeIndices()

class AG(Field):
  _basename='AG'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'8','su2':'1','lorL':'2l','lorR':'2l'},{'Y':frac(0,1)},'AG','AG',term)
    self.isSelfConjugate=True
    #self.makeIndices()

class G(Field):
  _basename='G'
  _vectorPotential=AG
  def __init__(self,term=None):
    Field.__init__(self,frac(2,1),{'su3':'8','su2':'1','lorL':'3','lorR':'1'},{'Y':frac(0,1)},'G','G',term)
    #self.isSelfConjugate=True
    #self.makeIndices()
#    self.funcDiffDict['AG']='self._FDNAFS(field,\'su3\')'

class AG4(Field):
  _basename='AG4'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su4':'15','su2':'1','lorL':'2l','lorR':'2l'},{'Y':frac(0,1)},'AG4','AG4',term)
    self.isSelfConjugate=True
class G4(Field):
  _basename='G4'
  _vectorPotential=AG4
  def __init__(self,term=None):
    Field.__init__(self,frac(2,1),{'su4':'15','su2':'1','lorL':'3','lorR':'1'},{'Y':frac(0,1)},'G4','G4',term)

class S3R(Field):
  _basename='S3R'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'3','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'S3R','S_3',term)
    self.isSelfConjugate=True
class S3C(Field):
  _basename='S3C'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'3','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'S3C','S_3',term)

class S1(Field):
  _basename='S1'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'1','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'S1','S_1',term)
    self.isSelfConjugate=True
    #self.makeIndices()
class S2(Field):
  _basename='S2'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'2','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'S2','S_2',term)
    #self.makeIndices()
class S3(Field):
  _basename='S3'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'3','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'S3','S_3',term)
    #self.makeIndices()
class S4(Field):
  _basename='S4'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'4','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'S4','S_4',term)
    #self.makeIndices()
class S5(Field):
  _basename='S5'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'5','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'S5','S_5',term)
    #self.makeIndices()
class S6(Field):
  _basename='S6'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'6','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'S6','S_6',term)
    #self.makeIndices()
class S7(Field):
  _basename='S7'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'7','lorL':'1','lorR':'1'},{'Y':frac(0,1)},'S7','S_7',term)
    #self.makeIndices()

class V1(Field):
  _basename='V1'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'1','lorL':'2','lorR':'2'},{'Y':frac(0,1)},'V1','V_1',term)
    #self.makeIndices()
class V2(Field):
  _basename='V2'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'2','lorL':'2','lorR':'2'},{'Y':frac(0,1)},'V2','V_2',term)
    #self.makeIndices()
class V3(Field):
  _basename='V3'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'3','lorL':'2','lorR':'2'},{'Y':frac(0,1)},'V3','V_3',term)
    #self.makeIndices()
class V4(Field):
  _basename='V4'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'4','lorL':'2','lorR':'2'},{'Y':frac(0,1)},'V4','V_4',term)
    #self.makeIndices()
class V5(Field):
  _basename='V5'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'5','lorL':'2','lorR':'2'},{'Y':frac(0,1)},'V5','V_5',term)
    #self.makeIndices()
class V6(Field):
  _basename='V6'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'6','lorL':'2','lorR':'2'},{'Y':frac(0,1)},'V6','V_6',term)
    #self.makeIndices()
class V7(Field):
  _basename='V7'
  def __init__(self,term=None):
    Field.__init__(self,frac(1,1),{'su3':'1','su2':'7','lorL':'2','lorR':'2'},{'Y':frac(0,1)},'V7','V_7',term)
    #self.makeIndices()


