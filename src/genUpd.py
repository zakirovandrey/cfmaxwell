#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import time

dim=3
anisotropy=False 

def makeCombinations(lst, res=['']):
  if len(lst)<=0: return res
  else: return makeCombinations(lst[:-1], [c+n for n in res for c in lst[-1]])

xyz = (1,2,4);
rules    ={"D":"d+","S":"s+","L":"sd","R":"ds","I":"_s","Y":"s_","X":"dx","P":"dp"} 
rules_all={"D":"d+","S":"s+","L":"sd","R":"ds","I":"_s","Y":"s_","X":"dx", "Z":"xx"} 
dek = lambda L1, L2 : [(L1_item + L2_item) for L2_item in L2 for L1_item in L1]
retdatas = lambda X : reduce(dek,map(lambda x: rules[x],X))
def cyclic(st, shift=0): return [ st[-3+shift],st[-2+shift],st[-1+shift] ]
Pname = lambda cell_mini, Cell_big: reduce( lambda p,q:p+q, map(lambda c,C: (c,C.lower())[c=='+'],cell_mini,Cell_big) )  # cl.replace("+","s" или "d")
#getkpml = lambda name, ax: 'Kpml%s['%('XYZ'[ax]) + (name.split('(')[-1][:-1] if ax==2 else name.split('.')[0]+'.pmli%s'%('xyz'[ax])) + ']'  # for vectorization
getkpml = lambda dat, ax: 'Kpml%s['%('XYZ'[ax]) + (dat.name.split('.')[0]+'.pmli%s*2'%('xyz'[ax])) + '+%d'%dat.coords[ax] + ']'
getzind = lambda name: name.split('(')[-1][:-1]

def LRv(str): 
  'makes left value with [] instead of ()'
  lvalue,rvalue = str.split('=')
  return lvalue.replace('(','[').replace(')',']')+'='+rvalue

pmlTmpl  = lambda d: LRv(r"  %(pml)s = %(pml)s*%(kpml)s.k1_%(gamma)s + %(diff)s*%(coff)s*%(kpml)s.k2_%(gamma)s  ;"%d)
difTmpl1 = lambda d: LRv(r"  %(pml)s+=                             %(diff)s*%(coff)s                                        ;"%d)
difTmpl2 = lambda d: LRv(r"  %(pml)s+=                             %(dif1)s*%(cof1)s + %(dif2)s*%(cof2)s                    ;"%d)
difTmpl3 = lambda d: LRv(r"  %(pml)s+=                             %(dif1)s*%(cof1)s + %(dif2)s*%(cof2)s + %(dif3)s*%(cof3)s;"%d)
averTmpl = lambda d: LRv(r"  %(Ex)s = %(Dx)s*%(depsXX)s + 0.5*(%(sumDym)s*%(depsXYm)s + %(sumDyp)s*%(depsXYp)s + \
                                                               %(sumDzm)s*%(depsXZm)s + %(sumDzp)s*%(depsXZp)s);"%d)
coffArrTmpl = r"CoffArr[%(fname)s.I[%(iz)s]]" if dim==2 else r"CoffArr[%(fname)s.I]"
coffTmpl = r"k%(type)s[%(i0)s][%(i1)s]"
coffDict = lambda d,ax='',i='': {'fname':d.indname,'iz':getzind(d.name),'i0':ax,'i1':i,'type':d.type[0]}
matTmpl  = r"%(p0)s + %(p1)s + %(p2)s + %(p3)s + %(p4)s"
diffsDeclared=[]; CoffArrDeclared=[]; aversDeclared=[]
def declareDiff(diffs, dfnot=None):
  dfRet=[]
  for df in diffs: 
    newdiff = "diff%d"%(len(diffsDeclared))
    if df not in diffsDeclared: 
      if df!=dfnot: print "  MainType %s=%s;"%(newdiff,df); 
      diffsDeclared.append(df); dfRet.append(newdiff) 
    else: dfRet.append("diff%d"%diffsDeclared.index(df))
  return dfRet
def declareCoffArr(coff):
  coffRet=[]
  newcoff = "CfArr%d"%(len(CoffArrDeclared))
  if coff not in CoffArrDeclared: print "  CoffStruct& %s=%s;"%(newcoff,coff); CoffArrDeclared.append(coff); coffRet = newcoff
  else: coffRet = "CfArr%d"%CoffArrDeclared.index(coff)
  return coffRet
def declare(vals, Declared=None, pref='', text='', exclude=None):
  valRet=[]
  for val in vals: 
    newVal = pref+"%d"%(len(Declared))
    if val not in Declared: 
      if val!=exclude: print text+" %s=%s;"%(newVal,val); 
      Declared.append(val); valRet.append(newVal) 
    else: valRet.append(pref+"%d"%Declared.index(val))
  return valRet


#########################################Main Class###################################################
class data:                            # данные в точке
  def __init__(self, type='', n='', r=True, C=None): self.real=r; self.name=n; self.type=type; self.Cell=C
  name=''; fname=''; real=True                   # имя значения; может быть реальное и виртуальное(заграничное)
  type='Wi';                           # Ei или Ti
  near = [ [None,None], [None,None], [None,None] ]          # соседи (x.p, x.m)              __
  coords=[0,0,0]                       # координаты внутри ConeFolda  ( 0-0.5-1-1.5 // XOXO /_/ )
  def main(self, ax):
    finiteDiff = cyclic( map( lambda i: "(%s-(%s))"%(i[0].name,i[1].name), self.near ), ax )
    diff = declare(finiteDiff, exclude=finiteDiff[0], Declared=diffsDeclared, pref='diff', text="  MainType")
    coffArr = declare([coffArrTmpl%coffDict(self)], Declared=CoffArrDeclared, pref='CfArr', text="  CoffStruct&")[0]
    coff = map( lambda i: coffArr+'.'+coffTmpl%coffDict(self,ax,i), (0,1))
#    if self.type=='Di': print "  Diold = %s;"%self.name
    print difTmpl2({'pml':self.name, 'dif1':diff[1],'cof1':coff[0], 'dif2':diff[2],'cof2':coff[1]})
#    if self.type=='Di': print "  %s = %s-Diold;"%(self.name.replace('Di','rot'),self.name)
  def pml1d(self, ax, pmlax):
  # в 1d pmlе : Ei_pml0 --- вдоль pml , Ei_pml1 --- перпендикулярно pml  // H--- аналогично
    finiteDiff = map( lambda i: "(%s-(%s))"%(i[0].name,i[1].name), self.near )
    diff = declare(finiteDiff, exclude=finiteDiff[ax], Declared=diffsDeclared, pref='diff', text="  MainType")
    coffArr = declare([coffArrTmpl%coffDict(self)], Declared=CoffArrDeclared, pref='CfArr', text="  CoffStruct&")[0]
    coff = map( lambda i: coffArr+'.'+coffTmpl%coffDict(self,ax,i), (0,1))
#    if self.type=='Di': print "  Diold = %s;"%self.name
    if ax==pmlax:
      print difTmpl2({'pml':self.name, 'dif1':diff[pmlax-2],'cof1':coff[0], 'dif2':diff[pmlax-1],'cof2':coff[1]})
    else        : 
      print pmlTmpl ({'pml':self.pmlnames[0], 'diff':diff[pmlax     ], 'coff':coff[(1,0)[(ax+1)%3==pmlax]], 'kpml': getkpml(self,pmlax), 'gamma':'H'})
      print difTmpl1({'pml':self.pmlnames[1], 'diff':diff[3-ax-pmlax], 'coff':coff[(0,1)[(ax+1)%3==pmlax]]})
      print LRv("  %s = %s + %s;"%(self.name, self.pmlnames[0], self.pmlnames[1]))
#    if self.type=='Di': print "  %s = %s-Diold;"%(self.name.replace('Di','rot'),self.name)
  def pml2d(self, ax, pmlax1, pmlax2):
  # в 2d pmlе : Ei_pml0 --- вдоль pml0, Ei_pml1 --- вдоль pml1 (или может получится наоборот если pmlx+pmlz)
  #                                либо
  #             Ei_pml0 --- вдоль pml,  Ei_pml1 --- вдоль не-pml
    finiteDiff = map( lambda i: "(%s-(%s))"%(i[0].name,i[1].name), self.near )
    diff = declare(finiteDiff, exclude=finiteDiff[ax], Declared=diffsDeclared, pref='diff', text="  MainType")
    coffArr = declare([coffArrTmpl%coffDict(self)], Declared=CoffArrDeclared, pref='CfArr', text="  CoffStruct&")[0]
    coff = map( lambda i: coffArr+'.'+coffTmpl%coffDict(self,ax,i), (0,1))
#    if self.type=='Di': print "  Diold = %s;"%self.name
    if ax!=pmlax1 and ax!=pmlax2: 
      print pmlTmpl ({'pml':self.pmlnames[0], 'diff':diff[pmlax1], 'coff':coff[(1,0)[(ax+1)%3==pmlax1]], 'kpml':getkpml(self,pmlax1), 'gamma':'H'})
      print pmlTmpl ({'pml':self.pmlnames[1], 'diff':diff[pmlax2], 'coff':coff[(0,1)[(ax+1)%3==pmlax1]], 'kpml':getkpml(self,pmlax2), 'gamma':'H'})
      print LRv("  %s = %s + %s;"%(self.name, self.pmlnames[0], self.pmlnames[1]))
    else:
      print pmlTmpl ({'pml':self.pmlnames[0], 'diff':diff[{pmlax1:pmlax2,pmlax2:pmlax1}[ax]], 'coff':coff[(0,1)[(ax+1)%3==3-pmlax1-pmlax2]], 'kpml':getkpml(self,{pmlax1:pmlax2,pmlax2:pmlax1}[ax]), 'gamma':'H'})
      print difTmpl1({'pml':self.pmlnames[1], 'diff':diff[3-pmlax1-pmlax2                  ], 'coff':coff[(1,0)[(ax+1)%3==3-pmlax1-pmlax2]], 'kpml':getkpml(self,3-pmlax1-pmlax2)})
      print LRv("  %s = %s + %s;"%(self.name, self.pmlnames[0], self.pmlnames[1]))
#    if self.type=='Di': print "  %s = %s-Diold;"%(self.name.replace('Di','rot'),self.name)
  def pml3d(self, ax):
    finiteDiff = cyclic( map( lambda i: "(%s-(%s))"%(i[0].name,i[1].name), self.near ), ax )
    diff = declare(finiteDiff, exclude=finiteDiff[0], Declared=diffsDeclared, pref='diff', text="  MainType")
    coffArr = declare([coffArrTmpl%coffDict(self)], Declared=CoffArrDeclared, pref='CfArr', text="  CoffStruct&")[0]
    coff =         map( lambda i: coffArr+'.'+coffTmpl%coffDict(self,ax,i), (0,1)  )
#    if self.type=='Di': print "  Diold = %s;"%self.name
    print pmlTmpl({'pml':self.pmlnames[0], 'diff':diff[1], 'coff':coff[0], 'kpml':getkpml(self,(ax+1)%3), 'gamma':'H'})
    print pmlTmpl({'pml':self.pmlnames[1], 'diff':diff[2], 'coff':coff[1], 'kpml':getkpml(self,(ax+2)%3), 'gamma':'H'})
    print LRv("  %s = %s + %s;"%(self.name, self.pmlnames[0], self.pmlnames[1]))
#    if self.type=='Di': print "  %s = %s-Diold;"%(self.name.replace('Di','rot'),self.name)
  def make_avers(self, ax, coff_prepare, val):
    part0 = "%s*%s"%(coff_prepare[0],coords2name(self.coords, self.Cell, {'Ei':val, 'Hi':'Bi'}[self.type]))
    def nears(dir):
      _coord = map( lambda i: self.coords[i]+{ax:dir[0],(ax+1)%3:dir[1],(ax+2)%3:dir[2]}[i], (0,1,2) )
#      coord = map( lambda i: (_coord[i],1.0-_coord[i])[self.Cell[i]=='Y' and _coord[i]>0.5], (0,1,2) )
#      sign=0
#      for i in 0,1,2: sign = (sign + ( self.Cell[i]=='Y' and _coord[i]>0.5 and (i+1)%2 ))%2
      return (coords2name( _coord, self.Cell, {"Ei":val, "Hi":"Bi"}[self.type] ),"0.0")[True in map(lambda i:self.Cell[i]=='Y' and _coord[i]>0.5, (0,1,2))]
    Avers = [part0,];
    Avers.append( "%s1*(%s)"%(coff_prepare[1],"+".join([nears([-0.5,-0.5, 0  ]),nears([+0.5,-0.5, 0  ])])) )
    Avers.append( "%s2*(%s)"%(coff_prepare[1],"+".join([nears([-0.5,+0.5, 0  ]),nears([+0.5,+0.5, 0  ])])) )
    Avers.append( "%s1*(%s)"%(coff_prepare[2],"+".join([nears([-0.5, 0  ,-0.5]),nears([+0.5, 0  ,-0.5])])) )
    Avers.append( "%s2*(%s)"%(coff_prepare[2],"+".join([nears([-0.5, 0  ,+0.5]),nears([+0.5, 0  ,+0.5])])) )
    anisoParts = declare(Avers, Declared=aversDeclared, pref='aver', text='  MainType')
    return anisoParts
  def material_subpxl(self,ax):
    em = {"Ei":'deps',"Hi":'dmu'}[self.type]   # eps or mu
    coff_prepare = map( lambda i: "Mat.%s%s"%(em,'XYZ'[ax]+cyclic('XYZ',i)[ax]) , range(3))
    anisoParts = self.make_avers(ax, coff_prepare, "Di")
    print "  %s = %s;"%(self.name, "+".join(anisoParts))
  def material_simple(self,ax):
    coff0="Mat.d%s"%{"Ei":'eps',"Hi":'mu'}[self.type]
    print "  %s = %s*%s;"%(self.name,coff0,coords2name(self.coords, self.Cell, {'Ei':'Di', 'Hi':'Bi'}[self.type]))
  def material_dispersion(self,ax):
    return
    cell = self.name.split('.')[0]
    rot = coords2name(self.coords, self.Cell, {'Ei':'Di', 'Hi':'Bi'}[self.type]).replace("Di","rot") 
    if self.type=='Ei': print "  countEandJ(%s.Ei[%d],%s.Eim[%d],%s.Ji[%d],%s.Jim[%d], %s, DispCfArr, %d);"%(cell,ax, cell,ax, cell,ax, cell,ax, rot, ax)
    else: self.material_simple(ax)
  def material_subpixl_disp(self,ax):
    return
    if self.type=="Hi": self.material_subpxl(ax); return
    cell = self.indnameL.split('.')[0]
    coff_prepare = map( lambda i: "Mat.Proj.%s"%('XYZ'[ax]+cyclic('XYZ',i)[ax]) , range(3) )
    anisoParts = self.make_avers(ax, coff_prepare, "rot")
#    is_p = ("","_p")['Y' in self.Cell]
    is_p = ""

    #------------------------------------1-----------------------------
    print "  auxfld = &%s.auxfld%s[0]; "%(cell,is_p)
    print "  auxfld->rot[%d] = Mat.f1[%d]*(%s);"%(ax,ax,"+".join(anisoParts))
    print "  countEandJ(auxfld->Ei[%d], auxfld->Eim[%d], auxfld->Ji[%d], auxfld->Jim[%d], auxfld->rot[%d], DispArr[auxfld->I[%d]], %d);"%(ax, ax, ax, ax, ax, ax, ax)

    #------------------------------------2-----------------------------
    print "  auxfld = &%s.auxfld%s[1]; "%(cell,is_p)
    print "  auxfld->rot[%d] = Mat.f2[%d]*(%s);"%(ax,ax,"+".join(anisoParts))
    print "  countEandJ(auxfld->Ei[%d], auxfld->Eim[%d], auxfld->Ji[%d], auxfld->Jim[%d], auxfld->rot[%d], DispArr[auxfld->I[%d]], %d);"%(ax, ax, ax, ax, ax, ax, ax)
 
    #------------------------------------3-----------------------------
    print "  auxfld = &%s.auxfld%s[2]; "%(cell,is_p)
    print "  auxfld->rot[%d] = %s-(%s);"%(ax,coords2name(self.coords, self.Cell, "rot"),"+".join(anisoParts))
    print "  countEandJ(auxfld->Ei[%d], auxfld->Eim[%d], auxfld->Ji[%d], auxfld->Jim[%d], auxfld->rot[%d], DispArr[auxfld->I[%d]], %d, Md);"%(ax, ax, ax, ax, ax, ax, ax)

    #########################################################################
    print "  %s.Ei[%d] = %s.auxfld%s[0].Ei[%d] + %s.auxfld%s[1].Ei[%d] + %s.auxfld%s[2].Ei[%d];"%(self.name.split('.')[0],ax,cell,is_p,ax,cell,is_p,ax,cell,is_p,ax)
 
  def flux(self, ax): pass
#    if self.type=='Si':  
#      print "  %s+= 0.5*(%s+%s) - 0.5*(%s+%s);"%(self.name, 
#              self.near[(ax+2)%3][0].name.replace('EHi','EHi_p'),self.near[(ax+2)%3][1].name.replace('EHi','EHi_p'),
#              self.near[(ax+1)%3][0].name.replace('EHi','EHi_m'),self.near[(ax+1)%3][1].name.replace('EHi','EHi_m')  )
#    if self.type=='EHi': 
#      print "  MainType Eiav_%d = 0.5*(Eiold_%d+%s);"%(ax, ax, self.name.replace('EHi','Ei'))
#      print "  %s*= Eiav_%d;"%(self.name.replace('EHi','EHi_m'), ax)
#      print "  %s*= Eiav_%d;"%(self.name.replace('EHi','EHi_p'), ax)
#    if self.type=='avHi':
#      print "  %s = (%s+%s);    "%(self.name.replace('EHi','EHi_m'), self.near[(ax+1)%3][0].name, self.near[(ax+1)%3][1].name )
#      print "  %s = (%s+%s);    "%(self.name.replace('EHi','EHi_p'), self.near[(ax+2)%3][0].name, self.near[(ax+2)%3][1].name )
#      print "  %s = 0.5*(%s+%s);"%(self.name.replace('EHi','avH%s'%('xyz'[(ax+2)%3])), self.near[(ax+1)%3][0].name, self.near[(ax+1)%3][1].name )
#      print "  %s = 0.5*(%s+%s);"%(self.name.replace('EHi','avH%s'%('xyz'[(ax+1)%3])), self.near[(ax+2)%3][0].name, self.near[(ax+2)%3][1].name )

def getcellname(cl, Cbig):
  "возвращает имя ячейки"
  plus_pos=[]
  for pp in range(len(cl)):
    if cl[pp]=="+": plus_pos.append("i"+"xyz"[pp]+Pname(cl,Cbig).replace('p','d'))        # находим позиции где есть "+"
  array_num=("0","+".join(plus_pos))[len(plus_pos)>0]
  return "F%s[%s]"%(Pname(cl,Cbig),array_num) 
def coords2cell(coor, Cl):   
  cl = retdatas(Cl)
  idd=0; 
  for i in (0,1,2): idd=idd+(xyz[i],0)[coor[i]<1]
  return cl[idd]
def coords2name(coor, Cl, type=''):           # Используется только для задания Filds names, Ind names и Neighbours names
  "по координате возвращает имя"
  sum_coors=coor[0]+coor[1]+coor[2]
  cell,k='',''
  for c in coor: k=k+'01'[int(c)==c]
  nn = ('[o]',"[%d]"%k.find('1'),"[%d]"%k.find('0'),'[x]')[k.count('1')]
  if type!='Bi' and type!='Ei' and type!='Si': cell=coords2cell(coor,Cl)[:dim]
  else: cell=coords2cell(map(lambda c: (c,c-1)[int(c)==c],coor),Cl)[:dim]
  return "%s.%s"%(getcellname(cell, Cl[:dim]), type+nn)       # XXX убрать [:-1] при dim=3
def doUpdate(dat, axis=0, Cell=None, materialEqSubpxl=False, dispersion=False):
  cell = coords2cell(dat.coords,Cell)
  pml = [0,0,0]       # координаты по которым есть pml
  for i in (0,1,2): pml[i] = 1 if (cell[i]=="s" and dat.type!='Bi' and dat.type!='Si')  or (cell[i]=="+" and Cell[i]=="S") else 0
  if dat.real:
    if dat.type=='Di' or dat.type=='Bi':
      if pml.count(1)==0: dat.main (axis                                           )
      if pml.count(1)==1: dat.pml1d(axis, pml.index(1)                             )
      if pml.count(1)==2: dat.pml2d(axis, pml.index(1), pml.index(1,pml.index(1)+1))
      if pml.count(1)==3: dat.pml3d(axis                                           )
    elif dat.type=='Si' or dat.type=='EHi' or dat.type=='avHi':
      dat.flux(axis)
    else:
      if not ((dat.type=="Ei" or (dat.type=="Hi" and Cell.count('I')>1) or (dat.type=="Hi" and Cell[axis]!='I')) and 'I' in Cell) :
        if dat.type=='Ei': print "  MainType Eiold_%d = %s;"%(axis,dat.name)
        if materialEqSubpxl                       : dat.material_subpxl(axis)
        if dispersion and not materialEqSubpxl    : dat.material_dispersion(axis)
        if dispersion and materialEqSubpxl        : dat.material_subpixl_disp(axis)
        if not materialEqSubpxl and not dispersion: dat.material_simple(axis)
  if dat.real and dat.type=='Hi' and (Cell[(axis+1)%3]=='I' or Cell[(axis+2)%3]=='I'):
    Bcord = [dat.coords[0],dat.coords[1],dat.coords[2]]; Bcord[Cell.find('I')]+=1.0
    Bname=coords2name(Bcord,Cell,'Bi');
    if '_' not in Bname: print "  %s = %s%s;"%(dat.name,'+-'[Cell.find('I')%2],Bname)
  ##-----------Добавление источника---------------------------------------#
  field = dat.name.split('.')[0]
#  func  = dat.name.split('.')[1][:1]+'xyz'[int(dat.name.split('.')[1][3:4])]
  func  = dat.name.split('.')[1][:1]+'xyz'[int(dat.name[-2])]
  sourceAxeN = cell.find('x')
  if dim==2: x1, x2, it = "%s.i%s"%(field,"YXX"[sourceAxeN]), getzind(dat.name), "%s.it"%field
  if dim==3: x1, x2, it = "%s.i%s"%(field,"YXX"[sourceAxeN]), "%s.i%s"%(field,"ZZY"[sourceAxeN]), "%s.it"%field
#  if cell.count('x')==1 and dat.real and (dat.type=="Hi"):
#    print "  if(%s.isSource) %s += Src%s(%s,%s,%s);"%(field, dat.name.replace('(','[').replace(')',']'), func, x1, x2, it)
#    if "  %s++;"%it not in post_grow: post_grow.append("  %s++;"%it)
  if Cell[-1]=='R' and dat.real and dat.type=='Di' and abs(dat.coords[2])-0.5<0.1:
    srcfield = getcellname(coords2cell([dat.coords[0],dat.coords[1],dat.coords[2]+1],Cell),Cell)
    x1, x2, it = "%s.SrcInd.iX"%srcfield, "%s.SrcInd.iY"%srcfield, "%s.SrcInd.it"%srcfield
    print "  if(%s.SrcInd.isSource) %s += Src%s(%s,%s,%s,%s);"%(srcfield, dat.name, func, x1,x2,"0.5*pars.dz",it)
    if "  %s++;"%it not in post_grow: post_grow.append("  %s++;"%it)
#  Fddd[0].Ex+=SrcEx(Fdds[ixdds].SrcInd)
#  Fddd[0].Ey+=SrcEy(Fdds[iydds].SrcInd)
  if dat.real and dat.type=='Di' and abs(dat.coords[2])-0.5<0.1:
    srcfield = getcellname(coords2cell([dat.coords[0],dat.coords[1],dat.coords[2]],Cell),Cell)
#    x1, x2, x3, it = "%s.ind.x"%srcfield, "%s.ind.y"%srcfield, "%s.ind.z"%srcfield, "%s.ind.time"%srcfield
#    print "  if(isSourceInside(%s.ind)) SrcInside%s(%s,%s,%s,%s,%s);"%(srcfield, func, dat.name, x1,x2,x3,it)

def vecSum(coor1,coor2): return [coor1[i]+coor2[i] for i in range(len(coor1))]    # векторная сумма координат

def setneighbours(Dat, Cell, zindx):
  Dat.near = [[data(),data()],[data(),data()],[data(),data()]]
  for id in 0,1,2: 
    for mp in 0,1: Dat.near[id][mp].coords = [ Dat.coords[ic] for ic in 0,1,2 ]
  for dx in  (0,1,2):
    Dat.near[dx][0].coords[dx] = Dat.coords[dx]+0.5                  # координаты соседей
    Dat.near[dx][1].coords[dx] = Dat.coords[dx]-0.5                  #  
  #-----------------neighbours-------------------------------------------#
  revtypes={'Ei':'Di','Di':'Hi','Hi':'Bi','Bi':'Ei','Si':'EHi','EHi':'Hi','avHi':'Hi','':''}
  #задать соседей
  for dx in (0,1,2):
#    inadx = i if i==dx else 3-i-dx        ## если i==dx равно i, если i!=dx --- равно не i и не dx
    for di in 0,1:                    # 0 -- plus, 1 -- minus
      Dat.near[dx][di].name = "%s%s"%(coords2name(Dat.near[dx][di].coords,Cell,revtypes[Dat.type]), zindx[int(Dat.near[dx][di].coords[2])])  # xp, yp, zp
    #задать значения для соседей, если они виртуальные
    for di in 0,1:
      koord = Dat.near[dx][di].coords; kcel = coords2cell(koord,Cell) 
      if kcel.count("_")>0 or map(lambda ci: koord[ci]==1.5 and kcel[ci]=='x', (0,1,2)).count(True)>0 or map(lambda ci: koord[ci]>0.6 and Cell[ci]=='Y', (0,1,2)).count(True)>0 :
#        if kcel[2]!="_": Dat.near[dx][di].real=False; Dat.near[dx][di].name = "-"+Dat.near[dx][1-di].name   ####!! УБрать if kcel[2]!="_" при отсутствии векторизации!!!
        Dat.near[dx][di].real=False; Dat.near[dx][di].name = "+-"[dx%2]+Dat.near[dx][1-di].name   ####!! УБрать if kcel[2]!="_" при отсутствии векторизации!!!

def makeUpd(Cell, zindxDicti=None, d=dim, ConeType="ConeFold", h=0.5):     # ConeType пока только одномерный по Z # h --- высота ConeFolda
  cell = retdatas(Cell)      #  --- фактически список ячеек в правильном порядке
  D =  (data(type="Di"  ,C=Cell), data(type="Di"  ,C=Cell), data(type="Di"  ,C=Cell));
  E =  (data(type="Ei"  ,C=Cell), data(type="Ei"  ,C=Cell), data(type="Ei"  ,C=Cell));
  B =  (data(type="Bi"  ,C=Cell), data(type="Bi"  ,C=Cell), data(type="Bi"  ,C=Cell));
  H =  (data(type="Hi"  ,C=Cell), data(type="Hi"  ,C=Cell), data(type="Hi"  ,C=Cell));
  S =  (data(type="Si"  ,C=Cell), data(type="Si"  ,C=Cell), data(type="Si"  ,C=Cell));
  EH = (data(type="EHi" ,C=Cell), data(type="EHi" ,C=Cell), data(type="EHi" ,C=Cell));
  avH =(data(type="avHi",C=Cell), data(type="avHi",C=Cell), data(type="avHi",C=Cell));
  for i in (0,1,2):
    D[i].coords =  cyclic([0,.5,.5],3-i); B[i].coords = cyclic([ .5,1,1],3-i)  # B[i] хранится в левой ячейке
    E[i].coords =  cyclic([1,.5,.5],3-i); H[i].coords = cyclic([1.5,1,1],3-i)  # E[i] хранится в левой ячейке
    S[i].coords =  cyclic([ .5,1,1],3-i)  # S[i] хранится в левой ячейке просто по аналогии с B и E
    EH[i].coords = cyclic([1,.5,.5],3-i)  # 
    avH[i].coords =cyclic([0,.5,.5],3-i)  #
    for Fd in E[i],D[i],H[i],B[i],S[i],EH[i],avH[i]:
      kcel = coords2cell(Fd.coords, Cell)
      if kcel.count("_")>0 or map(lambda ci: Fd.coords[ci]==1.5 and kcel[ci]=='x', (0,1,2)).count(True)>0: Fd.real=False
      if (Fd.type=='Bi' or Fd.type=='Ei' or Fd.type=='EHi' or Fd.type=='Si') and 'I' in Cell: Fd.real=False
      if map(lambda ci: Fd.coords[ci]>0.6 and cell[ci]=='Y', (0,1,2)).count(True)>0: Fd.real=False
  zindx=('','','')
  if d==2: zindx = zindxDict[Cell[-1]]   # (iz, iz+1)
  #тут надо понимать, что если Cell[-1] не pml по Z, то pmlzindx может быть нужен для pml по X или Y
  #---------------------------Fields names-------------------------------#
  D[0].name,  D[1].name, D[2].name = map( lambda x: "%s%s"%(coords2name(D[x].coords,Cell,D[x].type),zindx[(0,0,0)[x]]), (0,1,2) )
  E[0].name,  E[1].name, E[2].name = map( lambda x: "%s%s"%(coords2name(E[x].coords,Cell,E[x].type),zindx[(0,0,1)[x]]), (0,1,2) )
  B[0].name,  B[1].name, B[2].name = map( lambda x: "%s%s"%(coords2name(B[x].coords,Cell,B[x].type),zindx[(1,1,0)[x]]), (0,1,2) )
  H[0].name,  H[1].name, H[2].name = map( lambda x: "%s%s"%(coords2name(H[x].coords,Cell,H[x].type),zindx[(1,1,1)[x]]), (0,1,2) )
  S[0].name,  S[1].name, S[2].name = map( lambda x: "%s%s"%(coords2name(S[x].coords,Cell,S[x].type),zindx[(1,1,0)[x]]), (0,1,2) )
  EH[0].name ,EH[1].name,EH[2].name= map( lambda x: "%s%s"%(coords2name(EH[x].coords,Cell,EH[x].type),zindx[(1,1,0)[x]]), (0,1,2) )
  avH[0].name,avH[1].name,avH[2].name= map( lambda x: "%s%s"%(coords2name(avH[x].coords,Cell,'EHi'),zindx[(1,1,0)[x]]), (0,1,2) )
  #-------Ind names------------------------------------------------------#
  for W in D+E+B+H+S: W.indname = coords2name(map(lambda c: (0,1)[c=='I'],Cell),Cell).split('.')[0] + ".ind"
  celln=7
  for i in 0,1,2: 
    if cell[7][i]=="_": celln-= {0:1,1:2,2:4}[i]
  for W in D+E+B+H+S: W.indname = getcellname(cell[celln],Cell) + ".ind"
#  for W in D+E+B+H+S: W.indname = getcellname(cell[celln],Cell) + ".ind" if '_' not in cell[7] else W.indname
  #Left index#
  for W in D+E+B+H+S: W.indnameL = coords2name(map(lambda c: (0,1)[c=='I'],Cell),Cell).split('.')[0] + ".ind"
  for W in D+E+B+H+S: W.indnameL = getcellname(cell[0],Cell) + ".ind" if '_' not in cell[0] else W.indname
  #-------PML names------------------------------------------------------#
  # в 1d pmlе : Si_pml0 --- вдоль pml , Si_pml1 --- перпендикулярно pml                  // V и T --- аналогично
  # в 2d pmlе : Si_pml0 --- вдоль pml0, Si_pml1 --- вдоль pml1, Si_pml2 --- вдоль не pml // 
  D[0].pmlnames, D[1].pmlnames, D[2].pmlnames = map( lambda x: map( lambda pmlcoord: "%s_pml%s_%s%s"%(coords2name(D[x].coords,Cell,D[x].type)[:-3],pmlcoord,x,zindx[(0,0,0)[x]]), (0,1,2) ), (0,1,2) )
  E[0].pmlnames, E[1].pmlnames, E[2].pmlnames = map( lambda x: map( lambda pmlcoord: "%s_pml%s_%s%s"%(coords2name(E[x].coords,Cell,E[x].type)[:-3],pmlcoord,x,zindx[(0,0,1)[x]]), (0,1,2) ), (0,1,2) )
  B[0].pmlnames, B[1].pmlnames, B[2].pmlnames = map( lambda x: map( lambda pmlcoord: "%s_pml%s_%s%s"%(coords2name(B[x].coords,Cell,B[x].type)[:-3],pmlcoord,x,zindx[(1,1,0)[x]]), (0,1,2) ), (0,1,2) )
  H[0].pmlnames, H[1].pmlnames, H[2].pmlnames = map( lambda x: map( lambda pmlcoord: "%s_pml%s_%s%s"%(coords2name(H[x].coords,Cell,H[x].type)[:-3],pmlcoord,x,zindx[(1,1,1)[x]]), (0,1,2) ), (0,1,2) )
  #----------------------------------------------------------------------#
  for Fd in E+D+H+B+S+EH+avH: setneighbours(Fd,Cell,zindx)
  for Fd in E+D+H+B+S+EH+avH:
      for Fnear in Fd.near: setneighbours(Fnear[0],Cell,zindx); setneighbours(Fnear[1],Cell,zindx)
#  for Fd in E+D+B+H: 
#    for Fnear in Fd.near:
#      for fn in Fnear:
#        for Fnearnear in fn.near: setneighbours(Fnearnear[0],Cell,zindx); setneighbours(Fnearnear[1],Cell,zindx)
  #----------------------------------------------------------------------#
  if Cell[-1]=='Z': 
    print "  MatStruct& Mat = MatArr[%s.I];"%(D[0].indname)
    for i in (0,1,2):  doUpdate(D[i], i, Cell, materialEqSubpxl=False)
    return
#  if Cell[-1]!='X': print "  MatStruct& Mat = MatArr[%s.I];"%(D[0].indname)
#  if Cell[-1]!='X': print "  DispStruct& DispCfArr = DispArr[%s.I];"%(D[0].indname)
  print "  MatStruct& Mat = MatArr[%s.I];"%(D[0].indname)
#  print "  AuxField* auxfld;"
  for Subpxl,dsp in (False,False), (True,False) :
    global diffsDeclared; diffsDeclared=[]  
    global CoffArrDeclared; CoffArrDeclared=[]
    global aversDeclared; aversDeclared=[]
    Celltest=Cell.replace('P','D').replace('S','D').replace('L','D').replace('R','D').replace('Y','D')
    if Celltest=='DDD':
      print "  if(%s.isBound%s) {"%(D[0].indname,('==0','!=0')[Subpxl])
    if 'Y' not in Cell or not Subpxl or True:
      if dsp: print "  DispStruct& DispCfArr = DispArr[%s.I];"%(D[0].indname)
      for i in (0,1,2): 
        if {"ConeFold":1, 'OF':1, 'OT':0, 'XF':0, 'XT':0}[ConeType]: doUpdate(D[i],  i, Cell, materialEqSubpxl=Subpxl, dispersion=dsp)   
      for i in (0,1,2): 
        if {"ConeFold":1, 'OF':1, 'OT':0, 'XF':0, 'XT':0}[ConeType]: doUpdate(avH[i],i, Cell, materialEqSubpxl=Subpxl, dispersion=dsp)   
      for i in (0,1,2): 
        if {"ConeFold":1, 'OF':0, 'OT':1, 'XF':0, 'XT':0}[ConeType]: doUpdate(E[i],  i, Cell, materialEqSubpxl=Subpxl, dispersion=dsp)   
      for i in (0,1,2): 
        if {"ConeFold":1, 'OF':0, 'OT':1, 'XF':0, 'XT':0}[ConeType]: doUpdate(EH[i], i, Cell, materialEqSubpxl=Subpxl, dispersion=dsp)   
      for i in (0,1,2): 
        if {"ConeFold":1, 'OF':0, 'OT':0, 'XF':1, 'XT':0}[ConeType]: doUpdate(S[i],  i, Cell, materialEqSubpxl=Subpxl, dispersion=dsp)
      for i in (0,1,2): 
        if {"ConeFold":1, 'OF':0, 'OT':0, 'XF':1, 'XT':0}[ConeType]: doUpdate(B[i],  i, Cell, materialEqSubpxl=Subpxl, dispersion=dsp)
      for i in (0,1,2): 
        if {"ConeFold":1, 'OF':0, 'OT':0, 'XF':0, 'XT':1}[ConeType]: doUpdate(H[i],  i, Cell, materialEqSubpxl=Subpxl, dispersion=dsp)
    if Celltest!='DDD': break
    print "  }"
  if Cell.replace('P','D')=='DDD':
    print "  SaveNearField(%s, %s, %s, %s, %s, %s, %s, %s, %s);"%tuple(map(lambda a: getcellname(cell[a],Cell),(0,1,2,3,4,5,6,7))+[E[0].indname])
  if 'I' not in Cell: print "  %s.time++;"%(E[0].indnameL)
  print "  string Sfield; int DropIndex;"
#  print "  int DropIndex;"
  for Fd in (E[0],E[1],E[2],H[0],H[1],H[2]):#+H:#+D+B#+S:
    if Fd.real:
      nearestLeftInd = coords2name(map(int,Fd.coords),Cell).split('.')[0] + ".ind"
      Sfield=Fd.type[0]+'xyz'[int(Fd.name[-2])]
      print '  DropIndex = %s.isDrop; if(DropIndex>=0 && find(DropP.Fields4Sensors.begin(),DropP.Fields4Sensors.end(),"%s")!=DropP.Fields4Sensors.end()) {'%(nearestLeftInd,Sfield)
      print '     FILE* tt; char fname[256]; sprintf(fname,"%%s/tt%s-x%%04dy%%04dz%%04d.dat",DropP.DropPath.c_str(),DropP.Sensors[DropIndex].x,DropP.Sensors[DropIndex].y,DropP.Sensors[DropIndex].z); tt=fopen(fname,"a"); fprintf(tt,"%%g\\n",%s); fclose(tt);'%(Sfield,Fd.name)
      print "  }"
  Sfield="Ef"
  if False in map(lambda f:f.real, E): return
  print '  DropIndex = %s.isDrop; if(DropIndex>=0 && find(DropP.Fields4Sensors.begin(),DropP.Fields4Sensors.end(),"%s")!=DropP.Fields4Sensors.end()) {'%(nearestLeftInd,Sfield)
  print '     FILE* tt; char fname[256]; sprintf(fname,"%%s/tt%s-x%%04dy%%04dz%%04d.dat",DropP.DropPath.c_str(),DropP.Sensors[DropIndex].x,DropP.Sensors[DropIndex].y,DropP.Sensors[DropIndex].z); tt=fopen(fname,"a");'%Sfield
  print '     fprintf(tt,"%%g\\n",%s*%s+%s*%s+%s*%s); fclose(tt);'%(E[0].name,E[0].name,E[1].name,E[1].name,E[2].name,E[2].name)
  print "  }"

incpath = sys.argv[0][:sys.argv[0].find(sys.argv[0].split('/')[-1])]
stdout, sys.stdout = sys.stdout,open(incpath+'Update.inc.hpp', 'w')
for stc in reduce(dek,[rules.keys()]*dim) :
  s=''; post_grow=[]
  print "__attribute__((always_inline)) inline void " + stc + "act%s ("%(('PML','')[stc=='DD']),
  for l in retdatas(stc): # retdatas(stc) --- фактически список ячеек в правильном порядке 
    if not l.find("_")+1: s+=["fields%s* const F%s, "%(l.upper().replace('P','D'),l),('',"const int i%s%s, "%('xyz'[l.find("+")], Pname(l,stc)))[l.find('p')==-1],'',''][l.count('+')] #строка агруметов функции
  s = s[:-2] + "){"
  print s
  diffsDeclared=[]  
  CoffArrDeclared=[]
  aversDeclared=[]
  if dim==2:
  # zindxDict1 --- for Conefold2 , zindxDict2 --- for ConeFold1
    zindxDict1 =   {"D":("(iz)"   , "(iz+1)"    , "IZ+2!", "IZ+3!"), "S":("(pmliz)"         , "(pmliz+1)"        , "IZ+2!", "IZ+3!"), 
                    "L":(".LC(-1)", "(2*PMLNzV)", "IZ+2!", "IZ+3!"), "R":("(2*PMLNzV+NzV-1)", ".RC(NzV)"         , "IZ+2!", "IZ+3!"),
                    "I":(".LS(-1)", "(0)"       , "IZ+2!", "IZ+3!"), "Y":("(2*PMLNzV-1)"    , ".RS(2*PMLNzV)"    , "IZ+2!", "IZ+3!")   } 
    zindxDict2 =   {"D":("(iz-1)" , "(iz)"      , "IZ+2!", "IZ+3!"), "S":("(pmliz-1)"       , "(pmliz)"          , "IZ+2!", "IZ+3!"), 
                    "L":(".LC(-2)", ".LC(-1)"   , "IZ+2!", "IZ+3!"), "R":("(2*PMLNzV+NzV-2)", "(2*PMLNzV+NzV-1)" , "IZ+2!", "IZ+3!"),
                    "I":(".LS(-2)", ".LS(-1)"   , "IZ+2!", "IZ+3!"), "Y":("(2*PMLNzV-2)"    , "(2*PMLNzV-1)"     , "IZ+2!", "IZ+3!")   } 

   #--------------------------OF rank=1-------------------#
    makeUpd(stc+"Y", zindxDict1, ConeType="OF");
    makeUpd(stc+"Y", zindxDict2, ConeType="OF");
    makeUpd(stc+"Y", zindxDict2, ConeType="OT");
    makeUpd(stc+"Y", zindxDict2, ConeType="XF");
    #------------------------------------------------------#
    print "  for(int pmliz=2*PMLNzV     -3; pmliz>=0       ; pmliz--){ "; makeUpd(stc+"S", zindxDict1); print "  }"   # PML
    #------------------------------------------------------#
    #--------------------------OF rank=1-------------------#
    makeUpd(stc+"R", zindxDict1, ConeType="OF");
    makeUpd(stc+"R", zindxDict2, ConeType="OF");
    makeUpd(stc+"R", zindxDict2, ConeType="OT");
    makeUpd(stc+"R", zindxDict2, ConeType="XF");
    #------------------------------------------------------#
    #--------------------------OT rank=1-------------------#
    makeUpd(stc+"I", zindxDict1, ConeType="OT");
    makeUpd(stc+"I", zindxDict1, ConeType="XF");
    makeUpd(stc+"I", zindxDict2, ConeType="XT");
    makeUpd(stc+"I", zindxDict1, ConeType="XT");
    #------------------------------------------------------#
    print "  for(int    iz=2*PMLNzV+NzV-3;    iz>=2*PMLNzV ;    iz--){ "; makeUpd(stc+"D", zindxDict1); print "  }"   # main region
    #--------------------------OT rank=1-------------------#
    makeUpd(stc+"L", zindxDict1, ConeType="OT"); 
    makeUpd(stc+"L", zindxDict1, ConeType="XF"); 
    makeUpd(stc+"L", zindxDict2, ConeType="XT"); 
    makeUpd(stc+"L", zindxDict1, ConeType="XT"); 
    #------------------------------------------------------#

    # iz   : | 0 ------pml1d----- PMLNz-1 | PMLNz ----------------main---------------- PMLNz+Nz-1 | PMLNz+Nz -----pml1d----- 2*PMLNz+Nz-1 |
  if dim==3:
    print "  MainType Diold;"
#    if stc[2]=='X': makeUpd(stc[:-1]+'Z')
    makeUpd(stc);
#  for pg in post_grow: print pg
  if post_grow: print post_grow[-1];    # печатаем только последний it++, так как в остальных еще не сработали все Src.
  print "}"

sys.stdout.close()
