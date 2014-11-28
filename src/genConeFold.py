#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from types import *
from operator import *

def makeCombinations(lst, res=['']):
  if len(lst)<=0: return res
  else: return makeCombinations(lst[:-1], [c+n for n in res for c in lst[-1]])

class CFact:
  def __init__(self, Gen, actN):
    self.Gen=Gen; self.actN=actN
    self.podDatas,self.podDatasShift = ['',],[0,]
    for s in xrange(Gen.dim):
      #datas = self.pars4actN[actN[s]]
      datas = Gen.Rules2act['pars'][actN[s]]
      for k in xrange(1<<s):
        self.podDatas.append(self.podDatas[k] + datas[1])
        self.podDatas[k] += datas[0]
        self.podDatasShift.append(self.podDatasShift[k] + (0,1<<s)[datas[1] is 'p'])
        pass
      pass
    pass
  def PodActList(self, tC='A'):
    tier0=makeCombinations([self.Gen.Rules2act['subacts'][a][:2] for a in self.actN])
    tier1=makeCombinations([self.Gen.Rules2act['subacts'][a][2:] for a in self.actN])
    #tier0=makeCombinations([self.decomp4actN[a][:2] for a in self.actN])
    #tier1=makeCombinations([self.decomp4actN[a][2:] for a in self.actN])
    tier2=reduce(add,zip(tier0,tier1))
    tier0=filter(lambda (a,n,t):'-' not in a, zip(tier0,xrange(len(tier0)),['B']*len(tier0))); tier0.reverse()
    tier1=filter(lambda (a,n,t):'-' not in a, zip(tier1,xrange(len(tier1)),['T']*len(tier1))); tier1.reverse()
    tier2=filter(lambda (a,n,t):'-' not in a, zip(tier2,xrange(len(tier2)),['B','T']*len(tier0))); tier2.reverse()
    if tC is 'B': return tier0
    if tC is 'T': return tier1
    if tC is 'BT': return tier2
    return tier0+tier1
  def PodActList_mp(self):
    tier=makeCombinations([self.Gen.Rules2act['subacts_mp'][a] for a in self.actN])
    tier=filter(lambda (a,n):'-' not in a, zip(tier,xrange(len(tier)))); tier.reverse()
    return tier
  def getIsh(self, datI, shI):
    return filter(len, [('',self.Gen.rulesShift[self.podDatas[datI]][s])[(shI&(1<<s))>0] for s in xrange(self.Gen.dim)])
  def getPodActPar(self, par, parI, cntI, tC):
    '''Возвращает параметр parI (0..2^d) для под-act-а с базовым datas-ом cntI (0..2^d), исходя из шаблона par'''
    if tC is 'B': datI,poddatI = cntI&parI,cntI^parI
    if tC is 'T': datI,poddatI = cntI|parI,((1<<self.Gen.dim)-1)&(~(cntI^parI))
    if tC is 'X':
      parIl,cntIl=map(lambda v:(0,1,-1)[v],self.Gen.num2list(parI,3)),self.Gen.num2list(cntI)
      datI = self.Gen.list2num([(p+n)>=1 for (p,n) in zip(parIl,cntIl)])
      poddatI = self.Gen.list2num([1&~((p&1)^n) for (p,n) in zip(parIl,cntIl)])
    #datI,poddatI --- номера датаса и его поддатаса для параметра parI
    poddatIx = 0
    for s in xrange(self.Gen.dim-1,-1,-1):
      i2s = 1<<s
      if self.Gen.Rules2dim[self.podDatas[datI][s]]: poddatIx = 2*poddatIx+((poddatI&i2s)!=0)
    #poddatIx --- смещение параметра в массиве поддатасов для датаса datI
    #print par, parI, cntI, tC,self.podDatas[datI],"->",datI,poddatI,poddatIx
    if par.find('dat')>=0:
      dat_shift = self.podDatasShift[datI]
      return (self.Gen.datTmpl%datI+'->',self.Gen.datTmpl%(datI-dat_shift)+'[%s].'%('+'.join(self.getIsh(datI-dat_shift,dat_shift))))[dat_shift>0]+'datas'+('','+%d'%poddatIx)[poddatIx>0]
    if par.find('const int _I')>=0:
      sh = par[par.find('const int _I')+len('const int _I'):]
      sI = sh.find('p')
      if sI >= 0:
        if poddatI&(1<<sI): return '%d'%(1<<self.Gen.get_dim(sh[:sI]))
        return 'I'+sh
      sI = sh.find('m')
      if sI >= 0:
        if poddatI&(1<<sI): return '-I'+sh.replace('m','p')
        return '-%d'%(1<<self.Gen.get_dim(sh[:sI]))
      return '?'+sh
    return '======== not implemented yet,',par, n, t

class CFact_mp(CFact):
  def __init__(self, Gen, actN):
    CFact.__init__(self, Gen, actN)
    self.podDatas_mp,self.podDatasShift_mp = ['',],[0,]
    for s in xrange(Gen.dim):
      datas = Gen.Rules2act_mp['pars'][actN[s]]
      for k in xrange(3**s): self.podDatas_mp.append(self.podDatas_mp[k] + datas[2])
      for k in xrange(3**s): self.podDatas_mp.append(self.podDatas_mp[k] + datas[0])
      for k in xrange(3**s): self.podDatas_mp[k] += datas[1]
      for k in xrange(3**s): self.podDatasShift_mp.append(self.podDatasShift_mp[k] + (0,3**s)[datas[2] is 'p'])
      for k in xrange(3**s): self.podDatasShift_mp.append(self.podDatasShift_mp[k] + 2*(0,3**s)[datas[0] is 'm'])
      pass
    pass
  def PodActList_mp(self):
    #tier=makeCombinations([self.decomp4actNmp[a] for a in self.actN])
    tier=makeCombinations([self.Gen.Rules2act_mp['subacts'][a] for a in self.actN])
    tier=filter(lambda (a,n):'-' not in a, zip(tier,xrange(len(tier)))); tier.reverse()
    return tier
  def PodActList_mpPIC(self):
    #tier=makeCombinations([self.decomp4actNmpPIC[a] for a in self.actN])
    tier=makeCombinations([self.Gen.Rules2act_mp['subactsPIC'][a] for a in self.actN])
    tier=filter(lambda (a,n):'-' not in a, zip(tier,xrange(len(tier)))); tier.reverse()
    return tier
  def getIsh_mp(self, datI, shI):
    '''вычисляет сдвиг datas-а'''
    if '_' in self.podDatas_mp[datI]: return ('<no-data>',)
    #print datI, self.podDatas_mp[datI], self.Gen.rulesShift[self.podDatas_mp[datI]]
    shIl = self.Gen.num2list(shI,3)
    Ish=[]
    for s,sh in zip(xrange(self.Gen.dim), self.Gen.num2list(shI,3)):
      rul = self.Gen.rulesShift[self.podDatas_mp[datI]][s]
      Ish.append(('',rul, rul.replace('p','m'))[shIl[s]])
    return filter(len, Ish)
  def getPodActPar(self, par, parI, cntIl, tC):
    '''Возвращает параметр parI (0..3^d) для под-act-а с базовым datas-ом cntI (0..2^d), исходя из шаблона par'''
    if tC is 'X':
      parIl=map(lambda v:(0,1,-1)[v],self.Gen.num2list(parI,3))
      datIl = [(-1,-1,0,0,1,1)[p+n+2] for (p,n) in zip(parIl,cntIl)] # можно (p+n)/2 # номер (смещение) datas-а относительно базового cntI
      poddatIl = [(0,1,0,1,0,1)[p+n+2] for (p,n) in zip(parIl,cntIl)] # можно (p+n)%2 # номер poddatas-а в datas-е
      datI = self.Gen.list2num([((p+n)/2)%3 for (p,n) in zip(parIl,cntIl)], 3)
      poddatI = self.Gen.list2num([((p+n)%2)&3 for (p,n) in zip(parIl,cntIl)])
    poddatIx = 0
    for s in xrange(self.Gen.dim-1,-1,-1):
      i2s = 1<<s
      if self.Gen.Rules2dim[self.podDatas_mp[datI][s]]: poddatIx = 2*poddatIx+((poddatI&i2s)!=0)
    if par.find('dat')>=0:
      dat_shift = self.podDatasShift_mp[datI]
      return (self.Gen.datTmpl%datI+'->',self.Gen.datTmpl%(datI-dat_shift)+'[%s].'%('+'.join(self.getIsh_mp(datI-dat_shift,dat_shift))))[dat_shift!=0]+'datas'+('','+%d'%poddatIx)[poddatIx>0]
    if par.find('const int _I')>=0:
      sh = par[par.find('const int _I')+len('const int _I'):]
      sI = sh.find('p')
      if sI >= 0:
        if poddatIl[sI]==1: return '%d'%(1<<self.Gen.get_dim(sh[:sI]))
        if datIl[sI]<=0: return '-I'+sh.replace('p','m')
        return 'I'+sh
      sI = sh.find('m')
      if sI >= 0:
        if poddatIl[sI]==0: return '-%d'%(1<<self.Gen.get_dim(sh[:sI]))
        if datIl[sI]>=0: return '-I'+sh.replace('m','p')
        return 'I'+sh
      return '?'+sh
    return '======== not implemented yet,',par, n, t

class CFpodact(CFact):
  def __init__(self, act, actN):
    CFact.__init__(self, act.Gen, actN)
    self.nadact = act
    pass
  def getParsList(self, n, tC):
    pars_list, full_pars_list = self.Gen.get_pars(self.podDatas)
    fakt_pars = [self.nadact.getPodActPar(par, full_pars_list.index(par), n,tC) for par in pars_list]
    return fakt_pars

class CFpodact_mp(CFact_mp):
  def __init__(self, act, actN):
    CFact_mp.__init__(self, act.Gen, actN)
    self.nadact = act
  def getParsList(self, n, tC):
    pars_list, full_pars_list = self.Gen.get_pars_mp(self.podDatas_mp)
    fakt_pars = [self.nadact.getPodActPar(par, full_pars_list.index(par), n, tC) for par in pars_list]
    return fakt_pars

class Generator:
  def __init__(self, dim, types):
    self.Rules2act = {
      'pars': {'D':'dp', 'S':'dp', 'I':'_d', 'J':'dp', 'X':'dx', 'Y':'d_', 'P':'dd', 'Q':'dd', 'L':'ld', 'M':'dp', 'R':'dr'},
      'subacts':{'D':'DDDD', 'S':'SSSS', 'I':'-IID', 'X':'DXX-', 'Y':'DYY-', 'P':'DPPD', 'L':'-LLD', 'R':'DRR-'},
      'subacts_mp': {'D':'DD', 'I':'-J', 'J':'JD', 'X':'X-', 'Y':'Y-', 'P':'PQ', 'L':'-M', 'R':'R-'},
      'LR':{'L':'IL-', 'R':'-RY'},
    }
    self.Rules2act_mp = {
      #'pars': {'D':'mdp', 'S':'msp', 'J':'_dp', 'X':'mdp', 'P':'mdd', 'Q':'ddp', 'M':'ldp', 'R':'mdr'},
      'pars': {'D':'mdp', 'S':'msp', 'J':'_dp', 'X':'mdx', 'P':'mdd', 'Q':'ddp', 'M':'ldp', 'R':'mdr'},
      'subacts': {'D':'DD', 'I':'JD', 'J':'JD', 'X':'DX', 'P':'DP', 'Q': 'QD', 'L':'MD', 'M':'MD', 'R':'DR'},
      'subactsPIC': {'D':'DDDD', 'J':'-JDD', 'X':'DDX-', 'P':'DDPQ', 'Q': 'PQDD', 'M':'-MDD', 'R':'DDR-'},
    }
    self.Rules2rank = {}
    self.Rules2dim = {}
    for s in 'dmps': self.Rules2dim[s] = 1
    for s in 'xlr': self.Rules2dim[s] = 0
    self.dim=dim
    self.par_name_start=-(dim+2)-(dim==1)
    self.datTmpl='datas_____'[:2+(dim+1)/2:]+'%0'+'%dd'%((dim+1)/2)
    self.types=types
    print '// acts: %d'%len(makeCombinations([types]*dim))
    pass
  def num2list(self, num, bas=2):
    numLst = []
    for s in xrange(self.dim):
      numLst.append(num%bas)
      num /= bas
    return numLst
  def list2num(self, numLst, bas=2):
    numLstR = numLst[:]
    numLstR.reverse()
    return reduce(lambda r,v: r*bas+v, numLstR, 0)
  def get_dim(self, pd):
    return reduce(lambda r,s: r+self.Rules2dim[s], pd, 0)
  def add2rules(self, pd):
    self.rules[pd] = 'cubeLR<%d,T%%(Npd)d,%s>'%(self.get_dim(pd),self.Rules2rank.get(pd,self.rank)) + '* const '+self.datTmpl.replace('%','%(Npd)')
    self.add2rulesShift(pd)
    pass
  def add2rulesShift(self, pd):
    self.rulesShift[pd] = []
    #self.rulesShiftM[pd] = []
    for s in xrange(len(pd)):
      sh=pd[:s]+'p'+pd[s+1:]
      self.rules[sh] = 'const int _I'+sh
      self.rulesShift[pd].append('_I'+sh)
      shM=pd[:s]+'m'+pd[s+1:]
      self.rules[shM] = 'const int _I'+shM
      #self.rulesShiftM[pd].append('_I'+shM)
      pass
    pass
  def get_pars(self, datas, shift=0):
    full_pars_list = [self.rules.get(pd,'')%{'Npd':i+shift} for (i,pd) in zip(xrange(1<<self.dim),datas)]
    pars_list = filter(len, full_pars_list)
    pars_list = map(lambda i: pars_list[i], filter(lambda i: pars_list.index(pars_list[i])==i, xrange(len(pars_list))))
    return pars_list, full_pars_list
  def get_pars_mp(self, datas):
    full_pars_list = [self.rules.get(pd,'')%{'Npd':i} for (i,pd) in zip(xrange(3**self.dim),datas)]
    pars_list = filter(len, full_pars_list)
    pars_list = map(lambda i: pars_list[i], filter(lambda i: pars_list.index(pars_list[i])==i, xrange(len(pars_list))))
    return pars_list, full_pars_list
  def getTmplPars(self, formal_pars):
    template_pars = ','.join(map(lambda s: 'class '+s, filter(lambda s: s[0] is 'T', ''.join(formal_pars).split(','))))
    rank_par = ('','int %s'%self.rank)[self.rank is 'rank']
    return 'template <%s>'%(', '.join(filter(len,(rank_par,template_pars))))
  def makeAct(self, actN):
    act = CFact(self, actN)
    formal_pars = self.get_pars(act.podDatas)[0]
    shift_pars = map(lambda fp: fp[10:], filter(lambda fp: fp[:12] == 'const int _I', formal_pars))
    shift_line = ', '.join(map(lambda p: '%s=(%s<<%d)-%d'%(p[1:],p,self.get_dim(p[2:]),1<<self.get_dim(p[2:][:p[2:].index('p')])), shift_pars))
    #вычисление сдвигов поддатасов из имени сдвига датаса (например, из _Ixpd получаем Ixpd=(_Ixpd<<2)-1, где 2=dim(xpd), а 1=1<<dim(x).
    print self.getTmplPars(formal_pars)+' inline void %s(%s) {'%(self.actTmpl%actN,', '.join(formal_pars))
    if len(shift_line)>0:
      if ''.join(self.subactTmpl.keys()) in 'SF': print '//',
      print '  const int %s;'%shift_line
    if 'B' in self.subactTmpl.keys():
      for (a,n,tC) in act.PodActList('B'):
        print '  %s(%s);'%(self.subactTmpl['B']%a, ', '.join(CFpodact(act,a).getParsList(n,tC)))
    if 'F' in self.subactTmpl.keys():
      print '  %s(%s);'%(self.subactTmpl['F']%actN,', '.join(map(lambda p: p[self.par_name_start:], formal_pars)))
    if 'S' in self.subactTmpl.keys():
      tier=filter(lambda a:'-' not in a, makeCombinations([self.Rules2act['LR'].get(a, '-%c-'%a) for a in actN])); tier.reverse()
      for tactN in tier:
        tact = CFact(self, tactN)
        tformal_pars = self.get_pars(tact.podDatas, shift=self.list2num([{'I':-1,'Y':1}.get(c,0) for c in tactN]))[0]
        print '  %s(%s);'%(self.subactTmpl['S']%tactN,', '.join(map(lambda p: p[self.par_name_start:], tformal_pars)))
    if 'X' in self.subactTmpl.keys():
      for (a,n) in act.PodActList_mp():
        print '  %s(%s);'%(self.subactTmpl['X']%a, ', '.join(CFpodact_mp(act,a).getParsList(n,'X')))
    if 'T' in self.subactTmpl.keys():
      for (a,n,tC) in act.PodActList('T'): print '  %s(%s);'%(self.subactTmpl['T']%a, ', '.join(CFpodact(act,a).getParsList(n,tC)))
    print '}'
    pass
  def makeAct_mp(self, actN):
    actT = CFact_mp(self, 'D'*self.dim)
    act = CFact_mp(self, actN)
    formal_pars = self.get_pars_mp(act.podDatas_mp)[0]
    shift_pars = map(lambda fp: fp[10:], filter(lambda fp: fp[:12] == 'const int _I', formal_pars))
    shift_line = ', '.join(map(lambda p: '%s=(%s<<%d)%c%d'%(p[1:],p,self.get_dim(p[2:]),"-+-+"[p[2:].count('m')],1<<self.get_dim(p[2:][:p[2:].replace('m','p').index('p')])), shift_pars))
    #вычисление сдвигов поддатасов из имени сдвига датаса (например, из _Ixmd получаем Ixmd=(_Ixmd<<2)+1, где 2=dim(xmd), а 1=1<<dim(x).
    print self.getTmplPars(formal_pars)+' inline void %s(%s) {'%(self.actTmpl%actN,', '.join(formal_pars))
    if len(shift_line)>0: print '  const int %s;'%shift_line
    if 'B' in self.subactTmpl.keys():
      for (a,n,tC) in act.PodActList('B'): print '  %s(%s);'%(self.subactTmpl['B']%a, ', '.join(CFpodact(act,a).getParsList(n,tC)))
    #if 'F' in self.subactTmpl.keys():
    #  print '  %s(%s);'%(self.subactTmpl['F']%actN,', '.join(map(lambda p: p[self.par_name_start:], formal_pars)))
    #  #for (a,n,tC) in act.PodActList('F'): print '  %s(%s);'%(self.subactTmpl['F']%a, ', '.join(CFpodact(act,a).getParsList(n,tC)))
    if 'J' in self.subactTmpl.keys():
      caseNshift = self.list2num([1]*self.dim, 4)
      print '  for(int ic=0; ic<4; ic++) {'
      print 'if(dat0->datas[ic].Npts>NptsMax) {\n  T0& datT=dat0->datas[ic];\n  NptsMax = datT.Npts;\n  printf("'+'===%s:'%actN+' inc NptsMax to %d in xyt: %.3g %.3g %d\\n", NptsMax, datT.x, datT.y, datT.it);\n}'
      print '    int ip=dat0->datas[ic].Nexch;\n    while(ip < dat0->datas[ic].Npts) {\n      pts& pt=dat0->datas[ic].ptslist[ip];\n      double dstep=1.0;\n      do {\n       switch(pt.ix+4*pt.iy) {'
      for (a,n) in act.PodActList_mpPIC():
        nL = map(lambda nt: nt-1, self.num2list(n,4))
        parList = CFpodact_mp(act,a).getParsList(nL,'X')
        print '        case %d: dstep=pt.%s(dstep, %s); break;'%(n-caseNshift,self.subactTmpl['J']%a, ', '.join(parList))
      print '       }\n      } while(dstep<1.0);\n      if((pt.ix&2)|(pt.iy&2)) {\n        int swk=pt.ix+4*pt.iy;\n        if(pt.ix<0) pt.ix += 2; else if(pt.ix>1) pt.ix -= 2;\n        if(pt.iy<0) pt.iy += 2; else if(pt.iy>1) pt.iy -= 2;\n        switch(swk) {'
      for (a,n) in act.PodActList_mpPIC():
        nL = map(lambda nt: nt-1, self.num2list(n,4))
        if len(filter(lambda _n: _n in (-1,2), nL))==0: continue;
        datPtr = '(%s)'%CFpodact_mp(act,a).getParsList(nL,'X')[0]
        oldPtr = 'dat0->datas[ic]'
        Npts = datPtr+'->Npts'
        Nxch = datPtr+'->Nexch'
        case_dict = {'ptN': datPtr, 'ptO': oldPtr,'Np': Npts,'Nx':Nxch}
        print '          case %(n)d:'%{'n':n-caseNshift},
        #print 'if(%(ptO)s.it > %(ptN)s->it) printf("Illegal Exch!\\n"); else'%case_dict,
        print 'if(%(ptO)s.it > %(ptN)s->it) printf("Illegal Exch!\\n"); else if(%(ptO)s.it < %(ptN)s->it) %(ptN)s->ptslist[%(Np)s].copyfrom(pt); else { %(ptN)s->ptslist[%(Np)s].copyfrom(%(ptN)s->ptslist[%(Nx)s]); %(ptN)s->ptslist[%(Nx)s].copyfrom(pt); %(Nx)s++; } %(Np)s++; break;'%case_dict
        #print 'if(%(ptO)s.it < %(ptN)s->it) %(ptN)s->ptslist[%(Np)s].copyfrom(pt); else { %(ptN)s->ptslist[%(Np)s].copyfrom(%(ptN)s->ptslist[%(Nx)s]); %(ptN)s->ptslist[%(Nx)s].copyfrom(pt); %(Nx)s++; } %(Np)s++; break;'%case_dict
      print '        }\n        dat0->datas[ic].Npts--;\n        if(ip<dat0->datas[ic].Npts) pt.copyfrom(dat0->datas[ic].ptslist[dat0->datas[ic].Npts]);\n      } else ip++;\n    }\n    dat0->datas[ic].Nexch=0; dat0->datas[ic].it++;\n  }'
    if 'X' in self.subactTmpl.keys():
      for (a,n) in act.PodActList_mp():
        print '  for(int ip=0; ip<Nz; ip++) dat0->datas[%d].ptslist[ip].%s(1.0, %s);'%(n,self.subactTmpl['X']%a, ', '.join(CFpodact_mp(act,a).getParsList(n,'X')))
    if 'T' in self.subactTmpl.keys():
      for (a,n,tC) in act.PodActList('T'): print '  %s(%s);'%(self.subactTmpl['T']%a, ', '.join(CFpodact(act,a).getParsList(n,tC)))
    print '}'
    pass
  def genConeFold(self, rank='rank', actTmpl=r'%sactCF', subactTmpl=None, knot='p', exclude=[], acts4gen=[]):
    '''Печатает ConeFold заданного ранга, имени и типа:
    rank --- имя ранга (либо диапазона), строка;
    actTmpl --- шаблон имени, строка, на которую накатывается имя act-а;
    subactTmpl --- правила разбиения и имена подConeFold-ов (мЕньшего ранга), на которые разбивается ConeFold, словарь, ключи которого ---
    слои по времени B/T/F/X --- bottom/top/flat/flat с , уровни подConeFold-ов'''
    if len(exclude): print '// exclude up to: %d'%len(exclude)
    self.rank=rank
    self.actTmpl=actTmpl
    if subactTmpl is None: subactTmpl = { 'BT' : actTmpl }
    if type(subactTmpl) is str: subactTmpl = { subactTmpl : actTmpl }
    self.subactTmpl={}; map(lambda k: self.subactTmpl.update(dict(zip(k,(subactTmpl[k],)*len(k)))), subactTmpl.keys())
    self.rules,self.rulesShift={},{}
    datasTypes = reduce(lambda r,t: r+filter(lambda c: c not in r+'_mp', self.Rules2act['pars'][t]), self.types, '')
    for pd in makeCombinations([datasTypes]*self.dim): self.add2rules(pd)
    if len(acts4gen) == 0: acts4gen = makeCombinations([self.types]*self.dim)
    for a in acts4gen:
      if a in exclude: continue
      if knot == 'p': self.makeAct(a)
      elif knot == 'mp': self.makeAct_mp(a)
      pass
    pass

#stdout, sys.stdout = sys.stdout,open('Test.inc.hpp', 'w')
#gPJ = Generator(dim=2, types='JDX')
#gPJ.genConeFold(rank="FFRank-1", actTmpl=r'PIC2update%s', subactTmpl={'J': r'PIC2update%s'}, knot='mp')

dim=3

incpath = sys.argv[0][:sys.argv[0].find(sys.argv[0].split('/')[-1])]
stdout, sys.stdout = sys.stdout,open(incpath+'CF2Dpic.inc.hpp', 'w')
g = Generator(dim=dim, types='IDX')
g.genConeFold(actTmpl=r'update%s')
g.genConeFold(actTmpl=r'FLDupdate%s')
#g.genConeFold(actTmpl=r'FLDupdate%s', subactTmpl={'B': r'FLDupdate%s'})
sys.stdout.close(); sys.stdout = stdout

stdout, sys.stdout = sys.stdout,open(incpath+'CF2Dpic.inc.hpp', 'w')
print 'int NptsMax=0;'
gP = Generator(dim=dim, types='IDX')
gPJ = Generator(dim=dim, types='JDX')
g.genConeFold(actTmpl=r'picNfld%s', subactTmpl={'B': r'picNfld%s', 'T':r'FLDupdate%s'})
gP.genConeFold(rank="PicRank+1", actTmpl=r'update%s', subactTmpl={'B': r'picNfld%s', 'T':r'FLDupdate%s'})
gP.genConeFold(rank="FFRank-1", actTmpl=r'PIC1update%s', subactTmpl={'X': r'PIC1update%s'})
gPJ.genConeFold(rank="FFRank-1", actTmpl=r'PIC2update%s', subactTmpl={'J': r'PIC2update%s'}, knot='mp')
gP.genConeFold(rank="FFRank-1", actTmpl=r'PIC3update%s', subactTmpl={'X': r'PIC3update%s'})
gP.genConeFold(rank="FFRank", actTmpl=r'pic%s', subactTmpl={'B': r'PIC1update%s', 'X':r'PIC2update%s', 'T':r'PIC3update%s'})
gP.genConeFold(rank="FFRank+1", actTmpl=r'picNfld%s', subactTmpl={'B': r'pic%s', 'T':r'FLDupdate%s'})
sys.stdout.close(); sys.stdout = stdout

#---------------PML-------------------------------------------
stdout, sys.stdout = sys.stdout,open(incpath+'CF2Dpic.inc.hpp', 'w')

#============= Non-PML ConeFold
g = Generator(dim=dim, types='DX')
#acts = ['DD','DX']
acts = ['D'*dim]

#============= BC ConeFold for rank>PMLrank
gBC = Generator(dim=dim, types='LDRX')
gBC.Rules2rank['d'*dim] = 'rank+PMLrank'
#gBC.Rules2rank['dx'] = 'rank+PMLrank'
actsBC = makeCombinations(['LDR','LDR','LDX'])

#============= PML ConeFold for rank<=PMLrank
gPML = Generator(dim=dim, types='ILDSRYX')
gPML.Rules2act['subacts'].update({'I':'-IIS', 'Y':'SYY-', 'S':'SSSS', 'L':'SLLD', 'R':'DRRS'})
gPML.Rules2act['pars'].update({'S':'sp', 'L':'sd', 'R':'ds', 'Y':'s_', 'I':'_s', 'J':'sp'})
for s in 'lr': gPML.Rules2dim[s] = 1
actsPML = makeCombinations(['ILDSRY','ILDSRY','ILDSX'])

print '//===========any rank==============DX'
g.genConeFold(actTmpl=r'%sact', acts4gen=acts)

print '//=========rank>PMLrank============LDR/X'
gBC.genConeFold(actTmpl=r'%sact', acts4gen=actsBC, exclude=acts)

print '//=========rank<=PMLrank============DX'
g.genConeFold(actTmpl=r'%sactPML', subactTmpl={'F':r'%sact'}, acts4gen=acts)

print '//=========rank<PMLrank============ILDSRY'
gPML.genConeFold(actTmpl=r'%sactPML', acts4gen=actsPML, exclude=acts)

print '//=========rank=PMLrank============LDR/X'
gPML.genConeFold(rank='PMLrank', actTmpl=r'%sact', subactTmpl={'S':r'%sactPML'}, acts4gen=actsBC, exclude=acts)

sys.stdout.close(); sys.stdout = stdout
