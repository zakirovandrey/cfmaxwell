#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from types import *
from operator import *

def makeCombinations(lst, res=['']):
  if len(lst)<=0: return res
  else: return makeCombinations(lst[:-1], [c+n for n in res for c in lst[-1]])

def makeCombinationsList(lst, res=[[]]):
  if len(lst)<=0: return res
  else: return makeCombinationsList(lst[:-1], [[c]+n for n in res for c in lst[-1]])
  pass

def num_bits(i):
  if i==0: return 0
  return (i&1)+num_bits(i>>1)

class CFact:
  def __init__(self, Gen, actN):
    self.Gen=Gen; self.actN=actN
    self.podDatas,self.podDatasShift = ['',],[0,]
    for s in xrange(Gen.dim):
      #datas = self.pars4actN[actN[s]]
      datas = Gen.Rules4act['pars'][actN[s]]
      for k in xrange(1<<s):
        self.podDatas.append(self.podDatas[k] + datas[1])
        self.podDatas[k] += datas[0]
        self.podDatasShift.append(self.podDatasShift[k] + (0,1<<s)[datas[1] is 'p'])
        pass
      pass
    pass
  def PodActList(self, tC='#'):
    tier0=makeCombinations([self.Gen.Rules4act['subacts'][a][:2] for a in self.actN])
    tier1=makeCombinations([self.Gen.Rules4act['subacts'][a][2:] for a in self.actN])
    tier0=zip(tier0,xrange(len(tier0)),['_']*len(tier0))
    tier1=zip(tier1,xrange(len(tier1)),['~']*len(tier1))
    tier2=filter(lambda (a,n,t):'-' not in a, list(reduce(add,zip(tier1,tier0)))); tier2.reverse()
    tier0=filter(lambda (a,n,t):'-' not in a, tier0); tier0.reverse()
    tier1=filter(lambda (a,n,t):'-' not in a, tier1); tier1.reverse()
    if tC is '_': return tier0
    if tC is '~': return tier1
    if tC is '#': return tier2
    return tier0+tier1
  def PodActList_pm(self):
    tier=makeCombinations([self.Gen.Rules4act['subacts_pm'][a] for a in self.actN])
    tier=filter(lambda (a,n):'-' not in a, zip(tier,xrange(len(tier)))); tier.reverse()
    return tier
  def getIsh(self, datI, shI):
    return filter(len, [('',self.Gen.rulesShift[self.podDatas[datI]][s])[(shI&(1<<s))>0] for s in xrange(self.Gen.dim)])
  def getPodActPar(self, par, parI, cntI, tC):
    '''Возвращает параметр parI (0..2^d) для под-act-а с базовым datas-ом cntI (0..2^d), исходя из шаблона par'''
    if tC is '_': datI,poddatI = cntI&parI,cntI^parI
    if tC is '~': datI,poddatI = cntI|parI,((1<<self.Gen.dim)-1)&(~(cntI^parI))
    if tC is '$':
      parIl,cntIl=map(lambda v:(0,1,-1)[v],self.Gen.num2list(parI,3)),self.Gen.num2list(cntI)
      datI = self.Gen.list2num([(p+n)>=1 for (p,n) in zip(parIl,cntIl)])
      poddatI = self.Gen.list2num([1&~((p&1)^n) for (p,n) in zip(parIl,cntIl)])
    #datI,poddatI --- номера датаса и его поддатаса для параметра parI
    poddatIx = 0
    for s in xrange(self.Gen.dim-1,-1,-1):
      i2s = 1<<s
      if self.Gen.Rules4dim[self.podDatas[datI][s]]: poddatIx = 2*poddatIx+((poddatI&i2s)!=0)
    #poddatIx --- смещение параметра в массиве поддатасов для датаса datI
    #print par, parI, cntI, tC,self.podDatas[datI],"->",datI,poddatI,poddatIx
    if par.find('dat')>=0:
      dat_shift = self.podDatasShift[datI]
      return (self.Gen.Rules4cpp['dat']%datI+'->',self.Gen.Rules4cpp['dat']%(datI-dat_shift)+'[%s].'%('+'.join(self.getIsh(datI-dat_shift,dat_shift))))[dat_shift>0]+'datas'+('','+%d'%poddatIx)[poddatIx>0]
    if par.find(self.Gen.Rules4cpp['preI'])>=0:
      sh = par[par.find(self.Gen.Rules4cpp['preI'])+len(self.Gen.Rules4cpp['preI']):]
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

class CFact_pm(CFact):
  def __init__(self, Gen, actN):
    CFact.__init__(self, Gen, actN)
    self.podDatas_pm,self.podDatasShift_pm = ['',],[0,]
    for s in xrange(Gen.dim):
      datas = Gen.Rules4act_pm['pars'][actN[s]]
      for k in xrange(3**s): self.podDatas_pm.append(self.podDatas_pm[k] + datas[2])
      for k in xrange(3**s): self.podDatas_pm.append(self.podDatas_pm[k] + datas[0])
      for k in xrange(3**s): self.podDatas_pm[k] += datas[1]
      for k in xrange(3**s): self.podDatasShift_pm.append(self.podDatasShift_pm[k] + (0,3**s)[datas[2] is 'p'])
      for k in xrange(3**s): self.podDatasShift_pm.append(self.podDatasShift_pm[k] + 2*(0,3**s)[datas[0] is 'm'])
      pass
    pass
  def PodActList_pm(self):
    tier=makeCombinations([self.Gen.Rules4act_pm['subacts'][a] for a in self.actN])
    tier=filter(lambda (a,n):'-' not in a, zip(tier,xrange(len(tier)))); tier.reverse()
    return tier
  def PodActList_pmPIC(self):
    tier=makeCombinations([self.Gen.Rules4act_pm['subactsPIC'][a] for a in self.actN])
    tier=filter(lambda (a,n):'-' not in a, zip(tier,xrange(len(tier)))); tier.reverse()
    return tier
  def getIsh_pm(self, datI, shI):
    '''вычисляет сдвиг datas-а'''
    if '_' in self.podDatas_pm[datI]: return ('<no-data>',)
    shIl = self.Gen.num2list(shI,3)
    Ish=[]
    for s,sh in zip(xrange(self.Gen.dim), self.Gen.num2list(shI,3)):
      rul = self.Gen.rulesShift[self.podDatas_pm[datI]][s]
      Ish.append(('',rul, rul.replace('p','m'))[shIl[s]])
    return filter(len, Ish)
  def getPodActPar(self, par, parI, cntIl, tC):
    '''Возвращает параметр parI (0..3^d) для под-act-а с базовым datas-ом cntI (0..2^d), исходя из шаблона par'''
    if tC is '$':
      parIl=map(lambda v:(0,1,-1)[v],self.Gen.num2list(parI,3))
      datIl = [(-1,-1,0,0,1,1)[p+n+2] for (p,n) in zip(parIl,cntIl)] # можно (p+n)/2 # номер (смещение) datas-а относительно базового cntI
      poddatIl = [(0,1,0,1,0,1)[p+n+2] for (p,n) in zip(parIl,cntIl)] # можно (p+n)%2 # номер poddatas-а в datas-е
      datI = self.Gen.list2num([((p+n)/2)%3 for (p,n) in zip(parIl,cntIl)], 3)
      poddatI = self.Gen.list2num([((p+n)%2)&3 for (p,n) in zip(parIl,cntIl)])
    poddatIx = 0
    for s in xrange(self.Gen.dim-1,-1,-1):
      i2s = 1<<s
      if self.Gen.Rules4dim[self.podDatas_pm[datI][s]]: poddatIx = 2*poddatIx+((poddatI&i2s)!=0)
    if par.find('dat')>=0:
      dat_shift = self.podDatasShift_pm[datI]
      return (self.Gen.Rules4cpp['dat']%datI+'->',self.Gen.Rules4cpp['dat']%(datI-dat_shift)+'[%s].'%('+'.join(self.getIsh_pm(datI-dat_shift,dat_shift))))[dat_shift!=0]+'datas'+('','+%d'%poddatIx)[poddatIx>0]
    if par.find(self.Gen.Rules4cpp['preI'])>=0:
      sh = par[par.find(self.Gen.Rules4cpp['preI'])+len(self.Gen.Rules4cpp['preI']):]
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

class CFpodact_pm(CFact_pm):
  def __init__(self, act, actN):
    CFact_pm.__init__(self, act.Gen, actN)
    self.nadact = act
  def getParsList(self, n, tC):
    pars_list, full_pars_list = self.Gen.get_pars_pm(self.podDatas_pm)
    fakt_pars = [self.nadact.getPodActPar(par, full_pars_list.index(par), n, tC) for par in pars_list]
    return fakt_pars

class Generator:
  def __init__(self, dim, formaABC):
    '''Создает генератор ConeFold-ов.
    dim --- размерность декомпозиции
    formaABC --- список символов алфавита, входящих в генерируемые формы ConeFold-ов.
    после создания объекта генератора можно его настроить, задав правила (словарь словарей):
    Набор правил (разумные правила заданы по умолчанию) находится в двух словарях Rules4act и Rules4act_pm
    правила следующие:
    Rules4act,pm['pars'] --- правила определения формальных и фактических параметров функций-методов
    Rules4act['subacts']
    Rules4act['subacts_pm']
    Rules4act[]
    '''
    self.Rules4act = {
      'pars': {'D':'dp', 'S':'dp', 'I':'_d', 'J':'dp', 'X':'dx', 'Y':'d_', 'P':'dd', 'Q':'dd', 'L':'ld', 'M':'dp', 'R':'dr'},
      'subacts':{'D':'DDDD', 'S':'SSSS', 'I':'-IID', 'X':'DXX-', 'Y':'DYY-', 'P':'DPPD', 'L':'-LLD', 'R':'DRR-'},
      'subacts_pm': {'D':'DD', 'I':'-J', 'J':'JD', 'X':'X-', 'Y':'Y-', 'P':'PQ', 'L':'-M', 'R':'R-'},
      'LR':{'L':'IL-', 'R':'-RY'},
      'CompositSeq': '<_#=&$@~>',
    }
    self.Rules4act_pm = {
      #'pars': {'D':'mdp', 'S':'msp', 'J':'_dp', 'X':'mdp', 'P':'mdd', 'Q':'ddp', 'M':'ldp', 'R':'mdr'},
      'pars': {'D':'mdp', 'S':'msp', 'J':'_dp', 'X':'mdx', 'P':'mdd', 'Q':'ddp', 'M':'ldp', 'R':'mdr'},
      'subacts': {'D':'DD', 'I':'JD', 'J':'JD', 'X':'DX', 'P':'DP', 'Q': 'QD', 'L':'MD', 'M':'MD', 'R':'DR'},
      'subactsPIC': {'D':'DDDD', 'J':'-JDD', 'X':'DDX-', 'P':'DDPQ', 'Q': 'PQDD', 'M':'-MDD', 'R':'DDR-'},
      'CompositSeq': '<_#=&$@~>',
    }
    self.Rules4rank = {}
    self.Rules4dim = {}
    par_size = 2+dim+(dim==1)
    ind_size = 1+(dim-1)/3
    self.Rules4cpp = {
      'preact': 'template <{pars}> inline void ',
      'act': '{name}{form}({pars})',
      'preI': 'const int _I',
      'I': '_I',
      'dat': 'datas_____'[:par_size-ind_size:]+'%0'+'%dd'%ind_size,
      'add_funcs': {},
    }
    for s in 'dmpse': self.Rules4dim[s] = 1
    for s in 'xlr': self.Rules4dim[s] = 0
    self.dim=dim
    self.par_name_start=-par_size
    self.formaABC=formaABC
    print '// acts: %d'%len(makeCombinations([formaABC]*dim))
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
    return reduce(lambda r,s: r+self.Rules4dim[s], pd, 0)
  def add2rules(self, pd):
    self.rules[pd] = 'cubeLR<%d,T%%(Npd)d,%s>'%(self.get_dim(pd),self.Rules4rank.get(pd,self.rank)) + '* const '+self.Rules4cpp['dat'].replace('%','%(Npd)')
    self.add2rulesShift(pd)
    pass
  def add2rulesShift(self, pd):
    self.rulesShift[pd] = []
    #self.rulesShiftM[pd] = []
    for s in xrange(len(pd)):
      sh=pd[:s]+'p'+pd[s+1:]
      self.rules[sh] = self.Rules4cpp['preI']+sh
      self.rulesShift[pd].append(self.Rules4cpp['I']+sh)
      shM=pd[:s]+'m'+pd[s+1:]
      self.rules[shM] = self.Rules4cpp['preI']+shM
      #self.rulesShiftM[pd].append('_I'+shM)
      pass
    pass
  def add2rules_nLA(self, pd):
    pdd = pd.replace('p','d').replace('m','d').replace('e','d').replace('c','d').replace('x','r')
    self.rules[pd] = 'get_%s<%s,%%(Npd)d>(net,{%%(Npd)d})'%(pdd,self.Rules4rank.get(pd,self.rank))
    self.add2rulesShift_nLA(pd)
    pass
  def add2rulesShift_nLA(self, pd):
    self.rulesShift[pd] = []
    #self.rulesShiftM[pd] = []
    for s in filter(lambda s: pd[s]=='d', xrange(len(pd))):
      pddim = self.get_dim(pd)
      #pddim = self.get_dim(pd[:s]+'d'+pd[s+1:])
      sh=pd[:s]+'p'+pd[s+1:]
      if pddim == self.dim:
        self.rules[sh] = 'IM{0}-I{0}'.format(s)
        srep = '-I{0}'.format(s)
        self.rules4replace_nLA[s][srep] = '+(1<<({0}*({2})))*MultiBrickShift[{1}]'.format(pddim,s,self.Rules4rank.get(pd,self.rank))+srep
      else:
        nind = pd.find('d') if pddim==1 else len(pd)-len(pd.replace('d',' ').lstrip())
        IndStr = 'ind{dim}_{d}<{r},{nind}>'.format(dim=self.dim,d=pddim,r=self.Rules4rank.get(pd,self.rank),nind=nind)
        self.rules[sh] = '/*{c}*/{ind}(IM{s})-{ind}(I{s})'.format(ind=IndStr, s=s, c=pd)
        srep = '-{ind}(I{s})'.format(ind=IndStr, s=s)
        self.rules4replace_nLA[s][srep] = '+(1<<({0}*({2})))*MultiBrickShift[{1}]'.format(pddim,s,self.Rules4rank.get(pd,self.rank))+srep
        #print sh,self.rules[sh]
      #self.rules[sh] = '_%d_'%pddim+sh
      self.rulesShift[pd].append(self.Rules4cpp['I']+sh)
      shM=pd[:s]+'m'+pd[s+1:]
      self.rules[shM] = '_'+shM
      #self.rulesShiftM[pd].append('_I'+shM)
      pass
    pass
  def get_pars(self, datas, shift=0):
    full_pars_list = [self.rules.get(pd,'')%{'Npd':i+shift} for (i,pd) in zip(xrange(1<<self.dim),datas)]
    pars_list = filter(len, full_pars_list)
    pars_list = map(lambda i: pars_list[i], filter(lambda i: pars_list.index(pars_list[i])==i, xrange(len(pars_list))))
    return pars_list, full_pars_list
  def get_pars_pm(self, datas):
    full_pars_list = [self.rules.get(pd,'')%{'Npd':i} for (i,pd) in zip(xrange(3**self.dim),datas)]
    pars_list = filter(len, full_pars_list)
    pars_list = map(lambda i: pars_list[i], filter(lambda i: pars_list.index(pars_list[i])==i, xrange(len(pars_list))))
    return pars_list, full_pars_list
  def getTmplPars(self, formal_pars):
    template_pars = ','.join(map(lambda s: 'class '+s, filter(lambda s: s[0] is 'T', ''.join(formal_pars).split(','))))
    rank_par = ('','int %s'%self.rank)[self.rank is 'rank']
    return ', '.join(filter(len,(rank_par,template_pars)))
    #return 'template <%s>'%(', '.join(filter(len,(rank_par,template_pars))))
  def checkNmakeSubActsTier(self, ck, act, formal_pars):
    if ck not in self.composit.keys(): return
    if ck in '_~#':
      for (a,n,tC) in act.PodActList(ck):
        print '  '+self.Rules4cpp['act'].format(name=self.composit[ck], form=a, pars=', '.join(CFpodact(act,a).getParsList(n,tC)))+';'
    elif ck == '=':
      print '  '+self.Rules4cpp['act'].format(name=self.composit[ck], form=act.actN, pars=', '.join(map(lambda p: p[self.par_name_start:], formal_pars)))+';'
    elif ck == '&':
      tier=filter(lambda a:'-' not in a, makeCombinations([self.Rules4act['LR'].get(a, '-%c-'%a) for a in act.actN])); tier.reverse()
      for tactN in tier:
        tact = CFact(self, tactN)
        tformal_pars = self.get_pars(tact.podDatas, shift=self.list2num([{'I':-1,'Y':1}.get(c,0) for c in tactN]))[0]
        print '  '+self.Rules4cpp['act'].format(name=self.composit[ck], form=tactN, pars=', '.join(map(lambda p: p[self.par_name_start:], tformal_pars)))+';'
    elif ck == '$':
      for (a,n) in act.PodActList_pm():
        print '  '+self.Rules4cpp['act'].format(name=self.composit[ck], form=a, pars=', '.join(CFpodact_pm(act,a).getParsList(n,ck)))+';'
    elif ck == '@':
      apply(self.Rules4cpp['add_funcs'][self.composit[ck]], (self, act))
    pass
  def checkNmakeSubActsTier_pm(self, ck, act, formal_pars):
    if ck not in self.composit.keys(): return
    if ck in '_~#':
      for (a,n,tC) in act.PodActList(ck):
        print '  '+self.Rules4cpp['act'].format(name=self.composit[ck], form=a, pars=', '.join(CFpodact(act,a).getParsList(n,tC)))+';'
    elif ck == '$':
      for (a,n) in act.PodActList_pm():
        print '  for(int ip=0; ip<Nz; ip++) dat0->datas[%d].ptslist[ip].%s;'%(n,self.Rules4cpp['act'].format(name=self.composit[ck], form=a, pars='1.0, '+', '.join(CFpodact_pm(act,a).getParsList(n,ck))))
    elif ck == '=':
      pass
    #if 'F' in self.composit.keys():
    #  print '  %s(%s);'%(self.composit['F']%actN,', '.join(map(lambda p: p[self.par_name_start:], formal_pars)))
    #  #for (a,n,tC) in act.PodActList('F'): print '  %s(%s);'%(self.composit['F']%a, ', '.join(CFpodact(act,a).getParsList(n,tC)))
    elif ck == '@':
      apply(self.Rules4cpp['add_funcs'][self.composit[ck]], (self, act))
    pass
  def makeAct(self, actN):
    act = CFact(self, actN)
    formal_pars = self.get_pars(act.podDatas)[0]
    shift_pars = map(lambda fp: fp[10:], filter(lambda fp: fp[:12] == self.Rules4cpp['preI'], formal_pars))
    shift_line = ', '.join(map(lambda p: '%s=(%s<<%d)-%d'%(p[1:],p,self.get_dim(p[2:]),1<<self.get_dim(p[2:][:p[2:].index('p')])), shift_pars))
    #вычисление сдвигов поддатасов из имени сдвига датаса (например, из _Ixpd получаем Ixpd=(_Ixpd<<2)-1, где 2=dim(xpd), а 1=1<<dim(x).
    print self.Rules4cpp['preact'].format(pars=self.getTmplPars(formal_pars))+self.Rules4cpp['act'].format(name=self.actName, form=actN, pars=', '.join(formal_pars))+' {'
    if len(shift_line)>0:
      if ''.join(self.composit.keys()) in '&=': print '//',
      print '  const int %s;'%shift_line
    for ck in self.Rules4act['CompositSeq']: self.checkNmakeSubActsTier(ck, act, formal_pars)
    print '}'
  def makeAct_pm(self, actN):
    actT = CFact_pm(self, 'D'*self.dim)
    act = CFact_pm(self, actN)
    formal_pars = self.get_pars_pm(act.podDatas_pm)[0]
    shift_pars = map(lambda fp: fp[10:], filter(lambda fp: fp[:len(self.Rules4cpp['preI'])] == self.Rules4cpp['preI'], formal_pars))
    shift_line = ', '.join(map(lambda p: '%s=(%s<<%d)%c%d'%(p[1:],p,self.get_dim(p[2:]),"-+-+"[p[2:].count('m')],1<<self.get_dim(p[2:][:p[2:].replace('m','p').index('p')])), shift_pars))
    #вычисление сдвигов поддатасов из имени сдвига датаса (например, из _Ixmd получаем Ixmd=(_Ixmd<<2)+1, где 2=dim(xmd), а 1=1<<dim(x).
    print self.Rules4cpp['preact'].format(pars=self.getTmplPars(formal_pars))+self.Rules4cpp['act'].format(name=self.actName, form=actN, pars=', '.join(formal_pars))+' {'
    if len(shift_line)>0: print '  const int %s;'%shift_line
    for ck in self.Rules4act_pm['CompositSeq']: self.checkNmakeSubActsTier_pm(ck, act, formal_pars)
    print '}'
  def makeAct_pmP(self, actN):
    print '//Not implement yet'
    pass
  def genConeFold(self, rank='rank', actName='actCF', composit=None, knot='+', exclude=[], acts4gen=[]):
    '''Печатает ConeFold заданного ранга, имени и типа:
    rank --- имя ранга (либо диапазона), строка;
    actName --- уникальное имя действия;
    composit --- правила разбиения и имена подConeFold-ов (мЕньшего ранга), на которые разбивается ConeFold, словарь, ключи которого ---
    #/_/~/=/$/@'''
    if len(exclude): print '// exclude up to: %d'%len(exclude)
    self.rank=rank
    self.actName=actName
    if composit is None: composit = { '#' : actName }
    if type(composit) is str: composit = { composit : actName }
    self.composit={}; map(lambda k: self.composit.update(dict(zip(k,(composit[k],)*len(k)))), composit.keys())
    #if '#' in composit.keys(): self.composit.update({'_': composit['#'],'~': composit['#']})
    self.rules,self.rulesShift={},{}
    datasTypes = reduce(lambda r,t: r+filter(lambda c: c not in r+'_pm', self.Rules4act['pars'][t]), self.formaABC, '')
    for pd in makeCombinations([datasTypes]*self.dim): self.add2rules(pd)
    if len(acts4gen) == 0: acts4gen = makeCombinations([self.formaABC]*self.dim)
    for a in acts4gen:
      if a in exclude: continue
      if knot in ('+','p'): self.makeAct(a)
      elif knot in ('-+','+-','pm','mp'): self.makeAct_pm(a)
      elif knot in ('+-+','-++','pmP','mpP'): self.makeAct_pmP(a)
      else: print 'knot','"%s"'%a,'is not implement, variants: +,+-,+-+'
      pass
    pass
  def gen4nLAget_data(self, dataName='data', rank='rank', nLArank='MaxRank-rank', dataABC='DLR'):
    self.rank = rank
    self.dataName = dataName
    datas4gen = makeCombinations([dataABC]*self.dim)
    print '// datas:',len(datas4gen)
    for d in datas4gen:
      ND=d.count('D')
      dataN = d.lower().replace('x','r').replace('i','l')
      IndShift = ''
      if ND>0:
        IndShift = '+((1<<%d*%s)-1'%(ND,self.rank)
        if ND==self.dim: IndShift += '-I)'
        else:
          s = d.find('D') if ND==1 else len(d)-len(d.replace('D',' ').lstrip())
          IndShift += '-ind%d_%d<%s,%d>(I))'%(self.dim,ND, self.rank,s)
      print 'template <int {r},int ind> cubeLR<{ND},{dN}{d},{lr}>* get_{dn}(nLAnet<{dim}>* net, unsigned int I) {{'\
        ' return &(((brick<{dim},{ND},{dN}{d},{lr}>*)net->datas[ind])->data){IndShift}; }}'.format(
        r=self.rank, lr=self.Rules4rank.get(dataN,nLArank), dim=self.dim, dN=self.dataName, d=d, ND=ND, dn=dataN,IndShift=IndShift
        )
    pass
  def gen4nLA(self, rank='rank', actName='actCF', knot='+', LinkType='LCR', lcr=None):
    '''Печатает вставку в netnLA для вызова всех возможных ConeFold-ов:
    rank --- имя ранга, строка;
    actTmpl --- шаблон имени, строка, на которую накатывается имя act-а;'''
    self.LinkType = LinkType
    if lcr is None: self.lcr=LinkType
    else: self.lcr=lcr
    if type(self.lcr) is str: self.lcr = (self.lcr,)*3
    self.rank=rank
    self.actName=actName
    datasTypes = reduce(lambda r,t: r+filter(lambda c: c not in r+'_pm', self.Rules4act['pars'][t]), self.formaABC, '')
    Dcart = range(1<<self.dim)
    Dcart.sort(lambda a,b: num_bits(a)-num_bits(b))
    self.rules,self.rulesShift={},{}
    self.rules4replace_nLA = [{} for s in xrange(self.dim)]
    for pd in makeCombinations([datasTypes]*self.dim): self.add2rules_nLA(pd)
    #print '//',self.rules
    for iD in Dcart:
      if iD == 0: print 'if(iX == 0) {'
      elif iD == (1<<self.dim)-1: print '} else {'
      else: print "}} else if((iX&{0})==0) {{".format((1<<self.dim)-1-iD)
      self.print4nLA(''.join(['D?'[(iD>>s)&1] for s in xrange(self.dim)]))
      pass
    print '}'
    pass
  def print4nLA(self, actN, sh=0, spN=1, replace_rules=[]):
    #print '//===',actN, sh, spN, replace_rules
    rules_nLA = {
      3: ('I','I&~M |I0','I&~M1|I1','I0|I1|IM2','I&~M2|I2','I0|IM1|I2','IM0|I1|I2','I0|I1|I2'),
      2: ('I','I0|IM1','IM0|I1','I0|I1'), 1: ('I','I0')
    }
    if sh >= len(actN):
      act = CFact(self, actN)
      formal_pars = self.get_pars(act.podDatas)[0]
      #if actN == 'PX': print r'printf("I am PX, I=%d,%d,%d,%d; IM=%d,%d, I01=%d,%d\n",I,IM0-I0,IM0|I1,ind2_1<rank,0>(IM0)-ind2_1<rank,0>(I0), IM0,IM1, I0,I1);'
      #print '//',actN, act.podDatas, formal_pars
      #print '//',actN, self.get_pars(act.podDatas)
      ParsStr=(','.join(formal_pars)).format(*rules_nLA[self.dim])
      for rr in replace_rules: ParsStr=ParsStr.replace(rr[0],rr[1])
      print ' '*spN,self.Rules4cpp['act'].format(name=self.actName, form=actN, pars=ParsStr)+';'
      pass
    elif actN[sh] is 'D': self.print4nLA(actN, sh+1, spN, replace_rules)
    elif actN[sh] is '?':
      #print ' '*spN,'if(net->indL&%s) {'%(1<<sh)
      #self.print4nLA(actN[:sh]+self.lcr[sh][0]+actN[sh+1:], sh+1, spN+2)
      #print ' '*spN,'} else if(net->indR&%s) {'%(1<<sh)
      #self.print4nLA(actN[:sh]+self.lcr[sh][2]+actN[sh+1:], sh+1, spN+2)
      #print ' '*spN,'} else {'
      #self.print4nLA(actN[:sh]+self.lcr[sh][1]+actN[sh+1:], sh+1, spN+2)
      #print ' '*spN,'}'
      print ' '*spN,"if(net->linkType[%d]=='%c') {"%(sh,self.LinkType[1])
      if self.lcr[sh][1]=='D':
        #print ' '*spN,'  I%d -= 1<<(dim*rank);'%sh+'//'+actN
        #pd='d'*self.dim
        #self.rules[pd[:sh]+'p'+pd[sh+1:]] += '+(1<<(dim*rank))'
        self.print4nLA(actN[:sh]+self.lcr[sh][1]+actN[sh+1:], sh+1, spN+2, replace_rules + self.rules4replace_nLA[sh].items())
      else:
        self.print4nLA(actN[:sh]+self.lcr[sh][1]+actN[sh+1:], sh+1, spN+2, replace_rules)
      print ' '*spN,"} else if(net->linkType[%d]=='%c') {"%(sh,self.LinkType[0])
      self.print4nLA(actN[:sh]+self.lcr[sh][0]+actN[sh+1:], sh+1, spN+2, replace_rules)
      print ' '*spN,"} else if(net->linkType[%d]=='%c') {"%(sh,self.LinkType[2])
      self.print4nLA(actN[:sh]+self.lcr[sh][2]+actN[sh+1:], sh+1, spN+2, replace_rules)
      print ' '*spN,'}} else err_linkType(net->linkType,{0});'.format(sh)
    pass
