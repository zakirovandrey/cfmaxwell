#!/usr/bin/python
# -*- coding: utf-8 -*-

from ConeFold import *
import sys
#from types import *
#from operator import *

dim=3

#---------------PML-------------------------------------------
incpath = sys.argv[0][:sys.argv[0].find(sys.argv[0].split('/')[-1])]
stdout, sys.stdout = sys.stdout,open(incpath+'CF3Dpml.inc.hpp', 'w')

#============= Non-PML ConeFold
g = Generator(dim=dim, formaABC='DPX')
g.Rules4cpp.update({'act': '{form}{name}({pars})'})
#acts = ['DD','DX']
acts = makeCombinations(['DP']*dim)

#============= BC ConeFold for rank>PMLrank
gBC = Generator(dim=dim, formaABC='LDPRX')
gBC.Rules4cpp.update({'act': '{form}{name}({pars})'})
gBC.Rules4rank['ddd'] = 'rank+PMLrank'
gBC.Rules4rank['pdd'] = 'rank+PMLrank'
gBC.Rules4rank['dpd'] = 'rank+PMLrank'
gBC.Rules4rank['ddp'] = 'rank+PMLrank'
gBC.Rules4rank['dpp'] = 'rank+PMLrank'
gBC.Rules4rank['pdp'] = 'rank+PMLrank'
gBC.Rules4rank['ppd'] = 'rank+PMLrank'
gBC.Rules4rank['ppp'] = 'rank+PMLrank'
#gBC.Rules4rank['dx'] = 'rank+PMLrank'
#actsBC = makeCombinations(['LDPR','LDPR','LDPX'])
actsBC = makeCombinations(['LDPR']*3)

#============= PML ConeFold for rank<=PMLrank
gPML = Generator(dim=dim, formaABC='ILDPSRYX')
gPML.Rules4act['subacts'].update({'I':'-IIS', 'Y':'SYY-', 'S':'SSSS', 'L':'SLLD', 'R':'DRRS', 'P':'DPPD'})
gPML.Rules4act['pars'].update({'S':'sp', 'L':'sd', 'R':'ds', 'Y':'s_', 'I':'_s', 'J':'sp'})
gPML.Rules4cpp.update({'act': '{form}{name}({pars})'})
for s in 'lr': gPML.Rules4dim[s] = 1
#actsPML = makeCombinations(['ILDPSRY','ILDPSRY','ILDPSX'])
actsPML = makeCombinations(['ILDPSRY']*3)

print '//===========any rank==============DPX'
g.genConeFold(actName='act', acts4gen=acts)

g0=Generator(dim=dim, formaABC='DP')
g0.Rules4cpp.update({'act': '{form}{name}({pars})'})
#g0.genConeFold(actName='act_')
#g0.genConeFold(rank="1",actName='act_', composit={'#':"act"})

print '//=========rank>PMLrank============LDPR/X'
gBC.genConeFold(actName='act', acts4gen=actsBC)#, exclude=acts)

print '//=========rank<=PMLrank============DPX'
#g.genConeFold(actName='actPML', composit={'=':'act'}, acts4gen=acts)

print '//=========rank<PMLrank============ILDPSRY'
gPML.genConeFold(actName='actPML', acts4gen=actsPML)#, exclude=acts)

print '//=========rank=PMLrank============LDPR/X'
gPML.genConeFold(rank='PMLrank', actName='act', composit={'&':'actPML'}, acts4gen=actsBC)#, exclude=acts)

incpath = sys.argv[0][:sys.argv[0].find(sys.argv[0].split('/')[-1])]
stdout, sys.stdout = sys.stdout,open(incpath+'nLAnet.inc.hpp', 'w')
gnLA = Generator(dim=dim, formaABC='DPLRX')
gnLA.Rules4cpp.update({'act': '{form}{name}({pars})'})
gnLA.gen4nLA(actName='act', LinkType='LCR', lcr=('LPR','LPR','LPR'))
stdout, sys.stdout = sys.stdout,open(incpath+'nLAdat.inc.hpp', 'w')
gnLA.Rules4rank['ddd'] = 'MaxRank-rank'
gnLA.gen4nLAget_data(dataABC='LDR', nLArank='MaxRank-rank-PMLrank')

sys.stdout.close(); sys.stdout = stdout
