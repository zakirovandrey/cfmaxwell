import os, socket

env = Environment(
  tools = ['g++', 'swig', 'default'],
  CPPPATH = ['.', '/usr/src/linux/include'],
#  CPPFLAGS = ['-march=native', '-mtune=native', '-w', '-fPIC', '-time', '-O3'],
  CPPFLAGS = ['-w', '-fPIC', '-time', '-std=c++0x'],
  CPPDEFINES = [],
  LIBS = ['aiv'],
)

env['BUILDERS']['genIncs'] = Builder(action= 'python $SOURCES', src_suffix = '.py', suffix = '.inc.hpp')

#env.ParseConfig("pkg-config python --libs --cflags")
#env.ParseConfig("python-config --cflags")
env.ParseConfig("python-config --libs")
env.ParseConfig("python-config --includes")
conf = Configure(env)
for lib in ["libpthread", "libgomp", "libboost_serialization"]:
  if not conf.CheckLib(lib): print "You need %s to compile CFmaxwell"%lib; Exit(1)
if conf.CheckLibWithHeader("numa","numa.h",'c'): print "Compiling for NUMA"; conf.env.Append(CPPFLAGS = ['-DNUMA'])
env = conf.Finish()
env['LIBS'].reverse() # -libnuma goes first
env.Append(LIBS = ['pthread', 'gomp', 'boost_serialization'])
env.Append(CPPFLAGS = ['-fopenmp'])
#if(socket.gethostname()=='D'):        env.Append(CPPFLAGS = ['-march=native', '-mtune=native'])
if(socket.gethostname()=='D'):        env.Append(CPPFLAGS = ['-march=core-avx-i', '-mtune=core-avx-i'])
if(socket.gethostname()=='photon'):   env.Append(CPPFLAGS = ['-march=bdver2', '-mtune=bdver2'])
if(socket.gethostname()=='electron'): env.Append(CPPFLAGS = ['-march=amdfam10', '-mtune=amdfam10'])
if(socket.gethostname()=='positron'): env.Append(CPPFLAGS = ['-march=amdfam10', '-mtune=amdfam10'])
if(socket.gethostname()=='O48.cluster.local'): env.Append(CPPFLAGS = ['-march=amdfam10', '-mtune=amdfam10'])

vars = Variables('build-setup.conf')
vars.Add('MATERIALS','Set .cpp file with Material(...) function','materials/materials.cpp')
vars.Add('SIGNAL','Set .cpp file with Signal(...) function','signals/Signal.cpp')
vars.Update(env)

config_hpp_defines = {"MaxRank":"8", "nLArank":"4", "PMLrank":"1"}
#def config_hpp_build(target, source, env):
#  for a_target, a_source in zip(target, source):
#    config_hpp = file(str(a_target), "w")
#    build_conf = file(str(a_source), "r")
#    config_hpp.write('#ifndef CONFIG_HPP\n#define CONFIG_HPP\n')
#    for val in build_conf.read().split('\n'):
#      sval = val.upper().split("=")
#      if len(sval)==2: config_hpp_defines.update( dict((''.join(val.upper().split()).split("="),)) )
#    for val in ["MAXRANK","NLARANK","PMLRANK"]:
#      config_hpp.write( "const int %s = %s;\n"%({"MAXRANK":"MaxRank","NLARANK":"nLArank","PMLRANK":"PMLrank"}[val],config_hpp_defines[val]) )
#    config_hpp.write('#include "params.hpp"\n')
#    config_hpp.write("#endif")
#    config_hpp.close()
#    build_conf.close()
#env.AlwaysBuild(env.Command("build/config.hpp","build-setup.conf",config_hpp_build,duplicate=1))

print socket.gethostname()
if(socket.gethostname()=='positron'): 
  env.Replace(CXX = 'g++-4.6'); env.Append(CPPFLAGS = ['-Ofast'])
if(socket.gethostname() in ['photon','electron','D',"O48.cluster.local"]):
  env.Replace(CXX = 'g++-4.7'); env.Append(CPPFLAGS = ['-Ofast'])
else:
  env.Append(CPPFLAGS = ['-O3'])

if not os.path.exists("build/ConeFold.py"): os.system("cp src/ConeFold.py build/")
SConscript('src/SConscript', variant_dir='build', exports=['env'], duplicate=1)

