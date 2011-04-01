
# Convert project files between Visual Studio 2010 and 
# Visual C++ 2020 Express (use Python 2.6 or later).
#
# Copyright 2010, Dr B R Gladman

from __future__ import print_function

import os
import shutil
import string
import fileinput, sys

express = False

mkr = r'''  <Import Project="$(VCTargetsPath)\Microsoft.Cpp.Default.props" />''' + '\n'
pg1 = r'''  <PropertyGroup Condition="'$(Configuration)|$(Platform)'=='Release|x64'" Label="Configuration">''' '\n'
pg2 = r'''  <PropertyGroup Condition="'$(Configuration)|$(Platform)'=='Debug|x64'" Label="Configuration">''' '\n'
pgc = r'''    <ConfigurationType>Application</ConfigurationType>''' '\n'
pgt = r'''    <PlatformToolset>Windows7.1SDK</PlatformToolset>''' '\n'
pgf = r'''  </PropertyGroup>''' '\n'

def find_lines(l):
  i = i0 = i1 = i2 = 0
  while i < len(l):
    if l[i].strip() == mkr.strip():
      i0 = i + 1
    if (l[i].strip() == pg1.strip() and l[i + 1].strip() == pgc.strip() and 
        l[i + 2].strip() == pgt.strip() and l[i + 3].strip() == pgf.strip()):
      i1 = i
    if (l[i].strip() == pg2.strip() and l[i + 1].strip() == pgc.strip() and 
        l[i + 2].strip() == pgt.strip() and l[i + 3].strip() == pgf.strip()):
      i2 = i
    i = i + 1
  return (i0, i1, i2)
  
def add_it(fn):
  with open(fn, 'r') as f:
    l = f.readlines()
  t = find_lines(l)  
  if t[0] and not t[1] and not t[2]:
    print("Adding Express: ", fn)
    l[t[0] : t[0]] = [pg1, pgc, pgt, pgf, pg2, pgc, pgt, pgf]
    with open(fn, 'w') as f:
      f.writelines(l)

def remove_it(fn):
  with open(fn, 'r') as f:
    l = f.readlines()
  t = find_lines(l)
  if t[1] and t[2]:
    print("Removing Express: ", fn)
    l[t[2] : t[2] + 4] = []
    l[t[1] : t[1] + 4] = []
    with open(fn, 'w') as f:
      f.writelines(l)

if not os.getcwd().endswith("build.vc10"):
  print("This script must be run in the 'build.vc10' directory")
else:  
  for root, dirs, files in os.walk("./") :
    for file in files :
      if not (root.endswith("_p0") or root.endswith("_p3") or root.endswith("_p4")):
        if file.endswith(".vcxproj") :
          t = os.path.join(root, file)
          add_it(t) if express else remove_it(t)
