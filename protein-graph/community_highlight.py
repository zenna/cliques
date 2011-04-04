#!/usr/bin/env python
# community_highlight.py
# Copyright (c) 2010-11, Yun William Yu
# 
# Running this file in pymol, by "run /path/to/community_highlight.py",
# defines several new commands, which allow for easier visualisation of
# protein structural "communities".
#
#######################################################################
# Command descriptions:
#
# community_highlight(filename), where "filename" is the output from the
# "community" program of the stability package, will highlight each of
# the different communities in a different color. I have gone to a fair
# amount of trouble to ensure that the colors are as distinct as possible.
#
# multiple_pngs simply generates pngs from a directory full of the
# communities found by a stability analysis. These pngs can then be used
# to generate a move film using something like mencoder.
#
# file_select just selects all atom numbers specified by a text file.
#
#######################################################################
# Miscallaneous notes: 
# 
# Note that pymols "select" command takes the prefix "index" to select
# for atoms by its internal numbering scheme, while "id" selects based
# on an external numbering scheme. As you probably want the atom number-
# ing to match that of the PDB file, you should always use "id"!
#


from string import split 
from pymol import cmd 
import random
from colorsys import hsv_to_rgb

def community_highlight(filename):
  selections = cmd.get_names('selections')
  for i in range(len(selections)):
    cmd.delete(str( selections[i]))

  f = open(filename, 'r')
  a = []
  comm_id = 0
  for line in f :
    splitted = split(line)
    if (len(splitted) == 2):
      a.append([int(splitted[0]),int(splitted[1])])
      if (int(splitted[1])>comm_id): comm_id = int(splitted[1])
  for i in range(comm_id+1):
    cmd.select(i, 'none')

  b = []
  prev_node = a[0][0]
  prev_comm = a[0][1]
  length = 1
  for i in range(1,len(a)):
    if a[i][1] == prev_comm:
      length = length + 1
    else:
      b.append([a[i][0],length,prev_comm])
      length = 1
    prev_node = a[i][0]
    prev_comm = a[i][1]
  b.append([a[i][0],length,prev_comm])
  for i in range(len(b)):
    print str(b[i][2])+'|id '+str(b[i][0]-b[i][1]+1)+'-'+str(b[i][0])
    cmd.select(b[i][2],str(b[i][2])+'|id '+str(b[i][0]-b[i][1]+1)+'-'+str(b[i][0]))
  #for i in range(len(a)):
  #  cmd.select(a[i][1], str(a[i][1])+'|id '+str(a[i][0]+1))
  cmd.space('cmyk')
  #randomseed = random.randrange(1,100000)
  #random.seed(randomseed)
  #print str(randomseed)
  random.seed(10874)
  for i in range(comm_id+1):
    cmd.set_color("ccc"+str(i),hsv_to_rgb((((3.0*i)/13)+(1/13)*random.random())%1,0.5+0.5*random.random(),0.6+((2*i%5)*0.1)))
    #cmd.set_color("ccc" + str(i), [1-pow(random.random(),1.3),1-pow(random.random(),1.3),1-pow(random.random(),1.3)] )
    cmd.color("ccc" + str(i),i)
    #cmd.color(""+str(((i+1)*674)%997), i)
    #cmd.color("s"+format((random.randrange(0,127) + 5*128*i)%1000, "03d"),i)
  cmd.disable(str(comm_id))
cmd.extend('community_highlight',community_highlight)

import os
from string import find

def community_pngs(directory,filemask):
  directory_old = os.getcwd()
  os.chdir(directory)
  files = os.listdir('.')
  files.sort()
  for i in range(len(files)):
    if (find(files[i], filemask)!=-1):
      community_highlight(files[i])
      cmd.rotate('y',6)
      cmd.ray(800,600)
      cmd.png( 'png/'+files[i]+'.png')
  os.chdir(directory_old)
cmd.extend('community_pngs',community_pngs)

def file_select(filename):
  f = open(filename, 'r')
  cmd.select('selection', 'none')
  for line in f :
    cmd.select('selection', 'selection|id '+ line)
cmd.extend('file_select',file_select)
