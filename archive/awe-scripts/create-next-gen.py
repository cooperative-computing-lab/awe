#!/usr/bin/env python

import glob
import os
import random
import shutil
import sys

oldgendir = sys.argv[1]
newgendir = sys.argv[2]
run = int(sys.argv[3])
clone = int(sys.argv[4])
gen = int(sys.argv[5])
newgen = int(sys.argv[6])

me = sys.argv[0]

real_newgendir = reduce(lambda p, f: f(p), (newgendir, os.path.expanduser, os.path.expandvars, os.path.realpath, os.path.abspath))

if not os.path.exists(real_newgendir):
    os.makedirs(real_newgendir)

inp_files = glob.glob(os.path.join(oldgendir, "*.inp"))
psf_files = glob.glob(os.path.join(oldgendir, "*.psf"))

for inp_file in inp_files:
    shutil.copy(inp_file, real_newgendir)

for psf_file in psf_files:
    shutil.copy(psf_file, real_newgendir)

final_pdb = os.path.join(oldgendir, "final.pdb")
start_pdb = os.path.join(real_newgendir, "start.pdb")
shutil.copy(final_pdb, start_pdb)

conf_files = glob.glob(os.path.join(oldgendir, "*.conf"))

for conf_file in conf_files:
    basename = os.path.basename(conf_file)
    targetname = os.path.join(real_newgendir, basename)
    in_fl = open(conf_file)
    out_fl = open(targetname, "w")
    rand = random.randint(0, 100000)
    for ln in in_fl:
        seed_s = ln.find("seed")
        if seed_s != -1:
            seed_e = seed_s + 4 + 2
            ln = ln[:seed_e] + str(rand) + "\n"
        out_fl.write(ln)
    in_fl.close()
    out_fl.close()

os.chdir(oldgendir)
