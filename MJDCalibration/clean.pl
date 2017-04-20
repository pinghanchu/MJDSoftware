#!/usr/bin/perl
#This script reset everything.
system("make clean");
system("rm out.* err.* *.csh List/*_*/*.txt Plot/*_*/*.pdf *.root *~ reso_*.txt cov_*.txt calibration_*.txt");
