#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
#==============================================================================
                                Postprocessor
                Reference values from  DFG Benchmark FeatFlow 
-------------------------------------------------------------------------------
                             Flow around cylinder 
                     benchmark 2D, time-periodic case Re=100
                             implicit Euler scheme
#==============================================================================
"""
import numpy as np
import matplotlib.pyplot as plt
#==============================================================================
# load reference values
refValues = np.loadtxt("bdforces_lv6_bench2", usecols=(1,3,4))
refTime = refValues[:,0]    # reference time
refDrag = refValues[:,1]    # reference drag coefficents
refLift = refValues[:,2]    # reference lift coefficients
# load computed values
computed = np.loadtxt('coeff_comp_8s_re100_implicit_euler_final.gz', 
                      usecols=(1,2,3)) 
#==============================================================================
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
params = {'text.latex.preamble': [r'\usepackage{siunitx}', 
                                  r"\usepackage{amsmath}"]}
plt.rcParams.update({'font.size': 20})
plt.rcParams.update(params)
#==============================================================================
# plot comparison of drag coefficents t \in [0, 8] seconds
fig1 = plt.figure()
ax = fig1.add_subplot(1,1,1)
ax.plot(refTime, refDrag, color="red", label =r'Referenz [FeatFlow]')
ax.plot(computed[:,0], computed[:,1], "--", color="black", 
        label =r'implizites \textsc{Euler}-Verfahren')
ax.set_title(r'Widerstandsbeiwerte f\"ur $t \in [0\,\si{s}, 8\,\si{s}]$', 
                                                 fontsize = 22)
ax.set_xlabel(r'$t$', fontsize = 22)
ax.set_ylabel(r'$c_{\mathrm{d}}(t)$', fontsize = 22)
ax.grid(True)
ax.legend(loc='best', shadow=True, fontsize = 16)
fig1.savefig('plots/drag_coefficients_euler_8000_final.pdf', 
             bbox_inches='tight')
#==============================================================================
# plot comparison of lift coefficents t \in [0, 8] seconds
fig2  = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(refTime, refLift,  color="red",label =r'Referenz [FeatFlow]')
ax2.plot(computed[:,0], computed[:,2], "--", color="black", 
         label =r'implizites \textsc{Euler}-Verfahren')
ax2.set_title(r'Auftriebsbeiwerte f\"ur $t \in [0\,\si{s}, 8\,\si{s}]$', 
                                                fontsize = 22)
ax2.set_xlabel(r'$t$', fontsize = 22)
ax2.set_ylabel(r'$c_{\mathrm{l}}(t)$', fontsize = 22)
ax2.grid(True)
ax2.legend(loc='best', shadow=True, fontsize = 16)
fig2.savefig('plots/lift_coefficients_euler_8000_final.pdf', 
             bbox_inches='tight')
#==============================================================================
# print maximum computed values (see bachelor thesis)
time = computed[:,0]
cd = computed[:,1]
cl = computed[:,2]
cd_max = cd[np.argmax(cd)]
time_cd = time[np.argmax(cd)]
cl_max = cl[np.argmax(cl)]
time_cl = time[np.argmax(cl)]
print "The maximum drag coefficent cd_max = ", cd_max, " is at t = ", \
        time_cd, "s"
print "The maximum lift coefficent cl_max = ", cl_max, " is at t = ", \
        time_cl, "s"
