#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#==============================================================================
import numpy as np
import matplotlib.pyplot as plt
#==============================================================================
computed =\
np.loadtxt('saved_values_final/ray_3.4e+05_dt_1.0e-03_350s.gz')
iteration = computed[:,0]
time = computed[:,1]
T_1 = computed[:,2]
vx_1 = computed[:,3]
vy_1 = computed[:,4]
eps_12 =  computed[:,6] 
p_diff14 = computed[:,7]
p_diff51 = computed[:,8] 
p_diff35 = computed[:,9] 
vortex_metric = computed[:,10]  
velocity_metric = computed[:,11]
nu_left = computed[:,12]
delta_t = 1.0e-03
#==============================================================================
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
params = {'text.latex.preamble': [r'\usepackage{siunitx}', 
                                  r"\usepackage{amsmath}"]}
plt.rcParams.update({'font.size': 20})
plt.rcParams.update(params)
#==============================================================================
# plot T_1
fig1 = plt.figure()
ax = fig1.add_subplot(1,1,1)
ax.plot(time, T_1,  color="blue", 
        label =r'$T_1$')
ax.set_xlabel(r'$t$', fontsize = 22)
ax.set_ylabel(r'$T_1$', fontsize = 22)
ax.grid(True)
#ax.legend(loc='best', shadow=True, fontsize = 16)
fig1.savefig('plots/temperatur_point1.pdf', 
             bbox_inches='tight')
plt.close()
#==============================================================================
# plot v_x_1
fig2  = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(time, vx_1, color="blue", 
         label =r'$v_{x,1}$')
ax2.set_xlabel(r'$t$', fontsize = 22)
ax2.set_ylabel(r'$v_{x,1}$', fontsize = 22)
ax2.grid(True)
ax2.legend(loc='best', shadow=True, fontsize = 16)
fig2.savefig('plots/v_x_point1.pdf', 
             bbox_inches='tight')
plt.close()
#==============================================================================
# plot v_y_1
fig3  = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.plot(time, vy_1, color="blue", 
         label =r'$v_{y,1}$')
ax3.set_xlabel(r'$t$', fontsize = 22)
ax3.set_ylabel(r'$v_{y,1}$', fontsize = 22)
ax3.grid(True)
ax3.legend(loc='best', shadow=True, fontsize = 16)
fig3.savefig('plots/v_y_point1.pdf', 
             bbox_inches='tight')
plt.close()
#==============================================================================
# plot eps12
fig4  = plt.figure()
ax4 = fig4.add_subplot(111)
ax4.plot(time, eps_12, color="blue", 
         label =r'$\varepsilon_{12}$')
ax4.set_xlabel(r'$t$', fontsize = 22)
ax4.set_ylabel(r'$\varepsilon_{12}$', fontsize = 22)
ax4.grid(True)
#ax4.legend(loc='best', shadow=True, fontsize = 16)
fig4.savefig('plots/eps12.pdf', 
             bbox_inches='tight')
plt.close()
#==============================================================================
# plot p_diff14 
fig5  = plt.figure()
ax5 = fig5.add_subplot(111)
ax5.plot(time, p_diff14, color="blue", 
         label =r'$\Delta p_{14}$')
ax5.set_xlabel(r'$t$', fontsize = 22)
ax5.set_ylabel(r'$\Delta p_{14}$', fontsize = 22)
ax5.grid(True)
ax5.legend(loc='best', shadow=True, fontsize = 16)
fig5.savefig('plots/pdiff14.pdf', 
             bbox_inches='tight')
plt.close()
#==============================================================================
# plot p_diff51 
fig6  = plt.figure()
ax6 = fig6.add_subplot(111)
ax6.plot(time, p_diff51, color="blue", 
         label =r'$\Delta p_{51}$')
ax5.set_xlabel(r'$t$', fontsize = 22)
ax5.set_ylabel(r'$\Delta p_{51}$', fontsize = 22)
ax5.grid(True)
ax5.legend(loc='best', shadow=True, fontsize = 16)
fig5.savefig('plots/pdiff51.pdf', 
             bbox_inches='tight')
plt.close()
#==============================================================================
# plot p_diff35 
fig7 = plt.figure()
ax7 = fig7.add_subplot(111)
ax7.plot(time, p_diff35, color="blue", 
         label =r'$\Delta p_{35}$')
ax7.set_xlabel(r'$t$', fontsize = 22)
ax7.set_ylabel(r'$\Delta p_{35}$', fontsize = 22)
ax7.grid(True)
ax7.legend(loc='best', shadow=True, fontsize = 16)
fig7.savefig('plots/pdiff35.pdf', 
             bbox_inches='tight')
plt.close()
#==============================================================================
# plot vortex_metric = computed[:,10]   
fig8 = plt.figure()
ax8 = fig8.add_subplot(111)
ax8.plot(time, vortex_metric, color="blue", 
         label =r'$\hat{\omega}$')
ax8.set_xlabel(r'$t$', fontsize = 22)
ax8.set_ylabel(r'$\hat{\omega}$', fontsize = 22)
ax8.grid(True)
ax8.legend(loc='best', shadow=True, fontsize = 16)
fig8.savefig('plots/omega_metric.pdf', 
             bbox_inches='tight')
plt.close()
#==============================================================================
# plot velocity_metric = computed[:,11]   
fig9 = plt.figure()
ax9 = fig9.add_subplot(111)
ax9.plot(time, velocity_metric, color="blue", 
         label =r'$\hat{v}$')
ax9.set_xlabel(r'$t$', fontsize = 22)
ax9.set_ylabel(r'$\hat{v}$', fontsize = 22)
ax9.grid(True)
ax9.legend(loc='best', shadow=True, fontsize = 16)
fig9.savefig('plots/velocity_metric.pdf', 
             bbox_inches='tight')
plt.close()
#==============================================================================
# plot nu_left = computed[:,12]   
fig10 = plt.figure()
ax10 = fig10.add_subplot(111)
ax10.plot(time, nu_left, color="blue", 
         label =r'$Nu_{\mathrm{L.}}$')
ax10.set_xlabel(r'$t$', fontsize = 22)
ax10.set_ylabel(r'$Nu_{\mathrm{L.}}$', fontsize = 22)
ax10.grid(True)
ax10.legend(loc='best', shadow=True, fontsize = 16)
fig10.savefig('plots/nusselt_left.pdf', 
             bbox_inches='tight')
plt.close()
#==============================================================================
plt.close("all")
#==============================================================================
# calculate average and oscillatory components
"""
In this implementation the respective list (e.g. T_1, vx_1, ...) has been 
plotted in order to get the index of the maximum value of a
given list by "np.argmax()". This needs to be done 
for the first and the second maximum of one or more oscillations.

Better Alternative: FFT-Function
"""
# time period tau of T_1
T_1_max1 = T_1[330000 + np.argmax(T_1[330000:332000])]
T_1_max2 = T_1[347000 + np.argmax(T_1[347000:348000])]
tau_T_start = time[330000 + np.argmax(T_1[330000:332000])]
tau_T_ende = time[347000 + np.argmax(T_1[347000:348000])]
tau_T = tau_T_ende - tau_T_start
# average value of T_1
T_avg = np.trapz(
        T_1[330000 + np.argmax(T_1[330000:332000]):
            347000 + np.argmax(T_1[347000:348000])], 
                           dx=delta_t) / tau_T
# oscillatory component
T_osc  = T_1_max1 - T_avg
#==============================================================================
# time period tau of vx_1
vx_1_max1 = vx_1[331000 + np.argmax(vx_1[331000:332000])]
vx_1_max2 = vx_1[348000 + np.argmax(vx_1[348000:349000])]
tau_vx_start = time[331000 + np.argmax(vx_1[331000:332000])]
tau_vx_end = time[348000 + np.argmax(vx_1[348000:349000])]
tau_vx = tau_vx_end - tau_vx_start
# average value of vx_1
vx_avg = np.trapz(
        vx_1[331000 + np.argmax(vx_1[331000:332000]):
            348000 + np.argmax(vx_1[348000:349000])], 
                           dx=delta_t) / tau_vx
# oscillatory component
vx_osc = vx_1_max1 - vx_avg 
#==============================================================================
# time period tau of vy_1  (Attention: Does not appear in document or paper)
vy_1_max1 = vy_1[35100 + np.argmax(vy_1[35100:35200])]
vy_1_max2 = vy_1[35500 + np.argmax(vy_1[35500:35600])]
tau_vy_start = time[35100 + np.argmax(vy_1[35100:35200])]
tau_vy_end = time[35500 + np.argmax(vy_1[35500:35600])]
tau_vy = tau_vy_end - tau_vy_start
# average value of vy
vy_avg = np.trapz(
        vy_1[35100 + np.argmax(vy_1[35100:35200]):
            35500 + np.argmax(vy_1[35500:35600])], 
                           dx=delta_t) / tau_vy
# oscillatory component
vy_osc = vy_1_max1 - vy_avg 
#==============================================================================
# time period tau of nu(left) 
nu_l_max1 = nu_left[331000 + np.argmax(nu_left[331000:332000])]
nu_l_max2 = nu_left[348000 + np.argmax(nu_left[348000:349000])]
tau_nu_l_start = time[331000 + np.argmax(nu_left[331000:332000])]
tau_nu_l_end = time[348000 + np.argmax(nu_left[348000:349000])]
tau_nu_l = tau_nu_l_end - tau_nu_l_start
# average value of nu(left)
nu_left_avg = np.trapz(
        nu_left[331000 + np.argmax(nu_left[331000:332000]):
            348000 + np.argmax(nu_left[348000:349000])], 
                           dx=delta_t) / tau_nu_l
# oscillatory component
nu_left_osc = nu_l_max1 - nu_left_avg
#==============================================================================
# time period tau of Delta p14
p14_max1 = p_diff14[330000 + np.argmax(p_diff14[330000:331000])]
p14_max2 = p_diff14[347000 + np.argmax(p_diff14[347000:348000])]
tau_p14_start = time[330000 + np.argmax(p_diff14[330000:331000])]
tau_p14_end = time[347000 + np.argmax(p_diff14[347000:348000])]
tau_p14 = tau_p14_end - tau_p14_start
# average value of Delta p14
p14_avg = np.trapz(
        p_diff14[330000 + np.argmax(p_diff14[330000:331000]):
            347000 + np.argmax(p_diff14[347000:348000])], 
                           dx=delta_t) / tau_p14
# oscillatory component
p14_osc = p14_max1 - p14_avg 
#==============================================================================
omega_max1 = vortex_metric[330000 + np.argmax(vortex_metric[330000:331000])]
omega_max2 = vortex_metric[348000 + np.argmax(vortex_metric[348000:349000])]
tau_omega_start = time[330000 + np.argmax(vortex_metric[330000:331000])]
tau_omega_end = time[348000 + np.argmax(vortex_metric[348000:349000])]
tau_omega_metric = tau_omega_end - tau_omega_start
# average value of \hat{\omega}
omega_metric_avg = np.trapz(
        vortex_metric[330000 + np.argmax(vortex_metric[330000:331000]):
            348000 + np.argmax(vortex_metric[348000:349000])], 
                           dx=delta_t) / tau_omega_metric
# oscillatory component
omega_metric_osc = omega_max1 - omega_metric_avg
#==============================================================================
velo_metric_max1 = velocity_metric[29800 
                                   + np.argmax(velocity_metric[29800:29900])]
velo_metric_max2 = velocity_metric[30100 
                                   + np.argmax(velocity_metric[30100:30200])]
tau_velocity_metric_start = time[29800 
                                 + np.argmax(velocity_metric[29800:29900])]
tau_velocity_metric_end = time[30100 
                               + np.argmax(velocity_metric[30100:30200])]
tau_velocity_metric = tau_velocity_metric_end - tau_velocity_metric_start
# average value of \hat{\omega}
velocity_metric_avg = np.trapz(
        velocity_metric[29800 + np.argmax(velocity_metric[29800:29900]):
            30100 + np.argmax(velocity_metric[30100:30200])], 
                           dx=delta_t) / tau_velocity_metric
# oscillatory component
velocity_metric_osc = velo_metric_max1 - velocity_metric_avg