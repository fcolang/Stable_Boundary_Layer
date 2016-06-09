import numpy as np
import scipy.io as sio
from matplotlib.ticker import ScalarFormatter, MultipleLocator
import matplotlib.mlab as mlab
import scipy as sp
from scipy.interpolate import UnivariateSpline
import scipy.interpolate as si
from scipy.interpolate import interp1d
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

base_dir = os.path.expanduser('~')
path_data=base_dir+'/Dropbox/Monash_Uni/SBL/Dir_shifts/Dumosa/00Graficos_Paper/'
#*****************************************************************************\
colorh='#2929a3'
colorb='blue'
figsize3=(7.5, 5)
figsize6=(11, 5)
figsize2=(7.5, 7.5)
fsize0=12
fsize1=14
fsize2=16
#*****************************************************************************\
#Figure 2
#*****************************************************************************\
matb1= sio.loadmat(path_data+'mat/fig2_20160418.mat')

ncases=matb1['A2'][:]
n, bin_edges =np.histogram(ncases)

fig, ax1 = plt.subplots(figsize=(7, 5))
# Set the bar width
bar_width = 0.7

# positions of the left bar-boundaries
bar_l = [i+1 for i in range(len(n))]

# positions of the x-axis ticks (center of the bars as bar labels)
tick_pos = [i+(bar_width/2) for i in bar_l]

# Create the plot bars in x position
ax1.bar(bar_l,
        # using the pre_score data
        n,
        # set the width
        width=bar_width,
        # with alpha 0.5
        alpha=0.5,
        # with color
        color=colorh)

# set the x ticks with names
plt.xticks(tick_pos, range(1,len(n)+1),fontsize = 11)
plt.yticks(np.arange(0,60,5),fontsize = 11)
# add a grid
plt.grid()
# set axes labels and title
ax1.set_xlabel('Number of events')
ax1.set_ylabel('Ocurrences')
#plt.title('Mean Scores For Each Test')
# Set a buffer around the edge
plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])
plt.ylim([0,55])
fig.savefig('fig/fig2.eps', format='eps', dpi=1200)

# #****************************************************************************\
# # Figure 3
# #*****************************************************************************\
# matb2= sio.loadmat(path_data+'mat/fig3_20160418.mat')
# xdsh=matb2['shift_e'][:]
# xdt=matb2['xDTs'][:]
# xdwsp=matb2['xDwsp'][:]
# xdu=matb2['xDu'][:]
# xdv=matb2['xDv'][:]
# xdw=matb2['xDw'][:]

# #*****************************************************************************\
# #Wind Direction
# x=xdsh[0,:]
# bx=np.arange(-180,190,30)
# by=np.arange(0,25,5)
# x1=0
# x2=20
# y1=-180
# y2=180
# bins = np.arange(-180, 190, 10)

# fig, ax=plt.subplots(figsize=(7.5, 5))
# ax.hist(x, bins=bins, alpha=0.5, color=colorh)
# plt.xticks(bx,fontsize = 12)
# plt.yticks(by,fontsize = 12)
# plt.grid()
# plt.ylim([x1,x2])
# plt.xlim([y1,y2])
# plt.xlabel('Wind dir. shift (deg.)',fontsize = 12)
# plt.ylabel('Ocurrences',fontsize = 12)

# fig.savefig('fig/fig3_a.eps', format='eps', dpi=1200)
# #*****************************************************************************\
# #Wind Speed
# x=xdwsp[0,:]
# x1=0
# x2=25
# y1=-1.5
# y2=3.5
# by=np.arange(x1,x2+5,5)
# bx=np.arange(y1,y2+5,0.5)
# bins = np.arange(y1, y2+5, 0.1)

# fig, ax=plt.subplots(figsize=(7.5, 5))
# ax.hist(x, bins=bins, alpha=0.5, color=colorh)
# plt.xticks(bx,fontsize = 12)
# plt.yticks(by,fontsize = 12)
# plt.grid()
# plt.ylim([x1,x2])
# plt.xlim([y1,y2])
# plt.xlabel('Wind speed. shift (ms$^{-1}$)',fontsize = 12)
# plt.ylabel('Ocurrences',fontsize = 12)

# fig.savefig('fig/fig3_b.eps', format='eps', dpi=1200)
# #*****************************************************************************\
# #v-component
# x=xdv[0,:]
# x1=0
# x2=30
# y1=-3
# y2=3
# by=np.arange(x1,x2+5,5)
# bx=np.arange(y1,y2+5,0.5)
# bins = np.arange(y1, y2+5, 0.25)

# fig, ax=plt.subplots(figsize=figsize3)
# ax.hist(x, bins=bins, alpha=0.5, color=colorh)
# plt.xticks(bx,fontsize = 12)
# plt.yticks(by,fontsize = 12)
# plt.grid()
# plt.ylim([x1,x2])
# plt.xlim([y1,y2])
# plt.xlabel('v-comp. shift (ms$^{-1}$)',fontsize = 12)
# plt.ylabel('Ocurrences',fontsize = 12)

# fig.savefig('fig/fig3_c.eps', format='eps', dpi=1200)
# #*****************************************************************************\
# #u-component
# x=xdu[0,:]
# x1=0
# x2=30
# y1=-3
# y2=4.5
# by=np.arange(x1,x2+5,5)
# bx=np.arange(y1,y2+5,0.5)
# bins = np.arange(y1, y2+5, 0.25)

# fig, ax=plt.subplots(figsize=figsize3)
# ax.hist(x, bins=bins, alpha=0.5, color=colorh)
# plt.xticks(bx,fontsize = 12)
# plt.yticks(by,fontsize = 12)
# plt.grid()
# plt.ylim([x1,x2])
# plt.xlim([y1,y2])
# plt.xlabel('u-comp. shift (ms$^{-1}$)',fontsize = 12)
# plt.ylabel('Ocurrences',fontsize = 12)

# fig.savefig('fig/fig3_e.eps', format='eps', dpi=1200)
# #*****************************************************************************\
# #w-component
# x=xdw[0,:]
# x1=0
# x2=30
# y1=-0.2
# y2=0.2
# by=np.arange(x1,x2+5,5)
# bx=np.arange(y1,y2+5,0.04)
# bins = np.arange(y1, y2+5, 0.01)

# fig, ax=plt.subplots(figsize=figsize3)
# ax.hist(x, bins=bins, alpha=0.5, color=colorh)
# plt.xticks(bx,fontsize = 12)
# plt.yticks(by,fontsize = 12)
# plt.grid()
# plt.ylim([x1,x2])
# plt.xlim([y1,y2])
# plt.xlabel('w-comp. shift (ms$^{-1}$)',fontsize = 12)
# plt.ylabel('Ocurrences',fontsize = 12)

# fig.savefig('fig/fig3_d.eps', format='eps', dpi=1200)
# #*****************************************************************************\
# #Temperature
# x=xdt[0,:]
# x1=0
# x2=20
# y1=-3.6
# y2=1.6
# by=np.arange(x1,x2+5,5)
# bx=np.arange(y1,y2+5,0.4)
# bins = np.arange(y1, y2+5, 0.1)

# fig, ax=plt.subplots(figsize=figsize3)
# ax.hist(x, bins=bins, alpha=0.5, color=colorh)
# plt.xticks(bx,fontsize = 12)
# plt.yticks(by,fontsize = 12)
# plt.grid()
# plt.ylim([x1,x2])
# plt.xlim([y1,y2])
# plt.xlabel('Temperature shift (K)',fontsize = 12)
# plt.ylabel('Ocurrences',fontsize = 12)
# fig.savefig('fig/fig3_f.eps', format='eps', dpi=1200)
# #*****************************************************************************\
# #Figure 4
# #*****************************************************************************\
# matb3= sio.loadmat(path_data+'mat/fig4_20160418.mat')
# temp_p=matb3['Trmean_p'][:]
# temp_n=matb3['Trmean_n'][:]
# w_p=matb3['w_pt'][:]
# w_n=matb3['w_nt'][:]
# u_p=matb3['u_pt'][:]
# u_n=matb3['u_nt'][:]
# v_p=matb3['v_pt'][:]
# v_n=matb3['v_nt'][:]
# x=matb3['x'][:]
# x=x[0,:]
# #*****************************************************************************\
# #Temperature
# yp=temp_p[0,:]
# yn=temp_n[0,:]
# bx=np.arange(0,22,2)
# #bx=np.arange(y1,y2+5,0.4)

# fig, ax=plt.subplots(figsize=figsize3)
# ax.plot(x, yp, linewidth=2,marker='o',ms=5,ls='dotted',label='clockwise')
# ax.plot(x, yn, linewidth=2,marker='o',ms=5,ls='dotted',color='red',label='counter-clockwise')
# plt.xticks(bx,fontsize = 12)
# plt.ylabel('Temperature (K)',fontsize = 12)
# plt.xlabel('Time (min)',fontsize = 12)
# plt.legend(loc='lower left')
# plt.grid()

# fig.savefig('fig/fig4_a.eps', format='eps', dpi=1200)
# #*****************************************************************************\
# #w-component
# yp=w_p[0,:]
# yn=w_n[0,:]
# bx=np.arange(0,22,2)

# fig, ax=plt.subplots(figsize=figsize3)
# ax.plot(x, yp, linewidth=2,marker='o',ms=5,ls='dotted',label='clockwise')
# ax.plot(x, yn, linewidth=2,marker='o',ms=5,ls='dotted',color='red',label='counter-clockwise')
# plt.xticks(bx,fontsize = 12)
# plt.ylabel('w-comp. (ms$^{-1}$)',fontsize = 12)
# plt.xlabel('Time (min)',fontsize = 12)
# plt.grid()

# fig.savefig('fig/fig4_b.eps', format='eps', dpi=1200)
# #*****************************************************************************\
# #u-component
# yp=u_p[0,:]
# yn=u_n[0,:]
# bx=np.arange(0,22,2)


# fig, ax=plt.subplots(figsize=figsize3)
# ax.plot(x, yp, linewidth=2,marker='o',ms=5,ls='dotted',label='clockwise')
# ax.plot(x, yn, linewidth=2,marker='o',ms=5,ls='dotted',color='red',label='counter-clockwise')
# plt.xticks(bx,fontsize = 12)
# plt.ylabel('u-comp. (ms$^{-1}$)',fontsize = 12)
# plt.xlabel('Time (min)',fontsize = 12)
# plt.grid()

# fig.savefig('fig/fig4_c.eps', format='eps', dpi=1200)
# #*****************************************************************************\
# #v-component
# yp=v_p[0,:]
# yn=v_n[0,:]
# bx=np.arange(0,22,2)


# fig, ax=plt.subplots(figsize=figsize3)
# ax.plot(x, yp, linewidth=2,marker='o',ms=5,ls='dotted',label='clockwise')
# ax.plot(x, yn, linewidth=2,marker='o',ms=5,ls='dotted',color='red',label='counter-clockwise')
# plt.xticks(bx,fontsize = 12)
# plt.ylabel('v-comp. (ms$^{-1}$)',fontsize = 12)
# plt.xlabel('Time (min)',fontsize = 12)
# plt.grid()

# fig.savefig('fig/fig4_d.eps', format='eps', dpi=1200)

# #*****************************************************************************\
#Figure 2
#*****************************************************************************\
matb4= sio.loadmat(path_data+'mat/fig5_20160418.mat')
u=matb4['u_1s'][:]
u=u[0,:]
v=matb4['v_1s'][:]
v=v[0,:]
xdw=matb4['xdw'][:]
xdw=xdw[0,:]

u1=matb4['u1_1s'][:]
u1=u1[0,:]
v1=matb4['v1_1s'][:]
v1=v1[0,:]
xdw1=matb4['xdw1'][:]
xdw1=xdw1[0,:]

u2=matb4['u2_1s'][:]
u2=u2[0,:]
v2=matb4['v2_1s'][:]
v2=v2[0,:]
xdw2=matb4['xdw2'][:]
xdw2=xdw2[0,:]

u3=matb4['u3_1s'][:]
u3=u3[0,:]
v3=matb4['v3_1s'][:]
v3=v3[0,:]
xdw3=matb4['xdw3'][:]
xdw3=xdw3[0,:]
#*****************************************************************************\
row = 4
column = 2
fig, axes = plt.subplots(row, column, facecolor='w', figsize=(12,16))
ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7= axes.flat
xcomp=np.arange(1,len(u)+1,1)
x=np.arange(1,len(xdw)+1,1)
by=np.arange(-180,180+5,60)


ax0.plot(x,xdw, linewidth=2,marker='o',ms=5,ls='dotted')
ax0.set_xlabel('Time (min)',fontsize = fsize1)
ax0.set_ylabel('Wind dir. shift (deg)',fontsize = fsize1)
ax0.set_title('Main Tower',fontsize = fsize1,fontweight='bold')
ax0.tick_params(axis='both', which='major', labelsize=14)
ax0.set_ylim([-180,180])
ax0.set_yticks(by)
ax0.grid()


ax1.plot(xcomp,u, linewidth=1,ls='solid',label='u-comp.')
ax1.plot(xcomp,v, linewidth=1,ls='solid',color='red',label='v-comp.')
ax1.set_xlabel('Time (s)',fontsize = fsize1)
ax1.set_ylabel('(ms$^{-1}$)',fontsize = fsize1)
ax1.set_title('Main Tower',fontsize = fsize1,fontweight='bold')
ax1.tick_params(axis='both', which='major', labelsize=fsize1)
ax1.legend(loc='upper left')
ax1.set_ylim([-3,4])
ax1.grid()


#******************************************************************************
ax2.plot(x,xdw1, linewidth=2,marker='o',ms=5,ls='dotted')
ax2.set_xlabel('Time (min)',fontsize = fsize1)
ax2.set_ylabel('Wind dir. shift (deg)',fontsize = fsize1)
ax2.set_title('Station 1',fontsize = fsize1,fontweight='bold')
ax2.tick_params(axis='both', which='major', labelsize=fsize1)
ax2.set_ylim([-180,180])
ax2.set_yticks(by)
ax2.grid()

ax3.plot(xcomp,u1, linewidth=1,ls='solid',label='u-comp.')
ax3.plot(xcomp,v1, linewidth=1,ls='solid',color='red',label='v-comp.')
ax3.set_xlabel('Time (s)',fontsize = fsize1)
ax3.set_ylabel('(ms$^{-1}$)',fontsize = fsize1)
ax3.set_title('Station 1',fontsize = fsize1,fontweight='bold')
ax3.tick_params(axis='both', which='major', labelsize=fsize1)
ax3.set_ylim([-3,4])
ax3.grid()
#******************************************************************************
ax4.plot(x,xdw2, linewidth=2,marker='o',ms=5,ls='dotted')
ax4.set_xlabel('Time (min)',fontsize = fsize1)
ax4.set_ylabel('Wind dir. shift (deg)',fontsize = fsize1)
ax4.set_title('Station 2',fontsize = fsize1,fontweight='bold')
ax4.tick_params(axis='both', which='major', labelsize=fsize1)
ax4.set_ylim([-180,180])
ax4.set_yticks(by)
ax4.grid()

ax5.plot(xcomp,u2, linewidth=1,ls='solid',label='u-comp.')
ax5.plot(xcomp,v2, linewidth=1,ls='solid',color='red',label='v-comp.')
ax5.set_xlabel('Time (s)',fontsize = fsize1)
ax5.set_ylabel('(ms$^{-1}$)',fontsize = fsize1)
ax5.set_title('Station 2',fontsize = fsize1,fontweight='bold')
ax5.tick_params(axis='both', which='major', labelsize=fsize1)
ax5.set_ylim([-3,4])
ax5.grid()
#******************************************************************************
ax6.plot(x,xdw3, linewidth=2,marker='o',ms=5,ls='dotted')
ax6.set_xlabel('Time (min)',fontsize = fsize1)
ax6.set_ylabel('Wind dir. shift (deg)',fontsize = fsize1)
ax6.set_title('Station 3',fontsize = fsize1,fontweight='bold')
ax6.tick_params(axis='both', which='major', labelsize=fsize1)
ax6.set_ylim([-180,180])
ax6.set_yticks(by)
ax6.grid()

ax7.plot(xcomp,u3, linewidth=1,ls='solid',label='u-comp.')
ax7.plot(xcomp,v3, linewidth=1,ls='solid',color='red',label='v-comp.')
ax7.set_xlabel('Time (s)',fontsize = fsize1)
ax7.set_ylabel('(ms$^{-1}$)',fontsize = fsize1)
ax7.set_title('Station 3',fontsize = fsize1,fontweight='bold')
ax7.tick_params(axis='both', which='major', labelsize=fsize1)
ax7.set_ylim([-3,4])
ax7.grid()
#******************************************************************************
plt.tight_layout(pad=1, w_pad=1, h_pad=1.0)
#plt.show()
fig.savefig('fig/fig2.eps', format='eps', dpi=1200)
#*****************************************************************************\
# #Figure 6
# #*****************************************************************************\
matb5= sio.loadmat(path_data+'mat/fig6_propag_20160418.mat')
AA=matb5['AA'][:]
AV=matb5['AV'][:]
# #*****************************************************************************\
# x=np.arange(1,len(AA)+1,1)
# by=np.arange(0,390,60)
# bx=np.arange(0,len(AA)+5,5)
# fig, ax=plt.subplots(facecolor='w',figsize=figsize6)

# ax.plot(x, AA[:,0],marker='o',ms=9,ls='None',label='St1-St2-St3',markeredgecolor='blue', markerfacecolor='None',markeredgewidth=1)
# ax.plot(x, AA[:,1],marker='x',ms=7,ls='None',label='MSt-St1-St2',markeredgecolor='red', markerfacecolor='None',markeredgewidth=1)
# ax.plot(x, AA[:,2],marker='d',ms=9,ls='None',label='MSt-St1-St3',markeredgecolor='green', markerfacecolor='None',markeredgewidth=1)
# ax.plot(x, AA[:,3],marker='*',ms=9,ls='None',label='MSt-St2-St3',markeredgecolor='orange', markerfacecolor='None',markeredgewidth=1)


# # ax.plot(x, AA[:,0],marker='o',ms=9,ls='None',label='St1-St2-St3',color='blue')
# # ax.plot(x, AA[:,1],marker='o',ms=7,ls='None',label='MSt-St1-St2',color='red')
# # ax.plot(x, AA[:,2],marker='o',ms=9,ls='None',label='MSt-St1-St3',color='green')
# # ax.plot(x, AA[:,3],marker='o',ms=9,ls='None',label='MSt-St2-St3',color='orange')


# plt.ylabel('Direction of propagation (deg)',fontsize = fsize1)
# plt.xlabel('Case',fontsize = fsize1)
# plt.yticks(by,fontsize = fsize0)
# plt.xticks(bx,fontsize = fsize0)
# ax.set_ylim([0,360])
# ax.set_xlim([0,130])
# plt.grid()

# fig.savefig('fig/fig6_a.eps', format='eps', dpi=1200)

# #*****************************************************************************
# by=np.arange(0,22,2)
# fig, ax=plt.subplots(facecolor='w',figsize=figsize6)

# ax.plot(x, AV[:,0],marker='o',ms=9,ls='None',label='St1-St2-St3',markeredgecolor='blue', markerfacecolor='None',markeredgewidth=1)
# ax.plot(x, AV[:,1],marker='x',ms=7,ls='None',label='MSt-St1-St2',markeredgecolor='red', markerfacecolor='None',markeredgewidth=1)
# ax.plot(x, AV[:,2],marker='d',ms=9,ls='None',label='MSt-St1-St3',markeredgecolor='green', markerfacecolor='None',markeredgewidth=1)
# ax.plot(x, AV[:,3],marker='*',ms=9,ls='None',label='MSt-St2-St3',markeredgecolor='orange', markerfacecolor='None',markeredgewidth=1)

# plt.ylabel('Speed of propagation (m s$^{-1}$)',fontsize = fsize1)
# plt.xlabel('Case',fontsize = fsize1)
# plt.yticks(by,fontsize = fsize0)
# plt.xticks(bx,fontsize = fsize0)
# ax.set_ylim([0,20])
# ax.set_xlim([0,130])
# plt.legend(loc='upper right',numpoints=1)
# plt.grid()
# fig.savefig('fig/fig6_b.eps', format='eps', dpi=1200)

# #*****************************************************************************\
# #Figure 7
# #*****************************************************************************\
MProp=matb5['MProp'][:]
n_cases=np.count_nonzero(~np.isnan(MProp[:,2]))
matb0= sio.loadmat(path_data+'mat/windprop_20160421.mat')
binDP=matb0['binDP'][:]
binDP=binDP[0,:]
binSP=matb0['binSP'][:]
binSP=binSP[0,:]


#*****************************************************************************\
# bins = np.arange(0, 370, 30)
# nperc=binDP/float(np.sum(binDP))*100

# fig, ax1 = plt.subplots(figsize=figsize3)
# bar_width = 0.7

# bar_l = [i+1 for i in range(len(nperc))]
# tick_pos = [i+(bar_width/2) for i in bar_l]

# ax1.bar(bar_l,
#         # using the pre_score data
#         nperc,
#         # set the width
#         width=bar_width,
#         # with alpha 0.5
#         alpha=0.5,
#         # with color
#         color=colorh)

# plt.xticks(tick_pos, bins,fontsize = fsize1)
# plt.yticks(np.arange(0,26,2),fontsize = fsize1)
# plt.grid()
# ax1.set_xlabel('Direction of propagation (deg)',fontsize = fsize1)
# ax1.set_ylabel('Percentage of cases',fontsize = fsize1)
# plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])
# plt.ylim([0,24])
# fig.savefig('fig/fig7_a.eps', format='eps', dpi=1200)

#*****************************************************************************\
bins = np.arange(0, 21, 1)
bx=np.arange(0, 21,1)
by=np.arange(0,40,5)


#n, bins2, _ =ax.hist(binSP, bins=bins, alpha=0.5, color=colorh, weights=np.zeros_like(binSP) + 100./float(np.sum(binSP)))

nperc=binSP/float(np.sum(binSP))*100

fig, ax = plt.subplots(figsize=figsize3)
bar_width = 1

bar_l = [i+1 for i in range(len(nperc))]
tick_pos = [i+(bar_width/2) for i in bar_l]

ax.bar(bar_l,
        # using the pre_score data
        nperc,
        # set the width
        width=bar_width,
        # with alpha 0.5
        alpha=0.5,
        # with color
        color=colorh)

plt.xticks(tick_pos, bins,fontsize = fsize1)
plt.yticks(by,fontsize = fsize1)
plt.grid()
ax.set_xlabel('Speed of propagation (ms$^{-1}$)',fontsize = fsize1)
ax.set_ylabel('Percentage of cases',fontsize = fsize1)
plt.ylim([0,35])
fig.tight_layout()
fig.savefig('fig/fig7_b.eps', format='eps', dpi=1200)

# #*****************************************************************************\
# # Figure 9 y 10
# #*****************************************************************************\
binDW=matb0['binDW'][:]
binDW=binDW[0,:]
binSPW=matb0['binSPW'][:]
binSPW=binSPW[0,:]
#*****************************************************************************\
# bins = np.arange(0, 370, 30)
# nperc1=binDW/float(np.sum(binDW))*100

# fig, ax1 = plt.subplots(figsize=figsize3)
# bar_width = 0.7

# bar_l = [i+1 for i in range(len(nperc1))]
# tick_pos = [i+(bar_width/2) for i in bar_l]

# ax1.bar(bar_l,
#         # using the pre_score data
#         nperc1,
#         # set the width
#         width=bar_width,
#         # with alpha 0.5
#         alpha=0.5,
#         # with color
#         color=colorh)

# plt.xticks(tick_pos, bins,fontsize = fsize1)
# plt.yticks(np.arange(0,26,2),fontsize = fsize1)
# plt.grid()
# ax1.set_xlabel('Wind direction(deg)',fontsize = fsize1)
# ax1.set_ylabel('Percentage of cases',fontsize = fsize1)
# plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])
# plt.ylim([0,24])
# fig.savefig('fig/fig9.eps', format='eps', dpi=1200)

#*****************************************************************************\
bins = np.arange(0, 21, 1)
bx=np.arange(-5,21,1)
by=np.arange(0,50,5)

#n, bins2, _ =ax.hist(binSP, bins=bins, alpha=0.5, color=colorh, weights=np.zeros_like(binSP) + 100./float(np.sum(binSP)))

nperc1=binSPW/float(np.sum(binSPW))*100

fig, ax = plt.subplots(figsize=figsize3)
bar_width = 1

bar_l = [i for i in range(len(nperc1))]
tick_pos = [i+(bar_width/2) for i in bar_l]

ax.bar(np.arange(-1,18,1),
        # using the pre_score data
        nperc1,
        # set the width
        width=bar_width,
        # with alpha 0.5
        alpha=0.5,
        # with color
        color=colorh)

#plt.xticks(tick_pos, bins,fontsize = fsize1)
plt.xticks(bx,fontsize = fsize0)
plt.yticks(by,fontsize = fsize1)
plt.grid()
ax.set_xlabel('Difference prop. speed-wind speed (ms$^{-1}$)',fontsize = fsize1)
ax.set_ylabel('Percentage of cases',fontsize = fsize1)
plt.ylim([0,40])
fig.tight_layout()
fig.savefig('fig/fig10.eps', format='eps', dpi=1200)

#*****************************************************************************\
# Figure 12
#*****************************************************************************
matb6= sio.loadmat(path_data+'mat/propag_all_smo_20160418.mat')
prop123=matb6['propag123'][:]
wind=matb6['wind_mean'][:]

#*****************************************************************************
prop123[:,2]

indsr = np.nonzero(prop123[:,2]>30)
prop123[indsr,:]=np.nan
wind[indsr,:]=np.nan

#*****************************************************************************
y=prop123[:,1]
x=wind[:,1]
fig, ax = plt.subplots(figsize=figsize2)
ax.scatter(x, y,alpha=0.5)

plt.yticks(np.arange(0,420,60),fontsize = fsize1)
plt.xticks(np.arange(0,420,60),fontsize = fsize1)
ax.set_xlabel('Wind direction (deg)',fontsize = fsize1)
ax.set_ylabel('Prop. direction (deg)',fontsize = fsize1)
plt.ylim([0,360])
plt.xlim([0,360])
plt.grid()
fig.savefig('fig/fig12_a.eps', format='eps', dpi=1200)
#*****************************************************************************
y=prop123[:,2]
x=wind[:,2]
fig, ax = plt.subplots(figsize=figsize2)
ax.scatter(x, y,alpha=0.5)

plt.yticks(np.arange(0,40,5),fontsize = fsize1)
plt.xticks(np.arange(0,40,0.5),fontsize = fsize1)
ax.set_xlabel('Wind speed (ms$^{-1}$)',fontsize = fsize1)
ax.set_ylabel('Prop. speed (ms$^{-1}$)',fontsize = fsize1)
plt.ylim([0,30])
plt.xlim([0,4])
plt.grid()
fig.savefig('fig/fig12_b.eps', format='eps', dpi=1200)
#plt.show()
