import numpy as np
import matplotlib.pyplot as plt
import cphmdanalysis as cphmd
from cphmdanalysis import compute_pkas as pka
from cphmdanalysis import HH_fitting as pka_fit
from cphmdanalysis import plotting_lambda_data as plotting
import glob as g

# This script is used to generate figures from CpHMD simulations of YiiP

# Suppl. Fig. 7
def gather_stages(n):
    stages=[]
    for i in range(n):
        files=g.glob('radius14/stage1.{}/min.phrex.log*'.format(i))
        keys = [int(a.split('_')[1]) for a in files]
        combine=zip(keys, files)
        sorted_zipped=sorted(combine)
        sorted_files=[b for a, b in sorted_zipped ]
        stages.append(sorted_files)
    return stages

stages=gather_stages(260)

logs = cphmd.log_analysis_charmm(stages)

ph_list = [1.5, 1.75, 2, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 8.75, 9, 9.5, 10, 10.25, 10.5, 11, 11.5]

logs.plot_replica_walk(xlabel='Frame', save_fig=True, output='repwalk.pdf')

# ----------------------------------------------------------------------------------------------------------------
# Suppl. Fig. 8
l_files=[]
for i, ph in enumerate(ph_list):
    l_files.append('radius14/stage1.0/step0.ph1.5.lamb_{}'.format(i))

ph_objects = [cphmd.lambda_data(x) for x in l_files]

for n in range(1,260):
    for i, ph in enumerate(ph_list):
        ph_objects[i].add_l_file(file_path='radius14/stage1.{}/step{}.ph1.5.lamb_{}'.format(n,n,i))

print('Number of pHs: {}'.format(len(ph_objects)))
print('Number of titra site: {}'.format(ph_objects[0].n_ititr))
print('Number of lambda vals: {}'.format(len(ph_objects[0].lambda_and_x_vals[0])))

important_resids = [47-6, 51-6, 70-6, 73-6, 77-6, 155-6, 159-6, 
                    47+276, 51+276, 70+276, 73+276, 77+276, 155+276, 159+276, ]
resids = ph_objects[0].find_residues(important_resids)
print(resids)  

for x in range(len(ph_objects)):
    print('--- pH {} ----------------------------'.format(ph_list[x]))
    ph_objects[x].compute_all_s_values(output=True)

for x in range(len(ph_list)):
    ph_objects[x].compute_all_running_s()

pkas = pka(ph_list, ph_objects)

titles = ['D47A', 'D51A', 'D70A', 'H73A', 'H77A', 'H155A', 'D159A', 
          'D47B', 'D51B', 'D70B', 'H73B', 'H77B', 'H155B', 'D159B']
plotting.plot_running_s(ph_objects, ph_list, resids, titles, xlabel='Time [ns]', steps_to_time_conversion=(1/1000), save_fig=True, output='s_conv.pdf')


# ----------------------------------------------------------------------------------------------------------------
# Suppl. Fig. 8
for n, pka in enumerate(pkas):
    if ph_objects[0].ires[n] in important_resids:
        print('Resid: {0:3} pKa: {1:>3.1f} Hill: {2:>3.2f}'.format(ph_objects[0].ires[n], pka[0], pka[1]))

axes = plotting.plot_titration_curves(ph_objects, ph_list, resids, titles, xrange=[1.5,11.5], save_fig=False)

popt, pcov = curve_fit(couple, np.array(ph_list), np.array(s_values[0])+np.array(s_values[1]), p0=[8, 8], maxfev=100000)

axes[14].plot(ph_list, np.array(s_values[0])+np.array(s_values[1]), 'o', color='Black')
# Plot fit 
local_xs = np.linspace(0, 14, 100)
axes[14].plot(local_xs, couple(local_xs, np.max(popt), np.min(popt)), color='Black')
axes[14].axvline(np.max(popt), color='red')
axes[14].axvline(np.min(popt), color='red')
# Plot Details 
axes[14].set_title('H73AH77A')
axes[14].set_xlim(1.5, 11.5)
axes[14].set_xlabel('pH', fontsize=15)
#axes[index].set_xticks([x for x in range(xrange[0], xrange[1]+1, 2)])
axes[14].set_ylim(0, 2)
axes[14].set_ylabel('Unprot. Fraction', fontsize=12)
axes[14].text(2, 1.4, u'pK$_{a1}$='+'{:1.2f}'.format(np.max(popt))+ '\npK$_{a2}$='+'{:1.2f}'.format(np.min(popt)))

popt, pcov = curve_fit(couple, np.array(ph_list), np.array(s_values[2])+np.array(s_values[3]), p0=[8, 8], maxfev=100000)

axes[15].plot(ph_list, np.array(s_values[2])+np.array(s_values[3]), 'o', color='Black')
# Plot fit 
local_xs = np.linspace(0, 14, 100)
axes[15].plot(local_xs, couple(local_xs, np.max(popt), np.min(popt)), color='Black')
axes[15].axvline(np.max(popt), color='red')
axes[15].axvline(np.min(popt), color='red')
# Plot Details 
axes[15].set_title('H73BH77B')
axes[15].set_xlim(1.5, 11.5)
axes[15].set_xlabel('pH', fontsize=15)
#axes[index].set_xticks([x for x in range(xrange[0], xrange[1]+1, 2)])
axes[15].set_ylim(0, 2)
axes[15].set_ylabel('Unprot. Fraction', fontsize=12)
axes[15].text(2, 1.4, u'pK$_{a1}$='+'{:1.2f}'.format(np.max(popt))+ '\npK$_{a2}$='+'{:1.2f}'.format(np.min(popt)))

plt.savefig('titration_curves_coupled.pdf', dpi=300, transparent=False)