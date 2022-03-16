## python3

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from numpy.polynomial.polynomial import polyfit

## set fontsize
import matplotlib
ft = 12
params = {'axes.labelsize': ft,'axes.titlesize':ft,'legend.fontsize': ft, 'xtick.labelsize': ft, 'ytick.labelsize': ft}
matplotlib.rcParams.update(params)


## read original InSAR time series 
time = []
los = []	
f = open('./ts_example.txt','r')	
while True:
	line = f.readline().strip()
	if line=='':
		#print ('break')
		break
	raw = line.split('	')
	time.append(float(raw[0]))
	los.append(float(raw[1]))		
f.close()

###### normalize time series
time_min = min(time)
time_max = max(time)
los_min = min(los)
los_max = max(los)

for i in range(len(time)):
	time[i] = (time[i]-time_min)/(time_max-time_min)
	los[i] = (los[i]-los_min)/(los_max - los_min)



## get the slope of the original time series using best-fit line
b, s = polyfit(time, los, 1)
#print (b)
#print (s)
s_degree = np.arctan(s)*180/np.pi
print ("original slope: ", s_degree)		

## rotation angle
range_step = 60. # # degree, a_open, the angle constraining the final slope range of the rotated time series
theta_min = int(s_degree-range_step) 
theta_max = int(s_degree+range_step)

## get the information of peak distance and difference
r_degree = []
peak_dis_res_ave = []
peak_height_ave = []	
f = open('./clusters_distances_height.txt','r')
line = f.readline().strip()	
while True:
	line = f.readline().strip()
	if line=='':
		#print ('break')
		break
	raw = line.split('	')
	if raw[1] != '-999':
		r_degree.append(float(raw[0]))
		peak_dis_res_ave.append(float(raw[1]))	
		peak_height_ave.append(float(raw[2]))	
f.close()	


if len(r_degree) == 0:
	r_best = -999

else:
	## normalize them
	peak_dis_res_ave_norm = [0.0]*len(peak_dis_res_ave)
	peak_height_ave_norm = [0.0]*len(peak_height_ave)

	peak_dis_res_ave_min = min(peak_dis_res_ave)
	peak_dis_res_ave_max = max(peak_dis_res_ave)
	peak_height_ave_min = min(peak_height_ave)
	peak_height_ave_max = max(peak_height_ave)
	for i in range(len(peak_dis_res_ave_norm)):
		peak_dis_res_ave_norm[i] = (peak_dis_res_ave[i]-peak_dis_res_ave_min)/(peak_dis_res_ave_max-peak_dis_res_ave_min)
		peak_height_ave_norm[i] = (peak_height_ave[i]-peak_height_ave_min)/(peak_height_ave_max - peak_height_ave_min)

	## find the point closest to point (0, 0)
	min_dis = 1.0
	r_best = -999
	peak_dis_best = peak_dis_res_ave[0]
	peak_height_best = peak_height_ave[0]
	peak_dis_norm_best = peak_dis_res_ave_norm[0]
	peak_height_norm_best = peak_height_ave_norm[0]
	for i in range(len(peak_dis_res_ave_norm)):
		tmp = (peak_dis_res_ave_norm[i]**2 + (peak_height_ave_norm[i]-1.0)**2)**(0.5)
		if tmp < min_dis:
			min_dis = tmp
			r_best = r_degree[i]
			peak_dis_best = peak_dis_res_ave[i]
			peak_height_best = peak_height_ave[i]
			peak_dis_norm_best = peak_dis_res_ave_norm[i] 
			peak_height_norm_best = peak_height_ave_norm[i]
				
	##### start of plotting
	fig = plt.figure(figsize=(10,4))
	gs = gridspec.GridSpec(1,2)
	gs.update(wspace=0.25, hspace=0.2)

	cm = plt.cm.get_cmap('jet')
	for i in range(2):
		ax = fig.add_subplot(gs[i])
		
		theta_min_tmp = -50 
		theta_max_tmp = 70
		
		if i ==0:
			sc = ax.scatter(peak_dis_res_ave, peak_height_ave, c=r_degree, vmin=theta_min_tmp, vmax=theta_max_tmp, cmap=cm)				
			ellipse = Ellipse((peak_dis_best, peak_height_best), 0.5, 0.04, color='red', ls='--', lw = 1, fill = False)			
			ax.set_xlabel('Averaged peak distance residual (mm)')
			ax.set_ylabel('Averaged peak height')
		
		 
		if i==1:
			sc = ax.scatter(peak_dis_res_ave_norm, peak_height_ave_norm, c=r_degree, vmin=theta_min_tmp, vmax=theta_max_tmp, cmap=cm) #vmin=theta_min, vmax=theta_max
			ellipse = Ellipse((peak_dis_norm_best, peak_height_norm_best), 0.05, 0.05, color='red', ls='--', lw = 1, fill = False)	
			
			cbar = plt.colorbar(sc, fraction=0.046, pad=0.04, orientation='vertical', ticks=np.arange(theta_min_tmp, theta_max_tmp, 20)) #ticks=np.arange(theta_min, theta_max+0.1, 15)
			cbar.ax.tick_params(labelsize=12)  # colorbar tick
			cbar.set_label('Rotation angle (degree)', fontsize=14, labelpad = 6) # colorbar label
			ax.set_xlabel('Averaged peak distance residual (normalized)')
			ax.set_ylabel('Averaged peak height (normalized)')
			
		ax.add_patch(ellipse)
		#ax.set_ylim([ymin, ymax])
		#ax.set_xlim([xmin, xmax])
		plt.gca().invert_yaxis()
		 

	fig.savefig('peak_distance_vs_height.png', dpi=300, bbox_inches='tight')



f_out = open('best_rotation_degree.txt', 'w')
f_out.write("%f"%(r_best))
print ("Best rotaion angle is: ", r_best)

f_out.close()


	



