## python3

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from scipy.signal import find_peaks
import scipy
from numpy.polynomial.polynomial import polyfit

## set fontsize
import matplotlib
ft = 12
params = {'axes.labelsize': ft,'axes.titlesize':ft,'legend.fontsize': ft, 'xtick.labelsize': ft, 'ytick.labelsize': ft}
matplotlib.rcParams.update(params)


###### set the threshold to avoid over-correction
thred_ambiguity_res = 0.2


## read original InSAR time series 
time = []
los = []	
time_new = []
los_new = []
f = open('./ts_example.txt','r')	
while True:
	line = f.readline().strip()
	if line=='':
		#print ('break')
		break
	raw = line.split('	')
	time.append(float(raw[0]))
	los.append(float(raw[1]))	
	time_new.append(float(raw[0]))
	los_new.append(float(raw[1]))	
f.close()


###### normalize time series
time_min = min(time)
time_max = max(time)
los_min = min(los)
los_max = max(los)


time_n = [0]*len(time)
los_n = [0]*len(los)
for i in range(len(time)):
	time_n[i] = (time[i]-time_min)/(time_max-time_min)
	los_n[i] = (los[i]-los_min)/(los_max - los_min)


theta_degree = np.loadtxt("best_rotation_degree.txt")


theta = theta_degree*np.pi/180.
time_r = [0.]*len(time)
los_r = [0.]*len(los)
for i in range(len(time_r)):
	time_r[i] = time_n[i]*np.cos(theta)+los_n[i]*np.sin(theta)
	los_r[i] = -1.*time_n[i]*np.sin(theta)+los_n[i]*np.cos(theta)


time_r_shift = min(time_r)
los_r_shift = min(los_r)
for i in range(len(time_r)):
	time_r[i] = time_r[i]-time_r_shift
	los_r[i] = los_r[i]-los_r_shift


##### start of plotting
ymin = -0.1
ymax = 1.2
	
xmin = -0.1
xmax = 1.2

ymax2 = ymax*(los_max-los_min) + los_min
ymin2 = ymin*(los_max-los_min) + los_min

xmax2 = xmax*(time_max-time_min) + time_min
xmin2 = xmin*(time_max-time_min) + time_min


fig = plt.figure(figsize=(11,10))
gs = gridspec.GridSpec(2,2)
gs.update(wspace=0.3, hspace=0.3)

ax = fig.add_subplot(gs[0])
plt.title('Original')
ax.scatter(time, los, color='black') 		
ax.set_ylim([ymin2, ymax2])
ax.set_xlim([xmin2, xmax2])
ax.set_xlabel('Year')
ax.set_ylabel('LOS displacement (mm)')

ax2 = fig.add_subplot(gs[1])
plt.title('Classified')
ax3 = fig.add_subplot(gs[2])
plt.title('Corrected')
color_list = ['cyan', 'green', 'blue', 'yellow','purple'] ## please modify as needed, e.g., add more colors


### calculate gaussian Kernel density of the rotated time series
iqr = scipy.stats.iqr(los_r)
mean = sum(los_r)/len(los_r)
tmp = 0.0
for j in range(len(los_r)):
	tmp = tmp + (los_r[j] - mean)**2

std = (tmp/(len(los_r)-1.0))**(0.5)

if std<iqr:
	kde_bw = 0.9*std*((len(los_r)**(-0.2)))
else:
	kde_bw = 0.9*iqr*((len(los_r)**(-0.2)))

kde = KernelDensity(kernel='gaussian', bandwidth=kde_bw).fit(np.array(los_r).reshape(-1,1))

X_plot = np.arange(-0.1, 1.5, 0.002)[:, np.newaxis]
log_dens = kde.score_samples(X_plot)
dens = np.exp(log_dens)


## detect crests and troughs
peaks = find_peaks(np.array(dens))[0]	
#print (peaks)	
troughs = find_peaks(-np.array(dens))[0]

if len(peaks)>1:
	### get the peaks that are due to phaes jumps
	mark = [0]*len(peaks)		

	## get the main peak (the peak with highest density)		 
	main_peak = peaks[0]
	for i in range(len(peaks)):
		#print ("peak:", peaks[i],  X_plot[int(peaks[i]),0])	
		if dens[int(peaks[i])]>dens[int(main_peak)]:
			main_peak = peaks[i]		
	#print ("main peak:",  X_plot[int(main_peak),0])	
	
	
	## check the distance between main peak and other peaks
	for j in range(len(peaks)):			
		### note that, convert peak distance to the real LOS displacement
		peak_dis = abs(X_plot[int(peaks[j]),0] - X_plot[int(main_peak),0])*(los_max - los_min)
		ambiguity_int = (peak_dis/27.75)
		ambiguity_res = (peak_dis/27.75)%1
		
		if ambiguity_res>=(1-thred_ambiguity_res) or (ambiguity_res<=thred_ambiguity_res and ambiguity_int>=1):
			mark[j] = 1			

	if sum(mark) == 0:		
		print("warning: the distance between peaks do not match half wavelenght!!")


	### correct other peaks to the main peak
	for i in range(len(peaks)):
		if mark[i] == 1: #and peaks[i]!=main_peak:
			## get which points below to this peak group
			if i == 0:
				for j in range(len(troughs)):
					if (X_plot[int(troughs[j]),0]-X_plot[int(peaks[0]),0])*(X_plot[int(troughs[j]),0]-X_plot[int(peaks[1]),0])<0:		
						cluster_start = X_plot[0, 0]
						cluster_end = X_plot[int(troughs[j]), 0]
			
			
			elif i == len(peaks)-1:
				for j in range(len(troughs)):
					if (X_plot[int(troughs[j]),0]-X_plot[int(peaks[i]),0])*(X_plot[int(troughs[j]),0]-X_plot[int(peaks[i-1]),0])<0:
						cluster_start = X_plot[int(troughs[j]), 0]
						cluster_end = X_plot[-1, 0]
			
			else:
				for j in range(len(troughs)):
					if (X_plot[int(troughs[j]),0]-X_plot[int(peaks[i-1]),0])*(X_plot[int(troughs[j]),0]-X_plot[int(peaks[i]),0])<0:
						cluster_start = X_plot[int(troughs[j]), 0]
							
				for j in range(len(troughs)):
					if (X_plot[int(troughs[j]),0]-X_plot[int(peaks[i]),0])*(X_plot[int(troughs[j]),0]-X_plot[int(peaks[i+1]),0])<0:
						cluster_end = X_plot[int(troughs[j]), 0]
							
			### 
			time_tmp = []
			los_tmp = []				
			for j in range(len(los_r)):
				if los_r[j]>=cluster_start and los_r[j]<cluster_end:
					time_tmp.append(time[j])
					los_tmp.append(los[j])
			
			if peaks[i]==main_peak:
				cl = 'red'
				z = 1
			else:
				cl = color_list[i]
				z = 100
			ax2.scatter(time_tmp, los_tmp, color=cl)
			### correct the group
			## find the ambiguity
			peak_dis = (X_plot[int(peaks[i]),0] - X_plot[int(main_peak),0])*(los_max - los_min)
			ambiguity = round(peak_dis/27.75)
			for j in range(len(los_tmp)):
				los_tmp[j] = los_tmp[j] - ambiguity*27.75
				
			ax3.scatter(time_tmp, los_tmp, color=cl, zorder=z)
			
				
			## update the time series
			for j in range(len(time_tmp)):
				for k in range(len(time_new)):
					if time_tmp[j] == time_new[k]:
						los_new[k] = los_tmp[j]
						
		
		if mark[i] == 0:
			## get which points below to this peak group
			if i == 0:
				for j in range(len(troughs)):
					if (X_plot[int(troughs[j]),0]-X_plot[int(peaks[0]),0])*(X_plot[int(troughs[j]),0]-X_plot[int(peaks[1]),0])<0:		
						cluster_start = X_plot[0, 0]
						cluster_end = X_plot[int(troughs[j]), 0]
			
			
			elif i == len(peaks)-1:
				for j in range(len(troughs)):
						if (X_plot[int(troughs[j]),0]-X_plot[int(peaks[i]),0])*(X_plot[int(troughs[j]),0]-X_plot[int(peaks[i-1]),0])<0:
							cluster_start = X_plot[int(troughs[j]), 0]
							cluster_end = X_plot[-1, 0]
			
			else:
				for j in range(len(troughs)):
						if (X_plot[int(troughs[j]),0]-X_plot[int(peaks[i-1]),0])*(X_plot[int(troughs[j]),0]-X_plot[int(peaks[i]),0])<0:
							cluster_start = X_plot[int(troughs[j]), 0]
							
				for j in range(len(troughs)):
						if (X_plot[int(troughs[j]),0]-X_plot[int(peaks[i]),0])*(X_plot[int(troughs[j]),0]-X_plot[int(peaks[i+1]),0])<0:
							cluster_end = X_plot[int(troughs[j]), 0]
							
			### 
			time_tmp = []
			los_tmp = []				
			for j in range(len(los_r)):
				if los_r[j]>=cluster_start and los_r[j]<cluster_end:
					time_tmp.append(time[j])
					los_tmp.append(los[j])
			
			if peaks[i]==main_peak:
				cl = 'red'
				z = 1
			else:
				cl = color_list[i]
				z = 100				
	
			ax2.scatter(time_tmp, los_tmp, color=cl)			
			ax3.scatter(time_tmp, los_tmp, color=cl, zorder=z)			
				
					
			
ax2.set_ylim([ymin2, ymax2])
ax2.set_xlim([xmin2, xmax2])
ax2.set_xlabel('Year')
ax2.set_ylabel('LOS displacement (mm)')

ax3.set_ylim([ymin2, ymax2])
ax3.set_xlim([xmin2, xmax2])
ax3.set_xlabel('Year')
ax3.set_ylabel('LOS displacement (mm)')			
				
				
f_out = open('ts_example_corrected.txt', 'w')
  
for i in range(len(time_new)):
	f_out.write("%f\t%f\n"%(time_new[i], los_new[i]))	

f_out.close()

fig.savefig('ts_phase_jump_corrected.png', dpi=300, bbox_inches='tight')





	



