## python3

import numpy as np
#import matplotlib.gridspec as gridspec
#import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from scipy.signal import find_peaks
import scipy
from numpy.polynomial.polynomial import polyfit


## read time series
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


###### normalize data
time_min = min(time)
time_max = max(time)
los_min = min(los)
los_max = max(los)

for i in range(len(time)):
	time[i] = (time[i]-time_min)/(time_max-time_min)
	los[i] = (los[i]-los_min)/(los_max - los_min)


## get the slope of the original time series using best-fit line
b, s = polyfit(time, los, 1)
s_degree = np.arctan(s)*180/np.pi
print ("original slope: ", s_degree)	
	

## rotation the time series
range_step = 60. # degree, a_open, the angle constraining the final slope range of the rotated time series
theta_range = range(int(s_degree-range_step), int(s_degree+range_step+1), 1)


f_out = open('clusters_distances_height.txt', 'w')
f_out.write("rotation angle (unit: degrees)	average_peak_distance_residual d_res(unit: mm)	average_peak_height_hp\n")  


for theta_degree in theta_range:
	print ("rotation angle: ", theta_degree)
	theta = theta_degree*np.pi/180.
	time_r = [0.]*len(time)
	los_r = [0.]*len(los)
	for i in range(len(time_r)):
		time_r[i] = time[i]*np.cos(theta)+los[i]*np.sin(theta)
		los_r[i] = -1.*time[i]*np.sin(theta)+los[i]*np.cos(theta)

	#time_r_shift = time_r[0]
	time_r_shift = min(time_r)
	los_r_shift = min(los_r)
	for i in range(len(time_r)):
		time_r[i] = time_r[i]-time_r_shift
		los_r[i] = los_r[i]-los_r_shift

	
	ymin = -0.1	
	ymax = 1.5
	
	xmin = -0.1	
	xmax = 1.5

	'''
	##### start of plotting
	
	fig = plt.figure(figsize=(15,4))
	gs = gridspec.GridSpec(1,3)
	gs.update(wspace=0.4, hspace=0.0)

	####### plot original normalized time series
	ax = fig.add_subplot(gs[0])
	plt.title('Original normalized time series')
	ax.scatter(time, los, color='black') 			
	ax.set_ylim([ymin, ymax])
	ax.set_xlim([xmin, xmax])
	#ax.set_aspect('equal')
	ax.set_xlabel('Year (normalized)')
	ax.set_ylabel('LOS displacement (normalized)')
		
	x0 = 0.0
	y0 = 1.1
	length = 0.2
	#ax.text(x0+length*0.2,y0, '%d$^\circ$'%((theta_degree)),color='red') #fontweight='bold'	
	## twin axis
	ax2 = ax.twinx()
	ymax2 = ymax*(los_max-los_min) + los_min
	ymin2 = ymin*(los_max-los_min) + los_min
	ax2.set_ylim([ymin2, ymax2])
	ax2.set_ylabel('LOS displacement (mm)')
	
	
	####### plot rotated time series
	ax = fig.add_subplot(gs[1])	
	plt.title('Rotated time series')
	ax.scatter(time_r, los_r, color='blue') 	
	ax.set_ylim([ymin, ymax])
	ax.set_xlim([xmin, xmax])
	ax.set_xlabel('Year (normalized)')
	ax.set_ylabel('LOS displacement (normalized)')
	
	ax.text(x0+length*0.2,y0+0.2, 'Rotation angle: %d$^\circ$'%((theta_degree)),color='red') 
	
	## twin axis,
	ax2 = ax.twinx()		
	ymax2 = ymax*(los_max-los_min) + los_min
	ymin2 = ymin*(los_max-los_min) + los_min
	ax2.set_ylim([ymin2, ymax2])
	ax2.set_ylabel('LOS displacement (mm)')
	'''

	## histogram
	## # check the link below to find a opitmal bin width for the histogram
	#https://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule
	#https://stats.stackexchange.com/questions/798/calculating-optimal-number-of-bins-in-a-histogram
	iqr = scipy.stats.iqr(los_r)
	his_bin = 2.0*iqr/(len(los_r)**(1/3.0))
	
	
	####### plot histogram and density curve	
	'''
	ax = fig.add_subplot(gs[2])	
	plt.title('Histogram and density curve')
	ax.hist(los_r,bins=30,  edgecolor= 'black', orientation='horizontal', density=True) 
	'''
	
	### calculate gaussian Kernel density	
	# estimate the bandwidth based on the rule-of-thumb bandwidth estimator
	#https://en.wikipedia.org/wiki/Kernel_density_estimation
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
		
	X_plot = np.arange(ymin, ymax, 0.002)[:, np.newaxis]
	log_dens = kde.score_samples(X_plot)	
	dens = np.exp(log_dens)
	#ax.plot(dens, X_plot[:,0], color='black', lw=1.5, linestyle='-')
	


	## detect crests and troughs
	peaks = find_peaks(np.array(dens))[0]
	if len(peaks)>0:
		for j in range(len(peaks)):
			index = int(peaks[j])
			#ax.scatter(dens[index], X_plot[index,0], color='red',zorder=100)
			
	troughs = find_peaks(-np.array(dens))[0]
	if len(troughs)>0:
		for j in range(len(troughs)):
			index = int(troughs[j])
			#ax.scatter(dens[index], X_plot[index,0], color='lime',zorder=100)
	
	
	if len(peaks)>1:
		print ("Mutiple peaks.")
		### get the averaged (weighted) peak distance residual		
		## get the number of points in each peak/group
		n_points = [0]*len(peaks)		
						
		for i in range(len(peaks)):			
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
			n_points[i] = len(los_tmp)
			
		if sum(n_points) != len(los_r):
			print ("error! Total number of points do not match.")
			
		else:
			print (sum(n_points), " = ", len(los_r), ". Total number of points matches.")
						
		
		n_points_sum_non_main = 0.0		
		
		
		## get the main peak (the peak with highest density)		 
		main_peak = peaks[0]
		for j in range(len(peaks)):
			if dens[int(peaks[j])]>dens[int(main_peak)]:
				main_peak = peaks[j]		
		#print ("main peak:",  X_plot[int(main_peak),0])	
		
		n_points_sum_on_main = 0.0
		for j in range(len(n_points)):
			if peaks[j] != main_peak:
				n_points_sum_on_main = n_points_sum_on_main + n_points[j]
				
		#print("Total number of points in non-main groups: ", n_points_sum_on_main)
		
				
		peak_dis_res_ave = 0.0
		for j in range(len(peaks)):		
			if peaks[j] != main_peak:
				### note that, convert peak distance to the real LOS displacement
				peak_dis = abs(X_plot[int(peaks[j]),0] - X_plot[int(main_peak),0])*(los_max - los_min)				
				tmp = []
				amb = []
				for k in range(-10, 10):					
					tmp.append(abs(peak_dis-k*27.75))				

				peak_dis_res = min(tmp)
				
				peak_dis_res_ave = peak_dis_res_ave + peak_dis_res*n_points[j]/n_points_sum_on_main
			
				
		### get the averaged peak height (the height difference between a peak and its nearest trough(s))
		peak_h_ave = 0.0
		for j in range(len(peaks)):
			if j == 0:
				for k in range(len(troughs)):
					if (X_plot[int(troughs[k]),0]-X_plot[int(peaks[j]),0])*(X_plot[int(troughs[k]),0]-X_plot[int(peaks[j+1]),0])<0:		
						trough_mid = troughs[k]
				
				peak_h_ave = peak_h_ave + dens[int(peaks[j])] - dens[int(trough_mid)]
				
			elif j == len(peaks)-1:
				for k in range(len(troughs)):
					if (X_plot[int(troughs[k]),0]-X_plot[int(peaks[j]),0])*(X_plot[int(troughs[k]),0]-X_plot[int(peaks[j-1]),0])<0:
						trough_mid = troughs[k]
				peak_h_ave = peak_h_ave + dens[int(peaks[j])] - dens[int(trough_mid)]
			
			else:
				for k in range(len(troughs)):
					if (X_plot[int(troughs[k]),0]-X_plot[int(peaks[j]),0])*(X_plot[int(troughs[k]),0]-X_plot[int(peaks[j-1]),0])<0:
						trough_left = troughs[k]
							
				for k in range(len(troughs)):
					if (X_plot[int(troughs[k]),0]-X_plot[int(peaks[j]),0])*(X_plot[int(troughs[k]),0]-X_plot[int(peaks[j+1]),0])<0:
						trough_right = troughs[k]				
				
				peak_h_ave = peak_h_ave + dens[int(peaks[j])]-0.5*(dens[int(trough_left)]+dens[int(trough_right)])
		
		peak_h_ave = peak_h_ave/(len(peaks))
		
		f_out.write("%d\t%f\t%f\n"%(theta_degree, peak_dis_res_ave, peak_h_ave))  
		
	
	else: # when there is only a single peak
		print ("Single peak.")
		f_out.write("%d\t-999\t-999\n"%(theta_degree))
			
	'''
	ax.set_ylim([ymin, ymax])
	ax.set_xlabel('Frequency/Probability density (normalized)')
	ax.set_ylabel('LOS displacement (normalized)')
	ax.set_xlim([0, 3.0])

	## twin axis,
	ax2 = ax.twinx()
	## twin axis, two y axis, different y axis
	ymax2 = ymax*(los_max-los_min) + los_min
	ymin2 = ymin*(los_max-los_min) + los_min
	ax2.set_ylim([ymin2, ymax2])
	ax2.set_ylabel('LOS displacement (mm)')		
	
	fig.savefig('ts_rotation%d_histogram.png'%(theta_degree), dpi=300, bbox_inches='tight')
	'''









