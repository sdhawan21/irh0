#a series of functions to fit the lightcurves used in the H0 analysis.
#Dependencies: SNooPy (Burns et al. 2011)
#pypi installation should be sufficient

import numpy as np 
import matplotlib.pyplot as plt 
import sys
import glob
import snpy	#required for the pymc GP implementation
import matplotlib.gridspec as gridspec


calib = glob.glob('calib_lcs/*.dat')
hflow = glob.glob('hflow_lcs/*.dat')

#a dictionary with parameters to fit and write to a file the results for the calibrators and the 
##hubble flow
input_params = {'fithubble':True, 'writehubble':False, 'fitcalib':False, 'writecalib':False}

def fit_multiple_gps(gl = calib, printres=True, plotlc=True, out_direc='paper_plots/hflow/', out_direc_fits='hflow_lcs/fits/'):
	"""
	fits gaussian processes with a Matern kernel for a series of lightcurves
	Input: Lightcurves in a format that SNooPy can read
	Options:	
		Print the result
		Plot the light curve
	Output: 4 column array of name, z (not required for calibrators), mag, mag_err
	"""

	arr_out = []
	for i in gl:
		try:
			#load the SN
			ss = snpy.get_sn(i)

			#use the default GP fitting algorithm from Pymc written in SNooPy
			ss.J.template(method='gp', Nboot=500)

			#correct for E(B-V) MW
			jm_ebvcorr = ss.J.Mmax-0.81*ss.EBVgal

			if printres:
				print ss.name, ss.J.Mmax, ss.J.e_Mmax

			#use the plotting function to output the file for each SN that is in Appendix A
			if plotlc:
				do_gpfit_plot(ss, direc=out_direc)
			ss.save(out_direc_fits+ss.name+'.snpy')
			#redshift is superfluous for the calibrators but the code outputs it anyhow
			arr_out.append([ss.name,ss.z, round(jm_ebvcorr, 3), round(ss.J.e_Mmax, 3)])
		except:
			ss.name
	return np.array(arr_out)


def fit_tuned(arr, direc = 'hflow_lcs/fits/', plotlc=True):
	"""
	Fits the specific lightcurves of the 4 SNe that didnt have satisfactory plot 
	"""
	lf = snpy.get_sn('hflow_lcs/spl/SN2006lf.dat')
	lf.J.template(method='gp', scale=5., amp = lf.J.mag.std(), diff_degree=3, Nboot=500)
	arr.append([lf.name, lf.z, lf.J.Mmax - 0.81*lf.EBVgal, lf.J.e_Mmax])
	if plotlc:
		do_gpfit_plot(lf)		

	#lf.save(direc+lf.name+'.snpy')


	eq = snpy.get_sn('hflow_lcs/spl/SN2005eq.dat')
	eq.J.template(method='gp', scale=10., amp = eq.J.mag.std()/2., diff_degree=3., Nboot=500)
	arr.append([eq.name, eq.z, eq.J.Mmax - 0.81*eq.EBVgal, eq.J.e_Mmax])
	if plotlc:
		do_gpfit_plot(eq)
	#eq.save(direc+eq.name+'.snpy')

	#since Barone-Nugent et al. 2012 already correct for MW extinction, we dont correct these objects
	ufj = snpy.get_sn('hflow_lcs/spl/ptf10ufj.dat')
	ufj.J.template(method='gp', scale=5., amp = ufj.J.mag.std()*1.5, diff_degree=4, Nboot=500)
	arr.append([ufj.name, ufj.z, ufj.J.Mmax, ufj.J.e_Mmax])
	if plotlc:
		do_gpfit_plot(ufj)
	#ufj.save(direc+ufj.name+'.snpy')

	mwb = snpy.get_sn('hflow_lcs/spl/ptf10mwb.dat')
	mwb.J.template(method='gp', scale=10., amp = mwb.J.mag.std()/2., diff_degree=3, Nboot=500)
	arr.append([mwb.name, mwb.z, mwb.J.Mmax, mwb.J.e_Mmax])
	if plotlc:
		do_gpfit_plot(mwb)
	#mwb.save(direc+mwb.name+'.snpy')
	return arr

def do_gpfit_plot(sn, direc = "paper_plots/hflow/", n=100, xy=(.08, .1)):
	"""
	This function plots the lightcurves and the best fit gaussian process error region for 
	Input: SN from the SNooPy class
	Optional:
		direc: output directory to store the plots in 
		n: number of draws
	Output: PDF figure with the name "sn.name+.pdf"
	"""
	sn.Tmax = sn.J.Tmax
	#setup the plotting environment from the data and the fit
	gs = gridspec.GridSpec(2, 1,height_ratios=[4, 1])

	#plot for the lightcurve
	fig = plt.figure(figsize=(12,8))
	plt.subplot(gs[0]) 
	plt.ylabel("Magnitude (mag)", fontsize=20)
	plt.errorbar(sn.J.MJD - sn.Tmax, sn.J.mag, sn.J.e_mag, fmt='go', alpha=.5)
	#fix the y-axis limits based on the observations (to no "cutoff" the peak of the fit)
	plt.ylim(sn.J.mag.max() + 0.5, sn.J.mag.min() - 0.5)

	#annotate the name of the SN on the bottom left
	plt.annotate(sn.name, xy=xy, xycoords='axes fraction', fontsize=27)

	#get an evenly space time axis on which to use the interpolator
	t = np.linspace(sn.J.MJD.min(), sn.J.MJD.max(), 100)

	#interpolate the data
	mlc = sn.J.interp(t)[0]
    
	interp_mags = sn.J.interp(sn.J.MJD)[0]
    
	#draw n number of realisations from the fit and the error region (n is set to 100)
	ts = []
	for i in range(n):
		#"draw" function  creates a realisation
		sn.J.interp.draw()
		ts.append(sn.J.interp(t)[0])
	ts = np.array(ts)

	#calculate the standard deviation from the draws
	slc = np.std(ts, axis=0)

	
	#use fill_between for plotting
	plt.fill_between(t-sn.Tmax, mlc-slc, mlc+slc, alpha=0.5)
	frame1 = plt.gca()
	for xlabel_i in frame1.axes.get_xticklabels():
		xlabel_i.set_visible(False)
		xlabel_i.set_fontsize(0.0)
	plt.subplot(gs[1])
	res = sn.J.mag - interp_mags
	new_slc = sn.J.e_mag
	plt.errorbar(sn.J.MJD-sn.Tmax,res, new_slc, fmt='go', alpha=0.5)
	plt.plot(t-sn.Tmax, np.zeros_like(t))
	print "Res done"
	plt.ylim(min((res-new_slc))-0.1, max(res+new_slc)+0.1)
	plt.ylabel("Res", fontsize=20)
	plt.xlabel("Days from J-band maximum")
	gs.update(hspace=0.05)
	plt.savefig(direc+sn.name+'.pdf')


def do_kcorr(snfile, ext=False, Tmax=0, dm15=1.1):
	"""
	A function to calculate the K-corrections for an SN from the raw photometry
	Input: The SN file
	Optional: If there is a B-band light curve continue with ext = False
	else, input the Tmax and the dm15 for the SN (from other sources in the literature)
	"""
	sn = snpy.get_sn(snfile)
	sn.J.template(method='gp')
	if ext:
		sn.Tmax = tmax
		sn.dm15 = dm15
	else:
		sn.B.template(method='gp')
		sn.Tmax = sn.B.Tmax
		sn.dm15 = sn.B.dm15
	sn.kcorr(['J'])
	return sn.ks['J']

def main():
	"""
	Main function to do the fits and make output files
	"""
	#set the conditionals for the calibrators and the hubble flow 
	if input_params['fitcalib']:
		calib_fits = fit_multiple_gps(gl=calib, out_direc='paper_plots/calib/', out_direc_fits='calib_lcs/fits/')

	if input_params['fitcalib'] and input_params['writecalib']:
		np.savetxt('calib_fromscript.dat', calib_fits, fmt='%s', delimiter='\t')

	if input_params['fithubble']:
		#fit the hubble flow objects
		hflow_fits = list(fit_multiple_gps(gl=hflow))
		
		#fit the objects that need the parameters to be changed
		hflow_fits = fit_tuned([])
		hflow_fits = np.array(hflow_fits)

		#sort by the SN name
		hflow_fits = hflow_fits[hflow_fits[:,0].argsort()]

	#if the option to write is set, use savetxt to save the fits	
	if input_params['fithubble'] and input_params['writehubble']:
		np.savetxt('hubbleflow_fromscript.dat', hflow_fits, fmt='%s', delimiter='\t')

	

if __name__ == "__main__":
	main()
