import numpy as np
import pylab as plb
import matplotlib.pyplot as plt
import scipy
import scipy.signal
from scipy import optimize
from scipy import integrate
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.interpolate import make_interp_spline
from scipy import stats
import os
import os.path
import ast
import glob

for f in glob.glob("Results\Peak Data\*"):
    os.remove(f)
for f in glob.glob("Results\Peak Plots\*"):       # Delete any previously saved results, to avoid potential
    os.remove(f)                                  # confusion between different data sets (filenames are not unique)
for f in glob.glob("Results\Raw Data\*"):
    os.remove(f)

cur_dir = os.getcwd() # Finds the current directory of where the python file is installed
back_dir = cur_dir + "\Background\\" + os.listdir(cur_dir + "\Background\\")[0] # Finds directory of background reading
meas_dir = cur_dir + "\Measured\\" + os.listdir(cur_dir + "\Measured\\")[0]     # and measurement reading

backtime = int(open(back_dir, "r").readline().strip())
meastime = int(open(meas_dir, "r").readline().strip())  # Pull the live times from both files

with open(back_dir, "r") as b:
    with open(meas_dir, "r") as m:
        countssource = [int(i) for i in m] # create a list for both the source and background reading
        countsback = [int(j) for j in b]

del countssource[0]
del countssource[0]
del countsback[0]
del countsback[0]    # remove the timestamps from the start of the files

if os.stat("MCA Channel Calibration\calibration.txt").st_size != 0:
    with open("MCA Channel Calibration\calibration.txt", "r") as file:
        calidata = file.read().splitlines()
    unit = "(keV)"
    calidata = list(map(float, calidata)) # obtains calibration data for MCA to keV units conversion, if available.
else:
    calidata = []
    unit = "(MCA Channels)"
if len(calidata) != 4:
    if len(calidata) != 0:
        raise ValueError("Incorrect Calibration Data")

peakssource = scipy.signal.find_peaks(countssource,height= 50, distance= 100)[0] # find any peaks in Source array.
# This function will give the indices of the central points of the peaks

ranger = 18
peakxrange = []    # Zooms into each individual peak, using "ranger" to define the range
peakyrange = []
for i in range(len(peakssource)):
    rmin = peakssource[i] - ranger
    if rmin < 0:
        rmin = 0
    rmax = peakssource[i] + ranger
    if rmax > (len(countssource) - 1):
        rmax = len(countssource) - 1            # takes the counts of each peak, and matches with corresponding
    if scipy.stats.normaltest(countssource[rmin:(rmax+1)])[1] < 0.045:    # MCA channels
        peakyrange.append(countssource[rmin:(rmax+1)])    # Applies a statistical test, so that only normally
        peakxrange.append(range(rmin,rmax+1))             # distributed datasets are accepted

backxrange = peakxrange
backyrange = []
for i in backxrange:                                    # For each accepted peak in the source counts, the corresponding
    backyrange.append([countsback[j] for j in i])       # background counts is found


def gauss(x, a, xc, w, y0):
    return a*exp(-(x-xc)**2/(2*abs(w)**2))+y0 # Define Gaussian function to be plotted. Parameters will be estimated


def linear(x, m, c):
    return m*x + c   # Define linear function, just incase a Gaussian cannot be fitted


for i in range(len(peakxrange)):

    x = peakxrange[i]
    xc = peakxrange[i][scipy.signal.find_peaks(peakyrange[i],distance= 36)[0][0]]
    w = scipy.signal.peak_widths(peakyrange[i], scipy.signal.find_peaks(peakyrange[i],distance= 36)[0], 0.5)[0][0]
    FWTM = scipy.signal.peak_widths(peakyrange[i], scipy.signal.find_peaks(peakyrange[i], distance=36)[0], 0.05)[0][
        0]  # Define guesses for parameters of Gaussian

    wback = scipy.signal.peak_widths(backyrange[i], scipy.signal.find_peaks(backyrange[i],distance= 36)[0], 0.5)[0][0]
    FWTMback = scipy.signal.peak_widths(backyrange[i], scipy.signal.find_peaks(backyrange[i], distance=36)[0], 0.05)[0][0]

    xback = x
    xcback = xc

    prom = scipy.signal.peak_prominences(peakyrange[i], scipy.signal.find_peaks(peakyrange[i],distance= 36)[0])
    promback = scipy.signal.peak_prominences(backyrange[i], scipy.signal.find_peaks(
        scipy.signal.savgol_filter(backyrange[i],3,2),distance= 36)[0])  # finds the height of each peak, and the
                                                                         # positions of the peak bases
    h = prom[0][0]
    hback = promback[0][0]
    y = peakyrange[i]
    el = scipy.signal.peak_widths(peakyrange[i], scipy.signal.find_peaks(peakyrange[i],distance= 36)[0], 1)[1][0]
    yback = backyrange[i]
    elback = scipy.signal.peak_widths(backyrange[i], scipy.signal.find_peaks(peakyrange[i],distance= 36)[0], 1)[1][0]

    try:

        popt, pcov = optimize.curve_fit(gauss, x, y, [h, xc, w, el], method="trf")  # Fit a curve to each peak
        try:
            poptback, pcovback = optimize.curve_fit(gauss, xback, yback, [hback, xcback, wback, elback], method="trf")
            func = 0
            perrback = np.sqrt(np.diag(pcovback))
            netcountserrback = np.sqrt(np.pi * ((poptback[2] * perrback[0]) ** 2 + (poptback[0] * perrback[2]) ** 2 + 2
                                                * poptback[0] * poptback[2] * pcovback[0][2]))
            if netcountserrback > 2000:
                poptback, pcovback = optimize.curve_fit(linear, xback, yback, [0, elback], method="trf")
                func = 1
        except:
            poptback, pcovback = optimize.curve_fit(linear, xback, yback, [0, elback], method="trf")
            func = 1

        perr = np.sqrt(np.diag(pcov))
        perrback = np.sqrt(np.diag(pcovback))

        peakintegral = integrate.quad(lambda t: gauss(t, *popt), x[prom[1][0]], x[prom[2][0]])
        netcounts = peakintegral[0]
        netcountserr = np.sqrt(2 * np.pi * ((popt[2] * perr[0]) ** 2 + (popt[0] * perr[2]) ** 2 + 2 * popt[0] * popt[2]
                                           * pcov[0][2]))

        if func == 0:
            peakintegralback = integrate.quad(lambda t: gauss(t, *poptback), x[prom[1][0]], x[prom[2][0]])
            netcountserrback = np.sqrt(np.pi * ((poptback[2] * perrback[0]) ** 2 + (poptback[0] * perrback[2]) ** 2 + 2
                                                * poptback[0] * poptback[2] * pcovback[0][2]))
        else:
            peakintegralback = integrate.quad(lambda t: linear(t, *poptback), x[prom[1][0]], x[prom[2][0]])
            netcountserrback = np.sqrt((perrback[0]/2) ** 2 + (perrback[1] ** 2))

        netcountsback = peakintegralback[0]  # Integrate each peak to find the net counts

        yplot = [s - (l * meastime / backtime) for s, l in zip(peakyrange[i], backyrange[i])]
        FWHM = scipy.signal.peak_widths(yplot, scipy.signal.find_peaks(peakyrange[i],distance= 36)[0], 0.5)[0][0]

        with open("Results\Peak Data\peak " + str(i + 1) + ".txt", "w") as file:
            if len(calidata) != 0:  # background is normalised, and subtracted from the counts
                file.write("Net Counts:            " + str((netcounts - netcountsback*meastime/backtime))
                        + " +/- " + str(np.sqrt( netcountserr**2 + (netcountserrback*meastime/backtime)**2 )
                                        ) + " (Counts)")
                file.write("\nNet Count Rate:        " + str((netcounts - netcountsback * meastime / backtime)/meastime)
                           + " +/- " + str(np.sqrt(netcountserr ** 2 + (netcountserrback * meastime / backtime) ** 2)
                                           /meastime) + " (Counts/s)")
                file.write("\nPeak Height:           " + str(popt[0]) + " +/- " + str(perr[0]) + " (Counts)")
                file.write("\nPeak Central Position: " + str(popt[1]*calidata[0] + calidata[2]) + " +/- " + str(
                    np.sqrt(perr[1] ** 2 + calidata[3] ** 2)) + " " + unit)
                file.write("\nPeak FWHM:             " + str(FWHM*calidata[0]) + " +/- " + str(
                    np.sqrt(perr[2] ** 2 + calidata[3] ** 2)) + " " + unit)
            else:
                file.write("Net Counts:            " + str(netcounts - netcountsback*meastime/backtime) + " +/- " +
                        str(np.sqrt( netcountserr**2 + (netcountserrback*meastime/backtime)**2 )) + " (Counts)")
                file.write("\nNet Count Rate:        " + str((netcounts - netcountsback * meastime / backtime)/meastime)
                           + " +/- " + str(np.sqrt(netcountserr ** 2 + (netcountserrback * meastime / backtime) ** 2)
                                           /meastime) + " (Counts/s)")
                file.write("\nPeak Height:           " + str(popt[0]) + " +/- " + str(perr[0]) + " (Counts)")
                file.write("\nPeak Central Position: " + str(popt[1]) + " +/- " + str(perr[1]) + " " + unit)
                file.write("\nPeak FWHM:             " + str(FWHM) + " +/- " + str(perr[2]) + " " + unit)
        with open("Results\Raw Data\peak " + str(i + 1) + ".txt", "w") as file:
            for (x1, y1) in zip(x, yplot):
                file.write("{0}\t{1}\n".format(x1, y1))  # Saves all useful data from each peak

        if len(calidata) != 0:
            x = [j*calidata[0] + calidata[2] for j in x]
            xc = xc*calidata[0] + calidata[2]
            w = FWHM*calidata[0]

        poptplot, pcovplot = optimize.curve_fit(gauss, x, yplot, [h, xc, w, el], method="trf") # creates a new dataset,
        # with subtracted normalised background, so that a representation of each peak can be shown
        xsmooth = np.linspace(x[0], x[len(x) - 1], len(x) * 100)
        smooth = make_interp_spline(x, yplot, k=3)
        ysmooth = smooth(xsmooth) # Uses fitted function to create more data points, to plot a smooth line

        plt.plot(x, yplot, "b+", label="Counts")
        plt.plot(xsmooth, gauss(xsmooth, *poptplot), "k:", label="Best Fit")
        plt.legend()
        plt.xlabel("Energy " + unit)
        plt.ylabel("Counts")
        plt.savefig("Results\Peak Plots\Fig." + str(i + 1) + " - Full Energy Peak.png")
        plt.show() # Plots and saves peak

    except:
        pass  # If any of the fitting process fails due to "false peaks", the program will move on without stopping