# GammaRaySpectroscopy
Search and analysis software for gamma-ray spectroscopy

SET UP:

  Program requires a backround spectrum and a sample spectrum to work. Both files should be of .TKA or .txt file format, and both should have
  the same dimensions. The first two lines of each contain the live and real times respectively. The rest of the lines should encode the 
  gamma ray spectrum: each line contains the intensity for one bin, and each new line represents an increase of 1 MCA Channel, starting at
  line 3 as 1 MCA Channel.
  
  The background spectrum should be placed in the "Background" folder, and the sample spectrum in the "Measured" folder. If the detector has 
  been calibrated from MCA Channels to keV as follows:
    
                                                         keV = m*MCA + c
                                                         
  Where m and c are calculated with their errors, create a text file with the following on separate lines: m, m error, c, cerror, and store
  in the folder "MCA Channel Calibration".
  
RUNNING:

  Run "CosmicAnalyser2.py". This will detect, analyse, and plot photopeaks above background in the sample spectrum. A series of peaks will
  be shown to the user as they are found. After the program has run, the user can look through the "Results" folder to find all of the 
  photopeak data that was measured.
  
RESULTS:
 
  "Peak Data" contains a text file for each peak, each contains the net counts, net count rate, peak height, peak central position, and
  the FWHM of the peak, all with their uncertainties.
  
  "Peak Plots" contains a plot of each photopeak
  
  "Raw Data" contains a text file for each photopeak, containing the data points used to create the plots. If the user desires, this can
  be imported into software such as QTIplot, should they want to make their own plots. 
  
