# Description of "plot_3spectra.py"
The script displays a window containing simultaneously three consecutive spectral plots (SDSS .fits files) on the screen. After closing this window, in each subsequent one, the second and third spectrum from the previous window will be displayed, and the new one - the next spectrum in order in current working directory. If there are fewer than three spectra fits files in the working directory, all of them are displayed. Plots include most of the major emission and absorption lines.

# Description of "plotting_SDSS_bin_spectra.py"
Script which convert all .fits (binary) files to .asc (text files) in current directory. Extract main table from SDSS 1d spectra .fits files and save in asci text file. (Column of asci files are same as rows in .fits). 
Then, it displays every spectrum from the working directory to multiplot: without binning, binning bin1 and binning bin2 (by default bin1 = 3, bin2 = 5, you can change the values of these parameters in the script). After closing a given window, the next one is displayed until the last spectra from current working directory. Plots include most of the major emission and absorption lines.

# Required packages:
- python >=3.5
- numpy >=1.0.3
- astropy
- matplotlib

# To use it from the terminal:
Make sure you put the necessary script in the working directory containing your .fits spectra files, 
give it permission (chmod u + x <filename>) and run it ("<filename>" or "python <filename >").
You can try it - in 'spectra' directory I have placed several sample .fits files (WD / QSO / GALAXY).
  
# If You have question, please contact me: 
agn.gurgul@gmail.com
 
