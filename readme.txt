

Readme for scripts written by Kristi Webb (2016) for the purpose of analyzing the colours of stars seen in the NIR through observations of molecular cloud cores.

Procedure:

photometry
- determine optimal photometry parameters with a curve of growth
	phot_curves.py
- calculate the zeropoint offset comparing to a 2mass catalogue
	zeropoint.py
- plot colour magnitude diagrams, calculate colours
	colours.py
- determine the limiting magnitude by adding artificial stars 
	addstars.py (incomplete)


Note: the scripts phot_curves.py, zeropoint.py, and colours.py use the python modelue 'argparse' which reads in parameters given in the command line. The options for each script can be displayed easily with the command, for example:
python phot_curves.py --help
where parameters can be specified using the parameter name like so:
python phot_curves.py --image (or -img) image.fits

Note: the script addstars.py does not use the argparse module as the input is quite involved. The input is instead specified in the header of the script.

Note: Each script has the same comments in the header.


---------------
phot_curves.py
---------------

	input: fits image, 
	       list of aperture sizes, 
	       source extractor file
	       (optional) feducial radius
	output files: table of magnitudes of unsaturated unreddened stars
	output images: curve of growth, 
		       half light radius, 
		       mag aper and mag auto difference

   Preform an analysis on the initial photometry to determine the optimal parameters.
   Input consists of a fits image with corresponding source extractor photometric catalogue,
     and name of catalogue to be output.

   This script requires, at minimum, the source extractor values:
   'mag_aper', 'magerr_aper', 'kron_radius', 'x_image', 'y_image', 'fwhm_image' (if r_fed not specified)
   to use aper_auto_check: 'mag_auto', 'mag_err_auto'
   to use half_light: 'flux_radius'

   The apertures used with source extractor must be specified, or the default values are used:
   [5, 10, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60]

   To calculate the curve of growth:
   1. the photometry catalogue is parsed for stars with accurate measurements (i.e. mag < 99.)
   2. the location of each star is compared to the image, where stars that are saturated (have cores with 0's) are removed. 
   3. A curve of growth is then produced for each remaining star, and the final curve is the weighted average 
      (weight =  1/sigma^2). This curve is shown relative to a chosen fiducial radius, default 30, and fitted to an 
      exponential curve. The plot is saved if the parameter 'outcurvegrowth' is specified.

   Additional operations:
   aper_auto_check: calculation of the difference between the measurement of the magnitude from source
   extractors 'mag_aper' for the fiducial radius specified and the 'mag_auto'

   half_light: fits a linear function to the distribution of 'flux_radius' vs. 'mag_aper', removes objects
   outside of 1 standard deviation, fits a linear function again, and removes objects outside 3 standard
   deviations. The fitting is used twice as the large scatter of 'flux_radius' due to highly extended
   objects effects the linear fitting process.


---------------
zeropoint.py
---------------

	input: fits image, 
	       (optional) source extractor file (same file as phot_curves), 
	       (if not source extractor file) output table from phot_curves
	       2mass catalogue file,
	       aperture size selected from curve of growth
	output files: table of magnitudes of unsaturated unreddened stars
	output images: curve of growth, 
		       half light radius, 
		       mag aper and mag auto difference


   Determine the magnitude zeropoint of our CFHT observations by comparing bright unsaturated stars to the same stars
   in the 2mass catalog.
 
   Input includes either: image and corresponding OR the photometry file output from phot_curves.py
    as well as: a 2mass file listing photometric values for bright stars in the region of the image
                - this can be retreived at: http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd
             the aperture to calculate the zeropoint at. Any float may be specified, and the aperture
               closest to an aperture used with source extractor will be used.
               Default source extractor apertures are:
               [5, 10, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60]
             the band of the photometry, used to specify the magnitudes from the 2mass file for comparison to the
               measured magnitudes for the calculation of the zeropoint
    optional:  an outfile to print a few of the calculated values to

   To calculate the zeropoint correction to be entered into source extractor for the correct calculation of magnitudes,
    the process is as follows:
    1. if source extractor and file given:
         parse photometric cataloque to remove saturated stars
       else:
         read in photmetry file output from phot_curves, which has already been parsed
    2. use the cloud names, and shape paremeters specified in the header to remove stars which are in the region of the
        cloud core
    3. remove stars which do not have a star light PSF using the phot_curves half_light function.
    4. Identify the stars in our image remaining with those from the 2mass catalogue, uses a nearest neighbour approach
        with a maximum separation of 0.005 degrees (hardcoded into query_seperation function)
    5. subtract the magnitudes measured from those in the 2mass catalogue, calculate a weighted average


---------------
colours.py
---------------
	input: list of source extractor files for J,H,Ks bands, 
	       (optional) maximum separation (deg) between stars matched by nearest neighbour search, 
	       list of image files for J, H, Ks bands
	output files: table of magnitudes of unsaturated unreddened stars for all J,H,Ks bands
	output images: colour magnitude diagram

   Procedure:
   Match stars from J,H,Ks bands, plot colour magnitude diagrams

   Input consists of:
     a list of source extractor files corresponding to the J,H,K bands
     the fits images in the J,H,K bands
     the name of the photometry file to be output (this will be the input into NICEST)
     the name of the colour colour figure to be saved, by default not saved

     Note: if the output photometry file has already been produced, the script can be run without specifying
        the source extractor or images file, as long as the output photometry file is specified

   Optional input:
     maximum separation for nearest neighbour search, default 0.0001 degrees

   To produce the colour colour plots:
     1. remove saturated stars from each photometry catalogue using phot_curves function
     2. match the stars in each catalogue, only keep stars identified in all bands, uses nearest neighbour approach
     3. compare colours, and plot
     4. if a particular region is specifed (i.e. written directly into the script) the magnitude list can be cropped easily



