

Readme for scripts written by Kristi Webb for the purpose of analysing the colours of stars seen in the NIR through molecular cloud cores.

Procedure:

photometry
- determine optimal photometry parameters with a curve of growth
	phot_curves.py
- calculate the zeropoint offset comparing to a 2mass catalogue
	zeropoint.py
- determine the limiting magnitude by adding artificial stars 
	addstars.py (incomplete)
- plot colour magnitude diagrams, calculate colours
	colours.py



phot_curves.py
	input: fits image, 
	       list of aperture sizes, 
	       source extractor file
	       (optional) feducial radius
	output files: table of magnitudes of unsaturated unreddened stars
	output images: curve of growth, 
		       half light radius, 
		       mag aper and mag auto difference

zeropoint.py
	input: fits image, 
	       (optional) source extractor file (same file as phot_curves), 
	       (if not source extractor file) output table from phot_curves
	       2mass catalogue file,
	       aperture size selected from curve of growth
	output files: table of magnitudes of unsaturated unreddened stars
	output images: curve of growth, 
		       half light radius, 
		       mag aper and mag auto difference


colours.py
	input: list of source extractor files for J,H,Ks bands, 
	       (optional) maximum separation (deg) between stars matched by nearest neighbour search, 
	       list of image files for J, H, Ks bands
	output files: table of magnitudes of unsaturated unreddened stars for all J,H,Ks bands
	output images: colour magnitude diagram
