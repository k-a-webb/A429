

Readme for scripts written by Kristi Webb for the purpose of analysing the colours of stars seen in the NIR through molecular cloud cores.

Procedure:

photometry
- determine optimal photometry parameters with a curve of growth
	phot_curves.py
- calculate the zeropoint offset comparing to a 2mass catalogue
	zeropoint.py
- determine the limiting magnitude by adding artificial stars 
	incomplete

- plot colour magnitude diagrams
	plt_cmd.py
- calculate the colours of the stars and plot extinction map
	mapping.py