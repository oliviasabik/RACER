### RACER 1.0.0
This is the first release of RACER

## Features
# RACER has three major classes of functions
# (1) A formatting function. formatRACER; this function takes association data as an input and formats it to be compatible with downstream functions
# (2) A function for pulling linkage disequilibirum data. ldRACER; this function uses rs IDs from the input data to query the 1000G phase III data for a lead SNP, either specified or the most significant SNP in the input data and it's proxies.
# (3) Three functions for plotting your data:
	# (A) singlePlotRACER; this function makes a plot of one set of association data.
	# (B) scatterPlotRACER; this function makes a scatter plot, with p-values from two different studies on the x- and y-axes.
	# (C) mirrorPlotRACER; this function makes a mirror plot, composed of one inverted association plot and one standard association plot, stacked on top of one another, comparing the two associations.

