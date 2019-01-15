idx=find_idx(xi, xgrid)

Purpose: Fractional binning

Why:
- This function not only bin similarly to HISTC, but also returns
the fractional position of data points within the binning interval.
- It is equivalent to
	interp1(xgrid,(1:length(xgrid)), xi)
  but with the speed improvement up to 5 times.
- Few obvious examples of applications:
	+ binning step for more sophisticated interpolation schemes
	such as multi-dimensional spline or linear tensorial
	interpolations.
	+ Generate discrete random sequences with given probability
	distribution
