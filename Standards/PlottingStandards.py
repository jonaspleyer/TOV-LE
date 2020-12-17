from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'serif':['Computer Modern Roman']})

linestyles = [
	'solid',
	'dashed',
	'dotted',
	'dashdot',
	(0, (5, 10)),
	(0, (1, 10)),
	(0, (3, 10, 1, 10)),
	(0, (5, 1)),
	(0, (1, 1)),
	(0, (3, 1, 1, 1)),
	(0, (3, 5, 1, 5)),
	(0, (3, 5, 1, 5, 1, 5)),
	(0, (3, 10, 1, 10, 1, 10)),
	(0, (3, 1, 1, 1, 1, 1))]

linecolours = [
	'k',
	'r',
	'g',
	'b'
	]