# # # # EXTRACT/COMPARE DATA 

def make_mask( shp_fn, raster_fn ):
	shp = gpd.read_file(shp_fn)

	with rasterio.open(raster_fn) as rst:
		shp = shp.to_crs(rst.crs.to_dict())
		shapes = shp.geometry.tolist()
		mask = rasterize( shapes, out_shape=rst.shape, fill=0, transform=rst.transform, 
					all_touched=True, default_value=1 )
	shp = None					
	return mask

def extract(raster_fn, mask):
	with rasterio.open( raster_fn ) as rst:
		arr = rst.read(1).astype(np.float32)
		arr[arr == rst.nodata] = np.nan
	return np.nanmean(arr[mask == 1])

def run( x ):
	return extract(*x)

if __name__ == '__main__':
	import matplotlib
	matplotlib.use('agg')
	from matplotlib import pyplot as plt
	import os, rasterio, glob
	from rasterio.mask import mask
	from rasterio.features import rasterize
	import geopandas as gpd
	import numpy as np
	import pandas as pd
	from functools import partial
	import multiprocessing as mp

	out_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/climatologies/prism_climatology_assessment/plots'
	shp_fn = '/workspace/Shared/Tech_Projects/DOT/project_data/NOAA_Atlas14/shapefiles/AK_Extent_3338.shp'

	final = dict()
	for group in ['1961-1990','1971-2000','1981-2010']:
		files = glob.glob('/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/climatologies/prism_climatology_assessment/{}/GTiff/*.tif'.format(group))
		files = [fn for fn in files if not 'anl' in fn and '14' not in fn ]
		tmp_fn = files[0]

		# make mask
		mask = make_mask(shp_fn, tmp_fn)
		args = [ (fn, mask) for fn in files ]
		
		# run the averaging in parallel
		pool = mp.Pool( 32 )
		out = pool.map( run, args )
		pool.close()
		pool.join()

		# handle naming to common index conversion
		if group == '1961-1990':
			idx_tmp = [ ''.join(np.array(os.path.basename(fn).split('.')[0].split('_'))[[0,-3]]) for fn in files ]
			idx = []
			for i in idx_tmp:
				if 'tasmin' in i:
					i = i.replace('tasmin', 'tmin')
				elif 'tasmax' in i:
					i = i.replace('tasmax', 'tmax')
				elif 'tas' in i:
					i = i.replace('tas', 'tmean')
				elif 'pr' in i:
					i = i.replace('pr', 'ppt')
				idx = idx + [i]
		if group == '1971-2000':
			idx = [ os.path.basename(fn).split('.')[0] for fn in files ]
		if group == '1981-2010':
			idx = [ ''.join(np.array(os.path.basename(fn).split('.')[0].split('_'))[[1,-1]]) for fn in files ]

		# put it in the dict
		final[group] = dict(list(zip(idx,out)))

	# make a data frame with the results and the common index we added
	df = pd.DataFrame(final)

	# do this awkwardly, but effectively with a groupby
	grouper = [ i[:-2] for i in df.index ]
	for grp,sub_df in df.groupby(grouper):
		months = [ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
		sub_df.index = months
		out_fn = os.path.join(out_path, 'compare_prism_baseline_periods_{}.png'.format(grp))
		sub_df.to_csv( out_fn.replace('.png', '.csv') )
		title = 'Comparison of PRISM {} Monthly Climatologies\nAlaska Statewide Domain'.format(grp)
		ax = sub_df.plot( kind='line', title=title )
		ax.set_xticks(list(range(len(sub_df.index))))
		ax.set_xticklabels(sub_df.index)
		if grp == 'ppt':
			ylabel = 'Total Precipitation (mm)'
		else:
			ylabel = 'Temperature ($^\circ$C)'
		ax.set_ylabel(ylabel)
		plt.savefig( out_fn )
		plt.close()
		plt.cla()



