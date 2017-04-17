# this small script will take a root directory of xml metaadata files output by manually extracting a large number of
# metadata file, with particular mind to the REA work.  These output files have a deep directory structure with 
# non-unique filenames for the output metadata records (all named metadata.xml). This script will walk the subdirs in 
# this root directory and will generate filenames based on the information in the xml file and output the new file-
# name of the xml as determined by the file contents and dump it out to the output_path directory.
#
#  Written by: Michael Lindgren (malindgren@alaska.edu)
#  Purpose: REA Project Data Delivery Automation
#	
# [NOTE] this software requires an outside dependency named BeautifulSoup4
# 		This is installed simply with pip (setuptools) using: pip install beautifulsoup4
#
# # # # # 
def new_fn_from_meta( meta_fn, output_path ):
	'''
	take a SNAP Geonetwork created xml metadata file for the REA
	projects, extract the title, work some magic to generate a 
	proper output filename, write out the new file with the new
	name to the output_path.

	arguments:
		meta_fn = [str] path to the input xml file 
		output_path = [str] path to output directory

	returns:
		output filename as a string if successful, with the 
		side-effect of generating a file on disk in the 
		output_path location.

	'''
	import os, sys, re, glob, lxml
	from bs4 import BeautifulSoup

	# open using the lxml parser
	soup = BeautifulSoup( open( meta_fn ).read(), 'lxml-xml' )
	soup_select, = soup.find_all( 'identificationInfo' )
	
	# string formatting area this is plug and play as long as we return a 
	# useable string file basename ( no path )
	title = soup_select.title.get_text().strip()
	title = re.sub( '[^A-Za-z0-9]+', ' ', title ).replace( ' ', '_' ) + '.xml'
	out_fn = os.path.join( output_path, title.lower() )

	# we can do some other lookup-table-ish stuff here to generate proper 
	# filenames.  There cant be too many combinations and that is wayyy easier
	# than clicking a thousand times.
	# -------
	# write it out
	output_filename = os.path.join( output_path, out_fn )
	with open( meta_fn ) as input_meta:
		with open( out_fn, 'w' ) as output_meta:
			output_meta.write( input_meta.read() )
	return output_filename
	
if __name__ == '__main__':
	'''
	install BeautifulSoup4 before running or it will fail
	this tool is run at the commandline using the following syntax

	- cd to the directory containing the script
	- python walk_rename_meta.py -root {path_to_root_meta_directory} -out {path_to_output_directory}
	- running the script like:
		 python walk_rename_meta.py --help 
		will show the help for the commandline tool in case that is needed.
	'''

	# read in the needed packages
	import os, sys, re, argparse, glob
	from bs4 import BeautifulSoup

	# # set up an argument parser
	parser = argparse.ArgumentParser( description='program to find and rename metadata files downloaded from SNAPs GeoNetwork for REA projects' )
	parser.add_argument( "-root", "--root_path", action='store', dest='root_directory', type=str, help="path to metadata root directory" )
	parser.add_argument( "-out", "--output_path", action='store', dest='output_path', type=str, help="path to output directory" )
	args = parser.parse_args()
	
	# unpack the args
	root_directory = args.root_directory
	output_path = args.output_path

	# walk the directory and output the metadata with a new filename in the output directory
	for dir_name, subdir_list, file_list in os.walk( root_directory ):
		files = [ os.path.join(dir_name, fn) for fn in file_list if 'metadata.xml' in fn ]
		# l = glob.glob( os.path.join( dir_name, 'metadata.xml' ) )
		if len( files ) > 0:
			print files[0]
			# run the function in a list comprehension
	 		[ new_fn_from_meta( meta_fn, output_path ) for meta_fn in files ]
	print 'complete! new files located in %s' % output_path


# # # TESTING STUFF BELOW # # # 
# output_path = '/Users/malindgren/Documents/downscale_epscor/nov_fix/metadata/meta_export'
# root_directory = '/Users/malindgren/Documents/downscale_epscor/nov_fix/metadata/meta_export'

# python /Users/malindgren/Documents/repos/downscale/snap_scripts/epscor_sc/walk_rename_meta.py -root /Users/malindgren/Documents/downscale_epscor/nov_fix/metadata/meta_export -out /Users/malindgren/Documents/downscale_epscor/nov_fix/metadata
# python walk_rename_meta.py --help
# # # # # # # # # # # # # # # #


