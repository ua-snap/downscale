from setuptools import setup, find_packages

dependencies_list = [ 'xarray','rasterio','pandas','numpy','rasterio','pathos','pynio' ]
scripts_list = []

classifiers = [
		# How mature is this project? Common values are
		#   3 - Alpha
		#   4 - Beta
		#   5 - Production/Stable
		'Development Status :: 3 - Alpha',

		'Intended Audience :: Users',
		'Topic :: Delta Downscaling :: CMIP5 - Integrated Ecosystem Model Inputs',

		# Pick your license as you wish (should match "license" above)
		'License :: OSI Approved :: MIT License',

		# Specify the Python versions you support here. In particular, ensure
		# that you indicate whether you support Python 2, Python 3 or both.
		'Programming Language :: Python :: 2.7'
		]

setup(	name='downscale',
		version='0.01',
		description='tool to downscale CMIP5 model outputs for use in regional climate modeling',
		url='https://github.com/ua-snap/downscale',
		author='Michael Lindgren',
		author_email='malindgren@alaska.edu',
		license='MIT',
		packages=find_packages(),
		install_requires=dependencies_list,
		zip_safe=False,
		include_package_data=True,
		dependency_links=['https://github.com/uqfoundation/pathos','https://github.com/nioinnovation/pynio'],
		scripts=scripts_list,
		classifiers=classifiers,
		test_suite='nose.collector',
		tests_require=['nose']
		)

# if the tests are not running try this:
# $ cd tests
# $ chmod -x *.py
# and run again

