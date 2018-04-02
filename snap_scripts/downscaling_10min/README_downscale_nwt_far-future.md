## HOW TO RUN THE DOWNSCALING FOR THE LONGER-TIMESERIES NWT DOMAIN.
#### WRITTEN BY: MICHAEL LINDGREN (malindgren@alaska.edu)
#### DATE: SEPTEMBER 2017
---
** reference folder should be the snap_scripts/downscale_10min directory in the github repository **

1. Download the data from the ESGF using the synda application installed via the RPM.
	I used Fedora23 on the Phobos (phobos.snap.uaf.edu) server and installed the RPM there.
	See the below link for more docs on how to get that installed and the right ( patched )
	version too.  :)
	https://github.com/Prodiguer/synda/blob/master/sdt/doc/rpm_install.md

	SYNDA Notes:
		- make sure to examine the sdt.conf file and update as necessary: sudo vi /etc/synda/sdt/sdt.conf
		- run with: sudo synda install -s nwt_variables.txt
		- run: sudo synda autoremove (this removes all older versions of datasets)

2. remove the CCSM4 duplicate files that were noticed in the raw downladed files: 
	rm /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5_nwt/data/cmip5/output1/NCAR/CCSM4/rcp60/mon/atmos/Amon/r1i1p1/v20160829/tas/tas_Amon_CCSM4_rcp60_r1i1p1_200501-210012.nc

	rm /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5_nwt/data/cmip5/output1/NCAR/CCSM4/rcp45/mon/atmos/Amon/r1i1p1/v20160829/tas/tas_Amon_CCSM4_rcp45_r1i1p1_200501-210012.nc

3. move the GFDL-CM3 RCP60 tas/pr files from an older download since we cannot seem to access these files from any of the available nodes.
	mkdir /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5_nwt_v2/cmip5_raw_restructure_V2/GFDL-CM3/rcp60

	cp -R /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/raw/GFDL-CM3/rcp60/tas /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5_nwt/cmip5_raw_restructure/GFDL-CM3/rcp60

	cp -R /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/raw/GFDL-CM3/rcp60/pr /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5_nwt/cmip5_raw_restructure/GFDL-CM3/rcp60

	** this should complete the entire set of data to be used in downscaling. **

4. Move the downloaded files from the multiple path locations to a single and navigable directory structure.
	`move_raw_cmip5_common_dir.py`

5. Stack the multi-file datasets for a single model/scenario/variable combination to a single file for the entire available series. and in a navigable directory structure.
--> we do this because some of the datasets have CF-compliant metadata regarding time units that do not mesh properly with the MFDataset class in netCDF4 python's current version. Since this cannot be trusted and we want things to be somewhat clean for downscaling production I am using the `ncrcat` operator from NCO (through python) to do the stacking properly.  Since NCO seems to 'just work'™ I have chosen this route as the best way forward.
	`stack_raw_cmip5_ncrcat.py`

6. Change the name of the CCSM4 directory to the proper NCAR-CCSM4
	cd /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5_nwt/cmip5_raw_ncrcat
	mv CCSM4 NCAR-CCSM4

7. Run downscaling script on these newly stacked files. -- WE ARE STILL STUCK IN PYTHON2 -- use virtenv --> ds/bin/activate (atlas.snap.uaf.edu)
	`wrap_downscaler_cmip5_slurm_nwt_far-futures.py`

8. Run 5ModelAvg script
	`wrap_slurm_compute_5ModelAvg.py`

9. Move anomalies to a common directory away from the downscaled outputs
	`move_anomalies_to_common_dir.py`

10. Run Derived Grids Post-processing
	`wrap_slurm_run_post_downscaling.py` # will run all the needed derived grids.

11. Run clip to NWT domain:
	`crop_clip_to_nwt.py`

	on both the folders holding the derived_grids, _and_ the downscaled files.

12. Run Extractions for MineSites and NWT-wide --> use Python3 v3/bin/activate
	`run_derived_mine_profiles_nwt.py`
	`extract_nwt_aoi_means.py` # ** move the NWT clipped downscaled data to a parent dir named 'grid' -- if that has not been fixed in a subsequent version and an error is noticed at that point. 
	--> `/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/NWT_DELIVERABLES/downscaled`

13.	Plot datasets extracted from the larger set. [ still in progress and evolving ]
	`plot_nwt_aoi_means_v2.py`
	`plot_minesites.py`
