#!~/bin/bash

python /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/compute_decadal_grids_epscor_sc.py -b /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled -o /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/monthly_decadals -m ts323 -s historical -p cru -v tasmin -am mean -nc 32
echo 'DONE!'
python /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/compute_seasonal_annual_grids_epscor_sc.py -b /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled -o /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/annual_seasonals -m ts323 -s historical -p cru -v tasmin -am mean -nc 32
echo 'DONE!'
python /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/compute_seasonal_decadal_grids_epscor_sc.py -b /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/annual_seasonals -o /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/decadal_seasonals -m ts323 -s historical -p cru -v tasmin -am mean -nc 32
echo 'DONE!'
python /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/compute_decadal_grids_epscor_sc.py -b /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled -o /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/monthly_decadals -m ts323 -s historical -p cru -v tasmax -am mean -nc 32
echo 'DONE!'
python /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/compute_seasonal_annual_grids_epscor_sc.py -b /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled -o /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/annual_seasonals -m ts323 -s historical -p cru -v tasmax -am mean -nc 32
echo 'DONE!'
python /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/compute_seasonal_decadal_grids_epscor_sc.py -b /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/annual_seasonals -o /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/decadal_seasonals -m ts323 -s historical -p cru -v tasmax -am mean -nc 32
echo 'DONE!'
python /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/compute_decadal_grids_epscor_sc.py -b /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled -o /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/monthly_decadals -m ts323 -s historical -p cru -v tas -am mean -nc 32
echo 'DONE!'
python /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/compute_seasonal_annual_grids_epscor_sc.py -b /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled -o /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/annual_seasonals -m ts323 -s historical -p cru -v tas -am mean -nc 32
echo 'DONE!'
python /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/compute_seasonal_decadal_grids_epscor_sc.py -b /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/annual_seasonals -o /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/decadal_seasonals -m ts323 -s historical -p cru -v tas -am mean -nc 32
echo 'DONE!'
python /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/swi_calc_decadal_epscor_sc.py -b /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/monthly_decadals -o /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/swi_decadals -m ts323 -s historical -p cru -v tas