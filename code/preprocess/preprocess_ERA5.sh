output_path=/net/xenon/climphys/lbloin/boost_proba

# finds all files in the months of june and july of 2021, which is when the 2021 PNW heatwave occurred
for month in 06 07
        do
        # detrended ERA5 climatology data
        # Convert from Kelvin to Celsius, Remap bilinearly to CESM2 grid
        cdo -b F32 -addc,-273.15 -selname,T2M_clim -setattribute,T2M_clim@units="°C" \
        -remapbil,../../inputs/CESM_atm_grid.txt \
         -cat /net/litho/atmosdyn/hotzb/data/era5/N_clim_runyr9_runhour504/2021/${month}/T2021${month}??_?? \
        ${output_path}/tmp_2021_${month}_temperature_detrended_clim_ERA5.nc
        # instantaneous ERA5 data
        ## TEMPERATURE
        cdo -b F32 -addc,-273.15 -selname,T2M -setattribute,T2M@units="°C" \
        -remapbil,../../inputs/CESM_atm_grid.txt \
         -cat /net/thermo/atmosdyn/era5/cdf/2021/${month}/N2021${month}??_?? \
        ${output_path}/tmp_2021_${month}_temperature_ERA5.nc
        #Z500
        # instantaneeous ERA5 data
        cdo -b F32  -selname,Z -sellevel,50000\
        -remapbil,../../inputs/CESM_atm_grid.txt \
         -cat /net/thermo/atmosdyn/era5/cdf/2021/${month}/Z2021${month}??_?? \
        ${output_path}/tmp_2021_${month}_z500_ERA5.nc
        done
# merge all years
cdo mergetime ${output_path}/tmp_*_temperature_detrended_clim_ERA5.nc ${output_path}/temperature_detrended_clim_ERA5_2021.nc
cdo mergetime ${output_path}/tmp_*_temperature_ERA5.nc ${output_path}/temperature_ERA5_2021.nc
cdo mergetime ${output_path}/tmp_*_z500_ERA5.nc ${output_path}/z500_ERA5_2021.nc

# clean up
rm ${output_path}/tmp_*_temperature_detrended_clim_ERA5.nc
rm ${output_path}/tmp_*_temperature_ERA5.nc
rm ${output_path}/tmp_*_z500_ERA5.nc