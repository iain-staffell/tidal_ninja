##################################################################
#                                                                #
#  BSD 3-Clause License                                          #
#  Copyright (C) 2020-2023  Iain Staffell  <staffell@gmail.com>  #
#  All rights reserved.                                          #
#                                                                #
##################################################################
####  ###
####  ###
##   ##      TIDAL.NINJA MODEL
#####
##
#
#
#    this model performs the 'core' tidal.ninja calculations.
#    it provides a soft-link between the predict_tide fortran
#    code that's part of osu tidal prediction software (otps) 
#    and simple models for tidal stream and range operations.
#    it can simulate hourly capacity factors for tidal energy
#    production (or sub-hourly if preferred) spanning decades
#    and for any location in the world.
#
#    REQUIREMENTS:
#
#     - A copy of the OTPS software, downloaded from:
#		https://www.tpxo.net/otps
#
#     - OTPS needs to be compiled into an executable that
#       can be called from R's system() command
#
#     - A copy of the TPXO 10 Atlas (v2), downloaded from:
#		https://www.tpxo.net/global/tpxo10-atlas
#
#
#    MAIN USER FUNCTIONS:
#
#     - Feed simulations into TPXO:
#		generate_otps_dates = function(start_date, end_date, interval)
#		generate_otps_coords = function(lon, lat)
#
#     - Pull results out of TPXO:
#		read_otps_results = function(file_path, time_path, type)
#
#     - Simulate output from tidal stream:
#		speed_to_power = function(speed, curve)
#
#     - Simulate output from tidal range:
#		solve_0d_rk4 = function(r, input)



	# prereqs
	suppressPackageStartupMessages(library(terra))


#####
## ##  INTERFACING WITH TPXO TO EXTRACT RESOURCE DATA
#####

	# these functions could be eliminated if we instead do i/o directly between R ←→ Fortran

	# provide ymd_hms plain-text string for the start and end date (YYYY-MM-DD format)
	# note that end date should be the period *after* you end (e.g. 1st hour of next year)
	# provide the interval (e.g. 'hour', '5 min', '30 sec')
	# returns the string to write for otps
	generate_otps_dates = function(start_date, end_date, interval)
	{
		# generate sequence of dates
		dates = seq(from = as.POSIXct(start_date), to = as.POSIXct(end_date), by = interval)

		# format their components
		year = format(dates, '%Y')
		month = as.integer(format(dates, '%m'))
		day = as.integer(format(dates, '%d'))
		hour = as.integer(format(dates, '%H'))
		minute = as.integer(format(dates, '%M'))
		second = as.integer(format(dates, '%S'))

		# apply the OTPS formatting
		sprintf('  %4s   %3d   %3d   %3d   %3d   %3d', year, month, day, hour, minute, second)
	}


	# provide a single or a vector of coordinates
	# returns the string to write for otps
	generate_otps_coords = function(lon, lat)
	{
		paste(" ", lat, lon)
	}

	# read an otps results file and return in a useful format
	# pass the results file (u or z type)
	# pass the time input file which generated this result (bleh)
	# pass 'u' or 'z' to say what it is
	read_otps_results = function(file_path, time_path, type)
	{
		# Read the file, skip the headers
		lines = readLines(file_path)
		lines = tail(lines, -6)
		N = length(lines)

		# There was an error in the model run
		if (N == 1)
		{
			warning(paste(lines, '\n'))
			return(NULL)
		}

		# Read the time file to know how many results to expect
		ntim = readLines(time_path)
		ntim = length(ntim) - 1

		# Expected result column names
		if (type == 'u') {
			col_names = c('date', 'time', 'U', 'V', 'u', 'v', 'd')
		} else if (type == 'z') {
			col_names = c('date', 'time', 'z', 'd')
		} else {
			stop('DED\n')
		}

		# Initialise results
		results = list()

		i = 1

		# loop through all coordinate blocks
		while (i < N)
		{
			# check this first line is a coordinate
 			is_coord = grepl("^\\s*-?\\d+\\.\\d+\\s+-?\\d+\\.\\d+", lines[i])
			if (!is_coord) stop('FUCT\n')

			# check we are going to have results
			if (lines[i] %contains% '***** Site is out of model grid OR land *****')
			{
				i = i + 1
				next
			}

			# Extract and format the coordinates
			coords = gsub("^\\s+|\\s+$", "", lines[i])
			coords = gsub("\\s+", " ", coords)

			# read the corresponding block of data
			data = read.table(text = lines[(i+1):(i+ntim)], header = FALSE, stringsAsFactors = FALSE)

			# Assign column names
			colnames(data) = col_names

			# Convert 'date' column to POSIXct
			data$date = as.POSIXct(paste(data$date, data$time), format="%m.%d.%Y %H:%M:%S")
			data$time = NULL

			# Save into our list
			results[[coords]] = data

			i = i + 1 + ntim
		}

		# Return
		return(results)
	}



#####
## ##  UNDERPINNING PHYSICAL STUFF
#####

	# calculate gravity at a given latitude using the grs80 formula
	gravity_grs80 = function(lat)
	{
		lat = deg2rad(lat)

		9.7803267715 * (1
			+ 0.0052790414 * sin(lat)^2
			+ 0.0000232718 * sin(lat)^4
			+ 0.0000001262 * sin(lat)^6
			+ 0.0000000007 * sin(lat)^8
		)
	}


	# calculate average sea water density at a given coordinate
	sea_density_ssd = function(lon, lat)
	{

		# try to extract our exact point
		ssd = rast('inputs/average_sea_density.tif')
		d = extract(ssd, cbind(lon, lat), cells=FALSE)[[1]]

		# if we missed and it's on land...
		if (is.na(d))
		{
			# create a bounding box around the point (±1 degree)
			buffer_extent = ext(lon - 1, lon + 1, lat - 1, lat + 1)

			# calculate the mean of a cropped raster
			cssd = crop(ssd, buffer_extent)
			d = global(cssd, mean, na.rm = TRUE)[[1]]
		}

		# if we still missed, i just don't give a fuck
		if (is.na(d))
		{
			d = 1024.68
		}

		d
	}




#####
## ##  TIDAL STREAM MODEL
#####

	# use a simple lookup to convert speed (m/s) into power production (kW)
	speed_to_power = function(speed, curve)
	{
		# convert to a power curve index
		index = round(100 * speed + 1)

		# constrain this to within acceptable limits
		index[ index < 1 ] = 1
		index[ index > length(curve) ] = length(curve)

		# return
		curve[index]
	}




#####
## ##  TIDAL RANGE MODEL
#####

	# identify and interpolate the high and low tides
	interpolate_low_tide = function(x)
	{
		# identify indices where the value is a local maximum or minimum
		is_peak = function(i)
		{
			(x[i] < x[i - 1] && x[i] < x[i + 1])
		}
		
		# only test indices that have both neighbors
		i = 2:(length(x) - 1)
		
		# find indices that are peaks
		i = i[sapply(i, is_peak)]

		# assemble input, padding at the ends to keep the spline in check
		peaks = data.frame(
			i = c(1, i, length(x)),
			x = c(x[head(i,1)], x[i], x[tail(i,1)])
		)
		
		# interpolate them to the full length
		spline(peaks$i, peaks$x, xout=1:length(x))$y
	}

	interpolate_high_tide = function(x)
	{
		-1 * interpolate_low_tide(-x)
	}

	# calculate the water height above and below which we 
	# operate our turbines, from the user-defined op_range
	identify_operating_times = function(r, input)
	{
		low_tide = interpolate_low_tide(r$z)
		high_tide = interpolate_high_tide(r$z)

		# generate the thresholds for pumping water in and out
		pump_out = low_tide + (high_tide - low_tide) * input$op_range
		pump_in = high_tide - (high_tide - low_tide) * input$op_range

		# and thus which timesteps we operate in
		operating = rep(0, nrow(r))
		operating[ r$z < pump_out ] = -1
		operating[ r$z > pump_in ] = 1

		return(operating)
	}


	# solve the 0d model with single timesteps
	solve_0d_ts1 = function(r, input)
	{
		# pre-calculate the spatially-varying physical constants
		gravity = gravity_grs80(input$lat)
		sea_density = sea_density_ssd(input$lon, input$lat)

		# pre-calculate anything else
		timestep = difftime(r$date[2], r$date[1], units='secs') |> as.numeric()

		# identify when our turbines will operate
		operating = identify_operating_times(r, input)

		# initialise our variables
		H = Q = P = r$z*0


		# run through every timestep
		for (i in 2:length(H))
		{
			# copy across the last operating state, to check if we should keep going
			if (i > 2 & operating[i] == 0)
				operating[i] = operating[i-1]

			# if we aren't operating, nothing happens
			if (operating[i] == 0)
			{
				H[i] = H[i-1]
				next
			}

			# what are the previous heights
			z = r$z[i-1]
			h = H[i-1]

			# stop operating if we have gone too far
			if (operating[i] == -1 & h <= z)
			{
				H[i] = H[i-1]
				operating[i] = 0
				next
			}
			if (operating[i] == 1 & h >= z)
			{
				H[i] = H[i-1]
				operating[i] = 0
				next
			}

			# height differential
			head = abs(z-h)

			# what is the maximum flow rate available to us
			q = input$turbine_cd * input$turbine_area * sqrt(2 * gravity * head)

			# what is the power output
			p = q * (sea_density * gravity * head / 10^6)

			# deal with hitting our rated power
			if (p > input$turbine_power / input$turbine_eff)
			{
				# cap the power
				p = input$turbine_power / input$turbine_eff

				# recalculate the flow
				q = p / (sea_density * gravity * head / 10^6)
			}

			# reduce power output by the turbine efficiency
			p = p * input$turbine_eff

			# what is the new upstream height
			dh_dt = (q / input$water_area / 10^6) * timestep * sign(z-h)
			h = h + dh_dt

			# update
			H[i] = h
			Q[i] = q
			P[i] = p
		}

		# assign these into our results data.frame
		r$op = operating
		r$h = H
		r$q = Q
		r$p = P

		return(r)
	}

	# calculate the change in upstream water height
	calculate_dh_dt = function(parms, z, h)
	{
		# height differential
		head = abs(z-h)

		# what is the maximum flow rate available to us
		q = parms$turbine_cd * parms$turbine_area * sqrt(2 * parms$gravity * head)

			# what is the resulting power output
			p = q * (parms$sea_density * parms$gravity * head / 10^6)

			# recalculate if we exceed our rated power
			if (p > parms$p_max)
			{
				q = parms$p_max / (parms$sea_density * parms$gravity * head / 10^6)
			}

		# what is the new upstream height
		dh_dt = (q / parms$water_area / 10^6) * sign(z-h)
		return(dh_dt)
	}

	# solve the 0d model with a 4th order rungle-kutta
	solve_0d_rk4 = function(r, input)
	{
		# pre-calculate the spatially-varying physical constants
		input$gravity = gravity_grs80(input$lat)
		input$sea_density = sea_density_ssd(input$lon, input$lat)

		# pre-calculate anything else
		timestep = difftime(r$date[2], r$date[1], units='secs') |> as.numeric()
		input$p_max = input$turbine_power / input$turbine_eff

		# identify when our turbines will operate
		operating = identify_operating_times(r, input)

		# initialise our variables
		H = Q = P = r$z*0


		# run through every timestep
		for (i in 2:length(H))
		{
			# copy across the last operating state, to check if we should keep going
			if (i > 2 & operating[i] == 0)
				operating[i] = operating[i-1]
			
			# if we aren't operating, nothing happens
			if (!operating[i])
			{
				H[i] = H[i-1]
				next
			}

			# what are the tide heights
			z = r$z[i-1]
			z_new = r$z[i]

			# what are the upstream heights
			h = H[i-1]

			# calculate the change in upstream height
			k1 = calculate_dh_dt(input, z, h)
			k2 = calculate_dh_dt(input, mean(z, z_new), h + (timestep / 2) * k1)
			k3 = calculate_dh_dt(input, mean(z, z_new), h + (timestep / 2) * k2)
			k4 = calculate_dh_dt(input, z_new, h + timestep * k3)
			dh_dt = (timestep / 6) * (k1 + 2*k2 + 2*k3 + k4)

			# what are the new heights
			z_new = r$z[i]
			h_new = h + dh_dt


			# stop operating if we have gone too far
			if (operating[i] == -1 & h_new <= z_new)
			{
				H[i] = H[i-1]
				operating[i] = 0
				next
			}
			if (operating[i] == 1 & h_new >= z_new)
			{
				H[i] = H[i-1]
				operating[i] = 0
				next
			}

			# back calculate the flow
			q = dh_dt * input$water_area * 10^6 / timestep / sign(z-h)

			# and the power output
			head = c( abs(z-h), abs(z_new-h_new) ) |> mean()
			p = q * (input$sea_density * input$gravity * head / 10^6)


			# recalculate if we exceed our rated power
			# we do this inside calculate dh_dt which checks sub-timesteps
			# but this catches edge cases which fall through
			if (p > input$p_max)
			{
				# cap the power
				p = input$p_max

				# recalculate the flow
				q = input$p_max / (input$sea_density * input$gravity * head / 10^6)

				# recalculate the upstream height
				dh_dt = (q / input$water_area / 10^6) * timestep * sign(z-h)
				h_new = h + dh_dt
			}


			# Store results
			H[i] = h_new
			Q[i] = q
			P[i] = p * input$turbine_eff
		}

		# assign these into our results data.frame
		r$op = operating
		r$h = H
		r$q = Q
		r$p = P

		return(r)
	}
