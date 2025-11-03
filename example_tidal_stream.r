
	# load the tidal.ninja model
	setwd('M:/WORK/~TMP~/~TIDAL~/online_ninja/')
	source('tidal_ninja.r')

	# set our inputs for the simulation
	input = list(

		# coordinates
		lat = 51.37,
		lon = -3.09,

		# timesteps
		start_date = '2024-01-01 00:00:00',
		end_date = '2025-01-01 00:00:00',
		interval = '10 min',

		# power in MW
		turbine_power = 2,

		# turbine model (must exist within `power_curves`)
		turbine_curve = 'MCT SeaGen-S 2MW'

	)



	# load up your databaes of power curves
	power_curves = read_data('inputs/turbine_stream_power_kW_ti10.csv')


	# generate the otps timesteps for the simulation
	times = generate_otps_dates(input$start_date, input$end_date, input$interval)
	writeLines(times, 'setup.time.txt')

	# generate the timesteps as an R vector (for plotting)
	times = seq(as.POSIXct(input$start), as.POSIXct(input$end_date), input$interval)


	# generate the otps coordinates for the simulation
	coords = generate_otps_coords(input$lon, input$lat)
	writeLines(coords, 'setup.coords.txt')


	# run the otps software
	shell('predict_tide -tsetup.time.txt < setup.u.txt > NUL 2>&1')

	# extract results from it
	r = read_otps_results('u.txt', 'setup.time.txt', 'u')


	# calculate speed from u and v
	s = naz(sqrt(r$u^2 + r$v^2)) / 100

	# convert speed to power
	p = speed_to_power(s, power_curves[[input$turbine_curve]])



	# plot
	plot(times, p, type='l')
	