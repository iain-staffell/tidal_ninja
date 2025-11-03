
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

		# physical parameters of your barrage
		water_area = 570,
		turbine_power = 8640,
		turbine_area = 216 * pi * (9/2)^2,

		# system operating and technical parameters
		op_range = 0.20,
		turbine_cd = 1.00,
		turbine_eff = 0.85
	
	)



	# generate the otps timesteps for the simulation
	times = generate_otps_dates(input$start_date, input$end_date, input$interval)
	writeLines(times, 'setup.time.txt')

	# generate the timesteps as an R vector (for plotting)
	times = seq(as.POSIXct(input$start), as.POSIXct(input$end_date), input$interval)


	# generate the otps coordinates for the simulation
	coords = generate_otps_coords(input$lon, input$lat)
	writeLines(coords, 'setup.coords.txt')


	# run the otps software
	shell('predict_tide -tsetup.time.txt < setup.z.txt > NUL 2>&1')

	# extract results from it
	r = read_otps_results('z.txt', 'setup.time.txt', 'z')


	# calculate operating times and power output
	p = solve_0d_rk4(r[[1]], input)



	# plot
	plot(times, p, type='l')
	