def generate_settings():	
	settings = {}
	settings['x_lower'] = 1
	settings['x_upper'] = 1000
	settings['y_lower'] = 1
	settings['y_upper'] = 1000000
	
	# Scales can be linear or log
	settings['x_scale'] = 'linear'
	settings['y_scale'] = 'log'
	settings['x_label'] = 'k-mer Coverage'
	settings['y_label'] = 'k-mer Count Frequency'
	return settings