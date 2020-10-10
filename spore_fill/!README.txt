2020-05-23 4pm CST
filled in spinach downy mildew spore trap dataset from Robins PhD


Columns
date - date
site - gonz = Gonzales, kcn = King City North, kcs = King City South, sol = Soledad, spen = Salinas (Spence farm)
mean_copy - the mean of the trapping dataset for that day, typically 2 traps, although sometimes just 1.
mean_copy_fill - linearly interpolated mean copy value
missing_data - whether or not that datapoint was interpolated or not
year - what year the data was from
month - what month (1 = January)
day - day of the month (not Julian calendar)
ihs_data - inverse hyperbolic sine of the mean_copy_fill, sort of like log(data) but without needing to 'modify' with the addition of a small unit
log_data - the log(data+1)
yday - Julian calendar date
site_year - the site and year together, makes it easier to separate datasets