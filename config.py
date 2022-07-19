# csv_check_work_tapes.py
extra_ground_stations = []
station_order = []
cumsum_test = 10  # not set here, but loops through in the code 
# from 10, 20, 90 and 180 - used in check_compliance(), which tries to fix
# any data which don't fit 
cumsum_final_test = 180 # minimun consecutive records to be exported as a 
# block - reducing the number (for example to 10) can reduce the number
# of dropped files
# lower lets more records through, higher has fewer segments 10 to 180 seem
# reasonalble
# time to let the loops run for 
low_tolerance=0.6009 # used in calculate_gaps() to determine whether 
# it is a single gap between the current frame and the previous one 
high_tolerance=0.6071# used in calculate_gaps() to determine whether 
# it is a single gap between the current frame and the previous one 
lower_tolerance=0.5538 # used in the final calculation for 
# whether a record is good, and therefore included in the output 
higher_tolerance=0.6538 # used in the final calculation for 
# whether a record is good, and therefore included in the output 
low_single_gap=0.6009 # controls behaviour of calculate_gaps() - which 
# helps check_compliance() check the data and remove any that don't fit
# correct data generally fall within a narrow range 
high_single_gap=0.6071 # controls behaviour of calculate_gaps() - which 
# helps check_compliance() check the data and remove any that don't fit
time_test = 5*60 # in seconds
# use this setting to read the data in and make some checks
initial=False
last_station=None
last_timestamp=None
last_orig_idx=None


# csv_join_work_tapes.py
combine_ground_stations=True
clean_spikes=True
view_corrected_traces = False
fix_clock_error=True
fix_jump_error=True
exclude_masked_sections=True
potentially_suspect = 0
rejected = 0 
valid_end_frame = 0
valid_end_timestamp = 0
valid_end_time_index = 0
valid_end_delta4 = 0
clock_end_frame = 0
clock_end_timestamp = 0
clock_end_time_index = 0
clock_end_delta4 = 0
# this gets used to check that the new day starts on the right
# station and right time
# if the list doesn't apply for the previous day, no action is taken, 
# which means that the list doesn't need resetting 
last_frame_S12=None
last_timestamp_S12=None
last_ground_station_S12=None
time_divergence_S12=None
mean_S12_MH1=0
mean_S12_MH2=0
mean_S12_MHZ=0
mean_S12_SHZ=0

last_frame_S14=None
last_timestamp_S14=None
last_ground_station_S14=None
time_divergence_S14=None
mean_S14_MH1=0
mean_S14_MH2=0
mean_S14_MHZ=0
mean_S14_SHZ=0

last_frame_S15=None
last_timestamp_S15=None
last_ground_station_S15=None
time_divergence_S15=None
mean_S15_MH1=0
mean_S15_MH2=0
mean_S15_MHZ=0
mean_S15_SHZ=0

last_frame_S16=None
last_timestamp_S16=None
last_ground_station_S16=None
time_divergence_S12=None
mean_S16_MH1=0
mean_S16_MH2=0
mean_S16_MHZ=0
mean_S16_SHZ=0