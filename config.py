# csv_check_work_tapes.py
extra_ground_stations = []
station_order = []
# how many cumulative records to use as a test
# it's not clear whether smaller or larger number gives better results 
# larger requires more adjustments, but may result in fewer deletions
cumsum_test = 10
# time to let the loops run for 
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