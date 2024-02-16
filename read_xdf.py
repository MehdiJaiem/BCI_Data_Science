# %%
import pyxdf
# %%
# XDF file name
file = 'sub-P001_ses-S001_task-Default_run-001_eeg.xdf'

# Whether markers should be appended as a new channel to the Raw object
with_stim = False

# Read the XDF file
streams, header = pyxdf.load_xdf(file)

# Initialize lists for the streams
marker_streams = []
data_stream = []
