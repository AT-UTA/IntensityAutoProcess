# IntensityAutoProcess
# Must also have bfmatlab to read .nd2 files.
# Files must be .nd2 file type and have the same resolution size (ie. if 1024x1024 then all images must be 1024x1024)
# Can process any number of files with any number of z-planes but must have 3 channels.
# How to use:
# 1. Make a master folder that will hold all groups of interest to process the fluorescent intensity of.
# 2. Within the master folder, make folders for each group named appropriately, and place the .nd2 file for each group in their respective folder.
# 3. Make a variable for the path to the master folder.
# 4. Call on the script process_cell_imagev2(PATH);
# 5. Fill out the staining target for each channel.
# 6. Either use the made graphs or copy data and remake in the preferred graphing program.
