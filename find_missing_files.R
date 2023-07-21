# find out which files weren't uploaded

# What I did to generate a vector of file names from the Data folder
# all_files <- fs::dir_ls('data', recurse = T, type = 'file')
# write(all_files, 'data_file_names.txt', ncolumns = 1)

# read in the text file as a vector
all_files <- readLines('data_file_names.txt')
head(all_files)

# get the paths of all the data files uploaded to Teams
Teams_files <- fs::dir_ls('data', recurse = T, type = 'file')

# find which files are missing and tell Lisa
missing_files <- all_files[!all_files %in% Teams_files]
missing_files
