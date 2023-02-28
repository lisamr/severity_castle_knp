# functions for importing data from the https://ftp server
require(RCurl)
require(XML)
require(glue)

f_check_slash <- function(url){
  if (substr(url, nchar(url), nchar(url)) != "/") {
    url <- glue('{url}/')
  }
  return(url)
}

# get subfolders from the URL. remove any that don't have fire data. 
f_get_subfolders <- function(url) {
  # check: does the URL end in /?
  url <- f_check_slash(url)
  
  # Get the contents of the directory
  html <- getURLContent(url)
  html_parsed <- XML::getHTMLLinks(html)
  
  #Use grepl to match only the element with all numbers
  subfolders <- html_parsed[grepl('^\\d+.?$', html_parsed)]
  
  subfolder_paths <- paste0(url, subfolders)
  
  return(subfolder_paths)
}


# find the gdb files in a given subfolder. 
f_get_file_paths <- function(subfolder, suffix, exclude = NA){
  # check: does the URL end in /?
  subfolder <- f_check_slash(subfolder)
  
  files <- getURL(subfolder) %>% XML::getHTMLLinks()
  pattern <- glue('(?<!{exclude})\\{suffix}$')
  gdbfile <- grep(pattern, files, value = T, perl = T)
  
  if(length(gdbfile) == 0) return(character(0)) else {
    return(glue("{subfolder}{gdbfile}"))
  }
}

# do it recursively when gdb file isn't immediately there. 
f_get_file_paths_r <- function(subfolder, suffix, check_folders, exclude = NA){
  
  # check: does the URL end in /?
  subfolder <- f_check_slash(subfolder)
  check_folders <- map_vec(check_folders, f_check_slash)
  
  files <- getURL(subfolder) %>% XML::getHTMLLinks()
  pattern <- glue('(?<!{exclude})\\{suffix}$')
  gdbfile <- grep(pattern, files, value = T, perl = T)
  
  # if gdbfile doesn't exist, check if it's within a designated subfolder
  if (length(gdbfile) == 0) {
    subfolders <- paste0(subfolder, files[files %in% check_folders])
    res <- map(subfolders, f_get_file_paths, suffix)
    return(res)
  }else{
    return(glue(subfolder, gdbfile))
  }
}

# unzip any kind of file
f_unzip <- function(local_file){
  subdirectory <- tools::file_path_sans_ext(local_file)
  if (!dir_exists(subdirectory)) dir_create(subdirectory, recurse = T)
  unzip(zipfile = local_file, exdir = subdirectory)
  try(fs::file_delete(local_file)) # wrap in try so it doesn't stop function in case can't remove
}

# download gbd file and optionally turn into an sf object
f_download <- function(url, dwn_dir, unzip = F){
  # download file
  local_file <- glue('{dwn_dir}/{basename(url)}')
  download.file(url, destfile = local_file, mode = "wb")
  
  # unzip and delete optionally
  if (unzip) f_unzip(local_file)
  
}


