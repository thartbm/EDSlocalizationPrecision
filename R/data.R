

downloadOSFdata <- function(demographics=TRUE, overwrite=TRUE, removezip=TRUE, getprocessed=TRUE, getraw=FALSE) {
  
  # if (overwrite) {
  #   conflicts = 'overwrite'
  # } else {
  #   conflicts = 'skip'
  # }
  
  filenames <- c()
  
  if (getraw) {
    filenames <- c(filenames, 'rawdata.zip')
  }
  
  if (getprocessed) {
    
    filenames <- c(filenames, 'processed.zip')
    
  }
  
  if (demographics) {
    filenames <- c(filenames, c('participants.csv', 'participants_files.csv'))
  }
  
  OSFnode <- osfr::osf_retrieve_node("t2jrs")
  
  # get a list of files for the year and semester that is requested:
  files <- osfr::osf_ls_files(OSFnode, path='data/', n_max=30)
  
  for (filename in filenames) {
    
    cat(sprintf('making sure we have: %s\n',filename))
    
    # find which line corresponds to the file:
    idx <- which(files$name == filename)
    
    # check that the file exists on OSF, and is unique:
    # if not: skip to next file
    if (length(idx) != 1) {
      next
    }
    
    # download the file:
    if (!file.exists(sprintf('data/%s',files$name[idx])) | overwrite) {
      osfr::osf_download(x = files[idx,], 
                         path = sprintf('data/%s', filename), 
                         overwrite = overwrite)
    }
    
    # check if it is a zip file:
    if (grepl('\\.zip$', filename)) {
      
      # then unzip it there:
      unzip(sprintf('data/%s',filename), exdir='data/')
      
    }
    
  }
  
  if (removezip) {
    
    for (filename in filenames) {
    
      # check if it is a zip file:
      if (grepl('\\.zip$', filename)) {
        
        # and remove the zip file, if that is wanted:
        file.remove(sprintf('data/%s',filename))
      }
      
    }
    
  }
  
}
