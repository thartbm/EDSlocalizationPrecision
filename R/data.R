
getPreprocessedData <- function(check=TRUE) {
  
  cat('Checking if all data is available locally.\n')
  
  data_files <- c('participants.csv'          = 'https://osf.io/6xpag/download',
                  'participants_files.csv'    = 'https://osf.io/xs6g5/download',
                  
                  'sEDS_curves.csv'           = 'https://osf.io/tv46u/download',
                  'sEDS_loc_AOV.csv'          = 'https://osf.io/3bxnh/download',
                  'sEDS_localization.csv'     = 'https://osf.io/6u9zg/download',
                  'sEDS_localization_tCI.csv' = 'https://osf.io/qzyut/download',
                  'sEDS_localization_var.csv' = 'https://osf.io/skv7u/download',
                  'sEDS_nocursors.csv'        = 'https://osf.io/shk29/download',
                  
                  'zEDS_curves.csv'           = 'https://osf.io/tbj6s/download',
                  'zEDS_loc_AOV.csv'          = 'https://osf.io/4rcdt/download',
                  'zEDS_localization.csv'     = 'https://osf.io/er97w/download',
                  'zEDS_localization_tCI.csv' = 'https://osf.io/muygq/download',
                  'zEDS_localization_var.csv' = 'https://osf.io/9tcp8/download',
                  'zEDS_nocursors.csv'        = 'https://osf.io/wgnfq/download'
    
  )
  
  for (filename in names(data_files)) {
    
    folderfilename <- sprintf('data/%s',filename)
    
    if (!check | !file.exists(folderfilename)) {
      
      url = as.character(data_files[filename])
      
      cat(sprintf("Downloading: '%s' from '%s'\n", filename, url))
      
      df <- read.csv(url(url),stringsAsFactors=FALSE)
      
      write.csv(df,folderfilename,row.names=FALSE,quote=FALSE)
      
    }
    
  }
  
}