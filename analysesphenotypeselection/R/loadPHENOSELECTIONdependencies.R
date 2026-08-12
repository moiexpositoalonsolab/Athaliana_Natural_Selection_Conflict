loadPHENOSELECTIONdependencies_<-function(){
  # general
  library(dplyr)
  library(tidyverse)
  # plotting
  library(ggplot2)
  library(ggpmisc)
  library(cowplot)
  theme_set(theme_cowplot())
  library(RColorBrewer)
  library(ggfortify)
  # custom?
  # library(genemaps)
  # analuyses
  library(caret)
}
loadPHENOSELECTIONdependencies<-function() suppressMessages(loadPHENOSELECTIONdependencies_())
