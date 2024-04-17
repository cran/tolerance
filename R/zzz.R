.onAttach <- function(lib, pkg){
  info <- packageDescription("tolerance")
  packageStartupMessage(
    paste('tolerance package, version ', info$Version, ', Released ', info$Date,
          '\n', 'This package is based upon work supported by the ',
          'Chan Zuckerberg Initiative: Essential Open Source Software for ',
          'Science (Grant No. 2020-255193).\n', sep="")
 )
}


