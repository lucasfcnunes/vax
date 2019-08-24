#NECESSARY:
#have the following folders 'output', 'utils', 'cache'

#LOAD NECESSARY PACKAGES
for (package in c('readstata13', 'foreach', 'bigstatsr')) {
    if (!require(package, character.only=T, quietly=T)) {
        install.packages(package, repos="http://cran.us.r-project.org")
        library(package, character.only=T)
    }
}

priorities <- c('USA', 'ROW', 'CHN', 'JPN', 'DEU', 'GBR', 'FRA', 'BRA', 'IND', 'ITA', 'RUS', 'CAN', 'AUS', 'KOR', 'ESP', 'MEX', 'IDN', 'NLD', 'TUR', 'CHE', 'TWN', 'SWE', 'POL', 'BEL', 'NOR', 'AUT', 'DNK', 'FIN', 'IRL', 'GRC', 'PRT', 'CZE', 'ROU', 'HUN', 'SVK', 'LUX', 'BGR', 'HRV', 'LTU', 'SVN', 'LVA', 'EST', 'CYP', 'MLT')
datapath <- 'data'
dd_str <- function(int) substr(sprintf("%04d", int), 3, 4) #double digit string
metadata <- readRDS('utils/metadata.rds')

country_count <- length(metadata$countries)
country_sectors <- length(metadata$sectors) #sectors in a region
countries <- seq(from=1, to=country_count, by=1)
country_slice <- function(country) ((country-1)*country_sectors+1):(country*country_sectors)
country_sector_slice <- function(country, sector){
  FUN <- function(sector) (country-1)*country_sectors+sector
  return(sapply(sector, FUN=FUN))
}
country_number <- function(country_name) match(country_name, metadata$countries)

selected_countries <- country_number(priorities[1:10])

gen_vax_info <- function(year, sector, choice, corrected_to=''){
  aes <- endsWith(choice, 'aes')
  ##
  corrected_to <- if (corrected_to=='') year else corrected_to
  dim_ <- rep(country_count, 2)
  dimnames_ <- rep(list(metadata$countries), 2)
  if (aes){
    dim_[[2]] <- 1
    dimnames_[2] <- list('aes')
  }
  
  result <- list(
    choice=choice,
    year=year,
    corrected_to=corrected_to,
    sector=sector,
    sector_description=metadata$sectors[sector],
    varname=paste0(choice, dd_str(year), dd_str(corrected_to), '_sector', dd_str(sector)),
    dim=dim_,
    dimnames=dimnames_,
    method='2019-08-24'
  )
  return(result)
}

custom_save <- function(obj, outputpath, fextension){
  f <- attributes(obj)
  outfile <- file.path(outputpath, paste0(f$varname, '.', fextension), fsep=.Platform$file.sep)
  message('Saving file to "', outfile, '"')
  if (tolower(fextension) == 'rds'){
    saveRDS(obj, file=outfile)
  }
  else if (tolower(fextension) == 'csv'){
    df <- data.frame(rownames(obj))
    colnames(df) <- c(attributes(obj)$varname)
    df <- data.frame(df, obj)
    write.table(df, file=outfile, sep=',', row.names=FALSE)
  }
  else message("Couldn't save file...")
}

vax <- function (years, choice, sector, pyp=FALSE, outputpath='cache', datapath='data'){  
  message('\nCHOICE: ', choice)
  
  dnames <- c('VAdiag', 'L', 'DDfin')
  if (choice %in% c('VAXC', 'VAXP')) {}#dnames <- c('L')
  else if (choice %in% c('VAXD', 'VAXDaes')){
    dnames <- c(dnames, 'A')
    I <- diag(country_sectors*country_count)
  }
  else stop('impossible choice ', choice, '...')
  
  for(year in years){
    message('\nYEAR: ', year)

    ##LOAD EXTERNAL DATA
    for(dname in dnames){
      message('\n:: ', dname)
      ##READ FROM DISK (only dta files!)
      fullpath <- file.path(datapath, paste0(dname, dd_str(year), if (pyp) 'PYP' else '', '.', 'dta'), fsep = .Platform$file.sep)
      ##TRANSFORM AND PRE-PROCESS
      message('Loading to RAM... (Reading \'', fullpath,'\' from HD)')
      aux <- read.dta13(fullpath)[-(1:3)]
      message('Read!')
      message('Transforming/Pre-processing to a R matrix...')
      aux <- as.matrix(aux)
      #aux <- apply(aux[, -(1:3)], 1, FUN=as.numeric)
      try(if(all(dim(aux) != rep(country_count*country_sectors, 2))) stop('wrong matrix dim...'))
      assign(dname, aux)
      message('Success!')
      message('Cleaning cache...')
      rm(aux)
    }

    ##INITIALIZE VARS
    message('VAdiag to VA1D', '')
    VA1D <- diag(VAdiag) #VAdiag #
    if (choice %in% c('VAXD', 'VAXDaes')) GDP <- sum(VAdiag%*%L%*%DDfin)
    #rm(VAdiag)
    #######------######
    
    info <- c(gen_vax_info(year=year, sector=sector, choice=choice, corrected_to=if (pyp) year-1 else year),
      list(computation=list(start=Sys.time(), delta=NULL)))
    result <- FBM(info$dim[[1]], info$dim[[2]]) #(diag(0, country_count)+NA)
    result[,] <- NA
                           
    message('Calculating \'', attributes(result)$varname, '\'...')
    # create progress bar
    pb <- txtProgressBar(min=0, max=country_count^2, style = 3)
    pb_total <- -1
                           
    #cl <- parallel::makeCluster(2)
    #doParallel::registerDoParallel(cl)
    tmp <- foreach(countryA=countries, .combine = c) %:%
      foreach(countryB=countries, .combine = c) %do% {
        #update progress bar
        setTxtProgressBar(pb, pb_total<-pb_total+1)
        
        slice_sector <- country_sector_slice(countryA, sector)
        sliceA <- country_slice(countryA)
        sliceB <- country_slice(countryB)

        assign(choice, NA)
        if ((choice %in% c('VAXD', 'VAXDaes')) & (countryA %in% selected_countries) & (countryB %in% selected_countries)){
          #aes == 'all except self'
          if (choice == 'VAXDaes'){
            #stop('not implemented!')
            if (countryA == countryB){
              VA1D_ <- VA1D + 0 #replicates #VAdiag_ <- VAdiag + 0 #replicates
              A_ <- A + 0 #replicates
              DDfin_ <- DDfin + 0 #replicates

              VA1D_[slice_sector] <- 0 #VAdiag_[slice_sector,] <- 0
              
              A_[slice_sector, -sliceB] <- 0
              L_ <- solve(I-A_)

              DDfin_[slice_sector, -sliceB] <- 0
              
              assign(choice, GDP - sum(VA1D_%*%L_%*%DDfin_))
              countryB <- 1 # has only one 'aggregated' column 'aes'...
            }
            else return(NULL)
          }
          else if (choice == 'VAXD'){
            VA1D_ <- VA1D + 0 #replicates #VAdiag_ <- VAdiag + 0 #replicates
            A_ <- A + 0 #replicates
            DDfin_ <- DDfin + 0 #replicates

            VA1D_[slice_sector] <- 0 #VAdiag_[slice_sector,] <- 0

            A_[slice_sector, sliceB] <- 0
            L_ <- solve(I-A_)

            DDfin_[slice_sector, sliceB] <- 0
            
            assign(choice, GDP - sum(VA1D_%*%L_%*%DDfin_))
          }
        }
        else if (choice == 'VAXP'){
          assign(choice, sum(VA1D[slice_sector]*(L[slice_sector, sliceB]%*%DDfin[sliceB,])))
        }
        else if (choice == 'VAXC'){
          assign(choice, sum(VA1D[slice_sector]*(L[slice_sector,]%*%DDfin[,sliceB])))
        }
        else return(NULL)

        result[countryA, countryB] <- get(choice)
        return(NULL)
    }
    #parallel::stopCluster(cl)
    result <- result[]
    #close progress bar
    close(pb)
    #cleaning result
    result[!is.finite(result)] <- 0
    #stats
    info$computation$delta <- (Sys.time() - info$computation$start)
    #attr
    attributes(result) <- c(attributes(result), info)
    #saving
    custom_save(result, outputpath=outputpath, fextension='rds')
    rm(result)
  }
}

corrected_vax <- function(choice, year, corrected_to, sector, dirpath='cache'){
  fextension <- 'rds'
  aux <- function(year, corrected_to){
    f <- gen_vax_info(year=year, sector=sector, choice=choice, corrected_to=corrected_to)
    return(readRDS(file.path(dirpath, paste0(f$varname, '.', 'rds'), fsep=.Platform$file.sep)))
  }
  result <- aux(year, year)
  for (i in seq(year+1, corrected_to)){
    result <- result*(aux(i,i)/aux(i,i-1))
  }
  attributes(result) <- c(
    attributes(result),
    gen_vax_info(year=year, sector=sector, choice=choice, corrected_to=corrected_to)
  )
  #cleaning result
  result[!is.finite(result)] <- 0
  #saving
  custom_save(result, outputpath=dirpath, fextension='rds')
  rm(result) #return(result)
}

choices <- c('VAXDaes')
sectors <- c(17, 39, 40)

tmp <- foreach(choice=choices, .combine = c) %:% foreach(sector=sectors, .combine = c) %do% {
  ##BASIC PROCESSING
  vax(years=2000:2014, choice=choice, sector=sector, pyp=FALSE)
  vax(years=2001:2014, choice=choice, sector=sector, pyp=TRUE)

  ##YEARLY CORRECTIONS
  corrected_vax(choice, 2000, 2001, sector=sector)
  corrected_vax(choice, 2000, 2002, sector=sector)
  corrected_vax(choice, 2000, 2014, sector=sector)
  NULL
}

##TO CSV
dirpath <- 'output'
fextension <- 'csv'
for (f in list.files('cache', full.names=TRUE, include.dirs=FALSE, pattern="*.rds", ignore.case=TRUE)){
  message('Reading \'', f, '\' to RAM...')
  m <- readRDS(f)
  custom_save(obj=m, outputpath=dirpath, fextension=fextension)
}
