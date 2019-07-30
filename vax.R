
#NECESSARY:
#have the following folders 'output', 'utils', 'cache'

#LOAD NECESSARY PACKAGES
for (package in c('readstata13')) {
    if (!require(package, character.only=T, quietly=T)) {
        install.packages(package, repos="http://cran.us.r-project.org")
        library(package, character.only=T)
    }
}

datapath <- 'data'
dd_str <- function(int) substr(sprintf("%04d", int), 3, 4) #double digit string
metadata <- readRDS('utils/metadata.rds')

gen_vax_info <- function(year, sector, choice, corrected_to=''){
  corrected_to <- if (corrected_to=='') year else corrected_to
  result <- list(
    choice=choice,
    year=year,
    corrected_to=corrected_to,
    sector=sector,
    sector_description=metadata$sectors[sector],
    varname=paste0(choice, dd_str(year), dd_str(corrected_to), '_sector', dd_str(sector)),
    dimnames=rep(list(metadata$countries), 2)
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

vax <- function (years, choice, sector=17, pyp=FALSE, outputpath='cache', datapath='data'){  
  message('\nCHOICE: ', choice)
  
  country_count <- length(metadata$countries)
  country_sectors <- length(metadata$sectors) #sectors in a region
  
  if (choice %in% c('VAXC', 'VAXP')) dnames <- c('L')
  else if ('VAXD' == choice){
    dnames <- c('A')
    I <- diag(country_sectors)
  }
  else stop('impossible choice ', choice, '...')
  dnames <- c(dnames, c('VAdiag', 'DDfin'))
  
  countries <- seq(from=1, to=country_count, by=1)
  country_slice <- function(country) ((country-1)*country_sectors+1):(country*country_sectors)
  country_sector_line <- function(country, sector) (country-1)*country_sectors+sector
  
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
    rm(VAdiag)
    #######------######

    result <- (diag(0, country_count)+NA)
    attributes(result) <- c(
      attributes(result),
      gen_vax_info(year=year, sector=sector, choice=choice, corrected_to=if (pyp) year-1 else year),
      list(newmethod=TRUE, computation=list(start=Sys.time(), delta=NULL))
                   )
      
    message('Calculating \'', attributes(result)$varname, '\'...')
    # create progress bar
    pb <- txtProgressBar(min=0, max=country_count^2, style = 3)

    for (countryB in countries){
      #message('(index, jndex) = (', index, jndex, ')')
      slice_j <- country_slice(countryB)
      #message('"Zeroing" DDfin...')
      #else stop('not a possible choice')

      for (countryA in countries){
        index <- country_sector_line(countryA, sector)
        slice_i <- country_slice(countryA)
        if (choice == 'VAXD'){
          #L <- if (countryA != countryB) I else solve(I-A[slice_i,slice_j])
          #assign(choice, sum(VA1D[index]*(L[sector:sector,]%*%DDfin[slice_i,slice_j])))
          
          S <- if (countryA == countryB) solve(I-A[slice_i,slice_j])[sector:sector,]%*%DDfin[slice_i,slice_j] else DDfin[index:index,slice_j]
          assign(choice, sum(VA1D[index]*S))
          
        }
        else if (choice == 'VAXP'){
          assign(choice, sum(VA1D[index]*(L[index:index,slice_j]%*%DDfin[slice_j,])))
        }
        else if (choice == 'VAXC'){
          assign(choice, sum(VA1D[index]*(L[index:index,]%*%DDfin[,slice_j])))
        }          
        #update progress bar
        setTxtProgressBar(pb, getTxtProgressBar(pb)+1)

        result[countryA, countryB] <- get(choice)
        #message('Success!\n', names(result), '\n', paste(result))
      }
    }
    #close progress bar
    close(pb)
    #stats
    attributes(result)$computation$delta <- (Sys.time() - attributes(result)$computation$start)
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
  custom_save(result, outputpath=dirpath, fextension='rds')
  #return(result)
}

##BASIC PROCESSING
choice <- 'VAXP'
system.time(vax(years=2000:2014, choice=choice, sector=17, pyp=FALSE))
system.time(vax(years=2001:2014, choice=choice, sector=17, pyp=TRUE))

##YEARLY CORRECTIONS
corrected_vax(choice, 2000, 2001, sector=17)
corrected_vax(choice, 2000, 2002, sector=17)
corrected_vax(choice, 2000, 2014, sector=17)

##TO CSV
dirpath <- 'output'
fextension <- 'csv'
for (f in list.files('cache', full.names=TRUE, include.dirs=FALSE, pattern="*.rds", ignore.case=TRUE)){
  message('Reading \'', f, '\' to RAM...')
  m <- readRDS(f)
  custom_save(obj=m, outputpath=dirpath, fextension=fextension)
}
