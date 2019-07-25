metadata <- list(
  countries=read.dta13("data\\VAdiag00.dta")[seq(1,2464,56),2],
  sectors=read.dta13("data\\VAdiag00.dta")[1:56,1]
  )