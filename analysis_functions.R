# functions to do analysis and organize output 
lmfun<-function(data,y,buffer,x, ndvi_measure){
  # format function call 
  xname = paste(c(paste0("ndvi_",as.character(buffer), "_", as.character(ndvi_measure),"_iqr"),x) , collapse = " + ")
  formula1<-as.formula(paste(y,"~",xname))
  lm.fit<-do.call("lm",list(data=quote(data),formula1))
  return(lm.fit)
}

lmfull_processed<-function(data, y, buffers, x, sig_figs = 3, ndvi_measure){
  # run models and store as output
  model = vector(mode="list")
  for(i in 1:length(buffers)){
    model[[i]] = lmfun(data = data, y = y, buffer = buffers[i], x = x, ndvi_measure = ndvi_measure)
  }
  # format and store summary output
  results = data.frame(matrix(data = NA, nrow = 1, ncol = length(buffers)))
  for(i in 1:length(buffers)){
    main <- as.data.frame(cbind(coef(model[[i]]), confint(model[[i]])))[-1,]
    main <- round(main, sig_figs)
    results[,i]<-paste0(main$V1[1]," (", main$`2.5 %`[1], ", ", main$`97.5 %`[1], ")")
    colnames(results)[i] = paste0("NDVI IQR ", buffers[i], "m")
  }
  rownames(results) = c("ndvi")
  return(results)
}

lmfull<-function(data, y, buffers, x, sig_figs = 3, ndvi_measure){
  # run models and store as output
  model = vector(mode="list")
  for(i in 1:length(buffers)){
    model[[i]] = lmfun(data = data, y = y, buffer = buffers[i], x = x, ndvi_measure = ndvi_measure)
  }
  results = data.frame(matrix(data = NA, nrow = length(buffers), ncol = 3))
  for(i in 1:length(buffers)){
    main <- as.data.frame(cbind(coef(model[[i]]), confint(model[[i]])))[-1,]
    main <- round(main, sig_figs)
    results[i,] = main[1,]
  }
  colnames(results) <- c("Mean", "Lower_CI", "Upper_CI")
  return(list(model = model, results = results))
}
