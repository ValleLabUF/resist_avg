
extract.covars.internal = function(dat, layers, state.col, which_cat, dyn_names, ind, p) {
  ## dat = data frame containing at least the id, coordinates (x,y), date-time (date), and
  ##      step length (step)
  ## layers = a raster object (Raster, RasterStack, RasterBrick) object containing environ covars
  ## state.col = character. The name of the column that contains behavioral states w/in
  ##             dat (if present)
  ## which_cat = vector of names or numeric positions of discrete raster layers; NULL by default 
  ## dyn_names = vector of names dynamic raster layers (in same order as layers); NULL by default
  ## ind = character/integer. The name or column position of the indicator column of dat to be
  ##       matched w/ names of a dynamic raster
  ## p = a stored 'progressr' object for creating progress bar
  
  
    #Subset and prep data
    tmp<- dat %>% 
      # dplyr::filter(id == ind[i]) %>% 
      dplyr::mutate(dt = difftime(date, dplyr::lag(date, 1), units = "secs")) %>% 
      dplyr::mutate_at("dt", {. %>% 
          as.numeric() %>%
          round()})
      tmp$dt<- c(purrr::discard(tmp$dt, is.na), NA)
      
      if (!is.null(dyn_names) & !is.factor(tmp[,ind])) stop("The `ind` column must be a factor.")
      
      
    
    #Identify levels of categorical layer (if available)
    if (!is.null(which_cat)) lev<- layers[[which_cat]]@data@attributes[[1]][,1]
    
    extr.covar<- data.frame()
    
    #Extract values from each line segment
    for (j in 2:nrow(tmp)) { 
      # print(j)
      segment<- tmp[(j-1):j, c("x","y")] %>%
        as.matrix() %>% 
        st_linestring() %>% 
        st_sfc(crs = projection(layers)) %>% 
        st_sf()
      
      tmp1<- raster::extract(layers, segment, along = TRUE, cellnumbers = FALSE) %>% 
        purrr::map(., ~matrix(., ncol = nlayers(layers))) %>%
        purrr::pluck(1) %>% 
        data.frame() %>% 
        purrr::set_names(names(layers))
      
      #subset to only include time-matched vars (by some indicator variable)
      cond<- tmp[j-1, ind]
      cond2<- levels(cond)[which(cond != levels(cond))]
      tmp1<- tmp1[,!stringr::str_detect(names(tmp1), paste(cond2, collapse="|")), drop=F]
      
      ind1<- stringr::str_which(names(tmp1), as.character(cond))
      names(tmp1)[ind1]<- dyn_names
      
      
      #calculate segment means if continuous and proportions spent in each class if categorical
      if (is.null(which_cat)) {
        covar.means<- data.frame(t(colMeans(tmp1)))
      } else {
        covar.means<- data.frame(t(colMeans(tmp1[,names(tmp1) != which_cat, drop = FALSE])))
        cat<- factor(tmp1[,which_cat], levels = lev)
        cat<- data.frame(unclass(t(table(cat)/length(cat))))
        
        covar.means<- cbind(covar.means, cat)  #Merge continuous and categorical vars
      }
      
      
      tmp2<- cbind(n = nrow(tmp1), dist = tmp$step[j-1], covar.means) %>% 
        dplyr::mutate(dt = as.numeric(tmp$dt[j-1]), id = unique(dat$id), date = tmp$date[j-1],
               state = ifelse(!is.null(state.col), tmp[j-1,state.col], NA)) #%>% 
        # dplyr::select(-cell)
      # tmp2<- tmp2[,!apply(is.na(tmp2), 2, any)]
      
      extr.covar<- rbind(extr.covar, tmp2)
    }
    
    p()  #plot progress bar
    extr.covar
}

#----------------------------
extract.covars = function(data, layers, state.col = NULL, which_cat = NULL, dyn_names = NULL,
                          ind) {
  ## data must be a data frame with "id" column, coords labeled "x" and "y" and datetime as POSIXct labeled "date"; optionally can have column that specifies behavioral state
  
  dat.list<- bayesmove::df_to_list(data, "id")
  
  progressr::with_progress({
    #set up progress bar
    p<- progressr::progressor(steps = length(dat.list))
    
    tictoc::tic()
    path<- furrr::future_map(dat.list,
                             ~extract.covars.internal(dat = .x, layers = layers,
                                                      state.col = state.col,
                                                      which_cat = which_cat,
                                                      dyn_names = dyn_names, ind = ind, p = p),
                             .options = furrr_options(seed = TRUE))
    tictoc::toc()
  })
  
  
  
  path<- dplyr::bind_rows(path)
  
  path
}
