library("tidyverse")
library("reshape2")
library('parallel')
library("readxl")
library('mlr3verse')
library('latex2exp')
library("mlr3extralearners")

# functions ========
#### prerequisite materials ####
pearson_heat <- function(df, corm = as.data.frame(-1)){
  defaultW <- getOption("warn")
  options(warn = -1) 
  if(corm[1, 1] == -1){
    corm = cor(df)
  }
  options(warn = defaultW)
  res <- reshape2::melt(corm) %>%
    ggplot(aes(Var1,Var2,fill = value)) + 
    geom_tile(color = "black",alpha = 1, 
              show.legend = T) + 
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 0.5,
                                     vjust = 0.5)) + 
    scale_fill_gradient2(low = 'black', high = 'white', mid = "red", midpoint = 0,guide = "colorbar", limit = c(-1,1)) + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) + 
    xlab(NULL) + 
    guides(fill = guide_colorbar(title = NULL)) +
    ylab(NULL)
  return(res)
}

# feature selecting by pearson correlation coefficient
fea_slc_bycor <- function(df, corr = 0.8){
  # df = features
  # corr = 0.7
  corm <- cor(df)
  name_of_features <- colnames(df)
  name_of_features_d <- name_of_features
  origin_fea_length <- length(name_of_features)
  for(q in 1:(origin_fea_length - 1)){
    # q = 1
    fea_t <- name_of_features_d[q]
    other_fea_t <- name_of_features_d[(q+1):length(name_of_features_d)]
    de_fea <- names(corm[fea_t, other_fea_t][abs(corm[fea_t, other_fea_t]) >= corr])
    name_of_features <- name_of_features[!(name_of_features %in% de_fea)]
  }
  res <- df[, colnames(df) %in% name_of_features]
  return(res)
}

# write filter to disk
create_flt <- function(f_task,f_fname, f_fp = NULL){
  # f_fp = NULL
  filter <- flt(f_fname)
  if(!is.null(f_fp)){
    filter$param_set$values <- f_fp
  }
  filter$calculate(task)
  flt_rank <- as.data.table(filter)
  system('mkdir results/filter_rank/')
  filename <- paste0(getwd(), "/results/filter_rank/", f_fname, ".csv")
  write_csv(as.data.frame(flt_rank), file = filename)
  return(flt_rank)
}

# function to write pridictioni file
to_write <- function(df, path){
  # df = st_t
  ind_seq <- c(seq(NN + 1, nrow(df), 500000), nrow(df))
  ind_seq_l <- lapply(1:(length(ind_seq) - 1), function(q){
    # q = length(ind_seq) - 1
    if(q != 1){
      return((ind_seq[q] + 1):ind_seq[q+1])
    }else{
      ind_seq[q]:ind_seq[q+1]
    }
  })
  pred <- lapply(1:length(ind_seq_l), function(q){
    # q = 1
    prediction2 <- learner$predict(task, row_ids = ind_seq_l[[q]])
    data.frame(prediction2$response)
  })
  preds <- do.call(rbind, pred)
  to_w <- st2[(NN+1): nrow(st2),]
  to_w[, ncol(to_w)] <- preds
  write_csv(to_w, file = path)
}

### scientific plot =========
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  # l = '0e+00'
  l <- gsub("0e\\+00","0",l)
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # l <- str_c('  ', l, collapse = '')
  # return this as an expression
  parse(text=l)
}

getname <- function(f_c){
  # f_c = 'CsCsCsCsGeBrBrBrBrII'
  # f_c = 'NaNaKCsGeFFClClClBr'
  sp1 <- strsplit(f_c, split = '')[[1]]
  nameVector <- rep('0', 11)
  id1 = 1
  for (s in 1:length(sp1)) {
    if(s == length(sp1)){
      nameVector[id1] <- sp1[s]
      id1 = 0
    }else if(sp1[s] %in% LETTERS) {
      if(sp1[s + 1] %in% letters){
        nameVector[id1] <- str_c(sp1[s], sp1[s + 1], collapse = '')
        id1 = id1 + 1
      }else{
        nameVector[id1] <- sp1[s]
        id1 = id1 + 1
      }
    }
  }
  names(nameVector) <- c('A1', 'A2', 'A3', 'A4', 'B', 
                         'X1', 'X2', 'X3', 'X4', 'X5', 'X6')
  nameVector <- nameVector[c('A1', 'A2', 'A3', 'A4', 'B', 'X1', 'X2', 'X3', 'X4', 'X5', 'X6')]
  return(nameVector)
}
# getname('NaNaKCsGeFFClClClBr')

getinfo <- function(path){
  # path = dir('rawdata/finish-update-0D/', full.names = T)[1]
  if("t.outmol" %in% dir(path)){
    outmol <- read_lines(str_c(path, '/t.outmol', collapse = ''))  
  }else{
    return(NA)
  }
  
  if(!("Message: DMol3 job finished successfully" %in% outmol)){
    return(NA) 
  }
  
  lind <- max(1, which("Message: DMol3 job finished successfully" == outmol) - 2000)
  out1 <- outmol[lind: which("Message: DMol3 job finished successfully" == outmol)]
  out2 <- out1
  out3 <- out1
  
  if(length(out2[str_detect(out2, " Yes ")]) == 0){
    return(NA)
  }
  
  #  detect yes or no
  while(length(out3[str_detect(out3, " Yes | No ")]) > 3){
    out3 <- out3[-1]
  }
  
  # extract
  while(length(out2[str_detect(out2, "LUMO")]) > 1){
    out2 <- out2[-1]
  }
  
  if(length(out3[str_detect(out3, " Yes ")]) == 3){
    homo <- as.numeric(gsub('*eV',  '', tail(strsplit(out2[which(str_detect(out2, "LUMO")) - 4], split = ' ')[[1]], 1)))
    lumo <- as.numeric(gsub('*eV',  '', tail(strsplit(out2[which(str_detect(out2, "LUMO")) - 3], split = ' ')[[1]], 1)))
    return(list(homo, lumo))
  }else{
    res1 <- sapply(1:length(out3[str_detect(out3, " No ")]), function(q){
      # q = 1
      t1 <- as.double(strsplit(out3[str_detect(out3, " No ")][q], split = "\\|")[[1]])
      t2 <- strsplit(as.character(formatC(t1[!is.na(t1)], format = "e")), split = "e")
      if(as.numeric(t2[[1]][2]) <= as.numeric(t2[[2]][2])){
        return(TRUE)
      }else{
        return(FALSE)
      }
    })
    
    if(sum(res1) == length(res1)){
      homo <- as.numeric(gsub('*eV',  '', tail(strsplit(out2[which(str_detect(out2, "LUMO")) - 4], split = ' ')[[1]], 1)))
      lumo <- as.numeric(gsub('*eV',  '', tail(strsplit(out2[which(str_detect(out2, "LUMO")) - 3], split = ' ')[[1]], 1)))
      return(list(homo, lumo))
    }else{
      return(NA)
    }
  }
}


get_solvent_information <- function(path){
  # path = dir('rawdata/finish-update-0D/', full.names = T)[1]
  inpuotFile <- read_lines(dir(path, full.names = T)[str_detect(dir(path), '.input')])
  sp1 <- str_split(inpuotFile[str_detect(inpuotFile, 'COSMO_Dielectric')], ' ')[[1]]
  as.numeric(sp1[sp1 != ''][2])
}


### scientific plot =========
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  # l = '0e+00'
  l <- gsub("0e\\+00","0",l)
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # l <- str_c('  ', l, collapse = '')
  # return this as an expression
  parse(text=l)
}

### prediction plot ==========
predplot <- function(prediction, title, scale = 10){
  # measure <- msr("regr.mse")
  # mse <- format(prediction$score(measure), digits = 2)
  measure <- msr("regr.mae")
  mae <- format(prediction$score(measure), digits = 2)
  # measure <- msr("regr.rmse")
  # rmse <- format(prediction$score(measure), digits = 2)
  preds <- prediction$response
  actual <- prediction$truth
  pr <- round(cor(cbind(preds, actual))[1,2],2)
  rss <- sum((preds - actual) ^ 2)  ## residual sum of squares
  tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
  rsq <- round(1 - rss/tss, 2)
  res1 <- autoplot(prediction) +
    theme_bw() + 
    scale_x_continuous(n.breaks = 5, labels=fancy_scientific) +
    scale_y_continuous(n.breaks = 5, labels=fancy_scientific) +
    labs(title = title) +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(hjust = .4),
          plot.margin = margin(10, 30, 10, 10)) +
    annotate("text",
             x = (0.3 * (range(prediction$response)[2] - range(prediction$response)[1]) + range(prediction$response)[1]),
             y = c(0.9 * (range(prediction$truth)[2] - range(prediction$truth)[1]) + range(prediction$truth)[1],
                   0.9 * (range(prediction$truth)[2] - range(prediction$truth)[1]) + range(prediction$truth)[1] - range(prediction$truth)[2]/scale,
                   0.9 * (range(prediction$truth)[2] - range(prediction$truth)[1]) + range(prediction$truth)[1] - 2*range(prediction$truth)[2]/scale), 
             label = c(TeX(paste0("Pearson's r: ",pr)), 
                       fancy_scientific(paste0("MAE: ",mae)),
                       fancy_scientific(paste0("R^2: ",rsq))),
             size = 5)
  return(res1)
}

### prediction plot ==========
predplot_classif <- function(task, learner, prediction, type, title, subtitle){
  rr <- resample(
    task = task,
    learner = learner,
    resampling = rsmp("cv", folds = 10)
  )
  if(subtitle == ""){
    subtitle = NULL
  }
  res1 <- autoplot(rr, type = type) +
    theme_bw() + 
    scale_x_continuous(n.breaks = 5, labels=fancy_scientific) +
    scale_y_continuous(n.breaks = 5, labels=fancy_scientific) +
    labs(title = title, subtitle = subtitle) +
    theme(text = element_text(size = 10),
          axis.text.x = element_text(hjust = .4),
          plot.margin = margin(10, 30, 10, 10))
  return(res1)
}

### axis latex ===========
fancy_axis <- function(l) {
  temp1 <- sapply(l, function(q){
    if(str_length(q) > 1){
      res1 <- str_c('$', substr(q, 1, 1),
                    '_{',
                    substr(q, 2, str_length(q)),
                    '}', '$'
      )
      return(TeX(res1))
    }else{
      TeX(q)
    }
  })
  return((temp1))
}

# correlation with goal ========
cwg <- function(fea_new_andC, lim = c(-1, 1), hlab = ''){
  
  if(hlab != ''){
    hlab <- str_c(" with ", hlab)
  }
  
  cor_for_fea_sel<- cor(fea_new_andC)
  cor_for_fea_sel_df <- cor_for_fea_sel[-ncol(fea_new_andC), ncol(cor_for_fea_sel)] %>%
    as.data.frame()
  colnames(cor_for_fea_sel_df) <- c("val")
  cor_for_fea_sel_df_pic <-
    cor_for_fea_sel_df %>%
    mutate(Names = rownames(cor_for_fea_sel_df)) %>%
    arrange((val)) %>%
    ggplot() +
    geom_col(aes(y = val, x = factor(Names, levels = (Names))), fill = "yellow") +
    theme_bw() +
    coord_flip() +
    scale_y_continuous(limits = lim) +
    xlab("Features") +
    ylab(str_c("Pearson Coefficient", hlab)) +
    theme(axis.text.y = element_text(size = 15,
                                     angle = 0,
                                     hjust = 0.1,
                                     vjust = 0.5),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    scale_x_discrete(labels = fancy_axis)
  return(cor_for_fea_sel_df_pic)
}

# rank for task type ============

write_rank <- function(task, wd, type){
  tmp1 <- as.data.table(mlr_filters)
  tmp2 <- sapply(tmp1$task_type, function(q){
    type %in% q
  })
  
  for (keyid in tmp1$key[tmp2]) {
    if(!(keyid %in% c("performance","permutation"))){
      filter <- flt(keyid)
      filter$calculate(task)
      (flt_correlation <- as.data.table(filter))
      fileName <- str_c(wd, keyid, ".csv")
      write_csv(as.data.frame(flt_correlation), file = fileName)
    }
  }
}

# symbolic regression ================
SR <- function(df){
  # ipd = fea_new
  # opd = homod
  ipd <- df[, 1:(ncol(df) - 1)]
  opd <- df[, ncol(df)]
  tn <- tail(colnames(df), 1)
  # opd <- homod[, ncol(homod)]
  da <- df
  fn <- colnames(ipd)
  da <- as.data.frame(da)
  da1 <- as.data.frame(matrix(0, nr = nrow(da), nc = ncol(da)))
  for (s in 1:nrow(da)) {
    da1[s, ] <- unlist(da[s, ])
  }
  da <- da1
  da <- as.data.frame(scale(da))
  colnames(da) <- c(paste("feature",1:(ncol(da)-1),sep = "_"),"tgt")
  colnames(da1) <- c(paste("feature",1:(ncol(da1)-1),sep = "_"),"tgt")
  attach(da)
  ruleDef <- list(expr = grule(op(expr, expr), func(expr), var),
                  func = grule(log, sqrt),
                  op = grule(`+`, `-`, `*`, `/`, `^`),
                  op = grule(`+`, `-`, `*`, `/`),
                  var = grule(feature_1^n, feature_2^n, feature_3^n,
                              feature_4^n, feature_5^n, feature_6^n, feature_7^n, 
                              feature_8^n, feature_9^n, feature_10^n, feature_11^n, 
                              feature_12^n, feature_13^n, feature_14^n, 
                              feature_15^n, feature_16^n, feature_17^n, 
                              feature_18^n, feature_19^n, feature_20^n),
                  n = grule(1, 2, 3, 0.5)
  )
  grammarDef <- CreateGrammar(ruleDef)
  # print(grammarDef)
  SymRegFitFunc <- function(expr) { 
    result <- eval(expr)
    if (any(is.nan(result)) | any(is.na(result))){
      return(Inf)
    }
    if(is.nan((cor(cbind(unlist(result), unlist(opd)))[1,2]))){
      return(Inf)
    }
    if(is.na((cor(cbind(unlist(result), unlist(opd)))[1,2]))){
      return(Inf)
    }
    if((cor(cbind(unlist(result), unlist(opd)))[1,2]) <= 0){
      return(Inf)
    }
    return(1/abs(cor(result, opd)))
  }
  ge <- GrammaticalEvolution(grammarDef = grammarDef, 
                             evalFunc = SymRegFitFunc, 
                             terminationCost = 1.25,
                             monitorFunc = print,
                             optimizer = "ga",
                             max.depth = 4,
                             seqLen = 50,
                             iterations = 1e5,
                             popSize = 1e2,
                             mutationChance = 0.3,
                             plapply = mclapply)
  print(colnames(df))
}

# name extractor ==========
name_extractor <- function(df){
  prediction_set <- read_csv(file = '../A4BX6/MS文件生成/Results/prediction.csv')
  abbs <- read_csv(file = '../A4BX6/MS文件生成/../rawdata/Abbs.csv')
  pred_set_positive <- pred_set[pred_set$bandgap == 1,]
  for(s in 1:(ncol(pred_set_positive) - 2)){
    colnames(pred_set_positive)[s] <- abbs$Fullname[abbs$Abbreviation == colnames(pred_set_positive)[s]]
  }
}





