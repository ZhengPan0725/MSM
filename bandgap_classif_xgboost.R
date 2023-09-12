rm(list = ls()) 
setwd("~/Desktop/ZhengPan/workspace/NUIST/John/胡文广/A4BX6/")  
source('JAP/code/dependencies.R')
cl <- mlr_learners$keys()[str_detect(mlr_learners$keys(), 'classif')]

#==== Utilities for this specific task ====
clean_up <- function(df, abbs){
  # df = st1
  df <- df %>%
    select(-contains('Symbol'), -contains('Phase'), -contains('Radioactive'), 
           -contains('Natural'), -contains('Metal'), -contains('Nonmetal'), 
           -contains('Metalloid'), -contains('Type'),-contains('Discoverer'),
           -contains('Group'), -contains('Year'), -contains('NumberofValence'))
  df <- df[!is.na(df$homo), ]
  colnames(df) <- abbs$Abbreviation
  df <- as.data.frame(df)
  st <- df[, -1]
  for(s in 1:ncol(st)){
    st[, s] <- as.numeric(st[, s])
  }
  return(st)
}
get_tune_expression <- function(learner){
  paras_to_tune <- learner$param_set$ids()[!(learner$param_set$class %in% c("Paramlgl", "ParamUty", "ParamFct"))]
  ps_all <- learner$param_set
  res1 <- sapply(1:length(paras_to_tune), function(q){
    # q = 1
    tdf <- as.data.frame(as.data.table(ps_all))[which(paras_to_tune[q] == learner$param_set$ids()), ]
    if(tdf$class == "ParamInt"){
      if(tdf$upper == "Inf"){
        if(tdf$lower == "-Inf"){
          rts <- str_c(tdf$id, " = p_int(lower = ", -1e4, ", upper = ", 1e4, ")")  
        }else{
          rts <- str_c(tdf$id, " = p_int(lower = ", tdf$lower, ", upper = ", 1e4, ")")
        }
      }else{
        if(tdf$lower == "-Inf"){
          rts <- str_c(tdf$id, " = p_int(lower = ", -1e4, ", upper = ", tdf$upper, ")")
        }else{
          rts <- str_c(tdf$id, " = p_int(lower = ", tdf$lower, ", upper = ", tdf$upper, ")")
        }
      }
    }else if(tdf$class == "ParamDbl"){
      if(tdf$upper == "Inf"){
        if(tdf$lower == "-Inf"){
          rts <- str_c(tdf$id, " = p_dbl(lower = ", -1e4, ", upper = ", 1e4, ")")  
        }else{
          rts <- str_c(tdf$id, " = p_dbl(lower = ", tdf$lower, ", upper = ", 1e4, ")")
        }
      }else{
        if(tdf$lower == "-Inf"){
          rts <- str_c(tdf$id, " = p_dbl(lower = ", -1e4, ", upper = ", tdf$upper, ")")
        }else{
          rts <- str_c(tdf$id, " = p_dbl(lower = ", tdf$lower, ", upper = ", tdf$upper, ")")
        }
      }
    }else{
      rts <- ''
    }
    return(rts)
  })
  res1 <- res1[res1 != '']
  for (tid in 1:length(res1)) {
    if(tid != length(res1)){
      res1[tid] <- str_c(res1[tid], ',', collapse = '')  
    }
  }
  res1 <- c("search_space = ps(",res1, ")")
  res1 <- str_c(res1, collapse = '')
  return(res1)
}
#==== Constants ===========
##==== not changable ====
st1 <- read_csv(file = 'rawdata/A4BX6train_and_test_fill1.csv')
abbs <- read_csv(file = 'rawdata/Abbs.csv')
st <- clean_up(df = st1, abbs = abbs)
features <- st[, -ncol(st):(-ncol(st) + 1)] 
fea_new <- fea_slc_bycor(features, corr = 0.95)
means <- apply(fea_new, 2, mean)
sds <- apply(fea_new, 2, sd)
fea_new <- scale(fea_new, center = means, scale = sds)
bandgapd <- cbind(fea_new, bandgap = st$lumo - st$homo) 
bandgapd <- as.data.frame(bandgapd)
bandgapd$bandgap <- as.factor(as.numeric(bandgapd$bandgap <= 3.4))
NN <- nrow(bandgapd)
##==== changable ====
patition_ratio <- 0.9
data_for_ML <- bandgapd
#==== Main work flow ===============
##==== temp constant ====
learner_to_use <- "classif.xgboost"
train_plot_title <- str_c(str_replace(learner_to_use, 
                                      pattern = 'regr.',
                                      replacement = ''),
                          ":train set", collapse = '')
test_plot_title <- str_c(str_replace(learner_to_use, 
                                     pattern = 'regr.', 
                                     replacement = ''),
                         ":test set", collapse = '')
##==== workflow ====
task <- TaskClassif$new(id = "data_for_ML", 
                     backend = data_for_ML, 
                     target = tail(colnames(data_for_ML), 1))
splits <- partition(task, ratio = patition_ratio)
learner <- lrn(learner_to_use,
               predict_type = 'prob', 
               nrounds = round(exp(4)),
               booster = 'dart',
               eta = 1, 
               gamma = exp(-2.5), 
               lambda = exp(1.535057), 
               alpha = exp(1.535057),
               subsample = 1.0,
               max_depth = round(exp(15)), 
               min_child_weight = 1, 
               colsample_bytree = 0.34,
               colsample_bylevel = 0.67, 
               rate_drop = 1, 
               skip_drop = 0.66667)
# set_threads(learner, n = 5)
learner$train(task = task, row_ids = splits$train)
predictions_train <- learner$predict(task, row_ids = splits$train)
predictions_test <- learner$predict(task, row_ids = splits$test)
autoplot(predictions_test)
predplot(prediction = predictions_train,
         title = learner_to_use,
         scale = 10)

instance = ti(
  task = task,
  learner = learner,
  resampling = rsmp("cv", folds = 10),
  measures = msr("classif.auc"),
  # terminator = trm("evals", n_evals = 200),
  terminator = trm("none"),
  search_space = lts("classif.xgboost.rbv2")
)
trn1 <- tnr("grid_search",resolution = 4, batch_size = 5) 
trn1$optimize(instance)
# learner$param_set$values <- instance$result_learner_param_vals


