#' Developing the optimal predictive model for the dichotomous variables with machine learning algorithms
#'
#' A function can be used to develop the predictive model for dichotomous variables with four machine learning algorithms.
#'
#' @param train_data The training data with the 'ID' and 'Var' as the first two columns. Starting in the third column are the variables used to construct the model. 'Var' is the target predictor variable for constructing the model. 'Var' contains only Y or N.
#' @param list_train_vali_Data A list containing the training data and the other validation data. All the validation data have the same data form as the training data.
#' @param candidate_genes The candidate variables used for constructing the predictive model.
#' @param methods There are four algorithms for developing the predictive model including 'nb', 'svmRadialWeights', 'rf', 'kknn'.
#'   'nb': Naive Bayes algorithm.
#'   'svmRadialWeights': Support Vector Machine (SVM).
#'   'rf': Random Forest.
#'   'kknn': K-nearest Neighbors.
#' @param seed The seed you can set as any positive number, for example, 5201314.
#' @param cores_for_parallel The cores you can choose for parallel operation. The default is 12. The bigger the better if the configuration allows it.
#'
#' @return A list containing the predictive model, the AUC, the ROC, and the candidate variables, all of which are developed by each single algorithm.
#' @export
#'
#' @examples
ML.Dev.Pred.Category.Sigs2 <- function(train_data, # cohort data used for training, the colnames of which inlcuding ID, Var, and the other candidate genes。
                                     # Var 是用于构建预测模型的目标变量，Y/N，
                                     list_train_vali_Data, # cohort data used for training, the colnames of which inlcuding ID, Var, and the other candidate genes。
                                     # Var 是用于构建预测模型的目标变量，Y/N
                                     candidate_genes = NULL,
                                     methods = NULL, # c('nb','svmRadialWeights','rf','kknn')
                                     seed = 5201314, # 5201314
                                     cores_for_parallel = 12 #
) {
  ### === loading packages ===###
  message("---loading the packages ---")
  
  
  if (T) {
    library(stringr)
    library(gridExtra)
    library(future)
    library(sva)
    library(e1071)
    library(pROC)
    library(ROCit)
    library(caret)
    library(doParallel)
    library(dplyr)
  }
  
  
  ###### loading the function #######
  
  
  model.Dev <- function(training, method, sig) {
    training <- training[, colnames(training) %in% c("Var", sig)]
    # 4 models adopted in this study as followings:
    #' nb': naive bayes
    #' svmRadialWeights': Support Vector Machines with Class Weights
    #' rf': random forest
    #' kknn': k-Nearest Neighbors
    
    # Grid search for parameter tuning
    Grid <- list(
      nb = expand.grid(fL = c(0, 0.5, 1, 1.5, 2.0), usekernel = TRUE, adjust = c(0.5, 0.75, 1, 1.25, 1.5)),
      svmRadialWeights = expand.grid(sigma = c(0.0005, 0.001, 0.005, 0.01, 0.05), C = c(1, 3, 5, 10, 20), Weight = c(0.1, 0.5, 1, 2, 3, 5, 10)),
      rf = expand.grid(mtry = c(2, 42, 83, 124, 165, 205, 246, 287, 328, 369)),
      kknn = expand.grid(kmax = c(5, 7, 9, 11, 13), distance = 2, kernel = "optimal")
    )
    TuneLength <- list(
      nb = nrow(Grid[["nb"]]),
      svmRadialWeights = nrow(Grid[["svmRadialWeights"]]),
      rf = nrow(Grid[["rf"]]),
      kknn = nrow(Grid[["kknn"]])
    )
    
    ls_model <- lapply(method, function(m) {
      f <- 5 # f folds resampling
      r <- 10 # r repeats
      n <- f * r
      
      # sets random seeds for parallel running for each single resampling f-folds and r-repeats cross-validation
      seeds <- vector(mode = "list", length = n + 1)
      # the number of tuning parameter
      for (i in 1:n) seeds[[i]] <- sample.int(n = 1000, TuneLength[[m]])
      
      # for the last model
      seeds[[n + 1]] <- sample.int(1000, 1)
      
      
      ctrl <- trainControl(
        method = "repeatedcv",
        number = f, ## 5-folds cv
        summaryFunction = twoClassSummary, # Use AUC to pick the best model
        classProbs = TRUE,
        repeats = r, ## 10-repeats cv,
        seeds = seeds
      )
      
      
      
      model.tune <- train(Var ~ .,
                          data = training,
                          method = m,
                          metric = "ROC",
                          trControl = ctrl,
                          tuneGrid = Grid[[m]]
      )
      
      print(m)
      return(model.tune)
    })
    
    
    names(ls_model) <- method
    
    return(ls_model)
  }
  
  
  
  cal.model.auc <- function(res.by.model.Dev, cohort.for.cal, sig) {
    library(dplyr)
    
    rownames(cohort.for.cal) <- cohort.for.cal$ID
    validation <- cohort.for.cal[, colnames(cohort.for.cal) %in% c("Var", sig)]
    validation$Var <- factor(validation$Var, levels = c("N", "Y"))
    
    ls_model <- res.by.model.Dev
    models <- names(ls_model)
    auc <- lapply(1:length(models), function(i) {
      model.tune <- ls_model[[i]]
      prob <- predict(model.tune, validation[, -1], type = "prob")
      pre <- predict(model.tune, validation[, -1])
      test_set <- data.frame(obs = validation$Var, N = prob[, "N"], Y = prob[, "Y"], pred = pre)
      auc <- twoClassSummary(test_set, lev = levels(test_set$obs))
      return(auc)
    }) %>% base::do.call(rbind, .)
    
    rownames(auc) <- names(ls_model)
    
    return(auc)
  }
  
  
  
  cal.model.roc <- function(res.by.model.Dev, cohort.for.cal, sig) {
    library(dplyr)
    
    rownames(cohort.for.cal) <- cohort.for.cal$ID
    validation <- cohort.for.cal[, colnames(cohort.for.cal) %in% c("Var", sig)]
    validation$Var <- factor(validation$Var, levels = c("N", "Y"))
    
    ls_model <- res.by.model.Dev
    models <- names(ls_model)
    roc <- lapply(1:length(models), function(i) {
      prob <- predict(ls_model[[models[i]]], validation[, -1], type = "prob") #
      pre <- predict(ls_model[[models[i]]], validation[, -1]) #
      test_set <- data.frame(obs = validation$Var, N = prob[, "N"], Y = prob[, "Y"], pred = pre)
      roc <- ROCit::rocit(
        score = test_set$N,
        class = test_set$obs,
        negref = "Y"
      )
      return(roc)
    })
    
    names(roc) <- models
    
    
    return(roc)
  }
  
  
  message("---loading the function---")
  
  
  list_train_vali_Data <- lapply(list_train_vali_Data,function(x){
    colnames(x) = gsub('-','.',colnames(x))
    return(x)})
  colnames(train_data) <- gsub("-", ".", colnames(train_data))
  candidate_genes <- gsub("-", ".", candidate_genes)
  common_feature <- c("ID", "Var", candidate_genes)
  
  for (i in names(list_train_vali_Data)) {
    common_feature <- intersect(common_feature, colnames(list_train_vali_Data[[i]]))
  }
  
  ##### parameters check #####
  
  allowed_methods <- c("nb", "svmRadialWeights", "rf", "kknn")
  if (is.null(methods)) methods <- allowed_methods
  
  if (
    all(is.element(methods, allowed_methods))
  ) {
    ####### data preparation ######
    
    
    
    list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
      x[, common_feature]
    })
    
    list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
      x[, -c(1:2)] <- apply(x[, -c(1:2)], 2, as.numeric)
      rownames(x) <- x$ID
      return(x)
    })
    
    list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
      x[, c(1:2)] <- apply(x[, c(1:2)], 2, as.factor)
      return(x)
    })
    
    
    list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
      x <- x[!is.na(x$Var) & !is.na(x$Var), ]
      return(x)
    })
    
    # use the mean replace the NA
    list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
      x[, -c(1:2)] <- apply(x[, -c(1:2)], 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = T)
        return(x)
      })
      
      
      return(x)
    })
    
    train_data <- train_data[, common_feature]
    train_data[, -c(1:2)] <- apply(train_data[, -c(1:2)], 2, as.numeric)
    train_data[, c(1:2)] <- apply(train_data[, c(1:2)], 2, as.factor)
    rownames(train_data) <- train_data$ID
    
    est_dd <- as.data.frame(train_data)[, common_feature[-1]]
    pre_var <- common_feature[-c(1:2)]
    
    
    print(paste0("There existing ", length(candidate_genes), " genes in candidate genes"))
    
    
    print("Intersetion of the candidate genes and the colnames of the provided data")
    
    print(paste0("There existing ", length(pre_var), " genes in candidate genes, colnames of training data, colnames of validation data"))
    
    
    
    # parallel processing
    
    
    cl <- makePSOCKcluster(cores_for_parallel)
    registerDoParallel(cl)
    
    res.model <- model.Dev(
      training = train_data,
      method = methods,
      sig = pre_var
    )
    stopCluster(cl)
    
    
    ml.auc <- lapply(list_train_vali_Data, function(x) {
      res.tmp <- cal.model.auc(res.by.model.Dev = res.model, cohort.for.cal = x, sig = pre_var)
      return(res.tmp)
    })
    names(ml.auc) <- names(list_train_vali_Data)
    
    ml.roc <- lapply(list_train_vali_Data, function(x) {
      res.tmp <- cal.model.roc(res.by.model.Dev = res.model, cohort.for.cal = x, sig = pre_var)
      return(res.tmp)
    })
    names(ml.roc) <- names(list_train_vali_Data)
    
    res <- list()
    res[["model"]] <- res.model
    res[["auc"]] <- ml.auc
    res[["roc"]] <- ml.roc
    res[["sig.gene"]] <- pre_var
    
    return(res)
  } else {
    print("Please provide the correct parameters")
  }
}
