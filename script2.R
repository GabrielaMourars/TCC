#-----------------------------------------------------------------------------
# Classificação do bandgap de TMDs em semimetal (0) e semicondutor (1) 
#-----------------------------------------------------------------------------
pacotes <- c('corrplot', 'caret', 'Hmisc','pROC', "tidyverse", 'gbm',
             'randomForest', "kableExtra", 'rpart', "rpart.plot", 'class',
             'naivebayes', 'tidyr')

options(rgl.debug = TRUE)

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}
#-----------------------------------------------------------------------------

# função para o cálculo das métricas de performance
calculo_metricas <- function(predict, holdout){
  cm <- confusionMatrix(predict,holdout) %>% print()
  df <- data.frame(cm$byClass)
  kbl(df) %>% kable_styling(bootstrap_options = "striped", full_width = F,
                            position = "float_left")
  return(df)
}

load('df2.Rda')

# gráfico mostrando a quatidade de tmds semimetal (bandgap = 0) e semicondutor (bandgap != 0)
  ggplot(df2, aes(x=dir_gap2, fill = dir_gap2)) + 
  geom_bar(position = "dodge") +
  labs(x = 'Bandgap Direto', y = 'Contagem', fill = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA),
        legend.position= 'none')

# distribuição para os calcogênios
ggplot(df2, aes(x=dir_gap2,color=C_AN, y = phase)) + 
  geom_jitter() +
  labs(x = 'Bandgap Direto', y = 'Contagem', fill = 'Calcogênio') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA))

# separar as variáveis que vão ser usadas para a análise
df3 <- df2 %>% select(hform, mass, C_AW, TM_AW, TM_AD, TM_FH, TM_SH, 
                      TM_MP, TM_BP, TM_FIE, TM_AN, TM_GN, dir_gap2)
# salvar o data frame
save(df3,file="df3.Rda")

# carregar o data frame
load("df3.Rda")

# separar a amostra em 20% para teste e 80% para treino
set.seed(123)
sample_size = floor(0.8*nrow(df3))
picked = sample(seq_len(nrow(df3)),size = sample_size)

training = df3[picked,]
holdout = df3[-picked,]

# transformar os níveis de dir_gap2 em caracteres (é preciso para calcular 
# as probabilidades das classes)
# o train do caret não funciona se não for assim
levels(training$dir_gap2)=c("Yes","No")
levels(holdout$dir_gap2) = c("Yes","No")

# análise de classificação

# aplicação de modelos de machine learning para prever se o bandgap é 0 ou 1
#--------------------------------------------------------------------------
#                               arvore de decisão
#--------------------------------------------------------------------------
fitControl <- caret::trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10,
  ## Estimate class probabilities
  classProbs = TRUE,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = twoClassSummary)

# grade para encontrar os melhores parâmetros
treeGrid <- expand.grid(cp = c(0.01, 0.1, 0.3, 0.5, 0.7))

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
tree.training <- caret::train(dir_gap2 ~., data = training, 
                            method = "rpart", 
                            trControl = fitControl, 
                            #verbose = FALSE, 
                            tuneGrid = treeGrid,
                            ## Specify which metric to optimize
                            metric = "ROC")

# encontrar os melhores valores dos parâmetros
# a função train escolhe o melhor modelo com o maior valor de performance
# a função tolerance é usada para encontrar o modelo menos complexo
whichTwoPct <- tolerance(tree.training$results, metric = "ROC", 
                         tol = 2, maximize = TRUE)  

tree.training$results[whichTwoPct,1:6]

# cv
plot(tree.training, method = 'cv')

# preditores mais importantes
varImp(tree.training, scale = F) %>% plot(top = 5)

# cálculo da predição e métricas de performance
tree.pred = predict(tree.training, holdout, type = 'raw')

# comparação entre os dados e a predição
#result = data.frame(holdout$dir_gap2, tree.pred)

#tree.metricas <- calculo_metricas(tree.pred,holdout$dir_gap2) %>% 
#  dplyr::rename(tree = cm.byClass) %>% rownames_to_column('row_names')

# cruva roc
#tree.pred <- as.numeric(tree.pred)
#tree.roc <- plot(roc(holdout[,14], tree.pred), print.auc = TRUE,
#               levels = c('Yes','No'),
#               max.auc.polygon = TRUE,
#               main = 'Curva ROC - Decision Tree')
#tree.metricas <- tree.metricas %>%  add_row(row_names = 'ROC', tree = tree.roc$auc[1])

#--------------------------------------------------------------------------
#                               random forest
#--------------------------------------------------------------------------
# 10 folds repeat 10 times
fitControl <- caret::trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10,
  ## Estimate class probabilities
  classProbs = TRUE,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = twoClassSummary)

# grade para encontrar os melhores parâmetros
rfGrid <- expand.grid(mtry = c(1:15))

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
rf.training <- caret::train(dir_gap2 ~., data = training, 
                        method = "rf", 
                        trControl = fitControl, 
                        verbose = FALSE, 
                        tuneGrid = rfGrid,
                        ## Specify which metric to optimize
                        metric = "ROC")

# encontrar os melhores valores dos parâmetros
# a função train escolhe o melhor modelo com o maior valor de performance
# a função tolerance é usada para encontrar o modelo menos complexo
#whichTwoPct <- tolerance(rf.training$results, metric = "ROC", tol = 2, maximize = TRUE)  
#rf.training$results[whichTwoPct,1:6]$mtry

plot(rf.training, method = 'cv')

# preditores mais importantes
varImp(rf.training, scale = F) %>% plot(top = 5)

# cálculo da predição e métricas de performance
#rf.pred = predict(rf.training, holdout, type = 'raw')

# comparação entre os dados e a predição
#result = data.frame(holdout$dir_gap2, rf.pred)
#print(result)

#rf.metricas <- calculo_metricas(rf.pred,holdout$dir_gap2) %>% 
#  dplyr::rename(rf = cm.byClass) %>% rownames_to_column('row_names')

# cruva roc
#rf.pred <- as.numeric(rf.pred)
#rf.roc <- plot(roc(holdout[,19], rf.pred), print.auc = TRUE,
#               levels = c('Yes','No'),
#     max.auc.polygon = TRUE,
#     main = 'Curva ROC - Random Florest')
#rf.metricas <- rf.metricas %>%  add_row(row_names = 'ROC', rf = rf.roc$auc[1])
#--------------------------------------------------------------------------
#                         gradient boosted modeling
#--------------------------------------------------------------------------
#library('gbm') # gradient boosted modeling

fitControl <- caret::trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10,
  ## Estimate class probabilities
  classProbs = TRUE,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = twoClassSummary)

# grade para encontrar os melhores parâmetros
gbmGrid <- expand.grid(interaction.depth = c(1, 5, 9), 
                        n.trees = (1:30)*50, 
                        shrinkage = c(0.1, 0.3, 0.001),
                        n.minobsinnode = c(10,15,20))

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
boost.training <- train(dir_gap2 ~., data = training, 
                 method = "gbm", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 tuneGrid = gbmGrid,
                 ## Specify which metric to optimize
                 metric = "ROC")


# encontrar os melhores valores dos parâmetros
# a função train escolhe o melhor modelo com o maior valor de performance
# a função tolerance é usada para encontrar o modelo menos complexo
#whichTwoPct <- tolerance(boost.training$results, metric = "ROC", 
#                         tol = 2, maximize = TRUE)  

#boost.training$results[whichTwoPct,1:6]
#boost.training$bestTune

# plotar o resampling
plot(boost.training,method = "cv")

# preditores mais importantes
varImp(boost.training, scale = F) %>% plot(top = 5)

# predição
boost.pred <- predict(boost.training,holdout)
#print(boost.pred)

# comparação entre os dados e a predição
#result = data.frame(holdout$dir_gap2, boost.pred)
#print(result)

# cálculo das métricas
#boost.metricas <- calculo_metricas(boost.pred,holdout$dir_gap2) %>% 
#  dplyr::rename(GBM = cm.byClass) %>% rownames_to_column('row_names')

# curva roc
#boost.pred <- as.numeric(boost.pred)
#boost.roc <- plot(roc(holdout[,14], boost.pred), print.auc = TRUE,
#               levels = c('Yes','No'),
#               max.auc.polygon = TRUE,
#               main = 'Curva ROC - Gradient Boosted Modeling')
#boost.metricas <- boost.metricas %>% add_row(row_names = 'ROC', GBM = boost.roc$auc[1])
#--------------------------------------------------------------------------
#                K-Nearest Neighbor Neighbors Classifier
#--------------------------------------------------------------------------

fitControl <- caret::trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10,
  ## Estimate class probabilities
  classProbs = TRUE,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = twoClassSummary)

# grade para encontrar os melhores parâmetros
knnGrid <- expand.grid(k = c(1:50))

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
knn.training <- caret::train(dir_gap2 ~., data = training, 
                        method = "knn", 
                        trControl = fitControl,
                        tuneGrid = knnGrid,
                        ## Specify which metric to optimize
                        metric = "ROC")

# encontrar os melhores valores dos parâmetros
# a função train escolhe o melhor modelo com o maior valor de performance
# a função tolerance é usada para encontrar o modelo menos complexo
#whichTwoPct <- tolerance(knn.training$results, metric = "ROC", 
#                         tol = 2, maximize = TRUE)  

#cat("best model within 2 pct of best:\n")
#knn.training$results[whichTwoPct,1:6] # k = 2 
#knn.training$bestTune

# plotar o resampling
plot(knn.training,method = "cv")

# preditores mais importantes
varImp(knn.training, scale = F) %>% plot(top = 5)

# predição
knn.pred <- predict(knn.training,holdout)
#print(knn.pred)
roc(holdout,knn.training)

# comparação entre os dados e a predição
#result = data.frame(holdout$dir_gap2, knn.pred)
#print(result)

# cálculo das métricas
#knn.metricas <- calculo_metricas(knn.pred,holdout$dir_gap2) %>% 
#  dplyr::rename(knn = cm.byClass) %>% rownames_to_column('row_names')

# curva roc
#knn.pred <- as.numeric(boost.pred)
#knn.roc <- plot(roc(holdout[,14], knn.pred), print.auc = TRUE,
#                  levels = c('Yes','No'),
#                  max.auc.polygon = TRUE,
#                  main = 'K-Nearest Neighbor Neighbors Classifier')
#knn.metricas <- knn.metricas %>%  add_row(row_names = 'ROC', knn = knn.roc$auc[1])
#--------------------------------------------------------------------------
#                   Gaussian Naive Bayes Classifier
#--------------------------------------------------------------------------

fitControl <- caret::trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10,
  ## Estimate class probabilities
  classProbs = TRUE,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = twoClassSummary,
  search = "random")

# grade para encontrar os melhores parâmetros
nvGrid <- expand.grid(usekernel = c(TRUE, FALSE),
                      laplace = c(0, 0.5, 1), 
                      adjust = c(0.75, 1, 1.25, 1.5)
                      )

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
nv.training <- train(dir_gap2 ~., data = training, 
                        method = "naive_bayes", 
                        trControl = fitControl, 
                        verbose = FALSE, 
                        tuneGrid = nvGrid,
                        ## Specify which metric to optimize
                        metric = "ROC")

  summary(nv.training)

# encontrar os melhores valores dos parâmetros
# a função train escolhe o melhor modelo com o maior valor de performance
# a função tolerance é usada para encontrar o modelo menos complexo
#whichTwoPct <- tolerance(nv.training$results, metric = "ROC", 
#                         tol = 2, maximize = TRUE)  

#nv.training$results[whichTwoPct,1:6]

# plotar o resampling
plot(nv.training,method = "cv")

# preditores mais importantes
varImp(nv.training, scale = F) %>% plot(top = 5)

# predição
#nv.pred <- predict(nv.training,holdout)
#print(nv.pred)

# comparação entre os dados e a predição
#result = data.frame(holdout$dir_gap2, nv.pred)
#print(result)

# cálculo das métricas
#nv.metricas <- calculo_metricas(nv.pred,holdout$dir_gap2) %>% 
#  dplyr::rename(nv = cm.byClass) %>% rownames_to_column('row_names')

# curva roc
#nv.pred <- as.numeric(nv.pred)
#nv.roc <- plot(roc(holdout[,19], nv.pred), print.auc = TRUE,
#                  levels = c('Yes','No'),
#                  max.auc.polygon = TRUE,
#                  main = 'Curva ROC - Gaussian Naive Bayes Classifier')
#nv.metricas <- nv.metricas %>% add_row(row_names = 'ROC', nv = nv.roc$auc[1])


#--------------------------------------------------------------------------
#                Comparação entre os modelos
#--------------------------------------------------------------------------
# ROC_AUC is similar to Balanced Accuracy, but there are some key differences: 
# Balanced Accuracy is calculated on predicted classes, and ROC_AUC is 
# cal ~]culated on predicted scores for each data which can't be obtained on the 
# confusion matrix.

avalia <- function(modelo){
  #nome <- as.character(modelo)
  #p_treino <- predict(modelo, training, type='prob') # Probabilidade predita
  c_treino <- predict(modelo, training)              # Classificação
  
  #Base de teste
  #p_teste <- predict(modelo, holdout, type='prob')
  c_teste <- predict(modelo, holdout)
  
  c_treino <- as.numeric(c_treino)
  c_teste <- as.numeric(c_teste)
  
  roc_c_treino <- roc(training[,13], c_treino)
  roc_c_teste <- roc(holdout[,13], c_teste)
  
  dados <- data.frame(teste = roc_c_teste$auc,
    treino= roc_c_treino$auc) %>% t()  %>% as.data.frame() %>%
    rownames_to_column(var = 'rowname') %>% dplyr::rename(ROC = 2, 'Base' = 1)
  return(dados)
}

tree <- avalia(tree.training)
gbm <- avalia(boost.training)
rf <- avalia(rf.training)
knn <- avalia(knn.training)
nv <- avalia(nv.training)

# comparação entre base de tete e base de treino
# o gráfico contém somente os valores da curva roc para a base de teste
avaliacao <- purrr::reduce(list(tree,gbm,rf,knn,nv), 
                           dplyr::left_join, by = "Base")  %>% 
  dplyr::rename('Tree' = 2, 'GBM' = 3, 'RF' = 4, 'K-NN' = 5, 'NV' = 6) %>%
  t() %>% as.data.frame() %>% rownames_to_column(var = 'rowname') %>%
  'colnames<-'(.[1,]) %>% .[-1, ] %>% mutate_at(c('teste','treino'),as.numeric) %>%
  ggplot(aes(x = Base, y = teste, fill = Base))+
  geom_col()+
  labs(x = 'Modelo', y = 'ROC', title = 'Curva ROC na base de teste')+
  geom_text(aes(label = round(teste,3)), vjust = 1.5, colour = "white")+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA),
        legend.position="none")
avaliacao
# gmb é o modelo com o maior valor AUC 

#variaveis importantes dos modelos com maior roc
var_rf <- ggplot(varImp(rf.training, scale = F), top = 5) +
  labs(x = 'Variável', y = 'Importância', title ='Random Forest') + 
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA)
        )
var_boost <- ggplot(varImp(boost.training, scale = F), top = 5) +
  labs(x = 'Variável', y = 'Importância', title ='GBM') + 
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA)
        #plot.title = element_text(hjust = 0.5)
  )

  
# figura final
grid.arrange(arrangeGrob(avaliacao,left = textGrob("a)", x = unit(1, "npc"), y = unit(.975, "npc"))),
             arrangeGrob(var_rf, left = textGrob("b)", x = unit(1, "npc"), 
                          y = unit(.975, "npc"))),
             #arrangeGrob(var_knn,left = textGrob("c)", x = unit(1, "npc"), 
            #              y = unit(.95, "npc"))),
             #arrangeGrob(rf.training$bestTune %>%  
            #               tableGrob(rows = NULL, theme = ttheme_default(base_size = 10))
             #            ),
             #arrangeGrob(knn.training$bestTune %>%  
            #               tableGrob(rows = NULL, theme = ttheme_default(base_size = 10)),
             #            top = 'KNN'),
             ncol = 2, 
             layout_matrix = cbind(c(1), c(2,2)), heights = c(0.6,0.6)
)
