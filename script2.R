#-----------------------------------------------------------------------------
# Classificação do bandgap de TMDs em semimetal (0) e semicondutor (1) 
#-----------------------------------------------------------------------------
pacotes <- c('corrplot', 'caret', 'Hmisc','pROC', "tidyverse", 'gbm',
             'randomForest', "kableExtra", 'rpart', "rpart.plot", 'class',
             'naivebayes')

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
library(PerformanceAnalytics)
library(see)
library(ggraph)
library(psych)
library(reshape2)
library("plotly")
library("knitr")

# função para o cálculo das métricas de performance
calculo_metricas <- function(predict, holdout){
  cm <- confusionMatrix(predict,holdout) %>% print()
  df <- data.frame(cm$byClass)
  kbl(df) %>% kable_styling(bootstrap_options = "striped", full_width = F,
                            position = "float_left")
  return(df)
}


# estimativa da densidade de kernel
kernel <- ggplot(df2, aes(x = C_AN, y = dir_gap, fill = phase)) + 
  geom_violin() +
  labs(x = 'Número Atômico', y = 'Bandgap Direto', fill = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA))
kernel
# para determinar se há algum tipo de associação entre o bandgap e as outras 
# variáveis do banco de dados será realizado um teste de v de cramer para
# as variáveis categóricas e a associação por coeficiente de pearson para
# as variáveis numéricas

# variável preditora: numérica
# variávels explicativas: numérica (correlação de pearson) categórica ()

#variáveis categóricas

# gráfico de correlação usando todas as variáveis como numéricas

##df2 %>% select(-phase,-C_AN,-TM_AN,-C_PN,-TM_PN,-TM_GN) %>%  cor() %>%
##  corrplot(method = 'color',order = 'alphabet', tl.cex = 0.5)

# determinar se a correlação é estatisticamente significante 

# os p-valores usando o correlation não são os mesmos de outras funções
##rcor <- rcorr(as.matrix(df2 %>% 
##        select(-phase,-C_AN,-TM_AN,-C_PN,-TM_PN,-TM_GN)), type = "pearson")

# tabela com os valores das correlções
  
##flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
##  library(tidyr)
##  library(tibble)
##  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
##  cor_r <- gather(cor_r, column, cor, -1)
##  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
##  cor_p <- gather(cor_p, column, p, -1)
##  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
##  cor_p_matrix
##}

##flat_cor_mat(rcor$r, rcor$P) %>% 
##  filter(row == 'dir_gap' & p < 0.05) %>% kable() %>%
##  kable_styling(bootstrap_options = "striped", 
##                full_width = F, 
##                font_size = 12)

##flat_cor_mat(rcor$r, rcor$P) %>% 
##  filter((row == 'dir_gap'| column == 'dir_gap') & p < 0.05) 

# variáveis categóricas
# o dir_gap deve ser normal 
##library(hrbrthemes) #temas para o gráfico
##library(nortest) #teste de normalidade

##sf.test(df2$dir_gap)

##ggplot(df2, aes(x = dir_gap)) +
##  geom_histogram()+
##  stat_function(fun = dnorm, 
##                args = list(mean = mean(df2$dir_gap),
##                            sd = sd(df2$dir_gap)),
##                aes(color = "Curva Normal Teórica"),
##                size = 2) +
##  scale_color_manual("Legenda:",
##                     values = "#FDE725FF") +
##  theme(legend.position = "bottom")

##library('stats') # teste anova

##AOV <- aov(dir_gap ~ phase+ TM_AN + C_AN + C_PN + TM_PN + TM_GN, data = df2) 
##summary(AOV)

# classificar como zero quando o bandgap é 0 e como 1 quando o bandgap é 
# diferente de zero

df2 <- df2 %>% mutate(dir_gap2 = ifelse(dir_gap != 0, 1,0))
df2$dir_gap2 <- as.factor(df2$dir_gap2)

# gráfico mostrando a quatidade de tmds semimetal (bandgap = 0) e semicondutor (bandgap != 0)
ggplot(df2, aes(x=dir_gap2,fill=phase)) + 
  geom_bar(position="stack") +
  labs(x = 'Bandgap Direto', y = 'Contagem', fill = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA))

# distribuição para os calcogênios
ggplot(df2, aes(x=dir_gap2,fill=C_AN)) + 
  geom_bar(position="stack") +
  labs(x = 'Bandgap Direto', y = 'Contagem', fill = 'Calcogênio') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA))

ggplot(df2, aes(x=dir_gap, y = hform, color=phase)) + 
  geom_point() +
  labs(x = 'Bandgap Direto', y = 'Energia de Formação', fill = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA))

ggplot(df2, aes(x=dir_gap, y = a, color=phase)) + 
  geom_point() +
  labs(x = 'Bandgap Direto', y = 'Energia de Formação', fill = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA))

# separar as variáveis que vão ser usadas para a análise
df3 <- df2 %>% select(magmom,hform, phase, C_AW,C_AD,TM_AR,C_AR,C_VDW,TM_AV,C_AV,C_COV,
                      C_MP,C_BP,TM_FIE,C_FIE,C_EN,TM_AN,C_AN,dir_gap2)
summary(df2$magmom)
# salvar o data frame
save(df3,file="df3.Rda")

# carregar o data frame

load("df3.Rda")
# separar a amostra em 20% para teste e 80% para treino
set.seed(123)
sample_size = floor(0.7*nrow(df3))
picked = sample(seq_len(nrow(df3)),size = sample_size)

training = df3[picked,]
holdout = df3[-picked,]

# transformar os níveis de dir_gap2 em caracteres (é preciso para calcular 
# as probabilidades das classes)
# o train do caret não funciona se não for assim
levels(training$dir_gap2)=c("Yes","No")
levels(holdout$dir_gap2) = c("Yes", "No")

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
treeGrid <- expand.grid(cp = c(0.01, 0.1, 0.3, 0.5))

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

cat("best model within 2 pct of best:\n")
tree.training$results[whichTwoPct,1:6]
# cp       ROC      Sens      Spec     ROCSD    SensSD
# 1 0.01 0.8200116 0.8794444 0.4916667 0.1355577 0.1218966

# # update do modelo com os melhores parametros mtry = 6
tree.training.update <- update(tree.training, param = list(cp = 0.01))
summary(tree.training.update)
cv.tree <- plot(tree.training.update, method = 'cv')

# preditores mais importantes
varImp(tree.training.update, scale = F) %>% plot()

# cálculo da predição e métricas de performance
tree.pred = predict(tree.training, holdout, type = 'raw')

# comparação entre os dados e a predição
result = data.frame(holdout$dir_gap2, tree.pred)
print(result)

# comparação entre os dados e a predição
result = data.frame(holdout$dir_gap2, rf.pred)
print(result)

tree.metricas <- calculo_metricas(tree.pred,holdout$dir_gap2) %>% 
  rename(tree = cm.byClass) %>% rownames_to_column('row_names')

# cruva roc
tree.pred <- as.numeric(tree.pred)
tree.roc <- plot(roc(holdout[,19], tree.pred), print.auc = TRUE,
               levels = c('Yes','No'),
               max.auc.polygon = TRUE,
               main = 'Curva ROC - Decision Tree')
tree.metricas <- tree.metricas %>%  add_row(row_names = 'ROC', tree = tree.roc$auc[1])

fitControl <- caret::trainControl(## 10-fold CV
  method = "cv",
  number = 10,
  ## repeated ten times
  #repeats = 10
  ## Estimate class probabilities
  classProbs = TRUE, # To calculate ROC
  ## Evaluate performance using 
  ## the following function
  summaryFunction = twoClassSummary
  )

# grade para encontrar os melhores parâmetros
treeGrid <- expand.grid(cp = c(0,0.01, 0.1, 0.3, 0.5))

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
tree.training <- caret::train(dir_gap2 ~., data = training, 
                            method = "rpart", 
                            trControl = fitControl, 
                            #verbose = FALSE, 
                            tuneGrid = treeGrid,
                            ## Specify which metric to optimize
                            metric = "ROC"
                            )
plot(tree.training)
# encontrar os melhores valores dos parâmetros
# a função train escolhe o melhor modelo com o maior valor de performance
# a função tolerance é usada para encontrar o modelo menos complexo
whichTwoPct <- tolerance(tree.training$results, metric = "ROC", 
                         tol = 2, maximize = TRUE)  

#cat("best model within 2 pct of best:\n")
tree.training$results[whichTwoPct,1:6]
#    cp       ROC      Sens      Spec     ROCSD     SensSD
#     0 0.8574074 0.9319444 0.4833333 0.1092738 0.07869009

# # update do modelo com os melhores parametros mtry = 6
tree.training.update <- update(tree.training, param = list(cp = 0))
summary(tree.training.update)
plot(tree.training.update, method = 'cv')

# preditores mais importantes
varImp(tree.training.update, scale = F) %>% plot()

# cálculo da predição e métricas de performance
tree.pred = predict(tree.training, holdout, type = 'raw')

# comparação entre os dados e a predição
result = data.frame(holdout$dir_gap2, tree.pred)
print(result)

tree.metricas <- calculo_metricas(tree.pred,holdout$dir_gap2) %>% 
  rename(tree = cm.byClass) %>% rownames_to_column('row_names')

# cruva roc
tree.pred <- as.numeric(tree.pred)
tree.roc <- plot(roc(holdout[,19], tree.pred), print.auc = TRUE,
               levels = c('Yes','No'),
               max.auc.polygon = TRUE,
               main = 'Curva ROC - Decision Tree')
tree.metricas <- tree.metricas %>%  add_row(row_names = 'ROC', tree = tree.roc$auc[1])
                                        
#--------------------------------------------------------------------------
#                               random forest
#--------------------------------------------------------------------------
# 10 folds repeat 3 times
control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3)

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
rfGrid <- expand.grid(mtry = c(1:10))

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
whichTwoPct <- tolerance(rf.training$results, metric = "ROC", 
                         tol = 2, maximize = TRUE)  

cat("best model within 2 pct of best:\n")
rf.training$results[whichTwoPct,1:6]$mtry

# # update do modelo com os melhores parametros mtry = 6
rf.training.update <- update(rf.training, param = list(mtry=rf.training$results[whichTwoPct,1:6]$mtry))
summary(rf.training.update)
plot(rf.training.update, method = 'cv')

# preditores mais importantes
varImp(rf.training.update, scale = F) %>% plot()

# cálculo da predição e métricas de performance
rf.pred = predict(rf.training, holdout, type = 'raw')

# comparação entre os dados e a predição
result = data.frame(holdout$dir_gap2, rf.pred)
print(result)

rf.metricas <- calculo_metricas(rf.pred,holdout$dir_gap2) %>% 
  rename(rf = cm.byClass) %>% rownames_to_column('row_names')

# cruva roc
rf.pred <- as.numeric(rf.pred)
rf.roc <- plot(roc(holdout[,19], rf.pred), print.auc = TRUE,
               levels = c('Yes','No'),
     max.auc.polygon = TRUE,
     main = 'Curva ROC - Random Florest')
rf.metricas <- rf.metricas %>%  add_row(row_names = 'ROC', rf = rf.roc$auc[1])
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
whichTwoPct <- tolerance(boost.training$results, metric = "ROC", 
                         tol = 2, maximize = TRUE)  

cat("best model within 2 pct of best:\n")
boost.training$results[whichTwoPct,1:6]
boost.training$bestTune

# plotar o resampling
trellis.par.set(caretTheme())
plot(boost.training,method = "cv")
summary(boost.training$finalModel)

# update do modelo com os melhores parametros
boost.training.update <- update(boost.training, 
                                param = list(n.trees = boost.training$results[whichTwoPct,1:6]$n.trees,
                                                             interaction.depth = boost.training$results[whichTwoPct,1:6]$interaction.depth,
                                                             shrinkage = boost.training$results[whichTwoPct,1:6]$shrinkage,
                                                             n.minobsinnode = boost.training$results[whichTwoPct,1:6]$n.minobsinnode))
summary(boost.training.update)

plot(boost.training.update, method = "cv")

# preditores mais importantes
varImp(boost.training.update, scale = F) %>% plot()

# predição
boost.pred <- predict(boost.training,holdout)
print(boost.pred)

# comparação entre os dados e a predição
result = data.frame(holdout$dir_gap2, boost.pred)
print(result)

# cálculo das métricas
boost.metricas <- calculo_metricas(boost.pred,holdout$dir_gap2) %>% 
  rename(GBM = cm.byClass) %>% rownames_to_column('row_names')

# curva roc
boost.pred <- as.numeric(boost.pred)
boost.roc <- plot(roc(holdout[,19], boost.pred), print.auc = TRUE,
               levels = c('Yes','No'),
               max.auc.polygon = TRUE,
               main = 'Curva ROC - Gradient Boosted Modeling')
boost.metricas <- boost.metricas %>% add_row(row_names = 'ROC', GBM = boost.roc$auc[1])
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
gbmGrid <- expand.grid(k = c(1:50))

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
knn.training <- caret::train(dir_gap2 ~., data = training, 
                        method = "knn", 
                        trControl = fitControl,
                        tuneGrid = gbmGrid,
                        ## Specify which metric to optimize
                        metric = "ROC")

# encontrar os melhores valores dos parâmetros
# a função train escolhe o melhor modelo com o maior valor de performance
# a função tolerance é usada para encontrar o modelo menos complexo
whichTwoPct <- tolerance(knn.training$results, metric = "ROC", 
                         tol = 2, maximize = TRUE)  

#cat("best model within 2 pct of best:\n")
knn.training$results[whichTwoPct,1:6] # k = 2 
knn.training$bestTune

# plotar o resampling
trellis.par.set(caretTheme())
plot(knn.training,method = "cv")
summary(knn.training$finalModel)

# update do modelo com os melhores parametros
knn.training.update <- update(knn.training, param = list(k = 2))
summary(knn.training.update)
plot(knn.training.update, method = "cv")

# preditores mais importantes
varImp(knn.training.update, scale = F) %>% plot()

# predição
knn.pred <- predict(knn.training,holdout)
print(knn.pred)
#roc(holdout[,-19],knn.training)

# comparação entre os dados e a predição
result = data.frame(holdout$dir_gap2, knn.pred)
print(result)

# cálculo das métricas
knn.metricas <- calculo_metricas(knn.pred,holdout$dir_gap2) %>% 
rename(knn = cm.byClass) %>% rownames_to_column('row_names')

# curva roc
knn.pred <- as.numeric(boost.pred)
knn.roc <- plot(roc(holdout[,19], knn.pred), print.auc = TRUE,
                  levels = c('Yes','No'),
                  max.auc.polygon = TRUE,
                  main = 'K-Nearest Neighbor Neighbors Classifier')
knn.metricas <- knn.metricas %>%  add_row(row_names = 'ROC', knn = knn.roc$auc[1])
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
whichTwoPct <- tolerance(nv.training$results, metric = "ROC", 
                         tol = 2, maximize = TRUE)  

cat("best model within 2 pct of best:\n")
nv.training$results[whichTwoPct,1:6]
nv.training$bestTune

# plotar o resampling
trellis.par.set(caretTheme())
plot(nv.training,method = "cv")
summary(nv.training$finalModel)

# update do modelo com os melhores parametros
nv.training.update <- update(nv.training, param = list(laplace = 1,
                                                       usekernel = TRUE,
                                                       adjust = 0.75))
summary(nv.training.update)

plot(nv.training.update, method = "cv")

# preditores mais importantes
varImp(nv.training.update, scale = F) %>% plot()

# predição
nv.pred <- predict(nv.training,holdout)
print(nv.pred)

# comparação entre os dados e a predição
result = data.frame(holdout$dir_gap2, nv.pred)
print(result)

# cálculo das métricas
nv.metricas <- calculo_metricas(nv.pred,holdout$dir_gap2) %>% 
  rename(nv = cm.byClass) %>% rownames_to_column('row_names')

# curva roc
nv.pred <- as.numeric(nv.pred)
nv.roc <- plot(roc(holdout[,19], nv.pred), print.auc = TRUE,
                  levels = c('Yes','No'),
                  max.auc.polygon = TRUE,
                  main = 'Curva ROC - Gaussian Naive Bayes Classifier')
nv.metricas <- nv.metricas %>% add_row(row_names = 'ROC', nv = nv.roc$auc[1])


#--------------------------------------------------------------------------
#                Comparação entre os modelos
#--------------------------------------------------------------------------
# ROC_AUC is similar to Balanced Accuracy, but there are some key differences: 
# Balanced Accuracy is calculated on predicted classes, and ROC_AUC is 
# calculated on predicted scores for each data which can't be obtained on the 
# confusion matrix.
tabela <- purrr::reduce(list(rf.metricas,knn.metricas, boost.metricas,
                             tree.metricas, nv.metricas), 
                        dplyr::left_join) %>% kable(caption = 'Comparação das predições', digits = 3) %>% 
                        kable_classic_2(full_width = F)
tabela

tabela_2 <- data.frame(training = c(knn.training$results[whichTwoPct,1:6][1,'ROC'],
                                    rf.training$results[whichTwoPct,1:6][1,'ROC'],
                                    boost.training$results[whichTwoPct,1:6][1,'ROC'],
                                    tree.training$results[whichTwoPct,1:6][1,'ROC'],
                                    nv.training$results[whichTwoPct,1:6][1,'ROC']),
                       holdout = c(knn.roc$auc, rf.roc$auc, boost.roc$auc, 
                                   tree.roc$auc, nv.roc$auc), 
                       row.names = c('knn', 'rf', 'boost', 'tree', 'nv')) %>% 
  kable(caption = 'Curva ROC', digits = 3) %>% 
  kable_classic_2(full_width = F)
tabela_2


  



