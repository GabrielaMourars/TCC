#-----------------------------------------------------------------------------
# Classificação do bandgap de heteroestruturasd de TMDs 
#-----------------------------------------------------------------------------
pacotes <- c('corrplot', 'caret', 'Hmisc','pROC', "tidyverse", 'gbm',
             'randomForest', "kableExtra", 'rpart', "rpart.plot", 'class',
             'naivebayes', 'gtools', 'corrgram','diffdf', 'MLmetrics', 
             'xgboost', 'plyr','scales', 'ggplot2', "patchwork")

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
#-------------------------------------------------------------------------------
# Construir o banco de dados com permutações de tmds
#-------------------------------------------------------------------------------
# novo banco de dados formado somente com TMDs não metálicos (bandgap =/ 0)
load('df1.Rda')

df4 <- df1 %>% drop_na() %>% filter(dir_gap != 0) %>%
  left_join(number_key_values %>% filter(key == "cbm", keep = TRUE)%>% 
              select(value,id)) %>% rename('cbm'=value)  
df4 <- df4 %>% left_join(number_key_values %>% filter(key == "vbm", keep = TRUE) %>% 
                           select(value,id)) %>% rename('vbm'=value)
df4 <- df4 %>% left_join(text_key_values %>% filter(key == 'name', keep = TRUE) %>% 
                           select(value,id)) %>% rename('name'=value)
df4 <- df4 %>% filter(C_AN != 8) # banco de dados somente com TMDs
df4$C_AN <- as.factor(df4$C_AN)

# transformar os parametros para o tipo correto

df4 <- df4 %>% mutate_at(c('TM_FH', 'TM_AR', 'C_AR'), as.numeric)

df4 <- df4 %>% mutate_at(c('phase', 'C_AN', 'TM_AN', 'C_PN', 'TM_PN', 'TM_GN'),
                         as.factor)

save(df4, file = 'df4.Rda')

load('df4.Rda')
            
# análise descritiva do banco de dados
summary(df4)

# fazer um box plot
  ggplot(df4, aes(x = C_AN, y = dir_gap, fill = phase)) +
  geom_boxplot() +
  labs(x = 'Número Atômico', y = 'Bandgap Direto', fill = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA))

  # gráfico kernel
ggplot(df4, aes(x = C_AN, y = dir_gap, fill = phase)) + geom_violin() +
    labs(x = 'Número Atômico', y = 'Bandgap Direto', fill = 'Fase') +
    theme(panel.background = element_rect("white"),
          panel.grid = element_line("grey95"),
          panel.border = element_rect(NA))

# relação entre a energia de formação e o bandgap para cada fase do tmd
ggplot(df4, aes(x=dir_gap, y = hform, color=phase)) + 
  geom_point() +
  labs(x = 'Bandgap Direto', y = 'Energia de Formação', fill = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA))

# relação entre a numero atomico mt e o bandgap para cada fase do tmd
ggplot(df4, aes(x=dir_gap, y = TM_AN, color=phase)) + 
  geom_point() +
  labs(x = 'Bandgap Direto', y = 'Número Atômico do Metal de Transição', fill = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA))

# relação entre o momento magético e o bandgap para cada fase do tmd
ggplot(df4, aes(x=dir_gap, y = magmom, color=phase)) + 
  geom_point() +
  labs(x = 'Bandgap Direto', y = 'Momento Magnético', fill = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA))

# relação entre a constante de rede e o bandgap para cada fase do tmd
ggplot(df4, aes(x=dir_gap, y = a, color=phase)) + 
  geom_point() +
  labs(x = 'Bandgap Direto', y = 'Constante de Rede', fill = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA))

# relação entre o numero atomico C e o bandgap para cada fase do tmd
ggplot(df4, aes(x=dir_gap, y = C_AN, color=phase)) + 
  geom_point() +
  labs(x = 'Bandgap Direto', y = 'Número Atômico do Calcogênio', fill = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA))

## Regra de Anderson
# tipo I se VBMa < VBMb < CBMb < CBMa
# tipo II se VBMa < VMBb < CMBa < CMBb
# tipo III se VBMa < CMBa < VBMb < CBMb

# montar uma matrix com os tipos de heterojunções de acordo com a regra de anderson

# todas as combinações possíveis com os tmds
options(expressions = 1e5)
data <- permutations(38,2,df4[,1])
data <- as.data.frame(data)

# adicionar vbm e cbm para cada e tmd

data <- data %>% left_join(df4 %>% select(vbm,id), by = join_by(V1 == id)) %>% dplyr::rename(vbma = vbm)
data <- data %>% left_join(df4 %>% select(vbm,id), by = join_by(V2 == id)) %>% dplyr::rename(vbmb = vbm)
data <- data %>% left_join(df4 %>% select(cbm,id), by = join_by(V1 == id)) %>% dplyr::rename(cbma = cbm)
data <- data %>% left_join(df4 %>% select(cbm,id), by = join_by(V2 == id)) %>% dplyr::rename(cbmb = cbm)
data <- data %>% left_join(df4 %>% select(id,name), by = join_by(V1 == id)) %>% dplyr::rename(V1_name = 'name')
data <- data %>% left_join(df4 %>% select(id,name), by = join_by(V2 == id)) %>% dplyr::rename(V2_name = 'name')
data <- data %>% left_join(df4 %>% select(id,TM_AN), by = join_by(V1 == id)) %>% dplyr::rename(TM_AN_V1 = 'TM_AN')
data <- data %>% left_join(df4 %>% select(id,TM_AN), by = join_by(V2 == id)) %>% dplyr::rename(TM_AN_V2 = 'TM_AN')
data <- data %>% left_join(df4 %>% select(id,C_AN), by = join_by(V1 == id)) %>% dplyr::rename(C_AN_V1 = 'C_AN')
data <- data %>% left_join(df4 %>% select(id,C_AN), by = join_by(V2 == id)) %>% dplyr::rename(C_AN_V2 = 'C_AN')

tipo <- data.frame()
for (i in 1:nrow(data)) {
  
  if((data[i,3] < data[i,4]) & (data[i,4] < data[i,6]) & (data[i,6]< data[i,5])){
    tipo <- rbind(tipo,'I')
  } else if((data[i,3] < data[i,4]) & (data[i,4] < data[i,5]) & (data[i,5]< data[i,6])) {
    tipo <- rbind(tipo,'II')
  }else if((data[i,3] < data[i,5]) & (data[i,5] < data[i,4]) & (data[i,4]< data[i,6])){
    tipo <- rbind(tipo,'III')
  } else if ((data[i,4] < data[i,3]) & (data[i,3] < data[i,5]) & (data[i,5]< data[i,6])) {
    tipo <- rbind(tipo,'I')
  } else if ((data[i,4] < data[i,3]) & (data[i,3] < data[i,6]) & (data[i,6]< data[i,5])){
    tipo <- rbind(tipo,'II')
  } else if ((data[i,4] < data[i,6]) & (data[i,6] < data[i,3]) & (data[i,3]< data[i,5])) {
    tipo <- rbind(tipo,'III')
  } else {
    tipo <- rbind(tipo,'x')
  }
}

data <- cbind(data,tipo) %>% dplyr::rename('tipo' = X.III.)
data$tipo <- as.factor(data$tipo)  

data$V1 <- as.character(data$V1)
data$V2 <- as.character(data$V2)

# análise descritiva do banco de dados
ggplot(data, aes(V1_name, V2_name, fill= tipo)) + 
  geom_tile() +
  scale_x_discrete(guide = guide_axis(angle = 90))+
  labs(x = 'TMD', y = 'TMD', fill = 'Tipo') 

data %>% count(C_AN_V1,tipo)
data %>% select(tipo, TM_AN_V1) %>% count(TM_AN_V1,tipo)
count(data$tipo)%>% kable() %>%
  kable_styling(bootstrap_options = "striped",
                full_width = F,
                font_size = 28)
count(data$C_AN_V1)
count(data$TM_AN_V1)

data %>% select(tipo, TM_AN_V1) %>% group_by(tipo) %>% count() %>% kable() %>%
  kable_styling(bootstrap_options = "striped",
                full_width = F,
                font_size = 28)

# gráfico com a distribuição dos tipos de junção para cada metal de transição
TM_juncao <- data %>% select(tipo, TM_AN_V1) %>%
  ggplot(aes(x = TM_AN_V1, group = tipo, fill = tipo)) + 
  geom_bar()+
  xlab('Metal de Transição') +
  ylab('Frequência')

# gráfico com a distribuição dos tipos de junção para cada calcogênio
C_juncao <- data %>% select(tipo, C_AN_V1) %>%
  ggplot(aes(x = C_AN_V1,group = tipo, fill = tipo)) + 
  geom_bar()+
  xlab('Calcogênio') +
  ylab('Frequência')

TM_juncao + C_juncao + plot_layout(guides = "collect")

# gráfico com a distribuição dos tipos de junção para cada calcogênio
#ggplot(data, aes(x=C_AN_V1,fill=tipo)) + 
#geom_bar(position="stack") +
#labs(x = 'Bandgap Direto', y = 'Contagem', fill = 'Tipo') +
#theme(panel.background = element_rect("white"),
#      panel.grid = element_line("grey95"),
#      panel.border = element_rect(NA))

#-------------------------------------------------------------------------------
# Contruir um novo banco de dados com metade das observações do banco anterior
#-------------------------------------------------------------------------------

df <- combinations(38,2,df4[,1])
df <- as.data.frame(df)

df$V1 <- as.character(df$V1)
df$V2 <- as.character(df$V2)

df<- df %>% left_join(data) %>% select(-c(V1,V2, V1_name, V2_name))

# salvar o banco de dados
save(df, file = 'df.Rda')

# carregar o banco de dados
load('df.Rda')

#-------------------------------------------------------------------------------
# Análise de classificação
# aplicação de métodos de machine learning para classificar o tipo de junção
#-------------------------------------------------------------------------------

# função para o cálculo de métricas importantes
calculo_metricas <- function(predict, holdout){
  cm <- confusionMatrix(predict,holdout) %>% print()
  df <- data.frame(cm$byClass)
  kbl(df) %>% kable_styling(bootstrap_options = "striped", full_width = F,
                            position = "float_left")
  return(df)
}
# a métrica escolhida para fazer a análise entre os diferentes modelos é a acurácia
# neste momento tem 3 classes para classificação, não dá pra usar ROC

# separar a amostra em 20% para teste e 80% para treino
set.seed(123)
sample_size = floor(0.7*nrow(df))
picked = sample(seq_len(nrow(df)),size = sample_size)

training = df[picked,]
holdout = df[-picked,]

# transformar os níveis de dir_gap2 em caracteres (é preciso para calcular 
# as probabilidades das classes)
# o train do caret não funciona se não for assim
levels(training$tipo)=c("I","II", "III")
levels(holdout$tipo) = c("I", "II", "III")

# novo banco de dados com os valores preditos pelos modelos
yhat <- holdout %>% select(tipo)
#-----------------------------------------------------------------------------
#                              Árvore de Decisão 
#-----------------------------------------------------------------------------
fitControl <- caret::trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10,
  ## Estimate class probabilities
  classProbs = TRUE,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = multiClassSummary
  )

# grade para encontrar os melhores parâmetros
treeGrid <- expand.grid(cp = c(0.01, 0.1, 0.3, 0.5))

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
tree.training <- caret::train(tipo ~., data = training, 
                              method = "rpart", 
                              trControl = fitControl, 
                              #verbose = FALSE, 
                              tuneGrid = treeGrid,
                              ## Specify which metric to optimize
                              metric = 'Accuracy'
                              )

# encontrar os melhores valores dos parâmetros
# a função train escolhe o melhor modelo com o maior valor de performance
# a função tolerance é usada para encontrar o modelo menos complexo
whichTwoPct <- tolerance(tree.training$results, metric = "Accuracy", 
                         tol = 2, maximize = TRUE)  

cat("best model within 2 pct of best:\n")
cp <- tree.training$results[whichTwoPct,1:6]

# # update do modelo com os melhores parametros mtry = 6
tree.training.update <- update(tree.training, param = list(cp = cp))
summary(tree.training.update)
cv.tree <- plot(tree.training.update, method = 'cv')

# preditores mais importantes
varImp(tree.training.update, scale = F) %>% plot()

# cálculo da predição e métricas de performance
tree.pred = predict(tree.training, holdout, type = 'raw')

# comparação entre os dados e a predição
result = data.frame(holdout$tipo, tree.pred)
print(result)
yhat <- yhat %>% cbind(tree.pred)

tree.metricas <- calculo_metricas(tree.pred,holdout$tipo) #%>% 
  #rename(tree = cm.byClass) %>% rownames_to_column('row_names')

# cruva roc
tree.pred <- as.numeric(tree.pred)
tree.roc <- plot(roc(holdout[,9], tree.pred), print.auc = TRUE,
                 levels = c('Yes','No'),
                 max.auc.polygon = TRUE,
                 main = 'Curva ROC - Decision Tree')
tree.metricas <- tree.metricas %>%  add_row(row_names = 'ROC', tree = tree.roc$auc[1])

#--------------------------------------------------------------------------
#                               random forest
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
  summaryFunction = multiClassSummary)

# grade para encontrar os melhores parâmetros
rfGrid <- expand.grid(mtry = c(1:10))

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
rf.training <- caret::train(tipo ~., data = training, 
                            method = "rf", 
                            trControl = fitControl, 
                            verbose = FALSE, 
                            tuneGrid = rfGrid,
                            ## Specify which metric to optimize
                            metric = "Accuracy")

# encontrar os melhores valores dos parâmetros
# a função train escolhe o melhor modelo com o maior valor de performance
# a função tolerance é usada para encontrar o modelo menos complexo
whichTwoPct <- tolerance(rf.training$results, metric = "Accuracy", 
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
result = data.frame(holdout$tipo, rf.pred)
print(result)
yhat <- yhat %>% cbind(rf.pred)

rf.metricas <- calculo_metricas(rf.pred,holdout$tipo)

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
  summaryFunction = multiClassSummary)

# grade para encontrar os melhores parâmetros
# method ='gbm'
gbmGrid <- expand.grid(interaction.depth = c(1, 5, 9), 
                       n.trees = (1:30)*50, 
                       shrinkage = c(0.1, 0.3, 0.001),
                       n.minobsinnode = c(10,15,20))
# method = "xgbTree"
gbmGrid <- expand.grid(
                       nrounds = c(50, 100, 150),
                       max_depth = c(2, 3),
                       gamma = c(0),
                       eta = c(0.1, 0.4),
                       colsample_bytree = c(0.6, 0.8),
                       min_child_weight = c(1),
                       subsample = c(0.75, 1)
)

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
boost.training <- train(tipo ~., data = training, 
                        method = "xgbTree", 
                        trControl = fitControl, 
                        verbose = FALSE, 
                        tuneGrid = gbmGrid,
                        ## Specify which metric to optimize
                        metric = "Accuracy")


# encontrar os melhores valores dos parâmetros
# a função train escolhe o melhor modelo com o maior valor de performance
# a função tolerance é usada para encontrar o modelo menos complexo
whichTwoPct <- tolerance(boost.training$results, metric = "Accuracy", 
                         tol = 2, maximize = TRUE)  

cat("best model within 2 pct of best:\n")
boost.training$results[whichTwoPct,1:8]
boost.training$bestTune

# plotar o resampling
trellis.par.set(caretTheme())
plot(boost.training,method = "cv")
summary(boost.training$finalModel)

# update do modelo com os melhores parametros
boost.training.update <- update(boost.training, 
                                param = list(eta = boost.training$results[whichTwoPct,1:7]$eta,
                                             max_depth = boost.training$results[whichTwoPct,1:7]$max_depth,
                                             gamma = boost.training$results[whichTwoPct,1:7]$gamma,
                                             colsample_bytree = boost.training$results[whichTwoPct,1:7]$colsample_bytree,
                                             min_child_weight = boost.training$results[whichTwoPct,1:7]$min_child_weight,
                                             subsample = boost.training$results[whichTwoPct,1:7]$subsample,
                                             nrounds = boost.training$results[whichTwoPct,1:7]$nrounds))
summary(boost.training.update)

plot(boost.training, method = "cv")

# preditores mais importantes
varImp(boost.training, scale = F) %>% plot()

# predição
boost.pred <- predict(boost.training,holdout)
print(boost.pred)

# comparação entre os dados e a predição
result = data.frame(holdout$dir_gap2, boost.pred)
print(result)
yhat <- yhat %>% cbind(boost.pred)

# cálculo das métricas
boost.metricas <- calculo_metricas(boost.pred,holdout$tipo)

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
  summaryFunction = multiClassSummary)

# grade para encontrar os melhores parâmetros
knnGrid <- expand.grid(k = c(1:50))

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
knn.training <- caret::train(tipo ~., data = training, 
                             method = "knn", 
                             trControl = fitControl,
                             tuneGrid = knnGrid,
                             ## Specify which metric to optimize
                             metric = "Accuracy")

# encontrar os melhores valores dos parâmetros
# a função train escolhe o melhor modelo com o maior valor de performance
# a função tolerance é usada para encontrar o modelo menos complexo
whichTwoPct <- tolerance(knn.training$results, metric = "Accuracy", 
                         tol = 2, maximize = TRUE)  

#cat("best model within 2 pct of best:\n")
knn.training$results[whichTwoPct,1:6] # k = 2 
knn.training$bestTune

# plotar o resampling
trellis.par.set(caretTheme())
plot(knn.training,method = "cv")
summary(knn.training$finalModel)

# update do modelo com os melhores parametros
knn.training.update <- update(knn.training, param = list(k = knn.training$results[whichTwoPct,1:6]$k))
summary(knn.training.update)
plot(knn.training.update, method = "cv")

# preditores mais importantes
varImp(knn.training.update, scale = F) %>% plot()

# predição
knn.pred <- predict(knn.training,holdout)
print(knn.pred)
yhat <- yhat %>% cbind(knn.pred)

# comparação entre os dados e a predição
result = data.frame(holdout$tipo, knn.pred)
print(result)

# cálculo das métricas
knn.metricas <- calculo_metricas(knn.pred,holdout$tipo) 

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
  summaryFunction = multiClassSummary,
  search = "random")


# grade para encontrar os melhores parâmetros
nvGrid <- expand.grid(usekernel = c(TRUE, FALSE),
                      laplace = c(0, 0.5, 1), 
                      adjust = c(0.75, 1, 1.25, 1.5)
)

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
nv.training <- train(tipo ~., data = training, 
                     method = "naive_bayes", 
                     trControl = fitControl, 
                     verbose = FALSE, 
                     tuneGrid = nvGrid,
                     ## Specify which metric to optimize
                     metric = "Accuracy")

summary(nv.training)

# encontrar os melhores valores dos parâmetros
# a função train escolhe o melhor modelo com o maior valor de performance
# a função tolerance é usada para encontrar o modelo menos complexo
whichTwoPct <- tolerance(nv.training$results, metric = "Accuracy", 
                         tol = 2, maximize = TRUE)  

cat("best model within 2 pct of best:\n")
nv.training$results[whichTwoPct,1:6]
nv.training$bestTune

# plotar o resampling
trellis.par.set(caretTheme())
plot(nv.training,method = "cv")
summary(nv.training$finalModel)

# update do modelo com os melhores parametros
nv.training.update <- update(nv.training, param = list(laplace = nv.training$results[whichTwoPct,1:6]$laplace,
                                                       usekernel = nv.training$results[whichTwoPct,1:6]$usekernel,
                                                       adjust = nv.training$results[whichTwoPct,1:6]$adjust))
summary(nv.training.update)

plot(nv.training.update, method = "cv")

# preditores mais importantes
varImp(nv.training.update, scale = F) %>% plot()

# predição
nv.pred <- predict(nv.training,holdout)
print(nv.pred)
yhat <- yhat %>% cbind(nv.pred)
# comparação entre os dados e a predição
result = data.frame(holdout$dir_gap2, nv.pred)
print(result)

# cálculo das métricas
nv.metricas <- calculo_metricas(nv.pred,holdout$tipo) 

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


performance <- data.frame(logloss = c(tree.training$results[whichTwoPct,1:8][1,'logLoss'],
                                      rf.training$results[whichTwoPct,1:8][1,'logLoss'],
                                      boost.training$results[whichTwoPct,1:8][1,'logLoss'],
                                      knn.training$results[whichTwoPct,1:8][1,'logLoss'],
                                      nv.training$results[whichTwoPct,1:8][1,'logLoss']),
                          acuracia = c(tree.training$results[whichTwoPct,1:8][1,'Accuracy'],
                                       rf.training$results[whichTwoPct,1:8][1,'Accuracy'],
                                       boost.training$results[whichTwoPct,1:11][1,'Accuracy'],
                                       knn.training$results[whichTwoPct,1:8][1,'Accuracy'],
                                       nv.training$results[whichTwoPct,1:8][1,'Accuracy']),
                          row.names = c('tree', 'rf', 'boost', 'knn', 'nv')) %>% 
  kable(caption = 'Performance', digits = 3) %>% 
  kable_classic_2(full_width = F)

# o modelo com a melhor performance é o gradiente boosting
###################################
# Avaliar o XGBoosting            #
###################################
avalia <- function(modelo){
  #nome <- as.character(modelo)
  p_treino <- predict(modelo, training, type='prob') # Probabilidade predita
  c_treino <- predict(modelo, training)              # Classificação
  
  #Base de teste
  p_teste <- predict(modelo, holdout, type='prob')
  c_teste <- predict(modelo, holdout)
  
  # Data frame de avaliação (Treino)
  aval_treino <- data.frame(obs=training$tipo, 
                            pred=c_treino,
                            Y = p_treino[,2],
                            N = 1-p_treino[,2]
  )
  
  # Data frame de avaliação (Teste)
  aval_teste <- data.frame(obs=holdout$tipo, 
                           pred=c_teste,
                           Y = p_teste[,2],
                           N = 1-p_teste[,2]
  )
  
  tcs_treino <- caret::multiClassSummary(aval_treino, 
                                       lev=levels(aval_treino$obs))
  tcs_teste <- caret::multiClassSummary(aval_teste, 
                                      lev=levels(aval_teste$obs))
  
  
  print('Avaliação base de treino')
  print(tcs_treino)
  print('Avaliação base de teste')
  print(tcs_teste)
  dados <- data.frame(teste = tcs_teste,
    treino= tcs_treino) %>%
    tibble::rownames_to_column("VALUE") %>% t() %>% as.data.frame() %>%
    `colnames<-`(.[1, ]) %>% select(Accuracy) %>% filter(!row_number() %in% c(1)) %>% 
    tibble::rownames_to_column("VALUE") %>% dplyr::rename('Base' = 1)
  return(dados)
}

tree <- avalia(tree.training) 
gbm <- avalia(boost.training)
rf <- avalia(rf.training) 
knn <- avalia(knn.training) 
nv <- avalia(nv.training) 


purrr::pmap(list(tree.training,boost.training),avalia)


avaliacao <- purrr::reduce(list(tree,gbm,rf,knn,nv), 
                           dplyr::left_join, by = "Base")  %>% 
  dplyr::rename('tree' = 2, 'gbm' = 3, 'rf' = 4, 'knn' = 5, 'nv' = 6) %>%
  mutate_at(c('tree','gbm','rf','knn','nv'), as.numeric)


avaliacao %>% gather(key = Accuracy, value = Value, tree:nv) %>%
  ggplot(aes(x = Accuracy,y = Value, fill = Base)) + 
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

# random forest apresentou melhor acurácia para a classificação das junções 
# formadas por heteroestruturas de TMD

#-------------------------------------------------------------------------------
# Análise de Predição
# aplicação de métodos de machine learning para predizer o tipo de junção
#-------------------------------------------------------------------------------
# Modelos escolhidos:
# Modelo Linear por Míminos Quadrados (OLS), Mínimo Quadrado Parcial (PLS), 
# o Operador de Encolhimento em Seleção do Menor Absoluto (LASSO), 
# Gradient Boosting Regression (GBR), Gaussian Process Regression (GPR) e 
# Random Forest Regression (RFR).

# como não há heteroestruturas de tmds o suficiente para fazer um banco de dados 
# o bandgap será calculado de acordo com o modelo de anderson

df <-df %>% mutate(bandgap = case_when((as.numeric(tipo) == 2 & cbma > cbmb) ~ (cbmb - vbma),
                                       (as.numeric(tipo) == 2 & cbma < cbmb) ~ (cbma - vbmb),
                                       (as.numeric(tipo) == 1 & (vbma - cbma > vbmb - cbmb)) ~ (vbmb - cbmb),
                                       (as.numeric(tipo) == 1 & (vbma - cbma < vbmb - cbmb)) ~ (vbma - cbma),
                                       as.numeric(tipo) == 3 ~ 0)) 

df$bandgap <- abs(df$bandgap)


summary(df)

# gráfico de dispersão
ggplot(df, aes(x = TM_AN_V1, y = bandgap)) +
  geom_boxplot(color = "#39568CFF", size = 1) +
  labs(x = "Metal de Transição",
       y = "Bandgap") +
  scale_color_manual("Legenda:",
                     values = "grey50") +
  theme_classic()

ggplot(df, aes(x = C_AN_V1, y = bandgap)) +
  geom_boxplot(color = "#39568CFF", size = 1) +
  labs(x = "Calcogênio",
       y = "Bandgapo") +
  scale_color_manual("Legenda:",
                     values = "grey50") +
  theme_classic()

ggplot(df, aes(x = tipo, y = bandgap)) +
  geom_boxplot(color = "#39568CFF", size = 1) +
  labs(x = "Tipo",
       y = "Bandgap") +
  scale_color_manual("Legenda:",
                     values = "grey50") +
  theme_classic()










