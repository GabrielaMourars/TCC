#-----------------------------------------------------------------------------
# Regressão do bandgap de heteroestruturasd de TMDs 
#-----------------------------------------------------------------------------
pacotes <- c('corrplot', 'caret', 'Hmisc','correlation', "tidyverse", 'gbm',
             'randomForest', "kableExtra", 'rpart', "rpart.plot", 'class',
             'gtools', 'corrgram','diffdf', 'MLmetrics', 'jtools',
             'xgboost', 'plyr','scales', 'ggplot2',"fastDummies",'nortest',
             'gridExtra','grid','cowplot')
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
# como não há heteroestruturas de tmds o suficiente para fazer um banco de dados 
# o bandgap será calculado de acordo com o modelo de anderson

# carregar o banco de dados
load('data.Rda')

data <- data %>% mutate(bandgap = case_when((as.numeric(tipo) == 2 & cbm.x > cbm.y) ~ (cbm.y - vbm.x),
                                            (as.numeric(tipo) == 2 & cbm.x < cbm.y) ~ (cbm.x - vbm.y),
                                            (as.numeric(tipo) == 1 & (vbm.x - cbm.x > vbm.y - cbm.y)) ~ (vbm.y - cbm.y),
                                            (as.numeric(tipo) == 1 & (vbm.x - cbm.x < vbm.y - cbm.y)) ~ (vbm.x - cbm.x),
                                            (as.numeric(tipo) == 3) ~ 0)) 

data$bandgap <- abs(data$bandgap)
data <- data %>% select(-tipo)

summary(data)

# transformar em dummy
data_dummy <- fastDummies::dummy_cols(data, select_columns = c('TM_AN_V1','TM_AN_V2','C_AN_V1','C_AN_V2','tipo'),
                                       remove_selected_columns = T,
                                       remove_most_frequent_dummy = T) #categoria de referência

# separar a amostra em 30% para teste e 70% para treino
set.seed(123)
sample_size = floor(0.7*nrow(data_dummy))
picked = sample(seq_len(nrow(data_dummy)),size = sample_size)

training = data_dummy[picked,]
holdout = data_dummy[-picked,]

sample_size = floor(0.7*nrow(data))
picked = sample(seq_len(nrow(data)),size = sample_size)

training = data[picked,]
holdout = data[-picked,]
#-------------------------------------------------------------------------------
# Análise de Predição
# aplicação de métodos de machine learning para predizer o tipo de junção
#-------------------------------------------------------------------------------
# Modelos escolhidos:
# - Modelo Linear por Míminos Quadrados (OLS)
# - Mínimo Quadrado Parcial (PLS)
# - Operador de Encolhimento em Seleção do Menor Absoluto (LASSO)
# - Gradient Boosting Regression (GBR)
# - XgBoost (XGBoost)
# - Gaussian Process Regression (GPR)
# - Random Forest Regression (RFR).

# a performance dos modelos será avaliada pela métrica RMSE

fitControl <- caret::trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 5,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = defaultSummary
)

#-------------------------------------------------------------------------------
#    Modelo linear ok
#-------------------------------------------------------------------------------
# modelo
set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
# linear model
lm_model <- caret::train(bandgap ~., data = training, 
                              method = "lm", 
                              trControl = fitControl, 
                              verbose = FALSE, 
                              #tuneGrid = treeGrid,
                              ## Specify which metric to optimize
                              metric = 'RMSE')

summary(lm_model)




# generalized linear model with step wise
library('MASS')
glmStepAIC_model <- caret::train(bandgap ~., data = training, 
                                 method = "glmStepAIC", 
                                 trControl = fitControl, 
                                 verbose = FALSE, 
                                 #tuneGrid = treeGrid,
                                 ## Specify which metric to optimize
                                 metric = 'RMSE')
summary(glmStepAIC_model)
glmStepAIC_model$results
# teste de aderência dos resíduos à normalidade
#Shapiro-Francia: n > 30
# p-valor > 0.05 -> a distribuição dos dados não é significativamente diferentes 
# de uma distribuição normal
# lm
sf.test(lm_model$finalModel$residuals) 
# glm com step wise
sf.test(glmStepAIC_model$finalModel$residuals)
# os resíduos são aderentes a normalidade
sf.test(data$bandgap)

#Histograma dos resíduos do modelo linear (não consegui fazer o gráfico)
lm_model %>%
  mutate(residuos = finalModel$residuals) %>%
  ggplot(aes(x = finalModel$residuals)) +
  geom_histogram(aes(y = ..density..), 
                 color = "grey50", 
                 fill = "grey90", 
                 bins = 30,
                 alpha = 0.6) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(modelo_linear$residuals),
                            sd = sd(modelo_linear$residuals)),
                aes(color = "Curva Normal Teórica"),
                size = 2) +
  scale_color_manual("Legenda:",
                     values = "#FDE725FF") +
  labs(x = "Resíduos",
       y = "Frequência") +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA),
        legend.position = "bottom")

# como os resídios não são aderentes a normalidade. é necessário fazer uma 
# transformação de box-cox

#Para calcular o lambda de Box-Cox
library("car")
data_box_cox <- data %>% mutate(bandgap = bandgap +1)
lambda_BC <- powerTransform(data_box_cox$bandgap)
lambda_BC

#Inserindo o lambda de Box-Cox na base de dados para a estimação de um novo modelo
data_box_cox$bandgap <- (((data_box_cox$bandgap ^ lambda_BC$lambda) - 1) / 
                           lambda_BC$lambda)

training_box = data_box_cox[picked,]
holdout_box = data_box_cox[-picked,]

#Estimando um novo modelo OLS com variável dependente transformada por Box-Cox
lm_model_box <- caret::train(bandgap ~., data = training_box, 
                         method = "lm", 
                         trControl = fitControl, 
                         verbose = FALSE, 
                         #tuneGrid = treeGrid,
                         ## Specify which metric to optimize
                         metric = 'RMSE') 

summary(lm_model_box)
#teste shapiro-francia
sf.test(lm_model_box$finalModel$residuals)
# continua não normal

# modelo linear com step-wise
glmStepAIC_model_box <- caret::train(bandgap ~., data = training_box, 
                                 method = "glmStepAIC", 
                                 trControl = fitControl, 
                                 verbose = FALSE, 
                                 #tuneGrid = treeGrid,
                                 ## Specify which metric to optimize
                                 metric = 'RMSE')
summary(glmStepAIC_model_box)
glmStepAIC_model_box$results
sf.test(glmStepAIC_model_box$finalModel$residuals)
#continua não normal

# teste Breusch-Pagan para heterocedasticidade 
training_box %>%
  mutate(residuos = lm_model_box$finalModel$residuals) %>% 
  ggplot(data = ., aes(y = residuos, x = bandgap)) + 
  geom_point() + 
  geom_abline(slope = 0) +
  theme_classic()

training %>%
  mutate(residuos = glmStepAIC_model$finalModel$residuals) %>% 
  ggplot(data = ., aes(y = residuos, x = bandgap)) + 
  geom_point() + 
  geom_abline(slope = 0) +
  theme_classic()

training %>%
  ggplot(aes(x = predict(glmStepAIC_model,training), y = bandgap))+
  geom_point()+
  scale_y_log10()


data %>% ggplot(aes(x = sqrt(bandgap)))+
  geom_histogram(binwidth = 0.05)

data %>% ggplot(aes(x = bandgap))+
  geom_histogram(binwidth = 0.05) 


ggplot(data = training, aes(sample = glmStepAIC_model$finalModel$residuals)) + 
  geom_qq() + 
  geom_qq_line() + 
  labs( x = 'Theoretical', y = 'Sample Quantiles - a')

var_func <- lm(residuos^2 ~ tx_homicidio, data = msp)
summary(var_func)


# intervalos
confint(lm_model, level = 0.95) # siginificância 5%

plot_summs(lm_model, colors = "#440154FF", scale = TRUE, plot.distributions = TRUE,
           inner_ci_level = .95) 

export_summs(lm_model,
             model.names = c("Modelo Linear"),
             scale = F, digits = 6)

# preditores mais importantes
varImp(lm_model, scale = F) %>% plot(top = 10)
varImp(glmStepAIC_model, scale = F) %>% plot(top = 10)

# resultados
lm_model$results
plot(lm_model_box)
lm_model_box$results
glmStepAIC_model$results
glmStepAIC_model_box$results

# comparação dos modelos
model_list <- list(lm = lm_model, glm = glmStepAIC_model)
res <- resamples(model_list)
summary(res)

# testar se um modelo é estatisticamente diferrente do outro
compare_models(lm_model, glmStepAIC_model)
# p-value = 0.5128, os modelos não são estatisticamente diferentes

# cálculo da predição e métricas de performance
lm_model_pred = predict(lm_model, holdout)
glmStepAIC_model_pred = predict(glmStepAIC_model,holdout)

#-------------------------------------------------------------------------------
#    Mínimo Quadrado Parcial (PLS) ok
#-------------------------------------------------------------------------------
#install.packages('pls')
library('pls')

# grade para encontrar os melhores parâmetros
set.seed(123)
plsGrid <- expand.grid(ncomp = c(1:10))

pls_model <- caret::train(bandgap ~., data = training, 
                               method = "pls", 
                               trControl = fitControl, 
                               verbose = FALSE, 
                               tuneGrid = plsGrid,
                               ## Specify which metric to optimize
                               metric = 'RMSE')
summary(pls_model)

plot(pls_model)
sf.test(pls_model)
pls_model$results
pls_model$bestTune
# variáveis importantes
plot(varImp(pls_model), top = 10)

# previsão
pls_pred <- predict(pls_model,holdout)

# melhor ajuste: ncomp = 6
pls_model$bestTune
#-------------------------------------------------------------------------------
#    Operador de Encolhimento em Seleção do Menor Absoluto (LASSO) ok
#-------------------------------------------------------------------------------
library('elasticnet')
lassoGrid <- expand.grid(fraction = seq(0,1,by=0.01))

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
lasso_model <- train(bandgap ~., data = training, 
                   method = "lasso", 
                   trControl = fitControl, 
                   #verbose = FALSE, 
                   #tuneGrid = lassoGrid,
                   ## Specify which metric to optimize
                   metric = "RMSE")

summary(lasso_model)
plot(lasso_model)
lasso_model$results
lasso_model$bestTune
# variáveis importantes
plot(varImp(lasso_model), top = 10)

# previsão
lasso_pred <- predict(lasso_model,holdout)

#-------------------------------------------------------------------------------
#    Gradient Boosting Regression (GBR) ok
#-------------------------------------------------------------------------------
gbmGrid <- expand.grid(interaction.depth = c(1, 5, 9), 
                       n.trees = (1:30)*50, 
                       shrinkage = c(0.1, 0.3, 0.001),
                       n.minobsinnode = c(10,15,20))
set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
gbm_model <- train(bandgap ~., data = training, 
                  method = "gbm", 
                  trControl = fitControl, 
                  verbose = FALSE, 
                  tuneGrid = gbmGrid,
                  ## Specify which metric to optimize
                  metric = "RMSE")
gbm_model
summary(gbm_model)
gbm_model$results
gbm_model$finalModel
gbm_model$bestTune
plot(gbm_model)

# variáveis importantes
plot(varImp(gbm_model), top = 10)

# previsão
gbm_pred <- predict(gbm_model,holdout)
#-------------------------------------------------------------------------------
#    XgBoost ok
#-------------------------------------------------------------------------------
XgGrid <- expand.grid(
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
Xg_model <- train(bandgap ~., data = training, 
                        method = "xgbTree", 
                        trControl = fitControl, 
                        verbose = FALSE, 
                        tuneGrid = XgGrid,
                        ## Specify which metric to optimize
                        metric = "RMSE")

summary(Xg_model)
Xg_model$results
Xg_model$bestTune
plot(Xg_model)

# variáveis importantes
plot(varImp(Xg_model), top = 10)

# previsão
Xg_pred <- predict(Xg_model,holdout)
#-------------------------------------------------------------------------------
#    Gaussian Process Regression (GPR) ok
#-------------------------------------------------------------------------------
library('kernlab')
# no grid
set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
gpr_model <- caret::train(bandgap ~., data = training, 
                         method = "gaussprLinear", 
                         trControl = fitControl, 
                         verbose = FALSE, 
                         #tuneGrid = rfGrid,
                         ## Specify which metric to optimize
                         metric = "RMSE")
summary(gpr_model)
gpr_model$results
gpr_model$finalModel
plot(gpr_model)

# variáveis importantes
plot(varImp(gpr_model), top = 10)

# previsão
gpr_pred <- predict(gpr_model,holdout)

#-------------------------------------------------------------------------------
#    Random Forest Regression (RFR) ok
#-------------------------------------------------------------------------------
# grade para encontrar os melhores parâmetros
rfGrid <- expand.grid(mtry = c(1:10))

set.seed(123)
# ajuste do modelo no data set de treino para encontrar os melhores parâmetros
rf_model <- caret::train(bandgap ~., data = training, 
                            method = "rf", 
                            trControl = fitControl, 
                            verbose = FALSE, 
                            tuneGrid = rfGrid,
                            ## Specify which metric to optimize
                            metric = "RMSE")
summary(rf_model)
rf_model$results
rf_model$finalModel
rf_model$bestTune
plot(rf_model)

# variáveis importantes
plot(varImp(rf_model), top = 10)

# previsão
rf_pred <- predict(rf_model,holdout)
#-------------------------------------------------------------------------------
# Avaliação dos modelos
#-------------------------------------------------------------------------------

# comparação entre base de teste e base de treino
nome_modelo <- function(modelo){
  nome <- modelo$method
  return(nome)
}

avalia <- function(modelo){
  #nome <- as.character(modelo)
  c_treino <- predict(modelo, training)              
  
  #Base de teste
  c_teste <- predict(modelo, holdout)
  
  # Data frame de avaliação (Treino)
  aval_treino <- data.frame(obs=training$bandgap, 
                            pred=c_treino
                            )
  
  # Data frame de avaliação (Teste)
  aval_teste <- data.frame(obs=holdout$bandgap, 
                           pred=c_teste
                           )
  
  tcs_treino <- caret::defaultSummary(aval_treino, 
                                         lev=levels(aval_treino$obs))
  tcs_teste <- caret::defaultSummary(aval_teste, 
                                        lev=levels(aval_teste$obs))
  
  
  print('Avaliação base de treino')
  print(tcs_treino)
  print('Avaliação base de teste')
  print(tcs_teste)
  nome <- nome_modelo(modelo)
  dados <- data.frame(#teste = tcs_teste
    treino = tcs_treino
  ) %>%`colnames<-`(nome) %>%
    t() %>% as.data.frame() %>% select(RMSE, Rsquared)
  return(dados)
}

avalia_box <- function(modelo){
  
  c_treino <- predict(modelo, training_box)              
  
  #Base de teste
  c_teste <- predict(modelo, holdout_box)
  
  # Data frame de avaliação (Treino)
  aval_treino <- data.frame(obs=training_box$bandgap, 
                            pred=c_treino
  )
  
  # Data frame de avaliação (Teste)
  aval_teste <- data.frame(obs=holdout_box$bandgap, 
                           pred=c_teste
  )
  
  tcs_treino <- caret::defaultSummary(aval_treino, 
                                      lev=levels(aval_treino$obs))
  tcs_teste <- caret::defaultSummary(aval_teste, 
                                     lev=levels(aval_teste$obs))
  
  
  print('Avaliação base de treino')
  print(tcs_treino)
  print('Avaliação base de teste')
  print(tcs_teste)
  nome <- nome_modelo(modelo)
  dados <- data.frame(#teste = tcs_teste
    treino = tcs_treino
  ) %>%`colnames<-`(nome) %>%
    t() %>% as.data.frame() %>% select(RMSE, Rsquared) 
  
  return(dados)
}

modelos <- list(gpr_model, glmStepAIC_model, lasso_model, lm_model,
             pls_model, rf_model, Xg_model)

gbm <- avalia(gbm_model)
gpr <- avalia(gpr_model)
glm_box <- avalia_box(glmStepAIC_model_box)
lasso <- avalia(lasso_model)
lm_box <- avalia_box(lm_model_box)
pls <- avalia(pls_model)
rf <- avalia(rf_model)
Xg <- avalia(Xg_model)

aval_regre <- purrr::reduce(list(gpr,glm_box,lasso,lm_box,pls,rf,Xg,gbm), 
                           rbind) %>% rownames_to_column(var = 'Modelo') %>%
  mutate(Modelo = ifelse(Modelo == "gaussprLinear", "GPR", 
                         ifelse(Modelo == 'glmStepAIC', 'GLM', 
                                ifelse(Modelo == 'xgbTree', 'XgBoost',Modelo)))) %>%
  mutate_at(c('Modelo'), toupper)

# rais do erro quadrático médio na base de teste
rmse <-  ggplot(data = aval_regre, aes(x = reorder(Modelo, -RMSE),
                                       y = RMSE, fill = Modelo))+
  geom_col() +
  geom_text(aes(label = round(RMSE,3)), vjust = 1.5, colour = "white")+
  labs(x = 'Modelo', y = 'RMSE', title = 'Raiz do Erro Quadrático Médio na Base de Teste')+
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA),
        legend.position="none")
# r-quadrado na base de teste
ggplot(data = aval_regre, aes(x = reorder(Modelo, -Rsquared), y = Rsquared, fill = Modelo))+
  geom_col() +
  geom_text(aes(label = round(Rsquared,3)), vjust = 1.5, colour = "white")+
  labs(x = 'Modelo', y = 'R-quadrado')+
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA),
        legend.position="none")

# GBM apresentou melhor RMSE para a regressão do bandgap  

# variáveis importantes do GBM
library('ggpmisc')
var_gbm <- ggplot(varImp(gbm_model, scale = F), top = 5  )+
  labs(x = 'Variável', y = 'Importância', title ='GBM') + 
  annotate(geom = "table", x = 0, y = 40, label = list(gbm_model$bestTune))
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA)
        #plot.title = element_text(hjust = 0.5)
  )

# Put the figure and table together:
final_figure <- cowplot::plot_grid(rmse,var_gbm, labels = c('a)', 'b)'),
                                   ncol = 2,
                                   rel_heights = c(1, 1))






