#-----------------------------------------------------------------------------
# Regressão do bandgap de heteroestruturasd de TMDs 
#-----------------------------------------------------------------------------
pacotes <- c('corrplot', 'caret', 'Hmisc','correlation', "tidyverse", 'gbm',
             'randomForest', "kableExtra", 'rpart', "rpart.plot", 'class',
             'gtools', 'corrgram','diffdf', 'MLmetrics', 'jtools',
             'xgboost', 'plyr','scales', 'ggplot2',"fastDummies",'nortest')

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
load('df.Rda')

# adicionar uma coluna com o bandgap das heterojunções calculado pela regra de anderson
df <-df %>% mutate(bandgap = case_when((as.numeric(tipo) == 2 & cbma > cbmb) ~ (cbmb - vbma),
                                       (as.numeric(tipo) == 2 & cbma < cbmb) ~ (cbma - vbmb),
                                       (as.numeric(tipo) == 1 & (vbma - cbma > vbmb - cbmb)) ~ (vbmb - cbmb),
                                       (as.numeric(tipo) == 1 & (vbma - cbma < vbmb - cbmb)) ~ (vbma - cbma),
                                       as.numeric(tipo) == 3 ~ 0)) 

df$bandgap <- abs(df$bandgap)

#-------------------------------------------------------------------------------
# Análise descritiva
#-------------------------------------------------------------------------------
summary(df)

# gráfico de dispersão (não ficou muito legal =/)
ggplot(df, aes(x = TM_AN_V1, y = bandgap)) +
  geom_boxplot(color = "#39568CFF", size = 1) +
  labs(x = "Metal de Transição",
       y = "Bandgap") +
  scale_color_manual("Legenda:",
                     values = "grey50") +
  theme_classic()

# boxplot da distribuição do bandgap para cada metal de transição 
# agora isso não fa muito sentido pois o que está sendo avaliado é a junção 
# dos tmds e não cada tmd em si
# a única variável que capta alguma característica da heterojunção é o tipo da junção
df %>% 
  ggplot(aes(x = TM_AN_V1,y = bandgap, fill = tipo)) + 
  geom_boxplot()+
  xlab('Metal de Transição') +
  ylab('Bandgap')

df %>% ggplot(aes(x = C_AN_V1,y = bandgap, fill = tipo)) + 
  geom_boxplot()+
  xlab('Calcogênio') +
  ylab('Bandgap')

# correlação 
# como o bandgap é criado a partir das bandas de valência e condução dos tmds
# nao faz muito sentido calcular a correlação entre eles
#df %>% correlation(method = 'pearson') %>% dplyr::filter(Parameter2 == 'bandgap') %>%
#  kable() %>%
#  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 28)

# teste anova para as variáveis categóricas
aov(bandgap ~ tipo + TM_AN_V1 + TM_AN_V2 + C_AN_V1 + C_AN_V2, df) %>% summary()
# todos as variáveis são estatisticamente significantes

# transformar em dummy
df_dummy <- fastDummies::dummy_cols(df, select_columns = c('TM_AN_V1','TM_AN_V2','C_AN_V1','C_AN_V2','tipo'),
                                       remove_selected_columns = T,
                                       remove_most_frequent_dummy = T) #categoria de referência

# separar a amostra em 30% para teste e 70% para treino
set.seed(123)
sample_size = floor(0.7*nrow(df_dummy))
picked = sample(seq_len(nrow(df_dummy)),size = sample_size)

training = df_dummy[picked,]
holdout = df_dummy[-picked,]
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

# generalized linear model (não deu certo)
glm_model <- caret::train(bandgap ~., data = training, 
                          method = "glm", 
                          trControl = fitControl, 
                          verbose = FALSE, 
                          #tuneGrid = treeGrid,
                          ## Specify which metric to optimize
                          metric = 'RMSE')

# generalized linear model with step wise
library('MASS')
glmStepAIC_model <- caret::train(bandgap ~., data = training, 
                                 method = "glmStepAIC", 
                                 trControl = fitControl, 
                                 verbose = FALSE, 
                                 #tuneGrid = treeGrid,
                                 ## Specify which metric to optimize
                                 metric = 'RMSE')

# Boosted Generalized Linear Model(não deu certo)
glmboost_model <- caret::train(bandgap ~., data = training, 
                               method = "glmboost", 
                               trControl = fitControl, 
                               verbose = FALSE, 
                               #tuneGrid = treeGrid,
                               ## Specify which metric to optimize
                               metric = 'RMSE')
library(lmtest)
bptest(modelo_linear)
# data:  modelo_linear BP = 105.8, df = 32, p-value = 7.985e-10

# teste de aderência dos resíduos à normalidade
#Shapiro-Francia: n > 30
# lm
sf.test(lm_model$finalModel$residuals)
# glm com step wise
sf.test(glmStepAIC_model$finalModel$residuals)
# os resíduos são aderentes a normalidade

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

# como os resídios são aderentes a normalidade não é necessário fazer uma 
# transformação de box-cox

# procedimendo step-wise já foi feito com o modelo glmStepAIC

# intervalos
confint(modelo_linear_step, level = 0.95) # siginificância 5%

plot_summs(modelo_linear_step, colors = "#440154FF", scale = TRUE, plot.distributions = TRUE,
           inner_ci_level = .95) 

export_summs(lm_model,
             model.names = c("Modelo Linear"),
             scale = F, digits = 6)

# preditores mais importantes
varImp(lm_model, scale = F) %>% plot()
varImp(glmStepAIC_model, scale = F) %>% plot()

# resultados
lm_model$results
glmStepAIC_model$results

# comparação dos modelos
model_list <- list(lm = lm_model, glm = glmStepAIC_model)
res <- resamples(model_list)
summary(res)

# testar de um modelo é estatisticamente diferrente do outro
compare_models(lm_model, glmStepAIC_model)
# p-value = 0.5128, os modelos não são estatisticamente diferentes

# cálculo da predição e métricas de performance
lm_model_pred = predict(lm_model, holdout,type = 'raw')
glmStepAIC_model_pred = predict(glmStepAIC_model,holdout)

#-------------------------------------------------------------------------------
#    Mínimo Quadrado Parcial (PLS) ok
#-------------------------------------------------------------------------------
#install.packages('pls')
library('pls')

# grade para encontrar os melhores parâmetros
seed(123)
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

# variáveis importantes
plot(varImp(pls_model))

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
lasso_model$results
lasso_model$finalModel
lasso_model$bestTune
plot(lasso_model)

# variáveis importantes
plot(varImp(lasso_model))

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
plot(varImp(gbm_model))

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
Xg_model$finalModel
Xg_model$bestTune
plot(Xg_model)

# variáveis importantes
plot(varImp(Xg_model))

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
#plot(gpr_model)

# variáveis importantes
plot(varImp(gpr_model))

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
plot(varImp(rf_model))

# previsão
rf_pred <- predict(rf_model,holdout)
#-------------------------------------------------------------------------------
# Avaliação dos modelos
#-------------------------------------------------------------------------------
pred <- list(glmStepAIC_model_pred, pls_pred, lasso_pred, gbm_pred, Xg_pred, gpr_pred,
     rf_pred)

avaliacao <- function(prediction){
  postResample(prediction,holdout$bandgap)
}

res = map(pred,avaliacao) %>% unnest()

data.frame(modelo = c('glm', 'pls', 'lasso', 'gmb', 'xb', 'gpr', 'rf'),
           res)

tabela <- data.frame(glm = res[1],
           pls = res[2],
           lasso = res[3],
           gbm = res[4],
           xb = res[5],
           gpr = res[6],
           rf = res[7]) %>% dplyr::rename(glm=1, pls=2, lasso=3, gbm=4, xb=5,
                                          gpr=6, rf=7) %>% rownames_to_column() %>% t()
 
colnames(tabela) <- tabela[1,]
tabela <- tabela[-1,]
r_quad <- tabela %>%as.data.frame() %>% rownames_to_column('modelo') %>% 
  ggplot(aes(x = modelo,y = Rsquared)) +
  geom_col()
rmse <- tabela %>%as.data.frame() %>% rownames_to_column('modelo') %>% 
  ggplot(aes(x = modelo,y = RMSE)) +
  geom_col()

rmse + r_quad + plot_layout(guides = "collect")

tabela %>%as.data.frame() %>% rownames_to_column('modelo') %>% kable() %>%
  kable_styling(bootstrap_options = "striped",
                full_width = F,
                font_size = 28)

# escrever a equação final bonitinha
library(equatiomatic)
extract_eq(modelo_linear_step, use_coefs = T) %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped",
                full_width = F,
                font_size = 28)




