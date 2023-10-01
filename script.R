#-----------------------------------------------------------------------------
#         construção do banco de dados - data wrangling
#-----------------------------------------------------------------------------

#pacotes utilizados

pacotes <- c('DBI', "tidyverse", "rvest", 'PeriodicTable', 'jtools', 
             "kableExtra", 'RSQLite', "magick", "webshot", 'ggplot2',
             'gridExtra', 'ggpmisc', 'cowplot', 'gtable')

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

#conectar com o banco de dados
sqlite.driver <- dbDriver('SQLite')
db <- dbConnect(SQLite(), 'c2dm.db')

data(periodicTable)

# ver a lista de tabelas no bando de dados
dbListTables(db)

keys<- dbReadTable(db, 'keys')
information <- dbReadTable(db, 'information')
species<- dbReadTable(db, 'species')
number_key_values<- dbReadTable(db, 'number_key_values')
sqlite_sequence<- dbReadTable(db, 'sqlite_sequence')
systems <- dbReadTable(db, 'systems')
text_key_values <- dbReadTable(db, 'text_key_values')

# contruindo o banco de dados geral df1 no formato tibble

# adicionando as colunas id, magmom e mass da tabela systems e transformando
# o paramentro hform da tabela number_key_values em coluna 
df1 <- systems %>% select(id,magmom, mass) %>% 
  left_join(number_key_values %>% filter(key == "hform") %>% select(value,id)) %>%
  rename("hform" = value)
  
# transformando os paramentros dir_gap, phase e a da tabela number_key_values 
# em coluna 
df1 <- df1 %>% left_join(number_key_values %>% filter(key == "dir_gap", keep = TRUE) %>% 
                           select(value,id)) %>% rename("dir_gap" = value)

df1 <- df1 %>% left_join(text_key_values %>% filter(key == "phase", keep = TRUE) %>% 
                           select(value,id)) %>% rename("phase" = value)

df1 <-df1 %>% left_join(number_key_values %>% filter(key == "a", keep = TRUE) %>% 
                           select(value,id)) %>% rename("a" = value)

# transformando os paramentros Z da tabela number_key_values nas colunas 
# TM_AN (número atômico do metal de transição) e C_AN (número atômico do 
# calcogênio)
df1 <- df1 %>% left_join(species %>% filter(Z != 8 ) %>% filter(Z != 34 ) %>% 
                           filter(Z != 52) %>% filter(Z != 16) %>% 
                           select(Z,id)) %>% rename("TM_AN" = Z)

df1 <- df1 %>% left_join(species %>% filter(Z %in% c(8,52,34, 16)) %>% 
                           select(Z,id)) %>% rename("C_AN" = Z)

# adicionando os valores AW (atomic weight) para metais de transição (TM)
# e calcogênios (C). 
url <- "https://www.degruyter.com/document/doi/10.1515/pac-2019-0603/html"
atomic_weight_url <- read_html(url)
atomic_weight_data <- 
  html_nodes(atomic_weight_url, "table") %>% html_table() %>% .[[3]]  
  
colnames(atomic_weight_data) <- as.character(atomic_weight_data[2,])
atomic_weight_data <- atomic_weight_data[-1:-3,] 
colnames(atomic_weight_data)[7] <- 'atomic_weight'
colnames(atomic_weight_data)[3] <- 'Z'
atomic_weight_data <- atomic_weight_data %>% select(Z, atomic_weight)
atomic_weight_data$Z <- as.numeric(atomic_weight_data$Z)
atomic_weight_data$atomic_weight <- as.numeric(atomic_weight_data$atomic_weight)

df1 <- df1 %>% left_join(atomic_weight_data %>% rename(TM_AN = "Z", TM_AW ="atomic_weight") %>% 
                           select(TM_AN,TM_AW))
df1 <- df1 %>% left_join(atomic_weight_data %>% rename(C_AN = "Z", C_AW ="atomic_weight") %>% 
                           select(C_AN, C_AW))

# group number (GN)
df1 <- df1 %>% mutate(C_GN = 16, .keep = "all")
df1 <- df1 %>% mutate(TM_GN = if_else(TM_AN %in% c(21,29),3,
                              if_else(TM_AN %in% c(22,40,72,104),4,
                              if_else(TM_AN %in% c(23,41,73,105),5,
                              if_else(TM_AN %in% c(24,42,74,106),6,
                              if_else(TM_AN %in% c(25,43,75,107),7,
                              if_else(TM_AN %in% c(26,44,76,107),8,
                              if_else(TM_AN %in% c(27,45,77),9,
                              if_else(TM_AN %in% c(28,46,78),10,
                              if_else(TM_AN %in% c(29,47,79),11,
                              if_else(TM_AN %in% c(30,48,80,112),12,
                              if_else(TM_AN %in% c(13,31,49,81),13,
                              if_else(TM_AN %in% c(50,82,32),14,0)))))))))))), .keep = "all")
#period number
df1 <- df1 %>% mutate(TM_PN = if_else(TM_AN %in% 19:36,4,
                                      if_else(TM_AN %in% 37:54,5,
                                              if_else(TM_AN %in% 55:86,6,7))),.keep = "all")

df1 <- df1 %>% mutate(C_PN = if_else(C_AN %in% 8,2,
                                     if_else(C_AN %in% 16,4,
                                             if_else(C_AN %in% 34,5,6))),.keep = "all")

# atomic density D = m / V
periodicTable <- periodicTable %>% rename(C_AN=numb)
df1 <- df1 %>% left_join(periodicTable %>% select(C_AN,density) %>% 
                           filter(C_AN %in% c(8,16,34,52))) %>% 
                           rename(C_AD = density)

periodicTable <- periodicTable%>%rename(TM_AN=C_AN)
df1 <- df1 %>% left_join(periodicTable %>% select(TM_AN,density)) %>% 
                           rename(TM_AD = density)


# atomic radius (AR)
url <- "https://www.chemeurope.com/en/encyclopedia/Atomic_radii_of_the_elements_%28data_page%29.html"
atomic_radius_url <- read_html(url)
atomic_radius_data <- 
  html_nodes(atomic_radius_url, "table") %>% html_table() %>% .[[3]]


df1 <- df1 %>% left_join(atomic_radius_data %>% 
                               rename(TM_AR = 'calculated', TM_AN = 'number') %>%
                             select(TM_AR,TM_AN))

df1 <- df1 %>% left_join(atomic_radius_data %>% 
                           rename(C_AR = 'calculated', C_AN = 'number') %>%
                           select(C_AR,C_AN))

# van der Waals radio (vDW_radio)
periodicTable <- periodicTable%>%rename(C_AN = TM_AN)

df1 <- df1 %>% left_join(periodicTable %>% select(C_AN,rvdw)) %>% rename(C_VDW = rvdw)

periodicTable <- periodicTable%>%rename(TM_AN=C_AN)
df1 <- df1 %>% left_join(periodicTable %>% select(TM_AN,rvdw)) %>% rename(TM_VDW = rvdw)

# atomic volume V = m / D
df1$TM_AW <- as.numeric(df1$TM_AW) # AW: atomic weight
df1$C_AW <- as.numeric(df1$C_AW)

df1 <- df1 %>% mutate(TM_AV= TM_AW / TM_AD, .keep = "all")
df1 <- df1 %>% mutate(C_AV= C_AW / C_AD, .keep = "all")

# raio covalente do calcogenio
periodicTable <- periodicTable%>%rename(C_AN = TM_AN)
df1 <- df1 %>% left_join(periodicTable %>% select(C_AN,rcov)) %>% rename(C_COV = rcov)

# TM fusion heat
url <- "https://www.schoolmykids.com/learn/periodic-table/heat-of-fusion-of-all-the-elements"
fusion_heat_url <- read_html(url)
fusion_heat_data <- html_element(fusion_heat_url, "table") %>%
  html_table()
df1 <- df1 %>% left_join(fusion_heat_data %>% rename(TM_AN = "Element Atomic Number", TM_FH ="Element Heat of Fusion") %>% 
  select(TM_AN,TM_FH))

# TM specifc heat
periodicTable <- periodicTable%>%rename(TM_AN=C_AN)
df1 <- df1 %>% left_join(periodicTable %>% select(TM_AN,C)) %>% rename(TM_SH = C)

# melting point
df1 <- df1 %>% left_join(periodicTable %>% select(TM_AN,melting)) %>% rename(TM_MP = melting)
periodicTable <- periodicTable%>%rename(C_AN = TM_AN)
df1 <- df1 %>% left_join(periodicTable %>% select(C_AN,melting)) %>% rename(C_MP = melting)

# boiling point
df1 <- df1 %>% left_join(periodicTable %>% select(C_AN,boiling)) %>% rename(C_BP = boiling)
periodicTable <- periodicTable%>%rename(TM_AN=C_AN)
df1 <- df1 %>% left_join(periodicTable %>% select(TM_AN,boiling)) %>% rename(TM_BP = boiling)

# first ionization energy
df1 <- df1 %>% left_join(periodicTable %>% select(TM_AN,IP)) %>% rename(TM_FIE = IP)
periodicTable <- periodicTable%>%rename(C_AN = TM_AN)
df1 <- df1 %>% left_join(periodicTable %>% select(C_AN,IP)) %>% rename(C_FIE = IP)

# second ionization energy


# Electronegativity
df1 <- df1 %>% left_join(periodicTable %>% select(C_AN,Eneg)) %>% rename(C_EN = Eneg)
periodicTable <- periodicTable%>%rename(TM_AN=C_AN)
df1 <- df1 %>% left_join(periodicTable %>% select(TM_AN,Eneg)) %>% rename(TM_EN = Eneg)

# salvar o database df1
save(df1, file = 'df1.Rda')

# carregar o banco de dados
load('df1.Rda')

# análise do banco de dados inicial
colnames(df1)[1]


# database final com os parametros importantes
# C_GN tem só um valor, ao calcular a correlação o desvio padrão é zero
df2 <- df1 %>%  select(-id, -C_GN) %>% drop_na() %>% filter(C_AN != 8) 

# transformar os parametros para o tipo correto

df2 <- df2 %>% mutate_at(c('TM_FH', 'TM_AR', 'C_AR'), as.numeric)

df2 <- df2 %>% mutate_at(c('phase', 'C_AN', 'TM_AN', 'C_PN', 'TM_PN', 'TM_GN'),
                         as.factor)

summary(df2)

#-------------------------------------------------------------------------------
# análise descritiva do banco de dados
#-------------------------------------------------------------------------------

# tabela resumo de algumas estatísticas descritivas de acordo com cada calcogênio
tb2 <- df2 %>% select(C_AN,phase,dir_gap) %>% 
  group_by(C_AN) %>% dplyr::rename('Calcogênio' = C_AN) %>% 
  dplyr::summarise(N = n(),Min = round(min(dir_gap),3), Máx = max(dir_gap),Média = mean(dir_gap), 
            DP = sd(dir_gap)) %>% kable() %>% kable_classic_2(full_width = F)

df2 %>% select(TM_AN,phase,dir_gap) %>% 
  group_by(TM_AN) %>% dplyr::rename('Calcogênio' = TM_AN) %>% 
  dplyr::summarise(N = n(),Min = round(min(dir_gap),3), Máx = max(dir_gap),Média = mean(dir_gap), 
            DP = sd(dir_gap)) %>% kable() %>% kable_classic_2(full_width = F)

# lista de todas as variáveis no banco de dados final
ls(df2)

# theme da tabela
tb <- sub(".*_", "",colnames(df2)) %>% as.data.frame(sigla = c()) %>% unique() %>% 
  mutate(Descrição = c('Momento Magnético', 'Massa', 
                       'Energia de Formação', 'Bandgap Direto', 'Fase', 
                       'Constante de Rede', 'Número Atômico',
                       'Peso Atômico', 'Número do Grupo', 'Número do Período',
                       'Densidade Atômica', 'Raio Atômico', 'Raio de van der Waals',
                       'Volume Atômico', 'Raio Covalente', 'Temperatura de Fusão', 
                       'Calor Específico', 'Ponto de Fusão','Ponto de Ebulição',
                       'Primeira Energia de Ionização', 'Eletronegatividade'))  %>% 
  tableGrob(rows = NULL, theme = ttheme_minimal(base_size = 11))

#kable(row.names = F) %>% kable_classic_2(full_width = F) 

# estimativa da densidade de kernel
kernel <- ggplot(df2, aes(x = C_AN, y = dir_gap, fill = phase)) + 
  geom_violin() +
  labs(x = 'Número Atômico', y = 'Bandgap Direto (eV)', fill = 'Fase') +
  theme(#text = element_text(size = 12),
        panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA),
        legend.position="none")

# histograma com a distribuição dos tmds pelo bandgap
# a grande maioria dos tmds possui bandgap nulo nesta base
hist <- ggplot(df2, aes(x = dir_gap, fill = phase )) + 
  geom_histogram(position = "dodge") +
  labs(x = 'Bandgap Direto (eV)', y = 'Contagem', fill = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA),
        legend.position="none")

# magmom x dir_gap 
# o momento margético parece ser um fator importante para separar 
# os tmds de acordo com o bandgap, independente o metal de transição
md <- ggplot(df2, aes(x=magmom, y = dir_gap, color = TM_GN, shape = phase)) + 
  geom_point() +
  labs(x = 'Momento Magnético', y = 'Bandgap (eV)', color = 'TM_GN', 
       shape = 'Fase') +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA),
        legend.position="none")

legend <- get_legend(hist + theme(legend.position = "left"))

p <- plot_grid(hist, kernel, legend, rel_widths = c(0.8,0.8,0.3),
          nrow = 1, labels = c('a', 'b'))
p2 <- plot_grid(p, md, ncol = 1, labels =  c('','c'))

png(file="C:/Users/Gabriela/Desktop/MBA - USP/TCC/Banco de dados/banco_de_dados.png",
    width=800 , height=450, units = "px")
grid.arrange(arrangeGrob(tb, p2, ncol = 2, widths=c(0.8,2.2)))
dev.off()


