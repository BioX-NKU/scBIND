library(scMET)
library(readr)

data <- read_table('intermediate.csv', col_names=TRUE, col_type='ccdd')
scmet_obj <- scmet(Y=data, seed=42)
hvf <- scmet_hvf(scmet_obj)
save_hvf <- hvf$hvf$summary$feature_name[hvf$hvf$summary$is_variable]
write.table(data.frame(save_hvf), file='hvf.csv', sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

