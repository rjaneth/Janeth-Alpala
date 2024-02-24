
load("./Data/results_data_Ottawa_cqv_512_5.Rdata")
x1 <- as.vector(cv_values)

ruta_del_archivo <- "./Data/Ottawa_cqv_5.csv"


write.csv(data.frame(x1), file = "./Data/Ottawa_cqv_5.csv", row.names = FALSE)
