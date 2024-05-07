
#load("./Data/results_data_Ottawa_cqv_512_5.Rdata")
#load("./Data/results_data_simulated_z_up_le_mnad_250_7.Rdata")
load("../Programs/Data/results_data_Frankfurt_512_mnad_7.Rdata")
x1 <- as.vector(cd_values_mnad)
#x1 <- as.vector(cv_values)
ruta_del_archivo <- "./Data/Frankfurt_512_mnad_7.csv"


write.csv(data.frame(x1), file = "./Data/Frankfurt_512_mnad_7.csv", row.names = FALSE)
