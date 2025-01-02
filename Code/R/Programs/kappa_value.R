#561
confusion_matrix <- table(Shannon = as.vector(binary_map_shannon), Rényi = as.vector(binary_map_renyi))
print(confusion_matrix)


library(irr)  
kappa_value <- kappa2(data.frame(binary_map_shannon = as.vector(binary_map_shannon), 
                                 binary_map_renyi = as.vector(binary_map_renyi)))
cat("Kappa Coefficient:", kappa_value$value, "\n")

# 
discrepancy_map <- binary_map_shannon != binary_map_renyi
image(discrepancy_map, main = "Discrepancy Map (Shannon vs Rényi)",
      col = c("white", "red"), xlab = "Columns", ylab = "Rows")
