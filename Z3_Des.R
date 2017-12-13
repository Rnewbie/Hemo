library(seqinr)
fasta <- read.fasta("Hb-alpha.fasta")
df <- rapply(fasta, c)
df <- data.frame(df)
df <- t(df)
df_capital <- toupper(df)

mutated_position <- c(1, 2, 6, 13, 14, 16, 24, 26, 31, 34, 40, 41, 44, 45, 46, 47, 49, 50, 53, 55, 56, 57, 58, 59, 61, 63, 64, 65, 68, 
                      74, 75, 76, 77, 80, 81, 84, 85, 86, 87, 88, 89, 90, 92, 94, 95, 96, 97, 110, 112, 119, 126, 130, 133, 134, 135, 138, 139, 140, 141)

alpha <- df_capital[mutated_position]
names(alpha) <- mutated_position
descs <- data.frame(AA1 = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"),
                    Z3_1 = c(0.24, 0.84, 3.98, 3.11, -4.22, 2.05, 2.47, -3.89, 2.29, -4.28, -2.85, 3.05, -1.66, 1.75, 3.52, 2.39, 0.75, -2.59, -4.36, -2.54),
                    Z3_2 = c(-2.32, -1.67, 0.93, 0.26, 1.94, 4.06, 1.95, -1.73, 0.89, -1.3, -0.22, 1.62, 0.27, 0.5, 2.5, -1.07, -2.18, -2.64, 3.94, 2.44),
                    Z3_3 = c(0.6, 3.71, 1.93, -0.11, 1.06, 0.36, 0.26, -1.71, -2.49, -1.49, 0.47, 1.04, 1.84, -1.44, -3.5, 1.15, -1.12, -1.54, 0.59, 0.43))

match_amino_acid <- function(AA) {
  d <- descs[descs$AA1 == AA, c(2, 3, 4)]
  return(d)
}


des <- sapply(alpha,match_amino_acid)
ColDesNames <- paste0("AA_Position", rep(colnames(des),  each = 3), "_",  rep(rownames(des), length(alpha)))
des <- unlist(des)
names(des) <- ColDesNames
des


alpha_z3 <- function(x) {
  fasta <- read.fasta("Hb-alpha.fasta")
  df <- rapply(fasta, c)
  df <- data.frame(df)
  df <- t(df)
  df_capital <- toupper(df)
  mutated_position <- c(1, 2, 6, 13, 14, 16, 24, 26, 31, 34, 40, 41, 44, 45, 46, 47, 49, 50, 53, 55, 56, 57, 58, 59, 61, 63, 64, 65, 68, 
                        74, 75, 76, 77, 80, 81, 84, 85, 86, 87, 88, 89, 90, 92, 94, 95, 96, 97, 110, 112, 119, 126, 130, 133, 134, 135, 138, 139, 140, 141)
  
  alpha <- df_capital[mutated_position]
  names(alpha) <- mutated_position
  descs <- data.frame(AA1 = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"),
                      Z3_1 = c(0.24, 0.84, 3.98, 3.11, -4.22, 2.05, 2.47, -3.89, 2.29, -4.28, -2.85, 3.05, -1.66, 1.75, 3.52, 2.39, 0.75, -2.59, -4.36, -2.54),
                      Z3_2 = c(-2.32, -1.67, 0.93, 0.26, 1.94, 4.06, 1.95, -1.73, 0.89, -1.3, -0.22, 1.62, 0.27, 0.5, 2.5, -1.07, -2.18, -2.64, 3.94, 2.44),
                      Z3_3 = c(0.6, 3.71, 1.93, -0.11, 1.06, 0.36, 0.26, -1.71, -2.49, -1.49, 0.47, 1.04, 1.84, -1.44, -3.5, 1.15, -1.12, -1.54, 0.59, 0.43))
  
  match_amino_acid <- function(AA) {
    d <- descs[descs$AA1 == AA, c(2, 3, 4)]
    return(d)
  }
  des <- sapply(alpha,match_amino_acid)
  ColDesNames <- paste0("AA_Position_", rep(colnames(des),  each = 3), "_",  rep(rownames(des), length(alpha)))
  des <- unlist(des)
  names(des) <- ColDesNames
  return(des)
}

alpha_z3("Hb-alpha.fasta")



results_alpha_one <- alpha_z3("Hb-alpha.fasta")
results_alpha_one <- t(results_alpha_one)
results_alpha_one <- data.frame(results_alpha_one)
results_alpha_one <- results_alpha_one[, order(names(results_alpha_one))]

df_alpha <- read.csv("data.csv")
descriptors <- df_alpha[, 2:ncol(df_alpha)]
descriptors <- descriptors[, order(names(descriptors))]
names(descriptors) <- names(results_alpha_one)
alpha_final <- data.frame(Affinity = df_alpha$Affinity, descriptors)
write.csv(alpha_final, file = "alpha_data_final.csv", row.names = FALSE)

