library(seqinr)
fasta <- read.fasta("Hb-beta.fasta")
df <- rapply(fasta, c)
df <- data.frame(df)
df <- t(df)
df_capital <- toupper(df)

mutated_position <- c(1, 2, 5, 7, 8, 9, 11, 15, 17, 20, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 40, 41, 42, 45, 47, 48, 51, 52, 54, 57, 59, 60,
                      61, 62, 63, 64, 65, 66, 67, 68, 70, 72, 73, 74, 75, 76, 79, 81, 82, 84, 85, 86, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 
                      105, 106, 107, 108, 109, 111, 113, 116, 117, 118, 119, 120, 121, 122, 123, 124, 126, 128, 129, 130, 131, 132, 135, 136, 138, 139, 140,
                      142, 143, 144, 145, 146)

beta <- df_capital[mutated_position]
names(beta) <- mutated_position
descs <- data.frame(AA1 = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"),
                    Z3_1 = c(0.24, 0.84, 3.98, 3.11, -4.22, 2.05, 2.47, -3.89, 2.29, -4.28, -2.85, 3.05, -1.66, 1.75, 3.52, 2.39, 0.75, -2.59, -4.36, -2.54),
                    Z3_2 = c(-2.32, -1.67, 0.93, 0.26, 1.94, 4.06, 1.95, -1.73, 0.89, -1.3, -0.22, 1.62, 0.27, 0.5, 2.5, -1.07, -2.18, -2.64, 3.94, 2.44),
                    Z3_3 = c(0.6, 3.71, 1.93, -0.11, 1.06, 0.36, 0.26, -1.71, -2.49, -1.49, 0.47, 1.04, 1.84, -1.44, -3.5, 1.15, -1.12, -1.54, 0.59, 0.43))

match_amino_acid <- function(AA) {
  d <- descs[descs$AA1 == AA, c(2, 3, 4)]
  return(d)
}


des <- sapply(beta,match_amino_acid)
ColDesNames <- paste0("AA_Position", rep(colnames(des),  each = 3), "_",  rep(rownames(des), length(alpha)))
des <- unlist(des)
names(des) <- ColDesNames
des


beta_z3 <- function(x) {
  fasta <- read.fasta("Hb-beta.fasta")
  df <- rapply(fasta, c)
  df <- data.frame(df)
  df <- t(df)
  df_capital <- toupper(df)
  
  mutated_position <- c(1, 2, 5, 7, 8, 9, 11, 15, 17, 20, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 40, 41, 42, 45, 47, 48, 51, 52, 54, 57, 59, 60,
                        61, 62, 63, 64, 65, 66, 67, 68, 70, 72, 73, 74, 75, 76, 79, 81, 82, 84, 85, 86, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 
                        105, 106, 107, 108, 109, 111, 113, 116, 117, 118, 119, 120, 121, 122, 123, 124, 126, 128, 129, 130, 131, 132, 135, 136, 138, 139, 140,
                        142, 143, 144, 145, 146)
  
  beta <- df_capital[mutated_position]
  names(beta) <- mutated_position
  descs <- data.frame(AA1 = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"),
                      Z3_1 = c(0.24, 0.84, 3.98, 3.11, -4.22, 2.05, 2.47, -3.89, 2.29, -4.28, -2.85, 3.05, -1.66, 1.75, 3.52, 2.39, 0.75, -2.59, -4.36, -2.54),
                      Z3_2 = c(-2.32, -1.67, 0.93, 0.26, 1.94, 4.06, 1.95, -1.73, 0.89, -1.3, -0.22, 1.62, 0.27, 0.5, 2.5, -1.07, -2.18, -2.64, 3.94, 2.44),
                      Z3_3 = c(0.6, 3.71, 1.93, -0.11, 1.06, 0.36, 0.26, -1.71, -2.49, -1.49, 0.47, 1.04, 1.84, -1.44, -3.5, 1.15, -1.12, -1.54, 0.59, 0.43))
  
  match_amino_acid <- function(AA) {
    d <- descs[descs$AA1 == AA, c(2, 3, 4)]
    return(d)
  }
  
  
  des <- sapply(beta,match_amino_acid)
  ColDesNames <- paste0("AA_Position_", rep(colnames(des),  each = 3), "_",  rep(rownames(des), length(alpha)))
  des <- unlist(des)
  names(des) <- ColDesNames
  return(des)
}



results_beta_one <- beta_z3("Hb-beta.fasta")
results_beta_one <- t(results_beta_one)
results_beta_one <- data.frame(results_beta_one)
results_beta_one <- results_beta_one[, order(names(results_beta_one))]

df_beta <- read.csv("beta_data.csv")
descriptors <- df_beta[, 2:ncol(df_beta)]
descriptors <- descriptors[, order(names(descriptors))]
names(descriptors) <- names(results_beta_one)
beta_final <- data.frame(Affinity = df_beta$Affinity, descriptors)
write.csv(beta_final, file = "beta_data_final.csv", row.names = FALSE)
