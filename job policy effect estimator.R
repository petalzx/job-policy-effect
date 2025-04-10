library(glmnet)
library(MASS)
library(dplyr)

jtpa.cq <- read.csv("mydata.csv", header = TRUE)[, -(1:2)]

# Pre-earning categories
income_bracket <- function(x) {
  case_when(x <= 1000 ~ 1, x > 1000 & x <= 4000 ~ 2, x > 4000 ~ 3)
}
jtpa.cq$preearn <- sapply(jtpa.cq$bfyrearn, income_bracket)

# Education categories
edu_bracket <- function(x) {
  case_when(x <= 8 ~ 1,
            x > 8 & x <= 10 ~ 2,
            x > 10 & x <= 11 ~ 3,
            x > 11 & x <= 12 ~ 4,
            x > 12 & x <= 14 ~ 5,
            x > 14 ~ 6)
}
jtpa.cq$preedu <- sapply(jtpa.cq$bfeduca, edu_bracket)

# Create 18-category W variable
jtpa.cq$W <- as.numeric(interaction(jtpa.cq$preearn, jtpa.cq$preedu))

# Set up variables
set.seed(210010)
Y <- jtpa.cq$earnings
D <- jtpa.cq$D
X <- model.matrix(~ (
  male + hsorged + black + hispanic + married + wkless13 +
    scale(age) + scale(bfeduca) + scale(bfyrearn)
)^2 - 1,
data = jtpa.cq)
W <- jtpa.cq$W

# Create complete dataset
analysis_data <- data.frame(Y, D, W) %>%
  bind_cols(as.data.frame(X)) %>%
  na.omit()

# Create folds
create_folds <- function(n, K) {
  indices <- sample(1:n)
  split(indices, cut(seq_along(indices), K, labels = FALSE))
}

# Gamma Estimation Function
estimate_gamma <- function(data, treatment) {
  sub_data <- data %>% filter(D == treatment)
  if (nrow(sub_data) == 0)
    return(NULL)
  
  X <- as.matrix(sub_data %>% select(-Y, -D, -W))
  Y <- sub_data$Y
  
  cv.glmnet(X, Y, alpha = 1, nfolds = 10)
}

# Compute G and P Matrices
compute_GP <- function(data, gamma1_pred, gamma0_pred) {
  X <- as.matrix(data %>% select(-Y, -D, -W))
  n <- nrow(data)
  
  list(
    G1 = (t(X) %*% (data$D * X)) / n,
    G0 = (t(X) %*% ((1 - data$D) * X)) / n,
    P = (2 * t(X) %*% (gamma1_pred - gamma0_pred)) / n
  )
}

compute_omega_cv <- function(train_data,
                             test_data,
                             lambda_grid,
                             treatment,
                             method = "catergory_weighted") {
  # Gamma estimation
  gamma1_model <- estimate_gamma(train_data, 1)
  gamma0_model <- estimate_gamma(train_data, 0)
  
  # Predictions
  X_train <- as.matrix(train_data %>% select(-Y, -D, -W))
  gamma1_train <- predict(gamma1_model, X_train, s = "lambda.min")
  gamma0_train <- predict(gamma0_model, X_train, s = "lambda.min")
  
  # GP matrices
  GP <- compute_GP(train_data, gamma1_train, gamma0_train)
  G <- if (treatment == 1) GP$G1 else GP$G0
  P <- GP$P
  
  # Cross-validation errors
  cv_errors <- sapply(lambda_grid, function(lambda) {
    a <- ginv(G %*% G + lambda * G) %*% G %*% P
    omega_test <- as.matrix(test_data %>% select(-Y, -D, -W)) %*% a
    
    gamma_diff <- predict(gamma1_model, as.matrix(test_data %>% select(-Y, -D, -W)), s = "lambda.min") -
      predict(gamma0_model, as.matrix(test_data %>% select(-Y, -D, -W)), s = "lambda.min")
    
    term <- if (treatment == 1) {
      test_data$D * omega_test * test_data$Y - 2 * gamma_diff *
        predict(gamma1_model, as.matrix(test_data %>% select(-Y, -D, -W)), s = "lambda.min")
    } else {
      (1 - test_data$D) * omega_test * test_data$Y - 2 * gamma_diff *
        predict(gamma0_model, as.matrix(test_data %>% select(-Y, -D, -W)), s = "lambda.min")
    }
    
    if (method == "catergory_weighted") {
      sum(sapply(1:18, function(s) {
        sum(abs(term * (test_data$W == s)))
      }))
    } else if (method == "aggregate") {
      sum(abs(term))
    }
  })
  
  cv_errors
}

# Estimate treatment effects and policy-relevant predictions 
estimate_treatment_effect <- function(data = analysis_data, K = 5, method = "category_weighted") {
  n <- nrow(data)
  lambda_grid <- seq(0.1, 3, by = 0.1)
  folds <- create_folds(n, K)
  results <- list()
  
  for (k in 1:K) {
    cat("Processing fold", k, "/", K, "\n")
    I_k <- folds[[k]]
    I_kc <- setdiff(1:n, I_k)
    
    data_k <- data[I_k, ]
    data_kc <- data[I_kc, ]
    
    # Estimate gamma
    gamma1_model <- estimate_gamma(data_kc, 1)
    gamma0_model <- estimate_gamma(data_kc, 0)
    
    # Predict gamma for current fold
    X_k <- as.matrix(data_k %>% select(-Y, -D, -W))
    gamma1_k <- predict(gamma1_model, newx = X_k, s = "lambda.min")
    gamma0_k <- predict(gamma0_model, newx = X_k, s = "lambda.min")
    
    # Omega estimation with inner CV
    inner_folds <- create_folds(nrow(data_kc), 10)
    total_cv1 <- total_cv0 <- numeric(length(lambda_grid))
    
    for (j in 1:10) {
      inner_I_j <- inner_folds[[j]]
      inner_train <- data_kc[-inner_I_j, ]
      inner_test <- data_kc[inner_I_j, ]
      
      # Omega1 CV
      cv1 <- compute_omega_cv(inner_train, inner_test, lambda_grid, 1, method)
      total_cv1 <- total_cv1 + cv1
      
      # Omega0 CV
      cv0 <- compute_omega_cv(inner_train, inner_test, lambda_grid, 0, method)
      total_cv0 <- total_cv0 + cv0
    }
    
    # Optimal lambdas
    lambda1 <- lambda_grid[which.min(total_cv1)]
    cat("lambda1: ", lambda1, "\n")
    lambda0 <- lambda_grid[which.min(total_cv0)]
    cat("lambda0: ", lambda0, "\n")
    
    # Final omega coefficients
    GP <- compute_GP(
      data_kc,
      predict(
        gamma1_model,
        newx = as.matrix(data_kc %>% select(-Y, -D, -W)),
        s = "lambda.min"
      ),
      predict(
        gamma0_model,
        newx = as.matrix(data_kc %>% select(-Y, -D, -W)),
        s = "lambda.min"
      )
    )
    
    # Compute omega1
    mat1 <- GP$G1 %*% GP$G1 + lambda1 * GP$G1
    a1 <- ginv(mat1) %*% GP$G1 %*% GP$P
    omega1_k <- X_k %*% a1
    
    # Compute omega0
    mat0 <- GP$G0 %*% GP$G0 + lambda0 * GP$G0
    a0 <- ginv(mat0) %*% GP$G0 %*% GP$P
    omega0_k <- X_k %*% a0
    
    results[[k]] <- data_k %>%
      mutate(
        gamma1 = as.numeric(gamma1_k),
        gamma0 = as.numeric(gamma0_k),
        omega1 = as.numeric(omega1_k),
        omega0 = as.numeric(omega0_k)
      )
  }
  
  # Combine results and calculate xi
  final_data <- bind_rows(results) %>%
    mutate(xi = (gamma1 - gamma0)^2 +
             (D * omega1 * (Y - gamma1) - (1 - D) * omega0 * (Y - gamma0)))
  
  return(final_data)
}

# Example Usage
results <- estimate_treatment_effect(method = "aggregate")

# Reformat results for step 2
reformat_results <- function(df = results) {
  formatted_data <- df %>% select(-omega1, -omega0) %>% 
    select(xi, Y, D, male:`scale(bfeduca):scale(bfyrearn)`, W, gamma1, gamma0)
  return(formatted_data)
}

# Compute proportion (aˆj) of individuals who should receive job training in each group
compute_optimal_props <- function(data) {
  a_hat <- numeric(18)
  
  for (j in 1:18) {
    group_data <- data %>% filter(W == j)
    
    objective <- function(a) {
      with(group_data, sum(xi * (as.numeric(gamma1 - gamma0 >= 0) - a) ^ 2))
    }
    
    opt_result <- optim(
      par = 0.5,
      fn = objective,
      method = "L-BFGS-B",
      lower = 0,
      upper = 1
    )
    
    a_hat[j] = opt_result$par
  }
  return(a_hat)
}

# Create heat map of optimal proportion of individuals
# in each education and income group who should receive job training

generate_heatmap <- function(proportions) {
  prop_matrix <- matrix(proportions, nrow = 3, ncol = 6)
  
  library(reshape)
  library(ggplot2)
  
  our_rule_plot <- melt(prop_matrix)
  colnames(our_rule_plot) <- c("pre-program income", "pre-program education", "decision rule")
  
  # Create the base plot and add all layers in one go
  heatmap <- ggplot(our_rule_plot, aes(`pre-program education`, `pre-program income`)) +
    geom_tile(aes(fill = `decision rule`)) +
    scale_fill_gradient(low = "white", high = "#d12168", limits = c(0, 1)) +
    labs(x = "pre-program education bracket", 
         y = "pre-program income bracket", 
         title = "Our proposed approach") +
    annotate("text", x = 6, y = 1, label = sprintf("%.2f", prop_matrix[1, 6])) +
    annotate("text", x = 3, y = 2, label = sprintf("%.2f", prop_matrix[2, 3])) +
    annotate("text", x = 3, y = 3, label = sprintf("%.2f", prop_matrix[3, 3]))
  
  return(heatmap)
  
  # Adding some example annotations (you can modify these based on actual values)
}

formatted_data <- reformat_results()
proportions <- compute_optimal_props(formatted_data)
heatmap_plot <- generate_heatmap(proportions)

print(heatmap_plot)
