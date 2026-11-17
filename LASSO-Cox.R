
getwd()
library(survival)
library(glmnet)
library(readr)



data <- read.csv("data.csv")

str(data)

time <- data$PFS_TIME
status <- data$PFS_status

gene_data <- data[, 8:ncol(data)]

x <- as.matrix(gene_data)
y <- Surv(time, status)

sum(is.na(x))
sum(is.na(time))
sum(is.na(status))

set.seed(123)

lasso_fit <- glmnet(x, y, family = "cox", alpha = 1)

plot(lasso_fit, xvar = "lambda", label = T)
abline(a = 0, b = 1, v = 2.29, col = "grey", lwd = 2, lty = 3)

library(ggplot2)

beta_matrix <- as.matrix(coef(lasso_fit))
lambda_vector <- lasso_fit$lambda

plot_data <- data.frame(
  lambda = rep(lambda_vector, each = nrow(beta_matrix)),
  coefficient = as.vector(beta_matrix),
  variable = rep(rownames(beta_matrix), times = ncol(beta_matrix))
)

p <- ggplot(plot_data, aes(x = log(lambda), y = coefficient, group = variable, color = variable)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = -2.198868527, lty = 3) +
  scale_x_reverse() +
  coord_cartesian(ylim = c(-0.16, 0.01)) +
  labs(x = "Log(Lambda)", y = "Coefficient") +
  theme_bw() +
  scale_color_discrete(name = "Variable") +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 16)
  )

p

ggsave("figure1.pdf", p, width = 8, height = 6.5)

cv_lasso <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)

print(cv_lasso)

plot(cv_lasso)

par(
  font.axis = 1,
  font.lab = 1,
  cex.axis = 1.5,
  cex.lab = 1.5,
  cex.main = 2.0,
  cex.sub = 1.2
)

plot(
  cv_lasso,
  cex.axis = 1.5,
  cex.lab = 1.5,
  font.axis = 1,
  font.lab = 1,
  family = "Arial",
  mgp = c(2.7, 0.9, 0)
)

lambda_min <- cv_lasso$lambda.min
lambda_1se <- cv_lasso$lambda.1se

cat("Lambda min:", lambda_min, "\n")
cat("Lambda 1se:", lambda_1se, "\n")

coef_min <- coef(cv_lasso, s = "lambda.min")
coef_1se <- coef(cv_lasso, s = "lambda.1se")

selected_genes_min <- rownames(coef_min)[which(coef_min != 0)]
selected_genes_1se <- rownames(coef_1se)[which(coef_1se != 0)]

cat("Genes selected using lambda.min:", length(selected_genes_min), "\n")
print(selected_genes_min)

cat("Genes selected using lambda.1se:", length(selected_genes_1se), "\n")
print(selected_genes_1se)

coef_min <- coef(cv_lasso, s = "lambda.min")

nonzero_coef <- coef_min[coef_min != 0]

print(nonzero_coef)

genes <- rownames(coef_min)[which(coef_min != 0)]

betas <- coef_min[which(coef_min != 0)]

risk_formula <- paste0(
  "Risk Score = ",
  paste0(round(betas, 4), " × ", genes, collapse = " + ")
)

cat(risk_formula, "\n")

selected_data <- data.frame(
  time = time,
  status = status,
  gene_data[, selected_genes_min]
)

cox_model <- coxph(Surv(time, status) ~ ., data = selected_data)

summary(cox_model)

library(survminer)

ggforest(cox_model, data = selected_data)

ggforest(cox_model, data = selected_data, fontsize = 1.2) +
  theme(
    axis.text.x = element_text(size = 121),
    axis.text.y = element_text(size = 166),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )


