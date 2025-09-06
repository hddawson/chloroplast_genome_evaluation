library(data.table)
library(ggplot2)
data <- fread("data/alignment_input_data.csv")
table(data$Order)

Order_counts <- data[, .N, by = Order]
hist(Order_counts$N)

table(data$FixationLabel)
table(data$photosyntheticPathway)


order_dt <- melt(
  data,
  id.vars = "Order",
  measure.vars = "pheno_Topt_site_p50",
  value.name = "Topt_site_p50")
order_dt$Topt_site_p50 <-order_dt$Topt_site_p50 / 100 


par(mfrow=c(1,2))
hist(order_dt$Topt_site_p50, main="Species median T_opt_site",xlab="T_opt_site (C)")
hist(order_dt_bio8$wc2.1_2.5m_bio_8_p50, main="WC_bio8", xlab="Mean Temp wettest quarter (C)")
sort(table(order_dt$Topt_site_p50),decreasing = TRUE)


p <- ggplot(order_dt, aes(x = Topt_site_p50)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~Order, scales = "free_y") +
  theme_minimal(base_size = 10) +
  labs(
    title = "Distribution of Species Topt_site by Order",
    x = "Topt_site (°C)",
    y = "Count"
  ) +
  # add sample size labels
  geom_text(
    data = Order_counts,
    aes(x = Inf, y = Inf, label = paste0("n=", N)),
    inherit.aes = FALSE,
    hjust = 1.1, vjust = 1.5, size = 3
  )

p
