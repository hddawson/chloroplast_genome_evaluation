library(tidyverse)
library(pROC)

df <- read_csv("../agent/data/pfam_experiment_results.csv")
head(df)

conds <- c("baseline","pfam","pfam_seq","bicor","pfam_seq_bicor")

rocs <- map(conds, ~roc(df$label, df[[paste0("pred_", .x)]], quiet=TRUE))
names(rocs) <- conds

plot(rocs[[1]])
walk(rocs[-1], ~plot(.x, add=TRUE))
legend("bottomright", legend=conds, lwd=2)

df_long <- df |>
  pivot_longer(starts_with("pred_"),
               names_to="condition",
               values_to="pred") |>
  mutate(condition=str_remove(condition,"pred_"),
         label=factor(label, levels=c(0,1), labels=c("Null","True")))

ggplot(df_long, aes(label, pred, fill=label)) +
  geom_boxplot(outlier.shape=NA) +
  facet_wrap(~condition, nrow=1) +
  coord_cartesian(ylim=c(0,1))


df_long |>
  group_by(condition, label) |>
  summarise(mean_pred=mean(pred, na.rm=TRUE)) |>
  pivot_wider(names_from=label, values_from=mean_pred) |>
  mutate(delta=True-Null)

ggplot(df, aes(as.factor(label), bicor)) + geom_boxplot()

cor(df$bicor, df$pred_bicor, use="complete.obs")
cor(df$bicor, df$pred_pfam_seq_bicor, use="complete.obs")

plot(df$bicor, df$pred_bicor)
plot(df$bicor, df$pred_pfam_seq_bicor)

ggplot(df, aes(bicor, pred_bicor)) + geom_point() + geom_smooth()

df |>
  filter(label == 0) |>
  summarise(
    bicor = mean(pred_bicor, na.rm=TRUE),
    pfam_seq_bicor = mean(pred_pfam_seq_bicor, na.rm=TRUE)
  )

df |>
  filter(bicor > 0.7, label == 0) |>
  summarise(
    bicor = mean(pred_bicor, na.rm=TRUE),
    pfam_seq_bicor = mean(pred_pfam_seq_bicor, na.rm=TRUE)
  )


conds <- c("baseline","pfam","pfam_seq","bicor","pfam_seq_bicor")

auc_df <- tibble(
  condition = conds,
  auc = map_dbl(conds, ~auc(roc(df$label, df[[paste0("pred_", .x)]], quiet=TRUE)))
)


ggplot(auc_df, aes(reorder(condition, auc), auc)) +
  geom_col() +
  coord_flip() +
  ylim(0,1)

ggplot(df, aes(bicor, pred_bicor, color=factor(label))) +
  geom_point(alpha=0.7) +
  geom_smooth(method="lm", se=FALSE) +
  ylim(0,1)

ggplot(df, aes(bicor, pred_pfam_seq_bicor, color=factor(label))) +
  geom_point(alpha=0.7) +
  geom_smooth(method="lm", se=FALSE) +
  ylim(0,1)


df_long <- df |>
  select(bicor, label, pred_bicor, pred_pfam_seq_bicor) |>
  pivot_longer(starts_with("pred_"),
               names_to="condition",
               values_to="pred")

ggplot(df_long, aes(bicor, pred, color=factor(label))) +
  geom_point(alpha=0.6) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~condition) +
  ylim(0,1)

cor(df$bicor, df$pred_bicor, use="complete.obs")
cor(df$bicor, df$pred_pfam_seq_bicor, use="complete.obs")

# ---- changing prompt 

library(tidyverse)

dfs <- list(
  minimal = read_csv("../agent/data/pfam_experiment_results.csv"),
  protein   = read_csv("../agent/data/1_pfam_experiment_results.csv"),
  plant_prune = read_csv("../agent/data/betterPrompt_pfam_experiment_results.csv")
)

df_long <- bind_rows(
  lapply(names(dfs), function(n)
    dfs[[n]] |> mutate(prompt=n))
)

ggplot(df_long,
       aes(bicor, pred_pfam_seq_bicor, color=factor(label))) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~prompt) +
  ylim(0,1)

df_long |>
  group_by(prompt) |>
  summarise(cor_bicor = cor(bicor, pred_pfam_seq_bicor, use="complete.obs"))
df_wide <- df_long |>
  select(gene_a, gene_b, prompt, pred_pfam_seq_bicor) |>
  pivot_wider(names_from=prompt, values_from=pred_pfam_seq_bicor)

df_wide |>
  mutate(sd_pred = apply(select(., -gene_a, -gene_b), 1, sd, na.rm=TRUE)) |>
  summarise(
    mean_sd = mean(sd_pred),
    max_sd = max(sd_pred)
  )
