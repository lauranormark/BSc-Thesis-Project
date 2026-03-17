getwd()

library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

#laddar ner datan
data_immunoabundance <- read.delim(
  "Data/merged_immunoabundance_clinical.txt",
  header = TRUE,
  check.names = FALSE
)

#violin plots females mot males
violin_plot <- data_immunoabundance %>%
  filter(!is.na(SEX), SEX %in% c("Female", "Male")) %>%  # ta bort NA och oväntade värden
  ggplot(aes(x = SEX, y = ABUNDANCE, fill = SEX)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_compare_means(
    method        = "t.test",
    label         = "p.format",   # visar faktiskt p-värde
    label.x       = 1.5,          # centrerar p-värdet över violinerna
    hide.ns       = FALSE
  ) +
  facet_wrap(~CELL_TYPE, scales = "free_y") +
  theme_classic() +
  labs(
    title = "Immune cell deconvolution in GBM by sex",
    x     = "Sex",
    y     = "Relative abundance"
  )

violin_plot

colnames(data_immunoabundance)

violin_plot2 <- data_immunoabundance %>%
  filter(!is.na(OS_STATUS), OS_STATUS %in% c("0:LIVING", "1:DECEASED")) %>%  # ta bort NA och oväntade värden
  ggplot(aes(x = OS_STATUS, y = ABUNDANCE, fill = OS_STATUS)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_compare_means(
    method        = "t.test",
    label         = "p.format",   # visar faktiskt p-värde
    label.x       = 1.5,          # centrerar p-värdet över violinerna
    hide.ns       = FALSE
  ) +
  facet_wrap(~CELL_TYPE, scales = "free_y") +
  theme_classic() +
  labs(
    title = "Immune cell deconvolution in GBM by overall survival status",
    x     = "Overalll survival status",
    y     = "Relative abundance"
  )

violin_plot2

# Sparar bilden på violin plot i "Plots"
ggsave(
  filename = "Plots/violin_plot_sex.pdf",
  plot     = violin_plot,
  width    = 10, height = 7
)

ggsave(
  filename = "Plots/violin_plot_OS.pdf",
  plot     = violin_plot2,
  width    = 10, height = 7
)

# t_tester 
t_results <- immune_long %>%
  filter(!is.na(SEX), SEX %in% c("Female", "Male")) %>%  # <-- anpassa till dina exakta värden
  group_by(CELL_TYPE) %>%
  summarise(
    n_male        = sum(SEX == "Male"),
    n_female      = sum(SEX == "Female"),
    mean_male     = mean(ABUNDANCE[SEX == "Male"],   na.rm = TRUE),
    mean_female   = mean(ABUNDANCE[SEX == "Female"], na.rm = TRUE),
    p_value       = t.test(ABUNDANCE ~ SEX)$p.value
  )

t_results

write.table(t_results, #male gene expression 
            file = "Data/t_result_sex_immunoabundance.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
