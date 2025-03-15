##Social ALAN
library(Ouyang)


setwd("~/UNR/Social/ALAN/Data/Behavior")

# Specify the directory containing the data files
data_dir <- "~/UNR/Social/ALAN/Data/Behavior/Raw Data"
#data_dir <- "~/UNR/Social/ALAN/Data/Behavior/Raw Data2"
#data_dir <- "~/UNR/Social/ALAN/Data/Behavior/Raw Data3"
#data_dir <- "~/UNR/Social/ALAN/Data/Behavior/Raw Data4"

# Get a list of all the data files in the directory
file_paths <- list.files(data_dir, full.names = TRUE)

# Reformat data files
data <- reformat_data_files(file_paths)

# Attach metadata and calculate phase
data <- attach_meta(data, "meta.csv")
#data <- attach_meta(data, "meta2.csv")
#data <- attach_meta(data, "meta3.csv")
#data <- attach_meta(data, "meta4.csv")

# Print the final data frame with metadata and phase
print(data)


#Add all datas together

#dat <- data
dat <- rbind(dat, data)





df_subset <- data[data$Cage != 13, ]
df_subset <- df_subset[df_subset$Cage != 3, ]
df_subset <- df_subset[df_subset$Cage != 4, ]
df_subset <- df_subset[df_subset$Cage != 5, ]
df_subset <- df_subset[df_subset$Cage != 9, ]

library(ggplot2)
act_plot <- ggplot(df_subset, aes(x = Hour, y = Date, fill = HopsPerMinute)) +
  geom_tile() +
  facet_wrap(~ Cage, scales = "free_y") +
  labs(x = "Hour", y = "Date", fill = "Hops") +
  theme_minimal()

act_plot


#reformat_for_Rethomics(data)



# Load the required libraries
library(ggplot2)
library(dplyr)



#data <- dat

# Convert 'timestamp' to POSIXct format for better handling
dat$timestamp <- as.POSIXct(dat$timestamp)

# Create a new column 'DayNight' indicating whether it's day or night
dat <- dat %>%
  mutate(DayNight = ifelse(Phase == "night", "Night", "Day"))



proportion_data2 <- dat %>% 
  group_by(Treatment, Hour) %>% 
  summarise(Proportion = mean(HopsPerMinute),
            SEM = sd(HopsPerMinute) / sqrt(n()),
            SD = sd(HopsPerMinute), 
            n = n())





# Subset the dataset to only include rows with Treatment "A" or "LD"
iso <- proportion_data2 %>%
  filter(Treatment %in% c("A", "LD"))

actfig_iso <- ggplot(iso, aes(x = Hour, y = Proportion, color = Treatment, fill = Treatment, group = Treatment)) +
  geom_rect(aes(xmin=-1, xmax=9, ymin=0, ymax=Inf), colour="transparent", fill="grey85", alpha=0.2) +
  geom_rect(aes(xmin=21, xmax=24, ymin=0, ymax=Inf), colour="transparent", fill="grey85", alpha=0.2) +
  geom_errorbar(aes(ymin = Proportion - SEM, ymax = Proportion + SEM), 
                position = position_dodge(width = 0.7), width = 0.2, size = 0.5) +
  geom_line(aes(linetype = Treatment, color = Treatment), position = position_dodge(width = 0.7), size = 0.8) +
  geom_point(position = position_dodge(width = 0.7), size = 3, aes(shape = Treatment, fill = Treatment)) +
  scale_shape_manual(values = c("A" = 21, "LD" = 21, "SA" = 15, "SLD" = 15)) +
  #scale_fill_manual(values = c("A" = "yellow", "LD" = "black", "SA" = "orange", "SLD" = "black")) +
  scale_fill_manual(values = c("A" = "#FFD700", "LD" = "black", "SA" = "orange", "SLD" = "black")) +
  scale_colour_manual(values = c("A" = "#FFD700", "LD" = "black", "SA" = "orange", "SLD" = "black")) +
  scale_y_continuous(limits = c(0, 6.5)) +
  scale_linetype_manual(values = c("A" = "solid", "LD" = "dashed")) +  # Specify line type for each Treatment
  labs(x = "Hour of the day", y = "Average Activity per Hour", color = "Treatment Group") +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(family="Times", size=16, colour = "black"),
        axis.text.y = element_text(family="Times", size=16, colour = "black"),
        axis.title.x = element_text(family="Times", size=24, colour = "black"),
        axis.title.y = element_text(family="Times", size=24, angle=90))


actfig_iso #saved as 630 and 455
# Save the figure as a high-quality image
ggsave("actfig_iso.png", plot = actfig_iso, width = 6.3, height = 4.55, dpi = 300, units = "in")


# Subset the dataset 
soc <- proportion_data2 %>%
  filter(Treatment %in% c("SA", "SLD"))

actfig_soc <- ggplot(soc, aes(x = Hour, y = Proportion, color = Treatment, fill = Treatment, group = Treatment)) +
  geom_rect(aes(xmin=-1, xmax=9, ymin=0, ymax=Inf), colour="transparent", fill="grey85", alpha=0.2) +
  geom_rect(aes(xmin=21, xmax=24, ymin=0, ymax=Inf), colour="transparent", fill="grey85", alpha=0.2) +
  geom_errorbar(aes(ymin = Proportion - SEM, ymax = Proportion + SEM), 
                position = position_dodge(width = 0.7), width = 0.2, size = 0.5) +
  geom_line(aes(linetype = Treatment), position = position_dodge(width = 0.7), size = 0.8) +
  geom_point(position = position_dodge(width = 0.7), size = 3, aes(shape = Treatment)) +
  scale_shape_manual(values = c("A" = 16, "LD" = 16, "SA" = 15, "SLD" = 15)) +
  #scale_fill_manual(values = c("A" = "yellow", "LD" = "black", "SA" = "orange", "SLD" = "black")) +
  scale_colour_manual(values = c("A" = "yellow", "LD" = "black", "SA" = "orange", "SLD" = "black")) +
  scale_y_continuous(limits = c(0, 6.5)) +
  scale_linetype_manual(values = c("SA" = "solid", "SLD" = "dashed")) +  # Specify line type for each Treatment
  labs(x = "Hour of the day", y = "Average Activity per Hour", color = "Treatment Group") +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(family="Times", size=16, colour = "black"),
        axis.text.y = element_text(family="Times", size=16, colour = "black"),
        axis.title.x = element_text(family="Times", size=24, colour = "black"),
        axis.title.y = element_text(family="Times", size=24, angle=90))


actfig_soc #saved as 630 and 455

ggsave("actfig_soc.png", plot = actfig_soc, width = 6.3, height = 4.55, dpi = 300, units = "in")

library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)

# Rename the Treatment variable
proportion_data2 <- proportion_data2 %>%
  mutate(Treatment = recode(Treatment,
                            "A" = "Isolated ALAN",
                            "LD" = "Isolated Control",
                            "SA" = "Social ALAN",
                            "SLD" = "Social Control"))

# Check the updated dataset
head(proportion_data2)

# Dummy plot to generate the legend
dummy_plot <- ggplot(data = proportion_data2, aes(x = Hour, y = Proportion, color = Treatment, shape = Treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Isolated ALAN" = "#FFD700", "Isolated Control" = "black", "Social ALAN" = "orange", "Social Control" = "black")) +
  scale_shape_manual(values = c("Isolated ALAN" = 16, "Isolated Control" = 16, "Social ALAN" = 15, "Social Control" = 15)) +
  theme_minimal() +
  theme(legend.position = "bottom")  # Position legend for extraction
dummy_plot
# Extract the legend
legend <- get_legend(dummy_plot)

# Create a figure with just the legend
legend_plot <- ggarrange(legend, ncol = 1, nrow = 1)
legend_plot
# Save the legend as a separate image
ggsave("figure_legend.png", legend_plot, width = 6.5, height = 1.5, dpi = 300)

actfig2 <- ggplot(proportion_data2, aes(x = Hour, y = Proportion, color = Treatment, fill = Treatment, group = Treatment)) +
  geom_rect(aes(xmin=-1, xmax=9, ymin=0, ymax=Inf), colour="transparent", fill="grey85", alpha=0.2) +
  geom_rect(aes(xmin=21, xmax=24, ymin=0, ymax=Inf), colour="transparent", fill="grey85", alpha=0.2) +
  geom_errorbar(aes(ymin = Proportion - SEM, ymax = Proportion + SEM), 
                position = position_dodge(width = 0.7), width = 0.2, size = 0.5) +
  geom_point(position = position_dodge(width = 0.7), size = 3, aes(shape = Treatment)) +
  scale_shape_manual(values = c("A" = 16, "LD" = 16, "SA" = 15, "SLD" = 15)) +
  #scale_fill_manual(values = c("A" = "yellow", "LD" = "black", "SA" = "orange", "SLD" = "black")) +
  scale_colour_manual(values = c("A" = "yellow", "LD" = "black", "SA" = "orange", "SLD" = "black")) +
  labs(x = "Hour of the day", y = "Average activty", color = "Treatment Group") +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(family="Times", size=16, colour = "black"),
        axis.text.y = element_text(family="Times", size=16, colour = "black"),
        axis.title.x = element_text(family="Times", size=24, colour = "black"),
        axis.title.y = element_text(family="Times", size=24, angle=90))


actfig2


#For presentation 
library(dplyr)

filtered_data <- proportion_data2 %>%
  filter(Treatment %in% c("A", "LD", "SLD")) # Filters to include only treatments "A" and "SA"



ggplot(filtered_data, aes(x = Hour, y = Proportion, color = Treatment, fill = Treatment, group = Treatment)) +
  geom_rect(aes(xmin=-1, xmax=9, ymin=0, ymax=Inf), colour="transparent", fill="grey85", alpha=0.2) +
  geom_rect(aes(xmin=21, xmax=24, ymin=0, ymax=Inf), colour="transparent", fill="grey85", alpha=0.2) +
  geom_errorbar(aes(ymin = Proportion - SEM, ymax = Proportion + SEM), 
                position = position_dodge(width = 0.7), width = 0.2, size = 0.5) +
  geom_point(position = position_dodge(width = 0.7), size = 3, aes(shape = Treatment)) +
  scale_shape_manual(values = c("A" = 16, "LD" = 16, "SA" = 15, "SLD" = 15)) +
  #scale_fill_manual(values = c("A" = "yellow", "LD" = "black", "SA" = "orange", "SLD" = "black")) +
  scale_colour_manual(values = c("A" = "yellow", "LD" = "black", "SA" = "orange", "SLD" = "black")) +
  labs(x = "Hour of the day", y = "Proportion active", color = "Treatment Group") +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(family="Times", size=16, colour = "black"),
        axis.text.y = element_text(family="Times", size=16, colour = "black"),
        axis.title.x = element_text(family="Times", size=24, colour = "black"),
        axis.title.y = element_text(family="Times", size=24, angle=90))



###Calculate noctural and activity

noc <-subset(dat, Phase == "night")

noc_sum <- noc %>% 
  group_by(Treatment, Cage) %>% 
  summarise(noc_activity = sum(HopsPerMinute),
            mean = mean(HopsPerMinute),
            SEM = sd(HopsPerMinute) / sqrt(sum(HopsPerMinute)),
            SD = sd(HopsPerMinute), 
            n = sum(HopsPerMinute))



# Separate data into LD:A and SLD:SA groups
LD_A <- filter(noc_sum, Treatment %in% c("LD", "A"))
SLD_SA <- filter(noc_sum, Treatment %in% c("SLD", "SA"))

# Perform t-test for LD:A group
t_test_LD_A <- t.test(noc_activity ~ Treatment, data = LD_A)

# Perform t-test for SLD:SA group
t_test_SLD_SA <- t.test(noc_activity ~ Treatment, data = SLD_SA)

# Print t-test results
print(t_test_LD_A)

# Welch Two Sample t-test
# 
# data:  noc_activity by Treatment
# t = 4.2134, df = 9.0639, p-value = 0.002226
# alternative hypothesis: true difference in means between group A and group LD is not equal to 0
# 95 percent confidence interval:
#   1517.074 5026.493
# sample estimates:
#   mean in group A mean in group LD
# 3372.2000         100.4167

print(t_test_SLD_SA)


# Welch Two Sample t-test
# 
# data:  noc_activity by Treatment
# t = 4.2134, df = 9.0639, p-value = 0.002226
# alternative hypothesis: true difference in means between group A and group LD is not equal to 0
# 95 percent confidence interval:
#   1517.074 5026.493
# sample estimates:
#   mean in group A mean in group LD 
# 3372.2000         100.4167


###########