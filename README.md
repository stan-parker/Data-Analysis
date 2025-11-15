# Data-Analysis
A repository for the data and code used to produce results described in manuscript "Determining the temporal relationships between change in possession and player movement in Australian Football".

Data can be found in the below Dropbox folder. 

https://www.dropbox.com/scl/fo/g2ffvbo8osxa6wqty20lv/AJVvEmHjKWxmhyKE2B35YyU?rlkey=tjtvy83e2ub1p88lnwscmike2&st=vmucy4fd&dl=0



## Data Analysis

# Load libraries
library(tidyverse)
library(data.table)
library(here)
library(nlme)
library(arrow)
library(progress)
library(doParallel)
library(foreach)
library(broom.mixed)

# Set up parallel processing
num_cores <- detectCores() - 1  # Reserve 1 core for the OS
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Pull in all GPS data from offensive turnovers
GPS.Files <- list.files(here("Study 2 Offensive Chain Data/"), full.names = TRUE)
#GPS.Files <- list.files(here("OffensiveChainGPSOutcome/"), full.names = TRUE)
df <- rbindlist(lapply(GPS.Files, fread), idcol = TRUE, fill = TRUE)

# Remove times when theres > 18 people on the field
subset <- df %>%
  mutate(trimmed_name = str_remove(ID, "_[^_]+$"))

count <- subset %>%
  group_by(trimmed_name) %>%
  summarise(unique_names = n_distinct(`Player Name`)) %>% filter(unique_names <= 18)

subset2 <- count %>% left_join(subset, by = c("trimmed_name"))

# Get rid of faulty chains
subset2 <- subset2 %>% mutate(trimmed_name = str_remove(ID, "_[^_]+$"))
subset2 <- subset2 %>% filter(!str_count(trimmed_name, "__") == 1)

table(subset2_unique$unique_names)

# Sample data to 1Hz

df.samples <- subset2 %>% mutate(TimeRelativeToIntercept = floor(TimeRelativeToIntercept)) %>% 
  group_by(`Player Name`, ID, TimeRelativeToIntercept,start.tye, start.zone,chain.outcome) %>% 
  summarise(velocity = mean(Velocity),
            Distance.to.Intercept.Location = mean(Distance.to.Intercept.Location),
            chain.meters.gained.lag = mean(chain.meters.gained.lag),
            previous.chain.duration = mean(previous.chain.duration)) %>%
  mutate(
    previous.chain.duration.bucket = case_when( # Assign previous chain duration bucks
      previous.chain.duration < 5 ~ "< 5 sec",
      previous.chain.duration < 10 ~ "5–10 sec",
      previous.chain.duration < 15 ~ "10–15 sec",
      previous.chain.duration < 20 ~ "15–20 sec",
      previous.chain.duration < 25 ~ "20–25 sec",
      previous.chain.duration >= 25 ~ "≥ 25 sec",
      TRUE ~ "NA")) 

df.samples <- df.samples %>% filter(previous.chain.duration.bucket != "NA")


df.samples$previous.chain.duration.bucket <- factor(
  df.samples$previous.chain.duration.bucket,
  levels = c("< 5 sec", "5–10 sec", "10–15 sec", "15–20 sec","20–25 sec","≥ 25 sec"))


df.samples <- df.samples %>% mutate(
  start.type = ifelse(start.tye == "Intercept", "Intercept", "Stoppage")
)

# Filter pre-turnover and make seconds positive
pre <- df.samples %>% 
  filter(TimeRelativeToIntercept < 0) %>%
  mutate(TimeRelativeToIntercept = abs(TimeRelativeToIntercept))


# Compute player's average velocity in each start zone
pre_player_zone_duration_avg <- pre %>%
  group_by(`Player Name`, start.zone, start.type, previous.chain.duration.bucket) %>%
  summarise(zone_mean_velocity = mean(velocity, na.rm = TRUE), .groups = "drop")


# Convert to data.table and join
setDT(pre)
setDT(pre_player_zone_duration_avg)
pre <- merge(pre, pre_player_zone_duration_avg, by = c("Player Name", "start.zone", "start.type", "previous.chain.duration.bucket"), all.x = TRUE)
pre <- pre %>% filter(previous.chain.duration.bucket != "NA")
unique(pre$previous.chain.duration.bucket)
unique(pre$start.type)
unique(pre$start.zone)
unique(pre$chain.outcome)


pre <- pre %>% 
  mutate(
    chain.outcome = case_when( # Assign previous chain duration bucks
      chain.outcome == "Behind" ~ "Score",
      chain.outcome == "Goal" ~ "Score",
      chain.outcome == "Rushed Behind" ~ "Score",
      TRUE ~ "No Score"))

unique(pre$chain.outcome)


names(pre)[1] <- "name"
unique(pre$start.zone)

# Ensure factors and apply sum-to-zero contrasts
pre$name <- as.factor(pre$name)
pre$start.zone <- factor(pre$start.zone, levels = c("D50", "DM", "AM", "F50"))
pre$start.type <- factor(pre$start.type, levels = c("Stoppage", "Intercept"))
pre$previous.chain.duration.bucket <- factor(pre$previous.chain.duration.bucket, levels = c("< 5 sec", "5–10 sec", "10–15 sec", "15–20 sec","20–25 sec","≥ 25 sec"))
pre$chain.outcome <- factor(pre$chain.outcome, levels = c("No Score", "Score"))

contrasts(pre$name) <- contr.sum(length(levels(pre$name)))
contrasts(pre$start.zone) <- contr.sum(length(levels(pre$start.zone)))
contrasts(pre$start.type) <- contr.sum(length(levels(pre$start.type)))
contrasts(pre$previous.chain.duration.bucket) <- contr.sum(length(levels(pre$previous.chain.duration.bucket)))
contrasts(pre$chain.outcome) <- contr.sum(length(levels(pre$chain.outcome)))

# Remove unused column
#pre <- pre %>% select(-chain.outcome)

# Time windows for model fitting
time_windows <- c(5,10,15,20,25,30)

# Function to fit GLS model
fit_gls_model <- function(window_size) {
  pre_window <- pre %>%
    filter(TimeRelativeToIntercept <= window_size)
  
  pre_window$name <- as.factor(pre_window$name)
  pre_window$start.zone <- factor(pre_window$start.zone, levels = c("D50", "DM", "AM", "F50"))
  pre_window$start.type <- as.factor(pre_window$start.type)
  pre_window$previous.chain.duration.bucket <- factor(pre_window$previous.chain.duration.bucket, levels = c("< 5 sec", "5–10 sec", "10–15 sec", "15–20 sec","20–25 sec","≥ 25 sec"))
  
  
  contrasts(pre_window$name) <- contr.sum(length(levels(pre_window$name)))
  contrasts(pre_window$start.zone) <- contr.sum(length(levels(pre_window$start.zone)))
  contrasts(pre_window$start.type) <- contr.sum(length(levels(pre_window$start.type)))
  contrasts(pre_window$previous.chain.duration.bucket) <- contr.sum(length(levels(pre_window$previous.chain.duration.bucket)))
  
  model <- gls(
    velocity ~ name + start.zone + start.type + previous.chain.duration.bucket,
    correlation = corAR1(form = ~ TimeRelativeToIntercept | ID),
    data = pre_window,
    na.action = na.exclude
  ) 
  
  ## bin single seconds, plot performance - hopefully not perfectly linear...
  
  tibble(
    TimeWindow = window_size,
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = as.numeric(logLik(model)),
    Model = list(model)
  )
}

# Run models in parallel using foreach
model_results <- foreach(w = time_windows, .combine = bind_rows, .packages = c("nlme", "dplyr", "broom.mixed")) %dopar% {
  fit_gls_model(w)
}

# Stop parallel cluster
stopCluster(cl)

# Extract model coefficients
coefficients_df <- model_results %>%
  mutate(Tidy = map(Model, ~ broom.mixed::tidy(.))) %>%
  select(TimeWindow, Tidy) %>%
  unnest(Tidy)

# Plot coefficients over time
ggplot(coefficients_df %>% filter(term != "(Intercept)"), 
       aes(x = TimeWindow, y = estimate, color = term)) +
  geom_line() +
  geom_point(size = 2) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title = "All Coefficients Across Time Windows",
    x = "Seconds Before Turnover",
    y = "Estimate",
    color = "Term"
  )

library(scales)

ggplot(model_results, aes(x = TimeWindow, y = logLik)) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_continuous(breaks = c(5, 10, 15, 20, 25, 30)) +  # exact x-axis breaks
  scale_y_continuous(labels = scales::comma, breaks = pretty(model_results$logLik))+
  theme_minimal() +
  labs(
    #   title = "BIC Across Different Time Windows",
    #  subtitle = "Model Fit across various window lengths",
    x = "Seconds Before Turnover",
    y = "logLik"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

pre_window <- droplevels(pre_window)



## Investigate model structure ----

# Data Structure
time_window_summary <- pre %>% 
  mutate(window = cut(TimeRelativeToIntercept, breaks = c(0, 5, 10, 15, 20, 25, 30))) %>%
  group_by(window) %>%
  summarise(n = n())

ggplot(time_window_summary, aes(x = window, y = n)) +
  geom_bar(stat = "identity") +
  labs(title = "Sample Size by Time Window")

# Predictor Distribution
pre %>%
  mutate(window = cut(TimeRelativeToIntercept, breaks = c(0, 5, 10, 15, 20, 25, 30))) %>%
  group_by(window, start.zone) %>%
  summarise(n = n()) %>%
  group_by(window) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = window, y = prop, fill = start.zone)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Start Zone Distribution by Time Window")


### Analysis at 5 second window ---- 

pre_window <- pre %>%
  filter(TimeRelativeToIntercept <= 5)

pre_window <- pre_window %>% filter(start.zone != "NA")

unique(pre_window$start.zone)
names(pre_window)[1] <- "name"

pre_window$name <- as.factor(pre_window$name)
pre_window$start.zone <- factor(pre_window$start.zone, levels = c("D50", "DM","AM", "F50"))
pre_window$start.type <- as.factor(pre_window$start.type)
pre_window$previous.chain.duration.bucket <- factor(pre_window$previous.chain.duration.bucket, levels = c("< 5 sec", "5–10 sec", "10–15 sec", "15–20 sec", "20–25 sec", "≥ 25 sec"))
#pre_window$chain.outcome <- factor(pre_window$chain.outcome, levels = c("No Score", "Score"))


contrasts(pre_window$name) <- contr.sum(length(levels(pre_window$name)))
contrasts(pre_window$start.zone) <- contr.sum(length(levels(pre_window$start.zone)))
contrasts(pre_window$start.type) <- contr.sum(length(levels(pre_window$start.type)))
contrasts(pre_window$previous.chain.duration.bucket) <- contr.sum(length(levels(pre_window$previous.chain.duration.bucket)))
#contrasts(pre_window$chain.outcome) <- contr.sum(length(levels(pre_window$chain.outcome)))


# Fit GLS model
model <- gls(
  velocity ~ name + start.zone * start.type * previous.chain.duration.bucket,
  correlation = corAR1(form = ~ TimeRelativeToIntercept | ID),
  data = pre_window,
  na.action = na.exclude
)

# There arent iterations of > 60 second previous chain durations for all start zonees. 
summary(model)

table(pre_window$start.zone, pre_window$start.type, pre_window$previous.chain.duration.bucket)

table(pre_window$start.zone, pre_window$start.type, pre_window$previous.chain.duration.bucket)

model_matrix <- model.matrix(~ name + start.zone * start.type * previous.chain.duration.bucket, data = pre_window_clean)
alias(lm(velocity ~ name + start.zone * start.type * previous.chain.duration.bucket, data = pre_window_clean))


### Coung rare combinations 

count <- pre_window %>%
  count(start.zone, start.type, previous.chain.duration.bucket) %>%
  arrange(n)


# plot of previous chain duration
library(scales)

ggplot(pre, aes(x = previous.chain.duration.bucket)) +
  geom_bar() +
  #geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  scale_y_continuous(labels = comma) +
  theme_minimal() +
  labs(
    title = "Distribution of Previous Chain Duration",
    x = "Previous Chain Duration",
    y = "Number of unique incidents (n)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

# Model with chain outcome

model <- gls(
  velocity ~ name + (start.zone * start.type * previous.chain.duration.bucket) + chain.outcome,
  correlation = corAR1(form = ~ TimeRelativeToIntercept | ID),
  data = pre_window,
  na.action = na.exclude
)


summary(model)

#### Start Zone Marginal Means ----

library(emmeans)

# Estimated marginal means
emm_start_zone <- emmeans(model, ~ start.zone, mode = "df.error")

# Grand mean contrasts with significance labels
contrasts_grand_mean <- contrast(emm_start_zone, method = "eff") %>%
  as.data.frame() %>%
  mutate(
    start.zone = gsub(" effect", "", contrast),
    Signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    )
  )

# Estimated marginal means to data frame
emm_df <- as.data.frame(emm_start_zone)

# Join only by start.zone
emm_df_labeled <- emm_df %>%
  left_join(contrasts_grand_mean %>% select(start.zone, Signif), by = "start.zone")

# Optional: set order for plot
emm_df_labeled$start.zone <- factor(emm_df_labeled$start.zone, levels = c("D50", "DM", "AM", "F50"))

# Plot with significance labels
ggplot(emm_df_labeled, aes(x = start.zone, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(aes(label = Signif, y = upper.CL + 0.05), size = 6) +
  #theme_minimal(base_size = 14) +
  labs(
    x = "Start Zone",
    y = expression("Estimated Velocity (m·s"^-1*")"),
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, face = "italic")
  )



# Merge estimated means with effect sizes and p-values
table_df <- emm_df %>%
  left_join(
    contrasts_grand_mean %>%
      select(start.zone, estimate, SE, p.value, Signif),
    by = c("start.zone")
  ) %>%
  mutate(
    `Estimated Velocity (m/s)` = round(emmean, 3),
    `Effect vs Grand Mean (m/s)` = round(estimate, 3),
    `Standard Error` = round(SE, 3),
    `95% CI` = paste0("(", round(lower.CL, 3), ", ", round(upper.CL, 3), ")"),
    `p-value` = format_custom_pval(p.value)
  ) %>%
  select(
    `Start Zone` = start.zone,
    `Estimated Velocity (m/s)`,
    `Effect vs Grand Mean (m/s)`,
    `Standard Error`,
    `95% CI`,
    `p-value`,
    Signif
  )

# New version for table

# Get estimated marginal means
emm_start_zone <- emmeans(model, ~ start.zone, mode = "df.error")
emm_df <- as.data.frame(emm_start_zone)

# Get contrasts vs grand mean and calculate significance
contrasts_grand_mean <- contrast(emm_start_zone, method = "eff") %>%
  as.data.frame() %>%
  rename(
    Effect = estimate,
    StdError = SE
  ) %>%
  mutate(
    start.zone = gsub(" effect", "", contrast),
    Signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    )
  )

# Custom p-value formatter
format_custom_pval <- function(p, digits = 7) {
  ifelse(p < 1e-7, "<0.0000001", sprintf(paste0("%.", digits, "f"), p))
}

# Merge EMMs with contrasts for final table
table_df <- emm_df %>%
  left_join(
    contrasts_grand_mean %>%
      select(start.zone, Effect, StdError, p.value, Signif),
    by = "start.zone"
  ) %>%
  mutate(
    Estimate = round(emmean, 3),
    `Effect Size` = round(Effect, 3),
    `Standard Error` = round(StdError, 3),
    `p-value` = format_custom_pval(p.value)
  ) %>%
  select(
    `Start Zone` = start.zone,
    Estimate,
    `Effect Size`,
    `Standard Error`,
    `p-value`,
    Signif
  )

### Box Plot Previous Chain Data 

pre_window <- pre_window %>% mutate(velocity = velocity/3.6)

# join significance labels to factor order
sig_labels <- contrasts_grand_mean %>%
  select(start.zone, Signif)

# Plot with annotations
ggplot(pre_window, aes(x = start.zone, y = velocity)) +
  geom_boxplot(aes(fill = start.zone), alpha = 0.5, outlier.shape = NA) +
  # geom_jitter(width = 0.15, alpha = 0.2, size = 1) +
  geom_text(
    data = sig_labels,
    aes(x = start.zone, y = 5,
        label = Signif),
    inherit.aes = FALSE,
    size = 6
  ) +
  coord_cartesian(ylim = c(0, 6)) +   # 👈 adjust limits here
  scale_fill_viridis_d() +  
  labs(
    x = "Start Zone",
    y = expression("Velocity (m·s"^-1*")")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "none"
  )


## Start Type analysis ----


# 1. Get estimated marginal means for start type
emm_start_type <- emmeans(model, ~ start.type, mode = "df.error")
emm_df_type <- as.data.frame(emm_start_type)

# 2. Compare each start type to the grand mean
contrasts_start_type <- contrast(emm_start_type, method = "eff") %>%
  as.data.frame() %>%
  rename(
    Effect = estimate,
    StdError = SE
  ) %>%
  mutate(
    start.type = gsub(" effect", "", contrast),
    Signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    )
  )

# 3. Prepare labeled data for plotting
emm_df_labeled_type <- emm_df_type %>%
  left_join(contrasts_start_type %>% select(start.type, Signif), by = "start.type")

# 4. Plot estimated marginal means with significance
ggplot(emm_df_labeled_type, aes(x = start.type, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(aes(label = Signif, y = upper.CL + 0.05), size = 6) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Start Type",
    y = expression("Estimated Velocity (m·s"^-1*")"),
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, face = "italic")
  )

# 5. Custom p-value formatter
format_custom_pval <- function(p, digits = 7) {
  ifelse(p < 1e-7, "<0.0000001", sprintf(paste0("%.", digits, "f"), p))
}

# 6. Final table for academic reporting
table_df_type <- emm_df_type %>%
  left_join(
    contrasts_start_type %>% select(start.type, Effect, StdError, p.value, Signif),
    by = "start.type"
  ) %>%
  mutate(
    Estimate = round(emmean, 3),
    `Effect Size` = round(Effect, 3),
    `Standard Error` = round(StdError, 3),
    `p-value` = format_custom_pval(p.value)
  ) %>%
  select(
    `Start Type` = start.type,
    Estimate,
    `Effect Size`,
    `Standard Error`,
    `p-value`,
    Signif
  )


sig_labels_type <- contrasts_start_type %>%
  select(start.type, Signif)

ggplot(pre_window, aes(x = start.type, y = velocity)) +
  geom_boxplot(aes(fill = start.type), alpha = 0.7, colour = "black", outlier.shape = NA) +
  # optional: add raw data points for transparency
  # geom_jitter(width = 0.15, alpha = 0.2, size = 1) +
  geom_text(
    data = sig_labels_type,
    aes(x = start.type, y = 6.5,
        label = Signif),
    inherit.aes = FALSE,
    size = 6
  ) +
  coord_cartesian(ylim = c(0, 7)) +   # 👈 adjust limits here
  scale_fill_viridis_d() +   # colour-blind & print-friendly
  labs(
    x = "Start Type",
    y = expression("Velocity (m·s"^-1*")")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Previous Chain Duration Analysis ----
# 1. Get estimated marginal means for previous chain duration bucket

emm_duration <- emmeans(model, ~ previous.chain.duration.bucket, mode = "df.error")
emm_df_duration <- as.data.frame(emm_duration)

# 2. Compare each duration bucket to the grand mean
contrasts_duration <- contrast(emm_duration, method = "eff") %>%
  as.data.frame() %>%
  rename(
    Effect = estimate,
    StdError = SE
  ) %>%
  mutate(
    previous.chain.duration.bucket = gsub(" effect", "", contrast),
    Signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    )
  )

# 3. Prepare data for plotting
emm_df_labeled_duration <- emm_df_duration %>%
  left_join(contrasts_duration %>% select(previous.chain.duration.bucket, Signif), 
            by = "previous.chain.duration.bucket")

emm_df_labeled_duration$previous.chain.duration.bucket <- factor(
  emm_df_labeled_duration$previous.chain.duration.bucket,
  levels = c("< 5 sec", "5–10 sec", "10–15 sec", "15–20 sec", "20–25 sec", "≥ 25 sec"))


# 4. Plot
ggplot(emm_df_labeled_duration, aes(x = previous.chain.duration.bucket, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(aes(label = Signif, y = upper.CL + 0.05), size = 6) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Previous Chain Duration Bucket",
    y =  expression("Estimated Velocity (m·s"^-1*")"),
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.grid = element_blank(),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 5. Custom p-value formatter
format_custom_pval <- function(p, digits = 7) {
  ifelse(p < 1e-7, "<0.0000001", sprintf(paste0("%.", digits, "f"), p))
}

# 6. Table for academic reporting
table_df_duration <- emm_df_duration %>%
  left_join(
    contrasts_duration %>% select(previous.chain.duration.bucket, Effect, StdError, p.value, Signif),
    by = "previous.chain.duration.bucket"
  ) %>%
  mutate(
    Estimate = round(emmean, 3),
    `Effect Size` = round(Effect, 3),
    `Standard Error` = round(StdError, 3),
    `p-value` = format_custom_pval(p.value)
  ) %>%
  select(
    `Previous Chain Duration` = previous.chain.duration.bucket,
    Estimate,
    `Effect Size`,
    `Standard Error`,
    `p-value`,
    Signif
  )

sig_labels_duration <- contrasts_duration %>%
  select(previous.chain.duration.bucket, Signif)

ggplot(pre_window, aes(x = previous.chain.duration.bucket, y = velocity)) +
  geom_boxplot(aes(fill = previous.chain.duration.bucket), 
               alpha = 0.7, colour = "black", outlier.shape = NA) +
  # Optional: raw points
  # geom_jitter(width = 0.15, alpha = 0.2, size = 1) +
  geom_text(
    data = sig_labels_duration,
    aes(x = previous.chain.duration.bucket,
        y = 5.5,
        label = Signif),
    inherit.aes = FALSE,
    size = 6
  ) +
  coord_cartesian(ylim = c(0, 6)) +   # 👈 adjust limits here
  scale_fill_viridis_d() +   # colour-blind safe + print-friendly
  labs(
    x = "Previous Chain Duration Bucket",
    y = expression("Velocity (m·s"^-1*")")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



# Fixed Effect Interactions Visuals ----

library(emmeans)

# 1. Get interaction marginal means
emm_interactions <- emmeans(model, ~ start.zone * start.type * previous.chain.duration.bucket, mode = "df.error")
emm_df <- as.data.frame(emm_interactions)

# 2. Compare each interaction level to the grand mean
contrasts_df <- contrast(emm_interactions, method = "eff") %>%
  as.data.frame() %>%
  mutate(
    Signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    )
  )


contrast_df_clean <- contrasts_df %>%
  filter(grepl(" effect$", contrast)) %>%  # Keep only rows ending in "effect"
  mutate(
    contrast_clean = str_remove(contrast, " effect$"),
    start.zone = word(contrast_clean, 1),
    start.type = word(contrast_clean, 2),
    previous.chain.duration.bucket = str_trim(str_remove(contrast_clean, paste(start.zone, start.type)))
  )

# 2. Join with emmeans dataframe
emm_df_labeled <- emm_df %>%
  left_join(
    contrast_df_clean %>% 
      select(start.zone, start.type, previous.chain.duration.bucket, Signif),
    by = c("start.zone", "start.type", "previous.chain.duration.bucket")
  )


emm_df_labeled$previous.chain.duration.bucket <- factor(
  emm_df_labeled$previous.chain.duration.bucket,
  levels = c("< 5 sec", "5–10 sec", "10–15 sec", "15–20 sec", "20–25 sec", "≥ 25 sec")
)

emm_df_labeled$start.zone <- factor(
  emm_df_labeled$start.zone,
  levels = c("D50", "DM", "AM", "F50")
)

emm_df_labeled <- emm_df_labeled %>% filter(previous.chain.duration.bucket != "NA")

unique(emm_df_labeled$previous.chain.duration.bucket)

# 4. Plot

ggplot(emm_df_labeled, aes(x = start.zone, y = emmean, shape = start.type, group = start.type)) +
  geom_point(position = position_dodge(width = 0.5), size = 3, color = "black") +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.2, position = position_dodge(width = 0.5),
    color = "black"
  ) +
  geom_text(
    aes(label = Signif, y = upper.CL + 0.05),
    position = position_dodge(width = 0.5), size = 4.5, fontface = "bold", color = "black"
  ) +
  facet_wrap(~ previous.chain.duration.bucket, ncol = 3) +
  scale_shape_manual(values = c(16, 17)) +  # Circle and triangle
  theme_classic(base_size = 14) +
  labs(
    # title = "Estimated Velocity Prior to Turnover",
    #  subtitle = "By Start Zone, Start Type and Previous Chain Duration",
    # caption = "* p < 0.05   ** p < 0.01   *** p < 0.001",
    x = "Start Zone",
    y = expression("Estimated Velocity (m·s"^-1*")"),
    shape = "Start Type"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, face = "italic", size = 10),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )


ggplot(emm_df_labeled, aes(x = start.zone, y = emmean, shape = start.type, group = start.type)) +
  geom_point(position = position_dodge(width = 0.5), size = 3, color = "black") +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.2, position = position_dodge(width = 0.5),
    color = "black"
  ) +
  geom_text(
    aes(label = Signif, y = upper.CL + 0.05),
    position = position_dodge(width = 0.5), size = 4.5, fontface = "bold", color = "black"
  ) +
  facet_wrap(~ previous.chain.duration.bucket, ncol = 3) +
  scale_shape_manual(values = c(16, 17)) +  # Circle and triangle
  theme_classic(base_size = 14) +
  labs(
    # title = "Estimated Velocity Prior to Turnover",
    #  subtitle = "By Start Zone, Start Type and Previous Chain Duration",
    # caption = "* p < 0.05   ** p < 0.01   *** p < 0.001",
    x = "Start Zone",
    y = expression("Estimated Velocity (m·s"^-1*")"),
    shape = "Start Type"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, face = "italic", size = 10),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )


ggplot(
  emm_df_labeled, 
  aes(x = start.zone, y = emmean, 
      shape = start.type, 
      group = interaction(start.type, previous.chain.duration.bucket),
      colour = previous.chain.duration.bucket)  # 👈 colour by chain duration
) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.2, position = position_dodge(width = 0.5)
  ) +
  geom_text(
    aes(label = Signif, y = upper.CL + 0.05),
    position = position_dodge(width = 0.5), 
    size = 4.5, fontface = "bold"
  ) +
  scale_shape_manual(values = c(16, 17)) +  # Circle and triangle
  theme_classic(base_size = 14) +
  labs(
    x = "Start Zone",
    y = expression("Estimated Velocity (m·s"^-1*")"),
    shape = "Start Type",
    colour = "Prev. Chain Duration"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, face = "italic", size = 10),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )











#### Chain Outcome Marginal Means ----

library(emmeans)

# Estimated marginal means
emm_chain_outcome <- emmeans(model, ~ chain.outcome, mode = "df.error")

# Grand mean contrasts with significance labels
contrasts_grand_mean <- contrast(emm_chain_outcome, method = "eff") %>%
  as.data.frame() %>%
  mutate(
    chain.outcome = gsub(" effect", "", contrast),
    Signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    )
  )

# Estimated marginal means to data frame
emm_df <- as.data.frame(emm_chain_outcome)

# Join only by start.zone
emm_df_labeled <- emm_df %>%
  left_join(contrasts_grand_mean %>% select(chain.outcome, Signif), by = "chain.outcome")

# Optional: set order for plot
emm_df_labeled$chain.outcome <- factor(emm_df_labeled$chain.outcome, levels = c("No Score", "Score"))

# Plot with significance labels
ggplot(emm_df_labeled, aes(x = chain.outcome, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(aes(label = Signif, y = emmean + 0.05), size = 6) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Chain Outcome",
    y = "Estimated Velocity (m/s^2)",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, face = "italic")
  )


# New version for table

# Get estimated marginal means
emm_chain_outcome <- emmeans(model, ~ chain.outcome, mode = "df.error")
emm_df <- as.data.frame(emm_chain_outcome)

# Get contrasts vs grand mean and calculate significance
contrasts_grand_mean <- contrast(emm_chain_outcome, method = "eff") %>%
  as.data.frame() %>%
  rename(
    Effect = estimate,
    StdError = SE
  ) %>%
  mutate(
    chain.outcome = gsub(" effect", "", contrast),
    Signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    )
  )

# Custom p-value formatter
format_custom_pval <- function(p, digits = 7) {
  ifelse(p < 1e-7, "<0.0000001", sprintf(paste0("%.", digits, "f"), p))
}

# Merge EMMs with contrasts for final table
table_df <- emm_df %>%
  left_join(
    contrasts_grand_mean %>%
      select(chain.outcome, Effect, StdError, p.value, Signif),
    by = "chain.outcome"
  ) %>%
  mutate(
    Estimate = round(emmean, 3),
    `Effect Size` = round(Effect, 3),
    `Standard Error` = round(StdError, 3),
    `p-value` = format_custom_pval(p.value)
  ) %>%
  select(
    `Start Zone` = chain.outcome,
    Estimate,
    `Effect Size`,
    `Standard Error`,
    `p-value`,
    Signif
  )

write_csv(visuals, "Turnover ID in Dataset.csv")





### Plots next to each other 

library(patchwork)


p1 <- ggplot(pre_window, aes(x = start.type, y = velocity)) +
  geom_boxplot(aes(fill = start.type), alpha = 0.7, colour = "black", outlier.shape = NA) +
  geom_text(
    data = sig_labels_type,
    aes(x = start.type, y = 6.5,
        label = Signif),
    inherit.aes = FALSE,
    size = 6
  ) +
  scale_fill_viridis_d() +
  labs(x = "Start Type", y = expression("Velocity (m·s"^-1*")")) +
  coord_cartesian(ylim = c(0, 7)) +   # 👈 adjust limits here
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

p2 <- ggplot(pre_window, aes(x = start.zone, y = velocity)) +
  geom_boxplot(aes(fill = start.zone), alpha = 0.5, outlier.shape = NA) +
  geom_text(
    data = sig_labels,
    aes(x = start.zone, y = 6.5,
        label = Signif),
    inherit.aes = FALSE,
    size = 6
  ) +
  scale_fill_viridis_d() +
  labs(x = "Start Zone", y = expression("Velocity (m·s"^-1*")")) +
  coord_cartesian(ylim = c(0, 7)) +   # 👈 same range for consistency
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Combine side by side (patchwork)
p1 + p2

