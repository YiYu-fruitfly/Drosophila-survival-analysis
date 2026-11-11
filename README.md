# Drosophila-survival-analysis
Statistical analysis of the Drosophila survival using stratified Cox proportional hazards model
# ==============================
# 0. Load packages
# ==============================

library(survival)
library(broom)
library(dplyr)
library(readxl)
library(openxlsx)

# ==============================
# 1. Read data
# ==============================

survival_name <- "survival_data"
dir_path <- "file destination"

flydata <- read_excel(
  paste0(dir_path, "file destination", survival_name, ".xlsx")
)

# ==============================
# 2. Clean data
# ==============================

flydata <- flydata %>%
  mutate(
    Event = as.numeric(Event),
    Time = as.numeric(Time),
    Treatment = factor(Treatment),
    Repeat = factor(Repeat)
  )

flydata$Treatment <- relevel(flydata$Treatment, ref = "control")

# ==============================
# 3. PRIMARY MODEL
# Stratified Cox (biological repeats)
# ==============================

cox_strata <- coxph(
  Surv(Time, Event) ~ Treatment + strata(Repeat),
  data = flydata
)

summary(cox_strata)

# ==============================
# 4. Extract HR table (Nature style)
# ==============================

cox_tidy <- tidy(
  cox_strata,
  exponentiate = TRUE,
  conf.int = TRUE
)

cox_table <- cox_tidy %>%
  filter(term != "strata(Repeat)") %>%
  mutate(
    Treatment = gsub("Treatment", "", term),
    HR = round(estimate, 2),
    CI = paste0(
      round(conf.low, 2), " – ", round(conf.high, 2)
    ),
    P_value = ifelse(
      p.value < 0.001,
      formatC(p.value, format = "e", digits = 2),
      round(p.value, 3)
    )
  ) %>%
  select(Treatment, HR, CI, P_value)

cox_table

# ==============================
# 5. PH assumption test
# ==============================

ph_test <- cox.zph(cox_strata)
print(ph_test)
plot(ph_test)

# ==============================
# 6. Export Excel report
# ==============================

wb <- createWorkbook()

addWorksheet(wb, "Stratified Cox HR")
writeData(wb, "Stratified Cox HR", cox_table)

addWorksheet(wb, "Model summary")
writeData(
  wb,
  "Model summary",
  data.frame(
    N = nrow(flydata),
    Events = sum(flydata$Event),
    Concordance = summary(cox_strata)$concordance[1],
    Wald_test_p = summary(cox_strata)$waldtest["pvalue"]
  )
)

saveWorkbook(
  wb,
  paste0(
    dir_path,
    "/Statistics/result_",
    survival_name,
    "_",
    Sys.Date(),
    ".xlsx"
  ),
  overwrite = TRUE
)

