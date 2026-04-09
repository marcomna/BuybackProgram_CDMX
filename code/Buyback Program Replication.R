## =============================================================================
## Synthetic Control Analysis: Gun Buy-Back Programme and COVID-19 Lockdown
## Effects on Firearm Crime in Mexico City
##
## This replication script estimates the effects of two distinct interventions
## on firearm crime incidence across Mexican states. The first treatment
## corresponds to the start of the gun buy-back programme in January 2019.
## The second corresponds to the onset of COVID-19 restrictions in March 2020.
##
## For each intervention, a synthetic control is constructed using pre-treatment
## data to match Mexico City's crime trajectory. Robustness checks include
## weight distributions, placebo tests, MSPE ratios, and augmented synthetic
## control (augsynth) estimates for each post-treatment phase.
##
## Sections:
##   1.  Load packages
##   2.  Data loading and processing
##   3.  Panel construction
##   4.  Main synthetic control (armas_conUD)
##   5.  Augmented synthetic control — phase-wise ATT (armas_conUD)
##   6.  Secondary outcomes: SC and SynthDiD by crime type
##       6.1  Firearm homicides (homi_af_rate)
##       6.2  Firearm injuries (lesi_dolosa_af)
##       6.3  Non-firearm homicides (homi_nofire_rate)
##       6.4  Non-firearm injuries (lesi_nofire_rate)
##   7.  Phase-wise ATT for all secondary outcomes
##   8.  Firearm buyback analyses
##       8.1  SC and SynthDiD — total buybacks (tasa_total)
##   9.  Firearm seizures — staggered DiD
##  10.  Validation against INEGI homicide records
##  11.  Descriptive buyback charts
##
## Authors: Marco Méndez Atienza and Aurora Ramírez Álvarez
## =============================================================================


## =============================================================================
## 1. Load packages
## =============================================================================

library(pacman)
p_load(
  tidyverse,    # Data manipulation and visualization
  lubridate,    # Date handling
  reshape2,     # Wide-to-long reshaping (dcast)
  zoo,          # Year-month date handling
  tidysynth,    # Synthetic control methods
  foreign,      # Read DBF files
  readr,        # Efficient CSV/text reading
  readxl,       # Excel files
  ggthemes,     # Clean ggplot themes
  rlang,
  gridExtra,
  grid,
  augsynth,     # Augmented synthetic control
  synthdid,     # Synthetic DiD (Arkhangelsky et al. 2021)
  fixest,       # Two-way FE and event-study DiD
  broom,        # Tidy model output
  kableExtra,   # Table formatting
  stringi,
  stringdist,
  fuzzyjoin,
  stringr,
  sf,
  wesanderson,
  scales,
  forcats
)


## =============================================================================
## 2. Data loading and processing
## =============================================================================

# NOTE: Update the DATA_DIR path to point to your local data folder.
DATA_DIR <- "path/to/your/data"   # <-- edit this


# -----------------------------------------------------------------------------
# 2.1  Firearm crime data (IDEFC) — SESNSP repository
#      Keep 2017–2023; recode long state names to short forms.
# -----------------------------------------------------------------------------
delitos_raw <- read_csv(
  file.path(DATA_DIR, "IDEFC_NM_abr24.csv"),
  locale = locale(encoding = "UTF-8")
) %>%
  filter(Año > 2016, Año < 2024) %>%
  mutate(Entidad = recode(Entidad,
    "Coahuila de Zaragoza"            = "Coahuila",
    "Michoacán de Ocampo"             = "Michoacán",
    "Veracruz de Ignacio de la Llave" = "Veracruz"
  ))

# Pivot wide monthly columns to long format and construct a Date per observation.
delitos <- delitos_raw %>%
  gather(mes, valor, -c(Año, Clave_Ent, Entidad,
                        `Bien jurídico afectado`, `Tipo de delito`,
                        `Subtipo de delito`, Modalidad)) %>%
  mutate(
    mes_num = recode(mes,
      "Enero" = 1, "Febrero" = 2, "Marzo" = 3, "Abril" = 4,
      "Mayo" = 5, "Junio" = 6, "Julio" = 7, "Agosto" = 8,
      "Septiembre" = 9, "Octubre" = 10, "Noviembre" = 11, "Diciembre" = 12
    ),
    Fecha = as.Date(paste(Año, mes_num, "01", sep = "-"))
  )


# -----------------------------------------------------------------------------
# 2.2  Population (2020 census) — used to compute rates per 100,000 inhabitants
# -----------------------------------------------------------------------------
poblacion <- read_csv(
  file.path(DATA_DIR, "poblacion.csv"),
  locale = locale(encoding = "UTF-8")
) %>%
  mutate(Entidad = recode(Entidad,
    "Coahuila de Zaragoza"            = "Coahuila",
    "Estado de México"                = "México",
    "Michoacán de Ocampo"             = "Michoacán",
    "Veracruz de Ignacio de la Llave" = "Veracruz"
  ))


# -----------------------------------------------------------------------------
# 2.3  Unclassified deaths ("Evento no especificado") from INEGI defunction files
#      DBF files are read for 2017–2023, transformed to keep state, cause,
#      month and year, then joined with the cause catalogue.
# -----------------------------------------------------------------------------
dbf_paths <- list(
  file.path(DATA_DIR, "Defunciones/2023/DEFUN23.dbf"),
  file.path(DATA_DIR, "Defunciones/2022/DEFUN22.dbf"),
  file.path(DATA_DIR, "Defunciones/2021/defun21.dbf"),
  file.path(DATA_DIR, "Defunciones/2020/defun20.dbf"),
  file.path(DATA_DIR, "Defunciones/2019/DEFUN19.dbf"),
  file.path(DATA_DIR, "Defunciones/2018/DEFUN18.dbf"),
  file.path(DATA_DIR, "Defunciones/2017/DEFUN17.dbf")
)

process_dbf <- function(fp) {
  read.dbf(fp, as.is = TRUE) %>%
    mutate(across(where(is.character), ~ iconv(., from = "latin1", to = "UTF-8"))) %>%
    select(ENT_REGIS, CAUSA_DEF, MES_REGIS, ANIO_REGIS)
}

defunciones_all <- lapply(dbf_paths, process_dbf) %>% bind_rows()

# Cause catalogue
categorias <- read.dbf(
  file.path(DATA_DIR, "Defunciones/2022/CATMINDE.dbf"),
  as.is = TRUE
) %>%
  mutate(across(where(is.character), ~ iconv(., from = "latin1", to = "UTF-8")))

# State code → name lookup
entidades_map <- c(
  "01" = "Aguascalientes",    "02" = "Baja California",   "03" = "Baja California Sur",
  "04" = "Campeche",          "05" = "Coahuila",          "06" = "Colima",
  "07" = "Chiapas",           "08" = "Chihuahua",         "09" = "Ciudad de México",
  "10" = "Durango",           "11" = "Guanajuato",        "12" = "Guerrero",
  "13" = "Hidalgo",           "14" = "Jalisco",           "15" = "México",
  "16" = "Michoacán",         "17" = "Morelos",           "18" = "Nayarit",
  "19" = "Nuevo León",        "20" = "Oaxaca",            "21" = "Puebla",
  "22" = "Querétaro",         "23" = "Quintana Roo",      "24" = "San Luis Potosí",
  "25" = "Sinaloa",           "26" = "Sonora",            "27" = "Tabasco",
  "28" = "Tamaulipas",        "29" = "Tlaxcala",          "30" = "Veracruz",
  "31" = "Yucatán",           "32" = "Zacatecas"
)

defunciones_all <- defunciones_all %>%
  left_join(categorias, by = c("CAUSA_DEF" = "CVE")) %>%
  mutate(
    Entidad = recode(ENT_REGIS, !!!entidades_map),
    Fecha   = as.Date(paste(ANIO_REGIS, MES_REGIS, "01", sep = "-"))
  )

defunciones_noespec <- defunciones_all %>%
  filter(str_detect(DESCRIP, "Evento no especificado")) %>%
  group_by(Entidad, Fecha) %>%
  summarise(Cuenta = n(), .groups = "drop")


# -----------------------------------------------------------------------------
# 2.4  Unemployment rate (quarterly → monthly expansion)
#      Source: Encuesta Nacional de Ocupación y Empleo (ENOE), INEGI
# -----------------------------------------------------------------------------
ruta_unemp  <- file.path(DATA_DIR, "Desoc_EF.CSV")
unemp_lines <- read_lines(ruta_unemp)
header_idx  <- which(str_detect(unemp_lines, "Entidad federativa"))[1]

unemp_data <- read_csv(
  ruta_unemp,
  skip   = header_idx - 1,
  na     = c("ND", "NA"),
  locale = locale(decimal_mark = ".")
) %>%
  filter(Indicador == "Total",
         `Entidad federativa` != "Estados Unidos Mexicanos")

unemp_cols <- c(
  "Entidad federativa",
  paste0(rep(2017:2023, each = 4), "/", sprintf("%02d", 1:4))
)

unemp_data <- unemp_data %>%
  select(any_of(unemp_cols)) %>%
  rename(Entidad = `Entidad federativa`) %>%
  mutate(
    Entidad = str_replace_all(Entidad, "Ciudad de M.*xico", "Ciudad de México"),
    across(-Entidad, as.numeric)
  )

unemp_monthly <- unemp_data %>%
  pivot_longer(cols = -Entidad, names_to = "period", values_to = "unemp_rate") %>%
  mutate(
    year        = as.integer(str_extract(period, "^\\d{4}")),
    quarter     = as.integer(str_extract(period, "(?<=/)\\d{2}")),
    month_start = map_int(quarter, ~ if (.x == 1) 1L else if (.x == 2) 4L
                          else if (.x == 3) 7L else 10L)
  ) %>%
  group_by(Entidad, year, unemp_rate, month_start) %>%
  reframe(month = month_start + 0:2) %>%
  mutate(Fecha = as.Date(paste(year, month, "01", sep = "-"))) %>%
  select(Entidad, Fecha, unemp_rate)


# -----------------------------------------------------------------------------
# 2.5  Industrial Production Index (IMAIEF) — monthly state-level index
#      Source: INEGI
# -----------------------------------------------------------------------------
ipi_raw <- read_csv(
  file.path(DATA_DIR, "Act. Ind. Estatal/IMAIEF VF.csv"),
  locale = locale(encoding = "UTF-8")
)

ipi_long <- ipi_raw %>%
  separate(col = Descriptores, into = c("descriptor", "Entidad"),
           sep = "\\|", extra = "merge", fill = "right") %>%
  filter(!is.na(Entidad)) %>%
  pivot_longer(cols = -c(descriptor, Entidad),
               names_to = "period", values_to = "value") %>%
  mutate(
    period    = str_remove(period, "<.*>"),
    is_annual = str_detect(period, "Anual"),
    period    = if_else(is_annual, NA_character_, period)
  ) %>%
  filter(!is.na(period)) %>%
  separate(period, into = c("year", "month_name"), sep = "\\|") %>%
  mutate(
    year       = as.integer(year),
    month_name = str_trim(month_name),
    month      = match(month_name,
                       c("Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio",
                         "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre",
                         "Diciembre")),
    Fecha      = as.Date(paste(year, month, "01", sep = "-")),
    value      = as.numeric(value)
  )

ipi_long$Entidad <- str_replace_all(ipi_long$Entidad, c(
  "Estados Unidos Mexicanos" = "Nacional",
  "Ciudad de M.*xico"        = "Ciudad de México"
))

ipi_values <- ipi_long %>%
  filter(str_detect(descriptor, "Índice")) %>%
  select(Entidad, Fecha, ipi = value)

ipi_pct <- ipi_long %>%
  filter(str_detect(descriptor, "Variación")) %>%
  select(Entidad, Fecha, ipi_pct_change_month = value)

ipi_monthly <- left_join(ipi_values, ipi_pct, by = c("Entidad", "Fecha"))


## =============================================================================
## 3. Panel construction
## =============================================================================

# Aggregate crime counts by state, modality and month; pivot to wide format;
# compute total crimes; then merge population, unclassified deaths,
# unemployment and IPI. Finally convert counts to rates per 100,000.

delitosCS <- delitos %>%
  group_by(Fecha, Entidad, Modalidad) %>%
  summarise(Total = sum(valor, na.rm = TRUE), .groups = "drop") %>%
  dcast(Fecha + Entidad ~ Modalidad) %>%
  mutate(
    TotalDelitos = rowSums(across(where(is.numeric)), na.rm = TRUE),
    armas_sinUD  = `Con arma de fuego`
  ) %>%
  left_join(poblacion,          by = "Entidad") %>%
  left_join(defunciones_noespec, by = c("Fecha", "Entidad")) %>%
  mutate(
    Cuenta      = coalesce(Cuenta, 0),
    armas_conUD = armas_sinUD + Cuenta
  ) %>%
  select(-Cuenta) %>%
  left_join(unemp_monthly, by = c("Entidad", "Fecha")) %>%
  left_join(ipi_monthly,   by = c("Entidad", "Fecha")) %>%
  mutate(
    unemp_rate                               = coalesce(unemp_rate, 0),
    ipi                                      = coalesce(ipi, 0),
    ipi_pct_change_month                     = coalesce(ipi_pct_change_month, 0),
    `Allanamiento de morada`                 = coalesce(`Allanamiento de morada`, 0),
    `Robo de coche de 4 ruedas Con violencia`= coalesce(`Robo de coche de 4 ruedas Con violencia`, 0),
    Narcomenudeo                             = coalesce(Narcomenudeo, 0),
    `Secuestro extorsivo`                    = coalesce(`Secuestro extorsivo`, 0),
    Amenazas                                 = coalesce(Amenazas, 0)
  )

# Convert counts to rates per 100,000 inhabitants
numeric_cols  <- delitosCS %>% select(where(is.numeric)) %>% colnames()
cols_to_scale <- setdiff(numeric_cols,
                         c("Poblacion2020", "unemp_rate", "ipi", "ipi_pct_change_month"))

delitosCS <- delitosCS %>%
  mutate(across(all_of(cols_to_scale), ~ . * 100000 / Poblacion2020)) %>%
  distinct(Entidad, Fecha, .keep_all = TRUE) %>%
  mutate(mes = factor(month(Fecha)))

# Donor pool: exclude neighboring states to avoid spatial spillovers
donor_exclude <- c("Morelos", "México")


## =============================================================================
## 4. Main synthetic control — total firearm crime (armas_conUD)
## =============================================================================

main_CS <- delitosCS %>%
  filter(!Entidad %in% donor_exclude) %>%
  synthetic_control(
    outcome           = armas_conUD,
    unit              = Entidad,
    time              = Fecha,
    i_unit            = "Ciudad de México",
    i_time            = as.Date("2019-01-01"),
    generate_placebos = TRUE
  ) %>%
  # --- Baseline predictors (2017–2018 pre-treatment average) ---
  generate_predictor(
    time_window = c(as.Date("2017-01-01"), as.Date("2018-12-01")),
    trespassing       = mean(`Allanamiento de morada`, na.rm = TRUE),
    grand_theft_auto  = mean(`Robo de coche de 4 ruedas Con violencia`, na.rm = TRUE),
    unemployment_rate = mean(unemp_rate, na.rm = TRUE),
    ipi               = mean(ipi, na.rm = TRUE),
    ipi_pct_change    = mean(ipi_pct_change_month, na.rm = TRUE),
    drug_dealing      = mean(Narcomenudeo, na.rm = TRUE),
    kidnapping        = mean(`Secuestro extorsivo`, na.rm = TRUE),
    weapon            = mean(`Con arma blanca`, na.rm = TRUE),
    tractors          = mean(`Robo de tractores Con violencia`, na.rm = TRUE),
    threats           = mean(Amenazas, na.rm = TRUE),
    family_violence   = mean(`Violencia familiar`, na.rm = TRUE),
    property_damage   = mean(`Daño a la propiedad`, na.rm = TRUE),
    bodily_harm       = mean(`Otros delitos que atentan contra la vida y la integridad corporal`,
                             na.rm = TRUE)
  ) %>%
  # --- Monthly snapshots within the pre-treatment window ---
  generate_predictor(time_window = as.Date("2017-01-01"), FC0117 = armas_conUD) %>%
  generate_predictor(time_window = as.Date("2017-05-01"), FC0517 = armas_conUD) %>%
  generate_predictor(time_window = as.Date("2018-09-01"), FC0917 = armas_conUD) %>%
  generate_predictor(time_window = as.Date("2018-01-01"), FC0118 = armas_conUD) %>%
  generate_predictor(time_window = as.Date("2018-05-01"), FC0518 = armas_conUD) %>%
  generate_predictor(time_window = as.Date("2018-09-01"), FC0918 = armas_conUD) %>%
  # --- Quarterly averages (2017–2018) ---
  generate_predictor(time_window = c(as.Date("2017-01-01"), as.Date("2017-03-01")),
                     q1_2017 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2017-04-01"), as.Date("2017-06-01")),
                     q2_2017 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2017-07-01"), as.Date("2017-09-01")),
                     q3_2017 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2017-10-01"), as.Date("2017-12-01")),
                     q4_2017 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2018-01-01"), as.Date("2018-03-01")),
                     q1_2018 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2018-04-01"), as.Date("2018-06-01")),
                     q2_2018 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2018-07-01"), as.Date("2018-09-01")),
                     q3_2018 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2018-10-01"), as.Date("2018-12-01")),
                     q4_2018 = mean(armas_conUD, na.rm = TRUE)) %>%
  # --- Weight optimization ---
  generate_weights(
    optimization_window = c(as.Date("2017-01-01"), as.Date("2018-12-01")),
    margin_ipop = 0.02, sigf_ipop = 7, bound_ipop = 5
  ) %>%
  generate_control()


# --- Extract observed vs. synthetic series for Mexico City ---
principalCSdelitos <- main_CS %>%
  unnest(cols = c(".synthetic_control")) %>%
  filter(.id == "Ciudad de México", .type == "treated") %>%
  select(Fecha = time_unit, real_y, synth_y) %>%
  mutate(diff = real_y - synth_y) %>%
  left_join(
    delitosCS %>% filter(Entidad == "Ciudad de México") %>% select(Fecha, armas_sinUD),
    by = "Fecha"
  )


# --- Plot 4.1: Observed vs. Synthetic firearm crime rate ---
ggplot(principalCSdelitos, aes(x = Fecha)) +
  geom_line(aes(y = real_y,  color = "Observed"),  linewidth = 1.2) +
  geom_line(aes(y = synth_y, color = "Synthetic"), linetype = "longdash", linewidth = 1.1) +
  geom_vline(xintercept = as.Date("2019-01-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-03-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-12-01"), linetype = "dotted", linewidth = 1) +
  annotate("text", x = as.Date("2018-12-15"), y = 0.6,
           label = "Program starts",   angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = as.Date("2020-02-15"), y = 0.6,
           label = "Lockdown starts",  angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = as.Date("2020-11-15"), y = 0.6,
           label = "Lockdown ends",    angle = 90, vjust = -0.5, size = 3) +
  scale_color_manual(values = c("Observed" = "gray30", "Synthetic" = "#449DD1")) +
  scale_x_date(date_breaks = "4 months", date_labels = "%m/%Y",
               expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Date", y = "Rate per 100,000 inhabitants", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title  = element_text(face = "bold"),
        legend.position = "bottom")


# --- Plot 4.2: Donor and predictor weights ---
unit_weights <- main_CS %>%
  filter(.id == "Ciudad de México", .type == "treated") %>%
  select(.unit_weights) %>%
  unnest(cols = c(.unit_weights)) %>%
  filter(weight > 0) %>%
  mutate(unit = factor(unit, levels = unit[order(weight)]))

predictor_weights <- main_CS %>%
  filter(.id == "Ciudad de México", .type == "treated") %>%
  select(.predictor_weights) %>%
  unnest(cols = c(.predictor_weights)) %>%
  filter(weight > 0) %>%
  mutate(variable = factor(variable, levels = variable[order(weight)]))

p_units <- ggplot(unit_weights, aes(x = unit, y = weight)) +
  geom_bar(stat = "identity", fill = "#08306B", color = "black", width = 1, linewidth = 0.4) +
  labs(x = "Donor State", y = "Weight") +
  theme_minimal(base_size = 12) +
  theme(axis.title = element_text(size = 12, face = "bold")) +
  coord_flip()

p_predictors <- ggplot(predictor_weights, aes(x = variable, y = weight)) +
  geom_bar(stat = "identity", fill = "#6BAED6", color = "black", width = 1, linewidth = 0.4) +
  labs(x = "Predictor", y = "Weight") +
  theme_minimal(base_size = 12) +
  theme(axis.title = element_text(size = 12, face = "bold")) +
  coord_flip()

grid.arrange(
  p_units, p_predictors, nrow = 1,
  top = textGrob(
    "Optimal weights: donor states and predictors — synthetic control for Mexico City",
    gp = gpar(fontsize = 14, fontface = "bold")
  )
)


# --- Table 4.1: Predictor balance ---
main_CS %>%
  grab_balance_table() %>%
  kable(
    format = "pandoc", booktabs = TRUE, digits = 2,
    col.names = c("Predictor", "Treated (CDMX)", "Synthetic Control", "Donor Pool"),
    caption = "Predictor balance between Mexico City and its synthetic control."
  ) %>%
  kable_styling(latex_options = c("HOLD_position", "scale_down"))


# --- Plot 4.3: Placebo test — gap paths (pruned donor pool) ---
MSPEcrimes <- main_CS %>%
  grab_significance() %>%
  mutate(RMSPEpre = sqrt(pre_mspe), RMSPEpost = sqrt(post_mspe))

cdmx_rmspe  <- MSPEcrimes %>% filter(unit_name == "Ciudad de México") %>% pull(RMSPEpre)
valid_units  <- MSPEcrimes %>% filter(RMSPEpre / cdmx_rmspe < 2.5) %>% pull(unit_name)

CSplacebos <- main_CS %>%
  filter(.type == "treated", .id %in% valid_units) %>%
  select(.id, .synthetic_control) %>%
  unnest(cols = c(.synthetic_control)) %>%
  rename(State = .id, Date = time_unit, Synthetic = synth_y, Observed = real_y) %>%
  mutate(gap = Observed - Synthetic)

ggplot(CSplacebos, aes(x = Date, y = gap, group = State)) +
  geom_line(data = filter(CSplacebos, State != "Ciudad de México"),
            color = "gray", linewidth = 1, alpha = 0.35) +
  geom_line(data = filter(CSplacebos, State == "Ciudad de México"),
            color = "gray30", linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2019-01-01"), linetype = "dotted", linewidth = 1) +
  annotate("text", x = as.Date("2021-01-01"), y = -3,
           label = "Mexico City", size = 4, fontface = "bold", color = "gray30") +
  scale_x_date(date_breaks = "4 months", date_labels = "%m/%Y") +
  labs(
    x = "Date", y = "Gap (Observed – Synthetic)",
    caption = "\nNote: Gray lines show donor-state placebo gaps (pre-RMSPE < 2.5 × CDMX RMSPE)."
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 65, hjust = 1),
        axis.title   = element_text(size = 12, face = "bold"),
        plot.caption = element_text(hjust = 0, size = 8))


# --- Plot 4.4: RMSPE ratio (Post/Pre) bar chart ---
MSPEcrimes <- MSPEcrimes %>%
  mutate(
    mspe_ratio = RMSPEpost / RMSPEpre,
    highlight  = if_else(unit_name == "Ciudad de México", "Mexico City", "Donor Units"),
    unit_label = if_else(unit_name == "Ciudad de México", "Mexico City", unit_name),
    unit_label = factor(unit_label, levels = unit_label[order(mspe_ratio)])
  )

ggplot(MSPEcrimes, aes(x = mspe_ratio, y = unit_label, fill = highlight)) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  scale_fill_manual(values = c("Mexico City" = "gray30", "Donor Units" = "#449DD1")) +
  labs(x = "RMSPE Ratio (Post/Pre)", y = "State") +
  theme_minimal(base_size = 12) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9),
        legend.position = "none")


# --- Table 4.2: RMSPE-based placebo p-values (pruned pool) ---
tbl_rmspe <- MSPEcrimes %>%
  filter(unit_name == "Ciudad de México" | (RMSPEpre / cdmx_rmspe < 2.5)) %>%
  transmute(
    Unit        = if_else(unit_name == "Ciudad de México", "Mexico City", unit_name),
    Type        = if_else(unit_name == "Ciudad de México", "Treated", "Donor"),
    Pre_RMSPE   = RMSPEpre,
    Post_RMSPE  = RMSPEpost,
    RMSPE_Ratio = mspe_ratio,
    Rank        = rank(-mspe_ratio, ties.method = "min")
  )

ratios   <- setNames(tbl_rmspe$RMSPE_Ratio, tbl_rmspe$Unit)
N        <- nrow(tbl_rmspe)
fisher_p <- sapply(seq_len(N), function(i) {
  others <- ratios[-i]
  (sum(others >= ratios[i]) + 1) / (length(others) + 1)
})
z_scores <- sapply(seq_len(N), function(i) {
  others <- ratios[-i]
  (ratios[i] - mean(others)) / sd(others)
})

tbl_rmspe %>%
  mutate(
    `Pre-RMSPE`      = round(Pre_RMSPE,   6),
    `Post-RMSPE`     = round(Post_RMSPE,  6),
    `RMSPE Ratio`    = round(RMSPE_Ratio, 6),
    `Fisher p-value` = round(fisher_p,    4),
    `Z-score`        = round(z_scores,    3)
  ) %>%
  select(Unit, Type, `Pre-RMSPE`, `Post-RMSPE`, `RMSPE Ratio`, Rank, `Fisher p-value`, `Z-score`) %>%
  arrange(Rank) %>%
  kable(caption = "RMSPE-based placebo p-values and z-scores (pruned donor pool)")


## =============================================================================
## 5. Augmented synthetic control — phase-wise ATT (armas_conUD)
## =============================================================================

# --- 5.1  Pre-treatment gap (2017–2018) via bootstrap ---
weights_main <- main_CS %>%
  grab_unit_weights() %>%
  filter(unit != "Ciudad de México")

panel_pre_main <- delitosCS %>%
  filter(Fecha >= as.Date("2017-01-01"), Fecha <= as.Date("2018-12-01")) %>%
  select(Entidad, Fecha, armas_conUD) %>%
  pivot_wider(names_from = Entidad, values_from = armas_conUD) %>%
  arrange(Fecha)

treated_y_pre  <- panel_pre_main$`Ciudad de México`
control_mat_pre <- panel_pre_main %>%
  select(all_of(weights_main$unit)) %>%
  as.matrix()

synth_y_pre <- control_mat_pre %*% weights_main$weight
gap_pre      <- treated_y_pre - synth_y_pre

set.seed(123)
boot_means_main <- replicate(1000, mean(sample(gap_pre, replace = TRUE), na.rm = TRUE))

avg_gap_pre <- mean(gap_pre, na.rm = TRUE)
ci_gap_pre  <- quantile(boot_means_main, probs = c(0.025, 0.975))

cat("Pre-treatment gap:", round(avg_gap_pre, 4),
    "| 95% CI:", round(ci_gap_pre[1], 4), "–", round(ci_gap_pre[2], 4), "\n")


# --- Helper: build augsynth panel and compute phase ATT ---
build_augsynth_phase <- function(data, outcome_var, t_start, t_end, t_treat) {
  panel <- data %>%
    mutate(
      id       = Entidad,
      date     = Fecha,
      treated  = if_else(id == "Ciudad de México" & date >= t_treat, 1, 0),
      y        = .data[[outcome_var]],
      time_num = year(date) * 100 + month(date)
    ) %>%
    filter(date >= t_start, date <= t_end) %>%
    select(id, date, time_num, treated, y)

  model <- augsynth(
    y ~ treated,
    unit       = id,
    time       = time_num,
    data       = panel,
    t_int      = year(t_treat) * 100 + month(t_treat),
    progfunc   = "none",
    fixedeff   = FALSE,
    lambda     = "none",
    inf_method = "conformal"
  )

  gap <- model$data$y[model$data$trt == 1] -
    as.numeric(t(model$weights) %*% model$data$y[model$data$trt == 0, ])

  att <- mean(gap, na.rm = TRUE)
  se  <- sd(gap, na.rm = TRUE) / sqrt(length(gap))
  ci  <- att + c(-1, 1) * qnorm(0.975) * se
  list(model = model, att = att, ci = ci)
}

# --- Phase estimates ---
a1_main <- build_augsynth_phase(delitosCS, "armas_conUD",
  as.Date("2017-01-01"), as.Date("2020-03-01"), as.Date("2019-01-01"))
a2_main <- build_augsynth_phase(delitosCS, "armas_conUD",
  as.Date("2019-01-01"), as.Date("2020-12-01"), as.Date("2020-03-01"))
a3_main <- build_augsynth_phase(delitosCS, "armas_conUD",
  as.Date("2020-03-01"), as.Date("2023-12-01"), as.Date("2020-12-01"))

# --- Table 5.1: ATT summary ---
att_segments <- tibble(
  start = as.Date(c("2017-01-01", "2019-01-01", "2020-03-01", "2020-12-01")),
  end   = as.Date(c("2018-12-01", "2020-03-01", "2020-12-01", "2023-12-01")),
  att   = c(avg_gap_pre,    a1_main$att, a2_main$att, a3_main$att),
  low   = c(ci_gap_pre[1],  a1_main$ci[1], a2_main$ci[1], a3_main$ci[1]),
  high  = c(ci_gap_pre[2],  a1_main$ci[2], a2_main$ci[2], a3_main$ci[2]),
  Phase = c("Pre-treatment", "Phase 1", "Phase 2", "Phase 3")
)

tibble(
  Period        = c("Pre-treatment (2017–2018)", "Phase 1 (2019–Mar 2020)",
                    "Phase 2 (Mar–Dec 2020)",    "Phase 3 (Dec 2020–Dec 2023)"),
  `Gap / ATT`   = att_segments$att,
  `Lower 95% CI`= att_segments$low,
  `Upper 95% CI`= att_segments$high
) %>%
  kable(
    format = "pandoc", digits = 3, align = c("l", "r", "r", "r"),
    caption = "ATT estimates (armas_conUD) by phase. Pre-treatment: bootstrap; Phases 1–3: augsynth conformal inference."
  ) %>%
  kable_styling(position = "center", full_width = FALSE)


# --- Plot 5.1: Gap plot with phase-wise ATT and CI ribbons ---
ggplot(principalCSdelitos, aes(x = Fecha, y = diff)) +
  geom_line(aes(color = "Monthly gap"), linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_rect(data = att_segments,
            aes(xmin = start, xmax = end, ymin = low, ymax = high,
                fill = "95% CI"),
            inherit.aes = FALSE, alpha = 0.25) +
  geom_segment(data = att_segments,
               aes(x = start, xend = end, y = att, yend = att,
                   color = "Phase ATT"),
               inherit.aes = FALSE, linewidth = 1.2) +
  geom_vline(xintercept = as.Date("2019-01-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-03-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-12-01"), linetype = "dotted", linewidth = 1) +
  annotate("text", x = as.Date("2018-12-15"), y = -2.5,
           label = "Program starts",  angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = as.Date("2020-02-15"), y = -2.48,
           label = "Lockdown starts", angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = as.Date("2020-11-15"), y = -2.48,
           label = "Lockdown ends",   angle = 90, vjust = -0.5, size = 3) +
  scale_color_manual(name = NULL,
    values = c("Monthly gap" = "gray30", "Phase ATT" = "#449DD1")) +
  scale_fill_manual(name = NULL, values = c("95% CI" = "#B3CDE3")) +
  scale_x_date(date_breaks = "4 months", date_labels = "%m/%Y",
               expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Date", y = "Gap (Observed – Synthetic)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title  = element_text(face = "bold"),
        legend.position = "bottom")


# --- Mean pre-treatment level for effect-size benchmarking ---
delitosCS %>%
  filter(Entidad == "Ciudad de México",
         Fecha >= as.Date("2017-01-01"), Fecha <= as.Date("2018-12-31")) %>%
  summarise(mean_pre = mean(armas_conUD, na.rm = TRUE))


# --- Temporal placebos: re-run SC with fake treatment dates ---

# Helper: build a temporal placebo SC with the same predictor logic
build_temporal_placebo <- function(data, i_time, opt_start, opt_end,
                                   extra_predictors = NULL) {
  sc <- data %>%
    filter(!Entidad %in% donor_exclude) %>%
    synthetic_control(
      outcome           = armas_conUD,
      unit              = Entidad,
      time              = Fecha,
      i_unit            = "Ciudad de México",
      i_time            = i_time,
      generate_placebos = FALSE
    ) %>%
    generate_predictor(
      time_window = c(opt_start, opt_end),
      trespassing       = mean(`Allanamiento de morada`, na.rm = TRUE),
      grand_theft_auto  = mean(`Robo de coche de 4 ruedas Con violencia`, na.rm = TRUE),
      unemployment_rate = mean(unemp_rate, na.rm = TRUE),
      ipi               = mean(ipi, na.rm = TRUE),
      ipi_pct_change    = mean(ipi_pct_change_month, na.rm = TRUE),
      drug_dealing      = mean(Narcomenudeo, na.rm = TRUE),
      kidnapping        = mean(`Secuestro extorsivo`, na.rm = TRUE),
      weapon            = mean(`Con arma blanca`, na.rm = TRUE),
      tractors          = mean(`Robo de tractores Con violencia`, na.rm = TRUE),
      threats           = mean(Amenazas, na.rm = TRUE),
      family_violence   = mean(`Violencia familiar`, na.rm = TRUE),
      property_damage   = mean(`Daño a la propiedad`, na.rm = TRUE),
      bodily_harm       = mean(`Otros delitos que atentan contra la vida y la integridad corporal`,
                               na.rm = TRUE)
    ) %>%
    generate_predictor(time_window = as.Date("2017-01-01"), FC0117 = armas_conUD) %>%
    generate_predictor(time_window = as.Date("2017-05-01"), FC0517 = armas_conUD) %>%
    generate_predictor(time_window = as.Date("2017-09-01"), FC0917 = armas_conUD) %>%
    generate_predictor(time_window = as.Date("2018-01-01"), FC0118 = armas_conUD) %>%
    generate_predictor(time_window = as.Date("2018-05-01"), FC0518 = armas_conUD) %>%
    generate_predictor(time_window = as.Date("2018-09-01"), FC0918 = armas_conUD) %>%
    generate_predictor(time_window = c(as.Date("2017-01-01"), as.Date("2017-03-01")),
                       q1_2017 = mean(armas_conUD, na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2017-04-01"), as.Date("2017-06-01")),
                       q2_2017 = mean(armas_conUD, na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2017-07-01"), as.Date("2017-09-01")),
                       q3_2017 = mean(armas_conUD, na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2017-10-01"), as.Date("2017-12-01")),
                       q4_2017 = mean(armas_conUD, na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2018-01-01"), as.Date("2018-03-01")),
                       q1_2018 = mean(armas_conUD, na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2018-04-01"), as.Date("2018-06-01")),
                       q2_2018 = mean(armas_conUD, na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2018-07-01"), as.Date("2018-09-01")),
                       q3_2018 = mean(armas_conUD, na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2018-10-01"), as.Date("2018-12-01")),
                       q4_2018 = mean(armas_conUD, na.rm = TRUE)) %>%
    generate_predictor(c(as.Date("2017-01-01"), as.Date("2017-12-01")),
                       average2017 = mean(armas_conUD, na.rm = TRUE)) %>%
    generate_predictor(c(as.Date("2018-01-01"), as.Date("2018-12-01")),
                       average2018 = mean(armas_conUD, na.rm = TRUE)) %>%
    generate_predictor(c(as.Date("2019-01-01"), as.Date("2019-12-01")),
                       average2019 = mean(armas_conUD, na.rm = TRUE))

  if (!is.null(extra_predictors)) sc <- extra_predictors(sc)

  sc %>%
    generate_weights(
      optimization_window = c(opt_start, opt_end),
      margin_ipop = 0.02, sigf_ipop = 7, bound_ipop = 5
    ) %>%
    generate_control()
}

placebo_mar2020 <- build_temporal_placebo(
  delitosCS,
  i_time    = as.Date("2020-03-01"),
  opt_start = as.Date("2017-01-01"),
  opt_end   = as.Date("2020-02-01")
)

placebo_dic2020 <- delitosCS %>%
  filter(!Entidad %in% donor_exclude) %>%
  synthetic_control(
    outcome           = armas_conUD,
    unit              = Entidad,
    time              = Fecha,
    i_unit            = "Ciudad de México",
    i_time            = as.Date("2020-12-01"),
    generate_placebos = FALSE
  ) %>%
  generate_predictor(
    time_window = c(as.Date("2017-01-01"), as.Date("2020-11-01")),
    trespassing       = mean(`Allanamiento de morada`, na.rm = TRUE),
    grand_theft_auto  = mean(`Robo de coche de 4 ruedas Con violencia`, na.rm = TRUE),
    unemployment_rate = mean(unemp_rate, na.rm = TRUE),
    ipi               = mean(ipi, na.rm = TRUE),
    ipi_pct_change    = mean(ipi_pct_change_month, na.rm = TRUE),
    drug_dealing      = mean(Narcomenudeo, na.rm = TRUE),
    kidnapping        = mean(`Secuestro extorsivo`, na.rm = TRUE),
    weapon            = mean(`Con arma blanca`, na.rm = TRUE),
    tractors          = mean(`Robo de tractores Con violencia`, na.rm = TRUE),
    threats           = mean(Amenazas, na.rm = TRUE),
    family_violence   = mean(`Violencia familiar`, na.rm = TRUE),
    property_damage   = mean(`Daño a la propiedad`, na.rm = TRUE),
    bodily_harm       = mean(`Otros delitos que atentan contra la vida y la integridad corporal`,
                             na.rm = TRUE)
  ) %>%
  generate_predictor(time_window = as.Date("2017-01-01"), FC0117 = armas_conUD) %>%
  generate_predictor(time_window = as.Date("2017-05-01"), FC0517 = armas_conUD) %>%
  generate_predictor(time_window = as.Date("2017-09-01"), FC0917 = armas_conUD) %>%
  generate_predictor(time_window = as.Date("2018-01-01"), FC0118 = armas_conUD) %>%
  generate_predictor(time_window = as.Date("2018-05-01"), FC0518 = armas_conUD) %>%
  generate_predictor(time_window = as.Date("2018-09-01"), FC0918 = armas_conUD) %>%
  generate_predictor(time_window = c(as.Date("2017-01-01"), as.Date("2017-03-01")),
                     q1_2017 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2017-04-01"), as.Date("2017-06-01")),
                     q2_2017 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2017-07-01"), as.Date("2017-09-01")),
                     q3_2017 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2017-10-01"), as.Date("2017-12-01")),
                     q4_2017 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2018-01-01"), as.Date("2018-03-01")),
                     q1_2018 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2018-04-01"), as.Date("2018-06-01")),
                     q2_2018 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2018-07-01"), as.Date("2018-09-01")),
                     q3_2018 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2018-10-01"), as.Date("2018-12-01")),
                     q4_2018 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(c(as.Date("2017-01-01"), as.Date("2017-12-01")),
                     average2017 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(c(as.Date("2018-01-01"), as.Date("2018-12-01")),
                     average2018 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_predictor(c(as.Date("2019-01-01"), as.Date("2019-12-01")),
                     average2019 = mean(armas_conUD, na.rm = TRUE)) %>%
  # Dec-2020 placebo only: additional pre-year average through Nov 2020
  generate_predictor(c(as.Date("2020-01-01"), as.Date("2020-11-01")),
                     average2020 = mean(armas_conUD, na.rm = TRUE)) %>%
  generate_weights(
    optimization_window = c(as.Date("2017-01-01"), as.Date("2020-11-01")),
    margin_ipop = 0.02, sigf_ipop = 7, bound_ipop = 5
  ) %>%
  generate_control()

# Extract and plot temporal placebos
extract_synth_series <- function(sc_obj, delitosCS) {
  sc_obj %>%
    unnest(cols = c(".synthetic_control")) %>%
    filter(.id == "Ciudad de México", .type == "treated") %>%
    select(Fecha = time_unit, real_y, synth_y) %>%
    mutate(diff = real_y - synth_y) %>%
    left_join(
      delitosCS %>% filter(Entidad == "Ciudad de México") %>% select(Fecha, armas_sinUD),
      by = "Fecha"
    )
}

plot_obs_vs_synth <- function(df, y_annot = 1.3) {
  ggplot(df, aes(x = Fecha)) +
    geom_line(aes(y = real_y,  color = "Observed"),  linewidth = 1.2) +
    geom_line(aes(y = synth_y, color = "Synthetic"), linetype = "longdash", linewidth = 1.1) +
    geom_vline(xintercept = as.Date("2019-01-01"), linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = as.Date("2020-03-01"), linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = as.Date("2020-12-01"), linetype = "dotted", linewidth = 1) +
    annotate("text", x = as.Date("2018-12-15"), y = y_annot,
             label = "Disarmament program starts", angle = 90, vjust = -0.5, size = 3) +
    annotate("text", x = as.Date("2020-02-15"), y = y_annot * 0.62,
             label = "Lockdown starts", angle = 90, vjust = -0.5, size = 3) +
    annotate("text", x = as.Date("2020-11-15"), y = y_annot * 0.62,
             label = "Lockdown ends",   angle = 90, vjust = -0.5, size = 3) +
    scale_color_manual(values = c("Observed" = "gray30", "Synthetic" = "#449DD1")) +
    scale_x_date(date_breaks = "4 months", date_labels = "%m/%Y",
                 expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Date", y = "Rate per 100,000 inhabitants", color = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title  = element_text(face = "bold"),
          legend.position = "bottom")
}

principalCS_mar2020 <- extract_synth_series(placebo_mar2020, delitosCS)
principalCS_dec2020 <- extract_synth_series(placebo_dic2020, delitosCS)

plot_obs_vs_synth(principalCS_mar2020)   # Plot 5.2: Temporal placebo — Mar 2020
plot_obs_vs_synth(principalCS_dec2020)   # Plot 5.3: Temporal placebo — Dec 2020


## =============================================================================
## 6. Secondary outcomes: SC and Synthetic DiD by crime type
## =============================================================================

# Helper: build outcome panel, run SC, extract series, and plot
build_outcome_panel <- function(subtype_filter, modality_filter, var_name, negate = FALSE) {
  filter_fn <- if (negate) {
    function(d) filter(d, `Subtipo de delito` == subtype_filter, Modalidad != modality_filter)
  } else {
    function(d) filter(d, `Subtipo de delito` == subtype_filter, Modalidad == modality_filter)
  }

  filter_fn(delitos) %>%
    group_by(Entidad, Fecha) %>%
    summarise(valor = sum(valor, na.rm = TRUE), .groups = "drop") %>%
    right_join(delitosCS %>% select(Entidad, Fecha, Poblacion2020),
               by = c("Entidad", "Fecha")) %>%
    mutate(valor = coalesce(valor, 0),
           !!var_name := (valor * 100000) / Poblacion2020) %>%
    select(Entidad, Fecha, !!sym(var_name))
}

run_sc_outcome <- function(panel_full, outcome_var,
                           donor_excl = donor_exclude,
                           opt_margin = 0.02) {
  panel_full %>%
    filter(!Entidad %in% donor_excl) %>%
    synthetic_control(
      outcome           = !!sym(outcome_var),
      unit              = Entidad,
      time              = Fecha,
      i_unit            = "Ciudad de México",
      i_time            = as.Date("2019-01-01"),
      generate_placebos = TRUE
    ) %>%
    generate_predictor(
      time_window = c(as.Date("2017-01-01"), as.Date("2018-12-01")),
      trespassing       = mean(`Allanamiento de morada`, na.rm = TRUE),
      grand_theft_auto  = mean(`Robo de coche de 4 ruedas Con violencia`, na.rm = TRUE),
      unemployment_rate = mean(unemp_rate, na.rm = TRUE),
      ipi               = mean(ipi, na.rm = TRUE),
      ipi_pct_change    = mean(ipi_pct_change_month, na.rm = TRUE),
      drug_dealing      = mean(Narcomenudeo, na.rm = TRUE),
      kidnapping        = mean(`Secuestro extorsivo`, na.rm = TRUE),
      weapon            = mean(`Con arma blanca`, na.rm = TRUE),
      tractors          = mean(`Robo de tractores Con violencia`, na.rm = TRUE),
      threats           = mean(Amenazas, na.rm = TRUE),
      family_violence   = mean(`Violencia familiar`, na.rm = TRUE),
      property_damage   = mean(`Daño a la propiedad`, na.rm = TRUE),
      bodily_harm       = mean(`Otros delitos que atentan contra la vida y la integridad corporal`,
                               na.rm = TRUE)
    ) %>%
    generate_predictor(time_window = as.Date("2017-01-01"),
                       FC0117 = !!sym(outcome_var)) %>%
    generate_predictor(time_window = as.Date("2017-05-01"),
                       FC0517 = !!sym(outcome_var)) %>%
    generate_predictor(time_window = as.Date("2017-09-01"),
                       FC0917 = !!sym(outcome_var)) %>%
    generate_predictor(time_window = as.Date("2018-01-01"),
                       FC0118 = !!sym(outcome_var)) %>%
    generate_predictor(time_window = as.Date("2018-05-01"),
                       FC0518 = !!sym(outcome_var)) %>%
    generate_predictor(time_window = as.Date("2018-09-01"),
                       FC0918 = !!sym(outcome_var)) %>%
    generate_predictor(time_window = c(as.Date("2017-01-01"), as.Date("2017-03-01")),
                       q1_2017 = mean(!!sym(outcome_var), na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2017-04-01"), as.Date("2017-06-01")),
                       q2_2017 = mean(!!sym(outcome_var), na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2017-07-01"), as.Date("2017-09-01")),
                       q3_2017 = mean(!!sym(outcome_var), na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2017-10-01"), as.Date("2017-12-01")),
                       q4_2017 = mean(!!sym(outcome_var), na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2018-01-01"), as.Date("2018-03-01")),
                       q1_2018 = mean(!!sym(outcome_var), na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2018-04-01"), as.Date("2018-06-01")),
                       q2_2018 = mean(!!sym(outcome_var), na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2018-07-01"), as.Date("2018-09-01")),
                       q3_2018 = mean(!!sym(outcome_var), na.rm = TRUE)) %>%
    generate_predictor(time_window = c(as.Date("2018-10-01"), as.Date("2018-12-01")),
                       q4_2018 = mean(!!sym(outcome_var), na.rm = TRUE)) %>%
    generate_predictor(c(as.Date("2017-01-01"), as.Date("2017-12-01")),
                       average2017 = mean(!!sym(outcome_var), na.rm = TRUE)) %>%
    generate_predictor(c(as.Date("2018-01-01"), as.Date("2018-12-01")),
                       average2018 = mean(!!sym(outcome_var), na.rm = TRUE)) %>%
    generate_weights(
      optimization_window = c(as.Date("2017-01-01"), as.Date("2018-12-01")),
      margin_ipop = opt_margin, sigf_ipop = 7, bound_ipop = 5
    ) %>%
    generate_control()
}

run_synthdid <- function(panel_full, outcome_var, donor_excl = donor_exclude,
                         treat_date = as.Date("2019-01-01"),
                         post_end   = as.Date("2020-03-01"),
                         n_boot     = 500) {
  set.seed(123)
  treated_unit  <- "Ciudad de México"

  panel_sd <- panel_full %>%
    filter(!Entidad %in% donor_excl) %>%
    group_by(Entidad, Fecha) %>%
    summarise(y = mean(.data[[outcome_var]], na.rm = TRUE), .groups = "drop")

  time_points   <- sort(unique(panel_sd$Fecha))
  control_units <- setdiff(unique(panel_sd$Entidad), treated_unit)
  units_order   <- c(control_units, treated_unit)

  Y <- panel_sd %>%
    mutate(Entidad = factor(Entidad, levels = units_order)) %>%
    arrange(Entidad, Fecha) %>%
    pivot_wider(names_from = Fecha, values_from = y, values_fill = 0) %>%
    column_to_rownames("Entidad") %>%
    as.matrix()

  keep_dates <- time_points[time_points <= post_end]
  Y  <- Y[, as.character(keep_dates)]
  N0 <- length(control_units)
  T0 <- which(keep_dates == treat_date) - 1

  tau  <- synthdid_estimate(Y, N0 = N0, T0 = T0)
  se_p <- sqrt(as.numeric(vcov(tau, method = "placebo")))
  ci_p <- as.numeric(tau) + c(-1, 1) * 1.96 * se_p
  se_b <- tryCatch(
    sqrt(as.numeric(vcov(tau, method = "bootstrap", reps = n_boot, seed = 123))),
    error = function(e) NA_real_
  )
  ci_b <- if (is.na(se_b)) c(NA_real_, NA_real_) else
    as.numeric(tau) + c(-1, 1) * 1.96 * se_b

  sc_hat  <- sc_estimate (Y, N0, T0)
  did_hat <- did_estimate(Y, N0, T0)

  list(tau = tau, se_placebo = se_p, ci_placebo = ci_p,
       se_boot = se_b, ci_boot = ci_b,
       sc = sc_hat, did = did_hat, N0 = N0, T0 = T0, Y = Y,
       keep_dates = keep_dates)
}

# Generic SC + SynthDiD result printer
print_synthdid_results <- function(res, label) {
  cat(sprintf("\n=== Synthetic DiD (%s) — Post: 2019-01 to 2020-03 ===\n", label))
  cat(sprintf("ATT (post average): %0.4f\n", as.numeric(res$tau)))
  cat(sprintf("SE (placebo):       %0.4f | 95%% CI: (%0.4f, %0.4f)\n",
              res$se_placebo, res$ci_placebo[1], res$ci_placebo[2]))
  cat(sprintf("SE (bootstrap):     %s    | 95%% CI: (%s, %s)\n",
              ifelse(is.na(res$se_boot), "NA", sprintf("%0.4f", res$se_boot)),
              ifelse(is.na(res$se_boot), "NA", sprintf("%0.4f", res$ci_boot[1])),
              ifelse(is.na(res$se_boot), "NA", sprintf("%0.4f", res$ci_boot[2]))))
  cat(sprintf("SC ATT: %0.4f  |  DiD ATT: %0.4f\n",
              as.numeric(res$sc), as.numeric(res$did)))
  print(summary(res$tau))
}


# -----------------------------------------------------------------------------
# 6.1  Firearm homicides (homi_af_rate)
# -----------------------------------------------------------------------------
homi_af_panel <- build_outcome_panel("Homicidio doloso", "Con arma de fuego", "homi_af_rate") %>%
  { left_join(delitosCS, ., by = c("Entidad", "Fecha")) } %>%
  mutate(homi_af_rate = coalesce(homi_af_rate, 0))

cs_homi_af <- run_sc_outcome(homi_af_panel, "homi_af_rate")

panel_synth_homi_af <- cs_homi_af %>%
  unnest(cols = c(.synthetic_control)) %>%
  filter(.id == "Ciudad de México", .type == "treated") %>%
  select(Fecha = time_unit, Observed = real_y, Synthetic = synth_y) %>%
  mutate(Gap = Observed - Synthetic)

# Plot 6.1: Observed vs. Synthetic — firearm homicides
ggplot(panel_synth_homi_af, aes(x = Fecha)) +
  geom_line(aes(y = Observed,  color = "Observed"),  linewidth = 1.2) +
  geom_line(aes(y = Synthetic, color = "Synthetic"), linetype = "longdash", linewidth = 1.1) +
  geom_vline(xintercept = as.Date("2019-01-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-03-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-12-01"), linetype = "dotted", linewidth = 1) +
  annotate("text", x = as.Date("2018-12-15"),
           y = max(panel_synth_homi_af$Observed, na.rm = TRUE) * 0.35,
           label = "Disarmament program starts", angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = as.Date("2020-02-15"),
           y = max(panel_synth_homi_af$Observed, na.rm = TRUE) * 0.35,
           label = "Lockdown starts", angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = as.Date("2020-11-15"),
           y = max(panel_synth_homi_af$Observed, na.rm = TRUE) * 0.35,
           label = "Lockdown ends",   angle = 90, vjust = -0.5, size = 3) +
  scale_color_manual(values = c("Observed" = "gray30", "Synthetic" = "#449DD1")) +
  scale_x_date(date_breaks = "4 months", date_labels = "%m/%Y") +
  labs(x = "Date", y = "Rate per 100,000 inhabitants", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title  = element_text(face = "bold"),
        legend.position = "bottom")

# Synthetic DiD — firearm homicides
sdid_homi_af <- run_synthdid(homi_af_panel, "homi_af_rate")
print_synthdid_results(sdid_homi_af, "homi_af_rate")

print(synthdid_plot(
  list("SDID" = sdid_homi_af$tau, "SC" = sdid_homi_af$sc, "DiD" = sdid_homi_af$did),
  facet = NULL, overlay = 1, ci.alpha = 0.01, line.width = 0.8,
  treated.name = "CDMX", control.name = "Synthetic CDMX"
))


# -----------------------------------------------------------------------------
# 6.2  Firearm injuries (lesi_dolosa_af)
# -----------------------------------------------------------------------------
lesi_af_panel <- build_outcome_panel("Lesiones dolosas", "Con arma de fuego", "lesi_dolosa_af") %>%
  { left_join(delitosCS, ., by = c("Entidad", "Fecha")) } %>%
  mutate(lesi_dolosa_af = coalesce(lesi_dolosa_af, 0))

cs_lesi_af <- run_sc_outcome(lesi_af_panel, "lesi_dolosa_af")

panel_synth_lesi_af <- cs_lesi_af %>%
  unnest(cols = c(.synthetic_control)) %>%
  filter(.id == "Ciudad de México", .type == "treated") %>%
  select(Fecha = time_unit, Observed = real_y, Synthetic = synth_y) %>%
  mutate(Gap = Observed - Synthetic)

# Plot 6.2: Observed vs. Synthetic — firearm injuries
ggplot(panel_synth_lesi_af, aes(x = Fecha)) +
  geom_line(aes(y = Observed,  color = "Observed"),  linewidth = 1.2) +
  geom_line(aes(y = Synthetic, color = "Synthetic"), linetype = "longdash", linewidth = 1.1) +
  geom_vline(xintercept = as.Date("2019-01-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-03-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-12-01"), linetype = "dotted", linewidth = 1) +
  annotate("text", x = as.Date("2018-12-15"),
           y = max(panel_synth_lesi_af$Observed, na.rm = TRUE) * 0.95,
           label = "Disarmament program starts", angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = as.Date("2020-02-15"),
           y = max(panel_synth_lesi_af$Observed, na.rm = TRUE) * 0.85,
           label = "Lockdown starts", angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = as.Date("2020-11-15"),
           y = max(panel_synth_lesi_af$Observed, na.rm = TRUE) * 0.85,
           label = "Lockdown ends",   angle = 90, vjust = -0.5, size = 3) +
  scale_color_manual(values = c("Observed" = "gray30", "Synthetic" = "#449DD1")) +
  scale_x_date(date_breaks = "4 months", date_labels = "%m/%Y") +
  labs(x = "Date", y = "Rate per 100,000 inhabitants", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title  = element_text(face = "bold"),
        legend.position = "bottom")

sdid_lesi_af <- run_synthdid(lesi_af_panel, "lesi_dolosa_af")
print_synthdid_results(sdid_lesi_af, "lesi_dolosa_af")

print(synthdid_plot(
  list("SDID" = sdid_lesi_af$tau, "SC" = sdid_lesi_af$sc, "DiD" = sdid_lesi_af$did),
  facet = NULL, overlay = 1, ci.alpha = 0.01, line.width = 1,
  treated.name = "CDMX", control.name = "Synthetic CDMX"
))


# -----------------------------------------------------------------------------
# 6.3  Non-firearm homicides (homi_nofire_rate)
#      Additional step: drop top 15% most volatile donors in pre-period.
# -----------------------------------------------------------------------------
homi_nofire_raw <- build_outcome_panel(
  "Homicidio doloso", "Con arma de fuego", "homi_nofire_rate", negate = TRUE
)

homi_nofire_panel <- delitosCS %>%
  left_join(homi_nofire_raw, by = c("Entidad", "Fecha")) %>%
  mutate(homi_nofire_rate = coalesce(homi_nofire_rate, 0))

# Identify high-volatility donors and extend exclusion list
volatiles_hn <- homi_nofire_panel %>%
  filter(Fecha >= as.Date("2017-01-01"), Fecha <= as.Date("2018-12-01")) %>%
  group_by(Entidad) %>%
  summarise(sd_pre = sd(homi_nofire_rate, na.rm = TRUE), .groups = "drop") %>%
  mutate(rank = percent_rank(sd_pre)) %>%
  filter(rank >= 0.85) %>%
  pull(Entidad)

donor_exclude_hn <- union(donor_exclude, volatiles_hn)

# NOTE: This model uses a bespoke predictor set designed for non-firearm
# homicides (modality proxies, coercion variables, firearm availability
# controls) and bimonthly snapshots — intentionally different from the
# standard run_sc_outcome() template used for the other outcomes.
cs_homi_nofire <- homi_nofire_panel %>%
  filter(!Entidad %in% donor_exclude_hn) %>%
  synthetic_control(
    outcome           = homi_nofire_rate,
    unit              = Entidad,
    time              = Fecha,
    i_unit            = "Ciudad de México",
    i_time            = as.Date("2019-01-01"),
    generate_placebos = TRUE
  ) %>%
  # --- Baseline predictors: theory-driven for non-firearm homicides ---
  generate_predictor(
    time_window = c(as.Date("2017-01-01"), as.Date("2018-12-01")),
    with_knife         = mean(`Con arma blanca`, na.rm = TRUE),
    with_other_element = mean(`Con otro elemento`, na.rm = TRUE),
    with_violence      = mean(`Con violencia`, na.rm = TRUE),
    without_violence   = mean(`Sin violencia`, na.rm = TRUE),
    threats            = mean(Amenazas, na.rm = TRUE),
    extortion          = mean(`Extorsión`, na.rm = TRUE),
    family_violence    = mean(`Violencia familiar`, na.rm = TRUE),
    bodily_harm        = mean(`Otros delitos que atentan contra la vida y la integridad corporal`,
                              na.rm = TRUE),
    drug_dealing       = mean(Narcomenudeo, na.rm = TRUE),
    unemployment_rate  = mean(unemp_rate, na.rm = TRUE),
    ipi_level          = mean(ipi, na.rm = TRUE),
    ipi_change         = mean(ipi_pct_change_month, na.rm = TRUE),
    armas_conUD_mean   = mean(armas_conUD, na.rm = TRUE),
    armas_sinUD_mean   = mean(armas_sinUD, na.rm = TRUE)
  ) %>%
  # --- Bimonthly snapshots (every 2 months, 2017–2018) ---
  generate_predictor(time_window = as.Date("2017-01-01"), s_2017_01 = homi_nofire_rate) %>%
  generate_predictor(time_window = as.Date("2017-03-01"), s_2017_03 = homi_nofire_rate) %>%
  generate_predictor(time_window = as.Date("2017-05-01"), s_2017_05 = homi_nofire_rate) %>%
  generate_predictor(time_window = as.Date("2017-07-01"), s_2017_07 = homi_nofire_rate) %>%
  generate_predictor(time_window = as.Date("2017-09-01"), s_2017_09 = homi_nofire_rate) %>%
  generate_predictor(time_window = as.Date("2017-11-01"), s_2017_11 = homi_nofire_rate) %>%
  generate_predictor(time_window = as.Date("2018-01-01"), s_2018_01 = homi_nofire_rate) %>%
  generate_predictor(time_window = as.Date("2018-03-01"), s_2018_03 = homi_nofire_rate) %>%
  generate_predictor(time_window = as.Date("2018-05-01"), s_2018_05 = homi_nofire_rate) %>%
  generate_predictor(time_window = as.Date("2018-07-01"), s_2018_07 = homi_nofire_rate) %>%
  generate_predictor(time_window = as.Date("2018-09-01"), s_2018_09 = homi_nofire_rate) %>%
  generate_predictor(time_window = as.Date("2018-11-01"), s_2018_11 = homi_nofire_rate) %>%
  # --- Quarterly averages (2017–2018) ---
  generate_predictor(time_window = c(as.Date("2017-01-01"), as.Date("2017-03-01")),
                     q1_2017 = mean(homi_nofire_rate, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2017-04-01"), as.Date("2017-06-01")),
                     q2_2017 = mean(homi_nofire_rate, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2017-07-01"), as.Date("2017-09-01")),
                     q3_2017 = mean(homi_nofire_rate, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2017-10-01"), as.Date("2017-12-01")),
                     q4_2017 = mean(homi_nofire_rate, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2018-01-01"), as.Date("2018-03-01")),
                     q1_2018 = mean(homi_nofire_rate, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2018-04-01"), as.Date("2018-06-01")),
                     q2_2018 = mean(homi_nofire_rate, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2018-07-01"), as.Date("2018-09-01")),
                     q3_2018 = mean(homi_nofire_rate, na.rm = TRUE)) %>%
  generate_predictor(time_window = c(as.Date("2018-10-01"), as.Date("2018-12-01")),
                     q4_2018 = mean(homi_nofire_rate, na.rm = TRUE)) %>%
  # --- Annual level anchors ---
  generate_predictor(c(as.Date("2017-01-01"), as.Date("2017-12-01")),
                     average2017 = mean(homi_nofire_rate, na.rm = TRUE)) %>%
  generate_predictor(c(as.Date("2018-01-01"), as.Date("2018-12-01")),
                     average2018 = mean(homi_nofire_rate, na.rm = TRUE)) %>%
  # --- Weight optimization (tighter margin for better pre-fit) ---
  generate_weights(
    optimization_window = c(as.Date("2017-01-01"), as.Date("2018-12-01")),
    margin_ipop = 0.01, sigf_ipop = 7, bound_ipop = 5
  ) %>%
  generate_control()

panel_synth_homi_nofire <- cs_homi_nofire %>%
  unnest(cols = c(.synthetic_control)) %>%
  filter(.id == "Ciudad de México", .type == "treated") %>%
  select(Fecha = time_unit, Observed = real_y, Synthetic = synth_y) %>%
  mutate(Gap = Observed - Synthetic)

# Plot 6.3
ggplot(panel_synth_homi_nofire, aes(x = Fecha)) +
  geom_line(aes(y = Observed,  color = "Observed"),  linewidth = 1.2) +
  geom_line(aes(y = Synthetic, color = "Synthetic"), linetype = "longdash", linewidth = 1.1) +
  geom_vline(xintercept = as.Date("2019-01-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-03-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-12-01"), linetype = "dotted", linewidth = 1) +
  annotate("text", x = as.Date("2018-12-15"),
           y = max(panel_synth_homi_nofire$Observed, na.rm = TRUE) * 0.95,
           label = "Disarmament program starts", angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = as.Date("2020-02-15"),
           y = max(panel_synth_homi_nofire$Observed, na.rm = TRUE) * 0.85,
           label = "Lockdown starts", angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = as.Date("2020-11-15"),
           y = max(panel_synth_homi_nofire$Observed, na.rm = TRUE) * 0.85,
           label = "Lockdown ends",   angle = 90, vjust = -0.5, size = 3) +
  scale_color_manual(values = c("Observed" = "gray30", "Synthetic" = "#449DD1")) +
  scale_x_date(date_breaks = "4 months", date_labels = "%m/%Y") +
  labs(x = "Date", y = "Rate per 100,000 inhabitants", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title  = element_text(face = "bold"),
        legend.position = "bottom")

sdid_homi_nofire <- run_synthdid(homi_nofire_panel, "homi_nofire_rate",
                                 donor_excl = donor_exclude_hn)
print_synthdid_results(sdid_homi_nofire, "homi_nofire_rate")

print(synthdid_plot(
  list("SDID" = sdid_homi_nofire$tau, "SC" = sdid_homi_nofire$sc,
       "DiD" = sdid_homi_nofire$did),
  facet = NULL, overlay = 1, ci.alpha = 0.01, line.width = 1,
  treated.name = "CDMX", control.name = "Synthetic CDMX"
))


# -----------------------------------------------------------------------------
# 6.4  Non-firearm injuries (lesi_nofire_rate)
# -----------------------------------------------------------------------------
lesi_nofire_raw <- build_outcome_panel(
  "Lesiones dolosas", "Con arma de fuego", "lesi_nofire_rate", negate = TRUE
)

lesi_nofire_panel <- delitosCS %>%
  left_join(lesi_nofire_raw, by = c("Entidad", "Fecha")) %>%
  mutate(lesi_nofire_rate = coalesce(lesi_nofire_rate, 0))

cs_lesi_nofire <- run_sc_outcome(lesi_nofire_panel, "lesi_nofire_rate")

panel_synth_lesi_nofire <- cs_lesi_nofire %>%
  unnest(cols = c(.synthetic_control)) %>%
  filter(.id == "Ciudad de México", .type == "treated") %>%
  select(Fecha = time_unit, Observed = real_y, Synthetic = synth_y) %>%
  mutate(Gap = Observed - Synthetic)

# Plot 6.4
ggplot(panel_synth_lesi_nofire, aes(x = Fecha)) +
  geom_line(aes(y = Observed,  color = "Observed"),  linewidth = 1.2) +
  geom_line(aes(y = Synthetic, color = "Synthetic"), linetype = "longdash", linewidth = 1.1) +
  geom_vline(xintercept = as.Date("2019-01-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-03-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-12-01"), linetype = "dotted", linewidth = 1) +
  annotate("text", x = as.Date("2018-12-15"),
           y = max(panel_synth_lesi_nofire$Observed, na.rm = TRUE) * 0.95,
           label = "Disarmament program starts", angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = as.Date("2020-02-15"),
           y = max(panel_synth_lesi_nofire$Observed, na.rm = TRUE) * 0.85,
           label = "Lockdown starts", angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = as.Date("2020-11-15"),
           y = max(panel_synth_lesi_nofire$Observed, na.rm = TRUE) * 0.85,
           label = "Lockdown ends",   angle = 90, vjust = -0.5, size = 3) +
  scale_color_manual(values = c("Observed" = "gray30", "Synthetic" = "#449DD1")) +
  scale_x_date(date_breaks = "4 months", date_labels = "%m/%Y") +
  labs(x = "Date", y = "Rate per 100,000 inhabitants", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title  = element_text(face = "bold"),
        legend.position = "bottom")

sdid_lesi_nofire <- run_synthdid(lesi_nofire_panel, "lesi_nofire_rate")
print_synthdid_results(sdid_lesi_nofire, "lesi_nofire_rate")

print(synthdid_plot(
  list("SDID" = sdid_lesi_nofire$tau, "SC" = sdid_lesi_nofire$sc,
       "DiD" = sdid_lesi_nofire$did),
  facet = NULL, overlay = 1, ci.alpha = 0.01, line.width = 0.8,
  treated.name = "CDMX", control.name = "Synthetic CDMX"
))


## =============================================================================
## 7. Phase-wise ATT for all secondary outcomes
## =============================================================================

# Helper: pre-treatment bootstrap gap using equal-weight naive synthetic
pre_gap_bootstrap <- function(panel, outcome_var, n_boot = 1000, seed = 123) {
  wide <- panel %>%
    filter(Fecha >= as.Date("2017-01-01"), Fecha <= as.Date("2018-12-01")) %>%
    select(Entidad, Fecha, !!sym(outcome_var)) %>%
    pivot_wider(names_from = Entidad, values_from = !!sym(outcome_var)) %>%
    arrange(Fecha)

  treated_y <- wide$`Ciudad de México`
  ctrl_mat  <- wide %>% select(-Fecha, -`Ciudad de México`) %>% as.matrix()
  w         <- rep(1 / ncol(ctrl_mat), ncol(ctrl_mat))
  synth_y   <- ctrl_mat %*% w
  gap_pre   <- treated_y - synth_y

  set.seed(seed)
  boots <- replicate(n_boot, mean(sample(gap_pre, replace = TRUE), na.rm = TRUE))
  list(
    avg = mean(gap_pre, na.rm = TRUE),
    ci  = quantile(boots, probs = c(0.025, 0.975), na.rm = TRUE)
  )
}

# Compute all phases for each outcome variable
outcomes_list <- list(
  list(panel = homi_af_panel,     var = "homi_af_rate",     label = "Homicide with firearm"),
  list(panel = homi_nofire_panel, var = "homi_nofire_rate", label = "Homicide without firearm"),
  list(panel = lesi_af_panel,     var = "lesi_dolosa_af",   label = "Injury with firearm"),
  list(panel = lesi_nofire_panel, var = "lesi_nofire_rate", label = "Injury without firearm")
)

phases <- list(
  list(start = as.Date("2017-01-01"), end = as.Date("2020-03-01"),
       treat = as.Date("2019-01-01"), code = 201901L),
  list(start = as.Date("2019-01-01"), end = as.Date("2020-12-01"),
       treat = as.Date("2020-03-01"), code = 202003L),
  list(start = as.Date("2020-03-01"), end = as.Date("2023-12-01"),
       treat = as.Date("2020-12-01"), code = 202012L)
)

att_all <- purrr::map_dfr(outcomes_list, function(o) {
  pre  <- pre_gap_bootstrap(o$panel, o$var)
  rows <- list(tibble(Outcome = o$label, Phase = "Pre",
                      ATT = pre$avg, CI_Low = pre$ci[1], CI_High = pre$ci[2]))

  for (ph in phases) {
    res <- build_augsynth_phase(o$panel, o$var, ph$start, ph$end, ph$treat)
    rows <- c(rows, list(tibble(
      Outcome = o$label,
      Phase   = as.character(ph$code),
      ATT     = res$att, CI_Low = res$ci[1], CI_High = res$ci[2]
    )))
  }
  bind_rows(rows)
})

att_all %>%
  mutate(Phase = recode(Phase,
    "Pre"    = "Pre-treatment (2017–2018)",
    "201901" = "Phase 1 (2019–Mar 2020)",
    "202003" = "Phase 2 (Mar–Dec 2020)",
    "202012" = "Phase 3 (Dec 2020–Dec 2023)"
  ),
  Phase = factor(Phase, levels = c("Pre-treatment (2017–2018)", "Phase 1 (2019–Mar 2020)",
                                    "Phase 2 (Mar–Dec 2020)", "Phase 3 (Dec 2020–Dec 2023)"))
  ) %>%
  arrange(Outcome, Phase) %>%
  kable(
    format = "pandoc", digits = 3, align = c("l", "c", "r", "r", "r"),
    caption = "ATT estimates and 95% CIs for all outcomes across four intervention phases."
  ) %>%
  kable_styling(position = "center", full_width = FALSE, font_size = 10)


## =============================================================================
## 8. Firearm buyback analyses
## =============================================================================

# State name standardization dictionary (abbreviations → official names)
estado_dict_es <- c(
  "AGS." = "Aguascalientes", "Aguascalientes" = "Aguascalientes",
  "B.C." = "Baja California", "Baja California" = "Baja California",
  "B.C.S." = "Baja California Sur", "Baja California Sur" = "Baja California Sur",
  "CAMP." = "Campeche", "Campeche" = "Campeche",
  "CD. MÉX." = "Ciudad de México", "D.F." = "Ciudad de México",
  "Cd. de México" = "Ciudad de México",
  "CHIH." = "Chihuahua", "Chihuahua" = "Chihuahua",
  "CHIS." = "Chiapas",   "Chiapas"   = "Chiapas",
  "COAH." = "Coahuila",  "Coahuila"  = "Coahuila",
  "COL."  = "Colima",    "Colima"    = "Colima",
  "DGO."  = "Durango",   "Durango"   = "Durango",
  "GTO."  = "Guanajuato","Guanajuato"= "Guanajuato",
  "GRO."  = "Guerrero",  "Guerrero"  = "Guerrero",
  "HGO."  = "Hidalgo",   "Hidalgo"   = "Hidalgo",
  "JAL."  = "Jalisco",   "Jalisco"   = "Jalisco",
  "MEX."  = "México", "MÉX." = "México", "México" = "México",
  "MICH." = "Michoacán", "Michoacán" = "Michoacán",
  "MOR."  = "Morelos",   "Morelos"   = "Morelos",
  "N.L."  = "Nuevo León","Nvo. León" = "Nuevo León","Nuevo León" = "Nuevo León",
  "NAY."  = "Nayarit",   "Nayarit"   = "Nayarit",
  "OAX."  = "Oaxaca",    "Oaxaca"    = "Oaxaca",
  "PUE."  = "Puebla",    "Puebla"    = "Puebla",
  "QRO."  = "Querétaro", "Querétaro" = "Querétaro",
  "Q. ROO"= "Quintana Roo","Q. ROO." = "Quintana Roo","Quintana Roo" = "Quintana Roo",
  "S.L.P."= "San Luis Potosí","San Luis Potosí" = "San Luis Potosí",
  "SIN."  = "Sinaloa",   "Sinaloa"   = "Sinaloa",
  "SON."  = "Sonora",    "Sonora"    = "Sonora",
  "TAB."  = "Tabasco",   "Tabasco"   = "Tabasco",
  "TAMPS."= "Tamaulipas","Tamaulipas"= "Tamaulipas",
  "TLAX." = "Tlaxcala",  "Tlaxcala"  = "Tlaxcala",
  "VER."  = "Veracruz",  "Veracruz"  = "Veracruz",
  "YUC."  = "Yucatán",   "Yucatán"   = "Yucatán",
  "ZAC."  = "Zacatecas", "Zacatecas" = "Zacatecas"
)

# Load raw buyback data
canjes_raw <- read_csv(
  file.path(DATA_DIR, "ANEXO FOLIO 330026424002055.csv")
)

canjes <- canjes_raw %>%
  mutate(ESTADO = estado_dict_es[ESTADO])

# Load municipality catalogue
catalogo <- read_csv(
  file.path(DATA_DIR, "AGEEML_20249161546763.csv"),
  locale = locale(encoding = "ISO-8859-1")
) %>%
  mutate(NOM_ENT = recode(NOM_ENT,
    "Coahuila de Zaragoza"            = "Coahuila",
    "Michoacán de Ocampo"             = "Michoacán",
    "Veracruz de Ignacio de la Llave" = "Veracruz"
  ))

# Validate state name alignment
stopifnot(length(setdiff(unique(canjes$ESTADO), unique(catalogo$NOM_ENT))) == 0)

# Monthly panel by state
canjes <- canjes %>%
  mutate(`FECHA EVENTO` = as.Date(`FECHA EVENTO`, format = "%d/%m/%Y"))

canjes_por_estado_mes <- canjes %>%
  mutate(mes = floor_date(`FECHA EVENTO`, unit = "month")) %>%
  group_by(ESTADO, mes) %>%
  summarise(TOTAL = sum(TOTAL, na.rm = TRUE),
            CORTA = sum(CORTA, na.rm = TRUE),
            LARGA = sum(LARGA, na.rm = TRUE),
            .groups = "drop") %>%
  rename(Entidad = ESTADO) %>%
  left_join(poblacion, by = "Entidad") %>%
  mutate(
    tasa_total = (TOTAL / Poblacion2020) * 100000,
    tasa_corta = (CORTA / Poblacion2020) * 100000,
    tasa_larga = (LARGA / Poblacion2020) * 100000
  )

# Merge with main panel
panel_buybacks_total <- delitosCS %>%
  left_join(
    canjes_por_estado_mes %>% select(Entidad, Fecha = mes,
                                     tasa_total, tasa_corta, tasa_larga),
    by = c("Entidad", "Fecha")
  ) %>%
  mutate(across(c(tasa_total, tasa_corta, tasa_larga), ~ coalesce(., 0)))


# -----------------------------------------------------------------------------
# 8.1  SC — total buybacks (tasa_total), peak-anchored predictors
# -----------------------------------------------------------------------------
main_CS_buybacks <- panel_buybacks_total %>%
  filter(!Entidad %in% donor_exclude | Entidad == "Nayarit") %>%
  synthetic_control(
    outcome           = tasa_total,
    unit              = Entidad,
    time              = Fecha,
    i_unit            = "Ciudad de México",
    i_time            = as.Date("2019-01-01"),
    generate_placebos = TRUE
  ) %>%
  generate_predictor(
    time_window = c(as.Date("2017-01-01"), as.Date("2018-12-01")),
    trespassing       = mean(`Allanamiento de morada`, na.rm = TRUE),
    grand_theft_auto  = mean(`Robo de coche de 4 ruedas Con violencia`, na.rm = TRUE),
    unemployment_rate = mean(unemp_rate, na.rm = TRUE),
    ipi               = mean(ipi, na.rm = TRUE),
    ipi_pct_change    = mean(ipi_pct_change_month, na.rm = TRUE),
    drug_dealing      = mean(Narcomenudeo, na.rm = TRUE),
    threats           = mean(Amenazas, na.rm = TRUE),
    firearm           = mean(`Con arma de fuego`, na.rm = TRUE),
    knife             = mean(`Con arma blanca`, na.rm = TRUE),
    family_violence   = mean(`Violencia familiar`, na.rm = TRUE)
  ) %>%
  generate_predictor(time_window = as.Date("2017-05-01"),
                     snapshot_2017_05 = tasa_total) %>%
  generate_predictor(time_window = as.Date("2018-01-01"),
                     snapshot_2018_01 = tasa_total) %>%
  generate_predictor(time_window = as.Date("2018-08-01"),
                     snapshot_2018_08 = tasa_total) %>%
  generate_predictor(c(as.Date("2017-01-01"), as.Date("2017-12-01")),
                     average2017 = mean(tasa_total, na.rm = TRUE)) %>%
  generate_predictor(c(as.Date("2018-01-01"), as.Date("2018-12-01")),
                     average2018 = mean(tasa_total, na.rm = TRUE)) %>%
  generate_weights(
    optimization_window = c(as.Date("2017-01-01"), as.Date("2018-12-01")),
    margin_ipop = 0.01, sigf_ipop = 7, bound_ipop = 5
  ) %>%
  generate_control()

principalCS_buybacks <- main_CS_buybacks %>%
  unnest(cols = c(.synthetic_control)) %>%
  filter(.id == "Ciudad de México", .type == "treated") %>%
  select(Fecha = time_unit, real_y, synth_y) %>%
  mutate(diff = real_y - synth_y)

# Plot 8.1: Observed vs. Synthetic — total buybacks
max_y_bb <- max(pmax(principalCS_buybacks$real_y, principalCS_buybacks$synth_y), na.rm = TRUE)

ggplot(principalCS_buybacks, aes(x = Fecha)) +
  geom_line(aes(y = real_y,  color = "Observed"),  linewidth = 0.9) +
  geom_line(aes(y = synth_y, color = "Synthetic"), linewidth = 0.9, linetype = "dashed") +
  geom_vline(xintercept = as.Date("2019-01-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-03-01"), linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = as.Date("2020-12-01"), linetype = "dotted", linewidth = 1) +
  annotate("text", x = as.Date("2019-01-01"), y = max_y_bb * 0.80,
           label = "Program starts",  angle = 90, vjust = -0.5, size = 3.5) +
  annotate("text", x = as.Date("2020-03-01"), y = max_y_bb * 0.80,
           label = "Lockdown starts", angle = 90, vjust = -0.5, size = 3.5) +
  annotate("text", x = as.Date("2020-12-01"), y = max_y_bb * 0.80,
           label = "Lockdown ends",   angle = 90, vjust = -0.5, size = 3.5) +
  scale_color_manual(values = c("Observed" = "#08306B", "Synthetic" = "#6BAED6"), name = NULL) +
  scale_x_date(date_breaks = "3 months", date_labels = "%m/%Y") +
  labs(x = "Date", y = "Firearm Buyback Rate (per 100,000)") +
  theme_minimal(base_size = 12) +
  theme(axis.title   = element_text(face = "bold", size = 14),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


# --- Pre-treatment bootstrap gap (buybacks) ---
weights_bb <- main_CS_buybacks %>%
  grab_unit_weights() %>%
  filter(unit != "Ciudad de México")

panel_pre_bb <- panel_buybacks_total %>%
  filter(Fecha >= as.Date("2017-01-01"), Fecha <= as.Date("2018-12-01")) %>%
  select(Entidad, Fecha, tasa_total) %>%
  pivot_wider(names_from = Entidad, values_from = tasa_total) %>%
  arrange(Fecha)

treated_y_bb  <- panel_pre_bb$`Ciudad de México`
control_mat_bb <- panel_pre_bb %>% select(all_of(weights_bb$unit)) %>% as.matrix()
synth_y_bb     <- control_mat_bb %*% weights_bb$weight
gap_pre_bb     <- treated_y_bb - synth_y_bb

set.seed(777)
boots_bb    <- replicate(1000, mean(sample(gap_pre_bb, replace = TRUE), na.rm = TRUE))
avg_gap_bb  <- mean(gap_pre_bb, na.rm = TRUE)
ci_gap_bb   <- quantile(boots_bb, probs = c(0.025, 0.975))

# --- Phase ATT (buybacks) using asc_phase wrapper ---
a1_bb <- build_augsynth_phase(panel_buybacks_total, "tasa_total",
  as.Date("2017-01-01"), as.Date("2020-03-01"), as.Date("2019-01-01"))
a2_bb <- build_augsynth_phase(panel_buybacks_total, "tasa_total",
  as.Date("2019-01-01"), as.Date("2020-12-01"), as.Date("2020-03-01"))
a3_bb <- build_augsynth_phase(panel_buybacks_total, "tasa_total",
  as.Date("2020-03-01"), as.Date("2023-12-01"), as.Date("2020-12-01"))

# Table 8.1: ATT summary — buybacks
tibble(
  Period        = c("Pre-treatment (2017–2018)", "Phase 1 (2019–Mar 2020)",
                    "Phase 2 (Mar–Dec 2020)",    "Phase 3 (Dec 2020–Dec 2023)"),
  `Gap / ATT`   = c(avg_gap_bb,   a1_bb$att, a2_bb$att, a3_bb$att),
  `Lower 95% CI`= c(ci_gap_bb[1], a1_bb$ci[1], a2_bb$ci[1], a3_bb$ci[1]),
  `Upper 95% CI`= c(ci_gap_bb[2], a1_bb$ci[2], a2_bb$ci[2], a3_bb$ci[2])
) %>%
  kable(
    format = "pandoc", digits = 3, align = c("l", "r", "r", "r"),
    caption = "ATT estimates for firearm buybacks (tasa_total) by phase."
  ) %>%
  kable_styling(position = "center", full_width = FALSE)

# --- Synthetic DiD — buybacks ---
sdid_buybacks <- run_synthdid(panel_buybacks_total, "tasa_total")
print_synthdid_results(sdid_buybacks, "tasa_total")

print(synthdid_plot(
  list("SDID" = sdid_buybacks$tau, "SC" = sdid_buybacks$sc, "DiD" = sdid_buybacks$did),
  facet = NULL, overlay = 1, ci.alpha = 0.01, line.width = 0.8,
  treated.name = "CDMX", control.name = "Synthetic CDMX"
))

# --- Plot 8.2: Buybacks by gun type (stacked bar — Mexico City only) ---
cdmx_lc <- panel_buybacks_total %>%
  filter(Entidad == "Ciudad de México") %>%
  select(Fecha, tasa_larga, tasa_corta) %>%
  mutate(across(c(tasa_larga, tasa_corta), ~ replace_na(., 0))) %>%
  pivot_longer(cols = c(tasa_larga, tasa_corta), names_to = "type", values_to = "rate") %>%
  mutate(
    type = recode(type, "tasa_corta" = "Handguns (short)", "tasa_larga" = "Long guns"),
    type = factor(type, levels = c("Handguns (short)", "Long guns"))
  )

ggplot(cdmx_lc, aes(x = Fecha, y = rate, fill = type)) +
  geom_col(position = "stack", width = 25) +
  geom_vline(xintercept = as.Date("2019-01-01"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-03-01"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-12-01"), linetype = "dashed") +
  annotate("text", x = as.Date("2019-01-01"),
           y = max(cdmx_lc$rate, na.rm = TRUE),
           label = "Program starts", angle = 90, vjust = -0.5, hjust = 1, size = 3) +
  annotate("text", x = as.Date("2020-03-01"),
           y = max(cdmx_lc$rate, na.rm = TRUE),
           label = "Lockdown starts", angle = 90, vjust = -0.5, hjust = 1, size = 3) +
  annotate("text", x = as.Date("2020-12-01"),
           y = max(cdmx_lc$rate, na.rm = TRUE),
           label = "Lockdown ends", angle = 90, vjust = -0.5, hjust = 1, size = 3) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_date(date_breaks = "3 months", date_labels = "%m/%Y",
               expand = expansion(mult = c(0.01, 0.01))) +
  scale_fill_manual(values = c("Handguns (short)" = "#08306B", "Long guns" = "#6BAED6")) +
  labs(x = "Date", y = "Firearm Buyback Rate (per 100,000)", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title  = element_text(face = "bold"),
        legend.position = "bottom")


## =============================================================================
## 9. Firearm seizures — two-way FE staggered DiD
##    Source: CNSPF (Censo Nacional de Seguridad Pública Federal), INEGI
##    NOTE: Files are labelled with year N but contain data for year N-1.
## =============================================================================

# Load annual seizure files
d2017 <- read_csv(file.path(DATA_DIR, "Decomisos/2018/conjunto_de_datos/arms_ent_cnspf2018.csv"))
d2018 <- read_csv(file.path(DATA_DIR, "Decomisos/2019/conjunto_de_datos/arms_ent_cnspf2019.csv"))
d2019 <- read_csv(file.path(DATA_DIR, "Decomisos/2020/conjunto_de_datos/m2p1_23_cnspf2020.csv"))
d2020 <- read_csv(file.path(DATA_DIR, "Decomisos/2021/conjunto_de_datos/m2s1p23_cnspf2021.csv"))
d2021 <- read_csv(file.path(DATA_DIR, "Decomisos/2022/conjunto_de_datos/m2s4p4_cnspf2022.csv"))
d2022 <- read_csv(file.path(DATA_DIR, "Decomisos/2023/conjunto_de_datos/m2bs4p4_cnspf2023.csv"))
d2023 <- read_csv(file.path(DATA_DIR, "Decomisos/2024/conjunto_de_datos/m2bs4p4_cnspf2024.csv"))

# State name lookup for seizure files
entidades_seizures <- c(
  "01"="Aguascalientes","02"="Baja California","03"="Baja California Sur","04"="Campeche",
  "05"="Coahuila de Zaragoza","06"="Colima","07"="Chiapas","08"="Chihuahua",
  "09"="Ciudad de México","10"="Durango","11"="Guanajuato","12"="Guerrero",
  "13"="Hidalgo","14"="Jalisco","15"="México","16"="Michoacán de Ocampo",
  "17"="Morelos","18"="Nayarit","19"="Nuevo León","20"="Oaxaca","21"="Puebla",
  "22"="Querétaro","23"="Quintana Roo","24"="San Luis Potosí","25"="Sinaloa",
  "26"="Sonora","27"="Tabasco","28"="Tamaulipas","29"="Tlaxcala",
  "30"="Veracruz de Ignacio de la Llave","31"="Yucatán","32"="Zacatecas",
  "99"="No especificado"
)

# Build annual state-level dataset for each file format
base2017 <- d2017 %>% transmute(year = 2017L,
  state_code = sprintf("%02d", as.integer(entidad)),
  firearms_total = as.numeric(armexp1),
  firearms_long  = as.numeric(armexp2),
  firearms_short = as.numeric(armexp3))

base2018 <- d2018 %>% transmute(year = 2018L,
  state_code = sprintf("%02d", as.integer(entidad)),
  firearms_total = as.numeric(armexp1),
  firearms_long  = as.numeric(armexp2),
  firearms_short = as.numeric(armexp3))

base2019 <- d2019 %>% transmute(year = 2019L,
  state_code = sprintf("%02d", as.integer(entidad_a)),
  firearms_total = as.numeric(armatis1),
  firearms_long  = as.numeric(armati1),
  firearms_short = as.numeric(armati2))

base2020 <- d2020 %>% transmute(year = 2020L,
  state_code = sprintf("%02d", as.integer(entidad_a)),
  firearms_total = as.numeric(armatis1),
  firearms_long  = as.numeric(armati1),
  firearms_short = as.numeric(armati2))

agg_year_seizures <- function(df, y) {
  df %>%
    mutate(state_code = sprintf("%02d", as.integer(enticve1))) %>%
    group_by(state_code) %>%
    summarise(
      firearms_long    = sum(as.numeric(armasa1), na.rm = TRUE),
      firearms_short   = sum(as.numeric(armasa2), na.rm = TRUE),
      firearms_unknown = sum(as.numeric(armasa3), na.rm = TRUE),
      firearms_total   = firearms_long + firearms_short + firearms_unknown,
      .groups = "drop"
    ) %>%
    mutate(year = y) %>%
    select(year, state_code, firearms_total, firearms_long, firearms_short)
}

decomisos_estado <- bind_rows(
  base2017, base2018, base2019, base2020,
  agg_year_seizures(d2021, 2021L),
  agg_year_seizures(d2022, 2022L),
  agg_year_seizures(d2023, 2023L)
) %>%
  mutate(
    state_code = sprintf("%02d", as.integer(state_code)),
    state_name = unname(entidades_seizures[state_code])
  ) %>%
  filter(!is.na(state_code), state_code != "99", !is.na(state_name)) %>%
  mutate(firearms_total_h = firearms_long + firearms_short) %>%
  left_join(
    read_csv(file.path(DATA_DIR, "poblacion.csv")) %>%
      mutate(state_code = sprintf("%02d", as.integer(state_code))) %>%
      rename(pop2020 = Poblacion2020),
    by = "state_code"
  ) %>%
  mutate(
    rate_total_100k = 1e5 * firearms_total_h / pop2020,
    rate_long_100k  = 1e5 * firearms_long    / pop2020,
    rate_short_100k = 1e5 * firearms_short   / pop2020
  )

# DiD datasets
did_df <- decomisos_estado %>%
  mutate(treated = as.integer(state_code == "09"),
         post    = as.integer(year >= 2019L))

did_main    <- did_df
did_nospill <- did_df %>% filter(!(state_code %in% c("15", "17")))

# Helper: TWFE DiD + event-study
run_seizure_did <- function(data, outcome, label) {
  f_twfe <- as.formula(paste0(outcome, " ~ treated:post | state_code + year"))
  f_es   <- as.formula(paste0(outcome, " ~ i(year, treated, ref = 2018) | state_code + year"))
  list(
    model  = feols(f_twfe, data = data, cluster = ~ state_code),
    es     = feols(f_es,   data = data, cluster = ~ state_code),
    label  = label,
    outcome = outcome
  )
}

seizure_outcomes <- list(
  list(var = "rate_total_100k", lab = "Total (Long+Short) per 100k"),
  list(var = "rate_long_100k",  lab = "Long guns per 100k"),
  list(var = "rate_short_100k", lab = "Handguns per 100k")
)

res_main_sz    <- lapply(seizure_outcomes, function(x) run_seizure_did(did_main,    x$var, x$lab))
res_nospill_sz <- lapply(seizure_outcomes, function(x) run_seizure_did(did_nospill, x$var, x$lab))

# Table 9.1: ATT comparison — full vs. no-spillover sample
etable(
  res_main_sz[[1]]$model, res_nospill_sz[[1]]$model,
  res_main_sz[[2]]$model, res_nospill_sz[[2]]$model,
  res_main_sz[[3]]$model, res_nospill_sz[[3]]$model,
  headers = c("Total — Full", "Total — No spillovers",
              "Long — Full",  "Long — No spillovers",
              "Short — Full", "Short — No spillovers"),
  keep    = "%treated:post",
  dict    = c("treated:post" = "ATT (CDMX × Post-2019)"),
  se.below = TRUE,
  notes    = "SEs clustered by state; two-way FE (state, year)."
)

# Plot helper: unified theme
wd_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(axis.title.x    = element_text(margin = margin(t = 6)),
          axis.title.y    = element_text(margin = margin(r = 6)),
          legend.position = "top",
          panel.grid.minor= element_blank(),
          panel.grid.major= element_line(colour = "#D9D9D9", linewidth = 0.3))
}

# Plot 9.1: Event-study (pre-trends)
plot_event_study <- function(es_obj) {
  broom::tidy(es_obj, conf.int = TRUE) %>%
    filter(grepl("treated", term)) %>%
    mutate(year = as.integer(gsub(".*?(\\d{4}).*", "\\1", term))) %>%
    arrange(year) %>%
    ggplot(aes(x = year, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#7A7A7A") +
    geom_vline(xintercept = 2019,    linetype = "dotted", colour = "#7A7A7A") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#449DD1", alpha = 0.15) +
    geom_point(colour = "#449DD1", size = 2) +
    geom_line(colour  = "#449DD1", linewidth = 0.8) +
    labs(x = "Year", y = "Estimate (95% CI)") +
    wd_theme()
}

print(plot_event_study(res_main_sz[[1]]$es))  # Total
print(plot_event_study(res_main_sz[[2]]$es))  # Long guns
print(plot_event_study(res_main_sz[[3]]$es))  # Handguns

# Plot 9.2: CDMX vs. control means over time
plot_seizure_trend <- function(df, outcome) {
  df %>%
    mutate(group = if_else(state_code == "09", "CDMX", "Controls")) %>%
    group_by(group, year) %>%
    summarise(val = mean(.data[[outcome]], na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(x = year, y = val, group = group, colour = group)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_vline(xintercept = 2019, linetype = "dashed", colour = "#7A7A7A") +
    scale_color_manual(NULL, values = c("CDMX" = "#449DD1", "Controls" = "#7A7A7A")) +
    labs(x = "Year", y = "Rate per 100,000") +
    wd_theme()
}

print(plot_seizure_trend(did_main, "rate_total_100k"))
print(plot_seizure_trend(did_main, "rate_long_100k"))
print(plot_seizure_trend(did_main, "rate_short_100k"))


## =============================================================================
## 10. Validation against INEGI homicide records (ICD-10: X93–X95)
## =============================================================================

# Load INEGI mortality data
inegi_raw <- read_excel(
  file.path(DATA_DIR, "Homicidios INEGI/INEGI_exporta_20_8_2025_16_21_46.xlsx"),
  col_names = FALSE,
  skip      = 4,
  na        = c("ND", "NA", "", ".")
)

hdr1 <- as.character(inegi_raw[1, ])
hdr2 <- as.character(inegi_raw[2, ])
names(inegi_raw) <- c("year", "month",
                       paste(hdr1[-c(1,2)], hdr2[-c(1,2)], sep = " | "))
inegi_raw <- inegi_raw[-c(1,2), ]

inegi_tidy <- inegi_raw %>%
  mutate(
    year      = suppressWarnings(as.integer(as.character(year))),
    month     = str_squish(str_to_title(as.character(month))),
    month_num = case_when(
      month == "Enero" ~ 1L, month == "Febrero" ~ 2L, month == "Marzo"     ~ 3L,
      month == "Abril" ~ 4L, month == "Mayo"    ~ 5L, month == "Junio"      ~ 6L,
      month == "Julio" ~ 7L, month == "Agosto"  ~ 8L, month == "Septiembre" ~ 9L,
      month == "Octubre" ~ 10L, month == "Noviembre" ~ 11L, month == "Diciembre" ~ 12L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(year >= 2017, year <= 2023, !is.na(month_num)) %>%
  pivot_longer(
    cols = -c(year, month, month_num),
    names_to  = c("entidad", "causa"),
    names_sep = "\\s*\\|\\s*",
    values_to = "defunciones"
  ) %>%
  mutate(
    entidad     = str_squish(entidad),
    causa       = str_squish(causa),
    defunciones = suppressWarnings(parse_number(as.character(defunciones)))
  ) %>%
  filter(entidad != "Total", !causa %in% c("Total", "CIE-9", "CIE-10/2"))

# Filter to firearm deaths (X93–X95)
homi_af_inegi <- inegi_tidy %>%
  filter(str_detect(causa, regex("\\bX9[345]", ignore_case = TRUE))) %>%
  mutate(
    firearm_type = case_when(
      str_detect(causa, regex("\\bX93", ignore_case = TRUE)) ~ "Handgun",
      str_detect(causa, regex("\\bX94", ignore_case = TRUE)) ~ "Long gun",
      str_detect(causa, regex("\\bX95", ignore_case = TRUE)) ~ "Other/unspecified",
      TRUE ~ NA_character_
    )
  )

homi_af_collapsed <- homi_af_inegi %>%
  select(entidad, year, month, month_num, defunciones) %>%
  group_by(entidad, year, month, month_num) %>%
  summarise(defunciones = sum(defunciones, na.rm = TRUE), .groups = "drop") %>%
  mutate(fecha = make_date(year, month_num, 1L))

# IDEFC firearm homicides (Doloso + Culposo + Feminicidio, Con arma de fuego)
idefc_af <- delitos %>%
  filter(Año >= 2017, Año <= 2023) %>%
  mutate(Entidad = recode(Entidad,
    "Coahuila de Zaragoza"            = "Coahuila",
    "Michoacán de Ocampo"             = "Michoacán",
    "Veracruz de Ignacio de la Llave" = "Veracruz",
    "Estado de México"                = "México",
    "Ciudad de M??xico"               = "Ciudad de México"
  )) %>%
  filter(`Subtipo de delito` %in% c("Homicidio doloso", "Homicidio culposo", "Feminicidio"),
         Modalidad == "Con arma de fuego") %>%
  group_by(Entidad, Fecha) %>%
  summarise(crimes_af = sum(valor, na.rm = TRUE), .groups = "drop")

# Join both series for Mexico City
deaths_cdmx <- homi_af_collapsed %>%
  mutate(entidad = recode(entidad, "Distrito Federal" = "Ciudad de México")) %>%
  filter(entidad == "Ciudad de México") %>%
  transmute(Fecha = fecha, defunciones_af = defunciones)

combined_cdmx <- idefc_af %>%
  filter(Entidad == "Ciudad de México") %>%
  full_join(deaths_cdmx, by = "Fecha") %>%
  filter(Fecha >= as.Date("2017-01-01"), Fecha <= as.Date("2023-12-01")) %>%
  arrange(Fecha) %>%
  mutate(across(c(crimes_af, defunciones_af), ~ coalesce(., 0)))

# Plot 10.1: IDEFC vs. INEGI monthly counts (Mexico City)
combined_cdmx %>%
  pivot_longer(cols = c(crimes_af, defunciones_af), names_to = "series", values_to = "count") %>%
  mutate(series = recode(series,
    "crimes_af"     = "Firearm Homicides (SESNSP)",
    "defunciones_af"= "Firearm Homicides (INEGI)"
  )) %>%
  ggplot(aes(x = Fecha, y = count, color = series)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Firearm Homicides (SESNSP)" = "#08306B",
                                "Firearm Homicides (INEGI)"  = "#6BAED6")) +
  scale_x_date(date_breaks = "4 months", date_labels = "%m/%Y",
               expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Date", y = "Count", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.title  = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Synthetic DiD — INEGI firearm deaths (X93–X95)
# NOTE: Built manually (not via run_synthdid helper) because def_af_panel
# uses the raw count outcome rather than a rate and does not have a
# Poblacion2020 column required by the helper's panel structure.
set.seed(123)

def_af_panel <- homi_af_collapsed %>%
  transmute(
    Entidad = recode(entidad,
      "Coahuila de Zaragoza"            = "Coahuila",
      "Michoacán de Ocampo"             = "Michoacán",
      "Veracruz de Ignacio de la Llave" = "Veracruz",
      "Estado de México"                = "México",
      .default = entidad
    ),
    Fecha  = as.Date(fecha),
    def_af = as.numeric(defunciones)
  ) %>%
  filter(!is.na(Entidad), !is.na(Fecha)) %>%
  group_by(Entidad, Fecha) %>%
  summarise(def_af = sum(def_af, na.rm = TRUE), .groups = "drop")

treated_unit_inegi <- "Ciudad de México"
donors_exclude_inegi <- c("Morelos", "México")
treat_date_inegi     <- as.Date("2019-01-01")
post_end_inegi       <- as.Date("2020-03-01")

panel_sd_inegi <- def_af_panel %>%
  filter(!Entidad %in% donors_exclude_inegi) %>%
  arrange(Entidad, Fecha)

time_points_inegi   <- sort(unique(panel_sd_inegi$Fecha))
control_units_inegi <- setdiff(unique(panel_sd_inegi$Entidad), treated_unit_inegi)
units_order_inegi   <- c(control_units_inegi, treated_unit_inegi)

Y_inegi <- panel_sd_inegi %>%
  mutate(Entidad = factor(Entidad, levels = units_order_inegi)) %>%
  arrange(Entidad, Fecha) %>%
  pivot_wider(names_from = Fecha, values_from = def_af, values_fill = list(def_af = 0)) %>%
  column_to_rownames("Entidad") %>%
  as.matrix()

keep_dates_inegi <- time_points_inegi[time_points_inegi <= post_end_inegi]
Y_inegi <- Y_inegi[, as.character(keep_dates_inegi)]

N0_inegi <- length(control_units_inegi)
T0_inegi <- which(keep_dates_inegi == treat_date_inegi) - 1

tau_inegi  <- synthdid_estimate(Y_inegi, N0 = N0_inegi, T0 = T0_inegi)
se_p_inegi <- sqrt(as.numeric(vcov(tau_inegi, method = "placebo")))
ci_p_inegi <- as.numeric(tau_inegi) + c(-1, 1) * 1.96 * se_p_inegi

cat("\n=== Synthetic DiD (INEGI def_af — X93–X95) — Post: 2019-01 to 2020-03 ===\n")
cat(sprintf("ATT (post average): %0.4f\n", as.numeric(tau_inegi)))
cat(sprintf("SE (placebo):       %0.4f   | 95%% CI: (%0.4f, %0.4f)\n",
            se_p_inegi, ci_p_inegi[1], ci_p_inegi[2]))
print(summary(tau_inegi))

cat("\n--- Top donor controls (weights) ---\n")
print(head(synthdid_controls(tau_inegi), 15))

sc_inegi  <- sc_estimate (Y_inegi, N0_inegi, T0_inegi)
did_inegi <- did_estimate(Y_inegi, N0_inegi, T0_inegi)

cat(sprintf("\nSC ATT: %0.4f  |  DiD ATT: %0.4f\n",
            as.numeric(sc_inegi), as.numeric(did_inegi)))

print(synthdid_plot(
  list("SDID" = tau_inegi, "SC" = sc_inegi, "DiD" = did_inegi),
  facet = NULL, overlay = 1, ci.alpha = 0.01, line.width = 0.8,
  treated.name = "CDMX", control.name = "Synthetic CDMX"
))

# Dynamic effect curve
post_dates_inegi <- keep_dates_inegi[(T0_inegi + 1):length(keep_dates_inegi)]
eff_df_inegi <- tibble(
  Fecha  = post_dates_inegi,
  Efecto = as.numeric(synthdid_effect_curve(tau_inegi))
)

ggplot(eff_df_inegi, aes(x = Fecha, y = Efecto)) +
  geom_line(linewidth = 1, color = "#449DD1") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_date(date_breaks = "2 months", date_labels = "%m/%Y") +
  labs(x = "Date", y = "Effect (counts)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title  = element_text(face = "bold"))

cat("\n--- Design checks (INEGI def_af, trimmed panel) ---\n")
cat(sprintf("Number of controls (N0): %d\n", N0_inegi))
cat(sprintf("Pre length (T0):         %d months\n", T0_inegi))
cat(sprintf("Total months:            %d\n", ncol(Y_inegi)))
cat(sprintf("Treated unit last row?:  %s\n",
            ifelse(tail(rownames(Y_inegi), 1) != "Ciudad de México",
                   "NO (check ordering!)", "YES")))


## =============================================================================
## 11. Descriptive buyback charts
## =============================================================================

# English state name dictionary for descriptive figures
estado_dict_en <- c(
  "AGS." = "Aguascalientes",      "Aguascalientes" = "Aguascalientes",
  "B.C." = "Baja California",     "Baja California" = "Baja California",
  "B.C.S." = "Baja California Sur","Baja California Sur" = "Baja California Sur",
  "CAMP." = "Campeche",            "Campeche" = "Campeche",
  "CD. MÉX." = "Mexico City",      "D.F." = "Mexico City",
  "Cd. de México" = "Mexico City", "Ciudad de México" = "Mexico City",
  "CHIH." = "Chihuahua",  "Chihuahua" = "Chihuahua",
  "CHIS." = "Chiapas",    "Chiapas"   = "Chiapas",
  "COAH." = "Coahuila",   "Coahuila"  = "Coahuila",
  "COL."  = "Colima",     "Colima"    = "Colima",
  "DGO."  = "Durango",    "Durango"   = "Durango",
  "GTO."  = "Guanajuato", "Guanajuato"= "Guanajuato",
  "GRO."  = "Guerrero",   "Guerrero"  = "Guerrero",
  "HGO."  = "Hidalgo",    "Hidalgo"   = "Hidalgo",
  "JAL."  = "Jalisco",    "Jalisco"   = "Jalisco",
  "MEX."  = "Mexico State","MÉX." = "Mexico State","México" = "Mexico State",
  "MICH." = "Michoacan",  "Michoacán" = "Michoacan",
  "MOR."  = "Morelos",    "Morelos"   = "Morelos",
  "N.L."  = "Nuevo Leon", "Nvo. León" = "Nuevo Leon","Nuevo León" = "Nuevo Leon",
  "NAY."  = "Nayarit",    "Nayarit"   = "Nayarit",
  "OAX."  = "Oaxaca",     "Oaxaca"    = "Oaxaca",
  "PUE."  = "Puebla",     "Puebla"    = "Puebla",
  "QRO."  = "Queretaro",  "Querétaro" = "Queretaro",
  "Q. ROO"= "Quintana Roo","Q. ROO." = "Quintana Roo","Quintana Roo" = "Quintana Roo",
  "S.L.P."= "San Luis Potosi","San Luis Potosí" = "San Luis Potosi",
  "SIN."  = "Sinaloa",    "Sinaloa"   = "Sinaloa",
  "SON."  = "Sonora",     "Sonora"    = "Sonora",
  "TAB."  = "Tabasco",    "Tabasco"   = "Tabasco",
  "TAMPS."= "Tamaulipas", "Tamaulipas"= "Tamaulipas",
  "TLAX." = "Tlaxcala",   "Tlaxcala"  = "Tlaxcala",
  "VER."  = "Veracruz",   "Veracruz"  = "Veracruz",
  "YUC."  = "Yucatan",    "Yucatán"   = "Yucatan",
  "ZAC."  = "Zacatecas",  "Zacatecas" = "Zacatecas"
)

# Plot 11.1: Total buybacks 2017–2023 by state (horizontal bar chart)
canjes_raw %>%
  mutate(
    Fecha   = dmy(fecha_evento),
    Entidad = recode(estado, !!!estado_dict_en),
    total_cl = coalesce(corta, 0) + coalesce(larga, 0)
  ) %>%
  filter(Fecha >= as.Date("2017-01-01"), Fecha <= as.Date("2023-12-31")) %>%
  group_by(Entidad) %>%
  summarise(total_2017_2023 = sum(total_cl, na.rm = TRUE), .groups = "drop") %>%
  filter(!Entidad %in% c("Nacional", "United States of Mexico")) %>%
  ggplot(aes(y = total_2017_2023,
             x = fct_reorder(Entidad, total_2017_2023))) +
  geom_col(fill = "#08306B") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Number of Firearm Buybacks") +
  theme_minimal(base_size = 12) +
  theme(axis.title   = element_text(face = "bold", size = 14),
        legend.position = "none") +
  coord_flip()


# Plot 11.2: CDMX share of national buybacks over time
all_months <- tibble(Fecha = seq(as.Date("2017-01-01"), as.Date("2023-12-01"), by = "month"))

df_monthly <- canjes_raw %>%
  mutate(
    Fecha   = floor_date(dmy(fecha_evento), "month"),
    Entidad = recode(estado, !!!estado_dict_en, .default = estado),
    total_cl = coalesce(corta, 0) + coalesce(larga, 0)
  ) %>%
  filter(Fecha >= as.Date("2017-01-01"), Fecha <= as.Date("2023-12-31")) %>%
  group_by(Entidad, Fecha) %>%
  summarise(total_mes_ent = sum(total_cl, na.rm = TRUE), .groups = "drop") %>%
  group_by(Entidad) %>%
  complete(Fecha = all_months$Fecha, fill = list(total_mes_ent = 0)) %>%
  ungroup()

nat_month  <- df_monthly %>% group_by(Fecha) %>%
  summarise(nacional = sum(total_mes_ent, na.rm = TRUE), .groups = "drop")
cdmx_month <- df_monthly %>% filter(Entidad == "Mexico City") %>%
  select(Fecha, cdmx = total_mes_ent)

share_df <- all_months %>%
  left_join(nat_month,  by = "Fecha") %>%
  left_join(cdmx_month, by = "Fecha") %>%
  mutate(
    cdmx     = replace_na(cdmx, 0),
    nacional = replace_na(nacional, 0),
    share    = if_else(nacional > 0, cdmx / nacional, 0)
  )

cut_date <- as.Date("2019-01-01")
pre_avg  <- share_df %>% filter(Fecha <  cut_date) %>% summarise(m = mean(share, na.rm = TRUE)) %>% pull(m)
post_avg <- share_df %>% filter(Fecha >= cut_date) %>% summarise(m = mean(share, na.rm = TRUE)) %>% pull(m)

ggplot(share_df, aes(x = Fecha, y = share)) +
  geom_line(linewidth = 0.8, color = "#4292C6") +
  geom_vline(xintercept = cut_date, linetype = "dotted", linewidth = 0.7) +
  geom_segment(aes(x = min(share_df$Fecha), xend = cut_date,
                   y = pre_avg, yend = pre_avg),
               inherit.aes = FALSE, linewidth = 0.7, linetype = "dashed") +
  geom_segment(aes(x = cut_date, xend = max(share_df$Fecha),
                   y = post_avg, yend = post_avg),
               inherit.aes = FALSE, linewidth = 0.7, linetype = "dashed") +
  annotate("text", x = cut_date + 200, y = 0.01,
           label = "Program starts", color = "black", fontface = "bold", size = 3.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Date", y = "Share of National Buybacks") +
  theme_minimal(base_size = 12) +
  theme(axis.title   = element_text(face = "bold", size = 14),
        plot.margin  = margin(t = 10, r = 10, b = 30, l = 10)) +
  coord_cartesian(ylim = c(0, 1), clip = "off")


# Plot 11.3: Stacked share — Mexico City vs. rest of country (monthly %)
totals <- nat_month %>%
  left_join(cdmx_month, by = "Fecha") %>%
  mutate(
    cdmx     = replace_na(cdmx, 0),
    nacional = replace_na(nacional, 0),
    resto    = pmax(nacional - cdmx, 0),
    no_data  = nacional == 0
  )

stack_df <- bind_rows(
  totals %>% filter(!no_data) %>%
    transmute(Fecha, `Mexico City` = cdmx / nacional,
              `Rest of Country` = resto / nacional) %>%
    pivot_longer(-Fecha, names_to = "Categoria", values_to = "share"),
  totals %>% filter(no_data) %>%
    transmute(Fecha, Categoria = "No data", share = 1)
) %>%
  mutate(Categoria = factor(Categoria,
    levels = c("Mexico City", "Rest of Country", "No data")))

ggplot(stack_df, aes(x = Fecha, y = share, fill = Categoria)) +
  geom_col() +
  geom_vline(xintercept = as.Date("2019-01-01"),
             linetype = "dotted", linewidth = 1.8, color = "red") +
  annotate("text", x = as.Date("2019-08-01"), y = -0.02,
           label = "Program starts", color = "red", fontface = "bold", size = 3.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(-0.05, 1)) +
  scale_x_date(date_breaks = "3 months", date_labels = "%m/%Y") +
  scale_fill_manual(values = c(
    "Mexico City"     = "#08306B",
    "Rest of Country" = "#6BAED6",
    "No data"         = "#BDBDBD"
  )) +
  labs(x = "Date", y = "Percentage of National Buybacks", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.title   = element_text(face = "bold", size = 14),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.margin  = margin(t = 10, r = 10, b = 30, l = 10)) +
  coord_cartesian(ylim = c(0, 1), clip = "off")


# Table 11.1: Buybacks by phase and type (Mexico City only)
canjes_por_estado_mes %>%
  filter(Entidad == "Ciudad de México") %>%
  mutate(Fase = case_when(
    mes >= as.Date("2017-01-01") & mes <= as.Date("2018-12-31") ~ "Pre (Jan 2017–Dec 2018)",
    mes >= as.Date("2019-01-01") & mes <= as.Date("2020-02-29") ~ "Phase 1 (Jan 2019–Mar 2020)",
    mes >= as.Date("2020-03-01") & mes <= as.Date("2020-11-30") ~ "Phase 2 (Mar–Dec 2020)",
    mes >= as.Date("2020-12-01") & mes <= as.Date("2023-12-31") ~ "Phase 3 (Jan 2021–Dec 2023)",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Fase)) %>%
  group_by(Fase) %>%
  summarise(
    Total    = sum(TOTAL, na.rm = TRUE),
    Handguns = sum(CORTA, na.rm = TRUE),
    Long_guns = Total - Handguns,
    .groups  = "drop"
  ) %>%
  mutate(Fase = factor(Fase, levels = c(
    "Pre (Jan 2017–Dec 2018)", "Phase 1 (Jan 2019–Mar 2020)",
    "Phase 2 (Mar–Dec 2020)",  "Phase 3 (Jan 2021–Dec 2023)"
  ))) %>%
  arrange(Fase)


## =============================================================================
## End of replication script
## =============================================================================
