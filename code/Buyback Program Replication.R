# --- Minimal GitHub loader (CSV/RDS, processed/csv/rds, refs/heads/main or main) ---
suppressPackageStartupMessages({
  library(readr)
  library(rlang)
})

read_remote_dataset <- function(basename,
                                repo     = "marcomna/BuybackProgram_CDMX",
                                folders  = c("data/processed","data/csv","data/rds","data"),
                                exts     = c("csv","rds"),
                                branches = c("refs/heads/main","main")) {
  
  combos <- expand.grid(branch = branches, folder = folders, ext = exts, stringsAsFactors = FALSE)
  
  # Try processed CSV first (most common), then the rest
  ord <- with(combos, order(!(folder == "data/processed" & ext == "csv")))
  combos <- combos[ord, , drop = FALSE]
  
  urls <- with(combos, sprintf(
    "https://raw.githubusercontent.com/%s/%s/%s/%s.%s",
    repo, branch, folder, basename, ext
  ))
  
  tried <- character(0)
  for (u in urls) {
    tried <- c(tried, u)
    if (grepl("\\.csv$", u, ignore.case = TRUE)) {
      obj <- try(readr::read_csv(u, show_col_types = FALSE), silent = TRUE)
      if (!inherits(obj, "try-error")) { message("✓ Loaded: ", u); return(obj) }
    } else {
      con <- try(url(u, "rb"), silent = TRUE)
      if (inherits(con, "try-error")) next
      on.exit(try(close(con), silent = TRUE), add = TRUE)
      obj <- try(readRDS(con), silent = TRUE)
      if (!inherits(obj, "try-error")) { message("✓ Loaded: ", u); return(obj) }
    }
  }
  
  rlang::abort(paste0(
    "Could not fetch '", basename, "' from GitHub in any supported location.\n",
    "Tried:\n - ", paste(tried, collapse = "\n - "), "\n",
    "Check exact file names, folder, and branch."
  ))
}

# ---- Files expected on GitHub (names = objects you want in your session) ----
expected <- c(
  delitos             = "crime_long_2017_2023",
  poblacion_2020      = "population_2020",
  defunciones_noespec = "unclassified_deaths_monthly",
  unemp_monthly       = "unemployment_monthly",
  ipi_monthly         = "ipi_monthly",
  delitosCS           = "panel_crime_predictors"
)

# ---- Load and assign into the global environment (avoids purrr's map) ----
for (nm in names(expected)) {
  assign(nm, read_remote_dataset(expected[[nm]]), envir = .GlobalEnv)
}
