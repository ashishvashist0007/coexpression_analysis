# =============================================================================
# Run every step in order (instructor demo, or to catch up).
#   Rscript scripts/run_all.R        (from the workshop folder)
# =============================================================================
t0 <- Sys.time()
for (s in c("00_check_setup.R", "01_prepare_data.R", "02_network_modules.R",
            "03_module_traits_hubs.R", "04_genie3.R", "05_guilt_by_association.R")) {
  cat("\n================", s, "================\n")
  source(file.path("scripts", s), echo = FALSE, local = new.env())
}
cat("\nAll steps finished in", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
