library(stringr)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

mods <- c("mod_LPcommunity_outbreak_reduced2", "mod_LPcommunity_habitat_reduced",
          "mod_SFcommunity_outbreak_reduced2", "mod_SFcommunity_habitat_reduced")
out <- c("Model_summary_LPoutbrk.csv", "Model_summary_LPhabitat.csv",
         "Model_summary_SFoutbrk.csv", "Model_summary_SFhabitat.csv")

for(i in 1:length(mods)) {
  msum <- loadObject(mods[i])$summary
  rnams <- row.names(msum)
  for(s in spp.list) {
    sind <- which(spp.list == s)
    sind.bracketed <- str_c("[", sind, "]")
    par.ind <- which(str_sub(rnams, -nchar(sind.bracketed), -1) == sind.bracketed)
    repl.vals <- rnams[par.ind] %>% str_sub(1, -1*(nchar(sind.bracketed) + 1)) %>%
      str_c(".", s)
    rnams <- replace(rnams, par.ind, repl.vals)
  }
  row.names(msum) <- rnams
  write.csv(msum, out[i])
}
