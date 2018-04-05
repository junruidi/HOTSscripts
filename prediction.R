###################################################################
##                  Prediction Mortality                         ##
##                        04/05/2018                             ##
###################################################################

rm(list = ls())
load("data/survivaforall.rda")
source("HOTSscripts/fwd_select.R")


hosvd_sc = c(paste0("sc2_",c(1:12)),paste0("hosvd_sc3_",c(1:12)),paste0("hosvd_sc4_",c(1:12)))
hosvd_sc_ct = c(paste0("sc2_",c(1:12)),paste0("hosvd_sc_y_convert3_",c(1:12)),paste0("hosvd_sc_y_convert4_",c(1:12)))

rob_sc = c(paste0("sc2_",c(1:12)),paste0("rob_sc3_",c(1:12)),paste0("rob_sc4_",c(1:12)))
rob_sc_ct = c(paste0("sc2_",c(1:12)),paste0("rob_sc_y_convert3_",c(1:12)),paste0("rob_sc_y_convert4_",c(1:12)))


result_hosvd_sc = fwd_select(data = survival_all,ind_vars = hosvd_sc, nboot = 1000,seed.start = 1234)
warnings()
save(result_hosvd_sc,file = "result_hosvd_sc.rda")

result_hosvd_sc_ct = fwd_select(data = survival_all,ind_vars = hosvd_sc_ct, nboot = 1000,seed.start = 1234)
warnings()
save(result_hosvd_sc_ct,file = "result_hosvd_sc_ct.rda")

result_rob_sc = fwd_select(data = survival_all,ind_vars = rob_sc, nboot = 1000,seed.start = 1234)
warnings()
save(result_rob_sc,file = "result_rob_sc.rda")

result_rob_sc_ct = fwd_select(data = survival_all,ind_vars = rob_sc_ct, nboot = 1000,seed.start = 1234)
warnings()
save(result_rob_sc_ct,file = "result_rob_sc_ct.rda")