source(here::here("R/simulation.R"))
library(wacolors)

PAL = wa_pal("rainier")

ex_plot(x_conf=0.0, r_bias=0.0)

{
res_00_00 = run_sim(x_conf=0.0, r_bias=0.0)
res_00_00b = run_sim(x_conf=0.0, r_bias=0.0)
res_00_08 = run_sim(x_conf=0.0, r_bias=-0.8)
res_04_00 = run_sim(x_conf=0.4, r_bias=0.0)
res_02_05 = run_sim(x_conf=0.2, r_bias=0.5)
}

rmse(diff_act, diff_est, res_00_08)
rmse(diff_act, diff_est, res_04_00)

ggplot(NULL, aes(diff_est - diff_act)) +
    geom_vline(xintercept=0, lty="dashed") +
    geom_density(data=res_00_08, fill=alpha(PAL[1], 0.5), adjust=1.5) +
    geom_density(data=res_04_00, fill=alpha(PAL[3], 0.5), adjust=1.5) +
    annotate("text", x=-0.08, y=Inf, vjust=1.4, lineheight=1, size=4,
             label="X confouding = 0.4\nRMSE = 0.083") +
    annotate("text", x=0.02, y=Inf, vjust=1.4, lineheight=1, size=4,
             label="R bias = -0.8\nRMSE = 0.029") +
    scale_x_continuous("Absolute Bias", labels=\(x) number(100*x, 0.1, suffix="pp")) +
    scale_y_continuous(NULL, expand=expansion(mult=c(0, 0.1)))

ggsave("~/Desktop/ss.png", width=15, height=7)

res = bind_rows(`00_00` = res_00_00,
          `00_-08` = res_00_08,
          `04_00` = res_04_00,
          `02_05` = res_02_05,
          .id = "config") %>%
    separate(config, c("x_conf", "r_bias"), sep="_") %>%
    mutate(across(1:2, ~ as.numeric(.) / 10),
           label = str_glue("X conf. = {x_conf}\nR bias = {r_bias}"))

res %>%
    pivot_longer(diff_est:diff_corr, names_to="type",
                 values_to="diff", names_prefix="diff_") %>%
ggplot(aes(label, diff - diff_act, fill=type)) +
    geom_hline(yintercept=0, lty="dashed") +
    geom_boxplot() +
    scale_y_continuous("Absolute Bias", labels=\(x) number(100*x, 0.1, suffix="pp")) +
    scale_fill_wa_d("palouse", name="Estimate",
                    labels=c("Corrected", "Uncorrected")) +
    labs(x=NULL)

rmse(diff_act, diff_corr, res_00_00) / rmse(diff_act, diff_est, res_00_00)
rmse(diff_act, diff_corr, res_00_08) / rmse(diff_act, diff_est, res_00_08)
rmse(diff_act, diff_corr, res_04_00) / rmse(diff_act, diff_est, res_04_00)
rmse(diff_act, diff_corr, res_02_05) / rmse(diff_act, diff_est, res_02_05)

res %>%
    #filter(x_conf != 0.4) %>%
ggplot(aes(sample=indep_p, color=label)) +
    geom_qq(distribution=qunif, size=0.5) +
    scale_x_continuous("CDF", labels=percent) +
    scale_y_continuous("p-value for R_est") +
    coord_flip() +
    scale_color_manual(NULL, values=`names<-`(c("black", PAL[1:3]), NULL)) +
    guides(color=guide_legend(override.aes=list(size=4))) +
    theme(legend.key.height=unit(1.5, "cm"),
          legend.key.width=unit(0.9, "cm"),
          legend.text=element_text(size=12))

res2 = run_sim(r_bias=-1)
ex_plot(r_bias=-1)
