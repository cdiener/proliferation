#  analysis.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

library(ggplot2)
library(data.table)
library(survival)
load("combined.rda")

# Panel plot
means <- pred[, list(val=median(rates)), by=panel]
setkey(pred, panel)
pred <- pred[means[order(val), panel]]
pred[, panel := factor(panel, levels=unique(panel))]
panel_plot <- ggplot(pred, aes(x=panel, y=rates, color=tumor, shape=tumor)) + 
    geom_jitter(alpha=0.5, width=0.8, height=0, size=1) + theme_bw() + theme(axis.text.x = 
    element_text(angle = 45, vjust = 1, hjust=1), legend.position="none") + 
    xlab("") + ylab("proliferation rate [1/h]") +
    scale_colour_manual(values=c("royalblue", "red3"))
    
ggsave("panel_rates.png", panel_plot, width=10, height=4, dpi=300)

# Normal vs tumor test
wil <- wilcox.test(pred[!(tumor), rates], pred[(tumor), rates], conf.int=TRUE)
rp <- ggplot(pred, aes(x=c("normal", "tumor")[tumor+1], y=rates, fill=tumor)) + 
    geom_boxplot() + xlab("") + ylab("proliferation rate [1/h]") + theme_bw() + 
    theme(legend.position="none")
ggsave("rate_plot.png", rp, width=3, height=3, dpi=300)

# Survival analysis
comb <- comb[!is.na(vital) & tumor]
delta <- c(comb[vital=="Alive", days_to_contact], 
    comb[vital=="Dead", days_to_death])
status <- comb$vital == "Dead"
prolif <- vector(length=nrow(comb))
prolif[comb$rates > quantile(comb$rates, .75)] <- "high"
prolif[comb$rates < quantile(comb$rates, .25)] <- "low"
prolif <- factor(prolif, levels=c("low", "high"))


surv <- Surv(delta, status)
fit <- survfit(surv ~ prolif)

svg("surv.svg", width=4, height=4)
par(mar=c(3,3,1,1))
plot(fit, col=c("blue", "red"))
grid()
dev.off()

coxm <- coxph(Surv(delta, status) ~ comb$rates)


