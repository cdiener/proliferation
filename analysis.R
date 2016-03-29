#  analysis.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

library(ggplot2)
library(dplyr)
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

panels <- pred$panel
names(panels) <- pred$barcode

fluxes <- fread("fluxes.csv")
barcodes <- fluxes$V1
fluxes <- as.matrix(fluxes[, V1 := NULL])
rownames(fluxes) <- barcodes
info <- fread("flux_info.csv")


tissue_lfc <- function(v, map, extra) {
    p <- map[rownames(v)]
    
    out <- NULL
    
    out <- lapply(unique(p), function(pname) {
        ti <- p == pname
        in_t <- colMeans(v[ti,]) + 1e-6
        out_t <- colMeans(v[!ti,]) + 1e-6
        lfcs <- log(in_t, 2) - log(out_t, 2) 
        
        res <- data.table(rid=names(in_t), panel=pname, lfc=lfcs)
        cbind(res, extra)
    })
    
    return(rbindlist(out))
}

lfcs <- tissue_lfc(fluxes, panels, info[, list(subsystem)])
summ <- lfcs %>% group_by(rid) %>% summarize(mean=mean(lfc), sd=sd(lfc), 
    sys=unique(subsystem)) %>% mutate(cv=sd/mean) %>% arrange(desc(mean))

spec <- abs(summ$mean) > quantile(abs(summ$mean), 0.9)
core <- summ$cv < quantile(summ$cv, 0.10)

#cols <- viridis::viridis(256)
#s <- seq(-16, log(max(fluxes)+1e-16,2), length.out = 256)
#pheatmap(fluxes, breaks=c(-1e-6,2^s), col=cols, show_rownames=F, 
#    show_colnames=F, annotation_row=as.data.frame(panels), cluster_rows=F, 
#    file="fluxes.png", width=10, height=12)

lfc_plot <- ggplot(lfcs, aes(x=panel, y=lfc, col=panel)) + 
    geom_boxplot(outlier.colour=NA) + geom_jitter(width=0.5, alpha=0.25) + 
    theme_bw()

ggsave("lfcs.svg", lfc_plot, width=6, height=4)
