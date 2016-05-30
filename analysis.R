#  analysis.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

library(ggplot2)
library(dplyr)
library(data.table)
library(survival)
library(pheatmap)
load("combined.rda")

stage_map <- c("0"="0", "I"="1", "II"="2", "III"="3", "IV"="4")
days_per_year <- 365.25

# Panel plot
means <- pred[, list(val=median(rates)), by=panel]
setkey(pred, panel)
pred <- pred[means[order(val), panel]]
pred[, panel := factor(panel, levels=unique(panel))]
panel_plot <- ggplot(pred, aes(x=panel, y=rates, color=tumor, shape=tumor)) +
    geom_jitter(alpha=0.4, width=0.8, height=0, size=1) + theme_bw() + theme(axis.text.x =
    element_text(angle = 45, vjust = 1, hjust=1), legend.position="none") +
    xlab("") + ylab("proliferation rate [1/h]") +
    scale_colour_manual(values=c("royalblue", "red3"))

ggsave("images/fig3.png", panel_plot, width=180, height=70, units="mm", scale=1.2,  dpi=300)

# Normal vs tumor test
wil <- wilcox.test(pred[!(tumor), rates], pred[(tumor), rates], conf.int=TRUE)
rp <- ggplot(pred, aes(x=c("normal", "tumor")[tumor+1], y=rates, fill=tumor)) +
    geom_boxplot() + xlab("") + ylab("proliferation rate [1/h]") + theme_bw() +
    theme(legend.position="none")
ggsave("images/rate_plot.svg", rp, width=80, height=100, units="mm")

# Survival analysis
comb <- comb[!is.na(vital) & tumor]
delta <- c(comb[vital=="Alive", days_to_contact/days_per_year],
    comb[vital=="Dead", days_to_death/days_per_year])
status <- comb$vital == "Dead"
prolif <- vector(length=nrow(comb))
prolif[comb$rates > quantile(comb$rates, .75)] <- "high"
prolif[comb$rates < quantile(comb$rates, .25)] <- "low"
prolif <- factor(prolif, levels=c("low", "high"))


surv <- Surv(delta, status)
fit <- survfit(surv ~ prolif)

svglite::svglite("images/surv.svg", width=5, height=4.5)
par(mar=c(4,4,1,1))
plot(fit, col=c("blue", "red"), xlab="time [years]", ylab="survival")
grid()
dev.off()

coxm <- coxph(Surv(delta, status) ~ comb$rates)

#load("combined.rda")
#comb = comb[(tumor)]
# Clean up the panels
comb[, T := substr(T,2,2)]
comb[, N := substr(N,2,2)]
comb[, M := substr(M,2,2)]
comb[, stage := stage_map[gsub("Stage\\s+|[A-C][0-9]*", "", stage)]]

T_plot <- ggplot(comb, aes(x=T, y=rates)) + geom_boxplot(outlier.colour=NA) +
    geom_jitter(alpha=0.25, height=0, width=0.5, col="black", stroke=0) +
    ylim(-0.005, 0.03) + theme_bw() + ylab("proliferation rate [1/h]") +
    xlab("primary tumor [T stage]")
ggsave("images/T.svg", width=120, height=100, units="mm")

x <- melt(comb[, .(rates, N, M, stage)], id.vars="rates")
stage_plot <- ggplot(x, aes(x=value, y=rates, col=variable)) +
    geom_boxplot(outlier.colour=NA) +
    geom_jitter(alpha=0.5, height=0, width=0.5) + ylim(-0.005, 0.04) +
    facet_wrap(~ variable, scales="free_x", nrow=1) + theme_bw() +
    theme(legend.position="none") + ylab("proliferation rate [1/h]")
ggsave("images/stage.png", width=180, height=70, units="mm", scale=1.2, dpi=300)

panels <- pred$panel
names(panels) <- pred$barcode
panels <- sort(panels)

fluxes <- fread("fluxes.csv")
barcodes <- fluxes$V1
fluxes <- as.matrix(fluxes[, V1 := NULL])
rownames(fluxes) <- barcodes
info <- fread("flux_info.csv")


lfcs <- tissue_lfc(fluxes, panels, info[, list(subsystem)])
lfcs <- lfcs[order(-abs(lfc))]

pws <- lfcs$subsystem
enr <- sapply(unique(pws), NES, w=lfcs$lfc, pws=pws)
enr <- data.table(subsystem=colnames(enr), nes=enr[1,], p=enr[2,])
enr <- enr[order(nes)]
enr[, subsystem := factor(subsystem, levels=subsystem)]

cols <- viridis::viridis(256)
s <- seq(-16, log(max(fluxes)+1e-16,2), length.out = 256)
pheatmap(fluxes, breaks=c(-1e-6,2^s), col=cols, show_rownames=F,
    show_colnames=F, annotation_row=as.data.frame(panels), cluster_rows=F,
    file="images/fluxes.png", width=7, height=4)

lfc_plot <- ggplot(lfcs, aes(x=panel, y=lfc, col=panel)) +
    geom_boxplot(outlier.colour=NA) + geom_jitter(width=0.5, alpha=0.25) +
    theme_bw() + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1),
    legend.position="none") + ylab("specificity score")

enr_plot <- ggplot(enr, aes(x=nes, y=subsystem, col=p)) +
    geom_vline(xintercept=1, linetype="dashed") +
    geom_point() + scale_colour_continuous(low="red", high="black") +
    theme_bw() + scale_y_discrete(labels=shorten) +
    theme(legend.position=c(0.8,0.2)) + xlab("enrichment score") + ylab("")

ggsave("images/lfcs.png", lfc_plot, width=120, height=100, units="mm")
ggsave("images/ssea_over.pdf", enr_plot %+% enr[nes > 1], width=90, height=120, units="mm", dpi=300, scale=1.3)
ggsave("images/ssea_under.pdf", enr_plot %+% enr[nes <= 1], width=90, height=70, units="mm", dpi=300, scale=1.3)
