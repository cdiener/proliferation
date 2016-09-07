#  analysis.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

library(ggplot2)
library(data.table)
library(survival)
library(pheatmap)
library(prtools)

comb <- readRDS("combined.rds")
pred <- fread("pred_rates.csv")

days_per_year <- 365.25

# Panel plot
means <- pred[, list(val=median(rates)), by=panel]
setkey(pred, panel)
pred <- pred[means[order(val), panel]]
pred[, panel := factor(panel, levels=unique(panel))]
panel_plot <- ggplot(pred, aes(x=panel, y=rates, color=tumor, shape=tumor)) +
    geom_jitter(alpha=0.2, width=0.8, height=0, size=1) + theme_bw() + theme(axis.text.x =
    element_text(angle = 45, vjust = 1, hjust=1), legend.position="none") +
    xlab("") + ylab("proliferation rate [1/h]") +
    scale_colour_manual(values=c("royalblue", "red3"))

ggsave("images/fig3.png", panel_plot, width=180, height=70, units="mm", scale=1.2, dpi=300)
pred[, panel := as.character(panel)]

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
plot(fit, col=c("blue", "red"), xlab="time [years]", ylab="survival", lwd=4)
grid()
dev.off()

coxm <- coxph(Surv(delta, status) ~ comb$rates)

# Clean up the panels
comb[, T := gsub("[a-e][0-9]*$", "", T)]
comb[, N := gsub("[^0-4NX]", "", N)]
comb[, M := gsub("[^0-4MX]", "", M)]
comb[stage %in% c("I/II NOS", "IS"), stage := NA]
comb[, stage := gsub("[A-C]*$", "", stage)]

x <- melt(comb[, .(rates, panel, T, N, M, stage)], id.vars=c("rates", "panel"))
med_iqr <- function(x) median_hilow(x, conf.int=0.5)
stage_plot <- ggplot(x, aes(x=value, y=rates, col=variable)) + geom_violin() +
    stat_summary(geom="pointrange", fatten=2, fun.data=med_iqr) +
    facet_wrap(~ variable, scales="free_x", nrow=1) + theme_bw() +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1,
    hjust=1), legend.position="none", strip.text = element_blank()) +
    ylab("proliferation rate [1/h]") + xlab("")
ggsave("images/stage.png", width=180, height=60, units="mm", scale=1.2, dpi=300)

# Get statistics
kw_tests <- x[, kruskal.test(rates, factor(value)), by=variable]
cat("\nTNM association:\n----------------\n")
print(kw_tests)

panels <- pred$panel
names(panels) <- pred$patient_barcode
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
panels <- panels[!duplicated(names(panels))]
fluxes <- fluxes[order(panels[rownames(fluxes)], decreasing=TRUE), ]
in_fluxes <- names(panels) %in% rownames(fluxes)
annrow <- data.frame(panel=panels[in_fluxes], row.names=names(panels)[in_fluxes])
anncolors <- scales::hue_pal()(9)
names(anncolors) <- levels(annrow$panel)
anncolors <- list(panel=anncolors)

s <- seq(-16, log(max(fluxes)+1e-16,2), length.out = 256)
pheatmap(fluxes, breaks=c(-1e-6,2^s), col=cols, show_rownames=F, cellwidth=0.5,
    cellheight=0.08, show_colnames=F, annotation_row=annrow,
    annotation_colors=anncolors, cluster_rows=FALSE, border_color=NA,
    file="images/fluxes.png", width=10, height=6)

lfcs[, panel := factor(panel)]
lfc_plot <- ggplot(lfcs, aes(x=panel, y=lfc, col=panel)) +
    geom_boxplot(outlier.colour=NA) + geom_jitter(width=0.5, alpha=0.25) +
    theme_bw() + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1),
    legend.position="none") + scale_color_discrete() +
    ylab("specificity score")

enr_plot <- ggplot(enr, aes(x=nes, y=subsystem, col=p)) +
    geom_vline(xintercept=1, linetype="dashed") +
    geom_point() + scale_colour_continuous(low="red", high="black") +
    theme_bw() + scale_y_discrete(labels=shorten) +
    theme(legend.position=c(0.8,0.2)) + xlab("enrichment score") + ylab("")

ggsave("images/lfcs.png", lfc_plot, width=120, height=100, units="mm")
ggsave("images/ssea_over.pdf", enr_plot %+% enr[nes > 1], width=90, height=120, units="mm", dpi=300, scale=1.3)
ggsave("images/ssea_under.pdf", enr_plot %+% enr[nes <= 1], width=90, height=70, units="mm", dpi=300, scale=1.5)
