# =============================================================================
# maxSSMD overview figure: "Interpretable, Simple, and Broadly Applicable"
# Requires: ggplot2 (plus base 'grid').  Output: maxSSMD_figure.tiff (300 dpi, LZW)
# =============================================================================
library(grid)
library(ggplot2)

# ---- palette -----------------------------------------------------------------
navy   <- "#1B2A6B"; red    <- "#E0262B"
pal <- list(
  blue   = list(fill = "#EAF3FC", line = "#A9CBEE", head = "#1B2A6B", acc = "#2A6EBB"),
  green  = list(fill = "#EDF7EC", line = "#A6D5A3", head = "#1E6B33", acc = "#2E8B3E"),
  purple = list(fill = "#F2ECFB", line = "#C3AEE6", head = "#4B2390", acc = "#6A3DB8"),
  yellow = list(fill = "#FDF7E6", line = "#F1D58A", head = "#1B2A6B", acc = "#D9A21B")
)
txt_col <- "#222222"

# ---- small drawing helpers ----------------------------------------------------
box <- function(x0, y0, x1, y1, fill, col, r = 0.012, lwd = 1.2) {
  grid.roundrect(x = x0, y = y0, width = x1 - x0, height = y1 - y0,
                 just = c("left", "bottom"), r = unit(r, "snpc"),
                 gp = gpar(fill = fill, col = col, lwd = lwd))
}
check <- function(x, y, col, s = 0.011) {           # drawn check mark (no font needed)
  grid.lines(x = c(x, x + 0.35 * s, x + s), y = c(y, y - 0.9 * s, y + 0.9 * s),
             gp = gpar(col = col, lwd = 2.6, lineend = "round", linejoin = "round"))
}
bullets <- function(x, y, items, col, size = 11, dy = 0.052) {
  for (i in seq_along(items)) {
    yy <- y - (i - 1) * dy
    check(x, yy + 0.004, col)
    grid.text(items[[i]], x = x + 0.017, y = yy + 0.012, just = c("left", "top"),
              gp = gpar(fontsize = size, col = txt_col, lineheight = 0.95))
  }
}
draw_gg <- function(p, x0, y0, x1, y1) {
  pushViewport(viewport(x = x0, y = y0, width = x1 - x0, height = y1 - y0,
                        just = c("left", "bottom")))
  grid.draw(ggplotGrob(p)); popViewport()
}
elbow <- function(x, y, col) {                       # connector with arrow head
  grid.lines(x = x, y = y, gp = gpar(col = col, lwd = 2.2, fill = col),
             arrow = arrow(length = unit(0.10, "in"), type = "closed"))
}

theme_mini <- theme_classic(base_size = 9) +
  theme(plot.background = element_rect(fill = NA, colour = NA),
        panel.background = element_rect(fill = NA, colour = NA),
        legend.background = element_rect(fill = NA, colour = NA),
        axis.text = element_text(colour = txt_col), axis.title = element_text(colour = txt_col))

# ---- the four mini plots -------------------------------------------------------
# 1. Coefficients
coef_df <- data.frame(f = factor(c("Feature 1", "Feature 2", "Feature 3", "Feature 4", "...", " "),
                                 levels = rev(c("Feature 1", "Feature 2", "Feature 3", "Feature 4", "...", " "))),
                      b = c(0.92, -0.95, -0.95, -0.95, -0.72, -0.18))
p_coef <- ggplot(coef_df, aes(b, f)) +
  geom_col(fill = "#6F9BD8", colour = "#2F5E9E", width = 0.6, linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  scale_x_continuous(limits = c(-1.05, 1.05), breaks = c(-1, 0, 1)) +
  labs(x = "Coefficient value", y = NULL) + theme_mini +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 9.5))

# 2. Effect size (two densities)
xs <- seq(-3, 5.2, length.out = 400); m1 <- 0; m2 <- 2.3; s <- 0.75
dens <- rbind(data.frame(x = xs, d = dnorm(xs, m1, s), g = "Group 1"),
              data.frame(x = xs, d = dnorm(xs, m2, s), g = "Group 2"))
ytop <- max(dens$d)
p_eff <- ggplot(dens, aes(x, d, colour = g, fill = g)) +
  geom_area(alpha = 0.35, position = "identity", linewidth = 0.8, key_glyph = "path") +
  geom_segment(aes(x = m1, xend = m1, y = 0, yend = ytop), colour = "#1F5FAE",
               linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE) +
  geom_segment(aes(x = m2, xend = m2, y = 0, yend = ytop), colour = "#C8262C",
               linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE) +
  annotate("segment", x = m1, xend = m2, y = ytop * 1.12, yend = ytop * 1.12,
           arrow = arrow(ends = "both", length = unit(0.07, "in")), linewidth = 0.5) +
  annotate("text", x = (m1 + m2) / 2, y = ytop * 1.24, label = "SSMD", size = 3.4) +
  scale_colour_manual(values = c("Group 1" = "#1F5FAE", "Group 2" = "#D7262C")) +
  scale_fill_manual(values = c("Group 1" = "#6E9FDB", "Group 2" = "#F08A8D")) +
  scale_y_continuous(limits = c(0, ytop * 1.32), expand = c(0, 0)) +
  labs(x = "Response (linear combination)", y = "Density", colour = NULL, fill = NULL) +
  theme_mini +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = c(0.86, 0.78), legend.key.height = unit(0.32, "cm"),
        legend.key.width = unit(0.5, "cm"), legend.text = element_text(size = 9)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))

# 3. Confidence interval
p_ci <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  annotate("errorbarh", y = 0, xmin = -1, xmax = 1.35, height = 0.35,
           colour = pal$purple$acc, linewidth = 1.1) +
  annotate("point", x = 0.2, y = 0, size = 3.5, colour = "#4B2390") +
  scale_x_continuous(limits = c(-2.3, 2.3), breaks = -2:2) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "SSMD (95% CI)", y = NULL) + theme_mini +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# 4. ROC curve (binormal, AUC = 0.85)
a   <- sqrt(2) * qnorm(0.85)
fpr <- seq(0, 1, length.out = 300)
roc <- data.frame(fpr = fpr, tpr = pnorm(a + qnorm(fpr)))
p_roc <- ggplot(roc, aes(fpr, tpr)) +
  annotate("segment", x = 0, y = 0, xend = 1, yend = 1, linetype = "dashed", colour = "grey65") +
  geom_line(colour = "#D7262C", linewidth = 1.1) +
  annotate("text", x = 0.68, y = 0.2, label = "AUC = 0.85", size = 3.4) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.01, 0)) +
  coord_equal() + labs(x = "1 - Specificity", y = "Sensitivity") + theme_mini +
  theme(panel.border = element_rect(fill = NA, colour = txt_col, linewidth = 0.5))

# ---- icons -----------------------------------------------------------------------
icon_list <- function(x, y, s, col) {
  pushViewport(viewport(x, y, width = unit(s, "snpc"), height = unit(s, "snpc")))
  grid.roundrect(0.5, 0.5, width = 0.72, height = 0.95, r = unit(0.12, "snpc"),
                 gp = gpar(fill = NA, col = col, lwd = 2.5))
  for (k in 1:3) {
    yy <- 0.72 - (k - 1) * 0.22
    grid.circle(0.32, yy, r = 0.04, gp = gpar(fill = col, col = NA))
    grid.lines(c(0.42, 0.72), c(yy, yy), gp = gpar(col = col, lwd = 2.5, lineend = "round"))
  }
  popViewport()
}
icon_people <- function(x, y, s, col) {
  person <- function(px, py, sc) {
    grid.circle(px, py + 0.22 * sc, r = 0.13 * sc, gp = gpar(fill = col, col = NA))
    grid.polygon(px + sc * c(-0.24, -0.24, -0.16, 0.16, 0.24, 0.24),
                 py + sc * c(-0.30, -0.08, 0.04, 0.04, -0.08, -0.30),
                 gp = gpar(fill = col, col = NA))
  }
  pushViewport(viewport(x, y, width = unit(s, "snpc"), height = unit(s, "snpc")))
  person(0.22, 0.45, 0.75); person(0.78, 0.45, 0.75); person(0.5, 0.40, 1)
  popViewport()
}

# ---- compose the figure ------------------------------------------------------------
draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  # Title banner
  box(0.005, 0.905, 0.995, 0.99, fill = "#FDE4E6", col = "#F4B8BD", r = 0.06)
  t1 <- "maxSSMD:"; t2 <- " Interpretable, Simple, and Broadly Applicable"
  gp_t <- gpar(fontsize = 24, fontface = "bold")
  w1 <- convertWidth(grobWidth(textGrob(t1, gp = gp_t)), "npc", TRUE)
  w2 <- convertWidth(grobWidth(textGrob(t2, gp = gp_t)), "npc", TRUE)
  x0 <- 0.5 - (w1 + w2) / 2
  grid.text(t1, x0, 0.948, just = "left", gp = gpar(fontsize = 24, fontface = "bold", col = navy))
  grid.text(t2, x0 + w1, 0.948, just = "left", gp = gpar(fontsize = 24, fontface = "bold", col = red))

  # Equation box
  box(0.30, 0.815, 0.70, 0.89, fill = "#F2F2F2", col = "#D0D0D0", r = 0.08)
  grid.text(expression(italic(S) == beta[0] + beta[1] * italic(X)[1] + cdots + beta[p] * italic(X)[p]),
            0.5, 0.851, gp = gpar(fontsize = 24, col = txt_col))

  # Panels
  px  <- list(c(0.005, 0.245), c(0.255, 0.495), c(0.505, 0.745), c(0.755, 0.995))
  pc  <- list(pal$blue, pal$green, pal$purple, pal$yellow)
  ttl <- c("Coefficients", "Effect size", "Confidence interval", "AUROC")
  sub <- c("(interpretable feature contributions)", "(ratio of mean to variability)",
           "(uncertainty quantification)", "(ML performance key index)")
  bl  <- list(c("Quantifies contribution of each feature", "Facilitates biological interpretation"),
              c("Standardized measure (mean/variability)", "Enables comparison across assays\nand studies"),
              c("Rigorous statistical inference", "Quantifies uncertainty for decision\nmaking"),
              c("Equivalent to maximizing SSMD", "Optimizes overall sensitivity\nand specificity"))
  plots <- list(p_coef, p_eff, p_ci, p_roc)
  py0 <- 0.31; py1 <- 0.765
  for (i in 1:4) {
    a0 <- px[[i]][1]; a1 <- px[[i]][2]; cc <- pc[[i]]
    box(a0, py0, a1, py1, fill = cc$fill, col = cc$line)
    grid.text(ttl[i], (a0 + a1) / 2, 0.733, gp = gpar(fontsize = 15, fontface = "bold", col = cc$head))
    grid.text(sub[i], (a0 + a1) / 2, 0.697, gp = gpar(fontsize = 11, col = txt_col))
    pw <- if (i == 4) c(0.03, 0.03) else c(0.006, 0.006)
    draw_gg(plots[[i]], a0 + pw[1], 0.43, a1 - pw[2], 0.675)
    bullets(a0 + 0.012, 0.408, bl[[i]], col = cc$acc)
  }

  # Connectors from the equation to the panels
  elbow(c(0.30, 0.125, 0.125), c(0.83, 0.83, 0.768), pal$blue$acc)
  elbow(c(0.39, 0.39, 0.375, 0.375), c(0.815, 0.795, 0.795, 0.768), pal$green$acc)
  elbow(c(0.60, 0.60, 0.625, 0.625), c(0.815, 0.795, 0.795, 0.768), pal$purple$acc)
  elbow(c(0.70, 0.875, 0.875), c(0.83, 0.83, 0.768), pal$yellow$acc)

  # Bottom-left: Simple and Transparent
  box(0.005, 0.02, 0.44, 0.285, fill = "#EDF7EC", col = "#A6D5A3")
  grid.text("Simple and Transparent", 0.02, 0.252, just = "left",
            gp = gpar(fontsize = 15, fontface = "bold", col = pal$green$head))
  icon_list(0.06, 0.13, 0.16, "#2C7A8C")
  bullets(0.115, 0.205, c("A linear combination of imaging features",
                          "Easy to implement and interpret", "Not a black-box model"),
          col = pal$green$acc, dy = 0.055)

  # Bottom-right: Broadly Applicable
  box(0.45, 0.02, 0.995, 0.285, fill = "#EAF3FC", col = "#A9CBEE")
  grid.text("Broadly Applicable", 0.465, 0.252, just = "left",
            gp = gpar(fontsize = 15, fontface = "bold", col = pal$blue$head))
  icon_people(0.51, 0.12, 0.16, "#2A6EBB")
  bullets(0.555, 0.205, c("Applicable to both independent and correlated samples",
                          "Works for homoscedastic and heteroscedastic settings",
                          "Generalizable across assays, cell types, imaging platforms, and perturbation technologies"),
          col = pal$blue$acc, dy = 0.055)
}

# ---- write to files ----------------------------------------------------------------
W <- 16; H <- 6.65   # inches
tiff("maxSSMD_figure.tiff", width = W, height = H, units = "in", res = 300,
     compression = "lzw", type = "cairo")
draw_figure(); dev.off()
