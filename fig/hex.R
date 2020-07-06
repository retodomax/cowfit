library(hexSticker)
plot_cow <- function() {
  plot(1,1, type = "n", axes = FALSE, ann = FALSE)
  par(lheight=.35)
  text(0,1.2, labels = "  ^__^\n  (oo)\\_______\n  (__)\\       )\\/\\\n      ||----w |\n      ||     ||",
       adj =0, cex = 2.8, family = "mono")
}
s <- sticker(~plot_cow(),
             package="cowfit", p_size=20, s_x=.8, s_y=.6, s_width=1.4, s_height=1.2,
             h_fill="#F8EBC1", h_color="#3E3120", p_color = "#000000", p_x = 1.3,
             filename="fig/cowfit.png")
# plot(s)
