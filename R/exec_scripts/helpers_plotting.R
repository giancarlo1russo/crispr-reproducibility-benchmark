#################################
save_plot_pdf <- function(cfg, filename, plot_obj = NULL, width = 10, height = 10, code = NULL) {
  pdf(file.path(cfg$plot_dir, filename), width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  
  if (!is.null(code)) {
    plot(code)
  } else {
    print(plot_obj)
  }
}

#################################
combine_plots <- function(filename,cfg,width,height,p_l, main_layout_string,
                          plot_layout_w = NULL, plot_layout_h = NULL) {
  n_plots <- length(p_l)
  if (is.null(plot_layout_h)) {
    plot_layout_h <- rep(1,n_plots)
  }
  if (is.null(plot_layout_w)) {
    plot_layout_w <- rep(1,n_plots)
  }
  
wrapped_plot <- wrap_plots(p_l, design = main_layout_string) + 
    plot_annotation(tag_levels = 'A') + 
  plot_layout(widths = plot_layout_w,heights = plot_layout_h)
pdf(file.path(cfg$plot_dir,filename),height=height,width=width)
plot(wrapped_plot)
dev.off()
}
#################################
combine_plots_old <- function(filename,cfg,width,height,p_l, main_layout_string = NULL, 
                          plot_layout_w = c(1,1), plot_layout_h = c(1,1)) {
  # 1. Fallback: If no layout string is provided, use patchwork's auto-wrap
  if (is.null(main_layout_string)) {
    combined <- wrap_plots(p_l)
  } else {
    combined <- eval(parse(text = main_layout_string), envir = list(p = p_l))
  }
  
  # 3. Safely apply plot_layout() options passed via the ellipsis (...)
  # e.g., widths, heights, guides = "collect"
  combined <- combined + plot_layout(widths = plot_layout_w,
                                     heights = plot_layout_h) +
    plot_annotation(tag_levels = "A") +
    theme(plot.tag = element_text(face = "bold", size = 14))
  ggsave(file.path(cfg$plot_dir,filename), plot=combined, width = width, height = height)
}

