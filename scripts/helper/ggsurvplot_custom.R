
ggsurvplot_custom <- function(fit_object, mytitle="", legend.title, myXbreaks=c(0,1,2,3,4,5), legend.labs, 
                              plot_inlay=TRUE,
                              save=TRUE, plotname) {
  ## example: 
  ## legend.title = "ESC HCM Risk-SCD Score" OR ""
  ## legend.labs = c("Genotype Negative", "Genotype Positive")
  
  ## start with main plot ---------------------------------------
  plot1 <- ggsurvplot(fit_object, 
                      conf.int=FALSE, 
                      pval=FALSE, 
                      censor=FALSE,  # not wanted by BM
                      
                      # risk table
                      risk.table=TRUE, #xlim=c(0,20),
                      #risk.table = "abs_pct",
                      tables.theme = theme_cleantable(),
                      #risk.table.height=.15,
                      risk.table.fontsize = 4,
                      risk.table.title = "No. at Risk",
                      risk.table.y.text = FALSE,
                      risk.table.height = 0.15,
                      #tables.height= 0.2,
                      risk.table.pos = "out",  # default
                      
                      title=mytitle,
                     
                      # legend
                      legend.labs=legend.labs, 
                      legend.title=legend.title,

                      subtitle ="",
                      xlab = "Time in Years", 
                      ylab = "Survival propability\nAll-cause Death",
                      # xlim = c(0, xlimyears),
                      ylim = c(0,1),
                      font.x = c(12),
                      font.y = c(12),
                      font.tickslab = c(12), 
                      font.legend = c(12),
                      
                      #break.x.by = 10,  # affects both plot and table
                      # style
                      palette=c("nejm"),
                      #ggtheme = ggthemes::theme_few()
  )
  # customize further ---------------------------------------
  plot1$plot$theme$text <- element_text(family = "Arial")
  plot1$plot <- plot1$plot + 
    scale_x_continuous(breaks = myXbreaks) #,   labels = c(30,50))
  plot1$table$theme$plot.title <-  element_text(size = 12, color = "black", face = "bold", hjust = -0.1)
  plot1$plot <- plot1$plot + 
    theme_survminer(font.x = c(12, "bold", "black"),
                    font.y = c(12, "bold", "black")) 
  
  
  if (plot_inlay==TRUE) {
    ## add inlay plot ---------------------------------------
    subplot <- ggsurvplot(fit_object, 
                          conf.int=FALSE, 
                          pval=FALSE, 
                          censor=FALSE,
                          
                          # risk table
                          risk.table=FALSE, #xlim=c(0,20),
                          #risk.table = "abs_pct",
                          tables.theme = theme_cleantable(),
                          #risk.table.height=.15,
                          # legend
                          legend = "none",
                          legend.title="",
                          title="",
                          subtitle ="",
                          xlab = "", 
                          ylab = "",
                          # xlim = c(0, 18),
                          ylim = c(0.75,1),
                          font.x = c(12),
                          font.y = c(12),
                          font.tickslab = c(12), 
                          font.legend = c(12),
                          
                          # style
                          palette=c("nejm"),
                          #ggtheme = ggthemes::theme_few()
    )
    
    # make inlay plot transparent ---------------------------------------
    subplot$plot <- subplot$plot + 
      scale_y_continuous(breaks = c(0.7, 0.8, 0.9, 1), limits = c(0.6,1)) +
      scale_x_continuous(breaks = myXbreaks)+
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent")
      )
    
    ## add plots  ---------------------------------------
    plot1$plot <- plot1$plot +
      annotation_custom(grob = ggplotGrob(subplot$plot),   # dont get why smaller plot suddenly scales wrong...
                        # position of subplot
                        xmin = 0.2, xmax=2.3, ymin = -0.05 , ymax = 0.8)
  }
  
  if (save == TRUE) {
    #save
    ## https://github.com/kassambara/survminer/issues/152
    filename.p1 <- paste0("/mnt/users/reich/rockerprojects/bestageing2022/output/plots/survival/", Sys.Date(), "-survival-", plotname, ".png")
    # this works on my mac
    ggsave(
      filename = filename.p1,
      #plot = print(p1, newpage = FALSE),
      survminer:::.build_ggsurvplot(plot1),
      device = 'png',
      width = 10,
      height = 8
    )
    filename.p1 <- paste0("/mnt/users/reich/rockerprojects/bestageing2022/output/plots/survival/", Sys.Date(), "-survival-", plotname, ".svg")
    ggsave(
      filename = filename.p1,
      #plot = print(p1, newpage = FALSE),
      survminer:::.build_ggsurvplot(plot1),
      device = 'svg',
      width = 10,
      height = 7
    )
  }
  
  return(plot1)
}
