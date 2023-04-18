

## helper functions for radar plot
# https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/


# PREPARE PLOT DATA ------------------------------------------------------------

prepare_radar_plot_dat <- function(data, cluster_var) {
  library(dplyr)
  
  radar_plot_dat_prep <- data %>% 
    select({{cluster_var}}, 
           # demographics
           Age, Sex, 
           # ECG & Vitals
           ECG_normal, Heart_Rate, LV_dysfunction,
           # labs
           hsTrop, CRP, Creatinine, Hemoglobine) %>% 
    group_by(eval(as.symbol(cluster_var))) %>%  # group_by({{cluster_var}})   does not work
    summarise(#n=n(), 
      Age = mean(Age), Sex = sum(Sex == "0")/n(),   
      ECG_normal=sum(ECG_normal==1)/n(), Heart_Rate=mean(Heart_Rate), LV_dysfunction = sum(LV_dysfunction=="good")/n(),
      hsTrop= mean(hsTrop), CRP=mean(CRP), Creatinine=mean(Creatinine), Hemoglobine=mean(Hemoglobine)
    ) %>%    
    rename(Row_index=`eval(as.symbol(cluster_var))`)
  
  return(radar_plot_dat_prep)
}


# BEAUTIFUL PLOT ------------------------------------------------------------
# https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  
  library(fmsb)
  
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.3), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

# RADAR PLOT with AVERAGE Background  --------------------------------------------------
create_radar_plot <- function(radar_plot_dat) {
  # Load necessary libraries
  library(scales)
  library(fmsb)
  
  # Rescale each variable to range between 0 and 1
  df_scaled <- round(apply(radar_plot_dat, 2, scales::rescale), 2)
  df_scaled <- as.data.frame(df_scaled)
  
  # Prepare the data for creating the radar plot using the fmsb package
  col_max <- apply(df_scaled, 2, max)
  col_min <- apply(df_scaled, 2, min)
  col_summary <- t(data.frame(Max = col_max, Min = col_min))
  df_scaled2 <- as.data.frame(rbind(col_summary, df_scaled))
  
  # Produce radar plots showing both the average profile and the individual profile
  opar <- par()
  par(mar = rep(0.8, 4))
  par(mfrow = c(3, 4))
  
  for (i in 4:nrow(df_scaled2)) {
    radarchart(
      df_scaled2[c(1:3, i), ],
      pfcol = c("#99999980", NA),
      pcol = c(NA, 2), plty = 1, plwd = 2,
      title = row.names(df_scaled2)[i]
    )
  }
  
  # Restore the standard par() settings
  par(opar)
}