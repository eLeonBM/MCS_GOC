
# Loading packages --------------------------------------------------------

library(tidyverse)
library(dplyr) # A staple for modern data management in R
library(lubridate) # Useful functions for dealing with dates
library(ggplot2) # The preferred library for data visualisation
library(tidync) # For easily dealing with NetCDF data
library(rerddap) # For easily downloading subsets of data
library(doParallel) # For parallel processing
library(heatwaveR) # For detecting MHWs
library(sf)
library(rnaturalearth)

# Loading data ------------------------------------------------------------

OISST_data <- readRDS('data/OISST_GoC.RDS')

goc_pol <- sf::st_read('data/shp/GOC_ecoregions.shp') # GOC polygon to be used for intersect data points

spdf_mx <- st_transform(st_as_sf(ne_countries(country = 'mexico')), crs = 4326)


# Function to detect marine heatwaves -------------------------------------

event_only <- function(df){
        require(heatwaveR)
        # First calculate the climatologies
        clim <- ts2clm(data = df, climatologyPeriod = c("1982-01-01", "2021-12-31"), pctile = 10)
        # Then the events
        event <- detect_event(data = clim, coldSpells = T)
        # Return only the event metric dataframe of results
        return(event$event)
}



# Detecting marine heatwaves in the Mexican Pacific ----------------------------------------------

# NB: One should never use ALL available cores, save at least 1 for other essential tasks
# The computer I'm writing this vignette on has 8 cores, so I use 7 here
registerDoParallel(detectCores()-1)

system.time(
        MCS_result <- plyr::ddply(.data = OISST_data,
                                  .variables = c("longitude", "latitude"), 
                                  .fun = event_only, 
                                  .parallel = TRUE)
) # 46.81 seconds




# Summarizing data --------------------------------------------------------

# summarise the number of unique longitude, latitude and year combination for heatwaves events:
event_freq <- MCS_result %>% 
        mutate(year = lubridate::year(date_start)) %>% 
        group_by(longitude, latitude, year) %>% 
        summarise(n = n(), .groups = "drop")


# create complete grid for merging with:
sst_grid <- OISST_data %>% 
        select(longitude, latitude, t) %>% 
        mutate(t = lubridate::year(t)) %>% 
        dplyr::rename(year = t) %>% 
        distinct()

# and merge:
OISST_n <- left_join(sst_grid, event_freq, by = c("longitude", "latitude", "year")) %>% 
        mutate(n = ifelse(is.na(n), 0, n))

# Function to calculate slope and p-values

lin_fun <- function(ev) {
        mod1 <- glm(n ~ year, family = poisson(link = "log"), data = ev)
        # extract slope coefficient and its p-value
        tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
                         p = summary(mod1)$coefficients[2,4])
        return(tr)
}

# creating trend data 

OISST_nTrend <- OISST_n |> 
        plyr::ddply(c("longitude", "latitude"), lin_fun, .parallel = T) |> 
        mutate(pval=cut(p, breaks = c(0, 0.001, 0.01, 0.05, 1)))


# create a "outputs" folder to store files
dir.create("outputs/")

# saving data 

saveRDS(OISST_nTrend,"outputs/OISST_nTrend.Rds") # Mexican Pacific SST trends
saveRDS(OISST_n,  "outputs/OISST_n.Rds") # Mexican Pacific SST heatwaves n


# Gulf of California Summaries of SST --------------------------------------------

# summarize the number of unique longitude, latitude and year combination for SST mean:
sst_freq <- OISST_data %>% 
        ungroup() %>% 
        mutate(year = lubridate::year(t)) %>% 
        group_by(longitude, latitude, year) %>% 
        summarise(average = mean(temp, na.rm = T))

sst_freq_monthly <- OISST_data %>% 
        ungroup() %>% 
        mutate(year = lubridate::year(t), month = lubridate::month(t)) %>% 
        group_by(longitude, latitude, year, month) %>% 
        summarise(average = mean(temp, na.rm = T))


# create complete grid for merging with:
full_grid <- OISST_data %>% 
        ungroup() %>% 
        select(longitude, latitude, t) %>% 
        mutate(t = lubridate::year(t)) %>% 
        dplyr::rename(year = t) %>% 
        distinct()


month_grid <- OISST_data %>% 
        ungroup() %>% 
        select(longitude, latitude, t) %>% 
        mutate(year = lubridate::year(t), month = lubridate::month(t)) %>% 
        distinct()

full_grid_sf <- sf::st_as_sf(full_grid, coords = c("longitude", "latitude"), crs = 4326) # Goc SST data grid
month_grid_sf <- sf::st_as_sf(month_grid, coords = c("longitude", "latitude"), crs = 4326) # Goc SST data grid

# Intersection to the Gulf of California study area corresponding to rocky reef data
# goc_grid <- sf::st_intersection(full_grid_sf, goc_pol)
goc_idx <- st_intersects(full_grid_sf, goc_pol) |> 
        lengths()>0

goc_grid <- full_grid_sf[goc_idx, ]


# goc_month_grid <- sf::st_intersection(month_grid_sf, goc_pol)
month_idx <- st_intersects(month_grid_sf, goc_pol) |>  
        lengths()>0

goc_month_grid <- month_grid_sf[month_idx, ] 

# this recreate a d.f.
goc_grid_df  <-goc_grid |> 
        mutate(longitude = st_coordinates(geometry)[, "X"],
               latitude = st_coordinates(geometry)[, "Y"]) |> 
        as.data.frame() |> 
        select(-geometry) |> 
        rowid_to_column("id") |> 
        select(year, everything())
        


goc_month_grid_df  <- goc_month_grid |> 
        mutate(longitude = st_coordinates(geometry)[, "X"],
               latitude = st_coordinates(geometry)[, "Y"]) |> 
        as.data.frame() |> 
        select(-geometry) |> 
        rowid_to_column("id") |> 
        select(year, everything())

# grid of sst for the GOC
sst_goc <- left_join(goc_grid_df, sst_freq, by = c("longitude", "latitude", "year")) 
sst_goc_month <- left_join(goc_month_grid_df, sst_freq_monthly, by = c("longitude", "latitude", "year", "month")) 

# grid of MHWs for the GOC
mcs_goc <- left_join(goc_grid_df, event_freq, by = c("longitude", "latitude", "year")) %>% 
        mutate(n = ifelse(is.na(n), 0, n))

# joining sst and heatwaves in one single dataset
all_env_goc <- inner_join(sst_goc, mcs_goc, by = c("longitude", "latitude", "year", "id"))

# Summary by year 
all_env_goc_year <- all_env_goc %>% 
        mutate(Degree=round(latitude,0)) |> 
        group_by(year, Degree) %>% 
        replace(is.na(.), 0) %>% 
        summarise(average = mean(average), 
                  min = min(average), 
                  max = max(average),
                  MCSs = mean(n))
all_env_goc_year <- readRDS("outputs/all_env_goc_year.Rds")
#Heatmap
all_env_goc_year |> 
        mutate(year = as.Date(paste0(year,"-01","-01"))) |> 
        ggplot( aes(x = year, y = Degree, fill = MCSs)) +
        geom_tile() +
        scale_y_continuous(breaks=seq(21, 32, 1)) +
        scale_x_date(date_labels = "%y", date_breaks = "1 year") +
        scale_fill_gradient(high = "darkblue", low = "white", na.value = "white") +
        labs(title = "Marine Cold Spells",
             x = "Year",
             y = "Degree",
             fill = "Frequency") +
        theme_minimal()+
        theme(axis.text.x = element_text(angle=90),
              axis.text = element_text(size=18),
              plot.title = element_text(hjust=0.5, size=19, face="bold"))



#Summary by year and month
sst_goc_month <- sst_goc_month %>% 
        mutate(Degree=round(latitude,0)) |> 
        group_by(year, month, Degree) %>% 
        summarise(average = mean(average)) %>% 
        group_by(month) %>% 
        mutate(anomaly = average - mean(average)) 
     
      

# Saving GOC summary data
saveRDS(all_env_goc, "outputs/all_goc_env.Rds")
saveRDS(all_env_goc_year, "outputs/all_env_goc_year.Rds")
saveRDS(sst_goc_month, "outputs/monthly_env_goc_.Rds")



# End script --------------------------------------------------------------



t |> 
ggplot( aes(x = year, y = Degree, fill = anomaly)) +
        geom_tile() +
        scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
        scale_fill_gradient(low = "darkblue", high = "firebrick", na.value = "white") +
        labs(title = "SST anomalies",
             x = "Year",
             y = "Degree",
             fill = "Anomaly") +
        theme_minimal()






