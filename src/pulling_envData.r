#'##################################################################################################'
#'                                                                                                 #'
#' Project........ Environmental GWAS (eGWAS) for PandAnd -- Hackathon                             #'
#' Title.......... Environmental Data extraction for each occurence record                         #'
#' Created at..... 11-06-2023                                                                      #'
#' Updated at..... 04-30-2024 (Sheng-Kai Hsu)                                                      #'
#' Author: G.Costa-Neto, germano.cneto<at>gmail.com                                                #'
#'                                                                                                 #'
#'##################################################################################################'

require(tidyverse)
require(plyr)
require(reshape2)
require(terra)


#'------------------------------------------------------------------------------------------------------------
# (1) load geo data 
#'------------------------------------------------------------------------------------------------------------
data_clean = data.table::fread('/workdir/hdd29/chloroplast_genome_evaluation/data/combined_clean.csv')
#data_clean <- data_clean[1:5000, ]

head(data_clean)
dim(data_clean)

#'------------------------------------------------------------------------------------------------------------
# (2) extract env. features for each coordinate
#'------------------------------------------------------------------------------------------------------------
# creating a fake "environmental unit": species - sample
data_clean <-
  data_clean %>% 
  ddply(.(LUI),mutate,envScientificName = paste0('env_',LUI,'_',1:length(decimalLatitude)))
dim(data_clean)
head(data_clean)


data_clean <- data_clean %>% na.omit()

# table(data_clean$scientificName) %>% sort()
# table(data_clean$scientificName) %>% hist(breaks = seq(0,3500,10))

########### Bioclimatic variables
# check src_generating_FAO_GAEZ.R to see how to generate enviromeDB::WC_Bioclimate since the package is broken
source('https://raw.githubusercontent.com/gcostaneto/envirotypeR/main/R/get_spatial_fun.R')

#url = '/workdir/hdd29/env_data/data/WC_Bioclim.rds'
#tmp = readRDS(url)
#digital_raster <- terra::rast(digital_raster)
#rast_rds <- readRDS(url)
#names(rast_rds) #"bio08_Mean_Temperature_Wettest_Quarter" is the only layer I need
#desired_layer <- rast_rds[["bio08_Mean_Temperature_Wettest_Quarter"]]
#print(desired_layer)
#terra::ncell(desired_layer)  # Total number of cells
#terra::ext(desired_layer)   # Extent

# Path to your folder of .tif files
tif_dir <- "/local/workdir/hdd29/chloroplast_genome_evaluation/data/geoDataTifs"

# Make sure output base directory exists
base_out <- "data/geoDataOut"
dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

# List all the .tif files
tif_files <- list.files(tif_dir, pattern = "\\.tif$", full.names = TRUE)

# Loop through each file
for (tif_path in tif_files) {
  # load raster
  layer <- rast(tif_path)
  
  # optionally inspect
  message("Processing ", basename(tif_path), 
          ": ", ncell(layer), " cells; extent ", ext(layer))
  
  # extract spatial data
  geo_df <- get_spatial(
    env.dataframe  = data_clean,
    lat            = 'decimalLatitude',
    lng            = 'decimalLongitude',
    env.id         = 'envScientificName',
    digital.raster = layer
  )
  
  # clean up infinite or sentinel values
  geo_df[geo_df == Inf   ] <- NA
  geo_df[geo_df == -Inf  ] <- NA
  geo_df[geo_df == -99   ] <- NA
  geo_df[geo_df == -999  ] <- NA
  
  # create a per-layer output directory
  layer_name <- tools::file_path_sans_ext(basename(tif_path))
  out_dir    <- file.path(base_out, layer_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  #keep only the LUI and the climate variable, this is all we need
  geo_df <- geo_df %>%
    select(LUI, layer_name, decimalLatitude, decimalLongitude)
  
  # write out cleaned CSV
  out_file <- file.path(out_dir, paste0(layer_name, "_envData.csv"))
  write.csv(geo_df,
            file      = out_file,
            row.names = FALSE,
            quote     = FALSE)
  
  # print first few rows for sanity check
  print(head(geo_df))
}
saveRDS(geo_df, "data/topt_df.RDS")
#get the soil nitrogen data

########### Soil Features from GSDE
# the file needs to be converted (reading as brick)
# and this is a too big .rds file.
# so let's pull each layer per time
urlList = list.files('/workdir/hdd29/chloroplast_genome_evaluation/data/GIS_env_data/GSDE_raw_nc_files/',pattern = "*.nc",recursive = T,full.names = T)
print("getting soil data")
for(i in 1:length(urlList))
{
  # this takes a long time in relation to the previous raster files. Don't worry.
  soil_df = 
    get_spatial( env.dataframe = data_clean,
                 lat = 'decimalLatitude',
                 lng = 'decimalLongitude',
                 env.id = 'envScientificName',#which.raster.number = 1,
                 digital.raster = raster::brick(urlList[i]))#
}
print("got soil data")
# take the first layer only

soil_df[soil_df==Inf] = NA
soil_df[soil_df==-Inf] = NA
soil_df[soil_df==-99] = NA
soil_df[soil_df==-999] = NA
noNAIdx = apply(soil_df,2,function(x) !any(is.na(x)))
soil_df[,noNAIdx][soil_df[,noNAIdx]==156] = NA # for % data in GSDE, NA -> 156...

head(soil_df)
saveRDS(soil_df, "data/soil_df.RDS")

soil_df <- readRDS("data/soil_df.RDS")
topt_df <- readRDS("data/topt_df.RDS")

library(data.table)
topt_df$Topt_site <- topt_df$Topt_site * 0.01
soil_df$X4.5 <- soil_df$X4.5 * 0.01
soil_df$X9.10000038146973 <- soil_df$X9.10000038146973 * 0.01
soil_df$X16.6000003814697 <- soil_df$X16.6000003814697 * 0.01

setDT(topt_df) 
summary(topt_df)

setDT(soil_df) 
summary(soil_df)

sum(is.na(topt_df$Topt_site))
sum(is.na(soil_df$X4.5))
sum(is.na(soil_df$X9.10000038146973))
sum(is.na(soil_df$X16.6000003814697))

soil_df <- na.omit(soil_df)
topt_df <- na.omit(topt_df)

setkey(soil_df, LUI, decimalLatitude, decimalLongitude)
setkey(topt_df, LUI, decimalLatitude, decimalLongitude)

# Perform inner join
merged_df <- soil_df[topt_df, nomatch = 0]

#for each LUI, get the median T_Opt_site and median X4.5, and count the number of records
results <- merged_df[, .(
  Median_Topt_site = median(Topt_site, na.rm = TRUE),
  Median_X4.5 = median(X4.5, na.rm = TRUE),
  n = .N
), by = LUI]

# Save the results to a CSV file
write.csv(results, "data/geoDataOut/merged_results.csv", row.names = FALSE)

library(ggplot2)

p <- ggplot(merged_df, aes(x = Topt_site, y = X4.5)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "plasma", trans = "log", name = "Count") +
  labs(
    title = "Topt site vs TN1 @ min depth",
    x = "Topt_site (°C)",
    y = "% Soil N at 4.5cm depth"
  ) +
  theme_minimal(base_size = 14)

ggsave("topt_vs_soilN.pdf", plot = p, width = 5, height = 4, units = "in", dpi = 300)

median(merged_df$X4.5)

library(data.table)
merged_df$LUI <- as.factor(merged_df$LUI)

results <- results[n > 30]

results <- merged_df[!is.na(X4.5) & !is.na(Topt_site),
                     {
                       if (length(unique(X4.5)) > 1 && length(unique(Topt_site)) > 1) {
                         m <- tryCatch(lm(Topt_site ~ X4.5), error = function(e) NULL)
                         if (!is.null(m)) {
                           coef_m <- coef(m)
                           if ("X4.5" %in% names(coef_m)) {
                             pval <- summary(m)$coefficients["X4.5", "Pr(>|t|)"]
                             list(
                               slope = as.numeric(coef_m["X4.5"]),
                               pval = as.numeric(pval),
                               n = .N
                             )
                           } else {
                             list(slope = as.numeric(NA), pval = as.numeric(NA), n = .N)
                           }
                         } else {
                           list(slope = as.numeric(NA), pval = as.numeric(NA), n = .N)
                         }
                       } else {
                         list(slope = as.numeric(NA), pval = as.numeric(NA), n = .N)
                       }
                     },
                     by = LUI]



head(results)
results <- results[n > 10]
hist(results$slope)
plot(results$slope, -log10(results$pval))
ggplot(results, aes(x = slope)) + geom_histogram()

#plot(merged_df$Topt_site, merged_df$X4.5, main="Topt site vs TN1 @ min depth")

library(ggplot2)

ggplot(merged_df[sample(.N, 100000)], aes(x = Topt_site, y = X4.5)) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
  scale_fill_viridis_c() +
  labs(title = "Topt site vs TN1 @ min depth", x = "Topt_site", y = "X4.5") +
  theme_minimal()

library(ggrastr)  # for rasterized points (fast and memory-efficient)

ggplot(merged_df, aes(x = Topt_site, y = X4.5)) +
  rasterize(geom_point(alpha = 0.01, shape = 16, color = "darkblue"), dpi = 150) +
  labs(title = "Topt site vs TN1 @ min depth", x = "Topt_site", y = "X4.5") +
  theme_minimal(base_size = 14)

ggplot(merged_df, aes(x = Topt_site, y = X4.5)) +
  geom_hex(bins = 100) +  # hex bins for density
  scale_fill_viridis_c(option = "plasma", trans = "log", name = "Count") +  # better scale
  labs(
    title = "Topt site vs TN1 @ min depth",
    x = "Topt_site (°C)",
    y = "Soil N (X4.5)"
  ) +
  theme_minimal(base_size = 14)

# put NA on Weird results
# .ControlData <- function(x)
# {
#   if(isTRUE(x ==  Inf) ) x <- NA
#   if(isTRUE(is.nan(x)) ) x <- NA
#   if(isTRUE(x == -99)  ) x <- NA
#   if(isTRUE(x == -999)) x <- NA
#   if(isTRUE(x == -Inf) ) x <- NA
#   return(x)
# }


write.table(geographic_ranges_bien_filtered,
            "/workdir/sh2246/p_phyloGWAS/output/metadataFormalOut/formal_envData_20240820.txt",quote = F,sep = "\t")



quit(save="not", status=0)


bio8_tif_path <- "/workdir/hdd29/chloroThermometer/data/GIS_env_data/WorldClim_raw_2.5m_files/wc2.1_2.5m_bio_8.tif"
desired_layer <- terra::rast(bio8_tif_path)
terra::ncell(desired_layer)  # Total number of cells
terra::ext(desired_layer) 

geographic_ranges_bien  = 
  get_spatial( env.dataframe =data_clean,
               lat = 'decimalLatitude',
               lng = 'decimalLongitude',
               env.id = 'envScientificName',
               digital.raster = desired_layer)

geographic_ranges_bien_filtered <- geographic_ranges_bien
geographic_ranges_bien_filtered[geographic_ranges_bien_filtered==Inf] = NA
geographic_ranges_bien_filtered[geographic_ranges_bien_filtered==-Inf] = NA
geographic_ranges_bien_filtered[geographic_ranges_bien_filtered==-99] = NA
geographic_ranges_bien_filtered[geographic_ranges_bien_filtered==-999] = NA

head(geographic_ranges_bien)
write.csv(geographic_ranges_bien_filtered,
            "data/geoData/envData.csv", 
            row.names = FALSE, 
            quote = FALSE)


length(unique(geographic_ranges_bien_filtered$LUI))
# Select only the scientificName and the last column
envData_simple <- geographic_ranges_bien_filtered[, c("scientificName", "wc2.1_2.5m_bio_8")]

hist(envData_simple$wc2.1_2.5m_bio_8
# Save the simplified data to a CSV file
write.csv(envData_simple, 
          "data/geoData/envData_simple.csv", 
          row.names = FALSE, 
          quote = FALSE)

write.csv(geographic_ranges_bien_filtered,
            "data/geoData/envData_simple.csv", 
            row.names = FALSE, 
            quote = FALSE)

hist(geographic_ranges_bien_filtered$wc2.1_2.5m_bio_8)
head(geographic_ranges_bien_filtered)
########### Elevation
url = '/workdir/sh2246/p_evolBNI/data/GIS_env_data/WorldClim_raw_2.5m_files/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif' 
geographic_ranges_bien = 
  get_spatial( env.dataframe = geographic_ranges_bien,
               lat = 'decimalLatitude',
               lng = 'decimalLongitude',
               env.id = 'envScientificName',
               name.feature = 'Elevation_m',
               digital.raster = terra::rast(url))#envirotypeR::SRTM_elevation #terra::rast(url), # using a certain url)

########### Global Hydrologic Soil Groups
url = '/workdir/sh2246/p_evolBNI/data/GIS_env_data/Global_Hydrologic_Soil_Group_1566/Global_Hydrologic_Soil_Group_1566/data/HYSOGs250m.tif'
geographic_ranges_bien = 
  get_spatial( env.dataframe = geographic_ranges_bien,
               lat = 'decimalLatitude',
               lng = 'decimalLongitude',
               env.id = 'envScientificName',
               name.feature = 'HYSOGs',
               digital.raster = terra::rast(url) #terra::rast(url), # using a certain url)
  ) 


########### FAO-GAEZ 
url = '/workdir/sh2246/p_phyloGWAS/output/envData/GIS_raster/GAEZ_AEZ.rds'
geographic_ranges_bien = 
  get_spatial( env.dataframe = geographic_ranges_bien,
               lat = 'decimalLatitude',
               lng = 'decimalLongitude',
               env.id = 'envScientificName',
               digital.raster = readRDS(url))#

########### Soil Temperature 
url = '/workdir/sh2246/p_phyloGWAS/output/envData/GIS_raster/TEMP_soil.rds'
geographic_ranges_bien = 
  get_spatial( env.dataframe = geographic_ranges_bien,
               lat = 'decimalLatitude',
               lng = 'decimalLongitude',
               env.id = 'envScientificName',
               digital.raster = readRDS(url))#


########### Soil Features from GSDE
# the file needs to be converted (reading as brick)
# and this is a too big .rds file.
# so let's pull each layer per time
urlList = list.files('/workdir/sh2246/p_evolBNI/data/GIS_env_data/GSDE_raw_nc_files/',pattern = "*.nc",recursive = T,full.names = T)


for(i in 1:length(urlList))
{
  # this takes a long time in relation to the previous raster files. Don't worry.
  geographic_ranges_bien = 
    get_spatial( env.dataframe = geographic_ranges_bien,
                 lat = 'decimalLatitude',
                 lng = 'decimalLongitude',
                 env.id = 'envScientificName',#which.raster.number = 1,
                 digital.raster = raster::brick(urlList[i]))#
}

# name the GSDE variablees
GSDE_varNames = limma::strsplit2(list.dirs("/workdir/sh2246/p_evolBNI/data/GIS_env_data/GSDE_raw_nc_files",recursive = F),"/")[,8]

colnames(geographic_ranges_bien)[-c(1:71)] = paste(rep(GSDE_varNames,each = 4),c(5,15,30,200),sep = "_GSDE_")

# take the first layer only
rmIdx = grep('GSDE_15|_30|_200',colnames(geographic_ranges_bien))
geographic_ranges_bien_filtered = geographic_ranges_bien[,-rmIdx]

geographic_ranges_bien_filtered = geographic_ranges_bien_filtered[,c(6,4,5,7:104)]

geographic_ranges_bien_filtered[geographic_ranges_bien_filtered==Inf] = NA
geographic_ranges_bien_filtered[geographic_ranges_bien_filtered==-Inf] = NA
geographic_ranges_bien_filtered[geographic_ranges_bien_filtered==-99] = NA
geographic_ranges_bien_filtered[geographic_ranges_bien_filtered==-999] = NA
noNAIdx = apply(geographic_ranges_bien_filtered,2,function(x) !any(is.na(x)))
geographic_ranges_bien_filtered[,noNAIdx][geographic_ranges_bien_filtered[,noNAIdx]==156] = NA # for % data in GSDE, NA -> 156...


# put NA on Weird results
# .ControlData <- function(x)
# {
#   if(isTRUE(x ==  Inf) ) x <- NA
#   if(isTRUE(is.nan(x)) ) x <- NA
#   if(isTRUE(x == -99)  ) x <- NA
#   if(isTRUE(x == -999)) x <- NA
#   if(isTRUE(x == -Inf) ) x <- NA
#   return(x)
# }


write.table(geographic_ranges_bien_filtered,
            "/workdir/sh2246/p_phyloGWAS/output/metadataFormalOut/formal_envData_20240820.txt",quote = F,sep = "\t")