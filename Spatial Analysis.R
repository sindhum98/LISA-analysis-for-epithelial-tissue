# Authors: Tanishq Tejaswi and Tanya Chhabra

# Reference:
# https://mgimond.github.io/Spatial/spatial-autocorrelation-in-r.html
# https://mkram01.github.io/EPI563-SpatialEPI/spatial-structure-and-clustering-i-morans-i-and-lisa.html

# This code is meant to calculate Moran's index using the cell outline images
# obtained from Cellpose segmentation and the intensity values calculated using 
# ImageJ. Make sure the image you input has only the outlines, preferably
# on ablack background.

# The measurements file needs to have 3 columns marked X, Y and Mean (containing)
# the mean intensity values of the protein. Also include the area and perimeter
# of the cells if you wish to shape index, protein levels,
# forces etc. Make sure all the outliers such as anomalously segmented regions with
# minute areas etc have been removed already. The coordinate values should be in
# pixels, NOT microns.

# The output will give several plots: A global Moran's plot indicating the the correlation
# between the spatially lagged intensity of a cell (mean intensity of the neighbours)
# against the cell's intensity. Then, the correlogram will show how the Moran's Index 
# decreases as higher and higher order neighbours are considered, i.e., how it
# decreases with distance in terms of average cell diameters. There will be a
# heat map of the intensity values themselves. The local Moran's plot will give 
# the value of significant (p < 5%) local Moran's indices for each cell. The LISA
# plot will give the clusters based on the local Moran's I calculated above. Heat
# maps indicating the shape index of the cells and whether the cells are jammed
# or not will also be plotted. Note that all the plots will be flipped vertically
# (about X-axis) as the origin here is on the bottom-right, not top-right.


library(ggplot2);library(spdep); library(terra)
library(tidyterra); library(sf); library(stars)
library(svglite); library(stringr);library(dplyr)

source("E:/Sindhu/codes/Actin_segmentation_MoransOnly.R")

Main_Folder <- "E:/Sindhu/2025_analysis/Model_LISA"

MeasurementsFolder <- file.path(Main_Folder,"Measurements")
OutlinesFolder <- file.path(Main_Folder,"Bare Outlines")

MeasurementFiles_Paths <- file.path(MeasurementsFolder,list.files(MeasurementsFolder))
OutlineFiles_Paths <- file.path(OutlinesFolder,list.files(OutlinesFolder))

Main_OutputFolder <- file.path(Main_Folder,"Spatial Plots")
if (!dir.exists(Main_OutputFolder)) dir.create(Main_OutputFolder)
Moran_OutputFolder <- file.path(Main_OutputFolder,"MoranPlots")
if (!dir.exists(Moran_OutputFolder)) dir.create(Moran_OutputFolder)
Shape_OutputFolder <- file.path(Main_OutputFolder,"ShapeMaps")
if (!dir.exists(Shape_OutputFolder)) dir.create(Shape_OutputFolder)
Neighbours_OutputFolder <- file.path(Main_OutputFolder,"Neighbours")
if (!dir.exists(Neighbours_OutputFolder)) dir.create(Neighbours_OutputFolder)

input_folder <- file.path(Main_Folder, "Images")
output_folder <- file.path(Main_Folder, "CellposeOutput")

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

csv_file <- file.path(output_folder, "diameters.csv")

# Get all .tif file names from Images folder
ImageNames <- list.files(input_folder, pattern = "\\.tif$", full.names = FALSE)
ImageNames <- tools::file_path_sans_ext(ImageNames)
#ImageNames <- sub("_label$", "", ImageNames)

# Create a dataframe with placeholder diameter values
results <- data.frame(
  `Image Name` = ImageNames,
  Diameter = rep(5000, length(ImageNames)),  # replace with actual values later
  stringsAsFactors = FALSE
)

# Write CSV
write.csv(results, csv_file, row.names = FALSE)

message("CSV saved at: ", csv_file)
Cell_diameter_list <- read.csv(file.path(Main_Folder,"CellposeOutput","diameters.csv")) #Input average cell diameter here in pixels.
#Cell_diameter_list$Image.Name <- tools::file_path_sans_ext(Cell_diameter_list$Image.Name)
#Cell_diameter_list$Image.Name <- sub("_label$", "", Cell_diameter_list$Image.Name)
ImageNames <- Cell_diameter_list$Image.Name

cutoff_distance_cell <- 3 #Decide the distance cutoff in terms of cell diameters.

Area_cutoff <- 0 #Minimum cutoff-area in pixels below which cells are not considered.
Channel_of_interest <- "1" #"1" #Channel which has the protein of interest
Channel_name <- "Actin" #Name of protein of interest

nsim <- 9999 #Number of Monte-Carlo simulations to get p-values for Moran's analysis
Max_order_correlogram <- 6 #Highest order of neighbours up to which Moran's I is calculated

Scale_pixels_in_microns <- 1 #number of pixels per micron (9.6923 for 63x, 3.0769 for 20x)

Combined_GlobalMoranTable <- data.frame(row.names = ImageNames,MoranI=rep(0,length(ImageNames)),PVal=rep(1,length(ImageNames)))
Combined_CorrelogramTable <- data.frame(row.names = 1:Max_order_correlogram)
Combined_Polygon_intensity_data <- NULL

for (i in ImageNames) {
  
  #cell_diameter <- Cell_diameter_list$Diameter[Cell_diameter_list$Image.Name == i]
  MeasurementFile <- MeasurementFiles_Paths[grep(paste(i,"_measurements.csv",sep=""),MeasurementFiles_Paths)]
  #print(MeasurementFile)
  MeasurementsMain <- Load_Measurements(MeasurementFile = MeasurementFile, Area_cutoff = Area_cutoff)
  print(MeasurementsMain)
  
  #Channel_Measurements <- subset(MeasurementsMain, Ch == Channel_of_interest)
  Channel_Measurements <- MeasurementsMain
  OutlineFile <- OutlineFiles_Paths[grep(paste(i,"_outlines.png",sep=""),OutlineFiles_Paths)]
  #print(OutlineFiles_Paths)
  Image_height <- Load_Image_height(OutlineFile)
  Polygon_data <- Load_Outlines(OutlineFile = OutlineFile, Area_cutoff = Area_cutoff)
  Polygon_intensity_data <- Get_combined_measurements_and_polygons(PolyData = Polygon_data,Measurements = Channel_Measurements)
  
  Polygon_intensity_data$neighbours <- Get_neighbours(Polygon_intensity_data)
  
  Global_Moran_index <- Get_Global_Moran(Intensity_Poly_Data = Polygon_intensity_data, nsim = nsim, InputFileName = i, OutputFolder = Moran_OutputFolder, Channel_name = Channel_name)
  Combined_GlobalMoranTable[i,]$MoranI <- Global_Moran_index$statistic
  Combined_GlobalMoranTable[i,]$PVal <- Global_Moran_index$p.value
  Combined_CorrelogramTable[,i] <- Correlogram(Intensity_Poly_Data = Polygon_intensity_data,MaxOrder = Max_order_correlogram, Channel_name = Channel_name, InputFileName = i, OutputFolder = Moran_OutputFolder)
  Polygon_intensity_data <- Get_Local_Moran(Intensity_Poly_Data = Polygon_intensity_data, nsim = nsim)
  
  Polygon_intensity_data <- Get_number_of_neighbours(Intensity_Poly_Data = Polygon_intensity_data, InputFileName = i, OutputFolder = Neighbours_OutputFolder)
  Write_cell_data(Intensity_Poly_Data = Polygon_intensity_data, InputFileName = i, OutputFolder = Main_OutputFolder)
  
  if (is.null(Combined_Polygon_intensity_data)) {
    Combined_Polygon_intensity_data <- cbind(ImageName = i, Polygon_intensity_data)
  } else {Combined_Polygon_intensity_data <- rbind.data.frame(Combined_Polygon_intensity_data,cbind(ImageName = i, Polygon_intensity_data))}
  
  Cell_diameter_list$Number_of_cells_in_moran[Cell_diameter_list$Image.Name == i] <- nrow(Polygon_intensity_data)
  Cell_diameter_list$Image_area_in_sq_microns[Cell_diameter_list$Image.Name == i] <- (nrow(rast(OutlineFile))*Image_height)/(Scale_pixels_in_microns)^2
  
  print(i)
}

### Code for combined plots ###

Cell_diameter_list$Diameter_in_microns <- Cell_diameter_list$Diameter/Scale_pixels_in_microns
Cell_diameter_list$Cell_density_per_sq_mm <- (Cell_diameter_list$Number_of_cells_in_moran/Cell_diameter_list$Image_area_in_sq_microns)*10^6
write.csv(Cell_diameter_list, file.path(Main_OutputFolder,"Combined_Cell_Densities.csv"))

#Combined_DefectsTable$charge <- NULL
Combined_Polygon_intensity_data$ClosestPointFromCSV <- NULL
Combined_Polygon_intensity_data$geometry <- NULL
Combined_Polygon_intensity_data$neighbours <- NULL

write.csv(Combined_Polygon_intensity_data, file.path(Main_OutputFolder, "Combined_Polygon_intensity_data.csv"))

write.csv(Combined_GlobalMoranTable, file.path(Main_OutputFolder, "Combined_GlobalMoran_data.csv"))
write.csv(Combined_CorrelogramTable, file.path(Main_OutputFolder, "Combined_Correlogram_data.csv"))

dev.off()