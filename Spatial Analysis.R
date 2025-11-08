# Authors: Tanishq Tejaswi and Tanya Chhabra

# Reference:
# https://mgimond.github.io/Spatial/spatial-autocorrelation-in-r.html
# https://mkram01.github.io/EPI563-SpatialEPI/spatial-structure-and-clustering-i-morans-i-and-lisa.html

# This code was meant to calculate Moran's index using the cell outline images
# obtained from Cellpose segmentation and the intensity values calculated using 
# ImageJ. Now it is used to also calculate and plot a bunch of other stuff needed
# for our project. Make sure the image you input has only the outlines, preferably
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

# This code will also calculate the mutants location and which mutants are located
# near which clusters. The tagged mutant should be provided in another channel
# in the same file from labels to ROIs plugin, which will contain all the channels
# of the image.

# If a CSV file containing the location and charges of the defects is also provided,
# this code can also calculate the proximity of those defects to these clusters.
# If there is also a mutant channel, it will also calculate the proximity of the
# defects to the mutants.

library(ggplot2);library(spdep); library(terra)
library(tidyterra); library(sf); library(stars)
library(svglite); library(stringr);library(dplyr)

source("Clusters_withNeighbours_and_Defects_Functions.R")
source("Pressures_Functions.R")
source("Mutants_Functions.R")

Main_Folder <- "F:/Sindhu data analysed by Tanishq/New Data by month/2024-07/CC_40_1_dapi_ph594_betacat555_gfp/R2_4"

MeasurementsFolder <- file.path(Main_Folder,"Measurements")
OutlinesFolder <- file.path(Main_Folder,"Bare Outlines")
DefectsFolder <- file.path(Main_Folder,"DefectFiles")
PressuresFolder <- file.path(Main_Folder,"Pressures")

MeasurementFiles_Paths <- file.path(MeasurementsFolder,list.files(MeasurementsFolder))
OutlineFiles_Paths <- file.path(OutlinesFolder,list.files(OutlinesFolder))
DefectFiles_Paths <- file.path(DefectsFolder, list.files(DefectsFolder))
PressureFiles_Paths <- file.path(PressuresFolder, list.files(PressuresFolder))

Main_OutputFolder <- file.path(Main_Folder,"Spatial Plots")
if (!dir.exists(Main_OutputFolder)) dir.create(Main_OutputFolder)
Moran_OutputFolder <- file.path(Main_OutputFolder,"MoranPlots")
if (!dir.exists(Moran_OutputFolder)) dir.create(Moran_OutputFolder)
Shape_OutputFolder <- file.path(Main_OutputFolder,"ShapeMaps")
if (!dir.exists(Shape_OutputFolder)) dir.create(Shape_OutputFolder)
Neighbours_OutputFolder <- file.path(Main_OutputFolder,"Neighbours")
if (!dir.exists(Neighbours_OutputFolder)) dir.create(Neighbours_OutputFolder)
Defects_OutputFolder <- file.path(Main_OutputFolder,"Defects")
if (!dir.exists(Defects_OutputFolder)) dir.create(Defects_OutputFolder)
Pressures_OutputFolder <- file.path(Main_OutputFolder,"Pressures")
if (!dir.exists(Pressures_OutputFolder)) dir.create(Pressures_OutputFolder)

Cell_diameter_list <- read.csv(file.path(Main_Folder,"CellposeOutput","diameters.csv")) #Input average cell diameter here in pixels.
Cell_diameter_list$Image.Name <- tools::file_path_sans_ext(Cell_diameter_list$Image.Name)
ImageNames <- Cell_diameter_list$Image.Name

cutoff_distance_cell <- 3 #Decide the distance cutoff in terms of cell diameters.

Area_cutoff <- 100 #Minimum cutoff-area in pixels below which cells are not considered.
Channel_of_interest <- "2" #"1" #Channel which has the protein of interest
Channel_name <- "Actin" #Name of protein of interest

MutantChannel <- "4" # The reporter/mutant channel. Write "0" if no mutants in the images.
Mutant_cutoff <- 2.5 # Minimum Z-score cutoff for intensities in mutant channel to be considered as mutants.
Mutants_OutputFolder <- file.path(Main_OutputFolder,"Mutants")
if (!dir.exists(Mutants_OutputFolder) & !(MutantChannel == "0")) dir.create(Mutants_OutputFolder)

nsim <- 9999 #Number of Monte-Carlo simulations to get p-values for Moran's analysis
Max_order_correlogram <- 10 #Highest order of neighbours up to which Moran's I is calculated

Scale_pixels_in_microns <- 9.6923 #number of pixels per micron (9.6923 for 63x, 3.0769 for 20x)

Combined_GlobalMoranTable <- data.frame(row.names = ImageNames,MoranI=rep(0,length(ImageNames)),PVal=rep(1,length(ImageNames)))
Combined_CorrelogramTable <- data.frame(row.names = 1:10)
Combined_DefectsTable <- NULL
Combined_Polygon_intensity_data <- NULL
Combined_mutant_neighbour_data <- NULL

#set.seed(1810) #DO NOT USE THIS. MEANT FOR DEBUGGING PURPOSES ONLy.

for (i in ImageNames) {
  
  print(i)
  
  cell_diameter <- Cell_diameter_list$Diameter[Cell_diameter_list$Image.Name == i]
  MeasurementFile <- MeasurementFiles_Paths[grep(paste(i,"_measurements.csv",sep=""),MeasurementFiles_Paths)]
  MeasurementsMain <- Load_Measurements(MeasurementFile = MeasurementFile, Area_cutoff = Area_cutoff)
  
  print(2)
  
  Channel_Measurements <- subset(MeasurementsMain, Ch == Channel_of_interest)
  OutlineFile <- OutlineFiles_Paths[grep(paste(i,"_outlines.png",sep=""),OutlineFiles_Paths)]
  Image_height <- Load_Image_height(OutlineFile)
  Polygon_data <- Load_Outlines(OutlineFile = OutlineFile, Area_cutoff = Area_cutoff)
  Polygon_intensity_data <- Get_combined_measurements_and_polygons(PolyData = Polygon_data,Measurements = Channel_Measurements)
  
  print(3)
  
  PressureFile <- PressureFiles_Paths[grep(paste(i,"_pressures.txt",sep=""),PressureFiles_Paths)]
  Pressures <- Load_pressures(PressureFile)
  Pressures <- na.omit(Pressures)
  Polygon_intensity_data <- Get_combined_pressures_and_polygons(Polygon_intensity_data,Pressures)
  Polygon_intensity_data$neighbours <- Get_neighbours(Polygon_intensity_data)
  
  print(4)
  
  Global_Moran_index <- Get_Global_Moran(Intensity_Poly_Data = Polygon_intensity_data, nsim = nsim, InputFileName = i, OutputFolder = Moran_OutputFolder, Channel_name = Channel_name)
  Combined_GlobalMoranTable[i,]$MoranI <- Global_Moran_index$statistic
  Combined_GlobalMoranTable[i,]$PVal <- Global_Moran_index$p.value
  Combined_CorrelogramTable[,i] <- Correlogram(Intensity_Poly_Data = Polygon_intensity_data,MaxOrder = Max_order_correlogram, Channel_name = Channel_name, InputFileName = i, OutputFolder = Moran_OutputFolder)
  Polygon_intensity_data <- Get_Local_Moran(Intensity_Poly_Data = Polygon_intensity_data, nsim = nsim)

  print(5)
  
  DefectFile <- DefectFiles_Paths[grep(paste(i,"_defects.txt",sep=""),DefectFiles_Paths)]
  Defects <- Load_defects(DefectFile)
  Defects <- Get_cell_location_of_defect(Defects = Defects, Intensity_Polygon_Data = Polygon_intensity_data)
  Defects <- Get_defect_proximity_to_cluster(Defects = Defects, Intensity_Cluster_Data = Polygon_intensity_data, cutoff_distance_cell = cutoff_distance_cell, cell_diameter = cell_diameter)
  Print_defect_proximity(Defects = Defects, InputFileName = i, OutputFolder = Defects_OutputFolder)
  
  print(5)
  
  if (!(MutantChannel == "0")) {
    Polygon_intensity_data <- Identify_mutants(Measurements = MeasurementsMain, Intensity_Poly_Data = Polygon_intensity_data, MutantChannel = MutantChannel, Area_cutoff = Area_cutoff, Mutant_cutoff = Mutant_cutoff )
    Polygon_intensity_data$Comparison[Polygon_intensity_data$Mutant == "Mutant"] <- "Mutant"
    Mutant_only_neighbours <- Get_mutants_and_neighbours(Intensity_Poly_Data = Polygon_intensity_data)
    Mutant_only_neighbours <- Get_mutant_proximity_to_cluster(Mutant_neighbours_data = Mutant_only_neighbours, Intensity_cluster_data = Polygon_intensity_data)
    Print_mutant_proximity(Mutant_neighbours_data = Mutant_only_neighbours, Intensity_Poly_Data = Polygon_intensity_data, InputFileName = i, OutputFolder = Mutants_OutputFolder)
    if (is.null(Combined_mutant_neighbour_data)) {
      Combined_mutant_neighbour_data <- Mutant_only_neighbours
    } else {Combined_mutant_neighbour_data <- rbind.data.frame(Combined_mutant_neighbour_data,Mutant_only_neighbours)}
  } else {Polygon_intensity_data$Mutant <- "Wild-Type"}
  
  Polygon_intensity_data <- Get_defect_proximity_to_cell(Intensity_Poly_Data = Polygon_intensity_data, Defects = Defects, cutoff_distance_cell = cutoff_distance_cell, cell_diameter = cell_diameter)
  Write_cell_defect_proximity(Intensity_Poly_Data = Polygon_intensity_data, InputFileName = i, OutputFolder = Defects_OutputFolder)
  
  print(6)
  
  Polygon_intensity_data <- Get_number_of_neighbours(Intensity_Poly_Data = Polygon_intensity_data, InputFileName = i, OutputFolder = Neighbours_OutputFolder)
  
  print(7)
  Defects <- Get_pressure_around_defects(Defects = Defects, Intensity_poly_data = Polygon_intensity_data)
  Write_pressure_around_defects(Defects = Defects, InputFileName = i, OutputFolder = Pressures_OutputFolder)
  
  print(8)
  
  if (is.null(Combined_Polygon_intensity_data)) {
    Combined_Polygon_intensity_data <- cbind(ImageName = i, Polygon_intensity_data)
  } else {Combined_Polygon_intensity_data <- rbind.data.frame(Combined_Polygon_intensity_data,cbind(ImageName = i, Polygon_intensity_data))}
  
  if (is.null(Combined_DefectsTable)) {
    Combined_DefectsTable <- cbind(ImageName = i, Defects)
  } else {Combined_DefectsTable <- rbind.data.frame(Combined_DefectsTable,cbind(ImageName = i, Defects))}
  
  Cell_diameter_list$Number_of_cells_in_moran[Cell_diameter_list$Image.Name == i] <- nrow(Polygon_intensity_data)
  Cell_diameter_list$Number_of_cells_in_pressures[Cell_diameter_list$Image.Name == i] <- nrow(Pressures)
  Cell_diameter_list$Image_area_in_sq_microns[Cell_diameter_list$Image.Name == i] <- (nrow(rast(OutlineFile))*Image_height)/(Scale_pixels_in_microns)^2
  
  print(i)
}

### Code for combined plots ###

Cell_diameter_list$Diameter_in_microns <- Cell_diameter_list$Diameter/Scale_pixels_in_microns
Cell_diameter_list$Cell_density_per_sq_mm <- (Cell_diameter_list$Number_of_cells_in_moran/Cell_diameter_list$Image_area_in_sq_microns)*10^6
write.csv(Cell_diameter_list, file.path(Main_OutputFolder,"Combined_Cell_Densities.csv"))

Combined_DefectsTable$charge <- NULL
Combined_Polygon_intensity_data$ClosestPointFromCSV <- NULL
Combined_Polygon_intensity_data$geometry <- NULL
Combined_Polygon_intensity_data$ClosestDefect <- NULL

write.csv(Combined_Polygon_intensity_data, file.path(Main_OutputFolder, "Combined_Polygon_intensity_data.csv"))

Write_defect_proximity(Defects = Combined_DefectsTable, InputFileName = "Combined",OutputFolder = Main_OutputFolder)
Plott_pressure_around_defects(Defects = Combined_DefectsTable, InputFileName = "Combined",OutputFolder = Main_OutputFolder)

ggplot(Combined_Polygon_intensity_data, aes(x = as.numeric(ShapeIndex)))+
  geom_vline(xintercept = 3.813)+
  geom_histogram(bins = 30)+
  geom_density(kernel = "gaussian")+
  labs(x = "Shape Index", y = "Count") + theme_classic()
ggsave(file.path(Main_OutputFolder,"Combined Shape Histogram.svg"))

if (!(MutantChannel == "0")) {
  Combined_mutant_neighbour_data$Mutant_IDs <- NULL
  Combined_mutant_neighbour_data$ClosestPointFromCSV <- NULL
  Combined_mutant_neighbour_data$Centroid.x <- NULL
  Combined_mutant_neighbour_data$Centroid.y <- NULL
  Combined_mutant_neighbour_data$Mutant <- NULL
  Combined_mutant_neighbour_data$geometry <- NULL
  Combined_mutant_neighbour_data$Neighbours_List <- NULL
  
  Combined_MutantsProximityTable <- table(Combined_mutant_neighbour_data$Cluster)
  write.table(Combined_MutantsProximityTable,file.path(Main_OutputFolder,"Combined_MutantsProximityTable.txt"))
  
  #Plots a bar-graph of the above table.
  ggplot(data = Combined_mutant_neighbour_data, aes(x=Cluster)) + geom_bar()+ labs(y = "Number of Cells", x = "Near Cluster")+theme_classic()
  ggsave(file.path(Main_OutputFolder,"Combined_MutantsProximity.svg"))
}

dev.off()
