Data and R code accompanying:

**Effects of bleaching-associated mass coral mortality on reef structural complexity across a gradient of local disturbance**

Authors: Jennifer M.T. Magel, John H.R. Burns, Ruth D. Gates, Julia K. Baum

****

### KI_complexity_data.csv

This file contains the structural complexity, local human disturbance, and benthic composition data for each permanent photoquadrat (PPQ) surveyed in our study. Variables are defined as follows:

* ```Year``` = Year that the data was collected
* ```heat_stress``` = Indicates whether the data was collected before (2015) or after (2017) the 2015-2016 El Niño
* ```Month``` = Month that the data was collected
* ```hum_dist``` = Local human disturbance level, based on fishing pressure and proximity to local villages
* ```2D_Area``` = Two-dimensional area of the PPQ
* ```3D_Area``` = Three-dimensional surface area of the PPQ
* ```rugosity``` = Surface rugosity of the PPQ (3D_Area/2D_Area)
* ```terrain_rug``` = Mean terrain ruggedness of the PPQ
* ```curvature``` = Mean curvature of the PPQ
* ```abs_curv``` = Absolute value of curvature for the PPQ
* ```nbranching``` = Total abundance of branching corals within the PPQ
* ```nplating``` = Total abundance of plating corals within the PPQ
* ```nmassive``` = Total abundance of massive corals within the PPQ
* ```dbranching``` = Density of branching corals within the PPQ (nbranching/2D_Area)
* ```dplating``` = Density of plating corals within the PPQ (nplating/2D_Area)
* ```dmassive``` = Density of massive corals within the PPQ (nmassive/2D_Area)

****

### KI_habitat_volume_data.csv

This file contains parameter settings for the ICP registration and 2.5D volume computation processes in CloudCompare, and resulting data from the 3D point cloud comparison. Variables are defined as follows:

* ```ICP_Overlap``` = Estimated amount of overlap between 2015 and 2017 point clouds (set to 95% to account for slightly different placement of GCPs between years)
* ```ICP_RSL``` = Random sampling limit (i.e. maximum number of sub-sampled points) for the ICP registration
* ```ICP_RMS_Final``` = Square root of mean square errors resulting from the ICP registration
* ```Vol_Step``` = Grid cell size for the 2.5D volume computation (set to 1 cm to match the DEM cell size)
* ```Vol_Cell_Height``` = Specifies how the height of cells containing multiple points should be computed during the 2.5D volume computation
* ```Vol_Match_Cells``` = Percentage of cells present in both the 2015 and 2017 point clouds
* ```Vol_Surface``` = Two-dimensional area of the PPQ
* ```Vol_Change_Total``` = Net change in habitat volume within the PPQ from 2015 to 2017 (Vol_Change_Loss + Vol_Change_Gain)
* ```Vol_Change_Gain``` = Total gain in habitat volume within the PPQ from 2015 to 2017
* Vol_Change_Loss = Total loss in habitat volume within the PPQ from 2015 to 2017
