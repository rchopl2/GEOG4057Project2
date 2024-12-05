## This project is creating a script which allows the user to input their workspace, csv file, their output and what ESPG they will be using. 

You can use the script in ArcGIS to create a map which shows the points from you CSV and you can set the symbology to represent different elevations.

Before starting you will want to make sure your Arcpy_Clone is activated and you have Earth Engine set up with a project

## Step 1: Read the .csv file and convert it to a feature class in GIS
#open your VS code

# the following code will allow you to read the csv file
```
import csv
file = open('boundary.csv')
csv_reader = csv.reader(file)
for line in csv_reader:
    print(line)
```

# the following code will describe the data
```
import arcpy
desc = arcpy.da.Describe('flood_2class.tif')
sr =desc['spatialReference']
sr
```
# Now we are going to set the workspace, input & output paths and make sure our workspace has no files with the same name

```
import arcpy
import os

# Set workspace
arcpy.env.workspace = r'C:\Users\rylei\OneDrive\Desktop\Fall 2024\geog4057\Project2'

# Input and output paths
input = os.path.join(arcpy.env.workspace, 'boundary.csv')
output = os.path.join(arcpy.env.workspace, 'boundary_pnts.shp')

# Delete the existing output shapefile if it exists
if arcpy.Exists(output):
    arcpy.Delete_management(output)

# Run the XYTableToPoint tool
arcpy.management.XYTableToPoint(in_table=input, out_feature_class=output, x_field='X', y_field='Y', coordinate_system=32119)
```

# Once that is done you will want to create the geometrys
```
# this code was generated from AI via instructor video
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

# Load the CSV file
csv_file = 'boundary.csv'  # Replace with your CSV path
data = pd.read_csv(csv_file)

# Ensure the column names match your data (replace 'X' and 'Y' with actual column names)
if not {'X', 'Y'}.issubset(data.columns):
    raise ValueError("CSV must contain 'X' and 'Y' columns")

# Create geometry from X and Y columns
geometry = [Point(xy) for xy in zip(data['X'], data['Y'])]

# Create a GeoDataFrame
gdf = gpd.GeoDataFrame(data, geometry=geometry)

# Set the coordinate reference system (CRS)
gdf.set_crs(epsg=32119, inplace=True)  # WGS84 CRS (adjust as needed)

# export to shapefile
output_shapefile = "output_geopandas.shp"
gdf.to_file(output_shapefile, driver = 'ESRI Shapefile')

print(f"Shapefile saved to: {output_shapefile}")
```

## Step 2: Retrive data from 'ee'
```
import ee
ee.Initialize(project='ee-recmonkey')
# Tell the code which DEM from Google Earth Engine you wish to use
dem = ee.Image('USGS/3DEP/10m')
# Create a map in the code..
import geemap
map = geemap.Map()
map

map.addLayer(dem)

geometrys = [ee.Geometry.Point([x,y],'EPSG: 32119') for x,y in zip(data['X'], data['Y'])]

fc = ee.FeatureCollection(geometrys)
map.add_layer(fc)

origin_info = fc.getInfo()
origin_info

sampled_fc = dem.sampleRegions(
    collection=fc,
    scale=10, # Resolution of image
    geometries=True
)

sampled_info = sampled_fc.getInfo()
sampled_info

for ind, itm in enumerate(origin_info['features']):
    itm['properties'] = sampled_info['features'][ind]['properties']

origin_info['features']
```

## Create a feature class and add elevation data to the features

# Creating the feature class
```
import os
fcname = os.path.join(arcpy.env.workspace,'pnt_elev.shp')
if arcpy.Exists(fcname):
    arcpy.management.Delete(fcname)
arcpy.management.CreateFeatureclass(arcpy.env.workspace,'pnt_elev.shp',geometry_type='POINT',spatial_reference=32119)

arcpy.management.AddField(fcname,field_name='elevation',field_type='FLOAT')

with arcpy.da.InsertCursor(fcname,['SHAPE@','elevation']) as cursor:
    for feat in origin_info['features']:
        # get the coordinates and create a pointgeometry
        coords = feat['geometry']['coordinates']
        pnt = arcpy.PointGeometry(arcpy.Point(coords[0],coords[1]),spatial_reference=32119)
        # get the properties and write it to elevation
        elev = feat['properties']['elevation']
        cursor.insertRow([pnt,elev])

```

## Next we will create the python script
# Doing this allows the user to use the script so it isn't hard coded and they can add their own file names... the script is below
```
# Import all used packages
import arcpy
import os
import ee
import pandas as pd
import geopandas as gpd

"""
example usage:
python project2.py "C:\Users\rylei\OneDrive\Desktop\Fall 2024\geog4057\Project2" boundary.csv pnt_elev2.shp 32119
"""
def GetGeeElevation(workspace,csv_file,outfc_name,epsg=4326):

    """
    workspace: directory that contains input and output
    csv_file: input csv filename
    epsg: wkid code for the spatial reference, e.g. 4326 for WGS GCS
    """

    # Load the CSV file
    csv_file = os.path.join(workspace,csv_file)
    data = pd.read_csv(csv_file)
    dem = ee.Image('USGS/3DEP/10m')
    print(epsg)
    geometrys = [ee.Geometry.Point([x,y],f'EPSG:{epsg}') for x,y in zip(data['X'], data['Y'])]
    fc = ee.FeatureCollection(geometrys)
    origin_info = fc.getInfo()
    sampled_fc = dem.sampleRegions(
        collection=fc,
        scale=10, # Resolution of image
        geometries=True
    )
    sampled_info = sampled_fc.getInfo()
    for ind, itm in enumerate(origin_info['features']):
        itm['properties'] = sampled_info['features'][ind]['properties']



    fcname = os.path.join(workspace,outfc_name)
    if arcpy.Exists(fcname):
        arcpy.management.Delete(fcname)
    arcpy.management.CreateFeatureclass(workspace,outfc_name,geometry_type='POINT',spatial_reference=epsg)
    
    arcpy.management.AddField(fcname,field_name='elevation',field_type='FLOAT')
    
    with arcpy.da.InsertCursor(fcname,['SHAPE@','elevation']) as cursor:
        for feat in origin_info['features']:
            # get the coordinates and create a pointgeometry
            coords = feat['geometry']['coordinates']
            pnt = arcpy.PointGeometry(arcpy.Point(coords[0],coords[1]),spatial_reference=32119)
            # get the properties and write it to elevation
            elev = feat['properties']['elevation']
            cursor.insertRow([pnt,elev])
    
    
def main():
    import sys
    try:
        ee.Initialize(project='ee-recmonkey')
    except:
        ee.Authenticate()
        ee.Initialize(project='ee-recmonkey')
    workspace = sys.argv[1]
    csv_file = sys.argv[2]
    outfc_name = sys.argv[3]
    epsg = int(sys.argv[4])
    
    GetGeeElevation(workspace=workspace,csv_file=csv_file, outfc_name=outfc_name, epsg=epsg)
    



if __name__== '__main__':
    main()
```

## In the final step we will add the following code into the ArcGIS Script tool, found in your project's toolbox
```
"""
Script documentation

- Tool parameters are accessed using arcpy.GetParameter() or 
                                     arcpy.GetParameterAsText()
- Update derived parameter values using arcpy.SetParameter() or
                                        arcpy.SetParameterAsText()
"""
import arcpy
import os
import ee
import pandas as pd
import geopandas as gpd


def GetGeeElevation(workspace,csv_file,outfc_name,epsg=4326):

    ee.Initialize(project='ee-recmonkey')

    """
    workspace: directory that contains input and output
    csv_file: input csv filename
    epsg: wkid code for the spatial reference, e.g. 4326 for WGS GCS
    """

    # Load the CSV file
    csv_file = os.path.join(workspace,csv_file)
    data = pd.read_csv(csv_file)
    dem = ee.Image('USGS/3DEP/10m')
    print(epsg)
    geometrys = [ee.Geometry.Point([x,y],f'EPSG:{epsg}') for x,y in zip(data['X'], data['Y'])]
    fc = ee.FeatureCollection(geometrys)
    origin_info = fc.getInfo()
    sampled_fc = dem.sampleRegions(
        collection=fc,
        scale=10, # Resolution of image
        geometries=True
    )
    sampled_info = sampled_fc.getInfo()
    for ind, itm in enumerate(origin_info['features']):
        itm['properties'] = sampled_info['features'][ind]['properties']



    fcname = os.path.join(workspace,outfc_name)
    if arcpy.Exists(fcname):
        arcpy.management.Delete(fcname)
    arcpy.management.CreateFeatureclass(workspace,outfc_name,geometry_type='POINT',spatial_reference=epsg)
    
    arcpy.management.AddField(fcname,field_name='elevation',field_type='FLOAT')
    
    with arcpy.da.InsertCursor(fcname,['SHAPE@','elevation']) as cursor:
        for feat in origin_info['features']:
            # get the coordinates and create a pointgeometry
            coords = feat['geometry']['coordinates']
            pnt = arcpy.PointGeometry(arcpy.Point(coords[0],coords[1]),spatial_reference=32119)
            # get the properties and write it to elevation
            elev = feat['properties']['elevation']
            cursor.insertRow([pnt,elev])

if __name__ == "__main__":

    workspace = arcpy.GetParameterAsText(0)
    csv_file = arcpy.GetParameterAsText(1)
    outfc_name = arcpy.GetParameterAsText(2)
    epsg = int(arcpy.GetParameterAsText(3))

    GetGeeElevation(workspace,csv_file,outfc_name,epsg)
```
