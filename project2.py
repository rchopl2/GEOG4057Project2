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


