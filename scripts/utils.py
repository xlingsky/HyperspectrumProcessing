from osgeo import gdal

def combine_bands(input_files, output_file, raster_driver = "GTiff", dtype=None):
    """
    Combine multiple image files along the band dimension using GDAL
    
    Parameters:
    input_files -- List of input file paths (can be single or multi-band)
    output_file -- Path for output multi-band file
    dtype -- GDAL data type for output (default: Float32)
    
    Returns:
    None (writes output to file)
    """
    gdal.UseExceptions()
    # Open first file to get spatial reference
    ds = gdal.Open(input_files[0])
    if ds is None:
        raise ValueError(f"Could not open file: {input_files[0]}")
    
    # Get spatial parameters from first file
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    geotransform = ds.GetGeoTransform()
    projection = ds.GetProjection()
    if dtype is None:
        dtype = ds.GetRasterBand(1).DataType
    ds = None  # Close the dataset
    
    # Count total bands across all files
    total_bands = 0
    for file in input_files:
        temp_ds = gdal.Open(file)
        if temp_ds is None:
            raise ValueError(f"Could not open file: {file}")
        total_bands += temp_ds.RasterCount
        temp_ds = None
    
    # Create output file
    driver = gdal.GetDriverByName(raster_driver)
    out_ds = driver.Create(output_file, cols, rows, total_bands, dtype)
    out_ds.SetGeoTransform(geotransform)
    out_ds.SetProjection(projection)
    
    # Combine bands
    current_band = 1
    for file in input_files:
        temp_ds = gdal.Open(file)
        for band_num in range(1, temp_ds.RasterCount + 1):
            src_band = temp_ds.GetRasterBand(band_num)
            out_band = out_ds.GetRasterBand(current_band)
            
            # Read data and write to output
            data = src_band.ReadAsArray()
            try:
                out_band.WriteArray(data)
            except Exception as e:
                print(f"Error writing band {current_band} of {output_file}: {e}")
            
            current_band += 1
        
        temp_ds = None
    
    # Flush cache and close files
    out_ds.FlushCache()
    out_ds = None
    print(f"Successfully created {output_file} with {total_bands} bands")

if __name__ == "__main__":
    import sys
    
    with open(sys.argv[1], 'r') as f:
        input_files = [line.strip() for line in f.readlines()]
    output_file = sys.argv[2]

    combine_bands(input_files, output_file)