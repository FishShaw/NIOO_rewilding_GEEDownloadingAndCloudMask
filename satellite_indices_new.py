import ee
import geemap
import os
import datetime

# --- ROI Helper Function ---

def get_roi_from_shapefile(shapefile_path):
    """Loads ROI geometry from a shapefile."""
    try:
        # Attempt to convert shapefile to ee.FeatureCollection
        roi_fc = geemap.shp_to_ee(shapefile_path)
        # Check if the result is a non-empty FeatureCollection
        if isinstance(roi_fc, ee.FeatureCollection) and roi_fc.size().getInfo() > 0:
            # Return the geometry of the FeatureCollection
            return roi_fc.geometry()
        # Check if the result is directly a Geometry
        elif isinstance(roi_fc, ee.Geometry):
            return roi_fc
        else:
            print(f"Error: Could not identify a valid ee.Geometry or non-empty ee.FeatureCollection from {shapefile_path}.")
            return None
    except ee.EEException as e:
        print(f"Error: GEE error accessing shapefile ({shapefile_path}): {e}")
        return None
    except Exception as e:
        print(f"Error: Unexpected error loading ROI from shapefile ({shapefile_path}): {e}")
        return None

# --- Landsat Data Processing ---

def _mask_clouds_shadows_landsat(image):
    """Internal function: Masks clouds, shadows, and snow in Landsat C2 L2."""
    qa = image.select('QA_PIXEL')
    # Bits 3: Cloud, 4: Cloud Shadow, 5: Snow
    cloud_mask = qa.bitwiseAnd(1 << 3).eq(0)      # Not Cloud
    shadow_mask = qa.bitwiseAnd(1 << 4).eq(0)     # Not Cloud Shadow
    snow_mask = qa.bitwiseAnd(1 << 5).eq(0)       # Not Snow
    final_mask = cloud_mask.And(shadow_mask).And(snow_mask)
    # Apply mask and copy time properties
    return image.updateMask(final_mask).copyProperties(image, ['system:time_start'])

def _apply_landsat_scale_factors(image):
    """Internal function: Applies scaling factors to Landsat C2 L2 SR bands."""
    scale_factor = 0.0000275
    offset = -0.2
    # Select optical bands (starting with SR_B)
    optical_bands = image.select('SR_B.*')
    scaled_bands = optical_bands.multiply(scale_factor).add(offset)
    # Add scaled bands back, overwriting originals, keep other bands
    return image.addBands(scaled_bands, None, True).copyProperties(image, ['system:time_start'])

def calculate_landsat_indices(year, roi, month=None):
    """
    Calculates Landsat NDVI and NIRv for a given year/month and ROI (using median composite).

    Args:
        year (int): The year.
        roi (ee.Geometry): The region of interest.
        month (int, optional): The month (1-12). If None, processes growing season (May-Sep).

    Returns:
        dict: Contains 'NDVI', 'NIRv' ee.Image objects, or None if processing fails.
    """
    if roi is None:
        return None

    try:
        # Determine date range and time suffix for naming
        if month:
            start_date = f"{year}-{month:02d}-01"
            if month == 12:
                end_date = f"{year}-12-31"
            else:
                # Calculate the last day of the month
                end_date = (datetime.date(year, month + 1, 1) - datetime.timedelta(days=1)).strftime('%Y-%m-%d')
            time_suffix = f"{year}_{month:02d}"
            period_desc = f"{year}-{month:02d}"
        else:
            start_date = f"{year}-05-01"  # Growing season start
            end_date = f"{year}-09-30"    # Growing season end
            time_suffix = f"{year}"
            period_desc = f"{year} (Growing Season)"

        # Select appropriate Landsat collection(s)
        collection_l8 = 'LANDSAT/LC08/C02/T1_L2'
        collection_l9 = 'LANDSAT/LC09/C02/T1_L2'
        if year < 1999: # L5
            collections_to_process = ['LANDSAT/LT05/C02/T1_L2']
            cloud_prop = 'CLOUD_COVER'
        elif year < 2013: # L7
            collections_to_process = ['LANDSAT/LE07/C02/T1_L2']
            cloud_prop = 'CLOUD_COVER'
        elif year < 2022: # L8 only
            collections_to_process = [collection_l8]
            cloud_prop = 'CLOUD_COVER_LAND'
        else: # L8 + L9
            collections_to_process = [collection_l8, collection_l9]
            cloud_prop = 'CLOUD_COVER_LAND'

        # --- Data Acquisition and Preprocessing ---
        initial_thr, relaxed_thr, min_imgs = 20, 50, 1 # Min images needed for median

        landsat_merged = ee.ImageCollection([])
        # Initial filtering
        for coll_id in collections_to_process:
            try:
                landsat_merged = landsat_merged.merge(
                    ee.ImageCollection(coll_id)
                    .filterBounds(roi)
                    .filterDate(start_date, end_date)
                    .filter(ee.Filter.lt(cloud_prop, initial_thr))
                )
            except:
                continue

        # Relax cloud cover filter if too few images
        initial_size = landsat_merged.size().getInfo()
        if initial_size < min_imgs:
            landsat_merged = ee.ImageCollection([]) # Reset collection
            for coll_id in collections_to_process:
                try:
                    landsat_merged = landsat_merged.merge(
                        ee.ImageCollection(coll_id)
                        .filterBounds(roi)
                        .filterDate(start_date, end_date)
                        .filter(ee.Filter.lt(cloud_prop, relaxed_thr)) # Use relaxed filter
                    )
                except:
                    continue

        if year < 2013: # L5, L7
            l57_bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']
            std_bands = ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']
            
            landsat_processed = landsat_merged.map(_apply_landsat_scale_factors) \
                                             .map(_mask_clouds_shadows_landsat) \
                                             .select(l57_bands, std_bands)
        else: # L8, L9
            l89_bands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']
            std_bands = ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']
            
            landsat_processed = landsat_merged.map(_apply_landsat_scale_factors) \
                                             .map(_mask_clouds_shadows_landsat) \
                                             .select(l89_bands, std_bands)

        # Check if any images remain after processing
        if landsat_processed.size().getInfo() == 0:
            return None

        # --- Create Median Composite ---
        median_composite = landsat_processed.median().clip(roi)

        # --- Calculate Indices ---
        # Calculate NDVI: (NIR - Red) / (NIR + Red)
        # Use try-except to handle potential errors if bands are missing in the composite
        try:
            ndvi = median_composite.normalizedDifference(['NIR', 'Red']).rename(f'NDVI_{time_suffix}')

            # Calculate NIRv: NDVI * NIR
            nirv = ndvi.multiply(median_composite.select('NIR')).rename(f'NIRv_{time_suffix}')

            return {'NDVI': ndvi, 'NIRv': nirv}

        except:
            return None

    except:
        return None

# --- Sentinel-2 Data Processing ---

def _mask_clouds_s2_qa60(image):
    """Internal function: Masks clouds and cirrus in Sentinel-2 using QA60."""
    qa60 = image.select('QA60')
    # Bits 10: Cloud, Bits 11: Cirrus
    cloud_mask = qa60.bitwiseAnd(1 << 10).eq(0) # Not Cloud
    cirrus_mask = qa60.bitwiseAnd(1 << 11).eq(0) # Not Cirrus
    final_mask = cloud_mask.And(cirrus_mask)
    # Apply mask, select optical bands, scale, copy time
    optical_bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'] # Blue, Green, Red, NIR, SWIR1, SWIR2
    return image.select(optical_bands) \
                .updateMask(final_mask) \
                .divide(10000.0) \
                .copyProperties(image, ['system:time_start'])

def calculate_sentinel_indices(year, roi, month=None):
    """
    Calculates Sentinel-2 NDVI and NIRv for a given year/month and ROI (using median composite).

    Args:
        year (int): The year (>= 2017 recommended for S2_SR_HARMONIZED).
        roi (ee.Geometry): The region of interest.
        month (int, optional): The month (1-12). If None, processes growing season (May-Sep).

    Returns:
        dict: Contains 'NDVI', 'NIRv' ee.Image objects, or None if processing fails.
    """
    if roi is None or year < 2017:
        return None

    try:
        # Determine date range and time suffix
        if month:
            start_date = f"{year}-{month:02d}-01"
            if month == 12:
                end_date = f"{year}-12-31"
            else:
                end_date = (datetime.date(year, month + 1, 1) - datetime.timedelta(days=1)).strftime('%Y-%m-%d')
            time_suffix = f"{year}_{month:02d}"
        else:
            start_date = f"{year}-05-01"  # Growing season start
            end_date = f"{year}-09-30"    # Growing season end
            time_suffix = f"{year}"

        # --- Data Acquisition and Preprocessing ---
        collection_id = 'COPERNICUS/S2_SR_HARMONIZED'
        # Select optical bands + QA60 for masking
        bands_select_s2 = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'QA60']

        initial_thr, relaxed_thr = 20, 50

        # Initial filtering
        try:
            sentinel_collection = ee.ImageCollection(collection_id) \
                .filterBounds(roi) \
                .filterDate(start_date, end_date) \
                .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', initial_thr)) \
                .select(bands_select_s2)
        except:
            return None

        # Relax cloud cover filter if too few images
        if sentinel_collection.size().getInfo() < 1:
            try:
                sentinel_collection = ee.ImageCollection(collection_id) \
                    .filterBounds(roi) \
                    .filterDate(start_date, end_date) \
                    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', relaxed_thr)) \
                    .select(bands_select_s2)
            except:
                return None

        # Apply masking (using QA60 function)
        sentinel_masked = sentinel_collection.map(_mask_clouds_s2_qa60)

        # Check if any images remain after processing
        if sentinel_masked.size().getInfo() == 0:
            return None

        # --- Create Median Composite ---
        median_composite = sentinel_masked.median().clip(roi)

        # --- Calculate Indices ---
        # Calculate NDVI: (NIR - Red) / (NIR + Red) -> (B8 - B4) / (B8 + B4)
        try:
            ndvi = median_composite.normalizedDifference(['B8', 'B4']).rename(f'NDVI_{time_suffix}')

            # Calculate NIRv: NDVI * NIR -> NDVI * B8
            nirv = ndvi.multiply(median_composite.select('B8')).rename(f'NIRv_{time_suffix}')

            return {'NDVI': ndvi, 'NIRv': nirv}

        except:
            return None

    except:
        return None


# --- Data Download Main Function ---

def download_satellite_data(roi, start_year=1993, end_year=2024, parent_folder='GEEpreprocessing', monthly=False, months=None):
    """
    Downloads yearly (growing season) or monthly satellite data (NDVI/NIRv).
    Data is exported at 30m resolution in EPSG:28992 CRS.

    Args:
        roi (ee.Geometry): Region of interest.
        start_year (int): Start year.
        end_year (int): End year.
        parent_folder (str): Base folder for results.
        monthly (bool): If True, download monthly data.
        months (list[int], optional): List of months (1-12) for monthly download. Defaults to all 12.
    """
    if roi is None:
        raise ValueError("Error: ROI must be provided.")

    # Determine months and output folder structure
    if monthly:
        if months is None:
            months_to_process = list(range(1, 13))
        elif isinstance(months, list) and all(isinstance(m, int) and 1 <= m <= 12 for m in months):
            months_to_process = months
        else:
            raise TypeError("`months` must be a list of integers (1-12) or None.")
        period_type = "Monthly"
        output_folder = os.path.join(parent_folder, 'monthly')
    else:
        period_type = "Annual (Growing Season)"
        output_folder = os.path.join(parent_folder, 'annual_growing_season')
        months_to_process = [None] # Use None to signify annual processing

    print(f"--- Starting {period_type} Download ({start_year}-{end_year}) ---")
    print(f"--- Output Directory: {output_folder} ---")

    os.makedirs(output_folder, exist_ok=True)

    # Process each year and month (if applicable)
    for year in range(start_year, end_year + 1):
        for month in months_to_process:
            period_desc = f"{year}_{month:02d}" if month else f"{year}"
            print(f"\n--- Processing {period_desc} ---")

            try:
                indices = None
                # Try Sentinel-2 first for recent years
                if year >= 2017:
                    indices = calculate_sentinel_indices(year, roi, month=month)

                # Fallback to Landsat if needed
                if indices is None and year >= 1984:
                    indices = calculate_landsat_indices(year, roi, month=month)

                # Export results if available
                if indices and (indices.get('NDVI') is not None or indices.get('NIRv') is not None):
                    _export_data(indices, year, output_folder, roi, month=month)
                else:
                    print(f"⏩ Skipping export for {period_desc}: No valid data.")
            except:
                print(f"⏩ Skipping export for {period_desc}: Processing error.")

    print(f"\n✅ {period_type} Download Finished ({start_year}-{end_year})")


# --- Data Export Helper Function ---

def _export_data(indices, year, output_folder, roi, month=None):
    """Internal helper function: Exports NDVI and NIRv images."""
    if roi is None or not indices:
        return

    resolution = 30
    crs = 'EPSG:28992' # RD New coordinate system for the Netherlands
    time_suffix = f"{year}_{month:02d}" if month else f"{year}"

    # Export NDVI
    ndvi_image = indices.get('NDVI')
    if isinstance(ndvi_image, ee.Image):
        filename_ndvi = os.path.join(output_folder, f"NDVI_{time_suffix}.tif")
        try:
            print(f" Preparing to export NDVI: {os.path.basename(filename_ndvi)}")
            geemap.download_ee_image(
                ndvi_image,
                filename=filename_ndvi,
                region=roi,
                scale=resolution,
                crs=crs
            )
            print(f"✅ NDVI Exported: {os.path.basename(filename_ndvi)}")
        except Exception as e:
            print(f"❌ Error exporting NDVI: {e}")

    # Export NIRv
    nirv_image = indices.get('NIRv')
    if isinstance(nirv_image, ee.Image):
        filename_nirv = os.path.join(output_folder, f"NIRv_{time_suffix}.tif")
        try:
            print(f" Preparing to export NIRv: {os.path.basename(filename_nirv)}")
            geemap.download_ee_image(
                nirv_image,
                filename=filename_nirv,
                region=roi,
                scale=resolution,
                crs=crs
            )
            print(f"✅ NIRv Exported: {os.path.basename(filename_nirv)}")
        except Exception as e:
            print(f"❌ Error exporting NIRv: {e}")




