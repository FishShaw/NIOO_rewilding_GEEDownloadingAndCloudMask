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

def _mask_clouds_shadows_landsat(image, relaxed_mode=False):
    """Internal function: Masks clouds, shadows, and snow in Landsat Collection 2 Level 2."""
    qa = image.select('QA_PIXEL')
    
    if relaxed_mode:
        # Relaxed mode: only mask explicit clouds, preserve more pixels
        cloud_mask = qa.bitwiseAnd(1 << 3).eq(0)      # Not Cloud
        # Don't mask shadows and snow to preserve more data
        final_mask = cloud_mask
        print("      🔧 Using relaxed masking mode (clouds only)")
    else:
        # Standard mode
        cloud_mask = qa.bitwiseAnd(1 << 3).eq(0)      # Not Cloud
        shadow_mask = qa.bitwiseAnd(1 << 4).eq(0)     # Not Cloud Shadow
        snow_mask = qa.bitwiseAnd(1 << 5).eq(0)       # Not Snow
        final_mask = cloud_mask.And(shadow_mask).And(snow_mask)
    
    # Apply mask and copy time properties
    return image.updateMask(final_mask).copyProperties(image, ['system:time_start'])

def _apply_landsat_scale_factors(image):
    """Internal function: Applies scaling factors to Landsat Collection 2 Level 2 Surface Reflectance bands."""
    scale_factor = 0.0000275
    offset = -0.2
    # Select optical bands (starting with SR_B)
    optical_bands = image.select('SR_B.*')
    scaled_bands = optical_bands.multiply(scale_factor).add(offset)
    # Add scaled bands back, overwriting originals, keep other bands
    return image.addBands(scaled_bands, None, True).copyProperties(image, ['system:time_start'])

def _get_seasonal_params(month, max_coverage_mode=False):
    """Adjust quality control parameters based on season and mode"""
    if max_coverage_mode:
        # Maximum coverage mode: extremely relaxed parameters
        if month in [12, 1, 2]:  # Winter
            return {
                'initial_thr': 40, 'relaxed_thr': 95, 'final_thr': 99,
                'min_imgs': 1, 'time_buffers': [7, 15, 30],
                'use_relaxed_masking': True
            }
        elif month in [6, 7, 8]:  # Summer
            return {
                'initial_thr': 30, 'relaxed_thr': 80, 'final_thr': 95,
                'min_imgs': 1, 'time_buffers': [5, 10, 21],
                'use_relaxed_masking': False
            }
        else:  # Spring/Autumn
            return {
                'initial_thr': 35, 'relaxed_thr': 85, 'final_thr': 97,
                'min_imgs': 1, 'time_buffers': [7, 14, 25],
                'use_relaxed_masking': True
            }
    else:
        # Standard mode: current parameters
        return {
            'initial_thr': 20, 'relaxed_thr': 90, 'final_thr': 90,
            'min_imgs': 1, 'time_buffers': [0],
            'use_relaxed_masking': False
        }

def _process_landsat_collection(collections_to_process, roi, start_date, end_date, cloud_prop, 
                               initial_thr, relaxed_thr, final_thr=None, use_relaxed_masking=False):
    """Common function for processing Landsat image collections"""
    
    landsat_merged = ee.ImageCollection([])
    
    # First round: initial threshold
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
    
    # Second round: relaxed threshold
    if landsat_merged.size().getInfo() < 1:
        landsat_merged = ee.ImageCollection([])
        for coll_id in collections_to_process:
            try:
                landsat_merged = landsat_merged.merge(
                    ee.ImageCollection(coll_id)
                    .filterBounds(roi)
                    .filterDate(start_date, end_date)
                    .filter(ee.Filter.lt(cloud_prop, relaxed_thr))
                )
            except:
                continue
    
    # Third round: final threshold (only in max coverage mode)
    if final_thr and final_thr > relaxed_thr and landsat_merged.size().getInfo() < 1:
        landsat_merged = ee.ImageCollection([])
        for coll_id in collections_to_process:
            try:
                landsat_merged = landsat_merged.merge(
                    ee.ImageCollection(coll_id)
                    .filterBounds(roi)
                    .filterDate(start_date, end_date)
                    .filter(ee.Filter.lt(cloud_prop, final_thr))
                )
            except:
                continue
        print(f"      🚨 Using final threshold {final_thr}%")
    
    return landsat_merged

def calculate_landsat_indices(year, roi, month=None, max_coverage_mode=False):
    """
    Calculates Landsat NDVI and NIRv for a given year/month and ROI using median composite.
    
    Features:
    - Adaptive time window expansion for data-sparse periods
    - Maximum coverage mode for continuity over quality
    - Multi-sensor support (Landsat 5/7/8/9)
    - Intelligent cloud masking strategies
    """
    if roi is None:
        return None

    try:
        # Get seasonal parameters
        if month:
            params = _get_seasonal_params(month, max_coverage_mode)
        else:
            params = _get_seasonal_params(6, max_coverage_mode)
        
        # Determine date range and time suffix for naming
        if month:
            base_start_date = f"{year}-{month:02d}-01"
            if month == 12:
                base_end_date = f"{year}-12-31"
            else:
                base_end_date = (datetime.date(year, month + 1, 1) - datetime.timedelta(days=1)).strftime('%Y-%m-%d')
            time_suffix = f"{year}_{month:02d}"
            period_desc = f"{year}-{month:02d}"
        else:
            base_start_date = f"{year}-05-01"  # Growing season start
            base_end_date = f"{year}-09-30"    # Growing season end
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

        # --- Adaptive time window processing ---
        landsat_processed = None
        
        for attempt, time_buffer in enumerate(params['time_buffers']):
            # Calculate expanded time range
            if time_buffer > 0:
                start_date_obj = datetime.datetime.strptime(base_start_date, '%Y-%m-%d').date()
                end_date_obj = datetime.datetime.strptime(base_end_date, '%Y-%m-%d').date()
                
                expanded_start_date = (start_date_obj - datetime.timedelta(days=time_buffer)).strftime('%Y-%m-%d')
                expanded_end_date = (end_date_obj + datetime.timedelta(days=time_buffer)).strftime('%Y-%m-%d')
                
                print(f"      📅 Expanding time window ±{time_buffer} days: {expanded_start_date} to {expanded_end_date}")
                start_date, end_date = expanded_start_date, expanded_end_date
            else:
                start_date, end_date = base_start_date, base_end_date
            
            # Process image collection
            if attempt == 0:
                # First attempt: standard parameters
                landsat_merged = _process_landsat_collection(
                    collections_to_process, roi, start_date, end_date, cloud_prop,
                    params['initial_thr'], params['relaxed_thr']
                )
            else:
                # Subsequent attempts: use final threshold
                landsat_merged = _process_landsat_collection(
                    collections_to_process, roi, start_date, end_date, cloud_prop,
                    params['initial_thr'], params['relaxed_thr'], params['final_thr']
                )
            
            if landsat_merged.size().getInfo() >= params['min_imgs']:
                print(f"      ✅ Found {landsat_merged.size().getInfo()} images")
                break
            else:
                print(f"      ⚠️  Attempt {attempt + 1}: Only found {landsat_merged.size().getInfo()} images")
        
        # If still not enough images, return None
        if landsat_merged.size().getInfo() < params['min_imgs']:
            print(f"      ❌ All attempts failed")
            return None

        # Process images
        if year < 2013: # L5, L7
            l57_bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']
            std_bands = ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']
            
            landsat_processed = landsat_merged.map(_apply_landsat_scale_factors) \
                                             .map(lambda img: _mask_clouds_shadows_landsat(img, params['use_relaxed_masking'])) \
                                             .select(l57_bands, std_bands)
        else: # L8, L9
            l89_bands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']
            std_bands = ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']
            
            landsat_processed = landsat_merged.map(_apply_landsat_scale_factors) \
                                             .map(lambda img: _mask_clouds_shadows_landsat(img, params['use_relaxed_masking'])) \
                                             .select(l89_bands, std_bands)

        # Check if any images remain after processing
        if landsat_processed.size().getInfo() == 0:
            return None

        # --- Create Median Composite ---
        median_composite = landsat_processed.median().clip(roi)

        # --- Calculate Indices ---
        try:
            ndvi = median_composite.normalizedDifference(['NIR', 'Red']).rename(f'NDVI_{time_suffix}')
            nirv = ndvi.multiply(median_composite.select('NIR')).rename(f'NIRv_{time_suffix}')
            return {'NDVI': ndvi, 'NIRv': nirv}
        except:
            return None

    except:
        return None

# --- Sentinel-2 Data Processing ---

def _mask_clouds_s2_improved(image, relaxed_mode=False):
    """Improved Sentinel-2 cloud masking with relaxed mode option."""
    
    if relaxed_mode:
        # Relaxed mode: only mask most explicit clouds
        if image.bandNames().contains('SCL').getInfo():
            scl = image.select('SCL')
            # Only mask high and medium clouds, preserve more pixels
            clear_mask = scl.neq(8).And(scl.neq(9))  # Don't mask shadows and cirrus
            print("      🔧 Using relaxed masking mode (high/medium clouds only)")
        else:
            qa60 = image.select('QA60')
            clear_mask = qa60.bitwiseAnd(1 << 10).eq(0)  # Only mask clouds, not cirrus
        
        optical_bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']
        return image.select(optical_bands) \
                    .updateMask(clear_mask) \
                    .divide(10000.0) \
                    .copyProperties(image, ['system:time_start'])
    else:
        # Standard mode: use original strict masking
        if image.bandNames().contains('SCL').getInfo():
            scl = image.select('SCL')
            clear_mask = scl.neq(3).And(scl.neq(8)).And(scl.neq(9)).And(scl.neq(10))
            clear_mask = clear_mask.And(scl.neq(1)).And(scl.neq(2))
        else:
            qa60 = image.select('QA60')
            cloud_mask = qa60.bitwiseAnd(1 << 10).eq(0)
            cirrus_mask = qa60.bitwiseAnd(1 << 11).eq(0)
            clear_mask = cloud_mask.And(cirrus_mask)
        
        optical_bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']
        return image.select(optical_bands) \
                    .updateMask(clear_mask) \
                    .divide(10000.0) \
                    .copyProperties(image, ['system:time_start'])

def _process_sentinel_collection(roi, start_date, end_date, initial_thr, relaxed_thr, final_thr=None):
    """Common function for processing Sentinel-2 image collections"""
    
    collection_id = 'COPERNICUS/S2_SR_HARMONIZED'
    bands_select_s2 = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'SCL', 'QA60']
    
    # First round: initial threshold
    try:
        sentinel_collection = ee.ImageCollection(collection_id) \
            .filterBounds(roi) \
            .filterDate(start_date, end_date) \
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', initial_thr)) \
            .select(bands_select_s2)
    except:
        sentinel_collection = ee.ImageCollection([])
    
    # Second round: relaxed threshold
    if sentinel_collection.size().getInfo() < 1:
        try:
            sentinel_collection = ee.ImageCollection(collection_id) \
                .filterBounds(roi) \
                .filterDate(start_date, end_date) \
                .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', relaxed_thr)) \
                .select(bands_select_s2)
        except:
            sentinel_collection = ee.ImageCollection([])
    
    # Third round: final threshold
    if final_thr and final_thr > relaxed_thr and sentinel_collection.size().getInfo() < 1:
        try:
            sentinel_collection = ee.ImageCollection(collection_id) \
                .filterBounds(roi) \
                .filterDate(start_date, end_date) \
                .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', final_thr)) \
                .select(bands_select_s2)
            print(f"      🚨 Using final threshold {final_thr}%")
        except:
            sentinel_collection = ee.ImageCollection([])
    
    return sentinel_collection

def calculate_sentinel_indices(year, roi, month=None, max_coverage_mode=False):
    """
    Calculates Sentinel-2 NDVI and NIRv for a given year/month and ROI (using median composite).
    Enhanced version: supports adaptive time window and maximum coverage mode
    """
    if roi is None or year < 2017:
        return None

    try:
        # Get seasonal parameters
        if month:
            params = _get_seasonal_params(month, max_coverage_mode)
        else:
            params = _get_seasonal_params(6, max_coverage_mode)
            
        # Determine date range and time suffix
        if month:
            base_start_date = f"{year}-{month:02d}-01"
            if month == 12:
                base_end_date = f"{year}-12-31"
            else:
                base_end_date = (datetime.date(year, month + 1, 1) - datetime.timedelta(days=1)).strftime('%Y-%m-%d')
            time_suffix = f"{year}_{month:02d}"
        else:
            base_start_date = f"{year}-05-01"  # Growing season start
            base_end_date = f"{year}-09-30"    # Growing season end
            time_suffix = f"{year}"

        # --- Adaptive time window processing ---
        sentinel_masked = None
        
        for attempt, time_buffer in enumerate(params['time_buffers']):
            # Calculate expanded time range
            if time_buffer > 0:
                start_date_obj = datetime.datetime.strptime(base_start_date, '%Y-%m-%d').date()
                end_date_obj = datetime.datetime.strptime(base_end_date, '%Y-%m-%d').date()
                
                expanded_start_date = (start_date_obj - datetime.timedelta(days=time_buffer)).strftime('%Y-%m-%d')
                expanded_end_date = (end_date_obj + datetime.timedelta(days=time_buffer)).strftime('%Y-%m-%d')
                
                print(f"      📅 Expanding time window ±{time_buffer} days: {expanded_start_date} to {expanded_end_date}")
                start_date, end_date = expanded_start_date, expanded_end_date
            else:
                start_date, end_date = base_start_date, base_end_date
            
            # Process image collection
            if attempt == 0:
                sentinel_collection = _process_sentinel_collection(
                    roi, start_date, end_date, params['initial_thr'], params['relaxed_thr']
                )
            else:
                sentinel_collection = _process_sentinel_collection(
                    roi, start_date, end_date, params['initial_thr'], params['relaxed_thr'], params['final_thr']
                )
            
            if sentinel_collection.size().getInfo() >= params['min_imgs']:
                print(f"      ✅ Found {sentinel_collection.size().getInfo()} images")
                # Apply improved masking
                sentinel_masked = sentinel_collection.map(lambda img: _mask_clouds_s2_improved(img, params['use_relaxed_masking']))
                break
            else:
                print(f"      ⚠️  Attempt {attempt + 1}: Only found {sentinel_collection.size().getInfo()} images")
        
        # If still not enough images, return None
        if sentinel_masked is None or sentinel_masked.size().getInfo() < params['min_imgs']:
            print(f"      ❌ All attempts failed")
            return None

        # --- Create Median Composite ---
        median_composite = sentinel_masked.median().clip(roi)

        # --- Calculate Indices ---
        try:
            ndvi = median_composite.normalizedDifference(['B8', 'B4']).rename(f'NDVI_{time_suffix}')
            nirv = ndvi.multiply(median_composite.select('B8')).rename(f'NIRv_{time_suffix}')
            return {'NDVI': ndvi, 'NIRv': nirv}
        except:
            return None

    except:
        return None

# --- Data Download Main Function ---

def download_satellite_data(roi, start_year=1993, end_year=2024, parent_folder='GEEpreprocessing', 
                           monthly=False, months=None, max_coverage_mode=False):
    """
    Downloads yearly (growing season) or monthly satellite data (NDVI/NIRv) for the study area.
    Data is exported at 30m resolution in RD New (EPSG:28992) coordinate system.
    
    Key Features:
    - Multi-sensor integration: Landsat 5/7/8/9 + Sentinel-2
    - Adaptive quality control: Standard vs maximum coverage modes
    - Intelligent cloud masking with seasonal parameter adjustment
    - Comprehensive time series coverage (1993-2024)

    Args:
        roi (ee.Geometry): Region of interest geometry.
        start_year (int): Start year for data acquisition.
        end_year (int): End year for data acquisition.
        parent_folder (str): Base directory for exported GeoTIFF files.
        monthly (bool): If True, download monthly composites; if False, annual growing season.
        months (list[int], optional): Specific months (1-12) for monthly download. Defaults to all 12.
        max_coverage_mode (bool): If True, prioritize data continuity over strict quality thresholds.
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
        
        if max_coverage_mode:
            period_type = "Monthly (Max Coverage Mode)"
            output_folder = os.path.join(parent_folder, 'monthly_maxcov')
        else:
            period_type = "Monthly"
            output_folder = os.path.join(parent_folder, 'monthly')
    else:
        period_type = "Annual (Growing Season)"
        output_folder = os.path.join(parent_folder, 'annual_growing_season')
        months_to_process = [None] # Use None to signify annual processing

    print(f"--- Starting {period_type} Download ({start_year}-{end_year}) ---")
    if max_coverage_mode and monthly:
        print(f"🔧 Maximum coverage mode enabled: Prioritizing data continuity")
    print(f"--- Output Directory: {output_folder} ---")

    os.makedirs(output_folder, exist_ok=True)

    # Statistics tracking
    total_periods = (end_year - start_year + 1) * len(months_to_process)
    successful_exports = 0
    failed_periods = []

    # Process each year and month (if applicable)
    for year in range(start_year, end_year + 1):
        for month in months_to_process:
            period_desc = f"{year}_{month:02d}" if month else f"{year}"
            print(f"\n--- Processing {period_desc} ---")

            try:
                indices = None
                # Try Sentinel-2 first for recent years
                if year >= 2017:
                    print(f"    🛰️  Trying Sentinel-2...")
                    indices = calculate_sentinel_indices(year, roi, month=month, max_coverage_mode=max_coverage_mode)

                # Fallback to Landsat if needed
                if indices is None and year >= 1984:
                    print(f"    🛰️  Falling back to Landsat...")
                    indices = calculate_landsat_indices(year, roi, month=month, max_coverage_mode=max_coverage_mode)

                # Export results if available
                if indices and (indices.get('NDVI') is not None or indices.get('NIRv') is not None):
                    _export_data(indices, year, output_folder, roi, month=month)
                    successful_exports += 1
                    print(f"✅ {period_desc} successfully exported")
                else:
                    print(f"⏩ Skipping export for {period_desc}: No valid data.")
                    failed_periods.append(period_desc)
            except Exception as e:
                print(f"⏩ Skipping export for {period_desc}: Processing error - {str(e)}")
                failed_periods.append(period_desc)

    # Print final statistics
    print(f"\n{'='*50}")
    print(f"✅ {period_type} Download Finished ({start_year}-{end_year})")
    print(f"📊 Success rate: {successful_exports}/{total_periods} ({successful_exports/total_periods*100:.1f}%)")
    
    if failed_periods:
        print(f"❌ Failed periods: {len(failed_periods)} items")
        if len(failed_periods) <= 20:  # Only show first 20
            print(f"   Details: {', '.join(failed_periods)}")
        else:
            print(f"   Details: {', '.join(failed_periods[:20])}... (and {len(failed_periods)-20} more)")
    
    if max_coverage_mode and monthly:
        print(f"💡 Tip: If gaps still exist, consider:")
        print(f"   1. Using subsequent time series interpolation methods")
        print(f"   2. Checking specific failed months for image availability")
    print(f"{'='*50}")

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
            print(f"   📤 Preparing to export NDVI: {os.path.basename(filename_ndvi)}")
            geemap.download_ee_image(
                ndvi_image,
                filename=filename_ndvi,
                region=roi,
                scale=resolution,
                crs=crs
            )
            print(f"   ✅ NDVI export completed: {os.path.basename(filename_ndvi)}")
        except Exception as e:
            print(f"   ❌ NDVI export error: {e}")

    # Export NIRv
    nirv_image = indices.get('NIRv')
    if isinstance(nirv_image, ee.Image):
        filename_nirv = os.path.join(output_folder, f"NIRv_{time_suffix}.tif")
        try:
            print(f"   📤 Preparing to export NIRv: {os.path.basename(filename_nirv)}")
            geemap.download_ee_image(
                nirv_image,
                filename=filename_nirv,
                region=roi,
                scale=resolution,
                crs=crs
            )
            print(f"   ✅ NIRv export completed: {os.path.basename(filename_nirv)}")
        except Exception as e:
            print(f"   ❌ NIRv export error: {e}")




