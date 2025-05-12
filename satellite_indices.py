import ee
import geemap
import os
import datetime
import numpy as np

# 使用本地shapefile定义研究区域
def get_roi_from_shapefile(shapefile_path):
    """从shapefile获取研究区域的几何形状"""
    roi_asset = geemap.shp_to_ee(shapefile_path)
    if isinstance(roi_asset, ee.FeatureCollection):
        roi_geometry = roi_asset.geometry()
    else:
        roi_geometry = roi_asset
    return roi_geometry

# Landsat数据处理函数
def calculate_landsat_indices(year, month=None):
    """
    计算Landsat影像的NDVI和EVI指数
    
    参数:
        year: 年份
        month: 月份（如果提供，则处理特定月份的数据；否则处理整个生长季）
        
    返回:
        包含 'Original', 'NDVI', 'EVI' 键的字典，
        如果无可用影像则值为 None。
    """
    
    if month:
        # 处理特定月份的数据
        start_date = f"{year}-{month:02d}-01"
        # 计算当月最后一天
        if month == 12:
            next_year = year + 1
            next_month = 1
        else:
            next_year = year
            next_month = month + 1
        
        end_date = f"{next_year}-{next_month:02d}-01"
        # 将结束日期向前调整一天
        end_date = (datetime.datetime.strptime(end_date, "%Y-%m-%d") - datetime.timedelta(days=1)).strftime("%Y-%m-%d")
        time_suffix = f"{year}_{month:02d}"  # 用于波段命名
    else:
        # 处理整个生长季数据
        start_date = f"{year}-05-01"  # 生长季开始
        end_date = f"{year}-09-30"    # 生长季结束
        time_suffix = f"{year}"       # 用于波段命名
    
    # 选择适当的Landsat数据集
    if year < 1999:
        # Landsat 5
        collection = 'LANDSAT/LT05/C02/T1_L2'
        nir_band = 'SR_B4'
        red_band = 'SR_B3'
        blue_band = 'SR_B1'
        swir1_band = 'SR_B5'
        swir2_band = 'SR_B7'
        green_band = 'SR_B2'
        scale_factor = 0.0000275
        offset = -0.2
    elif year < 2013:
        # Landsat 7
        collection = 'LANDSAT/LE07/C02/T1_L2'
        nir_band = 'SR_B4'
        red_band = 'SR_B3'
        blue_band = 'SR_B1'
        swir1_band = 'SR_B5'
        swir2_band = 'SR_B7'
        green_band = 'SR_B2'
        scale_factor = 0.0000275
        offset = -0.2
    else:
        # Landsat 8/9
        # Note: Ensure roi is defined globally or passed as an argument
        if roi is None:
             raise ValueError("ROI must be set before calling calculate_landsat_indices")
             
        if year >= 2022:
             collection_l9 = 'LANDSAT/LC09/C02/T1_L2'
             collection_l8 = 'LANDSAT/LC08/C02/T1_L2'
             nir_band = 'SR_B5'
             red_band = 'SR_B4'
             blue_band = 'SR_B2'
             swir1_band = 'SR_B6'
             swir2_band = 'SR_B7'
             green_band = 'SR_B3'
             scale_factor = 0.0000275
             offset = -0.2
        else: # 2013 to 2021
             collection = 'LANDSAT/LC08/C02/T1_L2'
             nir_band = 'SR_B5'
             red_band = 'SR_B4'
             blue_band = 'SR_B2'
             swir1_band = 'SR_B6'
             swir2_band = 'SR_B7'
             green_band = 'SR_B3'
             scale_factor = 0.0000275
             offset = -0.2
             
    # 定义云掩膜函数
    def mask_clouds_landsat(image):
        # Bits 3 (Cloud Shadow), 5 (Cloud), 6 (Cirrus) should be cleared (set to 0).
        qa = image.select('QA_PIXEL')
        cloud_shadow_mask = (1 << 3)
        cloud_mask = (1 << 5)
        # Landsat 8/9 also have a cirrus mask (Bit 6)
        cirrus_mask = (1 << 6) 
        
        mask = qa.bitwiseAnd(cloud_shadow_mask).eq(0) \
                 .And(qa.bitwiseAnd(cloud_mask).eq(0)) \
                 .And(qa.bitwiseAnd(cirrus_mask).eq(0)) # Also check cirrus for L8/9
        
        # Also check QA_RADSAT for saturated pixels (Bands 1-7 for L8/9)
        radsat = image.select('QA_RADSAT')
        saturation_mask = radsat.bitwiseAnd(int('0000111', 2)).eq(0) # Check bits 0, 1, 2 (Bands 1, 2, 3)
                                                                    # Adjust if other bands need checking

        return image.updateMask(mask).updateMask(saturation_mask)

    # 应用比例因子和偏移量转换为反射率
    def apply_scale_factors(image):
        bands_to_scale = [red_band, nir_band, blue_band, green_band, swir1_band, swir2_band]
        optical_bands = image.select(bands_to_scale) \
            .multiply(scale_factor).add(offset)
        thermal_bands = image.select('ST_B.*').multiply(0.00341802).add(149.0) # Scale thermal bands if needed
        # Add scaled bands back, overwriting original optical bands but keeping thermal
        return image.addBands(optical_bands, overwrite=True).addBands(thermal_bands, overwrite=True)


    # 获取Landsat数据
    if year >= 2022:
        landsat_l9 = ee.ImageCollection(collection_l9) \
            .filterBounds(roi) \
            .filterDate(start_date, end_date) \
            .filter(ee.Filter.lt('CLOUD_COVER_LAND', 70)) # Use CLOUD_COVER_LAND for L8/9
            
        landsat_l8 = ee.ImageCollection(collection_l8) \
            .filterBounds(roi) \
            .filterDate(start_date, end_date) \
            .filter(ee.Filter.lt('CLOUD_COVER_LAND', 70))
            
        # Merge L8 and L9
        landsat_merged = ee.ImageCollection(landsat_l8.merge(landsat_l9))
        
        # Apply masking and scaling
        landsat = landsat_merged.map(mask_clouds_landsat).map(apply_scale_factors)
            
    else: # Before 2022 (L5, L7, or L8 only)
        landsat_col = ee.ImageCollection(collection) \
            .filterBounds(roi) \
            .filterDate(start_date, end_date) 
            
        # Adjust cloud filter property name
        cloud_property = 'CLOUD_COVER_LAND' if year >= 2013 else 'CLOUD_COVER'
        landsat_col = landsat_col.filter(ee.Filter.lt(cloud_property, 70)) # Higher initial threshold

        # Apply masking and scaling
        landsat = landsat_col.map(mask_clouds_landsat).map(apply_scale_factors)


    # --- 检查集合大小 ---
    collection_size = landsat.size().getInfo()
    if collection_size == 0:
        print(f"⚠️ {time_suffix}: 没有可用的 Landsat 影像，跳过此时间段。")
        return {'Original': None, 'NDVI': None, 'EVI': None}
    # --- 检查结束 ---
    
    # 选择近红外波段质量最好的影像
    best_image = landsat.qualityMosaic(nir_band).clip(roi)
    
    # 计算NDVI: (NIR - RED) / (NIR + RED)
    ndvi = best_image.normalizedDifference([nir_band, red_band]).rename(f'NDVI_{time_suffix}')
    
    # 计算EVI: 2.5 * ((NIR - RED) / (NIR + 6*RED - 7.5*BLUE + 1))
    # Check if necessary bands exist before calculating EVI
    required_bands = [nir_band, red_band, blue_band]
    has_bands = all(band in best_image.bandNames().getInfo() for band in required_bands)

    if has_bands:
        evi = best_image.expression(
            '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
            {
                'NIR': best_image.select(nir_band),
                'RED': best_image.select(red_band),
                'BLUE': best_image.select(blue_band)
            }
        ).rename(f'EVI_{time_suffix}')
    else:
         print(f"⚠️ {time_suffix}: 缺少计算 EVI 所需波段 ({required_bands}) in best_image. Setting EVI to None.")
         evi = None # Or handle as appropriate, maybe create an image with a specific nodata value


    
    # 重命名波段以保持一致性
    original_bands = [blue_band, green_band, red_band, nir_band, swir1_band, swir2_band]
    new_band_names = [f'Blue_{time_suffix}', f'Green_{time_suffix}', f'Red_{time_suffix}', 
                      f'NIR_{time_suffix}', f'SWIR1_{time_suffix}', f'SWIR2_{time_suffix}']
    
    # Select only bands that actually exist in best_image
    available_original_bands = [b for b in original_bands if b in best_image.bandNames().getInfo()]
    available_new_names = [new_band_names[original_bands.index(b)] for b in available_original_bands]

    if available_original_bands:
        renamed_bands = best_image.select(available_original_bands, available_new_names)
    else:
        renamed_bands = None # No bands to rename

    
    return {
        'Original': renamed_bands,
        'NDVI': ndvi,
        'EVI': evi
    }

# Sentinel-2数据处理函数
def calculate_sentinel_indices(year, month=None):
    """
    计算Sentinel-2影像的NDVI和EVI指数
    
    参数:
        year: 年份
        month: 月份（如果提供，则处理特定月份的数据；否则处理整个生长季）
        
    返回:
        包含 'Original', 'NDVI', 'EVI' 键的字典，
        如果无可用影像则值为 None。
    """
    if roi is None:
        raise ValueError("ROI must be set before calling calculate_sentinel_indices")

    if month:
        # 处理特定月份的数据
        start_date = f"{year}-{month:02d}-01"
        # 计算当月最后一天
        if month == 12:
            next_year = year + 1
            next_month = 1
        else:
            next_year = year
            next_month = month + 1
        
        end_date = f"{next_year}-{next_month:02d}-01"
        # 将结束日期向前调整一天
        end_date = (datetime.datetime.strptime(end_date, "%Y-%m-%d") - datetime.timedelta(days=1)).strftime("%Y-%m-%d")
        time_suffix = f"{year}_{month:02d}"  # 用于波段命名
    else:
        # 处理整个生长季数据
        start_date = f"{year}-05-01"  # 标准窗口
        end_date = f"{year}-09-30"
        time_suffix = f"{year}"       # 用于波段命名
    
    # 定义云掩膜函数 - Sentinel-2 (Using SCL band)
    def mask_clouds_s2(image):
        scl = image.select('SCL')
        # Values to mask: 1 (saturated), 3 (cloud shadow), 8 (cloud medium prob), 
        # 9 (cloud high prob), 10 (cirrus), 11 (snow/ice)
        # Keep: 2 (dark features), 4 (veg), 5 (bare soil), 6 (water), 7 (unclassified)
        undesired_pixels = scl.eq(1).Or(scl.eq(3)).Or(scl.eq(8)).Or(scl.eq(9)).Or(scl.eq(10)).Or(scl.eq(11))
        # Invert the mask: mask pixels that are undesired
        mask = undesired_pixels.Not() 
        
        # Also consider masking based on QA60 for extra safety if needed
        # qa60 = image.select('QA60')
        # cloud_mask_qa60 = qa60.bitwiseAnd(1 << 10).Or(qa60.bitwiseAnd(1 << 11)) # Cloud or Cirrus
        # combined_mask = mask.And(cloud_mask_qa60.Not())
        # return image.updateMask(combined_mask).divide(10000).copyProperties(image, ['system:time_start'])

        # Apply scale factor (divide by 10000 for SR) and copy time property
        return image.updateMask(mask).divide(10000).copyProperties(image, ['system:time_start'])

    
    # 获取Sentinel-2数据 (Level-2A preferred)
    try:
         # Try L2A first
         sentinel_col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
             .filterBounds(roi) \
             .filterDate(start_date, end_date) \
             .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 70)) # Higher initial threshold
    except ee.EEException as e:
         # Fallback to L1C if L2A is not available for the period (e.g., early dates)
         # Note: L1C requires atmospheric correction or careful interpretation
         print(f"⚠️ {time_suffix}: Sentinel-2 L2A SR collection might not be available. Error: {e}. Check date range or GEE catalog.")
         # Handle this case - maybe raise error, or try L1C (requires different processing)
         # For now, we'll return None as SR data isn't directly available
         return {'Original': None, 'NDVI': None, 'EVI': None}

    # Apply cloud masking
    sentinel = sentinel_col.map(mask_clouds_s2)

    # --- 检查集合大小 ---
    collection_size = sentinel.size().getInfo()
    if collection_size == 0:
        print(f"⚠️ {time_suffix}: 没有可用的 Sentinel-2 影像，跳过此时间段。")
        return {'Original': None, 'NDVI': None, 'EVI': None}
    # --- 检查结束 ---
    
    # 使用近红外波段(B8)质量最好的像素
    best_image = sentinel.qualityMosaic('B8').clip(roi) # B8 is NIR for Sentinel-2
    
    # 计算NDVI: (NIR - RED) / (NIR + RED) => (B8 - B4) / (B8 + B4)
    ndvi = best_image.normalizedDifference(['B8', 'B4']).rename(f'NDVI_{time_suffix}')
    
    # 计算EVI: 2.5 * ((NIR - RED) / (NIR + 6*RED - 7.5*BLUE + 1)) => 2.5 * ((B8 - B4) / (B8 + 6*B4 - 7.5*B2 + 1))
    # Check if necessary bands exist before calculating EVI
    required_bands = ['B8', 'B4', 'B2']
    has_bands = all(band in best_image.bandNames().getInfo() for band in required_bands)

    if has_bands:
        evi = best_image.expression(
            '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
            {
                'NIR': best_image.select('B8'),
                'RED': best_image.select('B4'),
                'BLUE': best_image.select('B2')
            }
        ).rename(f'EVI_{time_suffix}')
    else:
        print(f"⚠️ {time_suffix}: 缺少计算 EVI 所需波段 ({required_bands}) in best_image. Setting EVI to None.")
        evi = None

    
    # 重命名波段以保持一致性
    original_bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'] # Blue, Green, Red, NIR, SWIR1, SWIR2
    new_band_names = [f'Blue_{time_suffix}', f'Green_{time_suffix}', f'Red_{time_suffix}', 
                      f'NIR_{time_suffix}', f'SWIR1_{time_suffix}', f'SWIR2_{time_suffix}']

    # Select only bands that actually exist in best_image
    available_original_bands = [b for b in original_bands if b in best_image.bandNames().getInfo()]
    available_new_names = [new_band_names[original_bands.index(b)] for b in available_original_bands]

    if available_original_bands:
        renamed_bands = best_image.select(available_original_bands, available_new_names)
    else:
        renamed_bands = None # No bands to rename

    
    return {
        'Original': renamed_bands,
        'NDVI': ndvi,
        'EVI': evi
    }

# 初始化全局变量roi，在调用函数前需要设置
roi = None

def set_roi(shapefile_path):
    """设置全局研究区域变量"""
    global roi
    roi = get_roi_from_shapefile(shapefile_path)
    print("ROI set successfully.") # Add confirmation
    return roi

def download_satellite_data(start_year=1993, end_year=2024, parent_folder='GEEpreprocessing', monthly=False, months=None):
    """
    下载特定年份或月份的卫星数据，计算NDVI和EVI指数
    所有数据统一采用30m分辨率和RD New坐标系
    
    参数:
        start_year: 起始年份
        end_year: 结束年份
        parent_folder: 保存结果的文件夹
        monthly: 是否按月下载数据
        months: 需要下载的月份列表 (int list)，如果为None则下载所有12个月
    """
    if roi is None:
        raise ValueError("ROI must be set using set_roi() before downloading data.")

    if monthly:
        if months is None:
            months = list(range(1, 13))  # 默认下载所有12个月
        elif not isinstance(months, list) or not all(isinstance(m, int) for m in months):
             raise TypeError("`months` parameter must be a list of integers or None.")
        period_desc = f"按月({','.join([str(m) for m in months])})"
    else:
        period_desc = "按年度"
    
    print(f"开始{period_desc}处理{start_year}-{end_year}年的卫星数据...")
    
    # 确保输出文件夹存在
    if not os.path.exists(parent_folder):
        os.makedirs(parent_folder)
        print(f"创建目录: {parent_folder}")
    
    # 处理每一年的数据
    for year in range(start_year, end_year + 1):
        if monthly:
            # 按月处理数据
            for month in months:
                print(f"\n--- 处理 {year}年 {month}月 数据 ---")
                try:
                    # 根据年份选择不同的卫星数据
                    if year < 2013: # Landsat 5/7 era
                        print(f"使用 Landsat 5/7 数据...")
                        indices = calculate_landsat_indices(year, month)
                    elif year < 2017: # Landsat 8 era starts, Sentinel-2 may not have L2A
                        print(f"使用 Landsat 8 数据...")
                        indices = calculate_landsat_indices(year, month)
                        # Optional: try Sentinel-2 L1C if needed and implement its processing
                    else: # Sentinel-2 L2A likely available + Landsat 8/9
                        print(f"尝试 Sentinel-2 (优先) 或 Landsat 8/9 数据...")
                        indices = calculate_sentinel_indices(year, month)
                        # Fallback to Landsat if Sentinel fails or returns None
                        if indices is None or indices.get('NDVI') is None:
                             print(f"⚠️ Sentinel-2 在 {year}-{month:02d} 无有效数据或计算失败，尝试 Landsat...")
                             indices = calculate_landsat_indices(year, month)

                    # 检查计算结果是否有效
                    if indices and indices.get('NDVI') is not None:
                         _export_data(indices, year, parent_folder, month=month)
                    else:
                         print(f"⏩ 跳过导出 {year}年{month}月 的数据，因为没有有效影像或计算失败。")

                except ee.EEException as e:
                    print(f"❌ 处理 {year}年{month}月 数据时 GEE 出错: {e}")
                    print(f"   跳过此月份...")
                except Exception as e: # 捕获其他可能的非 GEE 错误
                    print(f"❌ 处理 {year}年{month}月 数据时发生意外错误: {e}")
                    print(f"   跳过此月份...")
        else:
            # 按年处理数据 (生长季)
            print(f"\n--- 处理 {year}年 生长季数据 ---")
            try:
                 # 根据年份选择不同的卫星数据
                if year < 2013: # Landsat 5/7 era
                    print(f"使用 Landsat 5/7 数据...")
                    indices = calculate_landsat_indices(year)
                elif year < 2017: # Landsat 8 era starts, Sentinel-2 may not have L2A
                    print(f"使用 Landsat 8 数据...")
                    indices = calculate_landsat_indices(year)
                else: # Sentinel-2 L2A likely available + Landsat 8/9
                    print(f"尝试 Sentinel-2 (优先) 或 Landsat 8/9 数据...")
                    indices = calculate_sentinel_indices(year)
                     # Fallback to Landsat if Sentinel fails or returns None
                    if indices is None or indices.get('NDVI') is None:
                        print(f"⚠️ Sentinel-2 在 {year}年生长季 无有效数据或计算失败，尝试 Landsat...")
                        indices = calculate_landsat_indices(year)

                # 检查计算结果是否有效
                if indices and indices.get('NDVI') is not None:
                     _export_data(indices, year, parent_folder)
                else:
                     print(f"⏩ 跳过导出 {year}年 的数据，因为没有有效影像或计算失败。")

            except ee.EEException as e:
                print(f"❌ 处理 {year}年 数据时 GEE 出错: {e}")
                print(f"   跳过此年份...")
            except Exception as e: # 捕获其他可能的非 GEE 错误
                print(f"❌ 处理 {year}年 数据时发生意外错误: {e}")
                print(f"   跳过此年份...")

    print(f"\n✅✅✅ 所有数据处理完成，文件已保存至 {parent_folder} 目录 ✅✅✅")


def _export_data(indices, year, parent_folder, month=None):
    """
    导出处理后的数据到本地 (内部辅助函数)
    
    参数:
        indices: 计算后的指数字典 (必须包含有效的 'NDVI' 和 'EVI')
        year: 年份
        parent_folder: 保存结果的文件夹
        month: 月份(int, 如果有)
    """
    if roi is None:
        print("❌ 导出错误: ROI 未设置。")
        return
        
    # 所有数据统一使用30m分辨率 (Landsat native resolution)
    # Sentinel-2 will be resampled during export
    resolution = 30
    
    # 构建时间后缀
    if month:
        time_suffix = f"{year}_{month:02d}"
    else:
        time_suffix = f"{year}"
    
    # 导出NDVI
    ndvi_image = indices.get('NDVI')
    if ndvi_image is not None:
        filename_ndvi = os.path.join(parent_folder, f"NDVI_{time_suffix}.tif")
        print(f"   准备导出 NDVI: {filename_ndvi}")
        try:
            geemap.download_ee_image(
                ndvi_image, 
                filename=filename_ndvi,
                region=roi,
                scale=resolution,
                crs='EPSG:28992'
            )
            print(f"✅ {time_suffix} NDVI 导出完成")
        except ee.EEException as e:
             print(f"❌ 导出 {time_suffix} NDVI 时 GEE 出错: {e}")
        except Exception as e:
             print(f"❌ 导出 {time_suffix} NDVI 时发生意外错误: {e}")
    else:
        print(f"ℹ️ {time_suffix}: 无有效 NDVI 数据可导出。")

    
    # 导出EVI
    evi_image = indices.get('EVI')
    if evi_image is not None:
        filename_evi = os.path.join(parent_folder, f"EVI_{time_suffix}.tif")
        print(f"   准备导出 EVI: {filename_evi}")
        try:
            geemap.download_ee_image(
                evi_image, 
                filename=filename_evi,
                region=roi,
                scale=resolution,
                crs='EPSG:28992'
            )
            print(f"✅ {time_suffix} EVI 导出完成")
        except ee.EEException as e:
             print(f"❌ 导出 {time_suffix} EVI 时 GEE 出错: {e}")
        except Exception as e:
             print(f"❌ 导出 {time_suffix} EVI 时发生意外错误: {e}")
    else:
        print(f"ℹ️ {time_suffix}: 无有效 EVI 数据可导出。") 