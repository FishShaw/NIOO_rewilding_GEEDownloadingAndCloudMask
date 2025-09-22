# 🌿 Quantifying Rewilding Impact: 32-Year Vegetation Productivity and Resilience Analysis

**Master's Research Internship | Netherlands Institute of Ecology (NIOO-KNAW) | Geo-Information Science**

Transforming decades of satellite data into actionable insights about ecosystem restoration effectiveness in the Gelderse Poort floodplain, Netherlands.

![Study Area Overview](Docs/Images/study_area.png)
![Land Cover Classification](Docs/Images/layouot_classification.png)
![Detailed Classification](Docs/Images/layout_classification_zoomin.png)

## 🎯 Research Innovation

Despite the implementation of extensive rewilding projects in the Gelderse Poort over the past three decades, **quantitative assessment of their impact on floodplain vegetation remains a knowledge gap**. This research addresses this critical need by analyzing vegetation productivity dynamics across a 32-year period (1993-2024) using advanced geospatial data processing and time series analysis.

### 🔬 Scientific Contribution
This interdisciplinary project bridges **Restoration Ecology** and **Geo-Information Science**, developing novel approaches to assess ecological resilience using remote sensing data. The framework provides critical insights into how rewilding modifies both baseline ecosystem functioning and functional resilience to climatic stressors.

![Research Workflow](Docs/Images/workflow.png)

## 🚀 Technical Innovation Highlights

### 🛰️ **Multi-Platform Satellite Data Integration**
- **32-year time series** (1993-2024) spanning multiple sensor generations
- **Multi-sensor fusion**: Landsat 5/7/8/9 + Sentinel-2 harmonization
- **Intelligent data continuity**: Adaptive quality control with seasonal parameter adjustment
- **Cloud-native processing**: Google Earth Engine for petabyte-scale analysis

### 🌦️ **Advanced Weather Data Processing**
- **Multi-station meteorological integration**: KNMI weather station networks
- **Extreme event identification**: Automated drought/flood detection algorithms
- **Climatological indices**: SPI/SPEI calculation for resilience assessment
- **Quality assurance**: Comprehensive data validation and gap-filling strategies

### 📊 **Innovative Resilience Metrics**
- **Three-dimensional resilience assessment**: Impact magnitude, relative recovery, absolute recovery
- **Seven independent extreme events**: 3 droughts, 4 floods for robust statistical analysis
- **Chronosequence approach**: Four rewilding sites with different establishment years
- **Multi-scale comparison**: Landscape-wide vs site-specific analysis


## 🔧 Technical Stack & Workflow

### **Core Technologies**
```
🐍 Python Ecosystem
├── Google Earth Engine API - Cloud-native satellite data processing
├── geemap - Interactive geospatial analysis
├── xarray/GDAL - Multi-dimensional data handling
├── pandas/numpy - Statistical computing
└── geopandas - Spatial data operations

📊 Statistical Analysis (Linked Repository)
├── R Statistical Computing
├── Time Series Analysis
├── GLMMs (Generalized Linear Mixed Models)
└── Spatial Autocorrelation Analysis
```

### **Data Processing Pipeline**

```mermaid
flowchart TD
    A[Data Acquisition] --> B[Multi-Sensor Integration]
    
    B --> C1[Historical Data\n1993-2015]
    B --> C2[Recent Data\n2015-2024]
    
    C1 --> D1[Landsat 5 TM\n1993-2011]
    C1 --> D2[Landsat 7 ETM+\n1999-2013]
    C1 --> D3[Landsat 8 OLI\n2013-2015]
    C2 --> D4[Sentinel-2 MSI\n2015-2024]
    
    D1 & D2 & D3 & D4 --> E[Spatial Filtering\nGelderse Poort Boundary]
    
    E --> F[Quality Control]
    
    F --> G1[Temporal Filtering\nSeasonal/Annual]
    F --> G2[Cloud Coverage\n<20% threshold]
    F --> G3[Adaptive Masking\nSeasonal Parameters]
    
    G1 & G2 & G3 --> H[Preprocessing]
    
    H --> I1[Radiometric Correction]
    H --> I2[Atmospheric Correction\nTOA to Surface Reflectance]
    H --> I3[Cloud/Shadow Detection]
    
    I1 & I2 & I3 --> J[Vegetation Index Calculation\nNDVI & NIRv]
    
    J --> K[Temporal Compositing]
    
    K --> L1[Long-term Products]
    K --> L2[Disturbance Analysis]
    
    L1 --> M1[Monthly Composites]
    L1 --> M2[Seasonal Aggregation]
    L1 --> M3[Annual Growing Season]
    
    L2 --> N1[Drought Events]
    L2 --> N2[Flood Events]
    
    M1 & M2 & M3 & N1 & N2 --> O[Data Export]
    
    O --> P1[GeoTIFF Format\nRD New Projection]
    O --> P2[Metadata Documentation\nJSON/CSV]
    
    P1 & P2 --> Q[Quality Assessment]
    
    Q --> R[Analysis Ready Data]
```

## 📈 Key Technical Achievements

### **1. Intelligent Data Continuity System**
```python
# Adaptive time window expansion for data-sparse periods
def _get_seasonal_params(month, max_coverage_mode=False):
    """Adjust quality control parameters based on season and coverage mode"""
    if max_coverage_mode:
        # Maximum coverage mode: prioritize data continuity
        return {
            'initial_thr': 30, 'relaxed_thr': 80, 'final_thr': 95,
            'min_imgs': 1, 'time_buffers': [5, 10, 21],
            'use_relaxed_masking': True
        }
```

### **2. Multi-Sensor Harmonization**
- **Seamless sensor transitions**: Automated detection and processing of optimal satellite data
- **Standardized band mapping**: Consistent spectral response across Landsat generations
- **Quality-first approach**: Hierarchical sensor selection (Sentinel-2 → Landsat 8/9 → Landsat 7/5)

## 🛠️ System Requirements & Setup

### **Core Dependencies**
```bash
pip install earthengine-api geemap numpy pandas jupyter
```

### **Getting Started**
```python
# Initialize Google Earth Engine
import ee
ee.Initialize(project='your-gee-project')

# Load study area and begin processing
from satellite_indices_new import get_roi_from_shapefile, download_satellite_data

roi = get_roi_from_shapefile('shapefile/GeldersePoort_cliped.shp')
download_satellite_data(roi=roi, start_year=1993, end_year=2024)
```

### **Data Structure**
```
📁 Project Structure
├── 🛰️ GEEpreprocessing/
│   ├── yearly/annual_growing_season/ (32 years × 2 indices)
│   └── monthly/monthly_maxcov/ (384 months × 2 indices)
├── 🌦️ weatherdata_postprocessing/
│   ├── processed_weather_data.csv
│   ├── drought_events.csv
│   └── flood_events.csv
├── 🗺️ shapefile/ (Study area boundaries)
└── 📊 Docs/Images/ (Workflow diagrams)
```

## 🔬 Technical Specifications

- **Spatial Resolution**: 30m × 30m pixels
- **Temporal Coverage**: 1993-2024 (32 years)
- **Coordinate System**: RD New (EPSG:28992)
- **Data Volume**: >1000 satellite images processed
- **Processing Environment**: Google Earth Engine cloud platform
- **Export Format**: GeoTIFF with comprehensive metadata

## 🤝 Collaboration & Integration
Subsequent analysis see:  https://github.com/FishShaw/NIOO_rewilding_TimeSeriesModellingAndStats

This preprocessing pipeline seamlessly integrates with advanced statistical analysis workflows:
- **R-based time series analysis** for trend detection
- **GLMM modeling** for multi-factor resilience assessment
- **Spatial autocorrelation analysis** for landscape-scale patterns
- **Visualization frameworks** for stakeholder communication

---
