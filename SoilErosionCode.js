var year=2010;
var ProvinceList = ['China_shp']

var table = ee.FeatureCollection("users/joeyqmf83/China_shp");
//Get NDVI
var landsat_ndvi_fvc = require('users/joeyqmf83/hhx:STLS/GetAnnualFVC');
var fvc =landsat_ndvi_fvc.FVC
//Get Factor R
var GPM_R = require('users/joeyqmf83/hhx:STLS/GetAnnualR');
var r =GPM_R.R


var threshold = 0.783
var Good_mask = fvc.gte(threshold)
var Bad_mask = fvc.lt(threshold).and(fvc.gt(0))
var worst_mask = fvc.lte(0)

// Cover management factor
function C_Bad(img) {
var fvc = img.select("FVC");
var c_factor = img.expression(
  "0.6508-(0.3436*FVC)",
  {
    "FVC": fvc.log(),
  }
).rename('C_Factor');
return c_factor;
}
var Bad = C_Bad(fvc)
var C_Factor = Bad.updateMask(Bad_mask.or(worst_mask)).unmask(0)
var C_Factor = C_Factor.updateMask(Good_mask.or(Bad_mask)).unmask(1)


//  *****Rainfall erosivity factor
var R_factor = r.select('R').sum().clip(table);


//  *****Soil erodibility factor
var fensha = ee.Image("users/joeyqmf83/fensha_500m"),
    nianli = ee.Image("users/joeyqmf83/nianli_500m"),
    xisha = ee.Image("users/joeyqmf83/xisha_500m"),
    youjizhi = ee.Image("users/joeyqmf83/youjizhi_500m");

var palettes = require('users/gena/packages:palettes');

var SAN = xisha
var SIL = fensha
var CLA = nianli
var C = youjizhi
var SN1 = (SAN.multiply(-1).add(1)).divide(100)


function cal_K(SAN, SIL, CLA, C, SN1) {
  var A = ((((SAN.multiply((SIL.divide(-100)).add(1))).multiply(-0.0256)).exp()).multiply(0.3)).add(0.2)
  var B = SIL.divide(CLA.add(SIL))
  var C_ = (C.multiply(0.25)).divide(C.add(((C.multiply(-2.95)).add(3.72)).exp()))
  var D = (SN1.multiply(0.7)).divide(SN1.add(((SN1.multiply(22.9)).add(-5.51)).exp()))
  
  var k_factor = SAN.expression(
    "A*(B**0.3)*(1-C)*(1-D)",
    {
      "A": A,
      "B": B,
      "C": C_,
      "D": D,
    }
  ).rename('K_Factor');
  return k_factor;
  }

var k_factor = cal_K(SAN, SIL, CLA, C, SN1).multiply(0.1317)


//  *****Slope length and steepness factor
var roi = table;
var dataset_dem = ee.Image('USGS/SRTMGL1_003').clip(roi);
var terrain = ee.Algorithms.Terrain(dataset_dem);
var slope = terrain.select('slope');

var Flow_accumulation = ee.Image("MERIT/Hydro/v1_0_1").select('upa').clip(table);

var threshold1 = 5
var threshold2 = 10
var Good_mask = slope.lt(threshold1)
var Mid_mask = slope.gte(threshold1).and(slope.lt(threshold2))
var Bad_mask = slope.gte(threshold2)

function S_Good(img) {
var slope = img.select("slope");
var s_factor = img.expression(
  "10.8*SLOPE+0.036",
  {
    "SLOPE": (slope.multiply(Math.PI).divide(180)).sin(),
  }
).rename('S_Factor');
return s_factor;
}
function S_Mid(img) {
var slope = img.select("slope");
var s_factor = img.expression(
  "16.8*SLOPE-0.5",
  {
    "SLOPE": (slope.multiply(Math.PI).divide(180)).sin(),
  }
).rename('S_Factor');
return s_factor;
}
function S_Bad(img) {
var slope = img.select("slope");
var s_factor = img.expression(
  "21.9*SLOPE-0.96",
  {
    "SLOPE": (slope.multiply(Math.PI).divide(180)).sin(),
  }
).rename('S_Factor');
return s_factor;
}
var Good = slope.updateMask(Mid_mask.or(Bad_mask)).unmask(S_Good(slope.updateMask(Good_mask)))
var Mid = Good.updateMask(Good_mask.or(Bad_mask)).unmask(S_Mid(slope.updateMask(Mid_mask)))
var S = Mid.updateMask(Good_mask.or(Mid_mask)).unmask(S_Bad(slope.updateMask(Bad_mask)))
// var Good = S_Good(slope)
// var Mid = Good.updateMask(Good_mask).unmask(S_Mid(slope.updateMask(Mid_mask)))
// var S = Good.updateMask(Good_mask.or(Mid_mask)).unmask(S_Bad(slope.updateMask(Bad_mask)))


function cal_slope_length(img) {
  var slope = img.select("slope");
  var slope_length = img.expression(
    "DEM/SLOPE",
    {
      "DEM": dataset_dem.select('elevation'),
      "SLOPE": (slope.multiply(Math.PI).divide(180)).sin()
    }
  ).rename('slope_length');
  return slope_length;
  }

var L_ = Flow_accumulation


function cal_B(img){
  var slope = img.select('slope')
var B = img.expression(
  "SLOPE/(3*(SLOPE**0.8)+0.56)",
  {
    "SLOPE": ((slope.multiply(Math.PI).divide(180)).sin()).divide(0.0896),
  }
).rename('B');
return B;
}
function cal_m(img){
  var b = img.select('B')
var m = img.expression(
  "B/(B+1)",
  {
    "B": b,
  }
).rename('m');
return m;
}

function cal_LS(img, m_img, slope_img){
  var SLope_length = img
  var m = m_img.select('m')
  var slope = slope_img.select('slope')
var ls = img.expression(
  "((lamda/22.13)**m)*slope",
  {
    "m": m,
    "lamda": SLope_length.abs(),
    "slope": slope,
    
  }
).rename('LS');
return ls;
}

var B = cal_B(terrain)  
var m = cal_m(B)
var LS = cal_LS(L_, m, S)


//  *****Support practice factor

var dataset = ee.ImageCollection('MODIS/006/MCD12Q1').filter(ee.Filter.date(year+'-01-01', year+'-12-31'));
var igbpLandCover = dataset.select('LC_Type1').first().clip(table);


var eq_one_mask = igbpLandCover.lte(9).or(igbpLandCover.eq(14))
var eq_zero_mask = igbpLandCover.eq(11).or(igbpLandCover.eq(15)).or(igbpLandCover.eq(17)).or(igbpLandCover.eq(13))
var eq_crop_mask = igbpLandCover.eq(12)
var _eq_1 = igbpLandCover.updateMask(eq_zero_mask.or(eq_crop_mask)).unmask(1).clip(table)
var _eq_0 = _eq_1.updateMask(eq_one_mask.or(eq_crop_mask)).unmask(0).clip(table)

var _eq_crop_1 = _eq_0.updateMask(eq_one_mask.or(eq_zero_mask).or(slope.lte(25))).unmask(0.8).clip(table)
var _eq_crop_2 = _eq_crop_1.updateMask(eq_one_mask.or(eq_zero_mask).or(slope.lte(20)).or(slope.gt(25))).unmask(0.705).clip(table)
var _eq_crop_3 = _eq_crop_2.updateMask(eq_one_mask.or(eq_zero_mask).or(slope.lte(15)).or(slope.gt(20))).unmask(0.575).clip(table)
var _eq_crop_4 = _eq_crop_3.updateMask(eq_one_mask.or(eq_zero_mask).or(slope.lte(10)).or(slope.gt(15))).unmask(0.305).clip(table)
var _eq_crop_5 = _eq_crop_4.updateMask(eq_one_mask.or(eq_zero_mask).or(slope.lte(5)).or(slope.gt(10))).unmask(0.221).clip(table)
var _eq_crop_6 = _eq_crop_5.updateMask(eq_one_mask.or(eq_zero_mask).or(slope.lte(0)).or(slope.gt(5))).unmask(0.1).clip(table)


var Soil = C_Factor.multiply(R_factor).multiply(k_factor).multiply(LS).multiply(_eq_crop_6).multiply(100)
var palettes = require('users/gena/packages:palettes');
var palette = palettes.crameri.lajolla[10];
Map.addLayer(Soil, palette, year)

