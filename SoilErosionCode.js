var table = ee.FeatureCollection("users/joeyqmf83/China");

var year = 2020   // define year

// ********** R Factor **********//

var year_num = 5   // five-year average
var R_sum = ee.Image(0).clip(table)
function cal_R(img) {
  var month = img.select("precipitation").clip(table)
  var r = month.multiply(730)
  var r_factor = img.expression(
    "-1.15527+1.792*r",
    {
      "r": r
    }
  ).rename('R');
  return r_factor;
}
for(var i=0;i<year_num;i++){
  var dataset = ee.ImageCollection('NASA/GPM_L3/IMERG_MONTHLY_V06').filterDate(year-i+'-01-01', year+'-12-31').map(cal_R);
  var Single_R_factor = dataset.select('R').sum().clip(table);
  R_sum = R_sum.add(Single_R_factor)
}
var R_factor = R_sum.divide(year_num)
Map.addLayer(R_factor, {min: 300, max: 20000, palette: ['a52508','ff3818','fbff18','25cdff','2f35ff','0b2dab']}, 'R Factor Map', 0);




// ********** K Factor **********//

var SIL = ee.Image("users/joeyqmf83/fensha_500m"),
    CLA = ee.Image("users/joeyqmf83/nianli_500m"),
    SAN = ee.Image("users/joeyqmf83/xisha_500m"),
    C = ee.Image("users/joeyqmf83/youjizhi_500m");
    
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
var K_factor = cal_K(SAN, SIL, CLA, C, SN1).multiply(0.1317).clip(table)
Map.addLayer(K_factor, {min: 0, max: 0.05, palette: ['a52508','ff3818','fbff18','25cdff','2f35ff','0b2dab']}, 'K Factor Map', 0);


// ********** LS Factor **********//

var dataset = ee.Image('USGS/SRTMGL1_003').clip(table);
var terrain = ee.Algorithms.Terrain(dataset);
var slope = terrain.select('slope');

var threshold1 = 5
var threshold2 = 10
var Good_mask = slope.lt(threshold1)
var Mid_mask = slope.gte(threshold1).and(slope.lt(threshold2))
var Bad_mask = slope.gte(threshold2)

function Level1(img) {
  var slope = img.select("slope");
  var s_factor = img.expression(
    "10.8*SLOPE+0.03",
    {
      "SLOPE": (slope.multiply(Math.PI).divide(180)).sin(),
    }
  ).rename('Level1');
  return s_factor;
}
function Level2(img) {
  var slope = img.select("slope");
  var s_factor = img.expression(
    "16.8*SLOPE-0.5",
    {
      "SLOPE": (slope.multiply(Math.PI).divide(180)).sin(),
    }
  ).rename('Level2');
  return s_factor;
}
function Level3(img) {
  var slope = img.select("slope");
  var s_factor = img.expression(
    "21.9*SLOPE-0.96",
    {
      "SLOPE": (slope.multiply(Math.PI).divide(180)).sin(),
    }
  ).rename('Level3');
  return s_factor;
}

var S = slope.addBands(Level1(slope))
var S = S.addBands(Level2(slope))
var S = S.addBands(Level3(slope))


var S_factor = S.expression(
      "(b('slope') < 5)? b('Level1')" +
      ": (b('slope') >= 5) and (b('slope')<10)? b('Level2')" +
      ": b('Level3') "
).rename('S_factor').clip(table);





function cal_slope_length(img) {
  var slope = img.select("slope");
  var slope_length = img.expression(
    "DEM/SLOPE",
    {
      "DEM": dataset.select('elevation'),
      "SLOPE": (slope.multiply(Math.PI).divide(180)).sin()
    }
  ).rename('slope_length');
  return slope_length;
  }

var L_ = cal_slope_length(slope)


function cal_Beta(img){
  var slope = img.select('slope')
var B = img.expression(
  "(SLOPE/0.3)/(3*(SLOPE**0.8)+0.56)",
  {
    "SLOPE": ((slope.multiply(Math.PI).divide(180)).sin()),
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

function cal_LS(img, m_img){
  var SLope_length = img
  var m = m_img.select('m')
  var ls = img.expression(
    "(lamda/22.13)**m",
    {
      "m": m,
      "lamda": SLope_length.abs(),
    }
  ).rename('LS');
  return ls;
}
var Beta = cal_Beta(terrain)  
var m = cal_m(Beta)
var LS_factor = cal_LS(L_, m).multiply(S_factor)

LS_factor = LS_factor.where(LS_factor.gt(45), 45)

Map.addLayer(LS_factor, {min: 0, max: 50, palette: ['a52508','ff3818','fbff18','25cdff','2f35ff','0b2dab']}, 'LS Factor Map', 0);




// ********** C Factor **********//

function FVC(img) {
 var ndvi = img.select("NDVI");
 var fvc = img.expression(
   "(NDVI/10000 - 0.2)/(0.65- 0.2)",   // It is relatively time-consuming to calculate the vegetation coverage of large areas by GEE, and this value is the pre-loaded multi-year average
   {
     "NDVI": ndvi,
   }
 ).rename('FVC');
 return fvc;
}

var dataset = ee.ImageCollection('MODIS/061/MOD13Q1').filter(ee.Filter.date(year+'-01-01', year+'-12-31'));
var ndvi = dataset.select('NDVI').median().clip(table);
var fvc = FVC(ndvi)  // Get FVC

function Cal(img) {
  var fvc = img.select("FVC");
  var c_factor = img.expression(
    "0.6508-(0.3436*a)",
    {
      "a": fvc.log(),
    }
  ).rename('C_Factor');
  return c_factor;
}

fvc = fvc.addBands(Cal(fvc))

var C_factor = fvc.expression(
      "(b('FVC') <= 0) ? 1" +
      ": (b('FVC') >= 0.783) ? 0" +
      ": (b('C_Factor')>=0) ? b('C_Factor')"+
      ": 0 "
).rename('C').clip(table);

Map.addLayer (C_factor, {min: 0, max: 1, palette:  ['000000','006600','009900', '33CC00', '996600', 'CC9900','CC9966','FFFFFF']}, 'C Factor Map',0);


// ********** P Factor **********//

var modis = ee.ImageCollection('MODIS/006/MCD12Q1')
var lulc = modis.filterDate(year+'-01-01', year+'-12-31').select('LC_Type1')
        .first().clip(table).rename('lulc');

var lulc_slope = lulc.addBands(slope)   // Combined LULC & slope in single image
// Create P Facor map using an expression
var P_factor = lulc_slope.expression(
      "(b('lulc') < 11) ? 0.6" +
      ": (b('lulc') == 11) ? 0" +
      ": (b('lulc') == 13) ? 0.2" +
      ": (b('lulc') == 15) ? 0.1" +
      ": (b('lulc') == 16) ? 1" +
      ": (b('lulc') == 17) ? 0" +
      ": (b('slope') < 2) and((b('lulc')==12) or (b('lulc')==14)) ? 0.4" +
      ": (b('slope') < 5) and((b('lulc')==12) or (b('lulc')==14)) ? 0.5" +
      ": (b('slope') < 8) and((b('lulc')==12) or (b('lulc')==14)) ? 0.5" +
      ": (b('slope') < 12) and((b('lulc')==12) or (b('lulc')==14)) ? 0.6" +
      ": (b('slope') < 16) and((b('lulc')==12) or (b('lulc')==14)) ? 0.7" +
      ": (b('slope') < 20) and((b('lulc')==12) or (b('lulc')==14)) ? 0.8" +
      ": (b('slope') > 20) and((b('lulc')==12) or (b('lulc')==14)) ? 0.9" +
      ": 1"
).rename('P').clip(table);
Map.addLayer (P_factor, {min: 0, max: 1, palette: ['FFFFFF','CC9966','CC9900', '996600', '33CC00', '009900','006600','000000']}, 'P Factor Map', 0)

// **********  Estimating Soil Erosion **********//
var SoilErosion = ee.Image(R_factor.multiply(K_factor).multiply(C_factor).multiply(P_factor)).multiply(LS_factor).rename("Soil Loss")

Map.addLayer (SoilErosion, {min: 0, max:25, palette: ['490eff','12f4ff','12ff50','e5ff12','ff4812']}, 'SoilErosion',1)


var BNU_2015 = ee.Image("projects/ee-joeyqmf83/assets/Soil_Erosion_Other_Result/BNU_2015")
var GloSEM_2012 = ee.Image("projects/ee-joeyqmf83/assets/Soil_Erosion_Other_Result/GloSEM_2012")
var Jialei_Li_2015 = ee.Image("projects/ee-joeyqmf83/assets/Soil_Erosion_Other_Result/Jialei_Li_2015")


// ********** Other Dataset **********//
var style = ['490eff','12f4ff','12ff50','e5ff12','ff4812']
Map.addLayer (BNU_2015.multiply(0.01), {min: 0, max:100, palette: style}, 'BNU_2015',0)
Map.addLayer (Jialei_Li_2015, {min: 0, max:100, palette: style}, 'Jialei_Li_2015',0)
Map.addLayer (GloSEM_2012, {min: 0, max:100, palette: style}, 'GloSEM_2012',0)


