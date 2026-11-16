// Poluição do ar na Baixada Santista por satélite (Sentinel-5P, MODIS MAIAC e ERA5-Land).
// Cole no Code Editor do Google Earth Engine e clique em Run. Na aba Tasks, rode as 9 exportações
// (as tabelas diárias podem levar de 10 a 40 minutos). Os arquivos vão para a pasta GEE_poluicao do Drive.

var area = ee.Geometry.Rectangle([-47.1, -24.45, -45.9, -23.4]);

// ---------- zonas de análise (círculos de 4 km)
function zona(nome, lon, lat) { return ee.Feature(ee.Geometry.Point([lon, lat]).buffer(4000), {zona: nome}); }
var zonas = ee.FeatureCollection([
  zona('Santos_Porto', -46.315, -23.950),
  zona('Cubatao_polo', -46.395, -23.865),
  zona('Sao_Vicente', -46.390, -23.960),
  zona('Guaruja', -46.250, -23.990),
  zona('Praia_Grande', -46.420, -24.010),
  zona('Mongagua', -46.620, -24.090),
  zona('Itanhaem', -46.790, -24.180),
  zona('Peruibe', -47.000, -24.310),
  zona('Oceano_fundo', -46.300, -24.350),
  zona('Sao_Paulo_centro', -46.633, -23.550)
]);

var INI = '2018-07-01', FIM = '2026-01-01';

// ---------- coleções Sentinel-5P (nível 3), com máscara de nuvem
var no2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2').filterBounds(area).filterDate(INI, FIM)
  .map(function (im) {
    return im.select('tropospheric_NO2_column_number_density')
      .updateMask(im.select('cloud_fraction').lt(0.3)).rename('NO2')
      .copyProperties(im, ['system:time_start']);
  });
var so2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_SO2').filterBounds(area).filterDate(INI, FIM)
  .map(function (im) {
    return im.select('SO2_column_number_density')
      .updateMask(im.select('cloud_fraction').lt(0.3)).rename('SO2')
      .copyProperties(im, ['system:time_start']);
  });
var co = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_CO').filterBounds(area).filterDate(INI, FIM)
  .select(['CO_column_number_density'], ['CO']);

// ---------- aerossóis MODIS MAIAC (550 nm), só pixels sem nuvem, desde 2003
var aod = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES').filterBounds(area).filterDate('2003-01-01', FIM)
  .map(function (im) {
    var limpo = im.select('AOD_QA').bitwiseAnd(7).eq(1);
    return im.select('Optical_Depth_055').multiply(0.001).updateMask(limpo).rename('AOD')
      .copyProperties(im, ['system:time_start']);
  });

// ---------- meteorologia diária ERA5-Land
var era = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR').filterDate(INI, FIM)
  .select(['u_component_of_wind_10m', 'v_component_of_wind_10m', 'total_precipitation_sum', 'temperature_2m'],
          ['u10', 'v10', 'chuva', 't2m']);

// ---------- composição diária e mensal
function diario(col, banda, ini, fim) {
  var n = ee.Date(fim).difference(ee.Date(ini), 'day');
  return ee.ImageCollection(ee.List.sequence(0, n.subtract(1)).map(function (d) {
    var t0 = ee.Date(ini).advance(d, 'day'), c = col.filterDate(t0, t0.advance(1, 'day'));
    var im = ee.Image(ee.Algorithms.If(c.size().gt(0), c.select(banda).mean(),
                                       ee.Image.constant(0).updateMask(0).rename(banda)));
    return im.set('system:time_start', t0.millis());
  }));
}
function mensal(col, banda, ini, fim) {
  var n = ee.Date(fim).difference(ee.Date(ini), 'month').round();
  return ee.ImageCollection(ee.List.sequence(0, n.subtract(1)).map(function (m) {
    var t0 = ee.Date(ini).advance(m, 'month'), c = col.filterDate(t0, t0.advance(1, 'month'));
    var im = ee.Image(ee.Algorithms.If(c.size().gt(0), c.select(banda).mean(),
                                       ee.Image.constant(0).updateMask(0).rename(banda)));
    return im.set('system:time_start', t0.millis());
  }));
}
function tabela(colecao, banda, escala) {
  return colecao.map(function (im) {
    var data = ee.Date(im.get('system:time_start')).format('YYYY-MM-dd');
    return im.reduceRegions({collection: zonas, reducer: ee.Reducer.mean().combine(ee.Reducer.count(), '', true),
                             scale: escala})
      .map(function (f) { return f.set('data', data).set('variavel', banda).setGeometry(null); });
  }).flatten();
}

var opT = {folder: 'GEE_poluicao', fileFormat: 'CSV'};
Export.table.toDrive(Object.assign({collection: tabela(diario(no2, 'NO2', INI, FIM), 'NO2', 1113.2),
  description: 'no2_diario_zonas', selectors: ['zona', 'data', 'variavel', 'mean', 'count']}, opT));
Export.table.toDrive(Object.assign({collection: tabela(mensal(so2, 'SO2', INI, FIM), 'SO2', 1113.2),
  description: 'so2_mensal_zonas', selectors: ['zona', 'data', 'variavel', 'mean', 'count']}, opT));
Export.table.toDrive(Object.assign({collection: tabela(mensal(co, 'CO', INI, FIM), 'CO', 1113.2),
  description: 'co_mensal_zonas', selectors: ['zona', 'data', 'variavel', 'mean', 'count']}, opT));
Export.table.toDrive(Object.assign({collection: tabela(mensal(aod, 'AOD', '2003-01-01', FIM), 'AOD', 1000),
  description: 'aod_mensal_zonas', selectors: ['zona', 'data', 'variavel', 'mean', 'count']}, opT));
Export.table.toDrive(Object.assign({collection: era.map(function (im) {
    var data = ee.Date(im.get('system:time_start')).format('YYYY-MM-dd');
    return im.reduceRegions({collection: zonas, reducer: ee.Reducer.mean(), scale: 11132})
      .map(function (f) { return f.set('data', data).setGeometry(null); });
  }).flatten(), description: 'era5_diario_zonas', selectors: ['zona', 'data', 'u10', 'v10', 'chuva', 't2m']}, opT));

// ---------- mapas médios (GeoTIFF)
var opI = {folder: 'GEE_poluicao', region: area, crs: 'EPSG:4326', maxPixels: 1e10};
var no2_med = no2.filterDate('2019-01-01', FIM).mean();
var abrmai = function (ano) { return no2.filterDate(ano + '-04-01', ano + '-06-01').mean(); };
var no2_pandemia = abrmai(2020).subtract(ee.ImageCollection([abrmai(2019), abrmai(2021), abrmai(2022)]).mean()).rename('dNO2');
Export.image.toDrive(Object.assign({image: no2_med.toFloat(), description: 'no2_media_2019_2025', scale: 1113.2}, opI));
Export.image.toDrive(Object.assign({image: no2_pandemia.toFloat(), description: 'no2_abrmai2020_menos_referencia', scale: 1113.2}, opI));
Export.image.toDrive(Object.assign({image: so2.filterDate('2019-01-01', FIM).mean().toFloat(), description: 'so2_media_2019_2025', scale: 1113.2}, opI));
Export.image.toDrive(Object.assign({image: aod.filterDate('2019-01-01', FIM).mean().toFloat(), description: 'aod_media_2019_2025', scale: 1000}, opI));

// visualização rápida
Map.centerObject(area, 9);
Map.addLayer(no2_med, {min: 0, max: 0.00015, palette: ['black', 'blue', 'purple', 'cyan', 'green', 'yellow', 'red']}, 'NO2 médio');
Map.addLayer(zonas.style({color: 'white', fillColor: '00000000'}), {}, 'Zonas');
