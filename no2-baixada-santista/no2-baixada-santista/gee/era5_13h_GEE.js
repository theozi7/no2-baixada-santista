// ERA5 horário no horário da passagem do TROPOMI (16 UTC ≈ 13h local): vento a 10 m e altura da camada limite (BLH)
// para as dez zonas do artigo, 2018-07 a 2025-12. Exporta um CSV para o Google Drive (pasta GEE_poluicao).
var zonas = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([-46.633, -23.550]).buffer(4000), {zona: 'Sao_Paulo_centro'}),
  ee.Feature(ee.Geometry.Point([-46.395, -23.865]).buffer(4000), {zona: 'Cubatao_polo'}),
  ee.Feature(ee.Geometry.Point([-46.390, -23.960]).buffer(4000), {zona: 'Sao_Vicente'}),
  ee.Feature(ee.Geometry.Point([-46.315, -23.950]).buffer(4000), {zona: 'Santos_Porto'}),
  ee.Feature(ee.Geometry.Point([-46.420, -24.010]).buffer(4000), {zona: 'Praia_Grande'}),
  ee.Feature(ee.Geometry.Point([-46.250, -23.990]).buffer(4000), {zona: 'Guaruja'}),
  ee.Feature(ee.Geometry.Point([-46.620, -24.090]).buffer(4000), {zona: 'Mongagua'}),
  ee.Feature(ee.Geometry.Point([-46.790, -24.180]).buffer(4000), {zona: 'Itanhaem'}),
  ee.Feature(ee.Geometry.Point([-47.000, -24.310]).buffer(4000), {zona: 'Peruibe'}),
  ee.Feature(ee.Geometry.Point([-46.300, -24.350]).buffer(4000), {zona: 'Oceano_fundo'})
]);
// Mesmas coordenadas do script original (poluicao_baixada_GEE.js).
var era5 = ee.ImageCollection('ECMWF/ERA5/HOURLY')   // ERA5 completo (31 km): tem dados sobre o mar e a BLH
  .filterDate('2018-07-01', '2026-01-01')
  .filter(ee.Filter.calendarRange(16, 16, 'hour'))
  .select(['u_component_of_wind_10m', 'v_component_of_wind_10m', 'boundary_layer_height', 'temperature_2m', 'total_precipitation']);
var tab = era5.map(function (im) {
  return im.reduceRegions({collection: zonas, reducer: ee.Reducer.mean(), scale: 27830})
    .map(function (f) { return f.set('data', im.date().format('YYYY-MM-dd')).setGeometry(null); });
}).flatten();
Export.table.toDrive({collection: tab, description: 'era5_16utc_zonas', folder: 'GEE_poluicao', fileFormat: 'CSV',
  selectors: ['zona', 'data', 'u_component_of_wind_10m', 'v_component_of_wind_10m', 'boundary_layer_height', 'temperature_2m', 'total_precipitation']});
