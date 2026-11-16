# NO₂ na Baixada Santista por satélite (2018–2025)

Dados processados e código do artigo **“Poluição do Ar na Baixada Santista Vista do Espaço: Dióxido de Nitrogênio entre o Polo Industrial de Cubatão, o Porto de Santos e o Litoral Sul de São Paulo (2018–2025)”**, de Theo Goulart Tavares de Lima Araujo (Universidade Estadual de Londrina).

O trabalho usa o NO₂ troposférico do TROPOMI (Sentinel-5P), a espessura óptica de aerossóis do MAIAC (MODIS), as reanálises ERA5 e ERA5-Land e as estações automáticas da CETESB (sistema QUALAR).

## Estrutura

| Pasta | Conteúdo |
|---|---|
| `gee/` | Scripts do Google Earth Engine que exportam as séries por zona e os mapas médios |
| `python/` | Análises estatísticas e figuras |
| `dados/gee/` | Arquivos exportados do Google Earth Engine (CSV e GeoTIFF) |
| `dados/processados/` | Tabelas e resultados usados no artigo |
| `dados/mapas/` | Limites municipais, países e relevo (SRTM) usados nos mapas |
| `dados/cetesb/` | Local para os dados horários da CETESB (não redistribuídos; ver abaixo) |
| `figuras/` | Figuras do artigo (300 dpi) |

## Como reproduzir

1. Rode `gee/poluicao_baixada_GEE.js` e `gee/era5_13h_GEE.js` no Code Editor do Google Earth Engine e salve as exportações em `dados/gee/`.
2. Baixe no QUALAR (https://qualar.cetesb.sp.gov.br) o NO₂ horário das estações Santos-Ponta da Praia e Cubatão-Centro (2018–2025) e salve os CSV em `dados/cetesb/`.
3. A partir da pasta `dados/processados/`, rode os scripts na ordem:

```
python ../../python/analise_pol.py ../gee .
python ../../python/figs_pol.py ../gee ../mapas .
python ../../python/analise_extra.py ../gee .
python ../../python/analise_era5.py ../gee/era5_16utc_zonas.csv .
python ../../python/valida_cetesb.py ../cetesb .
python ../../python/pandemia_estacoes.py . ../gee/era5_16utc_zonas.csv
python ../../python/revisao_stats.py
python ../../python/diagnostico_modelos.py
python ../../python/fig_cetesb.py
python ../../python/fig_fundeio.py
```

Dependências: Python 3.10 ou superior, `numpy`, `pandas`, `scipy`, `matplotlib` e `Pillow`.

## Licença

Código sob a licença MIT. Dados processados e figuras sob a licença Creative Commons Atribuição 4.0 Internacional (CC BY 4.0). Os dados de origem pertencem às respectivas instituições (Copernicus/ESA, ECMWF, NASA, CETESB, FAO e USGS/NASA SRTM).

## Como citar

Ver o arquivo `CITATION.cff` ou o DOI do Zenodo indicado na página do repositório.
