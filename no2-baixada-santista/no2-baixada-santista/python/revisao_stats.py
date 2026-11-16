import numpy as np, pandas as pd, json
from scipy import stats
from PIL import Image
R = {}
# (2) métricas da validação
P = pd.read_pickle('cetesb_pares.pkl')
for est, m in P.items():
    lr = stats.linregress(m.trop, m.st13); pred = lr.intercept + lr.slope * m.trop; res = m.st13 - pred
    mm = m[['trop', 'st13']].resample('MS').mean(); n = m.resample('MS').size(); mm = mm[n >= 5].dropna()
    lm = stats.linregress(mm.trop, mm.st13); rm = mm.st13 - (lm.intercept + lm.slope * mm.trop)
    # viés relativo: razão entre médias normalizadas não faz sentido; usar viés da amplitude relativa (CV)
    R[est] = dict(n=len(m), slope=lr.slope, slope_se=lr.stderr, intercept=lr.intercept, r=lr.rvalue, p=lr.pvalue,
                  rmse=np.sqrt((res ** 2).mean()), mae=res.abs().mean(), media_st=m.st13.mean(), media_trop=m.trop.mean(),
                  n_m=len(mm), slope_m=lm.slope, intercept_m=lm.intercept, r_m=lm.rvalue, rmse_m=np.sqrt((rm ** 2).mean()), mae_m=rm.abs().mean(),
                  cv_st=m.st13.std() / m.st13.mean(), cv_trop=m.trop.std() / m.trop.mean())
# (3) teste do raio no mapa médio
im = Image.open('../gee/no2_media_2019_2025.tif'); A = np.array(im, float) * 1e6
sx, sy = im.tag_v2[33550][:2]; x0, y0 = im.tag_v2[33922][3:5]
co = {'Sao_Paulo_centro': (-46.633, -23.550), 'Cubatao_polo': (-46.395, -23.865), 'Sao_Vicente': (-46.390, -23.960), 'Santos_Porto': (-46.315, -23.950),
      'Praia_Grande': (-46.420, -24.010), 'Guaruja': (-46.250, -23.990), 'Mongagua': (-46.620, -24.090), 'Itanhaem': (-46.790, -24.180), 'Peruibe': (-47.000, -24.310), 'Oceano_fundo': (-46.300, -24.350)}
ny, nx = A.shape; lon = x0 + (np.arange(nx) + 0.5) * sx; lat = y0 - (np.arange(ny) + 0.5) * sy; LO, LA = np.meshgrid(lon, lat)
rad = {}
for rk in (2, 4, 6, 8):
    v = {}
    for z, (a, b) in co.items():
        d = np.hypot((LO - a) * 111.32 * np.cos(np.radians(b)), (LA - b) * 111.32); v[z] = np.nanmean(A[d <= rk])
    rad[rk] = {z: v[z] / v['Oceano_fundo'] for z in v}
RAD = pd.DataFrame(rad); RAD.to_csv('t_raio.csv'); R['raio'] = RAD.round(2).to_dict(); R['pixel_m'] = [sx * 111320, sy * 111320]
print(RAD.round(2).to_string())
from scipy.stats import spearmanr
print('spearman ordem 2 vs 8 km', spearmanr(RAD[2], RAD[8])[0])
# (4) bootstrap dos dias de 2020 para os modelos A e B
e = pd.read_csv('../gee/era5_16utc_zonas.csv', parse_dates=['data']).rename(columns={'u_component_of_wind_10m': 'u13', 'v_component_of_wind_10m': 'v13', 'boundary_layer_height': 'blh', 'temperature_2m': 't13'})
e['vel13'] = np.hypot(e.u13, e.v13)
x = pd.read_csv('no2_diario_met.csv', parse_dates=['data']).merge(e[['zona', 'data', 'u13', 'v13', 'vel13', 'blh', 't13']], on=['zona', 'data'])
def X(d, anos, v):
    cols = [np.ones(len(d)), d.vel, d.vel ** 2, d.u10, d.v10, np.log1p(d.chuva_mm), d.t2m - 290] if v == 1 else \
           [np.ones(len(d)), d.vel13, d.vel13 ** 2, d.u13, d.v13, np.log(d.blh), np.log1p(d.chuva_mm), d.t13 - 290]
    cols += [(d.mes == m).astype(float) for m in range(2, 13)] + [(d.dow == k).astype(float) for k in range(1, 7)] + [(d.ano == a).astype(float) for a in anos]
    return np.column_stack(cols)
rng = np.random.default_rng(1); B = {}
for z in co:
    g = x[(x.zona == z) & (x.no2 > 0)].dropna(subset=['vel', 'vel13', 'blh']).copy(); g['l'] = np.log(g.no2)
    for v in (1, 2):
        anos = [a for a in range(2018, 2026) if a != 2020]; efe = anos[1:]; tr = g[g.ano.isin(anos)]
        b, *_ = np.linalg.lstsq(X(tr, efe, v), tr.l, rcond=None); s2 = (tr.l - X(tr, efe, v) @ b).var()
        te = g[(g.ano == 2020) & g.mes.isin([4, 5])]; Xt = X(te, efe, v)
        for j, a in enumerate(efe): Xt[:, -len(efe) + j] = 0.5 if a in (2019, 2021) else 0
        pr = np.exp(Xt @ b + s2 / 2); ob = te.no2.values
        bs = [(ob[i].mean() / pr[i].mean() - 1) * 100 for i in (rng.integers(0, len(ob), len(ob)) for _ in range(2000))]
        B[(z, v)] = dict(est=(ob.mean() / pr.mean() - 1) * 100, lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5))
B = pd.DataFrame(B).T; B.to_csv('t_pandemia_ic.csv'); print(B.round(1).unstack().to_string())
json.dump({k: v for k, v in R.items()}, open('resultados_revisao.json', 'w'), indent=1, default=float, ensure_ascii=False)
for k in P: print(k, {a: round(b, 3) for a, b in R[k].items()})
print(R['pixel_m'])
