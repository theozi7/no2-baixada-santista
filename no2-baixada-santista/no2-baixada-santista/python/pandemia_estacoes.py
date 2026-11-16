"""Normalização meteorológica do NO2 das estações (12h-14h) com ERA5 às 13h. Uso: python pandemia_estacoes.py pasta_pol era5.csv"""
import sys, numpy as np, pandas as pd
O, E = sys.argv[1].rstrip('/') + '/', sys.argv[2]
P = pd.read_pickle(O + 'cetesb_pares.pkl')
e = pd.read_csv(E, parse_dates=['data']).rename(columns={'u_component_of_wind_10m': 'u13', 'v_component_of_wind_10m': 'v13', 'boundary_layer_height': 'blh', 'temperature_2m': 't13'})
x = pd.read_csv(O + 'no2_diario_met.csv', parse_dates=['data'])
mapa = {'Santos-Ponta da Praia': 'Santos_Porto', 'Cubatão-Centro': 'Cubatao_polo'}
for est, z in mapa.items():
    s = P[est]  # tem só dias com TROPOMI; refazer com todos os dias da estação
for est, z in mapa.items():
    import glob
    pass
# séries completas das estações (todos os dias, não só os pareados)
h = []
for f in glob.glob('../cetesb/*.csv'):
    L = open(f, encoding='latin1').read().split('\n'); nm = L[4].split(';')[1].strip()
    d = pd.read_csv(f, sep=';', skiprows=8, header=None, names=['data', 'hora', 'v'], encoding='latin1'); d['est'] = nm; h.append(d)
h = pd.concat(h).dropna(subset=['v']).drop_duplicates(['est', 'data', 'hora'])
h['v'] = pd.to_numeric(h.v.astype(str).str.replace(',', '.'), errors='coerce'); h['hr'] = h.hora.str[:2].astype(int)
h['dia'] = pd.to_datetime(h.data, format='%d/%m/%Y')
p = h[h.hr.isin([13, 14])].groupby(['est', 'dia']).v.agg(['mean', 'count']); p = p[p['count'] == 2]['mean'].rename('no2').reset_index()
def X(d, anos):
    cols = [np.ones(len(d)), d.vel13, d.vel13 ** 2, d.u13, d.v13, np.log(d.blh), np.log1p(d.chuva_mm), d.t13 - 290]
    cols += [(d.mes == m).astype(float) for m in range(2, 13)] + [(d.dow == k).astype(float) for k in range(1, 7)]
    cols += [(d.ano == a).astype(float) for a in anos]
    return np.column_stack(cols)
for est, z in mapa.items():
    g = p[p.est == est].rename(columns={'dia': 'data'})
    met = x[x.zona == z][['data', 'chuva_mm']].merge(e[e.zona == z], on='data')
    g = g.merge(met, on='data'); g = g[g.no2 > 0].copy()
    g['vel13'] = np.hypot(g.u13, g.v13); g['ano'] = g.data.dt.year; g['mes'] = g.data.dt.month; g['dow'] = g.data.dt.dayofweek; g['l'] = np.log(g.no2)
    def aj(anos_tr, alvo, viz):
        tr = g[g.ano.isin(anos_tr)]; efe = anos_tr[1:]
        b, *_ = np.linalg.lstsq(X(tr, efe), tr.l, rcond=None); res = tr.l - X(tr, efe) @ b
        te = g[(g.ano == alvo) & g.mes.isin([4, 5])]; Xt = X(te, efe); vz = [a for a in viz if a in anos_tr]
        for j, a in enumerate(efe): Xt[:, -len(efe) + j] = (1 / len(vz)) if a in vz else 0
        return (te.no2.mean() / np.exp(Xt @ b + res.var() / 2).mean() - 1) * 100, 1 - res.var() / tr.l.var(), len(te)
    anos = sorted(g.ano.unique())
    pan, r2, n = aj([a for a in anos if a != 2020], 2020, (2019, 2021))
    errs = [aj([a for a in anos if a not in (2020, t)], t, (t - 1, t + 1))[0] for t in anos if t not in (2020,) and t in (2019, 2021, 2022, 2023, 2024)]
    am = g[g.mes.isin([4, 5])]; bruto = (am[am.ano == 2020].no2.mean() / am[am.ano.isin([2019, 2021, 2022])].no2.mean() - 1) * 100
    print(est, 'anos', anos, 'bruto %.1f normalizado %.1f R2 %.2f n2020 %d erro_mediano %.1f erros %s' % (bruto, pan, r2, n, np.median(np.abs(errs)), np.round(errs, 1)))
