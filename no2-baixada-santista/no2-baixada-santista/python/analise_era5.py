"""Normalização meteorológica com o ERA5 no horário da passagem (16 UTC): vento, camada limite (BLH), temperatura.
Uso: python analise_era5.py era5_16utc_zonas.csv pasta_pol"""
import sys, json, numpy as np, pandas as pd
E, O = sys.argv[1], sys.argv[2].rstrip('/') + '/'
e = pd.read_csv(E, parse_dates=['data']).rename(columns={'u_component_of_wind_10m': 'u13', 'v_component_of_wind_10m': 'v13',
    'boundary_layer_height': 'blh', 'temperature_2m': 't13'})
e['vel13'] = np.hypot(e.u13, e.v13); e['dir13'] = (np.degrees(np.arctan2(-e.u13, -e.v13)) + 360) % 360
x = pd.read_csv(O + 'no2_diario_met.csv', parse_dates=['data']).merge(e[['zona', 'data', 'u13', 'v13', 'vel13', 'dir13', 'blh', 't13']], on=['zona', 'data'])
Z = ['Sao_Paulo_centro', 'Cubatao_polo', 'Sao_Vicente', 'Santos_Porto', 'Praia_Grande', 'Guaruja', 'Mongagua', 'Itanhaem', 'Peruibe', 'Oceano_fundo']
def X(d, anos, v):
    if v == 1: cols = [np.ones(len(d)), d.vel, d.vel ** 2, d.u10, d.v10, np.log1p(d.chuva_mm), d.t2m - 290]
    else: cols = [np.ones(len(d)), d.vel13, d.vel13 ** 2, d.u13, d.v13, np.log(d.blh), np.log1p(d.chuva_mm), d.t13 - 290]
    cols += [(d.mes == m).astype(float) for m in range(2, 13)] + [(d.dow == k).astype(float) for k in range(1, 7)]
    cols += [(d.ano == a).astype(float) for a in anos]
    return np.column_stack(cols)
def ajusta(g, anos_tr, alvo, viz, v):
    tr = g[g.ano.isin(anos_tr)]; efe = anos_tr[1:]
    b, *_ = np.linalg.lstsq(X(tr, efe, v), tr.l, rcond=None); res = tr.l - X(tr, efe, v) @ b
    te = g[(g.ano == alvo) & g.mes.isin([4, 5])]; Xt = X(te, efe, v)
    vz = [a for a in viz if a in anos_tr]
    for j, a in enumerate(efe): Xt[:, -len(efe) + j] = (1 / len(vz)) if a in vz else 0
    pred = np.exp(Xt @ b + res.var() / 2)
    return (te.no2.mean() / pred.mean() - 1) * 100, 1 - res.var() / tr.l.var()
T = {}
for z in Z:
    g = x[(x.zona == z) & (x.no2 > 0)].dropna(subset=['vel', 'vel13', 'blh']).copy(); g['l'] = np.log(g.no2)
    for v in (1, 2):
        pan, r2 = ajusta(g, [a for a in range(2018, 2026) if a != 2020], 2020, (2019, 2021), v)
        errs = [ajusta(g, [a for a in range(2018, 2026) if a not in (2020, t)], t, (t - 1, t + 1), v)[0] for t in (2019, 2021, 2022, 2023, 2024)]
        T[(z, v)] = dict(normalizado=pan, r2=r2, erro_mediano=np.median(np.abs(errs)), erro_2021=errs[1], erros=[round(a, 1) for a in errs])
T = pd.DataFrame(T).T; print(T[['normalizado', 'r2', 'erro_mediano', 'erro_2021']].astype(float).round(2).unstack().to_string())
T.to_csv(O + 't_pandemia_era5.csv')
# camada limite e vento às 13h
R = {}
for z in ('Santos_Porto', 'Cubatao_polo', 'Sao_Paulo_centro', 'Praia_Grande'):
    g = x[x.zona == z].dropna(subset=['blh'])
    inv = g[g.mes.isin([6, 7, 8])]; ver = g[g.mes.isin([12, 1, 2])]
    cl = pd.cut(g.blh, [0, 500, 800, 1100, 1500, 5000]); 
    R[z] = dict(blh_inv=inv.blh.median(), blh_ver=ver.blh.median(),
                no2_blh={str(k): round(v, 1) for k, v in g.groupby(cl, observed=True).no2.mean().items()},
                n_blh={str(k): int(v) for k, v in g.groupby(cl, observed=True).size().items()},
                r_logno2_logblh=np.corrcoef(np.log(g.no2[g.no2 > 0]), np.log(g.blh[g.no2 > 0]))[0, 1],
                setor13={s: round(v, 1) for s, v in g.groupby(pd.cut(g.dir13, [0, 45, 90, 135, 180, 225, 270, 315, 360], labels=['NE', 'L', 'SE', 'S', 'SO', 'O', 'NO', 'N'])
                         , observed=True).no2.mean().items()},
                freq13={s: round(v * 100, 1) for s, v in g.groupby(pd.cut(g.dir13, [0, 45, 90, 135, 180, 225, 270, 315, 360], labels=['NE', 'L', 'SE', 'S', 'SO', 'O', 'NO', 'N'])
                         , observed=True).size().div(len(g)).items()},
                vel13_med=g.vel13.median(), r_vel_diaria_13h=np.corrcoef(g.vel, g.vel13)[0, 1])
    # quanto da razão inverno/verão some ao normalizar pela BLH (NO2*BLH ~ massa na camada)
    R[z]['razao_inv_ver'] = inv.no2.mean() / ver.no2.mean()
    R[z]['razao_inv_ver_blh_fixa'] = None
json.dump(R, open(O + 'resultados_era5.json', 'w'), indent=1, default=float, ensure_ascii=False)
for k, v in R.items(): print(k, json.dumps(v, default=float, ensure_ascii=False))
