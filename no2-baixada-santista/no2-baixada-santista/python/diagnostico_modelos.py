import numpy as np, pandas as pd
from scipy import stats
e = pd.read_csv('../gee/era5_16utc_zonas.csv', parse_dates=['data']).rename(columns={'u_component_of_wind_10m': 'u13', 'v_component_of_wind_10m': 'v13', 'boundary_layer_height': 'blh', 'temperature_2m': 't13'})
e['vel13'] = np.hypot(e.u13, e.v13)
x = pd.read_csv('no2_diario_met.csv', parse_dates=['data']).merge(e[['zona', 'data', 'u13', 'v13', 'vel13', 'blh', 't13']], on=['zona', 'data'])
MET = {1: lambda d: pd.DataFrame({'vel': d.vel, 'u': d.u10, 'v': d.v10, 'chuva': np.log1p(d.chuva_mm), 't': d.t2m}),
       2: lambda d: pd.DataFrame({'vel': d.vel13, 'u': d.u13, 'v': d.v13, 'lblh': np.log(d.blh), 'chuva': np.log1p(d.chuva_mm), 't': d.t13})}
def X(d, anos, v):
    m = MET[v](d); vc = m.vel - m.vel.mean()
    cols = [np.ones(len(d)), m.vel, vc ** 2] + [m[c] for c in m.columns if c != 'vel']
    cols += [(d.mes == k).astype(float) for k in range(2, 13)] + [(d.dow == k).astype(float) for k in range(1, 7)] + [(d.ano == a).astype(float) for a in anos]
    return np.column_stack(cols)
def vif(M):
    M = (M - M.mean()) / M.std(); out = {}
    for c in M.columns:
        o = M.drop(columns=c); b, *_ = np.linalg.lstsq(np.column_stack([np.ones(len(o)), o]), M[c], rcond=None)
        r2 = 1 - ((M[c] - np.column_stack([np.ones(len(o)), o]) @ b) ** 2).sum() / ((M[c] - M[c].mean()) ** 2).sum(); out[c] = 1 / (1 - r2)
    return out
Z = ['Sao_Paulo_centro', 'Cubatao_polo', 'Sao_Vicente', 'Santos_Porto', 'Praia_Grande', 'Guaruja', 'Mongagua', 'Itanhaem', 'Peruibe', 'Oceano_fundo']
rows = []
for z in Z:
    g = x[(x.zona == z) & (x.no2 > 0)].dropna(subset=['vel', 'vel13', 'blh']).sort_values('data').copy(); g['l'] = np.log(g.no2)
    for v in (1, 2):
        anos = [a for a in range(2018, 2026) if a != 2020]; tr = g[g.ano.isin(anos)]
        Xm = X(tr, anos[1:], v); b, *_ = np.linalg.lstsq(Xm, tr.l, rcond=None); fit = Xm @ b; res = tr.l - fit
        cons = tr.data.diff().dt.days.eq(1).values[1:]
        ac1 = np.corrcoef(res.values[:-1][cons], res.values[1:][cons])[0, 1]
        bp = stats.pearsonr(fit, res ** 2)[0]
        vv = vif(MET[v](tr))
        rows.append(dict(zona=z, modelo='AB'[v - 1], n=len(tr), r2=1 - res.var() / tr.l.var(), rmse_log=res.std(), ac1=ac1,
                         assim=stats.skew(res), curt=stats.kurtosis(res), het=bp, vif_max=max(vv.values()), vif_arg=max(vv, key=vv.get)))
D = pd.DataFrame(rows); D.to_csv('t_diagnostico.csv', index=False)
print(D.round(2).to_string())
# sensibilidade à retransformação: fator de Duan (smearing) x exp(s2/2)
print('zona modelo  red_lognormal  red_duan')
for z in Z:
    g = x[(x.zona == z) & (x.no2 > 0)].dropna(subset=['vel', 'vel13', 'blh']).copy(); g['l'] = np.log(g.no2)
    for v in (1, 2):
        anos = [a for a in range(2018, 2026) if a != 2020]; efe = anos[1:]; tr = g[g.ano.isin(anos)]
        Xm = X(tr, efe, v); b, *_ = np.linalg.lstsq(Xm, tr.l, rcond=None); res = tr.l - Xm @ b
        te = g[(g.ano == 2020) & g.mes.isin([4, 5])]; Xt = X(te, efe, v)
        for j, a in enumerate(efe): Xt[:, -len(efe) + j] = 0.5 if a in (2019, 2021) else 0
        p1 = np.exp(Xt @ b + res.var() / 2).mean(); p2 = (np.exp(Xt @ b) * np.exp(res).mean()).mean()
        print(z, 'AB'[v-1], round((te.no2.mean()/p1-1)*100,1), round((te.no2.mean()/p2-1)*100,1))
