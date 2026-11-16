"""Análises do artigo de poluição do ar na Baixada Santista (Sentinel-5P, MAIAC, ERA5-Land).
Uso: python analise_pol.py pasta_GEE pasta_saida"""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__))); from mk_hr import mk_hr
from scipy import stats
D, O = sys.argv[1].rstrip('/') + '/', sys.argv[2].rstrip('/') + '/'
R = {}
n = pd.read_csv(D + 'no2_diario_zonas.csv', parse_dates=['data']).rename(columns={'mean': 'v'})
n = n[n['count'] >= 3]                                   # ao menos 3 pixels válidos no círculo
n['no2'] = n.v * 1e6                                      # µmol/m²
e = pd.read_csv(D + 'era5_diario_zonas.csv', parse_dates=['data'])
# zona oceânica sem ERA5-Land: usa Guarujá (mais próximo)
# Guarujá e oceano sem ERA5-Land (pixels costeiros mascarados): usa a meteorologia de Santos_Porto
extra = []
for zz in ('Guaruja', 'Oceano_fundo'):
    eo = e[e.zona == 'Santos_Porto'].copy(); eo['zona'] = zz; extra.append(eo)
e = pd.concat([e[~e.zona.isin(['Guaruja', 'Oceano_fundo'])]] + extra)
e['vel'] = np.hypot(e.u10, e.v10); e['chuva_mm'] = e.chuva * 1000
e['dir'] = (np.degrees(np.arctan2(-e.u10, -e.v10)) + 360) % 360          # direção de onde vem o vento
x = n.merge(e[['zona', 'data', 'vel', 'u10', 'v10', 'chuva_mm', 'dir', 't2m']], on=['zona', 'data'], how='left')
x = x[x.data >= '2018-07-01']
x['ano'] = x.data.dt.year; x['mes'] = x.data.dt.month; x['dow'] = x.data.dt.dayofweek
x.to_csv(O + 'no2_diario_met.csv', index=False)
Z = ['Sao_Paulo_centro', 'Cubatao_polo', 'Sao_Vicente', 'Santos_Porto', 'Praia_Grande', 'Guaruja', 'Mongagua', 'Itanhaem', 'Peruibe', 'Oceano_fundo']

# 1. médias (2019–2025, anos completos) e dias válidos
y = x[(x.ano >= 2019) & (x.ano <= 2025)]
T1 = y.groupby('zona').no2.agg(media='mean', mediana='median', p90=lambda s: s.quantile(0.9), dias='size').reindex(Z)
T1['razao_fundo'] = T1.media / T1.loc['Oceano_fundo', 'media']
T1.round(2).to_csv(O + 't_medias.csv'); R['medias'] = T1.round(1).to_dict()

# 2. ciclo sazonal (média mensal 2019–2025) e razão inverno/verão
S = y.groupby(['zona', 'mes']).no2.mean().unstack(0)[Z]; S.to_csv(O + 'sazonal.csv')
R['inv_ver'] = (S.loc[[6, 7, 8]].mean() / S.loc[[12, 1, 2]].mean()).round(2).to_dict()

# 3. tendência: médias mensais dessazonalizadas (anomalia relativa) 2019–2025, MK com correção HR
M = y.groupby(['zona', pd.Grouper(key='data', freq='MS')]).no2.mean().unstack(0)[Z]
anom = M.groupby(M.index.month).transform(lambda s: s / s.mean() - 1) * 100
TR = {}
for z in Z:
    s = anom[z].dropna(); t = s.index.year + (s.index.month - 0.5) / 12
    sen, p0, p1 = mk_hr(s.values, t.values); TR[z] = dict(sen_pct_ano=sen / 10, p=p1)
TR = pd.DataFrame(TR).T; TR.to_csv(O + 't_tendencia.csv'); R['tendencia'] = TR.round(3).to_dict()
M.to_csv(O + 'no2_mensal.csv')
A = y.groupby(['zona', 'ano']).no2.mean().unstack(0)[Z]; A.to_csv(O + 'no2_anual.csv')

# 4. ciclo semanal: razão de cada dia pela média semanal (por zona), com IC 95% por bootstrap de semanas
W = {}
rng = np.random.default_rng(1)
for z in Z:
    g = y[y.zona == z].copy(); g['semana'] = g.data.dt.to_period('W')
    base = g.set_index('data').no2.rolling('28D', center=True, min_periods=8).mean().values
    g['rel'] = g.no2 / base
    dw = g.groupby('dow').rel.mean(); wk = dw.loc[0:4].mean()
    # bootstrap por semanas: soma e contagem de rel por (semana, dia útil/domingo)
    g['tipo'] = np.where(g.dow == 6, 'dom', np.where(g.dow <= 4, 'util', 'sab'))
    ag = g.groupby(['semana', 'tipo']).rel.agg(['sum', 'count']).unstack(fill_value=0)
    S_, C_ = ag['sum'], ag['count']; nW = len(ag); boots = []
    for _ in range(1000):
        k = rng.integers(0, nW, nW)
        dom = S_['dom'].values[k].sum() / max(C_['dom'].values[k].sum(), 1)
        ut = S_['util'].values[k].sum() / max(C_['util'].values[k].sum(), 1)
        boots.append(dom / ut)
    W[z] = dict(**{f'd{k}': dw.get(k, np.nan) / wk for k in range(7)}, dom_util=dw.get(6) / wk,
                sab_util=dw.get(5) / wk, ic_lo=np.nanpercentile(boots, 2.5), ic_hi=np.nanpercentile(boots, 97.5))
W = pd.DataFrame(W).T; W.to_csv(O + 't_semanal.csv'); R['semanal'] = W[['dom_util', 'sab_util', 'ic_lo', 'ic_hi']].round(2).to_dict()

# 5. vento: NO2 por classe de velocidade e por setor de direção (Santos_Porto e Cubatão)
for z in ('Santos_Porto', 'Cubatao_polo', 'Praia_Grande'):
    g = y[y.zona == z].dropna(subset=['vel'])
    g['cls'] = pd.cut(g.vel, [0, 1, 2, 3, 4, 10]); g['setor'] = pd.cut(g.dir, [0, 45, 90, 135, 180, 225, 270, 315, 360],
                                                                     labels=['N', 'NE', 'E', 'SE', 'S', 'SO', 'O', 'NO'], include_lowest=True)
    R[f'vento_{z}'] = dict(vel={str(k): float(v) for k, v in g.groupby('cls', observed=True).no2.mean().round(1).items()},
                           setor={str(k): [float(a), int(b)] for k, (a, b) in g.groupby('setor', observed=True).no2.agg(['mean', 'size']).round(1).iterrows()},
                           chuva_sim_nao=[round(g[g.chuva_mm >= 1].no2.mean(), 1), round(g[g.chuva_mm < 1].no2.mean(), 1)])
R['vento'] = {k: {str(a): b for a, b in v['vel'].items()} for k, v in R.items() if k.startswith('vento_')}

# 6. pandemia: modelo de normalização meteorológica (regressão log-linear treinada sem 2020)
import itertools
PAN = {}
for z in Z:
    g = x[(x.zona == z)].dropna(subset=['vel', 'no2']).copy(); g = g[g.no2 > 0]
    g['lno2'] = np.log(g.no2)
    def X(d):
        cols = [np.ones(len(d)), d.vel, d.vel ** 2, d.u10, d.v10, np.log1p(d.chuva_mm), d.t2m - 290]
        cols += [(d.mes == m).astype(float) for m in range(2, 13)] + [(d.dow == k).astype(float) for k in range(1, 7)]
        cols += [(d.ano == a).astype(float) for a in (2019, 2021, 2022, 2023, 2024, 2025)]
        return np.column_stack(cols)
    tr = g[g.ano != 2020]; b, *_ = np.linalg.lstsq(X(tr), tr.lno2, rcond=None)
    res = tr.lno2 - X(tr) @ b; sig2 = res.var()
    p = g[(g.ano == 2020) & g.mes.isin([4, 5])].copy()
    # previsão contrafactual para 2020: média dos efeitos de ano 2019 e 2021
    Xp = X(p); i19, i21 = Xp.shape[1] - 6, Xp.shape[1] - 5
    Xp[:, i19] = 0.5; Xp[:, i21] = 0.5
    pred = np.exp(Xp @ b + sig2 / 2)
    ref = g[(g.ano.isin([2019, 2021, 2022])) & g.mes.isin([4, 5])]
    PAN[z] = dict(bruto_pct=(p.no2.mean() / ref.no2.mean() - 1) * 100, normalizado_pct=(p.no2.mean() / pred.mean() - 1) * 100,
                  r2=1 - res.var() / tr.lno2.var(), n2020=len(p))
PAN = pd.DataFrame(PAN).T; PAN.to_csv(O + 't_pandemia.csv'); R['pandemia'] = PAN.round(2).to_dict()

# 7. AOD 2003–2025: média, sazonalidade e tendência
a = pd.read_csv(D + 'aod_mensal_zonas.csv', parse_dates=['data']).rename(columns={'mean': 'aod'})
AM = a.pivot(index='data', columns='zona', values='aod')[Z]; AM.to_csv(O + 'aod_mensal.csv')
R['aod_media'] = AM.mean().round(3).to_dict(); R['aod_mes'] = AM.groupby(AM.index.month).mean().mean(axis=1).round(3).to_dict()
TA = {}
for z in Z:
    an_ = AM[z].groupby(AM.index.year).mean(); s, p0, p1 = mk_hr(an_.values, an_.index.values); TA[z] = dict(sen_dec=s, p=p1)
R['aod_tend'] = pd.DataFrame(TA).T.round(4).to_dict()
# 8. SO2 e CO
for f in ('so2', 'co'):
    q = pd.read_csv(D + f + '_mensal_zonas.csv', parse_dates=['data'])
    R[f] = q.groupby('zona')['mean'].agg(['mean', 'std']).reindex(Z).round(6).to_dict()
json.dump(R, open(O + 'resultados.json', 'w'), indent=1, default=float)
print(json.dumps({k: R[k] for k in ('medias', 'inv_ver', 'tendencia', 'semanal', 'pandemia', 'vento', 'aod_media', 'aod_mes', 'aod_tend')}, indent=0, default=float)[:6000])
for k in ('vento_Santos_Porto', 'vento_Cubatao_polo'): print(k, R[k]['setor'], R[k]['chuva_sim_nao'])
