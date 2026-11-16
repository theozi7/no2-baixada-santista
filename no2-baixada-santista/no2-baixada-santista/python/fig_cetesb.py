import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Carlito', 'Liberation Sans'], 'font.size': 8.5,
                     'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 400, 'axes.linewidth': 0.6})
P = pd.read_pickle('cetesb_pares.pkl')
C = {'Santos-Ponta da Praia': '#e08e0b', 'Cubatão-Centro': '#c0392b'}
fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.6), gridspec_kw={'width_ratios': [1, 1.5, 1]})
for est, m in P.items():
    r = stats.pearsonr(m.trop, m.st13)[0]
    ax[0].scatter(m.trop, m.st13, s=3, alpha=0.25, color=C[est], lw=0, label=f'{est} (r = {r:.2f})'.replace('.', ','))
    mm = m[['trop', 'st13']].resample('MS').mean(); n = m.resample('MS').size(); mm = mm.where(n >= 5).asfreq('MS')
    ax[1].plot(mm.index, mm.st13, color=C[est], lw=1.1, label=f'{est}, estação')
    ax2 = ax[1].twinx() if est.startswith('Santos') else ax2
    ax2.plot(mm.index, mm.trop, color=C[est], lw=0.9, ls='--', label=f'{est}, TROPOMI')
    st = m.st24.dropna(); t = m.trop
    ax[2].plot(range(1, 13), st.groupby(st.index.month).mean() / st.mean(), color=C[est], lw=1.2, marker='o', ms=2.5)
    ax[2].plot(range(1, 13), t.groupby(t.index.month).mean() / t.mean(), color=C[est], lw=1, ls='--', marker='s', ms=2.5)
ax[0].set_xlabel('TROPOMI (µmol m$^{-2}$)'); ax[0].set_ylabel('Estação, 12h–14h (µg m$^{-3}$)'); ax[0].set_xlim(0, 500); ax[0].set_ylim(0, 150)
ax[0].legend(fontsize=6.5, markerscale=3, frameon=False, loc='upper left'); ax[0].set_title('(a)', loc='left', fontweight='bold')
ax[1].set_ylabel('Estação (µg m$^{-3}$), linha cheia'); ax2.set_ylabel('TROPOMI (µmol m$^{-2}$), tracejado'); ax2.spines['right'].set_visible(True)
ax[1].set_title('(b)', loc='left', fontweight='bold'); ax[1].tick_params(axis='x', labelsize=7)
ax[2].set_xticks(range(1, 13)); ax[2].set_xticklabels('JFMAMJJASOND'); ax[2].axhline(1, color='0.6', lw=0.5)
ax[2].set_ylabel('Média mensal / média anual'); ax[2].set_title('(c)', loc='left', fontweight='bold')
ax[2].plot([], [], color='k', lw=1.2, label='Estação (24 h)'); ax[2].plot([], [], color='k', ls='--', label='TROPOMI'); ax[2].legend(fontsize=6.5, frameon=False)
fig.tight_layout(); fig.savefig('fig9_cetesb.png')
