import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, UnivariateSpline, InterpolatedUnivariateSpline
from mectools.plotter import plotter

plt = plotter(backend='Agg')
import matplotlib as mpl
colors = mpl.rcParams['axes.prop_cycle'].by_key()['color']

##
## options
##

ptag = '1_2'
mtag = '15_m3'

##
## Tools
##

def read_split(fname, delim='='):
    with open(fname) as fid:
        pairs = [line.split(delim, 1) for line in fid if delim in line]
        return pd.Series({k.strip(): float(v.strip()) for k, v in pairs})

##
## general names
##

names_bri = [f'Basic Intensive $m = {i}$' if i < 8 else f'Basic Intensive $m >= {i}$' for i in range(1, 9)]
names_epb = [f'Basic Extensive $m = {i}$' if i < 8 else f'Basic Extensive $m >= {i}$' for i in range(1, 9)]

##
## Federal R&D Funding
##

params = read_split(f'params/params{ptag}.txt', delim=':')
eta = params['eta'] # 0.079363351762
alpha = params['alpha'] # 0.049674661161

def fill_stats(df):
    df['patents'] = (df['hot']*eta+(1-df['hot'])*alpha)*df['taua']
    df['applied'] = df['wage']*df['Lapp']
    df['pats_per_researcher'] = df['patents']/df['Lapp']
    df['pats_per_rnd'] = df['patents']/df['applied']
    return df

sim_varw = (100*pd.read_csv('robust/academic_gdp_varw.csv')).set_index('ac_gdp')
sim_fixw = (100*pd.read_csv('robust/academic_gdp_fixw.csv')).set_index('ac_gdp')
dat = pd.read_excel('robust/RD.GDP.xlsx', sheetname='raw', index_col='Year')

sim_varw = fill_stats(sim_varw)
sim_fixw = fill_stats(sim_fixw)
dat['total_rnd'] = 100*dat['Total_RD_Frac']

# interpolate
f_ppres = interp1d(sim_varw.index, sim_varw['pats_per_researcher'])
f_pprnd = interp1d(sim_varw.index, sim_varw['pats_per_rnd'])

total_rnd = 100*dat['Total_RD_Frac']
dat['pred_ppres'] = f_ppres(dat['total_rnd'])
dat['pred_pprnd'] = f_pprnd(dat['total_rnd'])

# plotting
fig, ax = plt.subplots(figsize=(5, 4))
sim_varw['pats_per_researcher'].plot(ax=ax)
ax.set_xlabel('Academic / GDP (%)')
ax.set_ylabel('Patents / Researcher')
fig.savefig('output/graphics/academic_gdp_theory.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
dat['total_rnd'].plot(ax=ax)
ax.set_xlabel('Year')
ax.set_ylabel('Federal R&D / GDP (%)')
fig.savefig('output/graphics/academic_gdp_data.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
dat['pred_ppres'].plot(ax=ax)
ax.set_xlabel('Year')
ax.set_ylabel('Patents / Researcher')
fig.savefig('output/graphics/academic_gdp_predict.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
dat['total_rnd'].plot(ax=ax, label='Federal R&D / GDP (%)')
dat['pred_ppres'].plot(ax=ax, label='Patents / Researcher')
ax.set_xlabel('Year')
ax.legend()
fig.savefig('output/graphics/academic_gdp_datpred.pdf', bbox_inches='tight')

##
## misclassification graph
##

targ_info = read_split('policy/targeted.txt')
unif_info = read_split('policy/uniform.txt')

targ_asubs, targ_bsubs, targ_welf = targ_info[['asubs', 'bsubs', 'alpha']]
unif_subs, unif_welf = read_split('policy/uniform.txt')[['asubs', 'alpha']]
zbar = 100*(targ_asubs/targ_bsubs)

zlist = [5, 25, 50, 75, 95]
basic = pd.concat({z: read_split(f'policy/basic_z{z:02d}.txt') for z in zlist}, axis=1)
basic[0] = pd.Series({'asubs': targ_asubs, 'bsubs': targ_bsubs, 'alpha': targ_welf})
basic[zbar] = pd.Series({'asubs': targ_asubs, 'bsubs': targ_bsubs, 'alpha': targ_welf})
basic[100] = pd.Series({'asubs': unif_subs, 'bsubs': unif_subs, 'alpha': unif_welf})
basic = basic.T.sort_index().rename_axis('z').reset_index()
basic['asubst'] = basic['z']*basic['asubs']
basic['bsubst'] = basic['bsubs']

zlist0 = np.r_[0, 5, zbar, 25, 50, 75, 95, 100]
zlist1 = np.linspace(0, 100, 101)
zinterp = lambda s: pd.Series(InterpolatedUnivariateSpline(zlist0, s)(zlist1), index=zlist1)
basic1 = basic.apply(zinterp)

# policy

fig, ax = plt.subplots(figsize=(5, 4))

ax.plot(basic1.loc[zbar:].index, basic1['bsubs'].loc[zbar:], label='True Basic, $s_b=\\tilde{s}_b$', color=colors[0])
ax.hlines(targ_bsubs, 0, zbar, color=colors[0])

ax.plot([0, zbar, 100], [targ_asubs, 0, 0], color=colors[2], label='True Applied, $s_a$')

ax.plot(basic1.loc[zbar:].index, basic1['asubs'].loc[zbar:], label='Effective Applied, $\\tilde{s}_a$', color=colors[1])
ax.hlines(targ_asubs, 0, zbar, color=colors[1])

ax.set_ylim(-5, 70)
ax.vlines(zbar, *ax.get_ylim(), linestyle='--', linewidth=1)

ax.set_xlabel('Misclassification Level $z$ (%)')
ax.set_ylabel('Optimal Subsidy Rate (%)')
ax.legend(loc='upper left', bbox_to_anchor=(0.6, 1.05))

fig.savefig('output/graphics/misclass_policy.pdf', bbox_inches='tight')

# welfare

fig, ax = plt.subplots(figsize=(5, 4))
(basic1['alpha']-100).loc[zbar:].plot(ax=ax, color=colors[0])
ax.hlines(targ_welf-100, 0, zbar, color=colors[0])

ax.set_xlim(0, 100)
ax.set_ylim(0, 3.5)
ax.vlines(zbar, *ax.get_ylim(), linestyle='--', linewidth=1)

ax.set_xlabel('Misclassification Level $z$ (%)')
ax.set_ylabel('Welfare Gain (%)')

fig.savefig('output/graphics/misclass_welfare.pdf', bbox_inches='tight')

##
## secrecy plots
##

slist = [12, 25, 37, 50, 62, 75]
secret = pd.concat({s: read_split(f'policy/bloss{s:02d}_targeted.txt') for s in slist}, axis=1)
secret[0] = pd.Series({'asubs': targ_asubs, 'bsubs': targ_bsubs, 'alpha': targ_welf})
secret = secret.T.sort_index().rename_axis('s')

slist0 = np.r_[0, 12, 25, 37, 50, 62, 75]
slist1 = np.linspace(0, 75, 76)
sinterp = lambda s: pd.Series(InterpolatedUnivariateSpline(slist0, s)(slist1), index=slist1)
secret1 = secret.apply(sinterp)

secret1.index = 100 - secret1.index
secret1 = secret1.sort_index().loc[35:]

fig, ax = plt.subplots(figsize=(5, 4))
secret1['asubs'].plot(ax=ax, label='Applied')
secret1['bsubs'].plot(ax=ax, label='Basic')
ax.legend()
ax.set_xlabel('Diffusion Level $d$ (%)')
ax.set_ylabel('Subsidy Rate (%)')
fig.savefig('output/graphics/secrecy_policy.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(5, 4))
(secret1['alpha']-100).plot(ax=ax)
ax.set_ylim(0, 3)
ax.set_xlabel('Diffusion Level $d$ (%)')
ax.set_ylabel('Welfare Gain (%)')
fig.savefig('output/graphics/secrecy_welfare.pdf', bbox_inches='tight')

##
## citation distribution plots
##

cmax = 40
L = 15
basic_dir = 'targets/citations'

fname_cite = f'{basic_dir}/freq_cites_frgranted_{L}yrs.csv'
cite_name = f'TotN_citations_{L}yrs'
freq_name = f'tot_freq_{L}yr'

cite_dat = pd.read_csv(fname_cite, skipfooter=2, engine='python')
cite_dat = cite_dat.rename(columns={cite_name: 'cites', freq_name: 'count'})
cite_dat = cite_dat.set_index(['cites', 'public'])
cite_dat = cite_dat.unstack()['count'].rename(columns={0: 'private', 1: 'public'})
cite_dat = cite_dat.rename_axis('type', axis=1)
cite_dat = cite_dat.fillna(0).reindex(np.arange(cmax+1), fill_value=0).astype(np.int)
cite_dat = cite_dat/cite_dat.sum()

cite_mod = pd.read_csv('output/citations.csv', index_col='cites')
cite_mod = cite_mod.rename_axis('type', axis=1)
cite_merge = pd.concat({'data': cite_dat, 'model': cite_mod}, axis=1, names=['source'])
cite_merge = cite_merge.reorder_levels([1,0], axis=1)

fig, ax = plt.subplots(figsize=(5, 4))
cite_merge['private', 'model'].loc[:40].plot(ax=ax, color=colors[0], label='Model')
cite_merge['private', 'data'].loc[:40].plot(ax=ax, color=colors[2], label='Data')
ax.set_ylim(0, 0.16)
ax.legend()
ax.set_xlabel('Citations')
ax.set_ylabel('Probability')
fig.savefig('output/graphics/citations_private.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(5, 4))
cite_merge['public', 'model'].loc[:40].plot(ax=ax, color=colors[0], label='Model')
cite_merge['public', 'data'].loc[:40].plot(ax=ax, color=colors[2], label='Data')
ax.set_ylim(0, 0.16)
ax.legend()
ax.set_xlabel('Citations')
ax.set_ylabel('Probability')
fig.savefig('output/graphics/citations_public.pdf', bbox_inches='tight')

##
## identify spillover plots
##

slast2 = lambda s: int(s[-2])
moments = read_split(f'targets/target{mtag}.txt', delim=':')
dat_bri = moments[names_bri].rename(slast2).rename('Data')
dat_epb = moments[names_epb].rename(slast2).rename('Data')

slast1 = lambda s: int(s[-1])
simulate = pd.read_csv('output/identify_spillover.csv', index_col='pb')
mod_bri = simulate[[f'bri{i}' for i in range(1, 9)]].T.rename(mapper=slast1).rename(mapper=lambda x: f'$p_b={x}$', axis=1)
mod_epb = simulate[[f'epb{i}' for i in range(1, 9)]].T.rename(mapper=slast1).rename(mapper=lambda x: f'$p_b={x}$', axis=1)

fig, ax = plt.subplots(figsize=(5, 4))
mod_bri.plot(ax=ax, c=colors[0])
dat_bri.plot(ax=ax, c=colors[2])

mod_patch = mpl.patches.Patch(color=colors[0], label='Model')
dat_patch = mpl.patches.Patch(color=colors[2], label='Data')
plt.legend(handles=[mod_patch, dat_patch])

xmin, xmax = ax.get_xlim()
for pb in mod_bri:
    ax.annotate(pb, xy=(xmax+0.1, mod_bri[pb].iloc[-1]-0.002), annotation_clip=False)
fig.subplots_adjust(right=0.8)

ax.set_xticks(range(1, 9))
ax.set_xlabel('Number of Industries')
ax.set_ylabel('Basic Research Intensity')
fig.savefig('output/graphics/identify_spillover.pdf', bbox_inches='tight')

##
## bri and epb match plots
##

compfs = pd.read_csv('output/moments_used.csv', index_col='Description')
compfs_bri = compfs.loc[names_bri].rename(mapper=slast2).rename_axis('Number of Industries')
compfs_epb = compfs.loc[names_epb].rename(mapper=slast2).rename_axis('Number of Industries')

fig, ax = plt.subplots(figsize=(5, 4))
compfs_bri['Model'].plot(ax=ax, color=colors[0])
compfs_bri['Data'].plot(ax=ax, color=colors[2])
ax.set_ylim(0, 0.15)
ax.legend()
ax.set_ylabel('Basic Research Intensity')
fig.savefig('output/graphics/bri_match.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(5, 4))
compfs_epb['Model'].plot(ax=ax, color=colors[0])
compfs_epb['Data'].plot(ax=ax, color=colors[2])
ax.set_ylim(0, 0.8)
ax.legend()
ax.set_ylabel('Fraction Positive Basic')
fig.savefig('output/graphics/epb_match.pdf', bbox_inches='tight')

##
## applied spillover objective
##

fig, ax = plt.subplots(figsize=(5, 4))
appobj = pd.read_csv('output/objective_applied_spill.csv', index_col='pa')
appobj['obj'].plot(ax=ax)
ax.set_xlabel('Applied Spillover')
ax.set_ylabel('Estimation Objective')
fig.savefig('output/graphics/objective_applied_spill.pdf', bbox_inches='tight')
