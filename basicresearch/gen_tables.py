import numpy as np
import pandas as pd
import mectools.report as rep

##
## templates
##

latex_template = """\\documentclass[10pt]{article}

\\usepackage{amsmath}
\\usepackage{amssymb}
\\usepackage{booktabs}
\\usepackage[left=1in,right=1in,top=1in,bottom=1in]{geometry}

\\setlength{\\parindent}{0cm}
\\setlength{\\parskip}{0.2cm}

\\begin{document}

%s

\\end{document}
"""

markdown_template = """# Quantitative

%s"""

##
## useful functions
##

def read_split(fname, delim='='):
    with open(fname) as fid:
        pairs = [line.split(delim, 1) for line in fid if delim in line]
        return pd.Series({k.strip(): float(v.strip()) for k, v in pairs})

##
## begin
##

tex = ''
md = ''

##
## parameters
##

colmap_par = [
    ('disc', '1.', 'Discount Rate', '$\\delta$'),
    ('sigma', '2.', 'CRRA Utility Parameter', '$\\gamma$'),
    ('epsilon', '3.', 'Elasticity of Subsitution', '$\\varepsilon$'),
    ('pb', '4.', 'Cross-industry Spillover', '$p$'),
    ('nu', '5.', 'Multi-spillover Distribution', '$\\nu$'),
    ('eta', '6.', 'Basic Step Size', '$\\eta$'),
    ('alpha', '7.', 'Applied Step Size', '$\\lambda$'),
    ('massout', '8.', 'Mass of Entrants', '$E$'),
    ('massac', '9.', 'Mass of Academic Labs', '$U$'),
    ('kappa', '10.', 'Exogenous Exit Rate', '$\\kappa$'),
    ('agamma', '11.', 'Applied Cost Curvature', '$\\nu_a$'),
    ('bgamma', '12.', 'Basic Cost Curvature', '$\\nu_b$'),
    ('asigma', '13.', 'Applied Cost Scale', '$\\xi_a$'),
    ('bsigma', '14.', 'Basic Cost Scale', '$\\xi_b$'),
    ('dmu', '15.', 'Basic Fixed Mean', '$\\bar{\\phi}$'),
    ('dsig', '16.', 'Basic Fixed Std. Dev.', '$\\sigma$'),
    ('zeta', '17.', 'Product Cooldown Rate', '$\\zeta$'),
    ('merge', '18.', 'Buyout Rate', '$\\iota$'),
    ('x', '19.', 'Citation Rate', '$x$')
]

def param_table(fname):
    par = pd.DataFrame(colmap_par, columns=['index', '\\#', 'Description', 'Sym']).set_index('index')
    par['Value'] = read_split(fname, delim=':')
    return par

##
## moments
##

mmt_colspec = [(0, 10), (10, 21), (21, 27), (27, 58), (58, 69), (69, 80)]
mmt_usecols = ['Description', 'Model', 'Data']

def moment_table(fname):
    mmt = pd.read_fwf(fname, mmt_colspec, usecols=mmt_usecols).dropna()
    mmt['Model'] = pd.to_numeric(mmt['Model'])
    mmt['Data'] = pd.to_numeric(mmt['Data'])
    return mmt

##
## optimal policy
##

polmap = {
    'baseline': 'Baseline',
    'socplan': 'Soc Plan',
    'uniform': 'Uniform',
    'targeted': 'Targeted',
    'academic': 'Academic',
    'uniform_academic': 'Unif Acad',
    'zero': 'Zero',
    'zero_academic': 'Zero Acad',
    'academic_zsub': 'Acad + Zsub',
    'bloss25_targeted': 'Secret 25',
    'bloss50_targeted': 'Secret 50',
    'bloss75_targeted': 'Secret 75',
    'bdole_academic': 'Academic+',
    'bdole_uniform_academic': 'Unif Acad+'
}
poltypes = list(polmap)

colmap_pol = {
    'asubs': '$\\psi_a$',
    'bsubs': '$\\psi_b$',
    'academic gdp frac': '$R/Z$',
    'abar': '$\\tau_a$',
    'bbar': '$\\tau_b$',
    'dbar': '$\\tau_u$',
    'Lprod': '$L_p$',
    'Lbas': '$L_b$',
    'Lac': '$L_u$',
    'Lapp_entrant': '$L_e$',
    'Lapp_incumbent': '$L_a$',
    'psi': '$\\alpha$',
    'g': '$g$',
    'alpha': '$\\beta$'
}
columns_pol = list(colmap_pol)

def policy_table(fname):
    return read_split(fname)[columns_pol].rename(colmap_pol)

##
## default behavior
##

if __name__ == '__main__':
    import os
    import argparse

    parser = argparse.ArgumentParser(description='Generate tables for back to basics.')
    parser.add_argument('--par', type=str, default='1_2', help='parameter set to use')
    parser.add_argument('--pol', type=str, default=None, help='policy optima to use')
    args = parser.parse_args()

    # parse files
    params = param_table(f'params/params{args.par}.txt')
    moments = moment_table(f'moments/moments{args.par}.txt')
    policies = {}
    realpols = []
    for p in poltypes:
        if args.pol is not None:
            fname = f'policy/{args.pol}_{p}.txt'
        else:
            fname = f'policy/{p}.txt'
        if os.path.exists(fname):
            policies[p] = policy_table(fname)
            realpols.append(p)
    policies = pd.concat(policies, axis=1)
    policies = policies[realpols].rename(polmap, axis=1)

    # write markdown
    md += '## Parameters\n\n'
    md += rep.to_markdown(params, fmt='%.3f', index=False) + '\n\n'
    md += '## Moments\n\n'
    md += rep.to_markdown(moments, fmt='%.4f', index=False) + '\n\n'
    md += '## Policy\n\n'
    md += rep.to_markdown(policies, fmt='%.4f', index=True) + '\n\n'
    with open('latex/tables.md', 'w+') as fid:
        fid.write(markdown_template % md)

    # write latex
    tex += '\\section{Parameters}\n\n'
    tex += params.to_latex(float_format='%.3f', escape=False, index=False) + '\n\n'
    tex += '\\newpage\n\n'
    tex += '\\section{Moments}\n\n'
    tex += moments.to_latex(float_format='%.4f', escape=False, index=False) + '\n\n'
    tex += '\\newpage\n\n'
    tex += '\\section{Policy}\n\n'
    tex += policies.T.to_latex(float_format='%.2f', escape=False, index=True) + '\n\n'
    with open('latex/tables.tex', 'w+') as fid:
        fid.write(latex_template % tex)
