"""
Piece of code used earlier to generate SN 1006 tables
embedded within models_all_exec.py

Not usable as is, just for archiving
Aaron Tran
2014 August 3
"""

def generate_SN1006_table(mu):
    """What it says.  Configured for single mu value.
    Move to a configuration script if any more SN1006 tables are generated."""

    #mu_vals = [0, 1./3, 1./2, 1, 1.5, 2]  # Following Sean
    #mu_vals = np.array([1./3, 1./2, 1])  # Kolmogorov, Kraichnan, Bohm
    mu_vals = [mu]
    eta2_vals = np.logspace(-2, 2, 50, base=10)
    eta2_vals = np.sort(np.append(eta2_vals, np.linspace(0, 10, 50)))
    n_B0 = 20  # In practice, you'll usually get ~1.5 to 2x as many points
               # as code tries to achieve good spacing

    snr = snrcat.make_SN1006()
    kevs = SN1006_KEVS
    data = np.array([SN1006_DATA[flmt][0] for flmt in [1,2,3,4,5]])
    data_min = np.amin(data/1.51, axis=0)
    data_max = np.amax(data*1.51, axis=0)

    fname = 'sn1006_grid-1-100-20_2014-07-28_mu-{:0.2f}_partest.pkl'.format(mu)

    tab = models.maketab(snr, kevs, data_min, data_max, mu_vals, eta2_vals,
                         n_B0, fname=fname)

    return tab
