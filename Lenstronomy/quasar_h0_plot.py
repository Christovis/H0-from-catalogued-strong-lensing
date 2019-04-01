import matplotlib as plt

cosmology = "FLCDM"

savedir = os.path.join(samplesdir, cosmology, "%ix%i" % (nwalkers, nsamples))
samples_list = [lu.readpickle(os.path.join(savedir, "%s_samples.pkl" % lens.name)) for lens in lenses]
samples_list.append(lu.readpickle(os.path.join(savedir, "%s_samples.pkl" 
% ("+".join([l.name for l in lenses])))))


### plot nice H0s histogram 
plt.figure(figsize=(5, 5), dpi=200)
plt.subplot(1, 1, 1)

H0s_list = [[s[0] for s in samples] for samples in samples_list]
oms_list = [[s[1] for s in samples] for samples in samples_list]

percentiles = [16, 50, 84]
colors = ["crimson", "royalblue", "lightslategrey", "seagreen", "black"]
nbins = 35
title = r"${{{0}}}_{{-{1}}}^{{+{2}}}$"
fmt = "{{0:{0}}}".format(".1f").format


# loop over the lenses, shaded histograms
for il, lens in enumerate(lenses):
    h, be = np.histogram(H0s_list[il], bins=nbins, density=True)
    xs = [(b+be[ind+1])/2. for ind, b in enumerate(be[:-1])]
    plt.plot(xs, h, alpha=0.5, color=colors[il], linewidth=0.0)
    plt.fill_between(
            xs, h, alpha=0.5, color=colors[il], linewidth=0.0,
            label=r'%s' % lens.longname)

    # add the values
    pcs = np.percentile(H0s_list[il], q=percentiles)
    txt = r'$H_{0}: $' + \
            title.format(fmt(pcs[1]), fmt(pcs[1]-pcs[0]), fmt(pcs[2]-pcs[1]))
    plt.annotate(
            txt, xy=(0.0, 0.0), xytext=(0.02, 0.9-0.1*il),
            xycoords='axes fraction', fontsize=18, color=colors[il])


# plot the combined result
h, be = np.histogram(H0s_list[-1], bins=nbins, density=True)
xs = [(b+be[ind+1])/2. for ind, b in enumerate(be[:-1])]
plt.plot(xs, h, alpha=1.0, color=colors[-1], linewidth=2.0, label=r'All')

# add the values
pcs = np.percentile(H0s_list[-1], q=percentiles)
txt = r'$H_{0}: $' + title.format(fmt(pcs[1]), fmt(pcs[1] - pcs[0]), fmt(pcs[2] - pcs[1]))
plt.annotate(txt, xy=(0.0, 0.0), xytext=(0.02, 0.9 - 0.1 * len(lenses)), xycoords='axes fraction',
fontsize=18, color=colors[-1])


# fine tune
if cosmology == "FLCDM":
title = r"$H_{\rm{0}} \rm{[0, 150]} \ \ \ \Omega _{\rm{m}} \rm{[0.05, 0.5]}$"
elif cosmology == "FwCDM":
title = r"$H_{\rm{0}} \rm{[0, 150]} \ \ \ \Omega _{\rm{m}} \rm{[0.05, 0.5]} \ \ \  w \rm{[-2.5, 0.5]}$"
elif cosmology == "oLCDM":
title = r"$H_{\rm{0}} \rm{[0, 150]} \ \ \ \Omega _{\rm{m}} \rm{[0.05, 0.5]} \ \ \ \Omega _{\rm{k}} \rm{[-0.5, 0.5]}$"

plt.title(title, fontsize=18)
plt.xlabel(r"$H_{\rm{0}}\rm{\ [km\,s^{-1}\,Mpc^{-1}]}$", fontsize=24)
plt.ylabel("probability density", fontsize=18)
plt.yticks(fontsize=14)
plt.xticks(fontsize=20)
plt.yticks([])

legend = plt.legend()
legend.get_frame().set_alpha(0.0)
if cosmology in ["FLCDM", "oLCDM"]:
plt.xlim([48, 99])
plt.ylim([max(h)*0.005, max(h)*1.1])
elif cosmology == "FwCDM":
plt.xlim([48, 109])
plt.ylim([max(h)*0.005, max(h)*2.2])
plt.tight_layout()

plt.savefig(os.path.join(savedir, "H0.png"))
plt.show()
