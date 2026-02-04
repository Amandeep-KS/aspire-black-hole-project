using Pyehtim
using PythonCall
plt = pyimport("matplotlib.pyplot")
builtins = pyimport("builtins")
n0 = 200 # replace it with the resolution being studied
im = Pyehtim.ehtim.image.load_txt("/home/amandeep-kaur/aspire-black-hole-project/Data_files/Picasso_2_ehtim/img_ehtim$(Int(n0)).dat")
eht = Pyehtim.ehtim.array.load_txt("/home/amandeep-kaur/aspire-black-hole-project/Data_files/EHT2017_m87.txt")
im.display()


target_flux = 0.61
current_flux = im.total_flux()
scale_factor = target_flux / current_flux
im.imvec *= scale_factor
zbl = im.total_flux() # total flux
println(zbl)


# Simulating observations using self defined parameters using im.observe

# simulation parameters
tint_sec = 30  # Integration time in seconds, 
tadv_sec = 300 # Advance time between scans
tstart_hr = 0  # GMST time of the start of the observation
tstop_hr = 24  # GMST time of the end of the observation
bw_hz = 2.e9   # Bandwidth in Hz

# generate the observation
obs = im.observe(eht, ttype="fast", tint_sec, tadv_sec, tstart_hr, tstop_hr, bw_hz, add_th_noise=true,
    sgrscat=false, ampcal=true, phasecal=true)

# Plots
# uv coverage
obs.plotall("u", "v", conj=true)
plt.savefig("uv coverage")

# Plot amplitude with baseline u-v distance.
obs.plotall("uvdist", "amp")
plt.savefig("amp vs baseline")

# Plot phase with baseline distance.
obs.plotall("uvdist", "phase")
plt.savefig("phase vs baseline")

# Simulate an observation based on real M87 data
# loading M87 observation
obs_M87_orig = Pyehtim.ehtim.obsdata.load_uvfits("/home/amandeep-kaur/Downloads/eht_imaging_tutorial-main/SR1_M87_2017_095_lo_hops_netcal_StokesI.uvfits")

# We have to change the image coordinates (ra, dec) and frequency (rf)
# to match the observed coordinates and frequency.
im_M87 = im.copy()
im_M87.ra = obs_M87_orig.ra
im_M87.dec = obs_M87_orig.dec
im_M87.rf = obs_M87_orig.rf
# Simulate an observation with the samee settings as those of the 
# M87 observation( assuming calibrated amplitudes and phases)
obs_M87_calib = im_M87.observe_same(
    obs_M87_orig,
    add_th_noise=true,
    phasecal=true,
    ampcal=true,
    ttype="fast")
# uv coverage for this simulation
obs_M87_calib.plotall("u", "v", conj=true)
plt.savefig("uv_calib coverage")

# Simulating more realistic M87-based observation


# Add scan info to Obsdata object so that per-scan stabilization and averaging
# can be done.
obs_M87_orig.add_scans()


# Std. dev. of the constant absolute gain of each telescope from a gain of 1:
gain_offset = builtins.dict([
    ("AA", 0.1),
    ("AP", 0.1),
    ("AZ", 0.1),
    ("LM", 0.6),
    ("PV", 0.1),
    ("SM", 0.1),
    ("JC", 0.1),
    ("SP", 0.1),
    ("SR", 0.0)
])
# Std. dev. of the time-varying gain differences:
gainp = builtins.dict([
    ("AA", 0.05),
    ("AP", 0.05),
    ("AZ", 0.05),
    ("LM", 0.5),
    ("PV", 0.05),
    ("SM", 0.05),
    ("JC", 0.05),
    ("SP", 0.15),
    ("SR", 0.0)
])


# Simulate an observation with phase and amplitude (gain) errors.
obs_M87 = im_M87.observe_same(
    obs_M87_orig,
    add_th_noise=true,
    phasecal=false,
    ampcal=false,
    stabilize_scan_phase=true,  # `True` makes phase error constant per scan
    stabilize_scan_amp=true,    # `True` makes gain error constant per scan
    gain_offset=gain_offset,    # use station-dependent median gain error
    gainp=gainp,                # use station-dependent variability
    ttype="fast")

# Compare phases with vs. without atmospheric noise.
Pyehtim.ehtim.plotting.comp_plots.plotall_obs_compare(
    [obs_M87_calib, obs_M87], "uvdist", "phase", ebar=false,
    legendlabels=["calibrated", "realistic"]);
plt.savefig("phase comparisons")

# Imaging


# Nominal resolution
res = obs_M87.res()
a = Pyehtim.ehtim.RADPERUAS
resolution = res / a
println("nominal resolution= $resolution muas")
npix = n0
fov = n0 * Pyehtim.ehtim.RADPERUAS
# Show the true image blurred to the nominal resolution.
blurry_im = im_M87.blur_circ(res)
blurry_im.display()
# dirty image
dim = obs_M87.dirtyimage(npix, fov)
dim.display()
# dirty beam
dbeam = obs_M87.dirtybeam(npix, fov)
dbeam.display()
# clean beam
cbeam = obs_M87.cleanbeam(npix, fov)
cbeam.display()

# Scan average observation data
# Make sure Obsdata objects have scan info.
obs_M87_calib.add_scans()

# Reduce each scan down to one average measurement.
obs_M87_calib_scan_avg = obs_M87_calib.avg_coherent(0, scan_avg=true)
obs_M87_scan_avg = obs_M87.avg_coherent(0, scan_avg=true)

# prior image
gauss_fwhm = 60 * Pyehtim.ehtim.RADPERUAS     # size of gaussian can be modified by replacing 60 
emptyim = Pyehtim.ehtim.image.make_square(obs_M87, npix, fov)
gaussim = emptyim.add_gauss(zbl, (gauss_fwhm, gauss_fwhm, 0, 0, 0))
# Make the initial image have the correct coordinates and reference frequency.
gaussim.source = "M87"
gaussim.ra = obs_M87.ra
gaussim.dec = obs_M87.dec
gaussim.rf = obs_M87.rf
gaussim.display()


# Imaging from calibrated observation

# We can use just the complex visibilities because they are well-calibrated.
data_term = data_term = Dict("vis" => 1)
# We don't need strong regularization because the complex visibilities are
# already well constraining.
reg_term = Dict("tv2" => 0.1)

# Set up imager.
imgr_calibrated = imgr = Pyehtim.ehtim.imager.Imager(
    obs_M87_calib_scan_avg, gaussim,
    prior_im=gaussim,
    flux=zbl,
    data_term=data_term,
    reg_term=reg_term,
    norm_reg=true,
    maxit=300,
    ttype="fast")

# Running imager for multiple times  with blurring to assure good convergence
function converge(imgr_calibrated; major=major, blur_frac=1.0)

    # First imaging pass
    imgr.make_image_I(show_updates=true)

    for repeat in 1:major
        init = imgr_calibrated.out_last().blur_circ(blur_frac * res)
        imgr.init_next = init
        imgr.make_image_I(show_updates=true)
    end

    return imgr_calibrated
end


# Show the results of the imaging round.
out_self_calibrated = imgr_calibrated.out_last()
plt.savefig("calibrated final plots")



# Imaging from realistic M87-based observation
# Inflate amplitude error bars with systematic noise, which accounts for
# gain errors, polarization leakage, etc.

# The following systematic noise estimates are taken from the fiducial
# eht-imaging pipeline.
systematic_noise = builtins.dict([
    ("AA", 0.01220),
    ("AP", 0.01339),
    ("AZ", 0.00860),
    ("LM", 0.17613),  # LMT had the most variability
    ("PV", 0.01220),
    ("SM", 0.01810),
    ("JC", 0.01693),
    ("SP", 0.00860)])

# The following values are taken from the fiducial eht-imaging pipeline.
data_term = Dict("amp" => 0.2, "cphase" => 1, "logcamp" => 1) # data term weights

reg_term = Dict("tv2" => 1, "simple" => 100, "tv" => 1, "flux" => 1e4)   # regularizer term weights
# To avoid gradient singularities associated in the first step:
gaussim = gaussim.add_gauss(
    zbl * 1e-3, (gauss_fwhm, gauss_fwhm, 0, gauss_fwhm, gauss_fwhm))
# set up the imager
imgr_realistic = imgr = Pyehtim.ehtim.imager.Imager(obs_M87_scan_avg, gaussim, prior_im=gaussim, flux=zbl,
    data_term=data_term, reg_term=reg_term,
    norm_reg=true, # this is very important!
    systematic_noise=systematic_noise,
    maxit=250, ttype="fast")

function converge_realistic(imgr_realistic; major=major, blur_frac=1.0)

    # First imaging pass
    imgr.make_image_I(show_updates=true)

    for repeat in 1:major
        init = imgr.out_last().blur_circ(blur_frac * res)
        imgr.init_next = init
        imgr.make_image_I(show_updates=true)
    end

    return imgr_realistic
end


# Show the results of the imaging rounds.
out_realistic = imgr_realistic.out_last()
plt.savefig("realistic_10_iter final plots")

#use the below command line to save reconstructed data file for further analysis
# out_realistic.save_txt("/home/amandeep-kaur/aspire-black-hole-project/Data_files/ehtim_reconstructed_data/ehtim_$(Int(n0)).txt")
