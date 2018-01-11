;#############################################################################
;
; This PPXF_POPULATION_GAS_EXAMPLE_SDSS routine shows how to study stellar population with
; the procedure PPXF, which implements the Penalized Pixel-Fitting (pPXF) method by
; Cappellari M., & Emsellem E., 2004, PASP, 116, 138.
;
; This example includes gas emission lines in the fit instead of masking them.
;
; MODIFICATION HISTORY:
;   V1.0.0: Adapted from PPXF_KINEMATICS_EXAMPLE.
;       Michele Cappellari, Oxford, 12 October 2011
;   V1.1.0: Made a separate routine for the construction of the templates
;       spectral library. MC, Vicenza, 11 October 2012
;   V1.1.1: Includes regul_error definition. MC, Oxford, 15 November 2012
;   V1.1.2: The MILES SSP models are now included in the PPXF distribution
;       with permission. MC, Oxford, 9 December 2013
;   V1.1.3: Illustrate how to plot the gas and stellar spectra separately.
;       MC, Oxford, 16 November 2015
;   V1.1.4: Use updated ppxf_emission_lines routine, with one extra parameter.
;       Fixed bug in printed fluxes. MC, Oxford, 25 January 2016
;   V1.1.5: Included treatment of the SDSS/MILES vacuum/air wavelength difference.
;       MC, Oxford, 10 August 2016
;
;#############################################################################
pro setup_spectral_library, velScale, FWHM_gal, templates, lamRange_temp, logLam_temp
compile_opt idl2

; Read the list of filenames from the Single Stellar Population library
; by Vazdekis et al. (2010, MNRAS, 404, 1639) http://miles.iac.es/.
;
; For this example I downloaded from the above website a set of
; model spectra with default linear sampling of 0.9A/pix and default
; spectral resolution of FWHM=2.51A. I selected a Salpeter IMF
; (slope 1.30) and a range of population parameters:
;
;     [M/H] = [-1.71, -1.31, -0.71, -0.40, 0.00, 0.22]
;     Age = range(1.0, 17.7828, 26, /LOG)
;
; This leads to a set of 156 model spectra with the file names like
;
;     Mun1.30Zm0.40T03.9811.fits
;
; IMPORTANT: the selected models form a rectangular grid in [M/H]
; and Age: for each Age the spectra sample the same set of [M/H].
;
; We assume below that the model spectra have been placed in the
; directory "miles_models" under the current directory.
;
vazdekis = file_search('Miles_Models/Mun1.30*.fits',COUNT=nfiles)
FWHM_tem = 2.51 ; Vazdekis+10 spectra have a resolution FWHM of 2.51A.

; Extract the wavelength range and logarithmically rebin one spectrum
; to the same velocity scale of the SAURON galaxy spectrum, to determine
; the size needed for the array which will contain the template spectra.
;
fits_read, vazdekis[0], ssp, h2
lamRange_temp = sxpar(h2,'CRVAL1') + [0d,sxpar(h2,'CDELT1')*(sxpar(h2,'NAXIS1')-1d)]
log_rebin, lamRange_temp, ssp, sspNew, logLam_temp, VELSCALE=velScale

; Create a three dimensional array to store the
; two dimensional grid of model spectra
;
nAges = 26
nMetal = 6
templates = dblarr(n_elements(sspNew),nAges,nMetal)

; Convolve the whole Vazdekis library of spectral templates
; with the quadratic difference between the SAURON and the
; Vazdekis instrumental resolution. Logarithmically rebin
; and store each template as a column in the array TEMPLATES.

; Quadratic sigma difference in pixels Vazdekis --> SAURON
; The formula below is rigorously valid if the shapes of the
; instrumental spectral profiles are well approximated by Gaussians.
;
FWHM_dif = SQRT(FWHM_gal^2 - FWHM_tem^2)
sigma = FWHM_dif/2.355/sxpar(h2,'CDELT1') ; Sigma difference in pixels

; IMPORTANT: To avoid spurious velocity offsets of the templates, the
; NPIXEL keyword in PSF_GAUSSIAN must be an odd integer as done below
;
lsf = psf_Gaussian(NPIXEL=2*ceil(4*sigma)+1, ST_DEV=sigma, /NORM, NDIM=1)

; Here we make sure the spectra are sorted in both [M/H]
; and Age along the two axes of the rectangular grid of templates.
; A simple alphabetical ordering of Vazdekis's naming convention
; does not sort the files by [M/H], so we do it explicitly below
;
metal = ['m1.71', 'm1.31', 'm0.71', 'm0.40', 'p0.00', 'p0.22']
for k=0,nMetal-1 do begin
    w = where(strpos(vazdekis, metal[k]) gt 0)
    for j=0,nAges-1 do begin
        fits_read, vazdekis[w[j]], ssp
        ssp = convol(ssp,lsf) ; Degrade template to SAURON resolution
        ; From IDL 8.1 one can use the following commented line instead of the
        ; above one, and the line with PSF_GAUSSIAN is not needed any more.
        ; ssp = gauss_smooth(ssp,sigma)
        log_rebin, lamRange_temp, ssp, sspNew, VELSCALE=velScale
        templates[*,j,k] = sspNew ; Templates are *not* normalized here
    endfor
endfor

end
;------------------------------------------------------------------------------
pro ppxf_population_gas_example_sdss
compile_opt idl2

; Read SDSS DR8 galaxy spectrum taken from here http://www.sdss3.org/dr8/
; The spectrum is *already* log rebinned by the SDSS DR8
; pipeline and log_rebin should not be used in this case.
;
t = mrdfits('spectra/NGC3522_SDSS.fits',1,h1,/SILENT)
z = float(sxpar(h1,"z")) ; SDSS redshift estimate

; Only use the wavelength range in common between galaxy and stellar library.
;
t = t[where(t.wavelength gt 3540 and t.wavelength lt 7409)]
galaxy = t.flux/median(t.flux)  ; Normalize spectrum to avoid numerical issues
wave = t.wavelength

; The SDSS wavelengths are in vacuum, while the MILES ones are in air.
; For a proper treatment, the SDSS vacuum wavelengths should be
; converted into air wavelengths and the spectra should be resampled.
; To avoid resampling, given that the wavelength dependence of the
; correction is very weak, I approximate it with a constant factor.
;
vactoair, wave, waveNew
wave *= median(waveNew/wave)

; The noise level is chosen to give Chi^2/DOF=1 without regularization (REGUL=0).
; A constant error is not a bad approximation in the fitted wavelength
; range and reduces the noise in the fit.
;
noise = galaxy*0 + 0.01528

; The velocity step was already chosen by the SDSS pipeline
; and we convert it below to km/s
;
c = 299792.458d ; speed of light in km/s
velScale = c*alog(wave[1]/wave[0])
FWHM_gal = 2.76 ; SDSS has an approximate instrumental resolution FWHM of 2.76A.

;------------------- Setup templates -----------------------

setup_spectral_library, velScale, FWHM_gal, $
    stars_templates, lamRange_temp, logLam_temp

; The stellar templates are reshaped into a 2-dim array with each spectrum
; as a column, however we save the original array dimensions, which are
; needed to specify the regularization dimensions
;
s1 = size(stars_templates,/DIM)
reg_dim = s1[1:*]
nstars = product(reg_dim) ; number of SSP stellar templates
stars_templates = reform(stars_templates,s1[0],nstars)

; See the pPXF documentation for the keyword REGUL,
; for an explanation of the following two lines
;
regul_err = 0.004 ; Desired regularization error
stars_templates /= median(stars_templates) ; Normalizes stellar templates by a scalar

; Construct a set of Gaussian emission line templates.
; Estimate the wavelength fitted range in the rest frame.
;
lamRange_gal = [min(wave), max(wave)]/(1 + z)
gas_templates = ppxf_emission_lines(logLam_temp, lamRange_gal, FWHM_gal, LINE_NAMES=line_names, LINE_WAVE=line_wave)
ngas = (size(gas_templates,/DIM))[1] ; number of gas emission line templates

; Combines the stellar and gaseous templates into a single array
; during the PPXF fit they will be assigned a different kinematic
; COMPONENT value
;
templates = [[stars_templates], [gas_templates]]

;-----------------------------------------------------------

; The galaxy and the template spectra do not have the same starting wavelength.
; For this reason an extra velocity shift DV has to be applied to the template
; to fit the galaxy spectrum. We remove this artificial shift by using the
; keyword VSYST in the call to PPXF below, so that all velocities are
; measured with respect to DV. This assume the redshift is negligible.
; In the case of a high-redshift galaxy one should de-redshift its
; wavelength to the rest frame before using the line below as described
; in PPXF_KINEMATICS_EXAMPLE_SAURON.
;
c = 299792.458d
dv = c*alog(lamRange_temp[0]/wave[0]) ; km/s
vel = c*alog(1 + z) ; Relation between redshift and velocity in pPXF

; Here the actual fit starts. The best fit is plotted on the screen.
;
; IMPORTANT: Ideally one would like not to use any polynomial in the fit
; as the continuum shape contains important information on the population.
; Unfortunately this is often not feasible, due to small calibration
; uncertainties in the spectral shape. To avoid affecting the line strength of
; the spectral features, we exclude additive polynomials (DEGREE=-1) and only use
; multiplicative ones (MDEGREE=10). This is only recommended for population, not
; for kinematic extraction, where additive polynomials are always recommended.
;
start = [vel, 180d] ; (km/s), starting guess for [V,sigma]
!P.MULTI=[0,1,2]

; Assign COMPONENT=0 to the stellar templates and
; COMPONENT=1 to the gas emission lines templates.
; One can easily assign different components to different gas species
; e.g. COMPONENT=1 for the Balmer series, COMPONENT=2 for the [OIII] doublet, ...)
; Input a negative MOMENTS value to keep fixed the LOSVD of a component.
;
component = [replicate(0,nstars), replicate(1,ngas)]
moments = [4, 2] ; fit (V,sig,h3,h4) for the stars and (V,sig) for the gas
start = [[start], [start]] ; adopt the same starting value for both gas and stars

ppxf, templates, galaxy, noise, velScale, start, sol, $
    GOODPIXELS=goodPixels, MOMENTS=moments, DEGREE=-1, MDEGREE=10, $
    VSYST=dv, REGUL=1./regul_err, WEIGHTS=weights, REG_DIM=reg_dim, $
    COMPONENT=component, MATRIX=matrix, BESTFIT=bestfit

; When the two numbers below are the same, the solution is the smoothest
; consistent with the observed spectrum.
;
print, 'Desired Delta Chi^2:', sqrt(2*n_elements(goodPixels))
print, 'Current Delta Chi^2:', (sol[6] - 1)*n_elements(goodPixels)

w = where(component eq 1)  ; Extract weights of gas emissions only
print, '++++++++++++++++++++++++++++++'
print, FORMAT='("Gas V=", G0.4, " and sigma=", G0.2, " km/s")', sol[0:1, 1]  ; component=1
print, 'Emission lines fluxes:'
for j=0, ngas-1 do $
    print, FORMAT='(A12, ": ", G0.3)', line_names[j], weights[w[j]]*total(matrix[*, w[j]])
print, '------------------------------'

; Plot best fit separately for the gas and stars
;
cgplot, wave, galaxy, color='black', XTITLE='Observed Wavelength A', $
    YTITLE='Relative Flux', /XSTYLE, YRANGE=[-0.1, 1.3]
cgoplot, wave, bestfit, color='blue'
cgOPlot, wave, galaxy - bestfit, PSYM=4, COLOR='limegreen'
stars = matrix[*, 0:nstars-1] # weights[0:nstars-1]
cgoplot, wave, stars, COLOR='red'  ; overplot stellar templates alone
gas = matrix[*, nstars:*] # weights[nstars:*]
cgoplot, wave, gas + 0.15, COLOR='blue'  ; overplot emission lines alone
cgoplot, wave, wave*0, LINE=2

; Plot mass weights, of the stellar components alone, as a function of age and metallicity
;
loadct, 3
s = size(templates)
ages = cap_range(1.0, 17.7828, s[2], /LOG)
metal = cap_range(-1.9, 0.45, s[3])  ; This axis is approximate
image_plot, reform(weights[0:nstars-1], reg_dim[0], reg_dim[1]), ages, metal, $
    XTITLE='Age (Gyr)', YTITLE='[M/H]', /XLOG, $
    AXISCOLOR='dark grey', TITLE='Mass Fraction'
!P.MULTI=0

end
;------------------------------------------------------------------------------
