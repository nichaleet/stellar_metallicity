
; This function is a calibration factor (1/thpt) for DEIMOS in KTRS Settings.
;
;	Converts e- counts to e- counts
;
FUNCTION deimos_correction_tkrs, $
	lambda	; Wavelength to calibrate

	thpt = -350.02941d + 0.26843026d * lambda - 8.5101313d-05 * lambda^2 $
		+ 1.4283871d-08 * lambda^3 -1.3388369d-12 * lambda^4 $
		+ 6.6463023d-17 * lambda^5 -1.3659080d-21 * lambda^6
	
	thpt[where(lambda lt 5363.5 or lambda gt 10000)] = 0

	return,1d/thpt
END
