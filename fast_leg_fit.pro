FUNCTION fast_leg_fit, x, y, DEGREE=degree, YFIT=y_fit
;+
; NAME:
;
; fast_leg_fit
;
; PURPOSE:
;
; Do a simple least squares fit to get the best-fitting legendre polynomials,
; up to degree 10, describing y as a function of x.
;
; CATEGORY:
;
; fitting
;
; CALLING SEQUENCE:
;
; coefficients = fast_leg_fit(x, y, degree=degree, y_fit = y_fit)
;
; INPUTS:
;
; x - vector of independent data
; y - vector of dependent data
;
; OPTIONAL INPUTS:
;
; degree - degree of the fit, default is 2
;
; KEYWORD PARAMETERS:
;
; none
;
; OUTPUTS:
;
; coefficients of the fit
;
; OPTIONAL OUTPUTS:
;
; y_fit - vector of fitted y values
;
; COMMON BLOCKS:
;
; none
;
; SIDE EFFECTS:
;
; none
;
; RESTRICTIONS:
;
; none
;
; PROCEDURE:
;
; Simple matrix inversion. If you are at all concerned about robustness, use
; SVDFIT, /LEGENDRE instead. This approach is *much* faster, though.
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
; 30Apr09 - documented, leroy@mpia.de
;
;-


    COMPILE_OPT idl2

    if n_elements(degree) eq 0 then degree=2

    nx = n_elements(x)

;   B ## coefs = y

    b = dblarr(degree+1,nx)

;   HARD-CODED LEGENDRE POLYNOMIALS UP TO DEGREE 10
    b[0,*] = 1.
    if (degree ge 1) then b[1,*] = x
    if (degree ge 2) then b[2,*] = 1./2.*(3.0*x^2-1)
    if (degree ge 3) then b[3,*] = 1./2.*(5.0*x^3-3.0*x)
    if (degree ge 4) then b[4,*] = 1./8.*(35.*x^4 - 30.*x^2 + 3.)
    if (degree ge 5) then b[5,*] = 1./8.*(63.*x^5 - 70.*x^3 + 15.*x)

    if (degree ge 6) then b[6,*] = $
       1./16.*(231.*x^6 - 315.*x^4 + 105.*x^2 - 5)
    if (degree ge 7) then b[7,*] = $
       1./16.*(429.*x^7 - 693.*x^5 + 315.*x^3 - 35.*x)
    if (degree ge 8) then b[8,*] = $
       1./128.*(6435*x^8 - 12012.*x^6 + 6930.*x^4 - 1260.*x^2 + 35.)
    if (degree ge 9) then b[9,*] = $
       1./128.*(12155.*x^9 - 25740.*x^7 + 18018.*x^5 - 4620.*x^3 + 315.*x)
    if (degree ge 10) then b[10,*] = $
       1./256.*(46189.*x^10 - 109395.*x^8 + 90090.*x^6 - 30030.*x^4. $
                + 3465.*x^4 - 63.)

    psi = invert(transpose(b) ## b)

    coefs = psi ## transpose(b) ## transpose(y)

    y_fit = reform(b ## coefs)

    return, coefs

END
