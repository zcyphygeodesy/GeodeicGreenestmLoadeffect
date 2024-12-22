## Fortran codes for landwater and load effect monitoring from heterogeneous variations by Green’s integral
https://www.zcyphygeodesy.com/en/h-nd-145.html
## [Algorithm purpose]
    Using various heterogeneous geodetic variation time series as the observations and the load Green's integral as the geodynamic constraints, estimate the regional surface load equivalent water height (EWH) and all-element load effect grid time series (usually employed to represent regional time-varying gravity field).
    It is technically required that the long wave parts of the load effects on geodetic variations should be removed to satisfy the local Green’s integral condition.
    The geodetic variations here can be one or more of the following five types of variations. (1) Height anomaly variations (mm) from GNSS-leveling monitoring network, (2) disturbance gravity variations (μGal) from GNSS-gravity monitoring network or CORS-gravity tide stations, (3) ground gravity variations (μGal) from gravity monitoring network or gravity tide stations, (4) ellipsoidal height variations (mm) from CORS network or GNSS monitoring network, and (5) normal or orthometric height variations (mm) from leveling monitoring network.
## [Geophysical models]
    The Green function file LoadGreen.txt of the load indirect effect on all-element geodetic variations.
## [Main program for test entrance]
    GeodeicGreenestmLoadeffect.f90
    Input parameters: obstmsqdfl - the heterogeneous geodetic variation record time series file name.
    The file header contains the time series length and the sampling epoch time arranged with time. Record format: ID (the site name / no), 
longitude, latitude, …, weight, variation type, …, variations arranged in time series length (default value is 9999.0000).
Variation type = 1 represents the height anomaly variation (mm), = 2 represents gravity disturbance variation (μGal), = 3 represents ground gravity variation (μGal), = 4 represents ground ellipsoidal height variation (mm), and = 5 represents normal or orthometric height variation (mm).
    Input parameters: dtmfl - The calculation surface height grid file name.
    The calculation surface height is the height of the calculation point relative to the ground surface. When calculating the ground load deformation field, enter the zero-value grid. The calculation surface height grid specification is employed to specify the latitude and longitude range and spatial resolution of the land water EWH grid to be estimated.
    The program requires that the grid range of the calculation surface height must be larger than the geodetic site distribution range to absorb the edge effect. The actual effective range of the land water EWH and its load deformation field grid to be estimated will be less than the coverage range of these geodetic sites.
    Input parameters: lnth,dr - the mean distance (m) between geodetic sites and Green's integral radius (m).
    Input parameters: lvb,itern - the Laplace operator weight and cumulative approach times.
    Input parameter: nfar - the edge effect suppression parameter.
    Input parameter: knd - the method of the solution of normal equation, knd=1 LU triangular decomposition method, =2 Cholesky decomposition, =3 least square qr decomposition, =4 Minimum norm singular value decomposition.
    Input parameter: hepch - the column ordinal number of the first epoch time in header.
    Input parameter: frow - the column ordinal number of the first variation in record.
    Input parameter: kndrow - the column ordinal number of the variation type in record.
    Input parameter: wghrow - the column ordinal number of the variation weight in record.
## (1) Module of landwater and load effect inversion with Green’s integral constraints
    Loadeffectestmgreen(dtmfl,obstmsqdfl,GF,dtrow,inp)
    Input parameters: GF(8000,9) - The Green functions of the load indirect effects on height anomaly (e-13), ground gravity (e-17), gravity disturbance (e-18), ground tilt (e-14), vertical deflection (e-19), horizontal displacement (e-12), ground radial displacement (e-11), radial gravity gradient (e-15) and horizontal gravity gradient (e-15). Where the integral distance of GF(i,1:9) is equal to 100i (m).
Input parameters: dtrow - the column ordinal number of the current variations in the geodetic variation record time series file record.
Output the land water EWH grid file ewh****.dat, residual geodetic variation file rnt***.txt and 10 kinds of load effect grid files in the following.
      Greengeoid***.dat is the load effect grid file on geoid or height anomaly (mm).
      Greenterrgrav***.dat is the load effect grid file on ground gravity (μGal).
      Greengravdist***.dat is the load effect grid file on gravity disturbance (μGal).
      Greengrndtilt***.dat is the load effect vector grid file on ground tilt (SW, to the south and to the west, mas).
      Greenvertdefl***.dat is the load effect vector grid file on vertical deflection (SW, to the south and to the west, mas).
      Greenhorzdisp***.dat is the load effect vector grid file on horizontal displacement (EN, to the east and to the north, mm).
      Greenelliphgt***.dat is the load effect grid file on ground radial displacement (mm).
      Greenorthohgt***.dat is the load effect grid file on ground normal or orthometric height (mm).
      Greengradient***.dat is the load effect grid file on radial gravity gradient (mE).
      Greenhorzgrad***.dat is the load effect vector grid file on horizontal gravity gradient (NW, to the north and to the west, mE).
    Here, *** is the sampling epoch time which is also saved as the last column attribute of the load effect grid file header.
## (2) Construction module for observation equation of variation by Green 's integral constraint
    CoefLoadobs(BLH,hd,BB,nx,GF,dr,GRS,kobs,sk,kk)
    Input parameters: BLH(3) - longitude (decimal degrees), latitude (decimal degrees), height (m) of the calculation point relative to the Earth’s surface. 
    Input parameters: dr, hd(6) - the Green’s integral radius (m) and grid specification parameters (minimum and maximum longitude, minimum and maximum latitude, longitude and latitude intervals of a cell grid) of EWH to be estimated.
    Input parameters: GRS(6) - gm, ae, j2, omega, 1/f, default value
    Input parameters: GF(8000,9).
    Input parameter: kobs - the type (1~5) of the observation variation.
    Return parameters: BB(nx) - the coefficient vector of observation equation, and nx = nlat×nlon.
    Return parameters: sk(kk) - kk (<nx) is the number of nonzeros in BB(nx) and sk(kk) is the ordinal number vector of nonzero coefficients in BB(nx).
## (3) Computation module for residual surface load effects by Green's Integral
    rntGreenintegral(BLH,ewh,hd,nlat,nlon,GF,direct,indrct,GRS,dr)
    Input parameters: ewh(nlat,nlon) - the regional residual equivalent water height (EWH) variation grid (cm).
    Return parameters: direct(10) - the residual load direct effect on the height anomaly (mm), ground gravity (μGal), gravity disturbance (μGal), ground tilt (SW, to the south and to the west, mas), vertical deflection (SW, to the south and to the west, mas), radial gravity gradient (1mE) and horizontal gravity gradient (NW, to the north and to the west, mE).
    Return parameters: indrct(14)  - the residual load indirect effect on on the height anomaly (mm), ground gravity (μGal), gravity disturbance (μGal), ground tilt (SW, to the south and to the west, mas), vertical deflection (SW, to the south and to the west, mas), horizontal displacement (EN, to the east and to the north, mm), ground radial displacement (mm), ground normal or orthometric height (mm), radial gravity gradient (1mE) and horizontal gravity gradient (NW, to the north and to the west, mE). Where Indrct(11) = default value for the normal or orthometric height variation.
## (4) Module for spatial Laplace filter constraints supplemented on unknown parameters
    Laplacesmth(BPB,s,nlat,nlon,lvb)
    Input parameter: lvb - the Laplace operator weight.
    Input parameters: s(nlat*nlon) - the double real work vector of length nlat*nlon.
    Input and output parameters: BPB(nlat*nlon,nlat*nlon) - return the normal equation coefficient matrix after spatial Laplace filter constraints supplemented.
## (5) Module for far-zone zero constraint supplemented on unknown parameters at grid edge
    Farzonexx(BPB,nlat,nlon,nfar,wgh)
    Input parameter: nfar - the edge effect suppression parameter.
    Input parameter: wgh - the zero constraint parameter. Let BPB(k,k) = BPB(k,k) +wgh when the kth unknown parameter at the edge of the grid.
    Input and output parameters: BPB(nlat*nlon,nlat*nlon) - return the normal equation coefficient matrix after zero constraint supplemented on unknown parameters at grid edge.
## (6) Reading module of the Green function of load indirect effect
    LGrnFunc(loadgrfl,GF)
    Input parameters: loadgrfl - The Green function file name LoadGreen.txt of the load indirect effect.
    Return parameters: GF(8000,9).
## (7) Calculation module for the normal gravity field
    GNormalfd(BLH,NFD,GRS)
    Return parameters: NFD(5) - the normal geopotential (m2/s2), normal gravity (mGal), normal gravity gradient (E), normal gravity line direction (', expressed by its north declination relative to the center of the Earth center of mass) or normal gravity gradient direction (', expressed by its north declination relative to the Earth center of mass).
## (8) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt2(pn,dp1,dp2,n,t) ! t=cos ψ
## (9) Large normal equation solution module package
    EqHestens(BPB,xx,nn,BPL); EqJordon(BPB,xx,nn,BPL)
    EqCholesky(BPB,XX,nn,BPL); EqueSVD(BPB,XX,nn,BPL)
    RidgeEstimate(BPB,xx,nn,BPL); Equsolve(BPB,xx,nn,BPL,knd,bf) 
## (10) Other auxiliary modules
    BLH_RLAT(GRS, BLH, RLAT); PickReclong(line, kln, rec, nn); Stat1d(dt,nn,rst)
    IntpGrnF(GF,dl,vfn); Gauss2D(lon,lat,dt,row,col,hd)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler. mkl_lapack95_ilp64.lib link library required.
DOS executable file and all input and output data.
