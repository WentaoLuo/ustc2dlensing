pro gauss_smooth_shear

;This is a program to calculate the shear from galaxy sample.
;Parameters;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;-----------------------------------------------------------
c     =300000.0                ; Speed of light in unit of km/s
pi    =3.14159265              ; Pi as circumference ratio
G     =6.754e-11               ; Gravitational constant in unit of N*m^2/kg
ckm   =3.240779e-17            ; 1 km equals 3.240779e-17 kpc/h
ckg   =5.027854e-31            ; 1 kg equals 5.027854e-31 solar mass
H0    =67.3                    ; Hubble constant at z=0 frim PLANK space mission in unit of km/s/Mpc
dH0   =c/H0                    ; Hubble distance in Mpc/h
fac   =(c*c*ckg)/(4.*pi*G*ckm) ; The value of c^2/(4.*pi*G)

;Meshgrids for 2D shear;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;------------------------------------------------------------
Nbin  =20
rlow  =-4.0
rstep =8./Nbin
Xg =fltarr(Nbin*Nbin)
Yg =fltarr(Nbin*Nbin)
for i=0,Nbin-1 do begin
	for j=0,Nbin-1 do begin
		ind=i*Nbin+j
		Xg[ind]=rlow+i*rstep
		Yg[ind]=rlow+j*rstep
	endfor
endfor
;Start to read files;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;-------------------------------------------------------------
print,'% Start to load lens file... ...'
;fname ='redmapper_z_angle.dat'
fname ='lens_rot_5.dat'
;fname ='e_0.1.dat'
nl    =file_lines(fname)
data  =dblarr(5,nl)
openr,lunr,fname,/get_lun
readf,lunr,data
free_lun,lunr

ral   =data(0,*) 
decl  =data(1,*) 
zl    =data(2,*) 
dl    =data(3,*) 
angle =data(4,*)

print,'% Finish loading lens file... ...'
;-------------------------------------------------------------

print,'% Start to load the source file... ...'
ns    =file_lines('../SDSSDr7Lensing_catalog.cat')
data  =dblarr(10,ns)
openr,lunr,'../SDSSDr7Lensing_catalog.cat',/get_lun
readf,lunr,data
free_lun,lunr

ras   =data(0,*) & decs  =data(1,*)
e1    =data(2,*) & e2    =data(3,*)
res   =data(4,*) & er1   =data(5,*)
er2   =data(6,*) & zs    =data(7,*)
ds    =data(8,*) & spa   =data(9,*)

delvarx,data
print,'% Finish loading the source file... ...'
;--------------------------------------------------------------
;Start to measure shear;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print,'% Start to create the list... ...'
print,'% Time to start create the list: ',systime(0)

Nmax  =25000
lists =ptrarr(nl,/allocate_heap)

for is=0l,nl-1 do begin
	sep1          =abs((ras-ral(is))*cos(decl(is)*pi/180.))
	sep2          =abs(decs-decl(is))
	idx           =where(sep1 le 1.5 and sep2 le 1.5 and zs gt zl(is)+0.1 and res gt 1.0/3.0,ct)
	;print,format='($%"\r%d")',is
        if n_elements(idx) gt 1 then begin
		*lists(is) =idx
	endif else begin
	        *lists(is) =[0]
	endelse
endfor	

print,'% End of creating the list: ',systime(0)
print,'% Finish creating the list... ...'

;--------------------------------------------------------------

print,'% Calculating the shear... ...'
gm1  =fltarr(Nbin,Nbin,nl)
gm2  =fltarr(Nbin,Nbin,nl)
wgt  =fltarr(Nbin,Nbin,nl)
gmer =fltarr(Nbin,Nbin,nl)
wgts =fltarr(Nbin,Nbin,nl)

resp =2.0*(1.0-variance(e1,/nan))
print,'% Starting time: ',systime(0)
	
for j=0l,nl-1 do begin
                
	ipt =*lists(j)
	;rab =ras[ipt] & decb=decs[ipt]
	ra1 =ras[ipt] & dec1=decs[ipt]
	me1 =e1[ipt]  & me2 =e2[ipt]
	rs  =res[ipt] & mer1=er1[ipt]
	mer2=er2[ipt] & zs1 =zs[ipt]
	dis1=ds[ipt]  & spa1=spa[ipt]
	;agg =atan((rab-ral(j))*cos(decl(i)*pi/180.)/(decb-decl(j)))
 	;rad =sqrt(((rab-ral(j))*cos(decl(i)*pi/180.))^2+(decb-decl(j))^2)	
	;ra1 =rab*cos(agg)+decb*sin(agg)
	;dec1=-rab*sin(agg)+decb*cos(agg)

	Sig=fac*dis1*4222.0/(dl(j)*(dis1-dl(j)))/dH0/(1.0+zl(j))/(1.0+zl(j))
        wt =1.0/(0.1648+mer1*mer1+mer2*mer2)/Sig/Sig
        ;wt =1.0/(0.1648+mer1*mer1+mer2*mer2)

        xm0=cos(pi/2.0-dec1*pi/180.)
	xm1=cos(pi/2.0-decl(j)*pi/180.)
	xm2=sin(pi/2.0-dec1*pi/180.)
	xm3=sin(pi/2.0-decl(j)*pi/180.)
	;xm4=cos(((ra1-ral(j))*cos(decl(j)*pi/180.))*pi/180.)
	xm4=cos((ra1-ral(j))*pi/180.)
	the=acos(xm0*xm1+xm2*xm3*xm4)

	dang1=dl(j)*sin(cos(dec1*pi/180.)*(ra1-ral(j))*pi/180.0)*(1.+zl(j))
	dang2=dl(j)*sin((dec1-decl(j))*pi/180.0)*(1.+zl(j))
	          
	;sph=(2.0*sin((ra1-ral(j))*pi/180.0)*sin((ra1-ral(j))*pi/180.0))/sin(the)/sin(the)-1.0
	;cph=(2.0*sin((ra1-ral(j))*pi/180.0)*sin((dec1-decl(j))*pi/180.0))/sin(the)/sin(the)
	tpsin=sin(((ra1-ral(j))*cos(dec1*pi/180.))*pi/180.)
	tpcin=sin((dec1-decl(j))*pi/180.)
	sph=(2.0*tpsin*tpsin)/sin(the)/sin(the)-1.0
	cph=(2.0*tpsin*tpcin)/sin(the)/sin(the)

        for iy=0,Nbin-1 do begin
        for iw=0,Nbin-1 do begin
		ipp=iy*Nbin+iw
		rr =(dang1-Xg[ipp])*(dang1-Xg[ipp])+(dang2-Yg[ipp])*(dang2-Yg[ipp])
		gss=exp(-0.5*(rr/0.16))	
		;gss=exp(-0.5*(rr/0.04))	
		ang=spa1*pi/90.0
		ep =cos(ang)*me1-sin(ang)*me2
		em =sin(ang)*me1+cos(ang)*me2
		
		e45=cph*ep+sph*em
		et =-sph*ep+cph*em
		;e45=em
		;et =ep
		gm1(iy,iw,j)=total(et*wt*Sig*gss/rs,/nan)
                gm2(iy,iw,j)=total(e45*wt*Sig*gss/rs,/nan)
                wgt(iy,iw,j)=total(wt*gss,/nan)
		gmer(iy,iw,j)=total(et*et*wt*wt*gss*gss*Sig*Sig/rs/rs,/nan)

		;gm1(iy,iw,j) =total(et*wt*gss/rs,/nan)
		;gm2(iy,iw,j) =total(e45*wt*gss/rs,/nan)
		;wgt(iy,iw,j) =total(wt*gss,/nan)
		;gmer(iy,iw,j)=total(et*et*wt*wt*gss*gss/rs/rs,/nan)

		wgts(iy,iw,j)=total(wt*wt*gss*gss,/nan)
	endfor
	endfor

endfor
delvarx,ral,decl,zl,dl,ras,decs,zs,ds
delvarx,e1,e2,res,er1,er2,spa
delvarx,lists 
print,'% Finish calculating the shear... ...'

;--------------------------------------------------------------

;Output;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
filecov='ESDnorot_400kpc_20_covar'
openw,luncov,filecov,/get_lun
;fileden='ESDnorot_400kpc_40'
;openw,lundes,fileden,/get_lun

Nboot=500

for ixx=0,Nbin-1 do begin
for iyy=0,Nbin-1 do begin
	imm  =ixx*Nbin+iyy
	gtan =fltarr(Nboot)
	for iboot=0,Nboot-1 do begin
		ibtx       =long(randomu(seed,nl)*(nl+1))
		gtan(iboot)=total(gm1(ixx,iyy,ibtx))/total(wgt(ixx,iyy,ibtx))/resp
		printf,luncov,ixx,iyy,gtan(iboot),$
			format='(I3,x,I3,x,D16.6)'
	endfor
;	err    =sqrt(variance(gtan)) 	
;	delsens=mean(gtan)
;	printf,lundes,Xg(imm),Yg(imm),deldens,err,$
;		format='(D16.6,x,D16.6,x,D16.6,x,D16.6)'
endfor
endfor
free_lun,luncov
;free_lun,lundes

;Make a covariance matrix image---------------------------------
;print,'% Start to make covariance matrix'
;readcol,filecov,iboot,id,gm1,format=('ID'),/silent
;covar=fltarr(Nbin*Nbin,Nbin*Nbin)

;for ix=0,Nbin*Nbin-1 do begin
;    for iy=0,Nbin*Nbin-1 do begin

 ;       idx =where(id eq ix,ctx)
 ;       idy =where(id eq iy,cty)
 ;       mux =mean(gm1[idx])
 ;       muy =mean(gm1[idy])
 ;       sum=0.0
 ;
;        for j=0,Nboot-1 do begin
;               inx=idx(j)
;               iny=idy(j)
;               sum=sum+(gm1(inx)-mux)*(gm1(iny)-muy)
;        endfor

;        covar(ix,iy)=sum/float(Nboot)
        ;print,sum/500.
;    endfor
;endfor

;writefits,'covariance.fits',covar

end


print,'% Finishing time: ',systime(0)
Print,'% Everything is done!'

end
