module parametros

	integer, parameter :: ikind=selected_real_kind(p=18)
	integer, parameter :: DIM=3
	!real(kind=ikind), parameter :: CG=0.0002959867255922384281755837575197813572649020894500715428756906
	!real(kind=ikind), parameter :: CG=6988.191152154 !tierra luna
	real(kind=ikind), parameter :: CG=0.0001184462068986038 !unidades astronomicas-years

	integer, parameter :: N=3
	
	integer, parameter :: nt=365000 !itearion number
	
	real(kind=ikind), parameter :: Rc=0.03270138 !radio de los planetas
	


	real(kind=ikind) :: vp

	
	real(kind=ikind), parameter :: dt=0.0000027 !itearation step
	real(kind=ikind), dimension(3,N,2) :: r,v,a,rret,rret1 !!posicion,velocidad,aceleracion
	real(kind=ikind), dimension(N) :: m !!masa
	real(kind=ikind) :: az,min
  real(kind=ikind) :: rtx1, rtx2

	character*80 :: datin
	integer :: dat

	
	integer :: l,i,j,o,vpg,jj,stp,ir
	real(kind=ikind) :: sp,sk,st !!energia potencial, cinetica
	real(kind=ikind), dimension (N,3):: p,k
    integer :: step,nn,aux
    real(kind=ikind) :: eka,epa,eta,AC,BC,CC,DC,EC,ret1,ret2,ret3
    real(kind=ikind), dimension(DIM) :: Sij,Rij,Rji,Dij
	real(kind=ikind) :: Rsqij,phij,dphij,Rsqji,phji,dphji,rtc,Rsq
	real(kind=ikind), dimension(DIM) :: real_vel

end module parametros