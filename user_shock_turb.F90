!
! User module
!
! This module contains functions that may be altered by the user to define a specific problem. 
! 
! These functions are called by tristanmainloop.F90 and should provide
! initialization of parameters specific to the problem being setup, 
! loading of initial plasma, injection of new plasma, and boundary conditions
! on particles and fields that are specific to user's problem.
! 
! Functions by name:
! read_input_user() -- parse the problem-specific section of the input file for user parameters
! get_external_fields() -- called by mover to provide externally imposed EM fields if external_fields is on in input file
! init_EMfields_user() -- initializes EM fields 
! init_particle_distribution_user() -- loads initial particle distribution
! inject_particles_user() -- injects new particles on every time step as needed
! field_bc_user() -- applies user-specific boundary conditions to fields (beyond periodic or radiation BCs)
! particle_bc_user() -- enforces user-specific BCs on particles (e.g., reflecting walls)
! shift_domain_user() -- needed for some shock problems to remove empty space; ignore. 
!
#ifdef twoD 

module m_user

	use m_globaldata
	use m_system
	use m_aux
	use m_communications
	use m_fields
	use m_particles
	use m_inputparser
	use m_fparser
	use m_domain
	
#else

module m_user_3d

	use m_globaldata_3d
	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_fields_3d
	use m_particles_3d
	use m_inputparser_3d 
	use m_fparser_3d
	use m_domain_3d

#endif

	implicit none
		
	private

!-------------------------------------------------------------------------------
!	PARAMETERS
!-------------------------------------------------------------------------------

	
!-------------------------------------------------------------------------------
!	TYPE DEFINITIONS
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!	VARIABLES
!-------------------------------------------------------------------------------

	real(sprec) :: temperature_ratio, sigma_ext, bz_ext0

	integer :: minTurbulentLambdaX, maxTurbulentLambdaX, minTurbulentLambdaY, &
			maxTurbulentLambdaY, minTurbulentLambdaZ, maxTurbulentLambdaZ
	real :: turbulenceFieldCorrection, slabFieldCorrection
	integer :: turbulenceSeed


!-------------------------------------------------------------------------------
!	INTERFACE DECLARATIONS
!-------------------------------------------------------------------------------
	
!-------------------------------------------------------------------------------
!	PUBLIC MODIFIERS
!-------------------------------------------------------------------------------

	public :: init_EMfields_user, init_particle_distribution_user, &
	inject_particles_user, read_input_user, field_bc_user, get_external_fields, &
	particle_bc_user, shift_domain_user
	public :: init_turbulent_field

	public :: evaluate_turbulent_b
	public :: evaluate_turbulent_b_slab, evaluate_turbulent_b_2d


	public :: minTurbulentLambdaX, maxTurbulentLambdaX, minTurbulentLambdaY, maxTurbulentLambdaY,&
			minTurbulentLambdaZ, maxTurbulentLambdaZ, turbulenceFieldCorrection, slabFieldCorrection



	!public :: evaluate_turbulent_b

!-------------------------------------------------------------------------------
!	MODULE PROCEDURES AND FUNCTIONS
!-------------------------------------------------------------------------------

!
! External variables read from the input file earlier; they are available for use here
!
! sigma (electron sigma),
! Binit (B field based on sigma); can be reset here
! me, mi, qe, qi (mass and charge of the species, abs(qe)/me = 1
! ppc0 (fiducial particles per cell for ions and electrons together)
! btheta, bhi -- inclination angles of initial B field, in degrees
! delgam -- k T/ m_e c^2 -- fiducial normalized temperature of electrons 
!
! gamma0, and its beta -- 4-velocity of the upstream flow. 
! mx0,my0,mz0 -- real number of cells in each direction, including the ghost zones. Active domain is from 3 to mx0-3, inclusive, and same for y and z. This is global size of the grid.
! each core knows also about mx, my, mz, which is the size of the local grid on this core. This size includes 5 ghost zones. 
!
! iglob, jglob, kglob -- functions that return global i,j,k based on local i,j,k
! xglob, yglob, zglob -- functions that return global x,y,z based on local x,y,z

	contains

!-------------------------------------------------------------------------------
! 						subroutine read_input_user		
!									
! Reads any variables related to (or needed by) this module from section "problem" in input file
! 							
!-------------------------------------------------------------------------------

subroutine read_input_user()

	implicit none
	integer :: lextflds, luserpartbcs, lwall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CHANGE THE NAME ON THIS LINE IF YOU ARE CREATING A NEW USER FILE
!This helps to identify which user file is being used. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(rank.eq.0) print *, "Using user file user_shock_turb.F90"

!inputpar_getd_def -- read real parameter; input_geti_def -- read integer parameter

	call inputpar_getd_def("problem", "temperature_ratio", 1._sprec, Temperature_ratio)
	print *, "inputpar_getd_def"
	call inputpar_getd_def("problem","sigma_ext",0._sprec,sigma_ext)
	print *, "inputpar_getd_def"
	call inputpar_geti_def("problem","external_fields",0,lextflds)
	print *, "inputpar_geti_def"
	if(lextflds .eq. 1) then 
	   external_fields =.true.
	else
	   external_fields =.false.
	endif

	if(external_fields) bz_ext0 = sqrt((gamma0)*.5*ppc0*c**2*me*(1+me/mi)*sigma_ext)

	call inputpar_geti_def("problem","user_part_bcs",0,luserpartbcs)
	print *, "inputpar_geti_def"
	if(luserpartbcs .eq. 1) then 
	   user_part_bcs = .true.
	else
	   user_part_bcs = .false.
	endif

	call inputpar_geti_def("problem", "wall", 0, lwall)
	print *, "inputpar_geti_def"
	if (lwall==1) then
		wall=.true.
	else
		wall=.false.
	endif
	
	if(wall) user_part_bcs=.true.

!read_input_user is called last, so any parameter can be overwritten here
	print *, "finish read_input_user"
end subroutine read_input_user

!-------------------------------------------------------------
!     Compute external fields to be added to the mover. 
!     These fields do not evolve via Maxwell Eqs, but can depend on time
!-------------------------------------------------------------
	subroutine get_external_fields(x,y,z,ex_ext, ey_ext, ez_ext, bx_ext,by_ext,bz_ext, qm)
	
	real,intent(inout):: bx_ext,by_ext,bz_ext, ex_ext, ey_ext, ez_ext
	real, intent(in):: x,y,z
	real, optional:: qm
	ex_ext=0.
	ey_ext=0.
	ez_ext=0.
	bx_ext=0.
	by_ext=0.
	bz_ext=0.
        ! can use bz_ext0 as fiducial field.
        !external field can be directed in any way; only invoked for external_fields = 1
        !x,y,z come in as local coordinates, convert to global by calling xglob(x), yglob(y), zglob(z)	
	end subroutine get_external_fields

!-------------------------------------------------------------------------------
! 						subroutine init_EMfields_user		 
!												
! Sets the electromagnetic fields of any specific user purpose
!							
!-------------------------------------------------------------------------------

subroutine init_EMfields_user()
	
	! local variables
	
	integer :: i, j, k
	real B0x, B0y, B0z, E0x, E0y, E0z
	real kw, v

! sqrt(sigma)=omega_c/omega_p = (qe B / gamma0 me c)/sqrt( qe 0.5 ppc0 (1+me/mi)/gamma0)  !4 pi=1, qe/me = 1 
! B = sqrt(gamma0*.5*ppc0*(1+me/mi)*c**2*(me)*sigma)
! (1+me/mi) term comes from omega_p^2 = omega_pe^2+omega_pi^2

       	Binit=sqrt(gamma0*ppc0*.5*c**2*(me*(1+me/mi))*sigma)

	!initialize B field to be set by Binit and the inclination angle 	
        !angles btheta and bphi are read earlier
	print *, 'periodic x = ', periodicx
	btheta=btheta/180.*pi
	bphi=bphi/180.*pi

	B0x=Binit*cos(btheta)
	B0y=Binit*sin(btheta)*sin(bphi)
	B0z=Binit*sin(btheta)*cos(bphi)

	E0x=0.
	E0y=(-beta)*B0z
	E0z=-(-beta)*B0y


	kw = 2*3.1415927/50;

	print *, 'init fields'

	do  k=1,mz
		do  j=1,my
			do  i=1,mx
! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)

				bx(i,j,k)=B0x;
				by(i,j,k)=B0y;
				bz(i,j,k)=B0z;

				ex(i,j,k)=E0x;
				ey(i,j,k)=E0y;
				ez(i,j,k)=E0z;
			enddo
		enddo
	enddo
#ifdef turbulence
	call init_turbulent_field()
#endif

end subroutine init_EMfields_user


!-------------------------------------------------------------------------------
! 						subroutine init_particle_distribution_user()	
!											
! Sets the particle distrubtion for a user defined case
!
!-------------------------------------------------------------------------------

subroutine init_particle_distribution_user()

	implicit none

	! local variables
	
	integer :: i, n, direction
	real    gamma_drift, delgam_i, delgam_e, ppc, weight

	real(sprec) :: x1,x2,y1,y2,z1,z2

	pcosthmult=0 !if 0 the Maxwellian distribution corresponding to 
	             ! temperature is initialized in 2D (z temperature = 0), 1 for 3D distribution (z temp nonzero). 
                     ! 1 is default

	!set initial injection points 
 
        leftwall=20.
	xinject=leftwall
	xinject2=(mx0-2.) 
        walloc = leftwall
	! ------------------------------------------

!initialize left going upstream

	gamma_drift=-gamma0 ! negative gamma_drift will send the plasma in the negative direction
	delgam_i=delgam  !delgam is read from input file in particles
	delgam_e=delgam*mi/me*Temperature_ratio

	x1=xinject  
	x2=xinject2 

	y1=3. !in global coordinates
	y2=my0-2.  
	z1=3.
	z2=mz0-2. !if 2D, it will reset automatically
	ppc=ppc0
	weight=1.

	direction=1  !drift along x, + or - determined by the sign of gamma_drift
 
	call inject_plasma_region(x1,x2,y1,y2,z1,z2,ppc,&
             gamma_drift,delgam_i,delgam_e,weight,direction)


	call check_overflow()
	call reorder_particles()
	
end subroutine init_particle_distribution_user


!-------------------------------------------------------------------------------
! 				subroutine inject_particles_user()					 
!										
! Injects new particles in the simulation on every step. To skip injection, set ppc=0 below
!
!-------------------------------------------------------------------------------

subroutine inject_particles_user()

	implicit none
	real(sprec) :: x1,x2,y1,y2,z1,z2
	real delgam_i, delgam_e, injector_speed, ppc, gamma_drift, weight

	injectedions=0 !set counter to 0 here, so that all possible other injections on this step add up
	injectedlecs=0

        x1=mx0-2.
        x2=x1
	y1=3. !in global coordinates
	y2=my0-2.  
	z1=3.
	z2=mz0-2. !if 2D, it will reset automatically
        weight=1
        injector_speed=0. !injector is not receding
        ppc=ppc0

	gamma_drift=-gamma0
	delgam_i=delgam
	delgam_e=delgam*mi/me
		
	call inject_from_wall(x1,x2,y1,y2,z1,z2,ppc,gamma_drift,&
	  delgam_i,delgam_e,injector_speed,weight)
!inject_from_wall sends a stream of e and ions along normal to the wall 
!which has x1=x2, or y1=y2, etc. This is how it finds which direction to send the stream
	
end subroutine inject_particles_user



!-------------------------------------------------------------------------------
! 				subroutine field_bc_user()	
!										
! Applies boundary conditions specific to user problem. 
! 
!-------------------------------------------------------------------------------

subroutine field_bc_user(time)
	implicit none

		real, intent(in) :: time
        real xmin, xmax
        integer i1, i2
 ! impose conductor behind the left reflecting wall
        
        xmin=1.
        xmax=leftwall-10. !global coordinates
        i1=iloc(int(xmin)) !convert to local index
        i2=iloc(int(xmax))

        if(i1 .ne. i2) then  !this means we are on the right CPU
           ey(i1:i2,:,:)=0.
           ez(i1:i2,:,:)=0.
        endif

        xmin=mx0-2.
        xmax=mx0
        i1=iloc(int(xmin)) !convert to local index
        i2=iloc(int(xmax))

        if(i1 .ne. i2) then 
           bx(i1:i2,:,:)=binit*cos(btheta) 
           by(i1:i2,:,:)=binit*sin(btheta)*sin(bphi)
           bz(i1:i2,:,:)=binit*sin(btheta)*cos(bphi)
           ex(i1:i2,:,:)=0.				
           ey(i1:i2,:,:)=(-beta)*bz(i1:i2,:,:) 
           ez(i1:i2,:,:)=-(-beta)*by(i1:i2,:,:)    
        endif

#ifdef turbulence
		!if(modulo(rank,sizex).eq.sizex-1)then
		!	call evaluate_turbulence_e_right_boundary(time)
		!	call evaluate_turbulence_b_right_boundary(time)
		!end if
#endif

	end subroutine field_bc_user


	real function evaluate_turbulent_b(kx, ky, kz)
		real turbulentE, Bamplitude, pi
		real kw, v
		real kx, ky, kz, kxy

#ifdef twoD
		kz = 0
#endif
		kw = sqrt(kx*kx + ky*ky + kz*kz);


#ifdef twoD
		Bamplitude = 1.0/sqrt((kw**(8.0/3.0)))
#else
		Bamplitude = 1.0/sqrt((kw**(11.0/3.0)))
#endif
		!print *, 'Binit', Binit
		!print *, 'Bamplitude', Bamplitude
		evaluate_turbulent_b = Bamplitude
	end function evaluate_turbulent_b

	real function evaluate_turbulent_b_slab(kx, ky, kz)
		real turbulentE, Bamplitude, pi
		real kw, v
		real kx, ky, kz, kxy

#ifdef twoD
		kz = 0
#endif
		kw = sqrt(kx*kx + ky*ky + kz*kz);

#ifdef twoD
		Bamplitude = 1.0/sqrt((kw**(8.0/3.0)))
#else
		Bamplitude = 1.0/sqrt((kw**(11.0/3.0)))
#endif
		!print *, 'Binit', Binit
		!print *, 'Bamplitude', Bamplitude
		evaluate_turbulent_b_slab = Bamplitude
	end function evaluate_turbulent_b_slab

	real function evaluate_turbulent_b_2d(kx, ky, kz)
		real turbulentE, Bamplitude, pi
		real kw, v
		real kx, ky, kz, kxy

#ifdef twoD
		kz = 0
#endif
		kw = sqrt(kx*kx + ky*ky + kz*kz);

#ifdef twoD
		Bamplitude = 1.0/sqrt((kw**(8.0/3.0)))
#else
		Bamplitude = 1.0/sqrt((kw**(11.0/3.0)))
#endif
		!print *, 'Binit', Binit
		!print *, 'Bamplitude', Bamplitude
		evaluate_turbulent_b_2d = Bamplitude
	end function evaluate_turbulent_b_2d

!------------------------------------------------------------------------------

	subroutine particle_bc_user()
	implicit none
	real invgam, walloc, gammawall, betawall, gamma, walloc0,t0,t1,q0
        real x0,y0,z0, tfrac, xcolis, ycolis, zcolis
	integer n1,i0,i1, iter
        logical in

	!loop over particles to check if they crossed special boundaries, like reflecting walls
	!outflow and periodic conditions are handled automatically in deposit_particles
	!
	!This routine is called after the mover and before deposit, thus allowing to avoid
        ! charge deposition behind walls. 
	
        if(wall) then 
           gammawall=1.
           betawall=0.
           walloc=leftwall

	   do iter=1,2
              if(iter.eq.1) then 
                  i0=1
                  i1=ions
                  q0=qi
	      else
		 i0=maxhlf+1
		 i1=maxhlf+lecs
                 q0=qe
	      endif
	      
	      do n1=i0,i1
              if(xglob(p(n1)%x) .lt. walloc) then 
                 gamma=sqrt(1+(p(n1)%u**2+p(n1)%v**2+p(n1)%w**2))
	      
	      !this algorithm ignores change in y and z coordinates
	      !during the scattering. Including it can result in rare
	      !conditions where a particle gets stuck in the ghost zones. 
	      !This can be improved. 

	      !unwind x location of particle
	      x0=p(n1)%x-p(n1)%u/gamma*c
	      y0=p(n1)%y !-p(n1)%v/gamma*c
	      z0=p(n1)%z !-p(n1)%w/gamma*c

	      !unwind wall location
	      walloc0=walloc-betawall*c-mxcum
	      
	      !where did they meet?
	      tfrac=abs((x0-walloc0)/(betawall*c-p(n1)%u/gamma*c))

		 xcolis=x0+p(n1)%u/gamma*c*tfrac
		 ycolis=y0	!+p(n1)%v/gamma*c*tfrac
		 zcolis=z0	!+p(n1)%w/gamma*c*tfrac

	      !deposit current upto intersection
		 q=p(n1)%ch*q0 
		 call zigzag(xcolis,ycolis,zcolis,x0,y0,z0,in)

	      !reset particle momentum, getting a kick from the wall
		 p(n1)%u=gammawall**2*gamma*(2*betawall - p(n1)%u/gamma*(1 + betawall**2))
		 gamma=sqrt(1+(p(n1)%u**2+p(n1)%v**2+p(n1)%w**2))
	      
		 tfrac=min(abs((p(n1)%x-xcolis)/max(abs(p(n1)%x-x0),1e-9)),1.)
              !move particle from the wall position with the new velocity
 
		 p(n1)%x = xcolis + p(n1)%u/gamma*c * tfrac
		 p(n1)%y = ycolis !+ p(n1)%v/gamma*c * tfrac
		 p(n1)%z = zcolis !+ p(n1)%w/gamma*c * tfrac
	     
!now clean up the piece of trajectory behind the wall, 
!that deposit_particles will be adding when it 
!unwinds the position of the particle by the full timestep. 
	      
		 q=-q
		 call zigzag(xcolis,ycolis,zcolis,p(n1)%x-p(n1)%u/gamma*c, & 
		 p(n1)%y-p(n1)%v/gamma*c,p(n1)%z-p(n1)%w/gamma*c,in)
               endif !xglob < walloc
	      enddo
	   enddo
    endif ! if(wall)

	end subroutine particle_bc_user

!-------------------------------------------------------------------------------
! 				subroutine shift_domain_user
!										
! shift fields and particles backward
! only used for some shock problems. Ignore but don't delete. 
!-------------------------------------------------------------------------------
 subroutine shift_domain_user
   implicit none

 end subroutine shift_domain_user

subroutine init_turbulent_field
	implicit none
	real B0x, B0y, B0z, E0x, E0y, E0z
	real pi;


	real turbulenceEnergyFraction

	integer, dimension(8) :: values
	integer ierr



	call date_and_time(VALUES=values)
	turbulenceSeed = values(8) + values(7)*1000
	print *, 'seed', turbulenceSeed
	call MPI_Bcast(turbulenceSeed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	print *, 'seed common', turbulenceSeed
	call srand(turbulenceSeed)
	!print *, rand(), rand(), rand()

	print *, 'start initializing turbulence'

	print *, mx0, my0, mz0

	maxTurbulentLambdaX = 2000;
	minTurbulentLambdaX = 500;
	maxTurbulentLambdaY = 2000;
	minTurbulentLambdaY = 500;
	maxTurbulentLambdaZ = 2000;
	minTurbulentLambdaZ = 500;
	turbulenceEnergyFraction = 0.5
	
	call init_turbulent_field_isotropic_plasma_frame(turbulenceEnergyFraction)

	print *, 'finish initializing turbulence'



	
end subroutine init_turbulent_field

subroutine init_turbulent_field_maltese_slab_plasma_frame(turbulenceEnergyFraction)
	real turbulenceEnergyFraction

	integer maxKx, maxKy, maxKz
	integer :: i, j, k, ki, kj, kk
	real kw
	real kx, ky, kz, kyz
	real phase1, phase2
	real cosTheta, sinTheta, cosPhi, sinPhi
	real Bturbulent
	real kmultr
	real localB1, localB2
	real turbulenceEnergy
	real slabFraction
	real slabEnergy
	real restEnergy
real maltAngle
real maltAngleSlab
real maltAngle2d

	real tempB, tempBx, tempBn, tempB0, tempB0x, tempB0n;
real Bxreg, Byreg, Bzreg, Exreg, Eyreg, Ezreg
	real Btheta0;
	real pi;

	!print *,maltAngle
	!print *, maltAngleSlab

	pi = 2.0*acos(0.0);

	tempB = 1.0;
	tempBx = tempB*cos(Btheta)
	tempBn = tempB*sin(Btheta);
	tempB0x = tempBx;
	tempB0n = tempBn/gamma0;
	Btheta0 = atan2(tempB0n, tempB0x);
	tempB0 = sqrt(tempB0x*tempB0x + tempB0n*tempB0n);
	print *,'theta0', Btheta0*180.0/pi;


	maxKx = maxTurbulentLambdaX/minTurbulentLambdaX;
	maxKy = maxTurbulentLambdaY/minTurbulentLambdaY;
	maxKz = maxTurbulentLambdaZ/minTurbulentLambdaZ;

	!maxKy = 0

	print *, maxKx, maxKy, maxKz

	!turulence energy fraction correction in plasma frame

	maltAngleSlab = 10.0*pi/180.0
	maltAngle2d = 10.0*pi/180.0

	slabFraction = 0.2;
	turbulenceEnergy = 0;
	slabEnergy = 0;
	restEnergy = 0;

#ifdef twoD
	maxKz = 1;
#endif


	do ki = 0, maxKx

		do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
			!do kk = 0, mz0-5
		do kk = 0, maxKz
#endif

			if ((ki + kj + kk) .ne. 0) then


				kx = ki*2*pi/maxTurbulentLambdaX;
				ky = kj*2*pi/maxTurbulentLambdaY;
				kz = kk*2*pi/maxTurbulentLambdaZ;

				maltAngle = acos(ky/sqrt(kx*kx + ky*ky + kz*kz))
				if(maltAngle < maltAngleSlab) then
					!print *, '2'
					Bturbulent = evaluate_turbulent_b_slab(kx, ky, kz);

					slabEnergy = slabEnergy + Bturbulent*Bturbulent;
				else if((pi/2.0) - maltAngle < maltAngle2d) then
					!print *, '3'
					Bturbulent = evaluate_turbulent_b_2d(kx, ky, kz);

					restEnergy = restEnergy + Bturbulent*Bturbulent;
				end if
			endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo

	print *,'energies', slabEnergy, restEnergy

	if(slabEnergy > 0) then
		slabFieldCorrection = sqrt(slabFraction*restEnergy/((1.0 - slabFraction)*slabEnergy))
	else
		slabFieldCorrection = 1.0;
	endif

	turbulenceEnergy = restEnergy + slabFieldCorrection*slabFieldCorrection*slabEnergy;



	if (turbulenceEnergy > 0) then
		turbulenceFieldCorrection = sqrt(turbulenceEnergyFraction*tempB0*tempB0/((1.0 - turbulenceEnergyFraction)*turbulenceEnergy));
	else
		turbulenceFieldCorrection = 1.0;
	endif

	print *, 'field corection updated', turbulenceFieldCorrection, slabFieldCorrection

	call srand(turbulenceSeed)

	!turbulence total energy correction in lab frame

	turbulenceEnergy = tempB*tempB;

#ifdef twoD
	maxKz = 1;
#endif


	do ki = 0, maxKx

		do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
			!do kk = 0, mz0-5
		do kk = 0, maxKz
#endif

			if ((ki + kj + kk) .ne. 0) then


				kx = ki*2*pi/maxTurbulentLambdaX;
				ky = kj*2*pi/maxTurbulentLambdaY;
				kz = kk*2*pi/maxTurbulentLambdaZ;


				kw = sqrt(kx*kx + ky*ky + kz*kz);
				kyz = sqrt(ky*ky + kz*kz);
				cosTheta = kx/kw;
				sinTheta = kyz/kw;
				if(kj + kk .ne. 0) then
					cosPhi = ky/kyz;
					sinPhi = kz/kyz;
				else
					cosPhi = 1.0
					sinPhi = 0.0
				endif

				maltAngle = acos(ky/sqrt(kx*kx + ky*ky + kz*kz))
				Bturbulent = 0
				if(maltAngle < maltAngleSlab) then
					Bturbulent = evaluate_turbulent_b_slab(kx, ky, kz)*slabFieldCorrection*turbulenceFieldCorrection;

				else if((pi/2.0) - maltAngle < maltAngle2d) then
					Bturbulent = evaluate_turbulent_b_2d(kx, ky, kz)*turbulenceFieldCorrection;

				end if

				turbulenceEnergy = turbulenceEnergy + &
	0.5*Bturbulent*Bturbulent*(gamma0*gamma0 + (sinTheta*sinTheta + cosTheta*cosTheta*gamma0*gamma0));
			endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo

	if (turbulenceEnergy > 0) then
		tempB = Binit/sqrt(turbulenceEnergy);
	else
		tempB = Binit;
	endif

	Binit = tempB

	!init fields

	Bxreg = tempB*cos(Btheta);
	Byreg = tempB*sin(Btheta)*sin(Bphi);
	Bzreg = tempB*sin(Btheta)*cos(Bphi);

	Exreg = 0;
	Eyreg = -beta*Bzreg;
	Ezreg = beta*Byreg;

	do  k=1,mz
		do  j=1,my
			do  i=1,mx
				! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
				!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

				bx(i,j,k)=Bxreg
				by(i,j,k)=Byreg
				bz(i,j,k)=Bzreg

				ex(i,j,k)=Exreg
				ey(i,j,k)=Eyreg
				ez(i,j,k)=Ezreg
			enddo
		enddo
	enddo

	do ki = 0, maxKx
		do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
				do kk = 0, maxKz
#endif
				!print *, ki, kj, kk

			if ((ki + kj + kk) .ne. 0) then

				phase1 = 2*pi*rand();
				phase2 = 2*pi*rand();


				kx = ki*2*pi/maxTurbulentLambdaX;
				ky = kj*2*pi/maxTurbulentLambdaY;
				kz = kk*2*pi/maxTurbulentLambdaZ;

				!print *,'k', kx, ky, kz


				kw = sqrt(kx*kx + ky*ky + kz*kz);
				kyz = sqrt(ky*ky + kz*kz);
				cosTheta = kx/kw;
				sinTheta = kyz/kw;
				if(kj + kk .ne. 0) then
					cosPhi = ky/kyz;
					sinPhi = kz/kyz;
				else
					cosPhi = 1.0
					sinPhi = 0.0
				endif

				maltAngle = acos(ky/sqrt(kx*kx + ky*ky + kz*kz))
				Bturbulent = 0
				if(maltAngle < maltAngleSlab) then
					Bturbulent = evaluate_turbulent_b_slab(kx, ky, kz)*slabFieldCorrection*turbulenceFieldCorrection*tempB;

				else if((pi/2.0) - maltAngle < maltAngle2d) then
					Bturbulent = evaluate_turbulent_b_2d(kx, ky, kz)*turbulenceFieldCorrection*tempB;

				end if


				!print *, 'Bturbulent', Bturbulent
				do  k=1,mz
					do  j=1,my
						do  i=1,mx
							! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
							!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

							kmultr = kx*gamma0*xglob(1.0*i) + ky*yglob(1.0*j) + kz*zglob(1.0*k)
							localB1 = Bturbulent*sin(kmultr + phase1);
							localB2 = Bturbulent*sin(kmultr + phase2);
							!localB2 = 0

							bx(i,j,k)=bx(i,j,k) - localB1*sinTheta
							by(i,j,k)=by(i,j,k) + (localB1*cosTheta*cosPhi - localB2*sinPhi)*gamma0
							bz(i,j,k)=bz(i,j,k) + (localB1*cosTheta*sinPhi + localB2*cosPhi)*gamma0

							ex(i,j,k)=ex(i,j,k)
							ey(i,j,k)=ey(i,j,k) - beta*gamma0*(localB1*cosTheta*sinPhi + localB2*cosPhi)
							ez(i,j,k)=ez(i,j,k) + beta*gamma0*(localB1*cosTheta*cosPhi - localB2*sinPhi)



							!bx(i,j,k)=bx(i,j,k) - localB1*cosTheta*cosPhi + localB2*sinPhi;
							!by(i,j,k)=by(i,j,k) - localB1*cosTheta*sinPhi - localB2*cosPhi;
							!bz(i,j,k)=bz(i,j,k) + localB1*sinTheta;

							!ex(i,j,k)=ex(i,j,k);
							!ey(i,j,k)=ey(i,j,k) - beta*(localB1*sinTheta);
							!ez(i,j,k)=ez(i,j,k) + beta*(- localB1*cosTheta*sinPhi - localB2*cosPhi);
						enddo
					enddo
				enddo
			endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo

end subroutine init_turbulent_field_maltese_slab_plasma_frame

subroutine init_turbulent_field_maltese_slab_lab_frame(turbulenceEnergyFraction)
			real turbulenceEnergyFraction

			integer maxKx, maxKy, maxKz
			integer :: i, j, k, ki, kj, kk
			real kw
			real kx, ky, kz, kyz
			real phase1, phase2
			real cosTheta, sinTheta, cosPhi, sinPhi
			real Bturbulent
			real kmultr
			real localB1, localB2
			real turbulenceEnergy
			real slabFraction
			real slabEnergy
			real restEnergy
			real maltAngle
			real maltAngleSlab
			real maltAngle2d

			real tempB

			real Bxreg, Byreg, Bzreg, Exreg, Eyreg, Ezreg
			real pi;

			!print *,maltAngle
			!print *, maltAngleSlab

			pi = 2.0*acos(0.0);

			tempB = 1.0;


			maxKx = maxTurbulentLambdaX/minTurbulentLambdaX;
			maxKy = maxTurbulentLambdaY/minTurbulentLambdaY;
			maxKz = maxTurbulentLambdaZ/minTurbulentLambdaZ;

			!maxKy = 0

			print *, maxKx, maxKy, maxKz

			!turulence energy fraction correction in plasma frame

			maltAngleSlab = 10.0*pi/180.0
			maltAngle2d = 10.0*pi/180.0

			slabFraction = 0.2;
			turbulenceEnergy = 0;
			slabEnergy = 0;
			restEnergy = 0;

#ifdef twoD
	maxKz = 1;
#endif


			do ki = 0, maxKx

				do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
			!do kk = 0, mz0-5
				do kk = 0, maxKz
#endif

					if ((ki + kj + kk) .ne. 0) then


						kx = ki*2*pi/maxTurbulentLambdaX;
						ky = kj*2*pi/maxTurbulentLambdaY;
						kz = kk*2*pi/maxTurbulentLambdaZ;

						maltAngle = acos(ky/sqrt(kx*kx + ky*ky + kz*kz))
						if(maltAngle < maltAngleSlab) then
							!print *, '2'
							Bturbulent = evaluate_turbulent_b_slab(kx, ky, kz);

							slabEnergy = slabEnergy + Bturbulent*Bturbulent;
						else if((pi/2.0) - maltAngle < maltAngle2d) then
							!print *, '3'
							Bturbulent = evaluate_turbulent_b_2d(kx, ky, kz);

							restEnergy = restEnergy + Bturbulent*Bturbulent;
						end if
					endif
#ifndef twoD
			enddo
#endif
		enddo
			enddo

			print *,'energies', slabEnergy, restEnergy

			if(slabEnergy > 0) then
				slabFieldCorrection = sqrt(slabFraction*restEnergy/((1.0 - slabFraction)*slabEnergy))
			else
				slabFieldCorrection = 1.0;
			endif

			turbulenceEnergy = restEnergy + slabFieldCorrection*slabFieldCorrection*slabEnergy;



			if (turbulenceEnergy > 0) then
				turbulenceFieldCorrection = sqrt(turbulenceEnergyFraction*tempB*tempB/((1.0 - turbulenceEnergyFraction)*turbulenceEnergy));
			else
				turbulenceFieldCorrection = 1.0;
			endif

			print *, 'field corection updated', turbulenceFieldCorrection, slabFieldCorrection

			call srand(turbulenceSeed)

			!turbulence total energy correction in lab frame

			turbulenceEnergy = tempB*tempB/(1.0 - turbulenceEnergyFraction);

			if (turbulenceEnergy > 0) then
				tempB = Binit/sqrt(turbulenceEnergy);
			else
				tempB = Binit;
			endif

			Binit = tempB

			!init fields

			Bxreg = tempB*cos(Btheta);
			Byreg = tempB*sin(Btheta)*sin(Bphi);
			Bzreg = tempB*sin(Btheta)*cos(Bphi);

			Exreg = 0;
			Eyreg = -beta*Bzreg;
			Ezreg = beta*Byreg;

			do  k=1,mz
				do  j=1,my
					do  i=1,mx
						! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
						!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

						bx(i,j,k)=Bxreg
						by(i,j,k)=Byreg
						bz(i,j,k)=Bzreg

						ex(i,j,k)=Exreg
						ey(i,j,k)=Eyreg
						ez(i,j,k)=Ezreg
					enddo
				enddo
			enddo

			do ki = 0, maxKx
				do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
				do kk = 0, maxKz
#endif
				!print *, ki, kj, kk

					if ((ki + kj + kk) .ne. 0) then

						phase1 = 2*pi*rand();
						phase2 = 2*pi*rand();


						kx = ki*2*pi/maxTurbulentLambdaX;
						ky = kj*2*pi/maxTurbulentLambdaY;
						kz = kk*2*pi/maxTurbulentLambdaZ;

						!print *,'k', kx, ky, kz


						kw = sqrt(kx*kx + ky*ky + kz*kz);
						kyz = sqrt(ky*ky + kz*kz);
						cosTheta = kx/kw;
						sinTheta = kyz/kw;
						if(kj + kk .ne. 0) then
							cosPhi = ky/kyz;
							sinPhi = kz/kyz;
						else
							cosPhi = 1.0
							sinPhi = 0.0
						endif

						maltAngle = acos(ky/sqrt(kx*kx + ky*ky + kz*kz))
						Bturbulent = 0
						if(maltAngle < maltAngleSlab) then
							Bturbulent = evaluate_turbulent_b_slab(kx, ky, kz)*slabFieldCorrection*turbulenceFieldCorrection*tempB;

						else if((pi/2.0) - maltAngle < maltAngle2d) then
							Bturbulent = evaluate_turbulent_b_2d(kx, ky, kz)*turbulenceFieldCorrection*tempB;

						end if


						!print *, 'Bturbulent', Bturbulent
						do  k=1,mz
							do  j=1,my
								do  i=1,mx
									! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
									!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

									kmultr = kx*xglob(1.0*i) + ky*yglob(1.0*j) + kz*zglob(1.0*k)
									localB1 = Bturbulent*sin(kmultr + phase1);
									localB2 = Bturbulent*sin(kmultr + phase2);
									!localB2 = 0

									bx(i,j,k)=bx(i,j,k) - localB1*sinTheta
									by(i,j,k)=by(i,j,k) + (localB1*cosTheta*cosPhi - localB2*sinPhi)
									bz(i,j,k)=bz(i,j,k) + (localB1*cosTheta*sinPhi + localB2*cosPhi)

									ex(i,j,k)=ex(i,j,k)
									ey(i,j,k)=ey(i,j,k) - beta*(localB1*cosTheta*sinPhi + localB2*cosPhi)
									ez(i,j,k)=ez(i,j,k) + beta*(localB1*cosTheta*cosPhi - localB2*sinPhi)



									!bx(i,j,k)=bx(i,j,k) - localB1*cosTheta*cosPhi + localB2*sinPhi;
									!by(i,j,k)=by(i,j,k) - localB1*cosTheta*sinPhi - localB2*cosPhi;
									!bz(i,j,k)=bz(i,j,k) + localB1*sinTheta;

									!ex(i,j,k)=ex(i,j,k);
									!ey(i,j,k)=ey(i,j,k) - beta*(localB1*sinTheta);
									!ez(i,j,k)=ez(i,j,k) + beta*(- localB1*cosTheta*sinPhi - localB2*cosPhi);
								enddo
							enddo
						enddo
					endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo
end subroutine init_turbulent_field_maltese_slab_lab_frame


subroutine init_turbulent_field_slab_plasma_frame(turbulenceEnergyFraction)
	real turbulenceEnergyFraction

	integer maxKx, maxKy, maxKz
	integer :: i, j, k, ki, kj, kk
	real kw
	real kx, ky, kz, kyz
	real phase1, phase2
	real cosTheta, sinTheta, cosPhi, sinPhi
	real Bturbulent
	real kmultr
	real localB1, localB2
	real turbulenceEnergy
	real slabFraction
	real slabEnergy
	real restEnergy

	real tempB, tempBx, tempBn, tempB0, tempB0x, tempB0n;
	real Bxreg, Byreg, Bzreg, Exreg, Eyreg, Ezreg
	real Btheta0;
	real pi;

	pi = 2.0*acos(0.0);

	tempB = 1.0;
	tempBx = tempB*cos(Btheta)
	tempBn = tempB*sin(Btheta);
	tempB0x = tempBx;
	tempB0n = tempBn/gamma0;
	Btheta0 = atan2(tempB0n, tempB0x);
	tempB0 = sqrt(tempB0x*tempB0x + tempB0n*tempB0n);
	print *,'theta0', Btheta0*180.0/pi;


	maxKx = maxTurbulentLambdaX/minTurbulentLambdaX;
	maxKy = maxTurbulentLambdaY/minTurbulentLambdaY;
	maxKz = maxTurbulentLambdaZ/minTurbulentLambdaZ;

	!maxKy = 0

	print *, maxKx, maxKy, maxKz

	!turulence energy fraction correction in plasma frame

	slabFraction = 0.2;
	turbulenceEnergy = 0;
	slabEnergy = 0;
	restEnergy = 0;

#ifdef twoD
	maxKz = 1;
#endif



	do ki = 0, maxKx

		do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
			!do kk = 0, mz0-5
		do kk = 0, maxKz
#endif

			if ((ki + kj + kk) .ne. 0) then


				kx = ki*2*pi/maxTurbulentLambdaX;
				ky = kj*2*pi/maxTurbulentLambdaY;
				kz = kk*2*pi/maxTurbulentLambdaZ;

				if(kk + ki .eq. 0) then
					!print *, '2'
					Bturbulent = evaluate_turbulent_b_slab(kx, ky, kz);

					slabEnergy = slabEnergy + Bturbulent*Bturbulent;
				else if(kj .eq. 0) then
					!print *, '3'
					Bturbulent = evaluate_turbulent_b_2d(kx, ky, kz);

					restEnergy = restEnergy + Bturbulent*Bturbulent;
				end if
			endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo

	print *,'energies', slabEnergy, restEnergy

	if(slabEnergy > 0) then
		slabFieldCorrection = sqrt(slabFraction*restEnergy/((1.0 - slabFraction)*slabEnergy))
	else
		slabFieldCorrection = 1.0;
	endif

	turbulenceEnergy = restEnergy + slabFieldCorrection*slabFieldCorrection*slabEnergy;



	if (turbulenceEnergy > 0) then
		turbulenceFieldCorrection = sqrt(turbulenceEnergyFraction*tempB0*tempB0/((1.0 - turbulenceEnergyFraction)*turbulenceEnergy));
	else
		turbulenceFieldCorrection = 1.0;
	endif

	print *, 'field corection updated', turbulenceFieldCorrection, slabFieldCorrection

	call srand(turbulenceSeed)

	!turbulence total energy correction in lab frame

	turbulenceEnergy = tempB*tempB;

#ifdef twoD
	maxKz = 1;
#endif


	do ki = 0, maxKx

		do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
			!do kk = 0, mz0-5
		do kk = 0, maxKz
#endif

			if ((ki + kj + kk) .ne. 0) then


				kx = ki*2*pi/maxTurbulentLambdaX;
				ky = kj*2*pi/maxTurbulentLambdaY;
				kz = kk*2*pi/maxTurbulentLambdaZ;


				kw = sqrt(kx*kx + ky*ky + kz*kz);
				kyz = sqrt(ky*ky + kz*kz);
				cosTheta = kx/kw;
				sinTheta = kyz/kw;
				if(kj + kk .ne. 0) then
					cosPhi = ky/kyz;
					sinPhi = kz/kyz;
				else
					cosPhi = 1.0
					sinPhi = 0.0
				endif

				Bturbulent = 0
				if(kk + ki .eq. 0) then
					!print *, '2'
					Bturbulent = evaluate_turbulent_b_slab(kx, ky, kz)*slabFieldCorrection*turbulenceFieldCorrection;

					slabEnergy = slabEnergy + Bturbulent*Bturbulent;
				else if(kj .eq. 0) then
					!print *, '3'
					Bturbulent = evaluate_turbulent_b_2d(kx, ky, kz)*turbulenceFieldCorrection;

					restEnergy = restEnergy + Bturbulent*Bturbulent;
				end if

				turbulenceEnergy = turbulenceEnergy + &
				0.5*Bturbulent*Bturbulent*(gamma0*gamma0 + (sinTheta*sinTheta + cosTheta*cosTheta*gamma0*gamma0));
			endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo

	if (turbulenceEnergy > 0) then
		tempB = Binit/sqrt(turbulenceEnergy);
	else
		tempB = Binit;
	endif

	Binit = tempB

	!init fields

	Bxreg = tempB*cos(Btheta);
	Byreg = tempB*sin(Btheta)*sin(Bphi);
	Bzreg = tempB*sin(Btheta)*cos(Bphi);

	Exreg = 0;
	Eyreg = -beta*Bzreg;
	Ezreg = beta*Byreg;

	do  k=1,mz
		do  j=1,my
			do  i=1,mx
				! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
				!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

				bx(i,j,k)=Bxreg
				by(i,j,k)=Byreg
				bz(i,j,k)=Bzreg

				ex(i,j,k)=Exreg
				ey(i,j,k)=Eyreg
				ez(i,j,k)=Ezreg
			enddo
		enddo
	enddo

	do ki = 0, maxKx
		do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
				do kk = 0, maxKz
#endif
				!print *, ki, kj, kk

			if ((ki + kj + kk) .ne. 0) then

				phase1 = 2*pi*rand();
				phase2 = 2*pi*rand();


				kx = ki*2*pi/maxTurbulentLambdaX;
				ky = kj*2*pi/maxTurbulentLambdaY;
				kz = kk*2*pi/maxTurbulentLambdaZ;

				!print *,'k', kx, ky, kz


				kw = sqrt(kx*kx + ky*ky + kz*kz);
				kyz = sqrt(ky*ky + kz*kz);
				cosTheta = kx/kw;
				sinTheta = kyz/kw;
				if(kj + kk .ne. 0) then
					cosPhi = ky/kyz;
					sinPhi = kz/kyz;
				else
					cosPhi = 1.0
					sinPhi = 0.0
				endif

				Bturbulent = 0
				if(kk + ki .eq. 0) then
					!print *, '2'
					Bturbulent = evaluate_turbulent_b_slab(kx, ky, kz)*slabFieldCorrection*turbulenceFieldCorrection*tempB;

					slabEnergy = slabEnergy + Bturbulent*Bturbulent;
				else if(kj .eq. 0) then
					!print *, '3'
					Bturbulent = evaluate_turbulent_b_2d(kx, ky, kz)*turbulenceFieldCorrection*tempB;

					restEnergy = restEnergy + Bturbulent*Bturbulent;
				end if


				!print *, 'Bturbulent', Bturbulent
				do  k=1,mz
					do  j=1,my
						do  i=1,mx
							! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
							!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

							kmultr = kx*gamma0*xglob(1.0*i) + ky*yglob(1.0*j) + kz*zglob(1.0*k)
							localB1 = Bturbulent*sin(kmultr + phase1);
							localB2 = Bturbulent*sin(kmultr + phase2);
							!localB2 = 0

							bx(i,j,k)=bx(i,j,k) - localB1*sinTheta
							by(i,j,k)=by(i,j,k) + (localB1*cosTheta*cosPhi - localB2*sinPhi)*gamma0
							bz(i,j,k)=bz(i,j,k) + (localB1*cosTheta*sinPhi + localB2*cosPhi)*gamma0

							ex(i,j,k)=ex(i,j,k)
							ey(i,j,k)=ey(i,j,k) - beta*gamma0*(localB1*cosTheta*sinPhi + localB2*cosPhi)
							ez(i,j,k)=ez(i,j,k) + beta*gamma0*(localB1*cosTheta*cosPhi - localB2*sinPhi)



							!bx(i,j,k)=bx(i,j,k) - localB1*cosTheta*cosPhi + localB2*sinPhi;
							!by(i,j,k)=by(i,j,k) - localB1*cosTheta*sinPhi - localB2*cosPhi;
							!bz(i,j,k)=bz(i,j,k) + localB1*sinTheta;

							!ex(i,j,k)=ex(i,j,k);
							!ey(i,j,k)=ey(i,j,k) - beta*(localB1*sinTheta);
							!ez(i,j,k)=ez(i,j,k) + beta*(- localB1*cosTheta*sinPhi - localB2*cosPhi);
						enddo
					enddo
				enddo
			endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo
end subroutine init_turbulent_field_slab_plasma_frame

subroutine init_turbulent_field_slab_lab_frame(turbulenceEnergyFraction)
			real turbulenceEnergyFraction

			integer maxKx, maxKy, maxKz
			integer :: i, j, k, ki, kj, kk
			real kw
			real kx, ky, kz, kyz
			real phase1, phase2
			real cosTheta, sinTheta, cosPhi, sinPhi
			real Bturbulent
			real kmultr
			real localB1, localB2
			real turbulenceEnergy
			real slabFraction
			real slabEnergy
			real restEnergy

			real tempB
			real Bxreg, Byreg, Bzreg, Exreg, Eyreg, Ezreg
			real Btheta0;
			real pi;

			pi = 2.0*acos(0.0);

			tempB = 1.0;



			maxKx = maxTurbulentLambdaX/minTurbulentLambdaX;
			maxKy = maxTurbulentLambdaY/minTurbulentLambdaY;
			maxKz = maxTurbulentLambdaZ/minTurbulentLambdaZ;

			!maxKy = 0

			print *, maxKx, maxKy, maxKz

			!turulence energy fraction correction in plasma frame

			slabFraction = 0.2;
			turbulenceEnergy = 0;
			slabEnergy = 0;
			restEnergy = 0;

#ifdef twoD
	maxKz = 1;
#endif


			do ki = 0, maxKx

				do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
			!do kk = 0, mz0-5
				do kk = 0, maxKz
#endif

					if ((ki + kj + kk) .ne. 0) then


						kx = ki*2*pi/maxTurbulentLambdaX;
						ky = kj*2*pi/maxTurbulentLambdaY;
						kz = kk*2*pi/maxTurbulentLambdaZ;

						if(kk + ki .eq. 0) then
							!print *, '2'
							Bturbulent = evaluate_turbulent_b_slab(kx, ky, kz);

							slabEnergy = slabEnergy + Bturbulent*Bturbulent;
						else if(kj .eq. 0) then
							!print *, '3'
							Bturbulent = evaluate_turbulent_b_2d(kx, ky, kz);

							restEnergy = restEnergy + Bturbulent*Bturbulent;
						end if
					endif
#ifndef twoD
			enddo
#endif
		enddo
			enddo

			print *,'energies', slabEnergy, restEnergy

			if(slabEnergy > 0) then
				slabFieldCorrection = sqrt(slabFraction*restEnergy/((1.0 - slabFraction)*slabEnergy))
			else
				slabFieldCorrection = 1.0;
			endif

			turbulenceEnergy = restEnergy + slabFieldCorrection*slabFieldCorrection*slabEnergy;



			if (turbulenceEnergy > 0) then
				turbulenceFieldCorrection = sqrt(turbulenceEnergyFraction*tempB*tempB/((1.0 - turbulenceEnergyFraction)*turbulenceEnergy));
			else
				turbulenceFieldCorrection = 1.0;
			endif

			print *, 'field corection updated', turbulenceFieldCorrection, slabFieldCorrection

			call srand(turbulenceSeed)

			!turbulence total energy correction in lab frame

			turbulenceEnergy = tempB*tempB/(1.0 - turbulenceEnergyFraction);

			if (turbulenceEnergy > 0) then
				tempB = Binit/sqrt(turbulenceEnergy);
			else
				tempB = Binit;
			endif

			Binit = tempB

			!init fields

			Bxreg = tempB*cos(Btheta);
			Byreg = tempB*sin(Btheta)*sin(Bphi);
			Bzreg = tempB*sin(Btheta)*cos(Bphi);

			Exreg = 0;
			Eyreg = -beta*Bzreg;
			Ezreg = beta*Byreg;

			do  k=1,mz
				do  j=1,my
					do  i=1,mx
						! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
						!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

						bx(i,j,k)=Bxreg
						by(i,j,k)=Byreg
						bz(i,j,k)=Bzreg

						ex(i,j,k)=Exreg
						ey(i,j,k)=Eyreg
						ez(i,j,k)=Ezreg
					enddo
				enddo
			enddo

			do ki = 0, maxKx
				do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
				do kk = 0, maxKz
#endif
				!print *, ki, kj, kk

					if ((ki + kj + kk) .ne. 0) then

						phase1 = 2*pi*rand();
						phase2 = 2*pi*rand();


						kx = ki*2*pi/maxTurbulentLambdaX;
						ky = kj*2*pi/maxTurbulentLambdaY;
						kz = kk*2*pi/maxTurbulentLambdaZ;

						!print *,'k', kx, ky, kz


						kw = sqrt(kx*kx + ky*ky + kz*kz);
						kyz = sqrt(ky*ky + kz*kz);
						cosTheta = kx/kw;
						sinTheta = kyz/kw;
						if(kj + kk .ne. 0) then
							cosPhi = ky/kyz;
							sinPhi = kz/kyz;
						else
							cosPhi = 1.0
							sinPhi = 0.0
						endif

						Bturbulent = 0
						if(kk + ki .eq. 0) then
							!print *, '2'
							Bturbulent = evaluate_turbulent_b_slab(kx, ky, kz)*slabFieldCorrection*turbulenceFieldCorrection*tempB;

							slabEnergy = slabEnergy + Bturbulent*Bturbulent;
						else if(kj .eq. 0) then
							!print *, '3'
							Bturbulent = evaluate_turbulent_b_2d(kx, ky, kz)*turbulenceFieldCorrection*tempB;

							restEnergy = restEnergy + Bturbulent*Bturbulent;
						end if


						!print *, 'Bturbulent', Bturbulent
						do  k=1,mz
							do  j=1,my
								do  i=1,mx
									! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
									!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

									kmultr = kx*xglob(1.0*i) + ky*yglob(1.0*j) + kz*zglob(1.0*k)
									localB1 = Bturbulent*sin(kmultr + phase1);
									localB2 = Bturbulent*sin(kmultr + phase2);
									!localB2 = 0

									bx(i,j,k)=bx(i,j,k) - localB1*sinTheta
									by(i,j,k)=by(i,j,k) + (localB1*cosTheta*cosPhi - localB2*sinPhi)
									bz(i,j,k)=bz(i,j,k) + (localB1*cosTheta*sinPhi + localB2*cosPhi)

									ex(i,j,k)=ex(i,j,k)
									ey(i,j,k)=ey(i,j,k) - beta*(localB1*cosTheta*sinPhi + localB2*cosPhi)
									ez(i,j,k)=ez(i,j,k) + beta*(localB1*cosTheta*cosPhi - localB2*sinPhi)



									!bx(i,j,k)=bx(i,j,k) - localB1*cosTheta*cosPhi + localB2*sinPhi;
									!by(i,j,k)=by(i,j,k) - localB1*cosTheta*sinPhi - localB2*cosPhi;
									!bz(i,j,k)=bz(i,j,k) + localB1*sinTheta;

									!ex(i,j,k)=ex(i,j,k);
									!ey(i,j,k)=ey(i,j,k) - beta*(localB1*sinTheta);
									!ez(i,j,k)=ez(i,j,k) + beta*(- localB1*cosTheta*sinPhi - localB2*cosPhi);
								enddo
							enddo
						enddo
					endif
#ifndef twoD
			enddo
#endif
		enddo
			enddo
end subroutine init_turbulent_field_slab_lab_frame


subroutine init_turbulent_field_isotropic_plasma_frame(turbulenceEnergyFraction)
	real turbulenceEnergyFraction

	integer maxKx, maxKy, maxKz
	integer :: i, j, k, ki, kj, kk
	real kw
	real kx, ky, kz, kyz
	real phase1, phase2
	real cosTheta, sinTheta, cosPhi, sinPhi
	real Bturbulent
	real kmultr
	real localB1, localB2
	real turbulenceEnergy
	real tempB, tempBx, tempBn, tempB0, tempB0x, tempB0n;
real Bxreg, Byreg, Bzreg, Exreg, Eyreg, Ezreg
	real Btheta0;
	real pi;

	pi = 2.0*acos(0.0);

	tempB = 1.0;
	tempBx = tempB*cos(Btheta)
	tempBn = tempB*sin(Btheta);
	tempB0x = tempBx;
	tempB0n = tempBn/gamma0;
	Btheta0 = atan2(tempB0n, tempB0x);
	tempB0 = sqrt(tempB0x*tempB0x + tempB0n*tempB0n);
	print *,'theta0', Btheta0*180.0/pi;


	maxKx = maxTurbulentLambdaX/minTurbulentLambdaX;
	maxKy = maxTurbulentLambdaY/minTurbulentLambdaY;
	maxKz = maxTurbulentLambdaZ/minTurbulentLambdaZ;

	!maxKy = 0

	print *, maxKx, maxKy, maxKz

	!turulence energy fraction correction in plasma frame

	turbulenceEnergy = 0;

#ifdef twoD
	maxKz = 1;
#endif


		do ki = 0, maxKx

			do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
			!do kk = 0, mz0-5
				do kk = 0, maxKz
#endif

					if ((ki + kj + kk) .ne. 0) then


						kx = ki*2*pi/maxTurbulentLambdaX;
						ky = kj*2*pi/maxTurbulentLambdaY;
						kz = kk*2*pi/maxTurbulentLambdaZ;


						Bturbulent = evaluate_turbulent_b(kx, ky, kz);

						turbulenceEnergy = turbulenceEnergy + Bturbulent*Bturbulent;
					endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo

	if (turbulenceEnergy > 0) then
		turbulenceFieldCorrection = sqrt(turbulenceEnergyFraction*tempB0*tempB0/((1.0 - turbulenceEnergyFraction)*turbulenceEnergy));
	else
		turbulenceFieldCorrection = 1.0;
	endif

	!turbulence total energy correction in lab frame

	turbulenceEnergy = tempB*tempB;

#ifdef twoD
	maxKz = 1;
#endif


	do ki = 0, maxKx

		do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
			!do kk = 0, mz0-5
		do kk = 0, maxKz
#endif

			if ((ki + kj + kk) .ne. 0) then


				kx = ki*2*pi/maxTurbulentLambdaX;
				ky = kj*2*pi/maxTurbulentLambdaY;
				kz = kk*2*pi/maxTurbulentLambdaZ;


				kw = sqrt(kx*kx + ky*ky + kz*kz);
				kyz = sqrt(ky*ky + kz*kz);
				cosTheta = kx/kw;
				sinTheta = kyz/kw;
				if(kj + kk .ne. 0) then
					cosPhi = ky/kyz;
					sinPhi = kz/kyz;
				else
					cosPhi = 1.0
					sinPhi = 0.0
				endif

				Bturbulent = evaluate_turbulent_b(kx, ky, kz)*turbulenceFieldCorrection;

				turbulenceEnergy = turbulenceEnergy + &
				0.5*Bturbulent*Bturbulent*(gamma0*gamma0 + (sinTheta*sinTheta + cosTheta*cosTheta*gamma0*gamma0));
			endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo

	if (turbulenceEnergy > 0) then
		tempB = Binit/sqrt(turbulenceEnergy);
	else
		tempB = Binit;
	endif

	Binit = tempB

	!init fields

	Bxreg = tempB*cos(Btheta);
	Byreg = tempB*sin(Btheta)*sin(Bphi);
	Bzreg = tempB*sin(Btheta)*cos(Bphi);

	Exreg = 0;
	Eyreg = -beta*Bzreg;
	Ezreg = beta*Byreg;

	do  k=1,mz
		do  j=1,my
			do  i=1,mx
				! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
				!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

				bx(i,j,k)=Bxreg
				by(i,j,k)=Byreg
				bz(i,j,k)=Bzreg

				ex(i,j,k)=Exreg
				ey(i,j,k)=Eyreg
				ez(i,j,k)=Ezreg
			enddo
		enddo
	enddo

	do ki = 0, maxKx
		do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
				do kk = 0, maxKz
#endif
				!print *, ki, kj, kk

					if ((ki + kj + kk) .ne. 0) then

						phase1 = 2*pi*rand();
						phase2 = 2*pi*rand();


						kx = ki*2*pi/maxTurbulentLambdaX;
						ky = kj*2*pi/maxTurbulentLambdaY;
						kz = kk*2*pi/maxTurbulentLambdaZ;

						!print *,'k', kx, ky, kz


						kw = sqrt(kx*kx + ky*ky + kz*kz);
						kyz = sqrt(ky*ky + kz*kz);
						cosTheta = kx/kw;
						sinTheta = kyz/kw;
						if(kj + kk .ne. 0) then
							cosPhi = ky/kyz;
							sinPhi = kz/kyz;
						else
							cosPhi = 1.0
							sinPhi = 0.0
						endif

						Bturbulent = evaluate_turbulent_b(kx, ky, kz)*tempB*turbulenceFieldCorrection;


						!print *, 'Bturbulent', Bturbulent
						do  k=1,mz
							do  j=1,my
								do  i=1,mx
									! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
									!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

									kmultr = kx*gamma0*xglob(1.0*i) + ky*yglob(1.0*j) + kz*zglob(1.0*k)
									localB1 = Bturbulent*sin(kmultr + phase1);
									localB2 = Bturbulent*sin(kmultr + phase2);
									!localB2 = 0

									bx(i,j,k)=bx(i,j,k) - localB1*sinTheta
									by(i,j,k)=by(i,j,k) + (localB1*cosTheta*cosPhi - localB2*sinPhi)*gamma0
									bz(i,j,k)=bz(i,j,k) + (localB1*cosTheta*sinPhi + localB2*cosPhi)*gamma0

									ex(i,j,k)=ex(i,j,k)
									ey(i,j,k)=ey(i,j,k) - beta*gamma0*(localB1*cosTheta*sinPhi + localB2*cosPhi)
									ez(i,j,k)=ez(i,j,k) + beta*gamma0*(localB1*cosTheta*cosPhi - localB2*sinPhi)



									!bx(i,j,k)=bx(i,j,k) - localB1*cosTheta*cosPhi + localB2*sinPhi;
									!by(i,j,k)=by(i,j,k) - localB1*cosTheta*sinPhi - localB2*cosPhi;
									!bz(i,j,k)=bz(i,j,k) + localB1*sinTheta;

									!ex(i,j,k)=ex(i,j,k);
									!ey(i,j,k)=ey(i,j,k) - beta*(localB1*sinTheta);
									!ez(i,j,k)=ez(i,j,k) + beta*(- localB1*cosTheta*sinPhi - localB2*cosPhi);
								enddo
							enddo
						enddo
					endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo

end subroutine init_turbulent_field_isotropic_plasma_frame

subroutine init_turbulent_field_isotropic_lab_frame(turbulenceEnergyFraction)
			real turbulenceEnergyFraction

			integer maxKx, maxKy, maxKz
			integer :: i, j, k, ki, kj, kk
			real kw
			real kx, ky, kz, kyz
			real phase1, phase2
			real cosTheta, sinTheta, cosPhi, sinPhi
			real Bturbulent
			real kmultr
			real localB1, localB2
			real turbulenceEnergy
			real tempB
			real Bxreg, Byreg, Bzreg, Exreg, Eyreg, Ezreg
			real Btheta0;
			real pi;

			pi = 2.0*acos(0.0);

			tempB = 1.0;


			maxKx = maxTurbulentLambdaX/minTurbulentLambdaX;
			maxKy = maxTurbulentLambdaY/minTurbulentLambdaY;
			maxKz = maxTurbulentLambdaZ/minTurbulentLambdaZ;

			!maxKy = 0

			print *, maxKx, maxKy, maxKz

			!turulence energy fraction correction in plasma frame

			turbulenceEnergy = 0;

#ifdef twoD
	maxKz = 1;
#endif


			do ki = 0, maxKx

				do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
			!do kk = 0, mz0-5
				do kk = 0, maxKz
#endif

					if ((ki + kj + kk) .ne. 0) then


						kx = ki*2*pi/maxTurbulentLambdaX;
						ky = kj*2*pi/maxTurbulentLambdaY;
						kz = kk*2*pi/maxTurbulentLambdaZ;


						Bturbulent = evaluate_turbulent_b(kx, ky, kz);

						turbulenceEnergy = turbulenceEnergy + Bturbulent*Bturbulent;
					endif
#ifndef twoD
			enddo
#endif
		enddo
			enddo

			if (turbulenceEnergy > 0) then
				turbulenceFieldCorrection = sqrt(turbulenceEnergyFraction*tempB*tempB/((1.0 - turbulenceEnergyFraction)*turbulenceEnergy));
			else
				turbulenceFieldCorrection = 1.0;
			endif

			!turbulence total energy correction in lab frame

			turbulenceEnergy = tempB*tempB/(1.0 - turbulenceEnergyFraction);


			if (turbulenceEnergy > 0) then
				tempB = Binit/sqrt(turbulenceEnergy);
			else
				tempB = Binit;
			endif

			Binit = tempB

			!init fields

			Bxreg = tempB*cos(Btheta);
			Byreg = tempB*sin(Btheta)*sin(Bphi);
			Bzreg = tempB*sin(Btheta)*cos(Bphi);

			Exreg = 0;
			Eyreg = -beta*Bzreg;
			Ezreg = beta*Byreg;

			do  k=1,mz
				do  j=1,my
					do  i=1,mx
						! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
						!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

						bx(i,j,k)=Bxreg
						by(i,j,k)=Byreg
						bz(i,j,k)=Bzreg

						ex(i,j,k)=Exreg
						ey(i,j,k)=Eyreg
						ez(i,j,k)=Ezreg
					enddo
				enddo
			enddo

			do ki = 0, maxKx
				do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
				do kk = 0, maxKz
#endif
				!print *, ki, kj, kk

					if ((ki + kj + kk) .ne. 0) then

						phase1 = 2*pi*rand();
						phase2 = 2*pi*rand();


						kx = ki*2*pi/maxTurbulentLambdaX;
						ky = kj*2*pi/maxTurbulentLambdaY;
						kz = kk*2*pi/maxTurbulentLambdaZ;

						!print *,'k', kx, ky, kz


						kw = sqrt(kx*kx + ky*ky + kz*kz);
						kyz = sqrt(ky*ky + kz*kz);
						cosTheta = kx/kw;
						sinTheta = kyz/kw;
						if(kj + kk .ne. 0) then
							cosPhi = ky/kyz;
							sinPhi = kz/kyz;
						else
							cosPhi = 1.0
							sinPhi = 0.0
						endif

						Bturbulent = evaluate_turbulent_b(kx, ky, kz)*tempB*turbulenceFieldCorrection;


						!print *, 'Bturbulent', Bturbulent
						do  k=1,mz
							do  j=1,my
								do  i=1,mx
									! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
									!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

									kmultr = kx*xglob(1.0*i) + ky*yglob(1.0*j) + kz*zglob(1.0*k)
									localB1 = Bturbulent*sin(kmultr + phase1);
									localB2 = Bturbulent*sin(kmultr + phase2);
									!localB2 = 0

									bx(i,j,k)=bx(i,j,k) - localB1*sinTheta
									by(i,j,k)=by(i,j,k) + (localB1*cosTheta*cosPhi - localB2*sinPhi)
									bz(i,j,k)=bz(i,j,k) + (localB1*cosTheta*sinPhi + localB2*cosPhi)

									ex(i,j,k)=ex(i,j,k)
									ey(i,j,k)=ey(i,j,k) - beta*(localB1*cosTheta*sinPhi + localB2*cosPhi)
									ez(i,j,k)=ez(i,j,k) + beta*(localB1*cosTheta*cosPhi - localB2*sinPhi)



									!bx(i,j,k)=bx(i,j,k) - localB1*cosTheta*cosPhi + localB2*sinPhi;
									!by(i,j,k)=by(i,j,k) - localB1*cosTheta*sinPhi - localB2*cosPhi;
									!bz(i,j,k)=bz(i,j,k) + localB1*sinTheta;

									!ex(i,j,k)=ex(i,j,k);
									!ey(i,j,k)=ey(i,j,k) - beta*(localB1*sinTheta);
									!ez(i,j,k)=ez(i,j,k) + beta*(- localB1*cosTheta*sinPhi - localB2*cosPhi);
								enddo
							enddo
						enddo
					endif
#ifndef twoD
			enddo
#endif
		enddo
			enddo

end subroutine init_turbulent_field_isotropic_lab_frame


subroutine init_turbulent_field_simple_anisotropic_plasma_frame(turbulenceEnergyFraction)
	real turbulenceEnergyFraction

	integer maxKx, maxKy, maxKz
	integer :: i, j, k, ki, kj, kk
	real kw
	real kx, ky, kz, kyz
	real phase1, phase2
	real cosTheta, sinTheta, cosPhi, sinPhi
	real Bturbulent
	real kmultr
	real localB1, localB2
	real turbulenceEnergy
	real slabFraction
	real slabEnergy
	real restEnergy

	real tempB, tempBx, tempBn, tempB0, tempB0x, tempB0n;
real Bxreg, Byreg, Bzreg, Exreg, Eyreg, Ezreg
	real Btheta0;
	real pi;

	pi = 2.0*acos(0.0);

	tempB = 1.0;
	tempBx = tempB*cos(Btheta)
	tempBn = tempB*sin(Btheta);
	tempB0x = tempBx;
	tempB0n = tempBn/gamma0;
	Btheta0 = atan2(tempB0n, tempB0x);
	tempB0 = sqrt(tempB0x*tempB0x + tempB0n*tempB0n);
	print *,'theta0', Btheta0*180.0/pi;


	maxKx = maxTurbulentLambdaX/minTurbulentLambdaX;
	maxKy = maxTurbulentLambdaY/minTurbulentLambdaY;
	maxKz = maxTurbulentLambdaZ/minTurbulentLambdaZ;

	!maxKy = 0

	print *, maxKx, maxKy, maxKz

	!turulence energy fraction correction in plasma frame

	turbulenceEnergy = 0;

#ifdef twoD
	maxKz = 1;
#endif


	do ki = 0, maxKx

		do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
			!do kk = 0, mz0-5
		do kk = 0, maxKz
#endif

			if (((ki + kj + kk) .ne. 0) .and. ((kj + kk).ne.0)) then


				kx = ki*2*pi/maxTurbulentLambdaX;
				ky = kj*2*pi/maxTurbulentLambdaY;
				kz = kk*2*pi/maxTurbulentLambdaZ;


				Bturbulent = evaluate_turbulent_b(kx, ky, kz);

				turbulenceEnergy = turbulenceEnergy + Bturbulent*Bturbulent;
			endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo

	if (turbulenceEnergy > 0) then
		turbulenceFieldCorrection = sqrt(turbulenceEnergyFraction*tempB0*tempB0/((1.0 - turbulenceEnergyFraction)*turbulenceEnergy));
	else
		turbulenceFieldCorrection = 1.0;
	endif

	!turbulence total energy correction in lab frame

	turbulenceEnergy = tempB*tempB;

#ifdef twoD
	maxKz = 1;
#endif


	do ki = 0, maxKx

		do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
			!do kk = 0, mz0-5
		do kk = 0, maxKz
#endif

			if (((ki + kj + kk) .ne. 0) .and. ((kj + kk).ne.0)) then


				kx = ki*2*pi/maxTurbulentLambdaX;
				ky = kj*2*pi/maxTurbulentLambdaY;
				kz = kk*2*pi/maxTurbulentLambdaZ;


				kw = sqrt(kx*kx + ky*ky + kz*kz);
				kyz = sqrt(ky*ky + kz*kz);
				cosTheta = kx/kw;
				sinTheta = kyz/kw;
				if(kj + kk .ne. 0) then
					cosPhi = ky/kyz;
					sinPhi = kz/kyz;
				else
					cosPhi = 1.0
					sinPhi = 0.0
				endif

				Bturbulent = evaluate_turbulent_b(kx, ky, kz)*turbulenceFieldCorrection;

				turbulenceEnergy = turbulenceEnergy + &
				0.5*Bturbulent*Bturbulent*(gamma0*gamma0 + (sinTheta*sinTheta + cosTheta*cosTheta*gamma0*gamma0));
			endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo

	if (turbulenceEnergy > 0) then
		tempB = Binit/sqrt(turbulenceEnergy);
	else
		tempB = Binit;
	endif

	Binit = tempB

	!init fields

	Bxreg = tempB*cos(Btheta);
	Byreg = tempB*sin(Btheta)*sin(Bphi);
	Bzreg = tempB*sin(Btheta)*cos(Bphi);

	Exreg = 0;
	Eyreg = -beta*Bzreg;
	Ezreg = beta*Byreg;

	do  k=1,mz
		do  j=1,my
			do  i=1,mx
				! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
				!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

				bx(i,j,k)=Bxreg
				by(i,j,k)=Byreg
				bz(i,j,k)=Bzreg

				ex(i,j,k)=Exreg
				ey(i,j,k)=Eyreg
				ez(i,j,k)=Ezreg
			enddo
		enddo
	enddo

	do ki = 0, maxKx
		do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
				do kk = 0, maxKz
#endif
				!print *, ki, kj, kk

			if (((ki + kj + kk) .ne. 0) .and. ((kj + kk).ne.0)) then

				phase1 = 2*pi*rand();
				phase2 = 2*pi*rand();


				kx = ki*2*pi/maxTurbulentLambdaX;
				ky = kj*2*pi/maxTurbulentLambdaY;
				kz = kk*2*pi/maxTurbulentLambdaZ;

				!print *,'k', kx, ky, kz


				kw = sqrt(kx*kx + ky*ky + kz*kz);
				kyz = sqrt(ky*ky + kz*kz);
				cosTheta = kx/kw;
				sinTheta = kyz/kw;
				if(kj + kk .ne. 0) then
					cosPhi = ky/kyz;
					sinPhi = kz/kyz;
				else
					cosPhi = 1.0
					sinPhi = 0.0
				endif

				Bturbulent = evaluate_turbulent_b(kx, ky, kz)*tempB*turbulenceFieldCorrection;


				!print *, 'Bturbulent', Bturbulent
				do  k=1,mz
					do  j=1,my
						do  i=1,mx
							! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
							!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

							kmultr = kx*gamma0*xglob(1.0*i) + ky*yglob(1.0*j) + kz*zglob(1.0*k)
							localB1 = Bturbulent*sin(kmultr + phase1);
							localB2 = Bturbulent*sin(kmultr + phase2);
							!localB2 = 0

							bx(i,j,k)=bx(i,j,k) - localB1*sinTheta
							by(i,j,k)=by(i,j,k) + (localB1*cosTheta*cosPhi - localB2*sinPhi)*gamma0
							bz(i,j,k)=bz(i,j,k) + (localB1*cosTheta*sinPhi + localB2*cosPhi)*gamma0

							ex(i,j,k)=ex(i,j,k)
							ey(i,j,k)=ey(i,j,k) - beta*gamma0*(localB1*cosTheta*sinPhi + localB2*cosPhi)
							ez(i,j,k)=ez(i,j,k) + beta*gamma0*(localB1*cosTheta*cosPhi - localB2*sinPhi)



							!bx(i,j,k)=bx(i,j,k) - localB1*cosTheta*cosPhi + localB2*sinPhi;
							!by(i,j,k)=by(i,j,k) - localB1*cosTheta*sinPhi - localB2*cosPhi;
							!bz(i,j,k)=bz(i,j,k) + localB1*sinTheta;

							!ex(i,j,k)=ex(i,j,k);
							!ey(i,j,k)=ey(i,j,k) - beta*(localB1*sinTheta);
							!ez(i,j,k)=ez(i,j,k) + beta*(- localB1*cosTheta*sinPhi - localB2*cosPhi);
						enddo
					enddo
				enddo
			endif
#ifndef twoD
			enddo
#endif
		enddo
	enddo
end subroutine init_turbulent_field_simple_anisotropic_plasma_frame

subroutine init_turbulent_field_simple_anisotropic_lab_frame(turbulenceEnergyFraction)
			real turbulenceEnergyFraction

			integer maxKx, maxKy, maxKz
			integer :: i, j, k, ki, kj, kk
			real kw
			real kx, ky, kz, kyz
			real phase1, phase2
			real cosTheta, sinTheta, cosPhi, sinPhi
			real Bturbulent
			real kmultr
			real localB1, localB2
			real turbulenceEnergy
			real slabFraction
			real slabEnergy
			real restEnergy

			real tempB
			real Bxreg, Byreg, Bzreg, Exreg, Eyreg, Ezreg
			real Btheta0;
			real pi;

			pi = 2.0*acos(0.0);

			tempB = 1.0;


			maxKx = maxTurbulentLambdaX/minTurbulentLambdaX;
			maxKy = maxTurbulentLambdaY/minTurbulentLambdaY;
			maxKz = maxTurbulentLambdaZ/minTurbulentLambdaZ;

			!maxKy = 0

			print *, maxKx, maxKy, maxKz

			!turulence energy fraction correction in plasma frame

			turbulenceEnergy = 0;

#ifdef twoD
	maxKz = 1;
#endif


			do ki = 0, maxKx

				do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
			!do kk = 0, mz0-5
				do kk = 0, maxKz
#endif

					if (((ki + kj + kk) .ne. 0) .and. ((kj + kk).ne.0)) then


						kx = ki*2*pi/maxTurbulentLambdaX;
						ky = kj*2*pi/maxTurbulentLambdaY;
						kz = kk*2*pi/maxTurbulentLambdaZ;


						Bturbulent = evaluate_turbulent_b(kx, ky, kz);

						turbulenceEnergy = turbulenceEnergy + Bturbulent*Bturbulent;
					endif
#ifndef twoD
			enddo
#endif
		enddo
			enddo

			if (turbulenceEnergy > 0) then
				turbulenceFieldCorrection = sqrt(turbulenceEnergyFraction*tempB*tempB/((1.0 - turbulenceEnergyFraction)*turbulenceEnergy));
			else
				turbulenceFieldCorrection = 1.0;
			endif

			!turbulence total energy correction in lab frame

			turbulenceEnergy = tempB*tempB/(1.0 - turbulenceEnergyFraction);

			if (turbulenceEnergy > 0) then
				tempB = Binit/sqrt(turbulenceEnergy);
			else
				tempB = Binit;
			endif

			Binit = tempB

			!init fields

			Bxreg = tempB*cos(Btheta);
			Byreg = tempB*sin(Btheta)*sin(Bphi);
			Bzreg = tempB*sin(Btheta)*cos(Bphi);

			Exreg = 0;
			Eyreg = -beta*Bzreg;
			Ezreg = beta*Byreg;

			do  k=1,mz
				do  j=1,my
					do  i=1,mx
						! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
						!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

						bx(i,j,k)=Bxreg
						by(i,j,k)=Byreg
						bz(i,j,k)=Bzreg

						ex(i,j,k)=Exreg
						ey(i,j,k)=Eyreg
						ez(i,j,k)=Ezreg
					enddo
				enddo
			enddo

			do ki = 0, maxKx
				do kj = 0, maxKy
#ifdef twoD
			kk = 0
#else
				do kk = 0, maxKz
#endif
				!print *, ki, kj, kk

					if (((ki + kj + kk) .ne. 0) .and. ((kj + kk).ne.0)) then

						phase1 = 2*pi*rand();
						phase2 = 2*pi*rand();


						kx = ki*2*pi/maxTurbulentLambdaX;
						ky = kj*2*pi/maxTurbulentLambdaY;
						kz = kk*2*pi/maxTurbulentLambdaZ;

						!print *,'k', kx, ky, kz


						kw = sqrt(kx*kx + ky*ky + kz*kz);
						kyz = sqrt(ky*ky + kz*kz);
						cosTheta = kx/kw;
						sinTheta = kyz/kw;
						if(kj + kk .ne. 0) then
							cosPhi = ky/kyz;
							sinPhi = kz/kyz;
						else
							cosPhi = 1.0
							sinPhi = 0.0
						endif

						Bturbulent = evaluate_turbulent_b(kx, ky, kz)*tempB*turbulenceFieldCorrection;


						!print *, 'Bturbulent', Bturbulent
						do  k=1,mz
							do  j=1,my
								do  i=1,mx
									! can have fields depend on xglob(i), yglob(j), zglob(j) or iglob(i), jglob(j), kglob(k)
									!print *, xglob(1.0*i), yglob(1.0*j),zglob(1.0*k)

									kmultr = kx*xglob(1.0*i) + ky*yglob(1.0*j) + kz*zglob(1.0*k)
									localB1 = Bturbulent*sin(kmultr + phase1);
									localB2 = Bturbulent*sin(kmultr + phase2);
									!localB2 = 0

									bx(i,j,k)=bx(i,j,k) - localB1*sinTheta
									by(i,j,k)=by(i,j,k) + (localB1*cosTheta*cosPhi - localB2*sinPhi)
									bz(i,j,k)=bz(i,j,k) + (localB1*cosTheta*sinPhi + localB2*cosPhi)

									ex(i,j,k)=ex(i,j,k)
									ey(i,j,k)=ey(i,j,k) - beta*(localB1*cosTheta*sinPhi + localB2*cosPhi)
									ez(i,j,k)=ez(i,j,k) + beta*(localB1*cosTheta*cosPhi - localB2*sinPhi)



									!bx(i,j,k)=bx(i,j,k) - localB1*cosTheta*cosPhi + localB2*sinPhi;
									!by(i,j,k)=by(i,j,k) - localB1*cosTheta*sinPhi - localB2*cosPhi;
									!bz(i,j,k)=bz(i,j,k) + localB1*sinTheta;

									!ex(i,j,k)=ex(i,j,k);
									!ey(i,j,k)=ey(i,j,k) - beta*(localB1*sinTheta);
									!ez(i,j,k)=ez(i,j,k) + beta*(- localB1*cosTheta*sinPhi - localB2*cosPhi);
								enddo
							enddo
						enddo
					endif
#ifndef twoD
			enddo
#endif
		enddo
			enddo
end subroutine init_turbulent_field_simple_anisotropic_lab_frame



#ifdef twoD
end module m_user
#else
end module m_user_3d
#endif

