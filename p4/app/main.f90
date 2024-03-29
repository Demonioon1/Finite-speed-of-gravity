PROGRAM ver
  use parametros
  USE constants_NSWC
  IMPLICIT NONE
 
  
  REAL (dp)    :: coef(0:4)
  COMPLEX (dp) :: root(4)
  
  INTERFACE
    SUBROUTINE qdcrt (coef, root)
      USE constants_NSWC
      IMPLICIT NONE
      REAL (dp), INTENT(IN)     :: coef(:)
      COMPLEX (dp), INTENT(OUT) :: root(:)
    END SUBROUTINE qdcrt
  
    SUBROUTINE cbcrt (coef, root)
      USE constants_NSWC
      IMPLICIT NONE
      REAL (dp), INTENT(IN)     :: coef(:)
      COMPLEX (dp), INTENT(OUT) :: root(:)
    END SUBROUTINE cbcrt
  
    SUBROUTINE qtcrt (coef, root)
      USE constants_NSWC
      IMPLICIT NONE
      REAL (dp), INTENT(IN)     :: coef(:)
      COMPLEX (dp), INTENT(OUT) :: root(:)
    END SUBROUTINE qtcrt
  END INTERFACE


  call cpu_time(start)

  read(*,'(a)') datin
  dat = len_trim(datin)
r=0.0q0
v=0.0q0
a=0.0q0
root=0.0q0

  ! m(1)=333000
  ! m(2)=1
  ! m(3)=0.0123

    m(1)=1.0q0
    m(2)=0.0123q0
  

! if(N.ge.3)then
!   do i=1,N
!     call random_seed()
!         call random_number(az)
!     m(i)=az
!     !write(*,*) m(i)/10
!   end do
! end if
!  m(1)=10000
! v(:,1,:)=0.d0
!  m(11)=10000
! v(:,11,:)=0.d0

    ! if(N.ge.3)then
    !   do i=1,N
    !     m(i)=1
    !   end do
    ! end if
    ! m(1)=333000

!vp=999999999
!vp=10000000
!vp=10000
!vp=100
!vp=10
vp=1
!vp=0.1q0
!vp=0.0017q0
!vp=0.001q0
!vp=0.0005
!vp=0.00032
!vp=0.0003
!vp=0.0001
!vp=0.000001
!vp=0.0000001
!vp=0.00000001
!vp=0.000000001
!vp=0.0000000001
!vp=0.00000000001
!vp=0.000000000001
!vp=0.0000000000001
!vp=0.0000000000000000001

vp=vp*63242

 
  open(unit=1,file=datin(1:dat),status='unknown')
     do l=1,N
      read (1,*) r(:,l,1),v(:,l,1),a(:,l,1)
     end do


      a=0
      aux=0
      step=1

    open(file='ener.dat',unit=77,status='unknown')
    open(file='dat.dat',unit=7,status='unknown')
    open(file='orb.dat',unit=9,status='unknown')
    open(file='coef.dat',unit=13,status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i = 1,N-1
         do j = i+1,N
            Rij = r(:,i,1) - r(:,j,1)
            Rsq= dot_product(Rij ,Rij)

           phij  = -(m(j)*CG*m(i))/sqrt(Rsq)
           dphij = -(m(j)*CG)/((sqrt(Rsq))**3)

           phji  = -(m(j)*CG*m(i))/sqrt(Rsq)
           dphji = -(CG*m(i))/((sqrt(Rsq))**3)

           p(i,1) = p(i,1) + 0.5q0*phij
           p(j,1) = p(j,1) + 0.5q0*phji

           a(:,i,1) = a(:,i,1) + dphij*Rij
           a(:,j,1) = a(:,j,1) - dphji*Rij


          enddo
        enddo



    do  i = 1,N-1
      do j = i+1,N
  
       AC=0.25*dot_product(a(:,j,1),a(:,j,1))
       BC=dot_product(v(:,j,1),a(:,j,1))
       CC=dot_product(v(:,j,1),v(:,j,1))-dot_product(r(:,i,1)-rret1(:,j,1),a(:,j,1))
       DC=-2*dot_product(r(:,i,1)-rret1(:,j,1),v(:,j,1))
       EC= dot_product(r(:,i,1)-rret1(:,j,1),r(:,i,1)-rret1(:,j,1))

       coef(0)= AC*dt**4 + BC*dt**3 + CC*dt**2 + DC*dt + EC
       coef(1)= -(4*AC*dt**3 +3*BC*dt**2 + 2*CC*dt + DC)
       coef(2)= 6*AC*dt**2 + 3*BC*dt + CC - vp**2
       coef(3)=-(4*AC*dt + BC)
       coef(4)=AC
  
        coef(4)=0
        CALL cbcrt(coef, root)
              
        do l=1,4
          if (aimag(root(l)).ne.0) then
            root(l)=0
          endif
          if (real(root(l)).lt.0) then
            root(l)=0
          endif
        enddo
          ret1=max(real(root(1)),real(root(2)),real(root(3)),real(root(4)))
        

        do l=1,4
          if (real(root(l)).ne.0) then
            if(real(root(l)).lt.ret1) then
            ret1=real(root(l))
            end if
          endif
        enddo


        AC=0.25*dot_product(a(:,i,1),a(:,i,1))
        BC=dot_product(v(:,i,1),a(:,i,1))
        CC=dot_product(v(:,i,1),v(:,i,1))-dot_product(r(:,j,1)-rret1(:,i,1),a(:,i,1))
        DC=-2*dot_product(r(:,j,1)-rret1(:,i,1),v(:,i,1))
        EC= dot_product(r(:,j,1)-rret1(:,i,1),r(:,j,1)-rret1(:,i,1))

        coef(0)= AC*dt**4 + BC*dt**3 + CC*dt**2 + DC*dt + EC
        coef(1)= -(4*AC*dt**3 +3*BC*dt**2 + 2*CC*dt + DC)
        coef(2)= 6*AC*dt**2 + 3*BC*dt + CC - vp**2
        coef(3)=-(4*AC*dt + BC)
        coef(4)=AC              

       
  
        coef(4)=0  
        CALL cbcrt(coef, root)

       
        
        
        do l=1,4
          if (aimag(root(l)).ne.0.0) then
            root(l)=0
          endif
          if (real(root(l)).lt.0) then
            root(l)=0
          endif
        enddo
          ret2=max(real(root(1)),real(root(2)),real(root(3)),real(root(4)))
        

        do l=1,4
         if (real(root(l)).ne.0) then
            if(real(root(l)).lt.ret2) then
            ret2=real(root(l))
            end if
          endif
        enddo
     
        

        rret(:,i,1) = r(:,i,1) + (dt-ret1)*v(:,i,1) + 0.5d0*((dt-ret1)**2)*a(:,i,1)
        rret(:,j,1) = r(:,j,1) + (dt-ret2)*v(:,j,1) + 0.5d0*((dt-ret2)**2)*a(:,j,1)

        Rij = r(:,i,1) - rret(:,j,1)
        Rji = r(:,j,1) - rret(:,i,1)

        Rsqij= dot_product(Rij ,Rij)
        Rsqji= dot_product(Rji ,Rji)

        phij  = -(m(j)*CG*m(i))/sqrt(Rsqij)
        dphij = -(m(j)*CG)/((sqrt(Rsqij))**3)

        phji  = -(m(j)*CG*m(i))/sqrt(Rsqji)
        dphji = -(CG*m(i))/((sqrt(Rsqji))**3)

        p(i,1) = p(i,1) + 0.5q0*phij
        p(j,1) = p(j,1) + 0.5q0*phji

        a(:,i,1) = a(:,i,1) + dphij*Rij
        a(:,j,1) = a(:,j,1) + dphji*Rji
  
      enddo
    enddo






















    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        do step=1,nt
    
          rret1=r
          v(:,:,1) = v(:,:,1) + 0.5q0*dt*a(:,:,1)
          r(:,:,1) = r(:,:,1) + dt*v(:,:,1) 
          !r(:,:,1) = r(:,:,1) + dt*v(:,:,1) + 0.5d0*(dt**2)*a(:,:,1)
            a(:,:,2)=0
            eka=0
            epa=0
            p=0
            eta=0

            do  i = 1,N-1
              do j = i+1,N
          
               AC=0.25*dot_product(a(:,j,1),a(:,j,1))
               BC=dot_product(v(:,j,1),a(:,j,1))
               CC=dot_product(v(:,j,1),v(:,j,1))-dot_product(r(:,i,1)-rret1(:,j,1),a(:,j,1))
               DC=-2*dot_product(r(:,i,1)-rret1(:,j,1),v(:,j,1))
               EC= dot_product(r(:,i,1)-rret1(:,j,1),r(:,i,1)-rret1(:,j,1))

               coef(0)= AC*dt**4 + BC*dt**3 + CC*dt**2 + DC*dt + EC
               coef(1)= -(4*AC*dt**3 +3*BC*dt**2 + 2*CC*dt + DC)
               coef(2)= 6*AC*dt**2 + 3*BC*dt + CC - vp**2
               coef(3)=-(4*AC*dt + BC)
               coef(4)=AC
          
                coef(4)=0
                CALL cbcrt(coef, root)
                      
                do l=1,4
                  if (aimag(root(l)).ne.0) then
                    root(l)=0
                  endif
                  if (real(root(l)).lt.0) then
                    root(l)=0
                  endif
                enddo
                  ret1=max(real(root(1)),real(root(2)),real(root(3)),real(root(4)))
                

                do l=1,4
                  if (real(root(l)).ne.0) then
                    if(real(root(l)).lt.ret1) then
                    ret1=real(root(l))
                    end if
                  endif
                enddo


                AC=0.25*dot_product(a(:,i,1),a(:,i,1))
                BC=dot_product(v(:,i,1),a(:,i,1))
                CC=dot_product(v(:,i,1),v(:,i,1))-dot_product(r(:,j,1)-rret1(:,i,1),a(:,i,1))
                DC=-2*dot_product(r(:,j,1)-rret1(:,i,1),v(:,i,1))
                EC= dot_product(r(:,j,1)-rret1(:,i,1),r(:,j,1)-rret1(:,i,1))

                coef(0)= AC*dt**4 + BC*dt**3 + CC*dt**2 + DC*dt + EC
                coef(1)= -(4*AC*dt**3 +3*BC*dt**2 + 2*CC*dt + DC)
                coef(2)= 6*AC*dt**2 + 3*BC*dt + CC - vp**2
                coef(3)=-(4*AC*dt + BC)
                coef(4)=AC              

               
          
                coef(4)=0  
                CALL cbcrt(coef, root)

               
                
                
                do l=1,4
                  if (aimag(root(l)).ne.0.0) then
                    root(l)=0
                  endif
                  if (real(root(l)).lt.0) then
                    root(l)=0
                  endif
                enddo
                  ret2=max(real(root(1)),real(root(2)),real(root(3)),real(root(4)))
                

                do l=1,4
                 if (real(root(l)).ne.0) then
                    if(real(root(l)).lt.ret2) then
                    ret2=real(root(l))
                    end if
                  endif
                enddo
             
                

                rret(:,i,1) = r(:,i,1) + (dt-ret1)*v(:,i,1) + 0.5d0*((dt-ret1)**2)*a(:,i,1)
                rret(:,j,1) = r(:,j,1) + (dt-ret2)*v(:,j,1) + 0.5d0*((dt-ret2)**2)*a(:,j,1)

                Rij = r(:,i,1) - rret(:,j,1)
                Rji = r(:,j,1) - rret(:,i,1)

                Rsqij= dot_product(Rij ,Rij)
                Rsqji= dot_product(Rji ,Rji)

                phij  = -(m(j)*CG*m(i))/sqrt(Rsqij)
                dphij = -(m(j)*CG)/((sqrt(Rsqij))**3)

                phji  = -(m(j)*CG*m(i))/sqrt(Rsqji)
                dphji = -(CG*m(i))/((sqrt(Rsqji))**3)

                p(i,1) = p(i,1) + 0.5q0*phij
                p(j,1) = p(j,1) + 0.5q0*phji

                a(:,i,2) = a(:,i,2) + dphij*Rij
                a(:,j,2) = a(:,j,2) + dphji*Rji
          
              enddo
            enddo

            v(:,:,1) = v(:,:,1) + 0.5q0*dt*(a(:,:,2))
            !v(:,:,1) = v(:,:,1) + 0.5d0*dt*(a(:,:,2)+a(:,:,1))
            a(:,:,1)=a(:,:,2)

            do i=1,N
              real_vel = v(:,i,1)
              eka = eka+0.5q0*m(i)*dot_product(real_vel,real_vel)
              epa=epa+p(i,1)
            enddo
  
            eta = eka+epa
          
            aux=aux+1

          if(aux.eq.100000)then


            do nn=1,N
              write(7,*) r(1,nn,1),r(2,nn,1),r(3,nn,1)
            end do
            write(7,*) ' '
            write(7,*) ' '


            write(9,*) r(:,2,1)-r(:,1,1)!,r(:,3,1)-r(:,2,1)
            write(77,*) step*dt, eka,epa,eta
            aux=0
          endif

          
          
        enddo
    close(13)
    close (9)
    close (7)
    close (77)
  close(1)
  

  call cpu_time(finish)
  write(*,*) finish-start
    !write(*,*)  r

  END PROGRAM ver