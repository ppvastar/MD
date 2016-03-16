!This is the main file 

!Declaration and initialization of parameters/variables 
module constants
  implicit none
  
  !Four 4x4 matrices $sigma_0\tau_z$,$\sigma_y\tau_z$, $\sigma_z\tau_0$, $\sigma_0\tau_x$ in Eqs.~(3) and (4)
  !in manuscript http://arxiv.org/abs/1601.01402.
  complex*16 s0tz(4,4),sytz(4,4),szt0(4,4),s0tx(4,4)
  
  !Site number of QD chain
  integer, parameter::N=1000
  
  !Superconducting pairing energy, set as the energy unit (equal to 1)
  real*8, parameter::delta=1.d0
  
  !Mean value of spin-conserving hopping
  real*8, parameter::t0=1.d0
  !Mean value of spin-flip hopping 
  real*8, parameter::tso0=0.5d0

  !Fluctuation magnitude of on-site chemical potential (in units of superconducting pairing energy, delta)
  real*8, parameter::delta_onsite=0.d0 !unit: delta 
  !Fluctuation magnitude of spin-conserving hopping (in units of its mean value, t0)
  real*8, parameter::delta_t=0.d0      !unit: t0
  !Fluctuation magnitude of spin-flip hopping (in units of superconducting pairing energy, delta)
  real*8, parameter::delta_tso=0.d0    !unit: delta

  !On-site chemical potentials along the chain with site number N
  real*8 onsite(N)   !unit:delta
  !Spin-conserving hopping energies along the chain with site number N
  real*8 t(N)        !unit:delta
  !Spin-flip hopping energies along the chain with site number N
  real*8 tso(N)      !unit:delta 

  !Discretization of the Zeeman energy in the region: zm_init+[0,zm_D]*zm_step 
  real*8, parameter::zm_init=0.d0
  real*8, parameter::zm_step=0.05d0
  integer, parameter::zm_D=100
  !Discretization of the mean value of on-site chemical potentials in the region: mu_init+[0,mu_D]*mu_step 
  real*8, parameter::mu_init=0.d0 
  real*8, parameter::mu_step=0.05d0
  integer, parameter::mu_D=100

  !Matrix of the varying Zeeman energy 
  real*8 zm(0:zm_D)
  !Matrix of the varying mean value of on-site chemical potentials
  real*8 mu(0:mu_D)
  
  !Topological charge against Zeeman energy and mean value of on-site chemical potentials
  real*8 q(0:zm_D,0:mu_D)
  
 
#ifdef __solve_eigenstates
  !Hamiltonian defined by Eq.~(2) in manuscript http://arxiv.org/abs/1601.01402.
  complex*16 Hami(4*N,4*N)
  !In the calculation of Fig.~2(b), the Zeeman energy is fixed at a value given by zm_constant
  real*8, parameter::zm_constant=2.d0  !unit:delta 
#endif
  
end module constants

!Main process of code 
program main
  use constants
  integer i,j

#ifdef __solve_eigenstates
  !Eigenenergies in the QD chain system
  real*8 eigenvalue(4*N)
#endif
  
  call initialize_constant

  !Randomize the QD chain 
  call random_onsite
  call random_t
  call random_tso


  !Solve the topological charge for different Zeeman energies and average chemical potentials 
#ifdef __solve_topological_charge
 
  do i=0,zm_D
     zm(i)=zm_init+zm_step*i
     do j=0,mu_D
        mu(j)=mu_init+mu_step*j
        call recursive(i,j)
     end do
  end do
 

  open(0,file='result.dat',access='append')
  do i=0,zm_D
     do j=0,mu_D
        write(0,'(3e14.6)') zm(i),mu(j),q(i,j)
     end do
     write(0,*)
  end do
  close(0)
#endif

  
  !Solve the eigen states of the QD chain  
#ifdef __solve_eigenstates
  open(1,file='eigen.dat',access='append')
  open(2,file='q.dat',access='append')
  zm(0)=zm_constant 
  do j=0,mu_D
     mu(j)=mu_init+mu_step*j
     call recursive(0,j)
     call diago(Hami,4*N,eigenvalue)
     write(1,'(7e14.6)') mu(j),eigenvalue(2*N-2),eigenvalue(2*N-1),eigenvalue(2*N),eigenvalue(2*N+1),eigenvalue(2*N+2),eigenvalue(2*N+3)
     write(2,'(2e14.6)') mu(j),q(0,j)
  end do
  close(1)
  close(2)
#endif
  
end program main



!Subroutine to initialize constant matrices
SUBROUTINE initialize_constant
  use constants
  
  s0tz=0.d0
  sytz=0.d0
  szt0=0.d0
  s0tx=0.d0
 
  !sigma_0*tau_z
  s0tz(1,1)=cmplx(1.d0,0.d0)
  s0tz(2,2)=cmplx(1.d0,0.d0)
  s0tz(3,3)=cmplx(-1.d0,0.d0)
  s0tz(4,4)=cmplx(-1.d0,0.d0)

  !sigma_y*tau_z 
  sytz(1,2)=cmplx(0.d0,-1.d0)
  sytz(2,1)=cmplx(0.d0,1.d0)
  sytz(3,4)=cmplx(0.d0,1.d0)
  sytz(4,3)=cmplx(0.d0,-1.d0)

  !sigma_z*tau_0
  szt0(1,1)=cmplx(1.d0,0.d0)
  szt0(2,2)=cmplx(-1.d0,0.d0)
  szt0(3,3)=cmplx(1.d0,0.d0)
  szt0(4,4)=cmplx(-1.d0,0.d0)
  
  !sigma_0*tau_x
  s0tx(1,3)=cmplx(1.d0,0.d0)
  s0tx(2,4)=cmplx(1.d0,0.d0)
  s0tx(3,1)=cmplx(1.d0,0.d0)
  s0tx(4,2)=cmplx(1.d0,0.d0)

end SUBROUTINE initialize_constant


!Subroutine to randomize on-site chemical potentials 
SUBROUTINE random_onsite
  use constants
  integer i
  real*8 x1
  
  call random_seed()
  do i=1,N
     call random_number(x1)
     onsite(i)=delta_onsite*(2.d0*x1-1)
  end do
  
end SUBROUTINE random_onsite

!Subroutine to randomize spin-conserving hopping energies
SUBROUTINE random_t
  use constants
  integer i
  real*8 x1

  call random_seed()
  do i=1,N
     call random_number(x1)
     t(i)=t0+delta_t*(2.d0*x1-1)*t0
  end do
    
end SUBROUTINE random_t

!Subroutine to randomize spin-flip hopping energies
SUBROUTINE random_tso
  use constants
  integer i
  real*8 x1
  
  call random_seed()
  do i=1,N
     call random_number(x1)
     tso(i)=tso0+delta_tso*(2.d0*x1-1.d0)
  end do
  
end SUBROUTINE random_tso

!Subroutine to solve the topological charge when the Zeeman energy is zm(i_tmp) and the average
!chemical potential is mu(j_tmp), following the recursive way given in the appendix of manuscript http://arxiv.org/abs/1601.01402.
SUBROUTINE recursive(i_tmp,j_tmp)
  use constants
  integer i_tmp,j_tmp
  integer i,a,b,j1,j2
  complex*16 h(4,4),tc(4,4),tc_inv(4,4),tc_hc(4,4)
  complex*16 M(8,8),M_tmp(8,8),M11(4,4),M12(4,4),M21(4,4),M22(4,4),M22_inv(4,4)
  complex*16 W(8,8),W_tmp(8,8),W11(4,4),W12(4,4),W21(4,4),W22(4,4)
  complex*16 R(4,4),T_tmp(4,4)
  complex*16 det

#ifdef __solve_eigenstates 
  Hami=0.d0
#endif
     
  do i=1,N
     h=(onsite(i)-mu(j_tmp))*s0tz+zm(i_tmp)*szt0+delta*s0tx
     tc=t(i)*s0tz+cmplx(0.d0,1.d0)*(tso(i)*sytz)

#ifdef __solve_eigenstates           
     do j1=1,4
        do j2=1,4
           Hami(4*(i-1)+j1,4*(i-1)+j2)=h(j1,j2)
           if (i.le.N-1) then            
              Hami(4*(i-1)+j1,4*i+j2)=tc(j1,j2)
           end if 
        end do
     end do 
#endif
     
     tc_inv=tc
     call inverse(tc_inv,4)
     tc_hc=conjg(transpose(tc))
    
     M11=0.5d0*(cmplx(0.d0,1.d0)*(tc_inv+tc_hc)-matmul(tc_inv,h))
     M12=0.5d0*(cmplx(0.d0,1.d0)*(tc_inv-tc_hc)+matmul(tc_inv,h))
     M21=0.5d0*(cmplx(0.d0,1.d0)*(-tc_inv+tc_hc)+matmul(tc_inv,h))
     M22=0.5d0*(cmplx(0.d0,1.d0)*(-tc_inv-tc_hc)-matmul(tc_inv,h))
     
     M22_inv=M22  
     call inverse(M22_inv,4)
     W11=-matmul(M22_inv,M21)
     W12=M22_inv
     W21=M11-matmul(M12,matmul(M22_inv,M21))
     W22=matmul(M12,M22_inv)

     do a=1,4
        do b=1,4 
           M_tmp(a,b)=M11(a,b)
           M_tmp(a,4+b)=M12(a,b)
           M_tmp(4+a,b)=M21(a,b)
           M_tmp(4+a,4+b)=M22(a,b)
           W_tmp(a,b)=W11(a,b)
           W_tmp(a,4+b)=W12(a,b)
           W_tmp(4+a,b)=W21(a,b)
           W_tmp(4+a,4+b)=W22(a,b)
        end do
     end do
     
     if(i.eq.1) then
        M=M_tmp
        W=W_tmp 
     else 
        M=matmul(M_tmp,M)
        call convolution(W_tmp,W)
     end if
  end do

  
  do a=1,4
     do b=1,4
        M22(a,b)=M(4+a,4+b)
        M21(a,b)=M(4+a,b)
        R(a,b)=W(a,b)
        T_tmp(a,b)=W(4+a,b)
     end do
  end do

  if(real(det(R,4))>0.d0) then
     q(i_tmp,j_tmp)=1.d0
  else
     q(i_tmp,j_tmp)=-1.d0
  end if
  
  
#ifdef __solve_eigenstates 
  do j1=1,4*N
     do j2=1,j1-1
        Hami(j1,j2)=conjg(Hami(j2,j1))
     end do
  end do
#endif
  
end SUBROUTINE recursive

!Subroutine to solve the inverse of matrix "A" with dimension "DIM"
SUBROUTINE inverse(A,DIM)
  implicit none
  integer DIM
  complex*16 A(1:DIM,1:DIM)
  integer IPIV(1:DIM)
  integer LWORK
  complex*16 WORK(1:DIM*DIM)
  integer INFO
  LWORK = DIM*DIM
  call zgetrf(DIM,DIM,A,DIM,IPIV,INFO)
  if(INFO.ne.0)then
     print*,"error1",INFO
     stop
  end if
  call zgetri(DIM,A,DIM,IPIV,WORK,LWORK,INFO)
  if(INFO.ne.0)then
     print*,"error2",INFO
     stop
  end if
end SUBROUTINE inverse


!Subroutine to solve the determinant of matrix "A" with dimension "DIM"
complex*16 function det(A,DIM)
  implicit none
  integer DIM
  complex*16 A(1:DIM,1:DIM)
  integer IPIV(1:DIM)
  integer INFO
  real*8 sgn
  integer i 
  call zgetrf(DIM,DIM,A,DIM,IPIV,INFO)
  if(INFO.ne.0)then
     print*,"error3",INFO
     stop
  end if
  
  det=dcmplx(1.d0,0.d0)
  sgn=1.d0
  
  do i=1, DIM
     det=det*A(i,i)
     if(IPIV(i).ne.i) then
        sgn=-sgn
     end if
  end do

  det=sgn*det   
end  function det


!Subroutine to calculate the "convolution" of matrics "A" and "B" defined by Eq.~(A6) in manuscript http://arxiv.org/abs/1601.01402.
SUBROUTINE convolution(A,B)
  use constants
  complex*16 A(8,8),B(8,8)
  complex*16 A11(4,4),A12(4,4),A21(4,4),A22(4,4),B11(4,4),B12(4,4),B21(4,4),B22(4,4)
  complex*16 u(4,4),mat_tmp1(4,4),mat_tmp2(4,4)
  integer i,j
  u=0.d0
  do i=1,4
     u(i,i)=cmplx(1.d0,0.d0)
  end do

  do i=1,4
     do j=1,4
        A11(i,j)=A(i,j)
        A12(i,j)=A(i,4+j)
        A21(i,j)=A(4+i,j)
        A22(i,j)=A(4+i,4+j)
        B11(i,j)=B(i,j)
        B12(i,j)=B(i,4+j)
        B21(i,j)=B(4+i,j)
        B22(i,j)=B(4+i,4+j)
     end do
  end do
         
  mat_tmp1=u-matmul(A11,B22)
  mat_tmp2=u-matmul(B22,A11)
  call inverse(mat_tmp1,4)
  call inverse(mat_tmp2,4)

  B11=B11+matmul(B12,matmul(mat_tmp1,matmul(A11,B21)))
  B12=matmul(B12,matmul(mat_tmp1,A12))
  B21=matmul(A21,matmul(mat_tmp2,B21))
  B22=A22+matmul(A21,matmul(mat_tmp2,matmul(B22,A12)))

  do i=1,4
     do j=1,4
        B(i,j)=B11(i,j)
        B(i,4+j)=B12(i,j)
        B(4+i,j)=B21(i,j)
        B(4+i,4+j)=B22(i,j)
     end do
  end do
end SUBROUTINE convolution


!Subroutine to diagonalize a Hermitian matrix "Ham" with dimension "order"
SUBROUTINE diago(Ham, order, Eigenvalue)
  implicit none
  integer order
  complex*16 Ham(1:order,1:order)
  real*8 Eigenvalue(1:order)
  complex*16 work(1:2*order)
  integer info,lwork
  real*8 rwork(1:3*order-2)
  
  lwork=2*order
  
  call zheev('V','U',order,Ham,order,Eigenvalue,work,lwork,rwork,info)      
  if(info.ne.0)then
     write(*,*)'zheev process is wrong'
     stop
  endif
  
end SUBROUTINE diago
