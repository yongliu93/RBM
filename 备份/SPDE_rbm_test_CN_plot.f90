     
      module comdata
      implicit none
      
      type vector
          integer , allocatable, dimension (:) :: index
      end type vector
      
      real, allocatable, dimension (:,:) :: u
      integer, allocatable, dimension(:,:) :: alpha,Ic
      type (vector), allocatable, dimension(:) :: s_eps,s_alpha
      
      type(vector), allocatable, dimension(:) :: I_cindex
      
      real, allocatable, dimension(:,:) ::u_num,A,B
      real, allocatable, dimension(:) ::x,D1,D2
      !real, allocatable, dimension (:,:) :: c
      
      real :: tn,endt,dt,pi
      
      integer :: NPC,K,M,numstep,Np
      
      
      real, allocatable,dimension(:) :: f_LHS,unum
      real, allocatable,dimension(:,:) :: r_Basis, C_basis, u_rbm, u_ori
      real, allocatable,dimension(:) :: u0
      
      real :: epsilon
      integer :: num_M
      integer , allocatable, dimension(:) :: basisnum,pickup
      !integer, allocatable, dimension(:):: sq
      real, allocatable,dimension(:,:) :: uAAu,uA,uAAu_i,AAu
      real,allocatable,dimension(:,:) :: AA
      
      real, allocatable,dimension(:) :: frb
      
      real , allocatable, dimension(:,:) :: Mk
      real, allocatable, dimension(:,:,:) :: Mu,uAMu!,uMAu
      
      real, allocatable,dimension(:,:) :: res_Q,res_R!,res_RTR
      
      integer :: maxbasisnum,nrank,ncolumn
      
      integer , allocatable,dimension(:):: ipiv
      
      integer :: ns,nend,num_Mup
      real :: rmineigen,rq1,error2
      
      real, allocatable, dimension(:,:) :: identity, idtA1t,idtA1,idtA2
      real, allocatable, dimension(:,:) :: idtA,idtAi,idtAiA,idtA1tA,idtAtA1
      
      real, allocatable, dimension(:) :: bu0
      
      integer :: numstep_np,numstep_np1
      
    contains
    
    recursive subroutine buildI(ibeg,iend,i_M, I_cindex)
      implicit none
      type(vector), allocatable, dimension (:) :: I_cindex,tmp
      
      integer :: ibeg,iend,i_M
      
      integer :: I_size,nchoosek,nsize
      
      integer :: i,j,length_tmp
      integer, dimension(0:i_M-1) :: itmp
      
      I_size=nchoosek(iend-ibeg+i_M,i_M)
      
      allocate(I_cindex(0:I_size-1))
      do i=0,I_size-1
      allocate(I_cindex(i)%index(0:i_M-1))
      end do
      
      
      
      if (i_M==1) then
          do i=ibeg,iend
              I_cindex(i-ibeg)%index(0)=i
          end do
      else
          nsize=0
          do i=ibeg, iend
              call buildI(i,iend,i_M-1,tmp)
              length_tmp=nchoosek(iend-i+i_M-1,i_M-1)
              do j=0, length_tmp-1
                  itmp(0:i_M-2)=tmp(j)%index(0:i_M-2)
                  itmp(i_M-1)=i
                  I_cindex(nsize+j)%index(:)=itmp(:)
              end do
              nsize=nsize+length_tmp
              deallocate(tmp)
          end do
          
          
      end if
      
          
      
      
      
      
      
    end subroutine buildI
    
    
    
      
    end module comdata
      
    
    !module try
    !implicit none
    !type vector
    !      integer , allocatable, dimension (:) :: index
    !end type vector
    !contains
    !recursive subroutine buildI(ibeg,iend,i_M, I_cindex)
    !  implicit none
    !  type(vector), allocatable, dimension (:) :: I_cindex,tmp
    !  
    !  integer :: ibeg,iend,i_M
    !  
    !  integer :: I_size,nchoosek
    !  
    !  integer :: i,j,length_tmp
    !  integer, dimension(0:i_M-1) :: itmp
    !  
    !  I_size=nchoosek(iend-ibeg+1,i_M)
    !  
    !  allocate(I_cindex(0:I_size-1))
    !  do i=0,I_size-1
    !  allocate(I_cindex(i)%index(0:i_M-1))
    !  end do
    !  
    !  
    !  if (i_M==1) then
    !      do i=ibeg,iend
    !          I_cindex(i-ibeg)%index(0)=i
    !      end do
    !  else
    !      do i=ibeg, iend
    !          call buildI(i,iend,i_M-1,tmp)
    !          length_tmp=nchoosek(iend-i+1,i_M-1)
    !          do j=0, length_tmp-1
    !              itmp(0:i_M-2)=tmp(j)%index(0:i_M-2)
    !              itmp(i_M-1)=i
    !              I_cindex((i-ibeg)*length_tmp+j)%index(:)=itmp(:)
    !          end do
    !      end do
    !      
    !      
    !  end if
    !  
    !      
    !  
    !  
    !  
    !  
    !  
    !  end subroutine buildI
    !
    !end module try
    
    
    
      
      
      program main
      use comdata
      
      numstep=5000
      
      call assignmemory
      call init
      
      
      !endt=5.0
      
      call CPU_TIME(t1)
      call run(numstep)
      call CPU_TIME(t2)
      
      !call comperror
      !pause
      write(*,*) 'time no rbm=',t2-t1
      write(2018,*) 'time no rbm=',t2-t1
      timenorbm=t2-t1
      !call comperror
      !do K=6,6
            !do M=K,K
                
        do i=-1,3
           !if(i<10) then
           !epsilon=100.0-10.*i
           !else
           !    epsilon=1.0-(i-10)/100.
           !end if!1.0E-04*(1-i/2.)
           do j=9,1,-1
           epsilon=10.**(-i)*j!1.0E-02-i*1.0E-03
        
        ! number of polynomials order and random variables 
        !K=3
        !M=3
        
        write(*,*) "(N,K)=",K,M
        write(*,*) "eps_tol=",epsilon
        write(2018,*) "(N,K)=",K,M
        write(2018,*) "eps_tol=",epsilon
        !maxbasisnum=
           
           
      
      
      !endt=1.0
      !
      !call CPU_TIME(t1)
      !call run(numstep)
      !call CPU_TIME(t2)
      !
      !write(*,*) 'time no rbm=',t2-t1
      !write(2018,*) 'time no rbm=',t2-t1
      
      !open(12,file='fort.12',status='old')
      !pickup=1
      !do j=1,15
      !    read(12,*) pos
      !    pickup(pos)=0
      !end do
      !close(12)
      
      
      call CPU_TIME(t1)
      call reduce_basis
      call CPU_TIME(t2)
      
      !call init_rbm
      
      !call assemble_mat
      write(*,*) 'time with rbm=',t2-t1
      write(*,*) 'basis num/ total num',num_M,'/ ', NPC
      
      
      
      write(2018,*) 'time with rbm=',t2-t1
      write(2018,*) 'basis num/ total num',num_M,'/ ', NPC
      
      
      !call run_ori(0)
      !do i=1,NPC-1
      !    call assemble_lhs(i)
      !    
      !    !call argminmize_coefficient(pos)
      !    call run_ori(i)
      !end do
      
      call comperror
      !open(2019,file="fort.2019",status='old')
      !write(*,"(5ES12.2)") -epsilon,error2,t2-t1,real(num_M)
      write(2019,"(5ES12.2)") epsilon,error2,t2-t1,real(num_M),timenorbm
      !call outpdata
      
      
      !pause
            !end do
        !end do
        
           end do
           
       end do
       call deletememory
      
    end
      
      
     
    
    
        subroutine assignmemory
        use comdata
        integer :: nchoosek
        
        K=4
        M=4
        
        Np=32
        
        NPC=nchoosek(K+M,K)
        maxbasisnum=500
        !maxbasisnum=NPC
        numstep_np=numstep*np
        numstep_np1=(numstep+1)*np-1

        
        allocate(Ic(0:NPC-1,0:K-1),alpha(0:NPC-1,0:M-1))
        allocate(s_eps(0:NPC-1),s_alpha(0:NPC-1))
        allocate(u_num(0:NPC-1,0:Np-1))
        
        !allocate(AA(1:numstep_np,1:numstep_np))
        
        allocate(D1(0:np-1),D2(0:np-1),x(0:np-1))
        allocate(A(0:np-1,0:np-1),B(0:np-1,0:np-1))
        
        allocate(identity(0:np-1,0:np-1),idtA1(0:np-1,0:np-1),idtA1t(0:np-1,0:np-1))
        
        allocate(idtA2(0:np-1,0:np-1),idtA(0:np-1,0:np-1),idtAi(0:np-1,0:np-1))
        allocate(idtAiA(0:np-1,0:np-1),idtA1tA(0:np-1,0:np-1),idtAtA1(0:np-1,0:np-1))
        allocate(bu0(0:np-1))
        
        
    do i=0,NPC-1
    
        allocate(s_eps(i)%index(0:0),s_alpha(i)%index(0:0))
        s_eps(i)%index(0)=0
        s_alpha(i)%index(0)=0
    end do
    Ic=0
    alpha=0
    u_num=0.
    !u1=0.
    
    
     allocate(C_basis(0:NPC-1,0:maxbasisnum-1))
     allocate(r_Basis(1:numstep_np,0:maxbasisnum-1))
     allocate(f_LHS(0:numstep_np-1),unum(0:numstep_np1+1))
     
     allocate(u_rbm(0:NPC-1,0:np-1),u0(0:np-1))
     
     allocate(u_ori(0:NPC-1,1:numstep_np))
     
     allocate(basisnum(1:maxbasisnum))
     allocate(pickup(0:NPC-1))
     
     allocate(uAAu(1:maxbasisnum,1:maxbasisnum),uA(1:maxbasisnum,1:numstep_np))
     allocate(uAMu(0:M-1,1:maxbasisnum,1:maxbasisnum),Mu(0:M-1,1:numstep_np,1:maxbasisnum))
     allocate(uAAu_i(1:maxbasisnum,1:maxbasisnum))
     
     allocate(frb(1:maxbasisnum))
     
     allocate(Mk(0:numstep,0:M-1))
     
     allocate(ipiv(1:maxbasisnum))
     
     allocate(res_Q(1:numstep_np,1:maxbasisnum*(M+1)),res_R(1:maxbasisnum*(M+1),1:maxbasisnum*(M+1)))
     
     !allocate(res_RTR(1:maxbasisnum*(M+1),1:maxbasisnum*(M+1)))
     
     end subroutine
      
      !subroutine readdata
      !use comdata
      !integer :: i,j,nsize
      !!character(len=20) :: buffer
      !open(unit=11,file='data.dat')
      !
      !!NPC=70
      !read (11,*) NPC
      !
      !K=4
      !allocate(u(0:NPC-1,1:10000),Ic(0:NPC-1,0:K-1),alpha(0:NPC-1,0:K-1))
      !allocate(s_eps(0:NPC-1),s_alpha(0:NPC-1))
      !
      !allocate(u_num(0:NPC-1),u1(0:NPC-1,1:numstep))
      
      !open(unit=11,file='data.dat')
      !
      !do i=0,NPC-1
      !    do j=1,10000
      !        read(11,*) u(i,j)
      !    end do
      !    !read(11,"(A20)") buffer
      !end do
      !
      !do i=0, NPC-1
      !    do j=0,K-1
      !        read(11,*) Ic(i,j)
      !    end do
      !end do
      !do i=0, NPC-1
      !    do j=0,K-1
      !        read(11,*) alpha(i,j)
      !    end do
      !end do
      !
      !
      !do i=0, NPC-1
      !    read(11,*) nsize
      !    allocate(s_eps(i)%index(0:nsize))
      !    s_eps(i)%index(0)=nsize
      !    do j=1, nsize
      !        read(11,*) s_eps(i)%index(j)
      !    end do
      !end do
      !
      !
      !do i=0, NPC-1
      !    read(11,*) nsize
      !    allocate(s_alpha(i)%index(0:nsize))
      !    s_alpha(i)%index(0)=nsize
      !    do j=1, nsize
      !        read(11,*) s_alpha(i)%index(j)
      !    end do
      !end do
      
      
      
      
      
      
      
      !end subroutine readdata
      
      
!      subroutine build_sequence
!      use comdata
!      integer, dimension(0:M,0:M-1):: eps
!      integer, dimension(0:M-1) :: alpha_tmp
!      
!      integer :: pos
!      
!      allocate(sq(0:K))
!      
!      eps=0
!      do i=1,M
!          eps(i,i-1)=1
!      end do
!      
!      alpha_tmp(:)=eps(0,:)
!      
!      call alphaToPos(alpha_tmp,M,K,pos)
!      
!      allocate(sq(0)%index(0:1))
!      
!      sq(0)%index(0)=1
!      sq(0)%index(1)=pos
!      
!      
!      do i=1,K
!          
!          
!          
!      end do
!      
!          
!          
!      
!      
!      
!      
!      end subroutine
      

    subroutine TriPoly(y,left,right,k,p)
    implicit none
    real :: y, left, right
    integer :: k
    real, dimension(0:k) :: p
    real :: x,pi
    integer :: i
    
    pi=4.*atan(1.0)
    
    x=pi*(y-left)/(right-left)
    p=0.
    
    p(0)=sqrt(1.0/(right-left))
    
    do i=1,k
        p(i)=sqrt(2.0/(right-left))*cos(i*x)
    end do
    
    
    
    
    end subroutine
      
      subroutine LegendrePoly(y,left,right,k,p)
      implicit none
      
      real :: y, left,right
      integer :: k
      real , dimension (0:k) :: p
      
      real :: x
      integer :: i
      p=0.
      
      x=2.0*(y-left)/(right-left)-1.0
      
      p(0)=1.0
      
      if (k>=1) then
          p(1)=x
      end if
      
          do i=2,k
              p(i)=((2.0*i-1.0)*x*p(i-1)-(i-1.0)*p(i-2))/i
          end do
          
          do i=0,k
              p(i)=p(i)*sqrt((2.0*i+1.0)/(right-left))
          end do
          
      

      end subroutine LegendrePoly
      
      
      
      subroutine getDiscretization(res)
      use comdata
      
      real, dimension (0:M-1) :: basis
      integer :: i,j,na
      
      real, intent(out),dimension (0:NPC-1,0:Np-1) ::res
      real, allocatable, dimension (:,:) :: tmp
      
      allocate(tmp(0:NPC-1,0:Np-1))
      
      !res=0.
      !tmp=0.
      
      call LegendrePoly(tn,0.,endt,M-1,basis)
      !call TriPoly(tn,0.,endt,M-1,basis)
      do na=0, NPC-1
          do i=0, NP-1
              res(na,i)=sum(A(i,0:np-1)*u_num(na,0:np-1))
          end do
      end do
      
      do na=0,NPC-1
          do i=0,Np-1
              tmp(na,i)=sum(B(i,0:np-1)*u_num(na,0:np-1))
          end do
      end do
      
      
      
      do na=1,NPC-1
          do i=0,Np-1
          do l=1,s_eps(na)%index(0)
          res(na,i)=res(na,i)+sqrt(real(alpha(na,s_eps(na)%index(l))))* &
          tmp(s_alpha(na)%index(l),i)*basis(s_eps(na)%index(l))
          end do
          end do
          
          
      end do
      
      
      deallocate(tmp)
      
      
      end subroutine getDiscretization
      
      
      subroutine init
      use comdata
      integer :: nchoosek,i,j,pos
      !integer ,allocatable,dimension(:) :: alpha_tmp
      real, dimension(0:M-1) :: basis
      real work(1:4*np)
      INTEGER info,IWORK(1:np)
      real, allocatable, dimension(:,:) :: tmp
      
      allocate(tmp(0:np-1,0:np-1))
      
      u_num=0.
      pi=4.*atan(1.0)
      do i=0,Np-1
          x(i)=2.0*pi/Np*i
          u_num(0,i)=cos(x(i))
          u0(i)=u_num(0,i)
      end do
      do i=0,np-1
        bu0(i)=sum(B(i,0:np-1)*u0(0:np-1))
      end do
      
      tn=0.
      endt=5.0
      
      !allocate(I_cindex(0:nchoosek(8,4)-1))
      dt=endt/numstep
      !K=4
      !M=4
      
      call buildI(0,M,K,I_cindex)
      
      !NPC=nchoosek(K+M,K)
      !alpha=0
      do i=0,NPC-1
          call Itoalpha(I_cindex(i)%index(:),alpha(i,:),M,K)
      end do
      !Ic=0
      do i=0,NPC-1
          call alphaToI(alpha(i,:),Ic(i,:),M,K)
      end do
      
     
      call buildSkorohod
      D1=0.
      D2=0.
      do i=1,Np-1
      
          if(mod(i,2) .eq. 0) then
          D1(i) = -0.5/tan(0.5*x(i))
          D2(i) = -0.5/sin(0.5*x(i))/sin(0.5*x(i))
          else
          D1(i) = 0.5/tan(0.5*x(i))
          D2(i) = 0.5/sin(0.5*x(i))/sin(0.5*x(i))
          end if
          
      end do
      
      D2(0)=-(Np*Np+2)/12.0
      
      do i=0, Np-1
          do j=0,i-1
              A(i,j)=0.145*D2(j+Np-i)+0.1*sin(x(i))*D1(j+Np-i)
              B(i,j)=0.5*D1(j+Np-i)
          end do
          
          do j=i,Np-1
              A(i,j)=0.145*D2(j-i)+0.1*sin(x(i))*D1(j-i)
              B(i,j)=0.5*D1(j-i)
          end do
          
      end do
      
      
    identity=0.
    do i=0,np-1
        identity(i,i)=1.0
    end do
    
    idtA2=matmul(identity+0.5*transpose(A)*dt,identity+0.5*A*dt)+matmul(identity-0.5*transpose(A)*dt,identity-0.5*A*dt)
    idtA1=identity+0.5*A*dt
    idtA1t=transpose(idtA1)
    idtA=identity-0.5*A*dt
    
    idtA1tA=matmul(idtA1t,idtA)
    
    idtAtA1=transpose(idtA1tA)
    
    do i=0,numstep
        call LegendrePoly(i*dt,0.,endt,M-1,basis)
        
        Mk(i,0:M-1)=basis(0:M-1)
        
    end do
    
    tmp=idtA
     call dgetrf(np,np,tmp(0:np-1,0:np-1),np,ipiv(1:np),info)
     
     call dgetri(np,tmp(0:np-1,0:np-1),np,ipiv(1:np),work,4*(np),info) 
    
    idtAi=tmp
    
    idtAiA=matmul(idtAi,idtA1)
    !call compute_mineigen
      rmineigen=1.989093198067378E-002
      !do i=0,70-1
      !    write(*,"(4I8)") I_cindex(i)%index(:)
      !end do
      
      deallocate(tmp)
      
      end subroutine init
      
      
      subroutine stepForward
      use comdata
      
      real, allocatable,dimension(:,:)  :: u00,k1,k2,k3,k4
      !real, allocatable,dimension(:,:) :: k1
      
      allocate(u00(0:NPC-1,0:np-1),k1(0:NPC-1,0:np-1),k2(0:NPC-1,0:np-1),k3(0:NPC-1,0:np-1),k4(0:NPC-1,0:np-1))
      !k1=0.0
      !u00=u_num
      call getDiscretization(k1)
      u_num=u_num+k1*dt
      tn=tn+dt
      
      !call getDiscretization(k1)
      !
      !u_num=u_num+0.5*k1*dt
      !
      !tn=tn+0.5*dt
      !
      !call getDiscretization(k2)
      !
      !u_num=u00+0.5*k2*dt
      !
      !call getDiscretization(k3)
      !
      !u_num=u00+k3*dt
      !
      !tn=tn+0.5*dt
      !
      !call getDiscretization(k4)
      !
      !u_num=u00+(0.5*k1+k2+k3+0.5*k4)*dt/3.0
      
      
      
      deallocate(u00,k1,k2,k3,k4)
      
      end subroutine stepForward
      
      
      
      subroutine run(Nstep)
      use comdata
      integer :: pos
      
      integer :: Nstep,step
      dt =endt/Nstep
      
      !do step=1,Nstep
      !    
      !    call stepForward
      !    !u_ori(0:NPC-1,(step-1)*np+1:step*np)=u_num(0:NPC-1,0:np-1)
      !    write(11,*) 'Step=',step,'t=',tn
      !    
      !    !u1(:,step)=u_num(:)
      !    !write(*,*) sum(abs(u_num(:)-u(:,step)))/NPC
      !end do
      
      do pos=0,NPC-1
          call run_ori(pos)
          u_ori(pos,1:numstep_np)=unum(np:numstep_np1)
          write(11,*) 'pos=',pos
      end do
      
      call run_ori(1)
    
     !write(*,*) maxval(abs(u_ori(1,(numstep-1)*np+1:numstep*np)-unum(numstep*np:numstep_np1)))
    
      
      
      
      end subroutine run
      
      
      
      subroutine deletememory
      use comdata
      
      
      deallocate(u_num,Ic,alpha,s_eps,s_alpha,I_cindex)
     
      
      deallocate(f_LHS,r_Basis,C_basis,u_rbm,unum,basisnum,u0,frb)
      deallocate(pickup)
      deallocate(uAAu,Mu,uA,uAMu)
      deallocate(uAAu_i)
      deallocate(ipiv,Mk)
      
      deallocate(res_Q,res_R)
      !deallocate(AA)
      
      deallocate(A,B,D1,D2,x)
      
      deallocate(identity, idtA1t,idtA1,idtA2,idtA,idtAi,idtAiA)
      deallocate(idtA1tA,idtAtA1)
      
      deallocate(bu0)
      
      deallocate(u_ori)
      
      end subroutine
      
      
      subroutine reduce_basis
      use comdata
      integer :: pos
      
      !dt =endt/numstep
      
      call init_rbm
      
      !call assemble_mat
      
      
      !ns=1
      
      !do i=1,K
      !    
      !    call pickupbais(i)
      !    
      !end do
      
      
      do pos=1,NPC-1
          
          !call assemble_lhs(pos)
          !call run_ori(pos)
          call argminmize_coefficient(pos)
          !u_rbm(pos,0:np-1)=0.
          !do j=0,num_M-1
          !u_rbm(pos,0:np-1)=u_rbm(pos,0:np-1)+C_basis(pos,j)*r_basis((numstep-1)*np+1:numstep*np,j)
          !end do
          
          !write(*,*) maxval(abs(u_ori(pos,(numstep-1)*np+1:numstep*np)-u_rbm(pos,0:np-1)))
          !write(*,*)
          !pause
          !write(*,*) abs(u_rbm(pos,0:np-1)-u_num(pos,0:np-1))
      end do
      
      
          
      
      
      
      end subroutine reduce_basis
      
      subroutine outpdata
      
      use comdata
      
      integer :: i,j
      
      !open(unit=12,file='my_numu.dat')
      !do i=0,NPC-1
      !    write(12,"(10000ES26.16)") (u_rbm(i,j),j=1,numstep)
      !end do
      
      !close(12)
      
      !open(unit=13,file='numu.dat')
      !do i=0,NPC-1
      !    write(13,"(10000ES26.16)") (u1(i,j),j=1,numstep)
      !end do
      !
      !close(13)
      
      
    end subroutine outpdata
      
    
    
      function nchoosek(n,k)
      implicit none
      integer :: nchoosek,n,k
      integer :: a,i,j,ibeg,iend
      
      integer, allocatable, dimension(:) :: T

      if (n-k>k) then
          a=k
      else
          a=n-k
      end if
      
      allocate(T(0:a))
      
      T=1
      
      do i=0, n
          
          if (i>n-a) then
              ibeg=i-n+a
          else
              ibeg=1
          end if
          
          if (i>a) then
              iend=a
          else
              iend=i-1
          end if
          
          do j=iend,ibeg,-1
              T(j)=T(j)+T(j-1)
          end do
      end do
      
      
      nchoosek=T(a)
      
      return
          
      
      
      
      
    end function
    
    
    
    subroutine Itoalpha(Ic_tmp,alpha_tmp,M,K)
    implicit none
    
    integer :: K,M,v
    integer , dimension(0:K-1) :: Ic_tmp
    integer, dimension(0:M-1) :: alpha_tmp
    
    alpha_tmp=0
    
    do v=0, K-1
        if (Ic_tmp(v)>0) then
        alpha_tmp(Ic_tmp(v)-1)=alpha_tmp(Ic_tmp(v)-1)+1
        end if
        
    end do

    
    end subroutine Itoalpha
    
    
    subroutine alphaToI(alpha,I,M,K)
    implicit none
    
    integer :: M,pos,s,j,K
    integer , dimension (0:M-1) :: alpha
    integer , dimension (0:K-1) :: I
    
    I=0
    
    pos=0
    do s=M-1,0,-1
        do j=0,alpha(s)-1
            I(pos)=s+1
            pos=pos+1
        end do
    end do
    
    
    
    end subroutine alphaToI
    
    
    subroutine IToPos(I,M,K,pos)
    implicit none
    
    integer :: pos,M,K
    integer, dimension(0:K-1) ::I
    
    integer :: j,v,nchoosek
    
    pos=0
    
    do j=K-1,0,-1
        if (j == K-1) then
            do v=0, I(j)-1
                pos=pos+nchoosek(M+j-v,j)
            end do
        else
            do v=I(j+1),I(j)-1
                pos=pos+nchoosek(M+j-v,j)
            end do
        end if
        
    end do
    
    return
        
    
    
    
    end subroutine IToPos
    
    
    subroutine alphaToPos(alpha,M,K,pos)
    implicit none
    
    integer :: M,K,pos
    integer , dimension (0:M-1) :: alpha
    integer, dimension(0:K-1) ::I
    
    call alphaToI(alpha,I,M,K)
    
    call IToPos(I,M,K,pos)
    
    return
    
    end subroutine alphaToPos
    
    
    
    subroutine buildSkorohod
    use comdata
    
    integer :: i,j,pos
    integer ,dimension (0:M-1):: alpha_eps
    integer ,allocatable, dimension(:) :: tmp
    
    !deallocate(s_eps,s_alpha)
    
    !allocate(s_eps(0:NPC-1),s_alpha(0:NPC-1))
    
    !do i=0,NPC-1
    !
    !    allocate(s_eps(i)%index(0:0),s_alpha(i)%index(0:0))
    !    s_eps(i)%index(0)=0
    !    s_alpha(i)%index(0)=0
    !end do
    

    
    
    do i=0,NPC-1
        do j=0,M-1
            if (alpha(i,j)>0) then
                alpha_eps(:)=alpha(i,:)
                alpha_eps(j)=alpha_eps(j)-1
                nsize=s_eps(i)%index(0)
                nsize=nsize+1
                allocate(tmp(1:nsize))
                if (nsize>1) then
                tmp(1:nsize-1)=s_eps(i)%index(1:nsize-1)
                tmp(nsize)=j
                else
                tmp(nsize)=j
                end if
                
                deallocate(s_eps(i)%index)
                allocate(s_eps(i)%index(0:nsize))
                s_eps(i)%index(0)=nsize
                s_eps(i)%index(1:nsize)=tmp(:)
                
                
                deallocate(tmp)
                nsize=s_alpha(i)%index(0)
                nsize=nsize+1
                allocate(tmp(1:nsize))
                if (nsize>1) then
                tmp(1:nsize-1)=s_alpha(i)%index(1:nsize-1)
                call alphaToPos(alpha_eps,M,K,pos)
                tmp(nsize)=pos
                else
                call alphaToPos(alpha_eps,M,K,pos)
                tmp(nsize)=pos
                end if
                
                deallocate(s_alpha(i)%index)
                allocate(s_alpha(i)%index(0:nsize))
                s_alpha(i)%index(0)=nsize
                s_alpha(i)%index(1:nsize)=tmp(:)
                deallocate(tmp)
            end if
        end do
    end do
    
                
                
                
                
    
    
    
    end subroutine buildSkorohod
    
    
    
    
    !subroutine assemble_mat
    !use comdata
    !
    !A_mass=0.
    !A_mass(1,1)=1.
    !do i=2,numstep
    !A_mass(i,i)=1.
    !A_mass(i,i-1)=-(1+dt)
    !end do
    !
    !
    !
    !
    !
    !end subroutine
    
    
    subroutine assemble_lhs_rbm(pos)
    use comdata
    
    integer :: pos
    real , dimension(0:np-1) :: tmp1,tmp2
     
     
    f_LHS=0.

        
        do l=1,s_eps(pos)%index(0)
            
            
            if (s_alpha(pos)%index(l) .eq. 0) then
                tmp1=0.
            do i=0,num_M-1
                tmp1(0:np-1)=tmp1(0:np-1)+C_basis(s_alpha(pos)%index(l),i)*r_basis(0*np+1:1*np,i)
            end do
            do j=0,np-1
          f_LHS(j)=f_LHS(j)+0.5*dt*sqrt(real(alpha(pos,s_eps(pos)%index(l))))*(bu0(j)*Mk(0,s_eps(pos)%index(l)) + &
          sum(B(j,0:np-1)*tmp1(0:np-1))*Mk(1,s_eps(pos)%index(l)))     
            end do
            end if
            
        end do
    
        do j=0,np-1
            f_LHS(j)=sum(idtAi(j,0:np-1)*f_LHS(0:np-1))
        end do
        

        
    do i=2,numstep

        do l=1,s_eps(pos)%index(0)
            tmp1=0.
            tmp2=0.
            do j=0,num_M-1
                tmp1(0:np-1)=tmp1(0:np-1)+C_basis(s_alpha(pos)%index(l),j)*r_basis((i-2)*np+1:(i-1)*np,j)
                tmp2(0:np-1)=tmp2(0:np-1)+C_basis(s_alpha(pos)%index(l),j)*r_basis((i-1)*np+1:(i)*np,j)
            end do
            
            do j=0,np-1
          f_LHS((i-1)*np+j)=f_LHS((i-1)*np+j)+sqrt(real(alpha(pos,s_eps(pos)%index(l))))* &
         (sum(B(j,0:np-1)*tmp1(0:np-1))*Mk(i-1,s_eps(pos)%index(l))+&
          sum(B(j,0:np-1)*tmp2(0:np-1))*Mk(i,s_eps(pos)%index(l)))
            end do
        end do
        f_LHS((i-1)*np:i*np-1)=0.5*dt*f_LHS((i-1)*np:i*np-1)
        
        
        do j=0,np-1
            f_LHS((i-1)*np+j)=sum(idtAi(j,0:np-1)*f_LHS((i-1)*np:i*np-1))
        end do
        
    end do
    
    
    
    end subroutine
       
    
    subroutine assemble_lhs(pos)
    use comdata
    
    integer :: pos
    !real, dimension (0:M-1) :: basis
     
     
     
    f_LHS=0.
    !tn=(0)*dt
        !call LegendrePoly(tn,0.,endt,M-1,basis)
        
        do l=1,s_eps(pos)%index(0)
            
        if (s_alpha(pos)%index(l) .eq. 0) then
        do j=0,np-1
          f_LHS(j)=f_LHS(j)+0.5*dt*sqrt(real(alpha(pos,s_eps(pos)%index(l))))*(bu0(j)*Mk(0,s_eps(pos)%index(l)) + &
          sum(B(j,0:np-1)*u_ori(0,1:np))*Mk(1,s_eps(pos)%index(l)))     
        end do
        end if
            
        end do
    
        do j=0,np-1
            f_LHS(j)=sum(idtAi(j,0:np-1)*f_LHS(0:np-1))
        end do
        

        
    do i=2,numstep
        !tn=(i-1)*dt
        !call LegendrePoly(tn,0.,endt,M-1,basis)
        do l=1,s_eps(pos)%index(0)
            do j=0,np-1
          f_LHS((i-1)*np+j)=f_LHS((i-1)*np+j)+sqrt(real(alpha(pos,s_eps(pos)%index(l))))* &
         (sum(B(j,0:np-1)*u_ori(s_alpha(pos)%index(l),(i-2)*np+1:(i-1)*np))*Mk(i-1,s_eps(pos)%index(l))+&
          sum(B(j,0:np-1)*u_ori(s_alpha(pos)%index(l),(i-1)*np+1:i*np))*Mk(i,s_eps(pos)%index(l)))
            end do
        end do
        f_LHS((i-1)*np:i*np-1)=0.5*dt*f_LHS((i-1)*np:i*np-1)
        
        
        do j=0,np-1
            f_LHS((i-1)*np+j)=sum(idtAi(j,0:np-1)*f_LHS((i-1)*np:i*np-1))
        end do
        
    end do
    
    
    !do i=2,numstep
    !    !tn=(i-1)*dt
    !    !call LegendrePoly(tn,0.,endt,M-1,basis)
    !    do l=1,s_eps(pos)%index(0)
    !        do j=0,np-1
    !      f_LHS((i-1)*np+j)=f_LHS((i-1)*np+j)+sqrt(real(alpha(pos,s_eps(pos)%index(l))))* &
    !     Mk(i-1,s_eps(pos)%index(l))*sum(B(j,0:np-1)*u_ori(s_alpha(pos)%index(l),(i-2)*np+1:(i-1)*np))
    !        end do
    !    end do
    !    f_LHS((i-1)*np:i*np-1)=dt*f_LHS((i-1)*np:i*np-1)
    !    
    !end do
    
    
    end subroutine
   
    subroutine assemble_frb(pos)
    use comdata
    
    integer :: pos
    
    !real, allocatable, dimension(:) :: tmp
    
    !allocate(tmp(1:num_M))
    
    !frb(1:num_M)=uA(1:num_M,1)*(u0(pos)+dt*u0(pos))
    frb=0.
    do l=1,s_eps(pos)%index(0)
        
        if (s_alpha(pos)%index(l) .eq. 0) then
            do i=1,num_M
        frb(i)=frb(i)+0.5*dt*sqrt(real(alpha(pos,s_eps(pos)%index(l))))*Mk(0,s_eps(pos)%index(l))*&
        sum(uA(i,1:np)*bu0(0:np-1))
            end do
        end if
        
        do i=1,num_M
          frb(i)=frb(i)+0.5*dt*sqrt(real(alpha(pos,s_eps(pos)%index(l))))* &
         sum(uAMu(s_eps(pos)%index(l),i,1:num_M)*C_basis(s_alpha(pos)%index(l),0:num_M-1))
        end do
    end do
    
    !do i=1,num_M
    !tmp(i)=sum(uA(i,1:numstep)*f_LHS(1:numstep))
    !end do
    
    !write(*,*) frb(1:num_M)-tmp(1:num_M)
    
    !deallocate(tmp)
    
    end subroutine
    
    
    ! compute the coefficient c_alpha
    subroutine argminmize_coefficient(pos)

    use comdata
    
    !real, allocatable, dimension(:,:) :: tmp,tmp2
    real work(1:4*num_M)
    INTEGER info,IWORK(1:num_M)
    real, allocatable, dimension(:,:) :: tmp_c
    integer pos
    
    !allocate(tmp(1:num_M,1:num_M),tmp_c(1:num_M))
    !allocate(tmp2(1:num_M,1:numstep))
    allocate(tmp_c(1:num_M,1:1))
    
    
     !call update_RBmatrix(pos)
     !tmp=uAAu
     !
     !call dgetrf(num_M,num_M,uAAu,num_M,ipiv,info)
     !call dgetri(num_M,tmp,num_M,ipiv,work,4*num_M,info)
     
     
!     ns=1
!        do i=2,K
!            nend=ns+nchoosek(M+i-1,i)
!            do pos=ns,nend-1
!                
!                call assemble_frb(pos)
!                
!                do j=1,num_M
!                tmp_c(j,1)=sum(uAAu_i(j,1:num_M)*frb(1:num_M))
!                end do
!                
!                resrb=sum(abs(tmp_c(1:num_M,1)))
!            
!                write(*,*) resrb
!                
!            end do
!            ns=nend
!        end do
        
        
     
     
     call assemble_frb(pos)
     !tmp_c(1:num_M,1)=frb(1:num_M)
     !call dgetrs('N',num_M,1,uAAu_i(1:num_M,1:num_M),num_M,ipiv(1:num_M),tmp_c(1:num_M,1),num_M,info)
     
     do i=1,num_M
     tmp_c(i,1)=sum(uAAu_i(i,1:num_M)*frb(1:num_M))
     end do
     
     !call assemble_lhs(pos)
     
     call residue(tmp_c(1:num_M,1),resrb,pos)
     
     !resrb=sum(abs(tmp_c(1:num_M,1)))
     
     !if (resrb/rmineigen<epsilon) then
     if (resrb/rmineigen<epsilon .and. pos > M) then
      !if(pickup(pos) .eq. 1) then 
          
         C_basis(pos,0:num_M-1)=tmp_c(1:num_M,1)
         write(13,*) 'pos=',pos,'resrb=',resrb
         
      else
     
         write(13,*) 'pos=',pos,'resrb=',resrb
         call run_ori_rbm(pos)
         
         call update_rbmatrix(pos)
         
         !num_M=num_M+1
         !C_basis(pos,num_M-1)=1.0
         
         !r_basis(numstep,num_M-1)=unum(numstep)
         
         if (num_Mup .eq. 1) then
         
         write(13,*) "bais number=", pos
         
         write(12,*) pos
         write(13,*)'alpha='
         write(13,'(100I2)')(alpha(pos,i),i=0,M-1)
         write(*,*) "bais number=", pos
         write(*,*)'alpha='
         write(*,'(100I2)')(alpha(pos,i),i=0,M-1)
         basisnum(num_M)=pos
         
         pickup(pos)=0
         
         else
         write(2020,*) pos
         write(2020,*)'alpha='
         write(2020,'(100I2)')(alpha(pos,i),i=0,M-1)
         
         end if
         
         
     endif
     
     
     
    
    
!    AU(1,1:num_M)=r_basis(1,0:num_M-1)
!    do i=2,numstep
!        do j=0,num_M-1
!            AU(i,j+1)=sum(A_mass(i,i-1:i)*r_basis(i-1:i,j))
!        end do
!    end do
!12    
!     tmp=Matmul(Transpose(AU),AU)
!     
!     call dgetrf(num_M,num_M,tmp,num_M,ipiv,info)
!     call dgetri(num_M,tmp,num_M,ipiv,work,4*num_M,info)
!     
!     tmp2=matmul(tmp,Transpose(AU))
!     do i=1,num_M
!     tmp_c(1,i)=sum(f_LHS(1:numstep)*tmp2(i,1:numstep))
!     end do
!     
!     do i=1,numstep
!     tmp2(1,i)=sum(AU(i,:)*tmp_c(1,:))
!     end do
!     
!     
!     r_min2=sqrt(sum(abs(f_LHS(1:numstep)-tmp2(1,1:numstep))**2))
!     
!     !write(*,*) r_min2
!     
!     if (r_min2<epsilon) then
!         
!         C_basis(pos,0:num_M-1)=tmp_c(1,1:num_M)
!         
!         do i=0,num_M-1
!             u_rbm(pos,1:numstep)=u_rbm(pos,1:numstep)+tmp_c(1,i+1)*r_basis(1:numstep,i)
!         end do
!         
!         
!         
!         
!     else
!         
!    
!         call run_ori(pos)
!         
!         num_M=num_M+1
!         
!         r_basis(1:numstep,num_M-1)=unum(1:numstep)
!         
!         C_basis(pos,num_M-1)=1.0
!         
!         u_rbm(pos,1:numstep)=unum(1:numstep)
!         
!         write(*,*) "bais number=", pos
!         write(*,*)'alpha='
!         write(*,'(100I2)')(alpha(pos,i),i=0,M-1)
!         
!         basisnum(num_M)=pos
!         
!         pickup(pos)=0
!     end if
     
    
    
    deallocate(tmp_c)
    end subroutine
!    
    
    
    subroutine pickupbais(n_alpha)
    use comdata
    
    real, allocatable, dimension(:,:) :: tmp_c
    integer :: nchoosek,n_alpha,pos
    integer :: pos_pick(1:M)
    real, allocatable, dimension(:) :: c_normL1
    real maxnorm(1:M)
    integer maxpos(1:M)
    
    real :: minnorm
    integer :: minlocal
    
    real work(1:4*(num_M+M))
    INTEGER info,IWORK(1:num_M+M)
    real, allocatable, dimension(:,:) :: tmp
    
    allocate(tmp(1:num_M+M,1:num_M+M))
    
    
    numalpha=nchoosek(M+n_alpha-1,n_alpha)
    nend=ns+numalpha
    
    allocate(tmp_c(1:num_M,1:numalpha))
    allocate(c_normL1(1:numalpha))
    
    
    do pos=ns,nend-1
        
        call assemble_frb(pos)
        
        do i=1,num_M
        tmp_c(i,pos-ns+1)=sum(uAAu_i(i,1:num_M)*frb(1:num_M))
        end do
        !c_normL1(pos-ns+1)=sum(abs(tmp_c(1:num_M,pos-ns+1)))
        
    end do
    
!    maxnorm(1:M)=c_normL1(1:M)
!    do i=1,M
!        maxpos(i)=ns+i-1
!    end do
!    minnorm=maxnorm(1)
!    !minpos=ns
!    minlocal=1
!    
!    do j=2,M
!        if (minnorm> maxnorm(j)) then
!            minnorm=maxnorm(j)
!            !minpos=ns+j-1
!            minlocal=j
!        end if
!    end do
!    
!        
!    
!    do i=ns+M,nend-1
!        
!     if (c_normL1(i-ns+1)>minnorm) then
!         
!        C_basis(maxpos(minlocal),0:num_M-1)=tmp_c(1:num_M,maxpos(minlocal)-ns+1)
!        maxnorm(minlocal)=c_normL1(i-ns+1)
!        maxpos(minlocal)=i
!            
!    minnorm=maxnorm(1)
!    !minpos=ns
!    minlocal=1
!    
!    do j=2,M
!        if (minnorm > maxnorm(j)) then
!            minnorm=maxnorm(j)
!            !minpos=ns+j-1
!            minlocal=j
!        end if
!    end do
!     
!     else
!         C_basis(i,0:num_M-1)=tmp_c(1:num_M,i-ns+1)
!         
!     end if
!        
!    end do
!
!    
!    
!    do j=1,M
!        pos=maxpos(j)
!        call run_ori(pos)
!        
!        call update_rbmatrix(pos)
!         
!        num_M=num_M+1
!        C_basis(pos,num_M-1)=1.0
!         
!         r_basis(1:numstep,num_M-1)=unum(1:numstep)
!         
!         write(*,*) "bais number=", pos
!         
!         !write(12,*) pos
!         !write(*,*)'alpha='
!         !write(*,'(100I2)')(alpha(pos,i),i=0,M-1)
!         
!         basisnum(num_M)=pos
!         
!         pickup(pos)=0
!         
!    end do
!    
!     tmp(1:num_M,1:num_M)=uAAu(1:num_M,1:num_M)
!     
!     
!     call dgetrf(num_M,num_M,tmp(1:num_M,1:num_M),num_M,ipiv(1:num_M),info)
!     call dgetri(num_M,tmp,num_M,ipiv(1:num_M),work,4*(num_M),info)
!     
!     
!     uAAu_i(1:num_M,1:num_M)=tmp(1:num_M,1:num_M)
!    
!    ns=nend
    
    
    num_Mold=num_M
    
     do pos=ns,nend-1
        
        if (alpha(pos,M-1) .eq. n_alpha-1 .or. pos .eq. nend-1) then
            

            call run_ori(pos)
        
            call update_rbmatrix(pos)
         
            num_M=num_M+1
            C_basis(pos,num_M-1)=1.0
         
            r_basis(numstep,num_M-1)=unum(numstep)
         
            !write(*,*) "bais number=", pos
            !write(*,*) "C_alpha_L1=",c_normL1(pos-ns+1)
            !write(*,*) "alpha="
            !write(*,'(100I2)')(alpha(pos,i),i=0,M-1)
            !write(*,*)
         
            basisnum(num_M)=pos
         
            pickup(pos)=0
        else
             C_basis(pos,0:num_Mold-1)=tmp_c(1:num_Mold,pos-ns+1)
             
        end if
        
     end do
     
        
        tmp(1:num_M,1:num_M)=uAAu(1:num_M,1:num_M)
     
     
        call dgetrf(num_M,num_M,tmp(1:num_M,1:num_M),num_M,ipiv(1:num_M),info)
        call dgetri(num_M,tmp,num_M,ipiv(1:num_M),work,4*(num_M),info)
     
     
        uAAu_i(1:num_M,1:num_M)=tmp(1:num_M,1:num_M)
    
        ns=nend
    
    
    
     deallocate(tmp,tmp_c,c_normL1)
    end subroutine
    
    
    ! residue term ||zf_alpha-AUc_alpha||
    subroutine residue(c_alpha,resrb,pos)
    use comdata
    real :: resrb
    real, dimension(1:num_M) :: c_alpha
    !real, dimension(1:num_M,1:1) :: tmp
    real, dimension(1:nrank) :: tmp2
    real :: tmp
    integer :: pos
    
    resrb=0.
    tmp2=0.
    
    do i=0,num_M-1
    do j=1,(i+1)*(M+1)
        tmp2(j)=tmp2(j)+c_alpha(i+1)*res_R(j,i*(M+1)+1)
        do l=1,s_eps(pos)%index(0)
            
            tmp2(j)=tmp2(j)-0.5*dt*sqrt(real(alpha(pos,s_eps(pos)%index(l))))*C_basis(s_alpha(pos)%index(l),i)*res_R(j,i*(M+1)+s_eps(pos)%index(l)+2)
        
        end do
        
    end do
    
    !tmp2(j)=tmp
    
    end do
    
    resrb=sqrt(sum(tmp2(1:nrank)**2))
    
    !f0=0.
    !
    !do l=1,s_eps(pos)%index(0)
    !        
    !        if (s_alpha(pos)%index(l) .eq. 0) then
    !        do j=0,np-1
    !      f0(j)=f0(j)+sqrt(real(alpha(pos,s_eps(pos)%index(l))))* &
    !     sum(u0(0:np-1)*B(j,0:np-1))*Mk(0,s_eps(pos)%index(l))*dt
    !        end do
    !        end if
    !        
    !end do
    !
    !do i=1,nrank
    !resrb=resrb+(tmp2(i)-sum(res_Q(1:np,i)*f0(0:np-1)))**2
    !end do
    
    
    !resrb=resrb+rq1*f0**2
    
    !resrb=sqrt(resrb)
    
    !tmp=0.
    !do i=0,num_M-1
    !    tmp=tmp+c_alpha(i+1)*tmp2(i*(M+1)+1)
    !    do l=1,s_eps(pos)%index(0)
    !        
    !        tmp=tmp+C_basis(s_alpha(pos)%index(l),i)*tmp2(i*(M+1)+s_alpha(pos)%index(l)+2)
    !    end do
    !    
    !end do
    
    
    !resrb=sqrt(sum(tmp2(1:nrank)**2))
    
    
    
    end subroutine
    
    
    !grow matrix uA,u^TA^TAu,...
    subroutine update_rbmatrix(pos)
    use comdata
    integer :: pos
    real work(1:4*(num_M+1))
    INTEGER info,IWORK(1:num_M+1)
    real, allocatable, dimension(:,:) :: tmp
    
    real, allocatable, dimension(:) :: tmpvec1,tmpvec2
    real :: tmp1,tmp2
    
    allocate(tmpvec1(1:numstep_np),tmpvec2(1:numstep_np))
    
    allocate(tmp(1:num_M+1,1:num_M+1))
    
    
    !do i=1, num_M
    !
    !    uAAu(num_M+1,i)=sum(unum(1:numstep)*AAu(1:numstep,i))
    !    uAAu(i,num_M+1)=uAAu(num_M+1,i)
    !end do
    !
    !uAAu(num_M+1,num_M+1)=(1+(1+dt)**2)*sum(unum(1:numstep-1)*unum(1:numstep-1))+&
    !    -2*(1+dt)*sum(unum(1:numstep-1)*unum(2:numstep))+unum(numstep)**2
    
    !AAu(1,num_M+1)=(1+(1+dt)**2)*unum(1)-(1+dt)*unum(2)
    !AAu(2:numstep-1,num_M+1)=(1+(1+dt)**2)*unum(2:numstep-1)-(1+dt)*unum(3:numstep)-(1+dt)*unum(1:numstep-2)
    !AAu(numstep,num_M+1)=unum(numstep)-(1+dt)*unum(numstep-1)
    
    
        C_basis(pos,num_M)=1.
        do j=1,num_M
            
            tmp1=sum(unum(np:numstep_np1)*r_basis(1:numstep_np,j-1))
            tmpvec2(1:numstep_np)=unum(np:numstep_np1)-tmp1*r_basis(1:numstep_np,j-1)
    
            tmp2=sqrt(sum(tmpvec2(1:numstep_np)**2))
            
            C_basis(pos,j-1)=tmp1*C_basis(pos,num_M)
            
            
            if (tmp2 < 1.0E-17) then
                
                num_Mup=0
                C_basis(pos,j:num_M)=0.
                exit
                
            else
                unum(np:numstep_np1)=tmpvec2(1:numstep_np)/tmp2
                num_Mup=1
                C_basis(pos,num_M)=tmp2*C_basis(pos,num_M)
                
            end if
            
        end do
        
    if (num_Mup .eq. 1) then
            num_M=num_M+1
            r_basis(1:numstep_np,num_M-1)=unum(np:numstep_np1)
      
            do i=1,np
                uA(num_M,i)=sum(idtA(i-1,0:np-1)*unum(np:2*np-1))
            end do
    do i=2,numstep
        do j=0,np-1
    uA(num_M,(i-1)*np+j+1)=sum(idtA(j,0:np-1)*unum(i*np:(i+1)*np-1))-sum(idtA1t(0:np-1,j)*unum((i-1)*np:i*np-1))
        end do
    end do
    
    do i=1,num_M
    uAAu(num_M,i)=sum(uA(num_M,1:numstep_np)*uA(i,1:numstep_np))
    uAAu(i,num_M)=uAAu(num_M,i)
    end do
    
    uAAu(num_M,num_M)=sum(uA(num_M,1:numstep_np)*uA(num_M,1:numstep_np))
    
    !uA(num_M,1:numstep)=Au(1:numstep,num_M)
    
    
    do i=0, M-1
        do j=0,np-1
        Mu(i,j+1,num_M)=Mk(1,i)*sum(B(j,0:np-1)*unum(np:2*np-1))
        end do
        
        do ns=2,numstep
            do j=0,np-1
            Mu(i,(ns-1)*np+j+1,num_M)=Mk(ns-1,i)*sum(B(j,0:np-1)*unum((ns-1)*np:ns*np-1))+mk(ns,i)*sum(B(j,0:np-1)*unum(ns*np:(ns+1)*np-1))
            end do
        end do
    
    do j=1,num_M
    uAMu(i,num_M,j)=sum(uA(num_M,1:numstep_np)*Mu(i,1:numstep_np,j))!sum(unum(2:numstep)*Mk(1:numstep-1,i)*r_basis(2:numstep,j-1))-(1+dt)*sum(r_basis(3:numstep,j-1)*unum(2:numstep-1)*Mk(1:numstep-2,i))
    !uMAu(i,j,num_M)=uAMu(i,num_M,j)
    uAMu(i,j,num_M)=sum(uA(j,1:numstep_np)*Mu(i,1:numstep_np,num_M))!sum(unum(2:numstep)*Mk(1:numstep-1,i)*r_basis(2:numstep,j-1))-(1+dt)*sum(r_basis(2:numstep-1,j-1)*unum(3:numstep)*Mk(1:numstep-2,i))
    !uMAu(i,num_M,j)=uAMu(i,j,num_M)
    end do
    
    uAMu(i,num_M,num_M)=sum(uA(num_M,1:numstep_np)*Mu(i,1:numstep_np,num_M))!sum(unum(2:numstep)*Mk(1:numstep-1,i)*unum(2:numstep))-(1+dt)*sum(unum(3:numstep)*unum(2:numstep-1)*Mk(1:numstep-2,i))
    
    
    end do
    
    
    
    
    
    
    do i=-1,M-1
        
        if (i .eq. -1) then
            tmpvec1(1:numstep_np)=uA(num_M,1:numstep_np)
        else
            
            tmpvec1(1:numstep_np)=Mu(i,1:numstep_np,num_M)
        end if
        
        res_R(nrank+1,ncolumn+1)=1.0
        do j=1,nrank
            
            tmp1=sum(tmpvec1(1:numstep_np)*res_Q(1:numstep_np,j))
            tmpvec2(1:numstep_np)=tmpvec1(1:numstep_np)-tmp1*res_Q(1:numstep_np,j)
    
            tmp2=sqrt(sum(tmpvec2(1:numstep_np)**2))
            
            res_R(j,ncolumn+1)=tmp1*res_R(nrank+1,ncolumn+1)
            
            
            if (tmp2 < 1.0E-17) then
                
                nrankup=0
                res_R(nrank+1,ncolumn+1)=0.
                exit
                
            else
                tmpvec1(1:numstep_np)=tmpvec2(1:numstep_np)/tmp2
                nrankup=1
                res_R(nrank+1,ncolumn+1)=tmp2*res_R(nrank+1,ncolumn+1)
                
            end if
            
        end do
        
        if (nrankup .eq. 1) then
            nrank=nrank+1
            
            res_Q(1:numstep_np,nrank)=tmpvec1(1:numstep_np)
            
            !do j=1,ncolumn
            !    res_RTR(ncolumn+1,j)=sum(res_R(1:nrank-1,ncolumn+1)*res_R(1:nrank-1,j))
            !    res_RTR(j,ncolumn+1)=res_RTR(ncolumn+1,j)
            !end do
            !
            !res_RTR(ncolumn+1,ncolumn+1)=sum(res_R(1:nrank,ncolumn+1)**2)
            
        !else
        !    
        !    do j=1,ncolumn
        !        res_RTR(ncolumn+1,j)=sum(res_R(1:nrank,ncolumn+1)*res_R(1:nrank,j))
        !        res_RTR(j,ncolumn+1)=res_RTR(ncolumn+1,j)
        !    end do
        !    res_RTR(ncolumn+1,ncolumn+1)=sum(res_R(1:nrank,ncolumn)**2)
            
        end if
        
        ncolumn=ncolumn+1
        
        
        
    end do
    
    
    
     tmp(1:num_M,1:num_M)=uAAu(1:num_M,1:num_M)
     
     call dgetrf(num_M,num_M,tmp(1:num_M,1:num_M),num_M,ipiv(1:num_M),info)
     
     call dgetri(num_M,tmp,num_M,ipiv(1:num_M),work,4*(num_M),info)
     
     
     uAAu_i(1:num_M,1:num_M)=tmp(1:num_M,1:num_M)
     
     
     
!     do i=0,np-1
!         do j=0,np-1
!     
!     rq1(i,j)=identity(i,j)-sum(res_Q(i,1:nrank)*res_Q(j,1:nrank))
!    
!    do i=2,numstep_np
!        rq1=rq1+sum(res_Q(i,1:nrank)*res_Q(1,1:nrank))**2
!    end do
    
    
    end if
        
     
     
     
     deallocate(tmp,tmpvec1,tmpvec2)
     
    end subroutine
    
    
    ! "truth" solution solver
    subroutine run_ori(pos)
    use comdata

      real, dimension (0:M-1) :: basis
      integer :: i,pos
      

      if (pos .eq. 0) then
      !res=unum(0)+1.0
      unum(0:np-1)=u0(0:np-1)
      do i=1,numstep
          do j=0,Np-1
              unum(i*Np+j)=sum(idtAiA(j,0:np-1)*unum((i-1)*np:i*np-1))
          end do
          
      end do
      
      !u1(pos,1:numstep)=unum(1:numstep)
      !u_rbm(pos,1:numstep)=unum(1:numstep)
      
      else
          
          call assemble_lhs(pos)
          
          unum(np:2*np-1)=f_LHS(0:np-1)
          do i=2,numstep
              do j=0,np-1
          unum(i*np+j)=f_LHS((i-1)*np+j)+sum(idtAiA(j,0:np-1)*unum((i-1)*np:i*np-1))
              end do
          end do
          
          !u_ori(pos,1:numstep_np)=unum(np:numstep_np1)
          !u1(pos,1:numstep)=unum(1:numstep)
          !u_rbm(pos,0:np-1)=unum(numstep_np:numstep_np1)
          !write(*,*) maxval(abs(u_rbm(pos,1:numstep)-u1(pos,1:numstep)))
      end if
      
      
      
    end subroutine
    
    
        subroutine run_ori_rbm(pos)
    use comdata

      real, dimension (0:M-1) :: basis
      integer :: i,pos
      

      if (pos .eq. 0) then
      !res=unum(0)+1.0
      do i=1,numstep
          do j=0,Np-1
              unum(i*Np+j)=sum(idtAiA(j,0:np-1)*unum((i-1)*np:i*np-1))
          end do
          
      end do
      
      !u1(pos,1:numstep)=unum(1:numstep)
      !u_rbm(pos,1:numstep)=unum(1:numstep)
      
      else
          
          call assemble_lhs_rbm(pos)
          
          unum(np:2*np-1)=f_LHS(0:np-1)
          do i=2,numstep
              do j=0,np-1
          unum(i*np+j)=f_LHS((i-1)*np+j)+sum(idtAiA(j,0:np-1)*unum((i-1)*np:i*np-1))
              end do
          end do
          
          !u_ori(pos,1:numstep_np)=unum(np:numstep_np1)
          !u1(pos,1:numstep)=unum(1:numstep)
          !u_rbm(pos,0:np-1)=unum(numstep_np:numstep_np1)
          !write(*,*) maxval(abs(u_rbm(pos,1:numstep)-u1(pos,1:numstep)))
      end if
      
      
      
    end subroutine
    
    
    subroutine compute_mineigen
    use comdata
    
    real , allocatable, dimension(:) :: w, work
    integer :: info,lwork,lda
    
    
    allocate(w(1:numstep_np),work(1:3*numstep_np-1))
    !epsilon=1.0E-04
    lda=numstep_np
    lwork=3*numstep_np-1
    info=0
    
    
    
    AA=0.
    do i=1,numstep-1
    AA((i-1)*np+1:i*np,(i-1)*np+1:i*np)=idtA2
    AA((i-1)*np+1:i*np,(i)*np+1:(i+1)*np)=-idtA1tA
    AA((i)*np+1:(i+1)*np,(i-1)*np+1:(i)*np)=-idtAtA1
    end do
    AA((numstep-1)*np+1:numstep*np,(numstep-2)*np+1:(numstep-1)*np)=-idtAtA1
    AA((numstep-1)*np+1:numstep*np,(numstep-1)*np+1:numstep*np)=matmul(transpose(idtA),idtA)
    
    
    call dsyev('N','U',numstep,AA,lda,w,work,lwork,info)
    
    rmineigen=sqrt(w(1))
    
    write(2018,*) 'min sqrt eigenvalue AtA=', rmineigen
    deallocate(w,work)
    
    end subroutine
    
    ! init reduced basis method
    subroutine init_rbm
    use comdata
    
    real, dimension(0:M-1) :: basis
    
    real, allocatable, dimension(:) :: tmpvec1,tmpvec2!,w,work
    real :: tmp1,tmp2
    integer :: pos,info,lwork,lda
    
    allocate(tmpvec1(1:numstep_np),tmpvec2(1:numstep_np))
    !allocate(w(1:numstep),work(1:3*numstep-1))
    !epsilon=1.0E-04
    lda=numstep_np
    lwork=3*numstep_np-1
    info=0
    
    
    unum=0.
    
    
    do i=0,NP-1
    unum(i)=u0(i)
    end do
    
    
    
    
    !call LegendrePoly(0.,0.,endt,M-1,basis0)
    
    
    num_M=1
    r_basis=0.
    
    call run_ori_rbm(0)
    
    !write(*,*) maxval(abs(u_ori(0,(numstep-1)*np+1:numstep*np)-unum(numstep*np:numstep_np1)))
    !
    !call run_ori(1)
    !
    !write(*,*) maxval(abs(u_ori(1,(numstep-1)*np+1:numstep*np)-unum(numstep*np:numstep_np1)))
    
    
    !u_ori(0,1:numstep_np)=unum(np:numstep_np1)
    
    rnorm=sqrt(sum(unum(np:numstep_np1)**2))
    unum(np:numstep_np1)=unum(np:numstep_np1)/rnorm
    
    !---------construct Arb
    
    !uAAu(1,1)=(1.+(1+dt)**2)*sum(unum(1:numstep-1)*unum(1:numstep-1))+unum(numstep)*unum(numstep)&
    !    -2*(1.+dt)*sum(unum(1:numstep-1)*unum(2:numstep))
    !
    !uAAu_i(1,1)=1./uAAu(1,1)
    
    
    
    
    !call compute_mineigen
    
    
    !do j=0,np-1
    !AAu(j+1,1)=sum(idtA2(j,0:np-1)*unum(np:2*np-1))-sum(idtA1t(j,0:np-1)*unum(2*np:3*np-1))
    !end do
    !do i=2,numstep-1
    !    do j=0,np-1
    !        
    !AAu((i-1)*np+j+1,1)=sum(idtA2(j,0:np-1)*unum(i*np:(i+1)*np-1))-sum(idtA1t(j,0:np-1)*unum((i+1)*np:(i+2)*np-1))-sum(idtA1(j,0:np-1)*unum((i-1)*np:(i)*np-1))
    !
    !    end do
    !end do
    !do j=0,np-1
    !AAu((numstep-1)*np+j+1,1)=unum(numstep*np+j)-sum(idtA1(j,0:np-1)*unum((numstep-1)*np:numstep*np-1))
    !end do
    
    !uAA(1,1:numstep)=AAu(1:numstep,1)
    
    uA(1,1:np)=matmul(idtA(0:np-1,0:np-1),unum(np:2*np-1))
    do i=2,numstep
        do j=0,np-1
    uA(1,(i-1)*np+j+1)=sum(idtA(j,0:np-1)*unum(i*np:(i+1)*np-1))-sum(idtA1t(0:np-1,j)*unum((i-1)*np:i*np-1))
        end do
    end do
    
    uAAu(1,1)=sum(uA(1,1:numstep_np)*uA(1,1:numstep_np))
    uAAu_i(1,1)=1./uAAu(1,1)
    
    
    
    
    
    !uA(1,1:numstep)=Au(1:numstep,1)
    
    !uAAu_i(1:num_M,1:num_M)=uAAu(1:num_M,1:num_M)
    ! 
    !call dgetrf(num_M,num_M,uAAu_i(1:num_M,1:num_M),num_M,ipiv(1:num_M),info)
     
    
    !do i=0,numstep-1
    !    call LegendrePoly(i*dt,0.,endt,M-1,basis)
    !    
    !    Mk(i,0:M-1)=basis(0:M-1)
    !    
    !end do
    
    Mu=0.
    uAMu=0.
    !uMAu=0.
    
    do i=0, M-1
        
        do j=0,np-1
        Mu(i,j+1,1)=Mk(1,i)*sum(B(j,0:np-1)*unum(np:2*np-1))
        end do
        
        do ns=2,numstep
            do j=0,np-1
            Mu(i,(ns-1)*np+j+1,1)=Mk(ns-1,i)*sum(B(j,0:np-1)*unum((ns-1)*np:ns*np-1))+mk(ns,i)*sum(B(j,0:np-1)*unum(ns*np:(ns+1)*np-1))
            end do
        end do
        
    !Mu(i,2:numstep,1)=Mk(1:numstep-1,i)*unum(1:numstep-1)
    
    uAMu(i,1,1)=sum(uA(1,1:numstep_np)*Mu(i,1:numstep_np,1))!sum(unum(2:numstep)**2*Mk(1:numstep-1,i))-(1+dt)*sum(unum(3:numstep)*unum(2:numstep-1)*Mk(1:numstep-2,i))
    !uMAu(i,1,1)=uAMu(i,1,1)
    
    end do
    
    
    res_Q=0.
    res_R=0.
    !res_RTR=0.
    
    
    
    tmp=sqrt(sum(uA(1,1:numstep_np)**2))
    res_Q(1:numstep_np,1)=uA(1,1:numstep_np)/tmp
    
    res_R(1,1)=tmp
    !res_RTR(1,1)=tmp*tmp
    nrank=1
    ncolumn=1
    do i=0,M-1
        
        tmpvec1(1:numstep_np)=Mu(i,1:numstep_np,num_M)
        res_R(nrank+1,ncolumn+1)=1.0
        do j=1,nrank
            
            tmp1=sum(tmpvec1(1:numstep_np)*res_Q(1:numstep_np,j))
            tmpvec2(1:numstep_np)=tmpvec1(1:numstep_np)-tmp1*res_Q(1:numstep_np,j)
    
            tmp2=sqrt(sum(tmpvec2(1:numstep_np)**2))
            
            res_R(j,ncolumn+1)=tmp1*res_R(nrank+1,ncolumn+1)
            
            
            if (tmp2 < 1.0E-17) then
                
                nrankup=0
                res_R(nrank+1,ncolumn+1)=0.
                exit
                
            else
                tmpvec1(1:numstep_np)=tmpvec2(1:numstep_np)/tmp2
                nrankup=1
                res_R(nrank+1,ncolumn+1)=tmp2*res_R(nrank+1,ncolumn+1)
                
            end if
            
        end do
        
        if (nrankup .eq. 1) then
            nrank=nrank+1
            
            res_Q(1:numstep_np,nrank)=tmpvec1(1:numstep_np)
            
            !do j=1,ncolumn
            !    res_RTR(ncolumn+1,j)=sum(res_R(1:nrank-1,ncolumn+1)*res_R(1:nrank-1,j))
            !    res_RTR(j,ncolumn+1)=res_RTR(ncolumn+1,j)
            !end do
            !
            !res_RTR(ncolumn+1,ncolumn+1)=sum(res_R(1:nrank,ncolumn+1)**2)
            
        !else
        !    
        !    do j=1,ncolumn
        !        res_RTR(ncolumn+1,j)=sum(res_R(1:nrank,ncolumn+1)*res_R(1:nrank,j))
        !        res_RTR(j,ncolumn+1)=res_RTR(ncolumn+1,j)
        !    end do
        !    res_RTR(ncolumn+1,ncolumn+1)=sum(res_R(1:nrank,ncolumn)**2)
            
        end if
        
        ncolumn=ncolumn+1
        
        
        
    end do
    
    
    !write(*,*) matmul(transpose(res_Q(1:numstep_np,1:nrank)),res_Q(1:numstep_np,1:nrank))
    
    !rq1=(1.-sum(res_Q(1,1:nrank)**2))**2
    !
    !do i=2,numstep_np
    !    rq1=rq1+sum(res_Q(i,1:nrank)*res_Q(1,1:nrank))**2
    !end do
    
                
                
    
    
    
    
    
    r_basis(1:numstep_np,0)=unum(np:numstep_np1)
    
    C_basis=0.
    C_basis(0,num_M-1)=rnorm
    
    basisnum(num_M)=0
    
    !u_rbm(0,0:np-1)=unum(numstep*np:(numstep+1)*np-1)
    
    pickup=1
    pickup(0)=0
    
    !do pos=1,M
    !     call run_ori(pos)
    !
    !     call update_rbmatrix(pos)
    !     
    !     num_M=num_M+1
    !     C_basis(pos,num_M-1)=1.0
    !     
    !     r_basis(numstep,num_M-1)=unum(numstep)
    !     
    !     write(*,*) "bais number=", pos
    !     
    !     !write(12,*) pos
    !     write(*,*)'alpha='
    !     write(*,'(100I2)')(alpha(pos,i),i=0,M-1)
    !     
    !     basisnum(num_M)=pos
    !     
    !     pickup(pos)=0
    !end do
    
    
    
    deallocate(tmpvec1,tmpvec2)
    
     !write(*,*) maxval(u_rbm(0,1:numstep)-u1(0,1:numstep))
    
    end subroutine
    
    
    ! compute the errror between rbm solution with "truth" solution
    subroutine comperror
    use comdata
    
    open(219,file="rbm.dat")
    do i=0,np-1
        write(219,*) 2.0*pi/np*i,u_ori(0,(numstep-1)*np+i+1),sum(u_ori(0:npc-1,(numstep-1)*np+i+1)*u_ori(0:npc-1,(numstep-1)*np+i+1))
    end do
    close(219)
    
    !open(220,file="data.dat")
    !do i=1,numstep_np
    !write(220,'(70ES26.16)') (u_ori(j,i),j=0,NPC-1)
    !end do
    
    
    
    error2=0.
    u_rbm=0.
    do i=0,NPC-1
        do j=0,num_M-1
        u_rbm(i,0:np-1)=u_rbm(i,0:np-1)+C_basis(i,j)*r_basis((numstep-1)*np+1:numstep*np,j)
        end do
        
    end do
    
    !open(unit=22,file='my_numu.dat')
    !do i=0,NPC-1
    !    write(22,"(160000ES26.16)") (u_ori(i,j),j=1,numstep_np)
    !end do
    
    
    !
    !
    !!do i=0,NPC-1
    !!    error=max(error,maxval(abs(u1(i,1:numstep)-u_rbm(i,1:numstep))))
    !!end do
    !
    !
    !
    !!write(13,*) 'max error=',error
    !!write(*,*) 'max error=',error
    !
    do j=0,np-1
    error2=error2+(sum(u_rbm(1:NPC-1,j)*u_rbm(1:NPC-1,j))-sum(u_ori(1:npc-1,(numstep-1)*np+1+j)*u_ori(1:npc-1,(numstep-1)*np+1+j)))**2
    end do
    error2=sqrt(error2)
    write(*,*) 'E[1] error=',error2
    write(*,*)
    write(2018,*) 'E[1] error=',error2
    write(2018,*)
    
    !
    !write(2018,*) 'E[0]=',u_rbm(0),'E[1]=',apem2
    !write(2018,*)
    !
    !
    !em2=sum(u_num(1:NPC-1)*u_num(1:NPC-1))
    !!em2=22.3982744180594         
    !error2=abs(apem2-em2)
    !write(*,*) 'E[0]_num=',u_rbm(0),'E[1]_num=',em2
    !write(*,*)
    !write(*,*) 'error0_truth=', abs(u_rbm(0)-u_rbm(0)),'errror1_truth=',error2
    !write(*,*)
    !
    !write(2018,*) 'E[0]_num=',u_num(0),'E[1]_num=',em2
    !write(2018,*)
    !write(2018,*) 'error0_truth=', abs(u_rbm(0)-u_num(0)),'errror1_truth=',abs(apem2-em2)
    !write(2018,*)
    !
    !em1=2.0*exp(1.0)-1.0
    !em2=(7.0*exp(3.0)-6.0*exp(1.0)+2.0)/3.0
    !em2=em2-em1*em1
    !write(*,*) 'em[0]=',em1,'em[1]=',em2
    !write(*,*)
    !write(*,*) 'error0_exact=', abs(u_rbm(0)-em1),'errror1_exact=',abs(apem2-em2)
    !
    !
    !write(2018,*) 'em[0]=',em1,'em[1]=',em2
    !write(2018,*)
    !write(2018,*) 'error0_exact=', abs(u_rbm(0)-em1),'errror1_exact=',abs(apem2-em2)
    !
    
    

!        write(17,*) 1.,0.,0
!        ns=1
!        do i=1,K
!            nend=ns+nchoosek(M+i-1,i)
!            do pos=ns,nend-1
!                if(pickup(pos) .eq. 0) then
!                write(17,*) real(pos-ns+1),real(i),pickup(pos)
!                end if
!                
!            end do
!            ns=nend
!        end do
!        
!    
!    write(15,"('$(',I2,',',I2,')$','&',ES12.2,'&',I4,'&',I4,'&',ES12.2,'\\')") K,M,epsilon,num_M,NPC,error
    
    end subroutine