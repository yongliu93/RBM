     
      module comdata
      implicit none
      
      type vector
          integer , allocatable, dimension (:) :: index
      end type vector
      
      real, allocatable, dimension (:,:) :: u
      integer, allocatable, dimension(:,:) :: alpha,Ic
      type (vector), allocatable, dimension(:) :: s_eps,s_alpha
      
      type(vector), allocatable, dimension(:) :: I_cindex
      
      real, allocatable, dimension(:) ::u_num
      
      !real, allocatable, dimension (:,:) :: c
      
      real :: tn,endt,dt
      
      integer :: NPC,K,M,numstep
      
      
      real, allocatable,dimension(:) :: f_LHS,unum
      real, allocatable,dimension(:,:) :: r_Basis,C_basis
      real, allocatable,dimension(:) :: u_rbm,u0
      
      real :: epsilon
      integer :: num_M
      integer , allocatable, dimension(:) :: basisnum,pickup
      !integer, allocatable, dimension(:):: sq
      real, allocatable,dimension(:,:) :: uAAu,uA,uAAu_i,AAu
      
      real, allocatable,dimension(:) :: frb
      
      real , allocatable, dimension(:,:) :: Mk
      real, allocatable, dimension(:,:,:) :: Mu,uAMu!,uMAu
      
      real, allocatable,dimension(:,:) :: res_Q,res_R!,res_RTR
      
      integer :: maxbasisnum,nrank,ncolumn
      
      integer , allocatable,dimension(:):: ipiv
      
      integer :: ns,nend
      
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
      
      numstep=10000
      
       do i=0,0
           epsilon=5.0E-4!1.0E-04*(1-i/2.)
           
        
        ! number of polynomials order and random variables 
        do K=8,11
            do M=K,K+1
        !K=3
        !M=3
        
        write(*,*) "(N,K)=",K,M
        write(2018,*) "(N,K)=",K,M
        !maxbasisnum=
           
           
      call assignmemory
      call init
      
      endt=1.0
      
      !call CPU_TIME(t1)
      !call run(numstep)
      !call CPU_TIME(t2)
      
      write(*,*) 'time no rbm=',t2-t1
      write(2018,*) 'time no rbm=',t2-t1
      
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
      !call outpdata
      call deletememory
      
      pause
            end do
        end do
        
      
       end do
       
      
    end
      
      
     
    
    
        subroutine assignmemory
        use comdata
        integer :: nchoosek
        
        
        
        NPC=nchoosek(K+M,K)
        maxbasisnum=30
        !maxbasisnum=NPC

        
        allocate(Ic(0:NPC-1,0:K-1),alpha(0:NPC-1,0:M-1))
        allocate(s_eps(0:NPC-1),s_alpha(0:NPC-1))
        allocate(u_num(0:NPC-1))
        
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
     allocate(r_Basis(1:numstep,0:maxbasisnum-1))
     allocate(f_LHS(1:numstep),unum(0:numstep))
     
     allocate(u_rbm(0:NPC-1),u0(0:NPC-1))
     
     allocate(basisnum(1:maxbasisnum))
     allocate(pickup(0:NPC-1))
     
     allocate(uAAu(1:maxbasisnum,1:maxbasisnum),uA(1:maxbasisnum,1:numstep),AAu(1:numstep,1:maxbasisnum))
     allocate(uAMu(0:M-1,1:maxbasisnum,1:maxbasisnum),Mu(0:M-1,1:numstep,1:maxbasisnum))
     allocate(uAAu_i(1:maxbasisnum,1:maxbasisnum))
     
     allocate(frb(1:maxbasisnum))
     
     allocate(Mk(0:numstep-1,0:M-1))
     
     allocate(ipiv(1:maxbasisnum))
     
     allocate(res_Q(1:numstep,1:maxbasisnum*(M+1)),res_R(1:maxbasisnum*(M+1),1:maxbasisnum*(M+1)))
     
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
          
          do i=2,k
              p(i)=((2.0*i-1.0)*x*p(i-1)-(i-1.0)*p(i-2))/i
          end do
          
          do i=0,k
              p(i)=p(i)*sqrt((2.0*i+1.0)/(right-left))
          end do
          
      end if

      end subroutine LegendrePoly
      
      
      
      subroutine getDiscretization(res)
      use comdata
      
      real, dimension (0:M-1) :: basis
      integer :: i
      
      real, intent(out),dimension (0:NPC-1) ::res
      
      res=0.
      res(0)=u_num(0)+1.0
      
      call LegendrePoly(tn,0.,endt,M-1,basis)
      
      do i=1,NPC-1
          res(i)=u_num(i)
          
          do l=1,s_eps(i)%index(0)
          res(i)=res(i)+sqrt(real(alpha(i,s_eps(i)%index(l))))* &
          u_num(s_alpha(i)%index(l))*basis(s_eps(i)%index(l))
          end do
          
      end do
      
      
      
      
      end subroutine getDiscretization
      
      
      subroutine init
      use comdata
      integer :: nchoosek,i,j,pos
      integer ,allocatable,dimension(:) :: alpha_tmp
      
      u_num=0.
      u_num(0)=1.0
      
      tn=0.
      !allocate(I_cindex(0:nchoosek(8,4)-1))
      
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
      
      !do i=0,NPC-1
      !    call IToPos(Ic(i,:),M,K,pos)
      !    write(*,*) pos
      !end do
      !
      !do i=0,NPC-1
      !    write(*,"(4I8)") alpha(i,:)
      !    !write(*,"(4I8)") Ic(i,:)
      !    write(*,*)
      !end do
      
      
      call buildSkorohod
      
      
      
      
      
      !do i=0,70-1
      !    write(*,"(4I8)") I_cindex(i)%index(:)
      !end do
      
      end subroutine init
      
      
      subroutine stepForward
      use comdata
      
      !real, dimension (0:NPC-1) :: u0,k1,k2,k3,k4
      real, dimension(0:NPC-1) :: k1
      !u0=u_num
      k1=0.0
      
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
      !u_num=u0+0.5*k2*dt
      !
      !call getDiscretization(k3)
      !
      !u_num=u0+k3*dt
      !
      !tn=tn+0.5*dt
      !
      !call getDiscretization(k4)
      !
      !u_num=u0+(0.5*k1+k2+k3+0.5*k4)*dt/3.0
      
      
      
      
      
      end subroutine stepForward
      
      
      
      subroutine run(Nstep)
      use comdata
      
      
      integer :: Nstep,step
      dt =endt/Nstep
      
      do step=1,Nstep
          
          call stepForward
          
          !write(*,*) 'Step=',step,'t=',tn
          
          !u1(:,step)=u_num(:)
          !write(*,*) sum(abs(u_num(:)-u(:,step)))/NPC
      end do
      
      
      end subroutine run
      
      
      
      subroutine deletememory
      use comdata
      
      
      deallocate(u_num,Ic,alpha,s_eps,s_alpha,I_cindex)
     
      
      deallocate(f_LHS,r_Basis,C_basis,u_rbm,unum,basisnum,u0,frb)
      deallocate(pickup)
      deallocate(uAAu,Mu,uA,uAMu,AAu)
      deallocate(uAAu_i)
      deallocate(ipiv,Mk)
      
      deallocate(res_Q,res_R)
      end subroutine
      
      
      subroutine reduce_basis
      use comdata
      integer :: pos
      
      dt =endt/numstep
      
      call init_rbm
      
      !call assemble_mat
      
      
      !ns=1
      
      !do i=1,K
      !    
      !    call pickupbais(i)
      !    
      !end do
      
      
      do pos=M+1,NPC-1
          
          !call assemble_lhs(pos)
          
          call argminmize_coefficient(pos)
          
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
    
    
    subroutine assemble_lhs(pos)
    use comdata
    
    integer :: pos
    !real, dimension (0:M-1) :: basis
     
     
     
     f_LHS=0.
    !tn=(0)*dt
        !call LegendrePoly(tn,0.,endt,M-1,basis)
        do l=1,s_eps(pos)%index(0)
          f_LHS(1)=f_LHS(1)+sqrt(real(alpha(pos,s_eps(pos)%index(l))))* &
         u0(s_alpha(pos)%index(l))*Mk(0,s_eps(pos)%index(l))
        end do
        f_LHS(1)=dt*f_LHS(1)+u0(pos)+dt*u0(pos)
    
        

        
    do i=2,numstep
        !tn=(i-1)*dt
        !call LegendrePoly(tn,0.,endt,M-1,basis)
        do l=1,s_eps(pos)%index(0)
          f_LHS(i)=f_LHS(i)+sqrt(real(alpha(pos,s_eps(pos)%index(l))))* &
         sum(Mu(s_eps(pos)%index(l),i,1:num_M)*C_basis(s_alpha(pos)%index(l),0:num_M-1))
        end do
        f_LHS(i)=dt*f_LHS(i)
        
    end do
    
    
    end subroutine
    
   
    subroutine assemble_frb(pos)
    use comdata
    
    integer :: pos
    
    !real, allocatable, dimension(:) :: tmp
    
    !allocate(tmp(1:num_M))
    
    frb(1:num_M)=uA(1:num_M,1)*(u0(pos)+dt*u0(pos))
    
    do l=1,s_eps(pos)%index(0)
        
        frb(1:num_M)=frb(1:num_M)+uA(1:num_M,1)*dt*sqrt(real(alpha(pos,s_eps(pos)%index(l))))*Mk(0,s_eps(pos)%index(l))*u0(s_alpha(pos)%index(l))
        do i=1,num_M
          frb(i)=frb(i)+dt*sqrt(real(alpha(pos,s_eps(pos)%index(l))))* &
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
     
     if (resrb<epsilon) then
     
      !if(pickup(pos) .eq. 1) then 
          
         C_basis(pos,0:num_M-1)=tmp_c(1:num_M,1)
         !write(*,*) resrb
         
      else
     
         !write(*,*) resrb
         call run_ori(pos)

         call update_rbmatrix(pos)
         
         num_M=num_M+1
         C_basis(pos,num_M-1)=1.0
         
         r_basis(numstep,num_M-1)=unum(numstep)
         
         write(*,*) "bais number=", pos
         
         write(12,*) pos
         write(*,*)'alpha='
         write(*,'(100I2)')(alpha(pos,i),i=0,M-1)
         
         basisnum(num_M)=pos
         
         pickup(pos)=0
         
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
            
            tmp2(j)=tmp2(j)-dt*sqrt(real(alpha(pos,s_eps(pos)%index(l))))*C_basis(s_alpha(pos)%index(l),i)*res_R(j,i*(M+1)+s_eps(pos)%index(l)+2)
        
        end do
        
    end do
    
    !tmp2(j)=tmp
    
    end do
    
    !tmp=0.
    !do i=0,num_M-1
    !    tmp=tmp+c_alpha(i+1)*tmp2(i*(M+1)+1)
    !    do l=1,s_eps(pos)%index(0)
    !        
    !        tmp=tmp+C_basis(s_alpha(pos)%index(l),i)*tmp2(i*(M+1)+s_alpha(pos)%index(l)+2)
    !    end do
    !    
    !end do
    
    
    resrb=sqrt(sum(tmp2(1:nrank)**2))
    
    
    
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
    
    allocate(tmpvec1(1:numstep),tmpvec2(1:numstep))
    
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
    
    uA(num_M+1,1)=unum(1)
    uA(num_M+1,2:numstep)=unum(2:numstep)-(1+dt)*unum(1:numstep-1)
    
    AAu(1:numstep-1,num_M+1)=uA(num_M+1,1:numstep-1)-(1+dt)*uA(num_M+1,2:numstep)
    AAu(numstep,num_M+1)=uA(num_M+1,numstep)
    
    do i=1,num_M
    uAAu(num_M+1,i)=sum(uA(num_M+1,1:numstep)*uA(i,1:numstep))
    uAAu(i,num_M+1)=uAAu(num_M+1,i)
    end do
    
    uAAu(num_M+1,num_M+1)=sum(uA(num_M+1,1:numstep)*uA(num_M+1,1:numstep))
    
    !uA(num_M+1,1:numstep)=Au(1:numstep,num_M+1)
    
    
    do i=0, M-1
    Mu(i,2:numstep,num_M+1)=Mk(1:numstep-1,i)*unum(1:numstep-1)
    
    do j=1,num_M
    uAMu(i,num_M+1,j)=sum(uA(num_M+1,1:numstep)*Mu(i,1:numstep,j))!sum(unum(2:numstep)*Mk(1:numstep-1,i)*r_basis(2:numstep,j-1))-(1+dt)*sum(r_basis(3:numstep,j-1)*unum(2:numstep-1)*Mk(1:numstep-2,i))
    !uMAu(i,j,num_M+1)=uAMu(i,num_M+1,j)
    uAMu(i,j,num_M+1)=sum(uA(j,1:numstep)*Mu(i,1:numstep,num_M+1))!sum(unum(2:numstep)*Mk(1:numstep-1,i)*r_basis(2:numstep,j-1))-(1+dt)*sum(r_basis(2:numstep-1,j-1)*unum(3:numstep)*Mk(1:numstep-2,i))
    !uMAu(i,num_M+1,j)=uAMu(i,j,num_M+1)
    end do
    
    uAMu(i,num_M+1,num_M+1)=sum(uA(num_M+1,1:numstep)*Mu(i,1:numstep,num_M+1))!sum(unum(2:numstep)*Mk(1:numstep-1,i)*unum(2:numstep))-(1+dt)*sum(unum(3:numstep)*unum(2:numstep-1)*Mk(1:numstep-2,i))
    
    
    end do
    
    
    
    
    
    
    do i=-1,M-1
        
        if (i .eq. -1) then
            tmpvec1(1:numstep)=uA(num_M+1,1:numstep)
        else
            
        tmpvec1(1:numstep)=Mu(i,1:numstep,num_M+1)
        end if
        
        res_R(nrank+1,ncolumn+1)=1.0
        do j=1,nrank
            
            tmp1=sum(tmpvec1(1:numstep)*res_Q(1:numstep,j))
            tmpvec2(1:numstep)=tmpvec1(1:numstep)-tmp1*res_Q(1:numstep,j)
    
            tmp2=sqrt(sum(tmpvec2(1:numstep)**2))
            
            res_R(j,ncolumn+1)=tmp1*res_R(nrank+1,ncolumn+1)
            
            
            if (tmp2 < 1.0E-17) then
                
                nrankup=0
                res_R(nrank+1,ncolumn+1)=0.
                exit
                
            else
                tmpvec1(1:numstep)=tmpvec2(1:numstep)/tmp2
                nrankup=1
                res_R(nrank+1,ncolumn+1)=tmp2*res_R(nrank+1,ncolumn+1)
                
            end if
            
        end do
        
        if (nrankup .eq. 1) then
            nrank=nrank+1
            
            res_Q(1:numstep,nrank)=tmpvec1(1:numstep)
            
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
    
    
    
     tmp(1:num_M+1,1:num_M+1)=uAAu(1:num_M+1,1:num_M+1)
     
     call dgetrf(num_M+1,num_M+1,tmp(1:num_M+1,1:num_M+1),num_M+1,ipiv(1:num_M+1),info)
     
     call dgetri(num_M+1,tmp,num_M+1,ipiv(1:num_M+1),work,4*(num_M+1),info)
     
     
     uAAu_i(1:num_M+1,1:num_M+1)=tmp(1:num_M+1,1:num_M+1)
     
     deallocate(tmp,tmpvec1,tmpvec2)
     
    end subroutine
    
    
    ! "truth" solution solver
    subroutine run_ori(pos)
    use comdata
    
    
      real, dimension (0:M-1) :: basis
      integer :: i,pos
      

      if (pos .eq. 0) then
      !res=unum(0)+1.0
      
      do i=1,numstep
          unum(i)=unum(i-1)+dt*(unum(i-1)+1.0)
      end do
      
      !u1(pos,1:numstep)=unum(1:numstep)
      !u_rbm(pos,1:numstep)=unum(1:numstep)
      
      else
          
          call assemble_lhs(pos)
          
          unum(1)=f_LHS(1)
          do i=2,numstep
          unum(i)=f_LHS(i)+(1.+dt)*unum(i-1)
          end do
          
          !u1(pos,1:numstep)=unum(1:numstep)
          !u_rbm(pos,1:numstep)=unum(1:numstep)
          !write(*,*) maxval(abs(u_rbm(pos,1:numstep)-u1(pos,1:numstep)))
      end if
      
      
      
    end subroutine
    
    
    ! init reduced basis method
    subroutine init_rbm
    use comdata
    
    real, dimension(0:M-1) :: basis
    
    real, allocatable, dimension(:) :: tmpvec1,tmpvec2
    real :: tmp1,tmp2
    integer :: pos
    
    allocate(tmpvec1(1:numstep),tmpvec2(1:numstep))
    !epsilon=1.0E-04
    
    u0=0.
    unum=0.
    u0(0)=1.0
    unum(0)=1.0
    
    !call LegendrePoly(0.,0.,endt,M-1,basis0)
    
    
    num_M=1
    r_basis=0.
    
    call run_ori(0)
    
    !---------construct Arb
    
    !uAAu(1,1)=(1.+(1+dt)**2)*sum(unum(1:numstep-1)*unum(1:numstep-1))+unum(numstep)*unum(numstep)&
    !    -2*(1.+dt)*sum(unum(1:numstep-1)*unum(2:numstep))
    !
    !uAAu_i(1,1)=1./uAAu(1,1)
    
    AAu(1,1)=(1+(1+dt)**2)*unum(1)-(1+dt)*unum(2)
    AAu(2:numstep-1,1)=(1+(1+dt)**2)*unum(2:numstep-1)-(1+dt)*unum(3:numstep)-(1+dt)*unum(1:numstep-2)
    AAu(numstep,1)=unum(numstep)-(1+dt)*unum(numstep-1)
    
    !uAA(1,1:numstep)=AAu(1:numstep,1)
    
    uA(1,1)=unum(1)
    uA(1,2:numstep)=unum(2:numstep)-(1+dt)*unum(1:numstep-1)
    uAAu(1,1)=sum(uA(1,1:numstep)*uA(1,1:numstep))
    uAAu_i(1,1)=1./uAAu(1,1)
    !uA(1,1:numstep)=Au(1:numstep,1)
    
    !uAAu_i(1:num_M,1:num_M)=uAAu(1:num_M,1:num_M)
    ! 
    !call dgetrf(num_M,num_M,uAAu_i(1:num_M,1:num_M),num_M,ipiv(1:num_M),info)
     
    
    do i=0,numstep-1
        call LegendrePoly(i*dt,0.,endt,M-1,basis)
        
        Mk(i,0:M-1)=basis(0:M-1)
        
    end do
    
    Mu=0.
    uAMu=0.
    !uMAu=0.
    
    do i=0, M-1
    Mu(i,2:numstep,1)=Mk(1:numstep-1,i)*unum(1:numstep-1)
    
    uAMu(i,1,1)=sum(uA(1,1:numstep)*Mu(i,1:numstep,1))!sum(unum(2:numstep)**2*Mk(1:numstep-1,i))-(1+dt)*sum(unum(3:numstep)*unum(2:numstep-1)*Mk(1:numstep-2,i))
    !uMAu(i,1,1)=uAMu(i,1,1)
    
    end do
    
    
    res_Q=0.
    res_R=0.
    !res_RTR=0.
    
    
    
    tmp=sqrt(sum(uA(1,1:numstep)**2))
    res_Q(1:numstep,1)=uA(1,1:numstep)/tmp
    
    res_R(1,1)=tmp
    !res_RTR(1,1)=tmp*tmp
    nrank=1
    ncolumn=1
    do i=0,M-1
        
        tmpvec1(1:numstep)=Mu(i,1:numstep,num_M)
        res_R(nrank+1,ncolumn+1)=1.0
        do j=1,nrank
            
            tmp1=sum(tmpvec1(1:numstep)*res_Q(1:numstep,j))
            tmpvec2(1:numstep)=tmpvec1(1:numstep)-tmp1*res_Q(1:numstep,j)
    
            tmp2=sqrt(sum(tmpvec2(1:numstep)**2))
            
            res_R(j,ncolumn+1)=tmp1*res_R(nrank+1,ncolumn+1)
            
            
            if (tmp2 < 1.0E-17) then
                
                nrankup=0
                res_R(nrank+1,ncolumn+1)=0.
                exit
                
            else
                tmpvec1(1:numstep)=tmpvec2(1:numstep)/tmp2
                nrankup=1
                res_R(nrank+1,ncolumn+1)=tmp2*res_R(nrank+1,ncolumn+1)
                
            end if
            
        end do
        
        if (nrankup .eq. 1) then
            nrank=nrank+1
            
            res_Q(1:numstep,nrank)=tmpvec1(1:numstep)
            
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
    
    
    !write(*,*) matmul(transpose(res_Q(1:numstep,1:nrank)),res_Q(1:numstep,1:nrank))
    
                
                
    
    
    
    
    
    r_basis(1:numstep,0)=unum(1:numstep)
    
    C_basis=0.
    C_basis(0,num_M-1)=1.
    
    basisnum(num_M)=0
    
    u_rbm(0)=unum(numstep)
    
    pickup=1
    pickup(0)=0
    
    do pos=1,M
         call run_ori(pos)

         call update_rbmatrix(pos)
         
         num_M=num_M+1
         C_basis(pos,num_M-1)=1.0
         
         r_basis(numstep,num_M-1)=unum(numstep)
         
         write(*,*) "bais number=", pos
         
         !write(12,*) pos
         write(*,*)'alpha='
         write(*,'(100I2)')(alpha(pos,i),i=0,M-1)
         
         basisnum(num_M)=pos
         
         pickup(pos)=0
    end do
    
    
    
    deallocate(tmpvec1,tmpvec2)
     !write(*,*) maxval(u_rbm(0,1:numstep)-u1(0,1:numstep))
    
    end subroutine
    
    
    ! compute the errror between rbm solution with "truth" solution
    subroutine comperror
    use comdata
    
    error=0.
    u_rbm=0.
    do i=0,NPC-1
        do j=0,num_M-1
        u_rbm(i)=u_rbm(i)+C_basis(i,j)*r_basis(numstep,j)
        end do
        
    end do
    
    
    !do i=0,NPC-1
    !    error=max(error,maxval(abs(u1(i,1:numstep)-u_rbm(i,1:numstep))))
    !end do
    
    
    
    !write(13,*) 'max error=',error
    !write(*,*) 'max error=',error
    
    apem2=sum(u_rbm(1:NPC-1)*u_rbm(1:NPC-1))
    write(*,*) 'E[0]=',u_rbm(0),'E[1]=',apem2
    write(*,*)
    
    write(2018,*) 'E[0]=',u_rbm(0),'E[1]=',apem2
    write(2018,*)
    
    
    em2=sum(u_num(1:NPC-1)*u_num(1:NPC-1))
    write(*,*) 'E[0]_num=',u_num(0),'E[1]_num=',em2
    write(*,*)
    write(*,*) 'error0_truth=', abs(u_rbm(0)-u_num(0)),'errror1_truth=',abs(apem2-em2)
    write(*,*)
    
    write(2018,*) 'E[0]_num=',u_num(0),'E[1]_num=',em2
    write(2018,*)
    write(2018,*) 'error0_truth=', abs(u_rbm(0)-u_num(0)),'errror1_truth=',abs(apem2-em2)
    write(2018,*)
    
    em1=2.0*exp(1.0)-1.0
    em2=(7.0*exp(3.0)-6.0*exp(1.0)+2.0)/3.0
    em2=em2-em1*em1
    write(*,*) 'em[0]=',em1,'em[1]=',em2
    write(*,*)
    write(*,*) 'error0_exact=', abs(u_rbm(0)-em1),'errror1_exact=',abs(apem2-em2)
    
    
    write(2018,*) 'em[0]=',em1,'em[1]=',em2
    write(2018,*)
    write(2018,*) 'error0_exact=', abs(u_rbm(0)-em1),'errror1_exact=',abs(apem2-em2)
    
    
    

        write(17,*) 1.,0.,0
        ns=1
        do i=1,K
            nend=ns+nchoosek(M+i-1,i)
            do pos=ns,nend-1
                if(pickup(pos) .eq. 0) then
                write(17,*) real(pos-ns+1),real(i),pickup(pos)
                end if
                
            end do
            ns=nend
        end do
        
    
    write(15,"('$(',I2,',',I2,')$','&',ES12.2,'&',I4,'&',I4,'&',ES12.2,'\\')") K,M,epsilon,num_M,NPC,error
    
    end subroutine