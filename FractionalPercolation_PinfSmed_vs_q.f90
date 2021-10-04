!This algorithm computes the fraction of nodes in the GC for a fractional percolation process
!and the mean finite cluster size as a function of q

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Input
!N_node: total number of nodes
!kmin:   minimum degree or connectivity
!kmax:   maximum degree or connectivity
!lambda: mean degree for an ER network, or the exponent for an SF network
!nrea  : number of realizations
!rfix  : probability that a node is disfunctional, given that the node is not fully functional

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Output
!Pinf.dat: relative size of the GC as a function of q
!PinfFF.dat: fraction of FF nodes that belong to the GC as a function of q
!PinfPF.dat: fraction of PF nodes that belong to the GC as a function of q
!SmedT.dat: mean finite cluster size as a function of q

module globals
  implicit none
  save

  integer nptosr,nptosq
  integer kmin,kmax
  integer N_node  
  integer cluster_number,max_mass,gc

  real(8) lambda
  real(8) r  !random number
  real(8) rfix

  integer,allocatable,dimension(:)::edge
  integer,allocatable,dimension(:)::node
  integer,allocatable,dimension(:)::kk
  integer,allocatable,dimension(:)::cluster,mass
  integer,allocatable,dimension(:)::State
  
  real(8),allocatable,dimension(:)::Pk
  real(8),allocatable,dimension(:)::fracQ,fracR
  real(8),allocatable,dimension(:,:)::Pinf,S2,Smed1,SmedT,Smed0,PinfPF,PinfFF
    
end module globals

module random
   save
   integer::q1=1,q2=104
   integer ir(670)
   integer::sem(670)
   real(8)::nmax=2*(2**30-1)+1
   real(8) sig
end module random


Program FractionalPercolation
  use globals
  use random

  implicit none

  integer i,i_r,i_q,rea,nrea
  real(8) product
  real(8) pmin,deltp


  !Parameters
  N_node   =100000

  kmin     =0         
  kmax     =20
  lambda   =4d0
  nrea     =100 
  rfix     =0.5d0
  
  
        
  nptosr   =1
  nptosq   =100
  pmin     =0d0
  deltp    =0.01d0


  allocate(kk(N_node),node(N_node),State(N_node))
  allocate(Pk(kmin:kmax))
  allocate(cluster(N_node),mass(N_node))
  allocate(Pinf(nptosq,nptosr),SmedT(nptosq,nptosr))
  allocate(PinfPF(nptosq,nptosr),PinfFF(nptosq,nptosr))
  allocate(fracQ(nptosq),fracR(nptosr)) 
  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Degree distribution
  Pk           =0d0
  Pk(0)        =exp(-lambda)
  product      =1d0
  do i=kmin,kmax
     if(i==0)cycle
     product=product*i         
     Pk(i)=1.d0*exp(-lambda)*lambda**(i)/dble(product)
  enddo

  Pk=Pk/sum(Pk)  !normalizing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fracQ=0d0
  do i=1,nptosq
     fracQ(i)=pmin+i*deltp
  enddo
  
  fracR=rfix
  
    
  call initialize_random

  Pinf   =0d0
  PinfPF =0d0
  PinfFF =0d0
  SmedT  =0d0  

  do rea=1,nrea
     print*,rea
     call ConfModel
     call FracPercSubroutine
     deallocate(edge)
     
     open(1,file='Pinf.dat')
     write(1,*)'#',rea
     do i_q=1,nptosq
        do i_r=1,nptosr
           !write(1,*) fracQ(i_q),fracR(i_r),Pinf(i_q,i_r)/dble(rea)/dble(N_node)
           write(1,*) fracQ(i_q),Pinf(i_q,i_r)/dble(rea)/dble(N_node)
        enddo
     enddo
     close(1) 
     
     open(1,file='PinfFF.dat')
     write(1,*)'#',rea
     do i_q=1,nptosq
        do i_r=1,nptosr
           !write(1,*) fracQ(i_q),fracR(i_r),PinfFF(i_q,i_r)/dble(rea)/dble(N_node)
           write(1,*) fracQ(i_q),PinfFF(i_q,i_r)/dble(rea)/dble(N_node)
        enddo
     enddo
     close(1)     
     open(1,file='PinfPF.dat')
     write(1,*)'#',rea
     do i_q=1,nptosq
        do i_r=1,nptosr
           !write(1,*) fracQ(i_q),fracR(i_r),PinfPF(i_q,i_r)/dble(rea)/dble(N_node)
           write(1,*) fracQ(i_q),PinfPF(i_q,i_r)/dble(rea)/dble(N_node)
        enddo
     enddo
     close(1)
        
     open(1,file='SmedT.dat')
     write(1,*)'#',rea
     do i_q=1,nptosq
        do i_r=1,nptosr
           !write(1,*) fracQ(i_q),fracR(i_r),SmedT(i_q,i_r)/dble(rea)/dble(N_node)
           write(1,*) fracQ(i_q),SmedT(i_q,i_r)/dble(rea)/dble(N_node)
        enddo
     enddo
     close(1) 
     
  enddo

end Program FractionalPercolation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  Fractional Percolation subroutine
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine FracPercSubroutine
  use globals
  use random
  implicit none
  integer i,j,k
  integer i_q,i_r
  integer N_neto
  
  real(8) q,rr
  !  1-q is the probability that a node is fully functional (State=0)
  !  q*(1-rr) is the probability that a node is partially functional (State=1)
  !  q*rr is the probability that a node is dysfunctional (State=-1)
  
  do i_q=nptosq,1,-1       
     q   =fracQ(i_q)
     do i_r=nptosr,1,-1 
        rr=fracR(i_r)
        State=0
        do i=1,N_node
           call rand
           if(r<q)then
              call rand
              if(r<rr)then
                 State(i)=-1
              else
                 State(i)=1
              endif
           else
              State(i)=0
           endif
           
        enddo
        call cluster_id
        Pinf(i_q,i_r)=Pinf(i_q,i_r)+max_mass 
        do i=1,N_node
           if(cluster(i)==gc)then
              if(State(i)==0)PinfFF(i_q,i_r)=PinfFF(i_q,i_r)+1d0
              if(State(i)==1)PinfPF(i_q,i_r)=PinfPF(i_q,i_r)+1d0
           else
              if(State(i)/=-1)then
                 SmedT(i_q,i_r)=SmedT(i_q,i_r)+mass(cluster(i))
              endif              
           endif
        enddo             
     enddo
  enddo

end subroutine FracPercSubroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Configuration Model
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This subroutine creates three arrays: kk, edge, and node, which contains all the network topology information.
!kk(i) is the degree of node "i".
!The subset edge(node(i):node(i)+kk(i)-1) contains the list of neighbors of node "i".

subroutine ConfModel
  use globals
  use random
  implicit none
  integer i,j,k
  integer m,n,stubmaxpos,stubminpos
  integer counting, nstubs
  integer nstubAux,stubpos1,stubpos2

  real(8) ws

  integer, allocatable, dimension(:)::listStub,kkAux

  real(8),dimension(kmin-1:kmax)::CumulativePk
  
  
  allocate(kkAux(n_node))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Cumulative distribution
  
  CumulativePk   =0d0
  ws             =0d0
  do i=kmin,kmax
     ws                = ws+Pk(i)
     CumulativePk(i)   = ws
  enddo
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Assigning Connectivity to all nodes
2 kk=0 
  do i=1,N_node
     call rand
     do j=kmin,kmax
        if(CumulativePk(j-1)<r.and.r<=CumulativePk(j)) then
           kk(i)    = j
           exit
        endif
     enddo
  enddo

  nstubs=sum(kk)
  
  if(mod(nstubs,2).ne.0) go to 2   !the total number of stubs must be an even number
  
  allocate(listStub(nstubs))  
  allocate(edge(nstubs))

  listStub      =0
  nstubAux      =0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!List of stubs

  do i=1,n_Node
     if(kk(i)==0) cycle
     do k=1,kk(i)
        nstubAux           =nstubAux+1
        listStub(nstubAux) =i
     enddo
  enddo

  node(1)       =1
  do i=1,n_node-1
     node(i+1)  = node(i)+kk(i)
  enddo

  edge        =0
  kkAux       =0
  counting    =0
  do while(nstubs>0)
3    if(counting==100) then
	deallocate(listStub,edge)
	go to 2
     end if
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Randomly choose two stubs
     
     call rand
     stubpos1    =r*nstubs+1
     call rand
     stubpos2    =r*nstubs+1

     m        =listStub(stubpos1)   !the first stub corresponds to node "m"
     n        =listStub(stubpos2)   !the second stub corresponds to node "n"

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Checking if these two stubs can be connected

     if(m.ne.n) then
	do k=1,kkAux(m)
           if(edge(node(m)+k-1)==n)then
              counting     =counting+1
              go to 3   !if nodes "m" and "n" are already connected, we have to choose another pair of stubs
           endif
	enddo
	edge(node(m)+kkAux(m))   =n
	edge(node(n)+kkAux(n))   =m

	kkAux(m)                 =kkAux(m)+1
	kkAux(n)                 =kkAux(n)+1
	stubmaxpos               =max0(stubpos1,stubpos2)
	stubminpos               =min0(stubpos1,stubpos2)
        !updating the list of stubs
	if(stubmaxpos==nstubs) then
           listStub(stubmaxpos)       =listStub(nstubs)
           listStub(stubminpos)       =listStub(nstubs-1)
	else 
           if (stubmaxpos==nstubs-1) then
              listStub(stubminpos)    =listStub(nstubs)
           else
              listStub(stubmaxpos)    =listStub(nstubs)
              listStub(stubminpos)    =listStub(nstubs-1)
           endif
	endif
	nstubs           =nstubs-2
	counting         =0
     else
        counting           =counting+1
        go to 3 !self-loops are forbidden
     endif
  enddo
  deallocate(kkAux,listStub)
end subroutine ConfModel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Cluster identification
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cluster_id
use Globals
implicit none

integer m,i,j,k
integer nburn,cs

integer,dimension(N_node)::ocupp
integer,dimension(0:N_node+1)::w

ocupp=1
cluster=0
cluster_number=0
mass=0
max_mass=0
do i=1,N_node
   if(State(i)==-1) cycle !dysfunctional nodes do not belong to the GC
   if(cluster(i) .eq. 0.and.ocupp(i)==1) then 
      cluster_number=cluster_number+1
      w(0)=i
      nburn=1
      cs=0
      do while(nburn.ne.0)
         nburn=nburn-1
         j=w(nburn)
         cluster(j)=cluster_number
         cs=cs+1
         do k=0,kk(j)-1
            m=edge(node(j)+k)
            if(State(m)==-1)cycle  !dysfunctional neighbors do not belong to the GC
            if(State(j)==State(m).and.State(j)==1)cycle  !The PF nodes are not connected to each other
            if(cluster(m) == 0.and.ocupp(m)==1) then
               ocupp(m)=0  
               w(nburn)=m
               nburn=nburn+1
            endif
         enddo
      enddo
      mass(cluster_number)=cs
      if(cs > max_mass) then
         max_mass=cs
         gc=cluster_number
      endif
    endif
enddo
return
end subroutine cluster_id


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Seed for the Pseudorandom number generator
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_random
  use globals
  use random
  implicit none

  integer i,see(33)
  integer hour

  CALL SYSTEM_CLOCK(hour)
  see=hour
  CALL RANDOM_SEED(put=see)		
  CALL RANDOM_SEED(get=see)

  do i=1,670
     call random_number(r)
     sem(i)   =r*nmax
     ir(i)    =i+1
  enddo
  ir(670)=1
  return
end subroutine initialize_random
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Pseudo-random number generator
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rand
  use globals
  use random
  implicit none

 1 q1=ir(q1)
  q2=ir(q2)
  sem(q1)= IEOR(sem(q1),sem(q2))
  r=dfloat(sem(q1))/nmax
  if(r==1d0) go to 1
  return
end subroutine rand


