!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
!  
!       Contributions to this routine:
!                         Daniela di Serafino    Second University of Naples
!                         Pasqua D'Ambra         ICAR-CNR
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
!   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   C                                                                      C
!   C  References:                                                         C
!   C          [1] Duff, I., Marrone, M., Radicati, G., and Vittoli, C.    C
!   C              Level 3 basic linear algebra subprograms for sparse     C
!   C              matrices: a user level interface                        C
!   C              ACM Trans. Math. Softw., 23(3), 379-401, 1997.          C
!   C                                                                      C
!   C                                                                      C
!   C         [2]  S. Filippone, M. Colajanni                              C
!   C              PSBLAS: A library for parallel linear algebra           C
!   C              computation on sparse matrices                          C
!   C              ACM Trans. on Math. Softw., 26(4), 527-550, Dec. 2000.  C
!   C                                                                      C
!   C         [3] M. Arioli, I. Duff, M. Ruiz                              C
!   C             Stopping criteria for iterative solvers                  C
!   C             SIAM J. Matrix Anal. Appl., Vol. 13, pp. 138-144, 1992   C
!   C                                                                      C
!   C                                                                      C
!   C         [4] R. Barrett et al                                         C
!   C             Templates for the solution of linear systems             C
!   C             SIAM, 1993                                               C
!   C                                                                      C
!   C                                                                      C
!   C         [5] G. Sleijpen, D. Fokkema                                  C
!   C             BICGSTAB(L) for linear equations involving unsymmetric   C
!   C             matrices with complex spectrum                           C
!   C             Electronic Trans. on Numer. Analysis, Vol. 1, pp. 11-32, C
!   C             Sep. 1993                                                C
!   C                                                                      C
!   C                                                                      C
!   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! File:  psb_dsgmres.f90
!
! Subroutine: psb_dsgmres
!    This subroutine implements the sketched restarted GMRES method with right
!    preconditioning.
!
! Arguments:
!
!    a      -  type(psb_dspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_dprec_type)       Input: preconditioner
!    b      -  real,dimension(:)       Input: vector containing the
!                                         right hand side B
!    x      -  real,dimension(:)       Input/Output: vector containing the
!                                         initial guess and final solution X.
!    eps    -  real                       Input: Stopping tolerance; the iteration is
!                                         stopped when the error estimate |err| <= eps
!    desc_a -  type(psb_desc_type).       Input: The communication descriptor.
!    info   -  integer.                   Output: Return code
!
!    itmax  -  integer(optional)          Input: maximum number of iterations to be
!                                         performed.
!    iter   -  integer(optional)          Output: how many iterations have been
!                                         performed.
!                                         performed.
!    err    -  real   (optional)          Output: error estimate on exit. If the
!                                         denominator of the estimate is exactly
!                                         0, it is changed into 1. 
!    itrace -  integer(optional)          Input: print an informational message
!                                         with the error estimate every itrace
!                                         iterations
!    istop  -  integer(optional)          Input: stopping criterion, or how
!                                         to estimate the error. 
!                                         1: err =  |r|/(|a||x|+|b|);  here the iteration is
!                                            stopped when  |r| <= eps * (|a||x|+|b|)
!                                         2: err =  |r|/|b|; here the iteration is
!                                            stopped when  |r| <= eps * |b|
!                                         where r is the (preconditioned, recursive
!                                         estimate of) residual. 
!    irst   -  integer(optional)          Input: restart parameter 
!
subroutine psb_dsgmres_vect(a,prec,b,x,eps,desc_a,info,&
     & itmax,iter,err,itrace,irst,istop)
  use psb_base_mod
  use psb_prec_mod
  use psb_d_linsolve_conv_mod
  use psb_linsolve_mod
  use psb_util_mod
  implicit none
  type(psb_dspmat_type), intent(in)    :: a
  Type(psb_desc_type), Intent(in)      :: desc_a
  class(psb_dprec_type), intent(inout) :: prec
  type(psb_d_vect_type), Intent(inout) :: b
  type(psb_d_vect_type), Intent(inout) :: x
  type(psb_d_multivect_type) :: SK
  Real(psb_dpk_), Intent(in)           :: eps
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, irst,istop
  integer(psb_ipk_), Optional, Intent(out)       :: iter
  Real(psb_dpk_), Optional, Intent(out) :: err
! =   local data
  real(psb_dpk_), allocatable   :: aux(:)
  real(psb_dpk_), allocatable   :: c(:), s(:), h(:,:), rs(:), rst(:)
  real(psb_dpk_), allocatable   :: SKAV(:,:), Sb(:), SKAV2(:,:), Sb2(:)
  type(psb_d_vect_type), allocatable :: v(:)
  type(psb_d_vect_type)              :: w, w1, xt
  real(psb_dpk_) :: tmp 
  real(psb_dpk_) :: scal, gm, rti, rti1
  integer(psb_ipk_) ::litmax, naux, it, k, itrace_,&
       & n_row, n_col, nl, nsketch, korth
  integer(psb_lpk_) :: mglob
  Logical, Parameter :: exchange=.True., noexchange=.False., use_srot=.true.
  integer(psb_ipk_), Parameter :: irmax = 8
  integer(psb_ipk_) :: itx, i, istop_, err_act
  integer(psb_ipk_) :: debug_level, debug_unit
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me  
  Real(psb_dpk_)     :: rni, xni, bni, ani,bn2, dt, r0n2
  real(psb_dpk_)     :: errnum, errden, deps, derr, dnrm2
  character(len=20)           :: name
  character(len=*), parameter :: methdname='SGMRES'

  info = psb_success_
  name = 'psb_dsgmres'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt = desc_a%get_context()
  Call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_ext_)&
       & write(debug_unit,*) me,' ',trim(name),': from psb_info',np
  if (.not.allocated(b%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  mglob = desc_a%get_global_rows()
  n_row = desc_a%get_local_rows()
  n_col = desc_a%get_local_cols()

  if (present(istop)) then 
    istop_ = istop 
  else
    istop_ = 2
  endif
!
!  ISTOP_ = 1:  Normwise backward error, infinity norm 
!  ISTOP_ = 2:  ||r||/||b||, 2-norm 
!

  if ((istop_ < 1 ).or.(istop_ > 2 ) ) then
    info=psb_err_invalid_istop_
    err=info
    call psb_errpush(info,name,i_err=(/istop_/))
    goto 9999
  endif

  if (present(itmax)) then 
    litmax = itmax
  else
    litmax = 1000
  endif

  if (present(itrace)) then
    itrace_ = itrace
  else
    itrace_ = 0
  end if
  
  if (present(irst)) then
    nl = irst
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' present: irst: ',irst,nl
  else
    nl = 10 
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' not present: irst: ',irst,nl
  endif
  if (nl <=0 ) then 
    info=psb_err_invalid_irst_
    err=info
    call psb_errpush(info,name,i_err=(/nl/))
    goto 9999
  endif

  call psb_chkvect(mglob,lone,x%get_nrows(),lone,lone,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_chkvect on X')
    goto 9999
  end if
  call psb_chkvect(mglob,lone,b%get_nrows(),lone,lone,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_    
    call psb_errpush(info,name,a_err='psb_chkvect on B')
    goto 9999
  end if


  naux=4*n_col 
  allocate(aux(naux),h(nl+1,nl+1),&
       &c(nl+1),s(nl+1),rs(nl+1), rst(nl+1),stat=info)

  if (info == psb_success_) call psb_geall(v,desc_a,info,n=nl+1)
  if (info == psb_success_) call psb_geall(w,desc_a,info)
  if (info == psb_success_) call psb_geall(w1,desc_a,info)
  if (info == psb_success_) call psb_geall(xt,desc_a,info)
  if (info == psb_success_) call psb_geasb(v,desc_a,info,mold=x%v)  
  if (info == psb_success_) call psb_geasb(w,desc_a,info,mold=x%v)  
  if (info == psb_success_) call psb_geasb(w1,desc_a,info,mold=x%v)  
  if (info == psb_success_) call psb_geasb(xt,desc_a,info,mold=x%v)

  ! Sketching: preallocate a Rademacher matrix that we will use for sketching.
  ! To ensure an epsilon-embedding, we select it twice as large as the maximum 
  ! number of iterations before a restart.
  nsketch = min(2 * (nl + 1), mglob)
  korth = 1
  if (info == psb_success_) call psb_geall(SK,desc_a,info,n=nsketch)
  call psb_dsgmres_vect_gen_sketch(SK, desc_a)

  allocate(Sb(nsketch), SKAV(nsketch, nl), Sb2(nsketch), SKAV2(nsketch, nl), stat=info)

  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_non_ 
    call psb_errpush(info,name)
    goto 9999
  end if
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Size of V,W,W1 ',v(1)%get_nrows(),size(v),&
       & w%get_nrows(),w1%get_nrows()


  if (istop_ == 1) then 
    ani = psb_spnrmi(a,desc_a,info)
    bni = psb_geamax(b,desc_a,info)
  else if (istop_ == 2) then 
    bn2 = psb_genrm2(b,desc_a,info)
  else if (istop_ == 3) then
    call psb_geaxpby(done,b,dzero,v(1),desc_a,info)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    
    call psb_spmm(-done,a,x,done,v(1),desc_a,info,work=aux)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    r0n2 = psb_genrm2(v(1),desc_a,info)
  endif
  errnum = dzero
  errden = done
  deps   = eps
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_non_ 
    call psb_errpush(info,name)
    goto 9999
  end if
  if ((itrace_ > 0).and.(me == 0)) call log_header(methdname)

  itx   = 0
  restart: do 
  
    ! compute r0 = b-ax0
    ! check convergence

    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' restart: ',itx,it
    it = 0      
    call psb_geaxpby(done,b,dzero,v(1),desc_a,info)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_spmm(-done,a,x,done,v(1),desc_a,info,work=aux)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if

    ! Sketch r0, store ||Sro|| as the first residual norm
    Sb = psb_gedot(SK, v(1), desc_a, info)
    rs(1) = dnrm2(nsketch, Sb, ione)
    rs(2:) = dzero
    scal = done/rs(1)
    ! scal = done / psb_genrm2(v(1), desc_a, info)

    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' on entry to amax: b: ',b%get_nrows(),rs(1),scal

    !
    ! check convergence
    !
    if (istop_ == 1) then 
      rni = psb_geamax(v(1),desc_a,info)
      xni = psb_geamax(x,desc_a,info)
      errnum = rni
      errden = (ani*xni+bni)
    else if (istop_ == 2) then 
      rni = psb_genrm2(v(1),desc_a,info)
      errnum = rni
      errden = bn2
    else if (istop_ == 3) then 
      rni = psb_genrm2(v(1),desc_a,info)
      errnum = rni
      errden = r0n2
    endif
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    
    if ((errnum <= eps*errden).or.(itx >= litmax)) exit restart  

    if (itrace_ > 0) &
         & call log_conv(methdname,me,itx,itrace_,errnum,errden,deps)
     
    call v(1)%scal(scal) !v(1) = v(1) * scal

    !
    ! inner iterations
    !
    inner:  Do i=1,nl
      itx  = itx + 1

      call prec%apply(v(i),w1,desc_a,info)
      call psb_spmm(done,a,w1,dzero,w,desc_a,info,work=aux)
      
      ! Sketch the action of the operator
      SKAV(:, i) = psb_gedot(SK, w, desc_a, info)

      ! Only partial reorthogonalization is done in the sketched variant
      do k = max(1, i - korth), i
        h(k,i) = psb_gedot(v(k),w,desc_a,info)
        call psb_geaxpby(-h(k,i),v(k),done,w,desc_a,info)
      end do

      ! LR: In principle the scaling here is not really needed
      h(i+1,i) = psb_genrm2(w,desc_a,info)
      scal=done/h(i+1,i)
      call psb_geaxpby(scal,w,dzero,v(i+1),desc_a,info)

      ! Build the solution for the sketched least square problem |SKAV*rst - Sb|_2 
      ! We rely on dgels that overwrite input data with the solution, so we need to make 
      ! a copy of the matrix SKAV and the RHS.
      SKAV2(:, 1:i) = SKAV(:, 1:i)
      Sb2(:) = Sb(:)

      call dgels('N', nsketch, i, 1, SKAV2, nsketch, Sb2, nsketch, aux, naux, info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_non_ 
        call psb_errpush(info,name)
        goto 9999
      end if

      rs(1:i) = Sb2(1:i)
      rst(1:i) = rs(1:i)
      rni = dnrm2(nsketch - i, Sb2(i+1 : nsketch), ione)
      
      if (istop_ == 1) then
        !
        ! build x and then compute the residual and its infinity norm
        !
        call w1%set(dzero)
        do k=1, i
          call psb_geaxpby(rst(k),v(k),done,xt,desc_a,info)
        end do
        call prec%apply(xt,desc_a,info)
        call psb_geaxpby(done,x,done,xt,desc_a,info)
        call psb_geaxpby(done,b,dzero,w1,desc_a,info)
        call psb_spmm(-done,a,xt,done,w1,desc_a,info,work=aux)
        rni = psb_geamax(w1,desc_a,info)
        xni = psb_geamax(xt,desc_a,info)
        errnum = rni
        errden = (ani*xni+bni)
        !

      else if (istop_ == 2) then 
        !
        ! compute the residual sketched 2-norm as byproduct of the solution
        ! procedure of the least-squares problem
        !
        errnum = rni
        errden = bn2
      else if (istop_ == 3) then 
        !
        ! compute the residual sketched 2-norm as byproduct of the solution
        ! procedure of the least-squares problem
        !
        errnum = rni
        errden = r0n2
      endif

      if (errnum <= eps*errden) then 

        if (istop_ == 1) then 
          call psb_geaxpby(done,xt,dzero,x,desc_a,info)
! =          x = xt 
        else if (istop_ == 2) then
          !
          ! build x
          !
          call w1%set(dzero)
          do k=1, i
            call psb_geaxpby(rs(k),v(k),done,w1,desc_a,info)
          end do
          call prec%apply(w1,w,desc_a,info)
          call psb_geaxpby(done,w,done,x,desc_a,info)
        end if

        if (itrace_ > 0) &
             & call log_conv(methdname,me,itx,ione,errnum,errden,deps)
        exit restart

      end if

      if (itrace_ > 0) &
           & call log_conv(methdname,me,itx,itrace_,errnum,errden,deps)

    end do inner

    if (istop_ == 1) then 
      call psb_geaxpby(done,xt,dzero,x,desc_a,info)!      x = xt 
    else if (istop_ == 2) then
      !
      ! build x
      !
      call w1%set(dzero)
      do k=1, nl
        call psb_geaxpby(rs(k),v(k),done,w1,desc_a,info)
      end do
      call prec%apply(w1,w,desc_a,info)
      call psb_geaxpby(done,w,done,x,desc_a,info)
    end if
    
    if (itx >= litmax) then 
      if (itrace_ > 0) then 
        if (mod(itx,itrace_)/=0) &
             &  call log_conv(methdname,me,itx,ione,errnum,errden,deps)
      end if
      exit restart
    end if
    
  end do restart

  call log_end(methdname,me,itx,itrace_,errnum,errden,deps,err=derr,iter=iter)
  if (present(err)) err = derr

  
  if (info == psb_success_) call psb_gefree(v,desc_a,info)
  if (info == psb_success_) call psb_gefree(w,desc_a,info)
  if (info == psb_success_) call psb_gefree(w1,desc_a,info)
  if (info == psb_success_) call psb_gefree(xt,desc_a,info)
  if (info == psb_success_) deallocate(aux,h,c,s,rs,rst,Sb,SKAV,Sb2,SKAV2,stat=info)
  if (info == psb_success_) call psb_gefree(SK,desc_a,info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_dsgmres_vect

subroutine psb_dsgmres_vect_gen_sketch(SK, desc_a)
  use psb_base_mod
  use psb_prec_mod
  use psb_d_linsolve_conv_mod
  use psb_linsolve_mod
  implicit none

  type(psb_d_multivect_type), intent(InOut) :: SK
  type(psb_desc_type), intent(In) :: desc_a

  ! functions
  real(psb_dpk_) :: dsqrt
 
  ! local variables
  integer(psb_ipk_) :: err_act
  character(len=20) :: name
  integer(psb_ipk_) :: nlr, i, j, info, n, nsketch
  real(psb_dpk_), allocatable :: val(:,:)

  integer(psb_lpk_), allocatable :: myidx(:)

  name = 'psb_gen_sketch'

  myidx = desc_a%get_global_indices()
  nlr = size(myidx)
  n = size(SK%v%v, 1)
  nsketch = size(SK%v%v, 2)

  allocate(val(n, nsketch), stat=info)
  if (info /= psb_success_) then
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  do j = 1, nsketch
    do i = 1, n
      call random_number(val(i,j))
        if (val(i,j) .gt. 0.5) then
          val(i,j) = done
        else
          val(i,j) = -done
        end if
    end do
  end do

  val = val / dsqrt(real(nsketch, psb_dpk_))

  call psb_geins(nlr, myidx, val, SK, desc_a, info)
  if (info /= psb_success_) then
    info = psb_err_wrong_ins_
    call psb_errpush(info,name)
    goto 9999
  end if

  deallocate(val)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end
