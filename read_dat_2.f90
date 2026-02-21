!
! Copyright (C) 2019- 
!
! This program is free software; you can redistribute it and/or modify 
! it under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2 of the License, or 
! any later version. See the file LICENSE in the root directory of the 
! present distribution, or http://www.gnu.org/copyleft/gpl.txt ,
! or contact developers via e-mail: cono.dipaola@gmail.com 
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! Author: Cono Di Paola
!-----------------------------------------------------------------------
module nabla
    use prec
    use params
    implicit none
!  local
    integer :: ik,itemp,ierror1,nbands1,adv ! the last one is the advance of ik
    integer :: ibands,jbands,k1bands,k2bands
    integer :: ikpoints,ispin
    integer :: icart,jcart,iss  ! index ispin
    real(q) :: deltae

!  global
    complex(q), allocatable :: psi_nabla_psi(:,:,:,:) ! <psi|nabla|psi> matrix elements
    complex(q), allocatable :: psi_nabla_psi_2(:,:,:,:) ! <psi|nabla|psi>*(0,-i) matrix elements
    real(q), allocatable :: en_2(:,:)  ! eigenvalues
    real(q), allocatable :: fw_2(:,:)
    real(q), allocatable :: en(:,:)  ! eigenvalues
    real(q), allocatable :: fw(:,:)
    real(q), allocatable :: kx(:),ky(:),kz(:),wk(:)
    logical :: sind  ! indicator whether we have spin polarized case
contains
!====================================================
subroutine read_nabla()
    IU2=7
    if (.not.n_bin) then
       open(unit=IU2,file=FNAME)
    else if (n_bin) then
       open(unit=IU2,file=FNAME,form='unformatted',access='sequential')
    endif
!
    allocate(kx(nkpoints),ky(nkpoints),kz(nkpoints))
    allocate(en(nbands,nkpoints),fw(nbands,nkpoints))
    allocate(wk(nkpoints))
    allocate(en_2(nbands,nkpoints),fw_2(nbands,nkpoints))  ! my Odyssey starts here
    en_2=0
    fw_2=0
    allocate(psi_nabla_psi(1:3,nbands,nbands,nkpoints),psi_nabla_psi_2(1:3,nbands,nbands,nkpoints))
    psi_nabla_psi_2=0
!    ik=1
    adv=0
    itemp=1
    sind=.false.
!    write(30,*) FNAME
!    DO WHILE (itemp == 1)
    DO ik=1,nkpoints
    ! may be risky if nkpoints loaded from INPUT file does not match that from nabla.dat
            do iss=1,niss
!
        if (.not.n_bin) then
           READ (IU2,*,IOSTAT=ierror1) kx(ik),ky(ik),kz(ik),wk(ik)
        else if (n_bin) then
           READ (IU2,IOSTAT=ierror1) kx(ik),ky(ik),kz(ik),wk(ik)
        endif
        !write(*,*) 'position 1, ik=',ik
!
        IF (ierror1 > 0) THEN
            write(*,*) 'nabla data file: something wrong!'
            STOP
        ELSE IF (ierror1 < 0) THEN
            write(*,*) 'nabla data file: end of file reached'
            EXIT
        ELSE
!
!            write(60,*) kx(ik),ky(ik),kz(ik),wk(ik)
!            write(30,*) ik,itemp
!
                write(*,*) 'position 1, iss=',iss
            DO icart=1,3
                !write(*,*) 'position 2, icart=',icart
                if (.not.n_bin) then
                   READ (IU2,*,IOSTAT=ierror1) ikpoints,ispin,jcart,nbands1
                   write(*,*) 'position 2.4, ikpoints=',ikpoints,'ispin=',ispin,'jcart=',jcart,'nbands1=',nbands1
                else if (n_bin) then
                   READ (IU2,IOSTAT=ierror1) ikpoints,ispin,jcart,nbands1
                   write(*,*) 'position 2.5, ikpoints=',ikpoints,'ispin=',ispin,'jcart=',jcart,'nbands1=',nbands1
                endif
!
!                WRITE(60,*) ikpoints,ispin,icart,jcart,nbands1
                IF (nbands1 /= nbands) THEN
                    write(*,*) 'Number of bands read from INPUT (',nbands &
                      ,' ) differs from the number read in data file (,',nbands1,' )'
                    STOP
                ENDIF
                ! add spin judgement here and then split the cases for
                ! ispin=1 and ispin=2
                if (ispin==1) then
                        do ibands=1,nbands
                            do jbands=1,nbands
        !                        write(*,*) 'position 3, ib,jb=',ibands,jbands
                                if (.not.n_bin) then
                                   READ (IU2,*,IOSTAT=ierror1) k1bands,en(k1bands,ik),fw(k1bands,ik)&
                                       ,k2bands,en(k2bands,ik),fw(k2bands,ik)
                               else if (n_bin) then
                                   READ (IU2,IOSTAT=ierror1) k1bands,en(k1bands,ik),fw(k1bands,ik)&
                                       ,k2bands,en(k2bands,ik),fw(k2bands,ik)
                               endif
        !
        !                        WRITE(60,*) k1bands,en(k1bands,ik),fw(k1bands,ik)&
        !                            ,k2bands,en(k2bands,ik),fw(k2bands,ik)
                                if (k1bands == k2bands) then
                                  deltae=dabs(en(k1bands,ik)-en(k2bands,ik))
                                    if (deltae > de) then
                                        write(*,*) 'intra-band energies with too different eigenvalues'
                                      STOP 
                                 endif
                                endif
                                if (.not.n_bin) then
                                   READ (IU2,*,IOSTAT=ierror1) psi_nabla_psi(icart,k1bands,k2bands,ik)
                                else if (n_bin) then
                                   READ (IU2,IOSTAT=ierror1) psi_nabla_psi(icart,k1bands,k2bands,ik)
                                endif
        !                        WRITE(60,*)  psi_nabla_psi_i(icart,k1bands,k2bands,ik)
        !                        WRITE(60,*)  psi_nabla_psi(icart,k1bands,k2bands,ik)
                            enddo  !do jbands
                        enddo  !do ibands
                 else if (ispin==2) then
                        sind=.TRUE.
                        do ibands=1,nbands
                            do jbands=1,nbands
        !                        write(*,*) 'position 3, ib,jb=',ibands,jbands
                                if (.not.n_bin) then
                                   READ (IU2,*,IOSTAT=ierror1) k1bands,en_2(k1bands,ik),fw_2(k1bands,ik)&
                                       ,k2bands,en_2(k2bands,ik),fw_2(k2bands,ik)
                               else if (n_bin) then
                                   READ (IU2,IOSTAT=ierror1) k1bands,en_2(k1bands,ik),fw_2(k1bands,ik)&
                                       ,k2bands,en_2(k2bands,ik),fw_2(k2bands,ik)
                               endif
        !
        !                        WRITE(60,*) k1bands,en(k1bands,ik),fw(k1bands,ik)&
        !                            ,k2bands,en(k2bands,ik),fw(k2bands,ik)
                                if (k1bands == k2bands) then
                                  deltae=dabs(en_2(k1bands,ik)-en_2(k2bands,ik))
                                    if (deltae > de) then
                                        write(*,*) 'intra-band energies with too different eigenvalues'
                                      STOP 
                                 endif
                                endif
                                if (.not.n_bin) then
                                   READ (IU2,*,IOSTAT=ierror1) psi_nabla_psi_2(icart,k1bands,k2bands,ik)
                                else if (n_bin) then
                                   READ (IU2,IOSTAT=ierror1) psi_nabla_psi_2(icart,k1bands,k2bands,ik)
                                endif
        !                        WRITE(60,*)  psi_nabla_psi_i(icart,k1bands,k2bands,ik)
        !                        WRITE(60,*)  psi_nabla_psi(icart,k1bands,k2bands,ik)
                            enddo  !do jbands
                        enddo  !do ibands
                 endif
!                write(80,*) icart,ikpoints,ik
            ENDDO ! cartesian components  DO icart=1,3
!
!            write(70,*) icart,ikpoints,ik,ibands,jbands
        ENDIF
            enddo !iss=1,niss
!            ik=ik+1
!
    ENDDO  !ik
!
!  active conversion only if matrix elements coming from quantum espresso
    if (.not.FROM_VASP) then
    !    en=en * RYDTOEV  !convert eigenvalues from eV to Ha
        en=en * AUTOEV
        psi_nabla_psi=psi_nabla_psi / BRTOANG
    endif
!
!    fw = fw /2.0_q
    !
!    write(20,*) ikpoints,ispin,icart,nbands
!    write(20,*) ik,nkpoints
!   
    CLOSE(IU2)
!
return
end subroutine read_nabla
!
end module

