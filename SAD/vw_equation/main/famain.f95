!c===============================c
        include 'table.mdu'
!c===============================c
        program lesann
!c-------------------------------c
        call ipd0fm
        call ipd1uf
        call ipdff0
        call augrid
	
        write(*,*)'*--------------------------------------------------*'
        write(*,*)'*                                                  *'
        write(*,*)'*    The program is to calculate 2D-rans           *'
        write(*,*)'*                                                  *'
        write(*,*)'*--------------------------------------------------*'
        write(*,*)'*                                                  *'
        call begins
        call nseqts
!        call update
!        call initia
        write(*,*)'*                                                  *'
        write(*,*)'*--------------------------------------------------*'
        write(*,*)'*                                                  *'
        write(*,*)'*         output the divergence-free field         *'
!c        call opd0fm
!c        call opd1fm
!c        call opdff0
!c        call opdhf0
	call shuchu
        write(*,*)'*                                                  *'
        write(*,*)'*         The----End----Of---The---Program         *'
        write(*,*)'*                                                  *'
        write(*,*)'*--------------------------------------------------*'
        stop
        end

!c======================================================================c
        subroutine nseqts 
!c======================================================================c
        include 'table.prc'
        write(*,*)'*--------------------------------------------------*'
        write(*,*)'*                                                  *'
        write(*,*)'*        solving  process for momentum Eqs.        *'

        write(*,*)'*                                                  *'
        call leinuo
        write(*,*)'*                                                  *'
        do k=1,10000
            write(*,*)'*                                                  *'
            write(*,611) k
            write(*,*)'*                                                  *'
            write(*,*)'*--------------------------------------------------*'
            write(*,*)'*        solving  momentum  equation  of  v        *'
            call initiv
            call equatv

            write(*,*)'*                                                  *'
            write(*,*)'*        solving  momentum  equation  of  w        *'
            call initiw
            call equatw
            write(*,*)'*                                                  *'
            write(*,*)'*--------end of handling momentum equations--------*'

            write(*,*)'*--------------------------------------------------*'
            write(*,*)'*                                                  *'
            write(*,*)'*        solving pressure-correct  equation        *'
            call initip
            call equatp
            call update
            write(*,*)'*                                                  *'
            write(*,*)'*--------------------------------------------------*'
        enddo
         
!        call equatu
611     format(' *',8x,'compute steps = ',i5,21x,'*')
	close(91) 

        return
        end

!c======================================================================c
        subroutine initia
!c======================================================================c
        call initiv
        call initiw
        return
        end
