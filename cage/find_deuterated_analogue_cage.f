        program find_deuterated_analogue_cage
        implicit real*8(a-h,o-z)
        parameter (nwavetot=50)
        parameter (nmax=200000)
        parameter (ndim=54)
        parameter (natoms=18)
        parameter (noxy=6)
        parameter (nbondtot=12)
        dimension n(nwavetot),n0(nwavetot),time(nwavetot),weight(nmax),
     1  psips(ndim,nmax),coord(3,natoms,nmax),n1(nwavetot),ibad(nmax),
     1  coord_water(3,3,noxy,nmax),tot_weight(nwavetot),icount(nmax),
     1  roh(2,noxy,noxy,nmax),roh_prim(2,noxy,nmax),itype(noxy,nmax),
     1  roh_sec(2,noxy,noxy,nmax),ndonate(2,noxy,noxy,nmax),nflip(nmax),
     1  naccept(2,noxy,noxy,nmax),roh_sec_min(2,noxy,nmax),nflip2(nmax),
     1  loc_roh_sec(2,noxy,nmax),isame(noxy,nmax),deut(noxy,nmax),
     1  roh_use(nbondtot,nmax),loc_from(nbondtot,nmax),icheck(nmax),
     1  loc_to(nbondtot,nmax),wt_zero(nmax),izero(nwavetot),
     1  psips_zero(ndim,nmax,nwavetot),iwave_check(nmax,nwavetot),
     1  id_cage(nmax),wt_cage(nmax),nfrom_cage(nmax),ndonor2(nmax,nmax),
     1  psips_cage(ndim,nmax),loc_from_cage(nbondtot,nmax),
     1  loc_to_cage(nbondtot,nmax),iwalk_check(nmax,nwavetot),
     1  roh_sec_cage(nbondtot,nmax),ndonor(noxy,nmax),wt_cage2(nmax),
     1  nacceptor(noxy,nmax),marker(noxy,nmax),id_cage2(nmax),
     1  nfrom_cage2(nmax),psips_cage2(ndim,nmax),marker2(noxy,nmax),
     1  loc_from_cage2(nbondtot,nmax),loc_to_cage2(nbondtot,nmax),
     1  roh_sec_cage2(nbondtot,nmax),nacceptor2(noxy,nmax),
     1  wt_final(nmax),itype_final(noxy,nmax),psips_check(ndim,nmax),
     1  iwave(nmax),iwalk(nmax),zero(nwavetot)
        open(unit=7,file='../../../../testing1.dat',status='old',
     1  form='unformatted')
        open(unit=8,file='../../../../coord1.dat',status='old')
        open(unit=10,file='unopt_walkers_h2o5_d2o_cage.dat',
     1  status='unknown')
C		Identify location of the deuterated water molecule (or un-deuterated
C		water molecule) in the water hexamer cage structure through 
C		understanding hydrogen bonding network. 
C		Hydrogen bond is defined as a donor to acceptor bond length of less
C		than 3 angstroms.
C		Inputs:
C		testing*.dat = wave functions from DMC calculaton 
C		Coordinates in units of bohrs
C		coord*.dat = descendant weights from DMC calculations
C		Outputs:
C		unopt_walkers_h2o5_d2o_cage = walkers that might need to optimized in order 
C		to classify better (coordinates in units of bohr)
        pi = dacos(-1.d0)
        read(7) nwave
        do k = 1,nwave
            tot_weight(k) = 0.d0
            izero(k) = 0
            zero(k) = 0.d0
        enddo
        ir = 0
        tot_type1 = 0.d0
        tot_type2 = 0.d0
        tot_type3 = 0.d0
        tot_type4 = 0.d0
        tot_type5 = 0.d0
        tot_type6 = 0.d0
        tot_zero = 0.d0
        iztot = 0
        do k = 1,nwave
            do j = 1,noxy
                deut(j,k) = 0.d0
            enddo
        enddo
        do k = 1,nwave
            read(7) n(k),n0(k),time(k)
            do i = 1,n(k)
                read(7) (psips(j,i),j=1,ndim)
            enddo
            read(8,*) n1(k)
            if (n(k).ne.n1(k)) then
                print *, 'n and n1 are not equal',n(k),n1(k)
                stop
            endif
            do i = 1,n1(k)
                read(8,*) weight(i)
                tot_weight(k) = tot_weight(k) + weight(i)
            enddo
c       convert to xyz
            do i = 1,n(k)
                ip = 0
                do j = 1,natoms
                    do l = 1,3
                        ip = ip + 1
                        coord(l,j,i) = psips(ip,i)
                    enddo
                enddo
            enddo
            do i = 1,n(k)
                it = 0
                do j = 1,noxy
                    do l = 1,3
                        it = it + 1
                        do m = 1,3
                            coord_water(m,l,j,i) = coord(m,it,i)
                        enddo
                    enddo
                enddo
            enddo
            do i = 1,n(k)
                do j = 1,noxy
                    itype(j,i) = 0
                enddo
            enddo
c           calculate all the d2o distances
            if (k.ge.14) then
                ir = ir + 1
                do i = 1,n(k)
                    do j = 1,noxy
                        do l = 1,noxy
                            roh(1,l,j,i) = sqrt(((coord_water(1,2,j,i)-
     1                      coord_water(1,1,l,i))**2)+
     1                      ((coord_water(2,2,j,i)-
     1                      coord_water(2,1,l,i))**2)+
     1                      ((coord_water(3,2,j,i)-
     1                      coord_water(3,1,l,i))**2))
                            roh(2,l,j,i) = sqrt(((coord_water(1,3,j,i)-
     1                      coord_water(1,1,l,i))**2)+
     1                      ((coord_water(2,3,j,i)-
     1                      coord_water(2,1,l,i))**2)+
     1                      ((coord_water(3,3,j,i)-
     1                      coord_water(3,1,l,i))**2))
                        enddo
                    enddo
                enddo
                do i = 1,n(k)
                    do j = 1,noxy
                        do l = 1,noxy
                            do m = 1,2
                                ndonate(m,l,j,i) = j
                                naccept(m,l,j,i) = l
                                roh_sec(m,l,j,i) = 9999.d0
                            enddo
                        enddo
                    enddo
                enddo
                do i = 1,n(k)
                    do j = 1,noxy
                        isame(j,i) = 0
                    enddo
                enddo
                do i = 1,n(k)
                    do j = 1,noxy
                        do l = 1,noxy
                            if (j.eq.l) then
                                do m = 1,2
                                    roh_prim(m,j,i) = roh(m,l,j,i)
                                enddo
                            else
                                do m = 1,2
                                    roh_sec(m,l,j,i) = roh(m,l,j,i)
                                enddo
                            endif
                        enddo
                    enddo
                enddo
                do i = 1,n(k)
                    do j = 1,2
                        do l = 1,noxy
                            roh_sec_min(j,l,i) = 
     1                      minval(roh_sec(j,:,l,i),1)
                            loc_roh_sec(j,l,i) = 
     1                      minloc(roh_sec(j,:,l,i),1)
                        enddo
                    enddo
                enddo
                do i = 1,n(k)
                    do j = 1,noxy
                        do l = 1,2
                            do m = 1,2
                                if ((loc_roh_sec(l,j,i).eq.
     1                          loc_roh_sec(m,j,i)).and.(l.ne.m)) then
                                    isame(j,i) = 1
                                    exit
                                endif
                            enddo
                        enddo
                    enddo
                enddo
c               count up all the donor acceptor pairs
                do i = 1,n(k)
                    icount(i) = 0
                    do j = 1,noxy
                        if (isame(j,i).ne.1) then
                            do l =1,2
                                if (roh_sec_min(l,j,i)*0.529177.lt.3.d0)
     1                           then
                                    icount(i) = icount(i) + 1
                                    roh_use(icount(i),i)=
     1                              roh_sec_min(l,j,i)
                                    loc_from(icount(i),i) = j
                                    loc_to(icount(i),i)=
     1                              loc_roh_sec(l,j,i)
                                endif
                            enddo
                        else if (isame(j,i).eq.1) then
                            icount(i) = icount(i) + 1
                            roh_use(icount(i),i) = 
     1                      minval(roh_sec_min([1,2],j,i),1)
                            loc_from(icount(i),i) = j
                            loc_to(icount(i),i) = loc_roh_sec(1,j,i)
                        endif
                    enddo
                enddo
                do i = 1,n(k)
                    ibad(i) = 0
                enddo
                do i = 1,n(k)
                    do j = 1,icount(i)
                        do l = j,icount(i)
                            if ((loc_from(j,i).eq.loc_to(l,i)).and.
     1                      (loc_from(l,i).eq.loc_to(j,i))) then
                                ibad(i) = 1
                                exit
                            endif
                        enddo
                    enddo
                enddo
                iz = 0
                icage = 0
                do i = 1,n(k)
                    if ((icount(i).gt.8).or.(ibad(i).eq.1)) then
                        iz = iz + 1
                        wt_zero(iz) = weight(i) 
                        izero(ir) = izero(ir) + 1
                        do j = 1,ndim
                            psips_zero(j,izero(ir),ir) = psips(j,i)
                        enddo
                        iwave_check(izero(ir),ir) = k
                        iwalk_check(izero(ir),ir) = i
                    else
                        icage = icage + 1
                        id_cage(icage) = icount(i)
                        wt_cage(icage) = weight(i)
                        nfrom_cage(icage) = i
                        do j = 1,ndim
                            psips_cage(j,icage) = psips(j,i)
                        enddo
                        do j = 1,icount(i)
                            loc_from_cage(j,icage) = loc_from(j,i)
                            loc_to_cage(j,icage)=loc_to(j,i)
                            roh_sec_cage(j,icage) = roh_use(j,i)
                        enddo
                    endif
                enddo
                do i = 1,icage
                    do j = 1,noxy
                        ndonor(j,i) = 0
                        nacceptor(j,i) = 0
                    enddo
                enddo
                do i = 1,icage
                    do j = 1,id_cage(i)
                        do l = 1,noxy
                            if (loc_from_cage(j,i).eq.l) then
                                ndonor(l,i) = ndonor(l,i) + 1
                            endif
                            if (loc_to_cage(j,i).eq.l) then
                                nacceptor(l,i) = nacceptor(l,i) + 1
                            endif
                        enddo
                    enddo
                enddo
c           find the single donor single acceptor
                do i = 1,icage
                    nflip(i) = 0
                enddo
                do i = 1,icage  
                    do j = 1,noxy
                        if ((ndonor(j,i).eq.1).and.(nacceptor(j,i)
     1                  .eq.1)) then
                            nflip(i) = nflip(i) + 1
                            marker(nflip(i),i) = j
                        endif
                    enddo
                enddo
                icage2 = 0
                do i = 1,icage 
                    if (nflip(i).ne.2) then
                        iz = iz + 1
                        wt_zero(iz) = wt_cage(i) 
                        izero(ir) = izero(ir) + 1
                        do j = 1,ndim
                            psips_zero(j,izero(ir),ir) = psips_cage(j,i)
                        enddo
                        iwave_check(izero(ir),ir) = k
                        iwalk_check(izero(ir),ir) = nfrom_cage(i)
                    else
                        icage2 = icage2 + 1
                        nflip2(icage2) = nflip(i)
                        do j = 1,nflip(i)
                            marker2(j,icage2) = marker(j,i)
                        enddo
                        id_cage2(icage2) = id_cage(i)
                        wt_cage2(icage2) = wt_cage(i)
                        nfrom_cage2(icage2) = nfrom_cage(i)
                        do j = 1,ndim
                            psips_cage2(j,icage2) = psips_cage(j,i)
                        enddo
                        do j = 1,id_cage(i)
                            loc_from_cage2(j,icage2)=loc_from_cage(j,i)
                            loc_to_cage2(j,icage2)=loc_to_cage(j,i)
                            roh_sec_cage2(j,icage2) = roh_sec_cage(j,i)
                        enddo
                        do j = 1,noxy
                            ndonor2(j,icage2) = ndonor(j,i)
                            nacceptor2(j,icage2) = nacceptor(j,i)
                        enddo
                    endif
                enddo
c               mark flipping one as 1 and bond to double acceptor
                do i = 1,icage2
                    do j = 1,nflip2(i)
                        do l = 1,id_cage2(i)
                            if (marker2(j,i).eq.loc_from_cage2(l,i)) 
     1                      then
                                if (nacceptor2(loc_to_cage2(l,i),i).eq.
     1                          2) then
                                    itype(marker2(j,i),i) = 1
                                    itype(loc_to_cage2(l,i),i) = 2
                                endif
                            endif
                        enddo
                    enddo
                enddo
c               find the double acceptor to double donor bond
                do i = 1,icage2
                    do j = 1,noxy
                        do l = 1,id_cage2(i)
                            if (itype(j,i).eq.2) then
                                if (loc_from_cage2(l,i).eq.j) then
                                    itype(loc_to_cage2(l,i),i) = 3
                                endif
                            endif
                        enddo
                    enddo
                enddo
c               find the first double donor that connects to other
c               flip
                do i = 1,icage2
                    do j = 1,id_cage2(i)
                        do l = 1,noxy
                            if ((itype(l,i).eq.3).and.
     1                      (loc_from_cage2(j,i).eq.l)) then
                                if ((ndonor2(loc_to_cage2(j,i),i).eq.1)
     1                          .and.(nacceptor2(loc_to_cage2(j,i),i)
     1                          .eq.1))then
                                    itype(loc_to_cage2(j,i),i) = 4
                                endif
                            exit
                            endif
                        enddo
                    enddo
                enddo
c               finding the double acceptor from the flip
                do i = 1,icage2
                    do j = 1,noxy
                        do l = 1,id_cage2(i)
                            if (itype(j,i).eq.4) then
                                if (loc_from_cage2(l,i).eq.j) then
                                    itype(loc_to_cage2(l,i),i) = 5
                                endif
                            endif
                        enddo
                    enddo
                enddo
c               find the last one which should be a double donor
                do i = 1,icage2
                    do j = 1,id_cage2(i)
                        do l = 1,noxy
                            if ((itype(l,i).eq.5).and.
     1                      (loc_from_cage2(j,i).eq.l)) then
                                if ((ndonor2(loc_to_cage2(j,i),i).eq.1)
     1                          .and.(nacceptor2(loc_to_cage2(j,i),i)
     1                          .eq.2).and.
     1                          (itype(loc_to_cage2(j,i),i).eq.0)) then
                                    itype(loc_to_cage2(j,i),i) = 6
                                endif
                            exit
                            endif
                        enddo
                    enddo
                enddo
c               check last one connects with the first
                do i = 1,icage2
                    do j = 1,noxy
                        do l = 1,id_cage2(i)
                            if (itype(j,i).eq.6) then
                                if (j.eq.loc_from_cage2(l,i)) then
                                    icheck(i) = loc_to_cage2(l,i)
                                endif
                            endif
                        enddo
                    enddo
                enddo
                nfinal = 0
                do i = 1,icage2
                    do j = 1,noxy
                        if (itype(j,i).eq.1) then
                            if (icheck(i).eq.j) then
                                nfinal = nfinal + 1
                                wt_final(nfinal) = wt_cage2(i)
                                do l = 1,noxy
                                    itype_final(l,nfinal) = itype(l,i)
                                enddo
                            else
                                iz = iz + 1
                                wt_zero(iz) = wt_cage2(i) 
                                izero(ir) = izero(ir) + 1
                                do l = 1,ndim
                                    psips_zero(l,izero(ir),ir) = 
     1                              psips_cage2(l,i)
                                enddo
                                iwave_check(izero(ir),ir) = k
                                iwalk_check(izero(ir),ir) = 
     1                          nfrom_cage2(i)
                            endif
                        endif
                    enddo
                enddo
                do i = 1,nfinal
                    if (itype_final(6,i).eq.1) then
                         deut(1,ir) = deut(1,ir) + (wt_final(i)/
     1                   tot_weight(k))
                     else if (itype_final(6,i).eq.2) then
                         deut(2,ir) = deut(2,ir) + (wt_final(i)/
     1                   tot_weight(k))
                     else if (itype_final(6,i).eq.3) then
                         deut(3,ir) = deut(3,ir) + (wt_final(i)/
     1                   tot_weight(k))
                     else if (itype_final(6,i).eq.4) then
                         deut(4,ir) = deut(4,ir) + (wt_final(i)/
     1                   tot_weight(k))
                     else if (itype_final(6,i).eq.5) then
                         deut(5,ir) = deut(5,ir) + (wt_final(i)/
     1                   tot_weight(k))
                     else if (itype_final(6,i).eq.6) then
                         deut(6,ir) = deut(6,ir) + (wt_final(i)/
     1                   tot_weight(k))
                     endif
                enddo
                iztot = iztot + iz
                do i = 1,iz
                    zero(ir) = zero(ir) + wt_zero(i)/tot_weight(k)
                enddo
                tot_type1 = tot_type1 + deut(1,ir)
                tot_type2 = tot_type2 + deut(2,ir)
                tot_type3 = tot_type3 + deut(3,ir)
                tot_type4 = tot_type4 + deut(4,ir)
                tot_type5 = tot_type5 + deut(5,ir)
                tot_type6 = tot_type6 + deut(6,ir)
                tot_zero = tot_zero + zero(ir)
            endif
        enddo
        dev_type1 = 0.d0
        dev_type2 = 0.d0
        dev_type3 = 0.d0
        dev_type4 = 0.d0
        dev_type5 = 0.d0
        dev_type6 = 0.d0
        dev_zero = 0.d0
        avg_type1 = tot_type1/dfloat(ir)
        avg_type2 = tot_type2/dfloat(ir)
        avg_type3 = tot_type3/dfloat(ir)
        avg_type4 = tot_type4/dfloat(ir)
        avg_type5 = tot_type5/dfloat(ir)
        avg_type6 = tot_type6/dfloat(ir)
        avg_zero = tot_zero/dfloat(ir)
        do k = 1,ir
            dev_type1 = dev_type1 + ((deut(1,k)-avg_type1)**2)
            dev_type2 = dev_type2 + ((deut(2,k)-avg_type2)**2)
            dev_type3 = dev_type3 + ((deut(3,k)-avg_type3)**2)
            dev_type4 = dev_type4 + ((deut(4,k)-avg_type4)**2)
            dev_type5 = dev_type5 + ((deut(5,k)-avg_type5)**2)
            dev_type6 = dev_type6 + ((deut(6,k)-avg_type6)**2)
            dev_zero = dev_zero + ((zero(k)-avg_zero)**2)
        enddo
        it = 0
        do k = 1,ir
            do i = 1,izero(k)
                it = it + 1
                do j = 1,ndim
                    psips_check(j,it) = psips_zero(j,i,k)
                enddo
                iwave(it) = iwave_check(i,k)
                iwalk(it) = iwalk_check(i,k)
            enddo
        enddo
        write(10,*) it
        do i = 1,it
            write(10,*) (psips_check(j,i),j=1,ndim),iwave(i),iwalk(i)
        enddo
        dev_type1 = sqrt((dev_type1)/dfloat(ir))
        dev_type2 = sqrt((dev_type2)/dfloat(ir))
        dev_type3 = sqrt((dev_type3)/dfloat(ir))
        dev_type4 = sqrt((dev_type4)/dfloat(ir))
        dev_type5 = sqrt((dev_type5)/dfloat(ir))
        dev_type6 = sqrt((dev_type6)/dfloat(ir))
        dev_zero = sqrt((dev_zero)/dfloat(ir))
        print *, avg_type1,'+/-',dev_type1,'type 1'
        print *, avg_type2,'+/-',dev_type2,'type 2'
        print *, avg_type3,'+/-',dev_type3,'type 3'
        print *, avg_type4,'+/-',dev_type4,'type 4'
        print *, avg_type5,'+/-',dev_type5,'type 5'
        print *, avg_type6,'+/-',dev_type6,'type 6'
        print *, avg_zero,'+/-',dev_zero,'zero'
        end program
