        program calc_roo
        implicit real*8(a-h,o-z)
        parameter (nwavetot=200)
        parameter (nmax=1000000)
        parameter (ndim=54)
        parameter (natoms=18)
        parameter (noxy=6)
        dimension n(nwavetot),n0(nwavetot),time(nwavetot),weight(nmax),
     1  psips(ndim,nmax),coord(3,natoms,nmax),n1(nwavetot),
     1  coord_water(3,3,noxy,nmax),tot_weight(nwavetot),id_prism(nmax),
     1  roh(2,noxy,noxy,nmax),idonor(nmax),roh_prim(2,noxy,nmax),
     1  roh_sec(2,noxy,noxy,nmax),ndonate(2,noxy,noxy,nmax),
     1  naccept(2,noxy,noxy,nmax),roh_sec_min(2,noxy,nmax),idon(nmax),
     1  loc_roh_sec(2,noxy,nmax),isame(noxy,nmax),loc_donor(noxy,nmax),
     1  loc_to(2,noxy,nmax),roh_donor_check(2,noxy,nmax),id(noxy,nmax),
     1  icount(noxy,nmax),loc_donor_prism(noxy,nmax),itype(noxy,nmax),
     1  loc_to_prism(2,noxy,nmax),icount_prism(noxy,nmax),ia(noxy,nmax),
     1  itype_prism(noxy,nmax),iterms(noxy,nmax),roh_keep(2,noxy,nmax),
     1  loc_roh_keep(2,noxy,nmax),ilower(nmax),itype_lower(noxy,nmax),
     1  roh_lower(2,noxy,nmax),loc_roh_lower(2,noxy,nmax),loc_min(nmax),
     1  loc_roh_lower_min(noxy,nmax),itop(2,noxy,nmax),roh_upper(2,nmax)
     1  ,loc_from_upper(2,nmax),loc_to_upper(2,nmax),roh_upper_min(nmax)
     1  ,loc_upper_min(nmax),roh_rest(2,noxy,nmax),icheck(2,noxy,nmax),
     1  loc_roh_upper(2,noxy,nmax),itype_upper(noxy,nmax),wt_lower(nmax)
     1  ,loc_in_upper(noxy,nmax),loc_from_upper2(noxy,nmax),iwalk(nmax),
     1  roh_check1(2,nmax),loc_to_check1(2,nmax),loc_from_check1(2,nmax)
     1  ,roh_check2(2,nmax),loc_to_check2(2,nmax),roh_check2_min(nmax),
     1  loc_from_check2(2,nmax),roh_check1_min(nmax),roh_check3(2,nmax),
     1  loc_roh_check1_min(nmax),loc_roh_check1(nmax),wt_upper(nmax),
     1  loc_roh_check2_min(nmax),loc_roh_check2(nmax),weight_final(nmax)
     1  ,loc_roh_check3_min(nmax),roh_check3_min(nmax),wt_prism(nmax),
     1  loc_roh_check3(nmax),loc_from_check3(2,nmax),iwave(nmax),
     1  loc_to_check3(2,nmax),deut(noxy,nwavetot),itype_final(noxy,nmax)
     1  ,psips_check(ndim,nmax),psips_zero(ndim,nmax,nwavetot),
     1  izero(nwavetot),wt_zero(nmax),iwave_check(nmax,nwavetot),
     1  iwalk_check(nmax,nwavetot),zero(nwavetot),nfrom_prism(nmax),
     1  psips_prism(ndim,nmax),nfrom_lower(nmax),psips_lower(ndim,nmax),
     1  nfrom_upper(nmax),psips_upper(ndim,nmax),psips_opt(ndim,nmax),
     1  iwave_bad(nmax),iwalk_bad(nmax),idonor2(nmax),iloc(noxy,nmax),
     1  loc_donor2(noxy,nmax),loc_to2(2,noxy,nmax),loc_bottom(nmax),
     1  roh_donor_check2(2,noxy,nmax),loc_from(12,nmax),loc_use(nmax),
     1  loc_acceptor(12,nmax),roh_sec_min2(2,noxy,nmax)
        open(unit=7,file='../../../../testing1.dat',status='old',
     1  form='unformatted')
        open(unit=8,file='../../../../coord1.dat',status='old')
        open(unit=11,file='unopt_walkers_h2o5_d2o_prism.dat',
     1  status='unknown')
C		Identify location of the deuterated water molecule (or un-deuterated
C		water molecule) in the water hexamer prism structure through 
C		understanding hydrogen bonding network. 
C		Hydrogen bond is defined as a donor to acceptor bond length of less
C		than 3 angstroms.
C		Inputs:
C		testing*.dat = wave functions from DMC calculaton 
C		Coordinates in units of bohrs
C		coord*.dat = descendant weights from DMC calculations
C		Outputs:
C		unopt_walkers_h2o5_d2o_prism = walkers that might need to optimized in order 
C		to classify better (coordinates in units of bohr)
        pi = dacos(-1.d0)
        read(7) nwave
        do k = 1,nwave
            tot_weight(k) = 0.d0
        enddo
        ir = 0
        tot_type1 = 0.d0
        tot_type2 = 0.d0
        tot_type3 = 0.d0
        tot_type4 = 0.d0
        tot_type5 = 0.d0
        tot_type6 = 0.d0
        tot_zero = 0.d0
        do k = 1,nwave
            do j = 1,noxy
                deut(j,k) = 0.d0
            enddo
            zero(k) = 0.d0
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
                    idonor(i) = 0
                enddo
c               find the double donors
                do i = 1,n(k)
                    do j = 1,noxy
                        if (loc_roh_sec(1,j,i).eq.
     1                  loc_roh_sec(2,j,i)) then
                            isame(j,i) = 1
                        endif
                    enddo
                enddo
                do i = 1,n(k)
                    do j = 1,noxy
                        if (all((roh_sec_min([1,2],j,i)*0.529177).lt.
     1                  3.0d0).and.(isame(j,i).ne.1)) then
                            idonor(i) = idonor(i) + 1
                            loc_donor(idonor(i),i) = j
                            do l = 1,2
                                loc_to(l,idonor(i),i)=loc_roh_sec(l,j,i)
                                roh_donor_check(l,idonor(i),i) = 
     1                          roh_sec_min(l,j,i)
                            enddo
                        endif
                    enddo
                enddo
c               determine where the double donor connecting a double
c               donor and double acceptor and label as icount(j,i) = 1
                do i = 1,n(k)
                    do j = 1,noxy
                        icount(j,i) = 0
                    enddo
                enddo
                do i = 1,n(k)
                    do j = 1,idonor(i)
                        do l = 1,2
                            do m = 1,idonor(i)
                                if (loc_donor(m,i).eq.loc_to(l,j,i)) 
     1                          then
                                    icount(j,i) = 1
                                    itype(loc_donor(j,i),i) = 1
                                endif
                            enddo
                        enddo
                    enddo
                enddo
                do i = 1,n(k)
                    do j = 1,idonor(i)
                        do l = 1,2
                            do m = 1,idonor(i)
                                do ni = 1,2
                                    if (loc_donor(j,i).eq.loc_to(ni,m,i)
     1                              .and.(loc_donor(m,i).eq.
     1                              (loc_to(l,j,i)))) then
                                        if (roh_sec_min(l,j,i).gt.
     1                                  roh_sec_min(ni,m,i)) then
                                            loc_to(l,j,i) = 0
                                            loc_donor(j,i) = 0
                                        else
                                            loc_to(ni,m,i) = 0
                                            loc_donor(m,i) = 0
                                        endif
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo
                enddo 
               do i = 1,n(k)
                    idonor2(i) = 0
                enddo
                do i = 1,n(k)
                    do j = 1,idonor(i)
                        if (loc_donor(j,i).ne.0) then
                            idonor2(i) = idonor2(i) + 1
                            loc_donor2(idonor2(i),i)=
     1                      loc_donor(j,i)
                            do l = 1,2
                                loc_to2(l,idonor2(i),i)=loc_to(l,j,i)
                                roh_donor_check2(l,idonor2(i),i) = 
     1                          roh_donor_check(l,j,i)
                            enddo
                        endif
                    enddo
                enddo
                iprism = 0
                ipb = 0
                itot = 0
                iz = 0
c               separate the ones that don't have double donors that
c               donate to double acceptors
                do i = 1,n(k)
                     if (idonor2(i).gt.3) then
                        iz = iz + 1
                        wt_zero(iz) = weight(i) 
                        izero(ir) = izero(ir) + 1
                        do j = 1,ndim
                            psips_zero(j,izero(ir),ir) = psips(j,i)
                        enddo
                        iwave_check(izero(ir),ir) = k
                        iwalk_check(izero(ir),ir) = i
                    else
                        iprism = iprism + 1
                        id_prism(iprism) = idonor2(i)
                        wt_prism(iprism) = weight(i)
                        nfrom_prism(iprism) = i
                        do j = 1,ndim
                            psips_prism(j,iprism) = psips(j,i)
                        enddo
                        do j = 1,idonor2(i)
                            loc_donor_prism(j,iprism) = loc_donor2(j,i)
                            icount_prism(j,iprism) = icount(j,i)
                            do l = 1,2
                                loc_to_prism(l,j,iprism)=loc_to2(l,j,i)
                            enddo
                        enddo
                        do j = 1,noxy
                            itype_prism(j,iprism) = itype(j,i)
                            do l = 1,2
                                roh_keep(l,j,iprism)=roh_sec_min(l,j,i)
                                loc_roh_keep(l,j,iprism) = 
     1                          loc_roh_sec(l,j,i)
                            enddo
                        enddo
                    endif
                enddo
c               find the other double donor
c               combine everything
                do i = 1,iprism
                    idon(i) = 0
                    do j = 1,id_prism(i)
                        do l = 1,2
                            idon(i) = idon(i) + 1
                            loc_from(idon(i),i) = loc_donor_prism(j,i)
                            loc_acceptor(idon(i),i)=loc_to_prism(l,j,i)
                        enddo
                    enddo
                enddo
                do i = 1,iprism
                    do j = 1,idon(i)
                        do l = 1,idon(i)
                            if ((loc_from(j,i).eq.loc_acceptor(l,i))
     1                      .and.(j.ne.l)) then
                                loc_use(i) = loc_from(l,i)
                            endif
                        enddo
                    enddo
                enddo
                do i = 1,iprism
                    do j = 1,idon(i)
                        if (loc_from(j,i).eq.loc_use(i)) then
                            do l = 1,idon(i)
                                if (loc_from(l,i).eq.
     1                          loc_acceptor(j,i)) then
                                    itype_prism(loc_acceptor(j,i),i) = 2
                                    exit
                                else
                                    itype_prism(loc_acceptor(j,i),i) = 3
                                endif
                            enddo
                        endif
                    enddo
                enddo
c               check if we get the lower triangle if ilower(i).eq.0
c               then dont continue further
                do i = 1,iprism
                    ilower(i) = 0
                enddo
                do i = 1,iprism
                    do j = 1,id_prism(i)
                        do l = 1,2
                            if (itype_prism(loc_to_prism(l,j,i),i).eq.3)
     1                      then
                                if (itype_prism(loc_donor_prism(j,i),i)
     1                          .eq.1) then
                                    ilower(i) = 1
                                endif
                            endif
                        enddo
                    enddo
                enddo
                nlower = 0
                do i = 1,iprism
                    do j = 1,noxy
                        iterms(j,i) = 0
                    enddo
                enddo
                do i = 1,iprism
                    do j = 1,noxy
                        do l = 1,noxy
                            if (itype_prism(j,i).eq.l) then
                                iterms(l,i) = iterms(l,i) + 1
                            endif
                        enddo
                    enddo
                enddo
                do i = 1,iprism
                    if (any(iterms([1,noxy],i).gt.1).or.(ilower(i)
     1              .eq.0)) then
                        iz = iz + 1
                        wt_zero(iz) = wt_prism(i) 
                        izero(ir) = izero(ir) + 1
                        do j = 1,ndim
                            psips_zero(j,izero(ir),ir)=psips_prism(j,i)
                        enddo
                        iwave_check(izero(ir),ir) = k
                        iwalk_check(izero(ir),ir) = nfrom_prism(i)
                    else
                        nlower = nlower + 1
                        nfrom_lower(nlower) = nfrom_prism(i)
                        do j = 1,ndim
                            psips_lower(j,nlower) = psips_prism(j,i)
                        enddo
                        wt_lower(nlower) = wt_prism(i)
                        do j = 1,noxy
                            itype_lower(j,nlower) = itype_prism(j,i)
                            do l = 1,2
                                roh_lower(l,j,nlower) = roh_keep(l,j,i)
                                loc_roh_lower(l,j,nlower) = 
     1                          loc_roh_keep(l,j,i)
                            enddo
                        enddo
                    endif
                enddo
c               find the connector to the bottom to the top of the prism
                do i = 1,nlower
                    do j = 1,noxy
                        do l = 1,2
                            itop(l,j,i) = 0
                        enddo
                    enddo
                enddo
                do i = 1,nlower 
                    do j = 1,noxy
                        do l = 1,2
                            if (itype_lower(j,i).eq.3) then
                                itop(l,j,i) = 1
                            endif
                        enddo
                    enddo
                enddo
                do i = 1,nlower
                    do j =1 ,noxy
                        do l =1,2
                            if (itop(l,j,i).eq.1) then
                                roh_upper(l,i) = roh_lower(l,j,i)
                                loc_to_upper(l,i) = loc_roh_lower(l,j,i)
                                loc_from_upper(l,i) = j
                            endif
                        enddo
                    enddo
                enddo
                do i = 1,nlower
                    roh_upper_min(i) = minval(roh_upper([1,2],i),1)
                    loc_min(i) = minloc(roh_upper([1,2],i),1)
                    loc_upper_min(i) = loc_to_upper(loc_min(i),i)
                enddo
                do i =1,nlower
                    do j = 1,noxy
                        itype_lower(loc_upper_min(i),i) = 4
                    enddo
                enddo
c               Make sure there are 4 things connected
                nupper = 0
                do i = 1,nlower
                    if (any(itype_lower(:,i).eq.4)) then
                        nupper = nupper + 1
                        nfrom_upper(nupper) = nfrom_lower(i)
                        do j = 1,ndim
                            psips_upper(j,nupper) = psips_lower(j,i)
                        enddo
                        wt_upper(nupper) = wt_lower(i)
                        do j = 1,noxy
                            do l = 1,2
                                roh_rest(l,j,nupper) = roh_lower(l,j,i)
                                loc_roh_upper(l,j,nupper) = 
     1                          loc_roh_lower(l,j,i)
                            enddo
                            itype_upper(j,nupper) = itype_lower(j,i)
                        enddo
                        do l = 1,2
                            loc_in_upper(l,nupper) = loc_to_upper(l,i)
                            loc_from_upper2(l,nupper) = 
     1                      loc_from_upper(l,i)
                        enddo
                    else
                        iz = iz + 1
                        wt_zero(iz) = wt_lower(i) 
                        izero(ir) = izero(ir) + 1
                        do j = 1,ndim
                            psips_zero(j,izero(ir),ir)=psips_lower(j,i)
                        enddo
                        iwave_check(izero(ir),ir) = k
                        iwalk_check(izero(ir),ir) = nfrom_lower(i)
                    endif
                enddo
c               find the 2 double acceptors in the top part
c               l = 2 is water 2 and l = 1 
c               connect water 3 to 2 then to 1 
c               check 1 connects back to 3
                do i = 1,nupper
                    do j = 1,noxy
                        if (itype_upper(j,i).eq.4) then
                            do l = 1,2
                                if (l.eq.2) then
                                    itype_upper(
     1                              loc_roh_upper(l,j,i),i) = 5
                                else if (l.eq.1) then
                                    itype_upper(
     1                              loc_roh_upper(l,j,i),i) = 6
                                endif
                            enddo
                        endif
                    enddo
                enddo
c               check connections from top down check for itype upper
c               5 and 6 to be connected to 2 and 3 respectively
                do i = 1,nupper
                    do j = 1,noxy
                        do l= 1,2
                            icheck(l,j,i) = 0
                        enddo
                    enddo
                enddo
                do i = 1,nupper
                    do j = 1,noxy
                        do l = 1,2
                            if (itype_upper(j,i).eq.5) then
                                icheck(l,j,i) = 1
                            else if (itype_upper(j,i).eq.6) then
                                icheck(l,j,i) = 2
                            else if (itype_upper(j,i).eq.2) then
                                icheck(l,j,i) = 3
                            endif
                        enddo
                    enddo
                enddo
                do i = 1,nupper
                    do j = 1,noxy
                        do l =1 ,2
                            if (icheck(l,j,i).eq.1) then
                                roh_check1(l,i) = roh_rest(l,j,i)
                                loc_to_check1(l,i)=loc_roh_upper(l,j,i)
                                loc_from_check1(l,i) = j
                            else if (icheck(l,j,i).eq.2) then
                                roh_check2(l,i) = roh_rest(l,j,i)
                                loc_to_check2(l,i)=loc_roh_upper(l,j,i)
                                loc_from_check2(l,i) = j
                            else if (icheck(l,j,i).eq.3) then
                                roh_check3(l,i) = roh_rest(l,j,i)
                                loc_to_check3(l,i)=loc_roh_upper(l,j,i)
                                loc_from_check3(l,i) = j
                            endif
                        enddo
                    enddo
                enddo
                do i = 1,nupper
                    roh_check1_min(i) = minval(roh_check1([1,2],i),1)
                    loc_roh_check1_min(i) =minloc(roh_check1([1,2],i),1)
                    loc_roh_check1(i) = loc_to_check1(
     1              loc_roh_check1_min(i),i)
                    roh_check2_min(i) = minval(roh_check2([1,2],i),1)
                    loc_roh_check2_min(i) =minloc(roh_check2([1,2],i),1)
                    loc_roh_check2(i) = loc_to_check2(
     1              loc_roh_check2_min(i),i)
                enddo
                nfinal = 0
                do i = 1,nupper
                    if (all(itype_upper(:,i).ne.0)) then
                        nfinal = nfinal + 1
                        weight_final(nfinal) = wt_upper(i)
                        do j = 1,noxy
                            itype_final(j,nfinal) = itype_upper(j,i)
                        enddo
                    else
                        iz = iz + 1
                        wt_zero(iz) = wt_upper(i) 
                        izero(ir) = izero(ir) + 1
                        do j = 1,ndim
                            psips_zero(j,izero(ir),ir)=psips_upper(j,i)
                        enddo
                        iwave_check(izero(ir),ir) = k
                        iwalk_check(izero(ir),ir) = nfrom_upper(i)
                    endif
                enddo
c               find atom 6
                do i = 1,nfinal
                    if (itype_final(6,i).eq.1) then
                         deut(1,ir) = deut(1,ir) + (weight_final(i)/
     1                   tot_weight(k))
                     else if (itype_final(6,i).eq.2) then
                         deut(2,ir) = deut(2,ir) + (weight_final(i)/
     1                   tot_weight(k))
                     else if (itype_final(6,i).eq.3) then
                         deut(3,ir) = deut(3,ir) + (weight_final(i)/
     1                   tot_weight(k))
                     else if (itype_final(6,i).eq.4) then
                         deut(4,ir) = deut(4,ir) + (weight_final(i)/
     1                   tot_weight(k))
                     else if (itype_final(6,i).eq.5) then
                         deut(5,ir) = deut(5,ir) + (weight_final(i)/
     1                   tot_weight(k))
                     else if (itype_final(6,i).eq.6) then
                         deut(6,ir) = deut(6,ir) + (weight_final(i)/
     1                   tot_weight(k))
                     endif
                enddo
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
                iwave_bad(it) = iwave_check(i,k)
                iwalk_bad(it) = iwalk_check(i,k)
            enddo
        enddo
        write(11,*) it
        do i = 1,it
            write(11,*) (psips_check(j,i),j=1,ndim),iwave_bad(i),
     1      iwalk_bad(i)
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
