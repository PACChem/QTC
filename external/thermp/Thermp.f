
      program thermp

c
c     This program reads the canonical partition functions from paper
c     and 
c     (i) converts them into a table of Cp, S, and H.
c     (ii) calls the Gordon and McBride code pac99
c     (iii) converts the polynomials from Nada new format/old polynomial to
c     ChemKin format 
c

      implicit double precision (a-h,o-z)
      integer dimel, dimeld
      parameter (maxt=1000, maxsp=200)
      parameter (dimel=30, dimeld=6)
      dimension tt(maxt),qtmp(maxsp*3)
      dimension qlpsi(maxt,maxsp),dlqpsi(maxt,maxsp),d2lqpsi(maxt,maxsp)
c     dimension qpsi(maxt,maxsp),dqpsi(maxt,maxsp),d2qpsi(maxt,maxsp),
c    $ dlqpsi(maxt,maxsp),d2lqpsi(maxt,maxsp)
      dimension qlps(maxt),dlqpdt(maxt),d2lqpdt2(maxt)
c     dimension qps(maxt),dqpsdt(maxt),dlqpdt(maxt),d2lqpdt2(maxt)
      dimension cps(maxt),ents(maxt),hints(maxt),acoeff(14)
      dimension hinc(dimel),hincd(dimeld)
      character*80 spname,prefix,infile
      character*24 stoich,anelem,stoichi
      character*80 today
      character*2 elemname,elab(dimel),elabd(dimeld)

      include 'data.fi'
      data hinc /6.197, 6.197, 6.197, 6.197, 6.197, 4.412, 5.360, 
     $ 1.050, 3.217, 4.636, 6.323, 6.870, 1.222, 4.540, 5.647,
     $ 6.247, 9.342, 5.004, 5.745, 4.824, 6.364, 6.35, 1.950, 
     $ 4.998, 5.736, 4.632, 6.460, 7.088, 7.489, 7.711 /
      data hincd /8.680, 8.468, 8.825, 9.181, 24.52, 8.670/
      data elab / 'He', 'Ne', 'Ar', 'Kr', 'Xe', 'S ', 'P ', 'C ',
     $ 'Si', 'Ge', 'Sn', 'Pb', 'B ', 'Al', 'Zn', 'Cd', 'Hg', 'Cu', 
     $ 'Ag', 'Ti', 'U ', 'Th', 'Be', 'Mg', 'Ca', 'Li', 'Na', 'K ',
     $ 'Rb', 'Cs' /
      data elabd /'O ', 'H ', 'F ', 'Cl', 'Br', 'N ' /

c     First evaluate the canonical partition function for the complex 
c     for tt = tt + dp and tt = tt - dp

      tp = dtderiv*dkbltz/cautoerg
      tm = -dtderiv*dkbltz/cautoerg
      tmin = dtderiv*dkbltz/cautoerg
      t0 = 0.0d0
      pstp = 760.0d0/pconv
      OPEN (UNIT=18,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     & STATUS='unknown',FILE='thermp.out')
      OPEN (UNIT=25,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     & STATUS='unknown',FILE='thermp.dat')
      read (25,*) nwell, nprod
      read (25,*) 
      nspecies = nwell + nprod
      read (25,*) nt
      read (25,*) 
      OPEN (UNIT=27,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     & STATUS='unknown',FILE='pf.dat')
c     read in the partition function data
      read(27,*)
      read(27,*)
      do it = 1 , nt + 1
         ntot = nspecies*3
         read (27,*) tt(it),(qtmp(itmp),itmp=1,ntot)
         iind = 1
         do ispecies = 1 , nspecies
            qlpsi(it,ispecies) = qtmp(iind)
            dlqpsi(it,ispecies) = qtmp(iind+1)
            d2lqpsi(it,ispecies) = qtmp(iind+2)
c           dqpsi(it,ispecies) = qtmp(iind+1)
c           d2qpsi(it,ispecies) = qtmp(iind+2)
            iind = iind+3
         enddo
      enddo
      do ispecies = 1 , nspecies
         qlpsi(nt+2,ispecies) = qlpsi(nt+1,ispecies)
         dlqpsi(nt+2,ispecies) = dlqpsi(nt+1,ispecies)
         d2lqpsi(nt+2,ispecies) = d2lqpsi(nt+1,ispecies)
c        qpsi(nt+2,ispecies) = qpsi(nt+1,ispecies)
c        dqpsi(nt+2,ispecies) = dqpsi(nt+1,ispecies)
c        d2qpsi(nt+2,ispecies) = d2qpsi(nt+1,ispecies)
      enddo
c     generate the logarithmic derivatives of the partition function
c     and convert units to au
      qconv = dlog(cautocm**3)
c     qconv = cautocm**3
c     write (18,*) 'qconv test',cautocm**3,qconv
      dlqconv = cautoerg/dkbltz
      d2lqconv = dlqconv*cautoerg/dkbltz
c     dqconv = qconv*cautoerg/dkbltz
c     dq2conv = dqconv*cautoerg/dkbltz
      do it = 1 , nt + 2
         tt(it) = tt(it)*dkbltz/cautoerg
         do ispecies = 1 , nspecies
            qlpsi(it,ispecies) = qlpsi(it,ispecies)+qconv
            dlqpsi(it,ispecies) = dlqpsi(it,ispecies)*dlqconv
            d2lqpsi(it,ispecies) = d2lqpsi(it,ispecies)*d2lqconv
c           qpsi(it,ispecies) = qpsi(it,ispecies)*qconv
c           dqpsi(it,ispecies) = dqpsi(it,ispecies)*dqconv
c           d2qpsi(it,ispecies) = d2qpsi(it,ispecies)*dq2conv
c           dlqpsi(it,ispecies) = dqpsi(it,ispecies)/qpsi(it,ispecies)
c           d2lqpsi(it,ispecies) = d2qpsi(it,ispecies)/
c    $       qpsi(it,ispecies) - 
c    $       (dqpsi(it,ispecies)/qpsi(it,ispecies))**2
         enddo
      enddo

      do 100 ispecies = 1 , nspecies
c
c Open .i97 for each species
c Put appropriate labels .i97 for NASA polynomials
C generate thermo data
c print thermo data in .i97
c repeat labels and thermo data for CHEMKIN polynomials
c do next two steps outside code
c call pac99
c read .ch97 and convert second polynomial to CHEMKIN style
c
         read (25,*) href,tref
         itref = idint(tref)
         if ((itref.ne.0).and.(itref.ne.298)) then
            write (18,*) 'error: Tref in fitdat must be', 
     $       'either 0 or 298',tref
            stop
         endif
c        OPEN (UNIT=26,ACCESS='SEQUENTIAL',FORM='FORMATTED',
c    &   STATUS='unknown',FILE='fitdat.dat')
         read (25,*) spname
c        write (26,8001) spname
c001     format (a80)
         prefix = spname
         ln = INDEX(prefix,' ') - 1
         infile = prefix(1:ln)//'.i97'
         OPEN (15,FILE=infile,STATUS='unknown',FORM='formatted')
         hincs = 0.0d0
         stoich = ''
c loop over elements in species
   10    continue
         read (25,*) elemname,dnelem
         if (elemname.eq.'**') go to 20
         nelem = idint(dnelem)
         write (26,8002) elemname,nelem
 8002    format (a2,i2)
         write(anelem,'(I0)') nelem

         ln1 = INDEX(elemname,' ') - 1
         ln2 = INDEX(anelem,' ') - 1
         stoichi = elemname(1:ln1)//anelem(1:ln2)
         ln3 = INDEX(stoich,' ') - 1
         ln4 = INDEX(stoichi,' ') - 1
         stoich = stoich(1:ln3)//stoichi(1:ln4)
c        write (15,*) 'stoich test',stoich,stoichi
c if input h is for zero K convert to 298 K
         if (itref.eq.0) then 
            iels = 0
            do iel = 1 , dimel
               if (elemname.eq.elab(iel)) then
                  iels = iel
                  hincs = hincs + dnelem*hinc(iels)
               endif
            enddo
            do ield = 1 , dimeld
               if (elemname.eq.elabd(ield)) then
                  iels = ield
                  hincs = hincs + dnelem*hincd(iels)/2.0d0
               endif
            enddo
            if (iels.eq.0) then
               write (18,*) 'error: element is not in element list',
     $          elemname
               stop
            endif
         endif
         go to 10
   20    continue
c        write (26,8002) elemname
         tbreak = dnelem

c        tref = tref*dkbltz/cautoerg
         t298 = 298.15d0*dkbltz/cautoerg
         tt(nt+1) = t298
         h298 = href - hincs/ckctokj
c presume h = 0 for T = 0

         write (18,*) 'thermo properties for species ',ispecies
         write (18,*) ' T(K) Cp(cal/(mol.K)) S(cal/(mol.K)) ',
     $    'H-H(1)(kcal/mol)  Q/V   H(kcal/mol) '
         do it = 1 , nt+1
            qlps(it) = qlpsi(it,ispecies)
c           qps(it) = qpsi(it,ispecies)
c           dqpsdt(it) = dqpsi(it,ispecies)
            dlqpdt(it) = dlqpsi(it,ispecies)
            d2lqpdt2(it) = d2lqpsi(it,ispecies)

            ent = qlps(it)+tt(it)*dlqpdt(it)
c           ent = dlog(qps(it))+tt(it)*dqpsdt(it)/qps(it)
            ent = ent + log(tt(it)*exp(1.0d0)/pstp)
            ent = ent*1.987d0
            uint = tt(it)**2*dlqpdt(it)
c           uint = tt(it)**2*dqpsdt(it)/qps(it)
            hint = uint + tt(it)
            hint = hint
            if (it.eq.1) hint1=hint
            cv = 2.0d0*tt(it)*dlqpdt(it) + tt(it)**2*d2lqpdt2(it) 
            cv = cv*1.987d0
            cp = cv + 1.987d0
            write (18,8000) tt(it)*cautoerg/dkbltz,cp,ent,
     $       (hint-hint1)*cautoerg/ckctoerg,exp(qlps(it))/cautocm**3,
     $       hint*cautoerg/ckctoerg
8000        format (1x,3g12.5,g13.5,2g12.5)
            cps(it) = cp
            ents(it) = ent
            hints(it) = hint*cautoerg/ckctoerg
         enddo
         write (18,*) 'thermodynamic fit results'
         if (itref.eq.0) h298 = h298+hints(nt+1)
         write (18,*) 'h298 final',h298
c        write (26,*) h298
         write (15,9000) 'NAME  ',spname
 9000    format(a6,a80)
         write (15,9001) stoich,'HF298',h298*ckctokj*1000,'JOULES'
 9001    format (a24,a5,F12.0,1x,a6)
         call date_and_time (today)
         write (15,9002) 'DATE  ',today
 9002    format (a6,a8)
         write (15,9003) 'REFN  ME'
 9003    format (a8)
         t1 = tt(1)*cautoerg/dkbltz
         t2 = 200.
         t3 = tbreak
         t4 = tt(nt)*cautoerg/dkbltz
         write (15,9005) 'LSTS  ','OLD','T',t1,'T',t2,'T',t3
 9005    format (a6,a3,15x,a1,f11.0,6x,a1,f11.0,6x,a1,f11.0)
         write (15,9006) 'LSTS  ','T',t4
 9006    format (a6,a1,f11.0)
         write (15,9007) 'OUTP  ','LSQS  '
 9007    format (2a6)
         write (15,9008) 'METH  ','READIN','KJOULE','BAR'
 9008    format (A6,A6,12x,a6,12x,a3)
c
c repeat for chemkin polynomials
c
         ttt = 0.0d0
         cpst = 0.0d0
         entst = 0.0
c        ht = -h298
         ht = -hints(nt+1)*ckctokj
         write (15,9009) 'T',ttt,'CP',cpst,'S',
     $    entst,'H-H2',ht
         do it = 1 , nt
            ttt = tt(it)*cautoerg/dkbltz
            ht = hints(it)*ckctokj
c           ht = (hints(it)-hints(nt+1))*ckctokj
            cpst = cps(it)*ckctokj
            entst = ents(it)*ckctokj
            write (15,9009) 'T',ttt,'CP',cpst,'S',
     $       entst,'H-H0',ht
 9009       format (6x,A1,f15.3,A4,f13.3,A4,f14.3,A7,f11.3)
            if ((tt(it+1).gt.tt(nt+1)).and.(tt(it).lt.tt(nt+1))) then
               ttt = tt(nt+1)*cautoerg/dkbltz
               ht = hints(nt+1)*ckctokj
c              ht = (hints(nt+1)-hints(nt+1))*ckctokj
               cpst = cps(nt+1)*ckctokj
               entst = ents(nt+1)*ckctokj
               write (15,9009) 'T',ttt,'CP',cpst,'S',
     $       entst,'H-H0',ht
            endif
         enddo
c        write (15,9010) 'FINISH'
         write (15,9010) 'FINISH'
9010     format (a6)
c
         write (15,9000) 'NAME  ',spname
         write (15,9001) stoich,'HF298',h298*ckctokj*1000,'JOULES'
         call date_and_time (today)
         write (15,9002) 'DATE  ',today
         write (15,9003) 'REFN  ME'
         write (15,9011) 'LSTS  ','T',t1,'T',t2,'T',t3
 9011    format (a6,a1,f11.0,6x,a1,f11.0,6x,a1,f11.0)
         write (15,9007) 'OUTP  ','LSQS  '
         write (15,9008) 'METH  ','READIN','KJOULE','BAR'
         ttt = 0.0d0
         cpst = 0.0d0
         entst = 0.0
c        ht = -h298
         ht = -hints(nt+1)*ckctokj
         
         write (15,9009) 'T',ttt,'CP',cpst,'S',
     $    entst,'H-H2',ht
         do it = 1 , nt
            ttt = tt(it)*cautoerg/dkbltz
            ht = hints(it)*ckctokj
c           ht = (hints(it)-hints(nt+1))*ckctokj
            cpst = cps(it)*ckctokj
            entst = ents(it)*ckctokj
            write (15,9009) 'T',ttt,'CP',cpst,'S',
     $       entst,'H-H0',ht
            if ((tt(it+1).gt.tt(nt+1)).and.(tt(it).lt.tt(nt+1))) then
               ttt = tt(nt+1)*cautoerg/dkbltz
               ht = hints(nt+1)*ckctokj
c              ht = (hints(nt+1)-hints(nt+1))*ckctokj
               cpst = cps(nt+1)*ckctokj
               entst = ents(nt+1)*ckctokj
               write (15,9009) 'T',ttt,'CP',cpst,'S',
     $       entst,'H-H0',ht
            endif
         enddo
         write (15,9010) 'FINISH'
         close (UNIT=15,status='keep')

c 
c run pac99
c 
c     call pac99run(spname)

c
c convert polynomial from NASA format to 
c ChemKin format
c

 100  continue

c     close (unit=18,status='keep')

      stop
      end


c     **********************************************************************


