C This program computes the equilibrium of the stochastic-beta economy in
C "Income and Wealth Heterogeneity in the Macroeconomy" (co-authored with
C Per Krusell of the University of Rochester).


      implicit real*8 (a-h,o-z)

      parameter (nkpts=132,durug=1.5D+00,nmupts=4,
     *           delta=0.025D+00,alpha=0.36D+00,unempg=0.04D+00,
     *           mxloop=12,sfac=0.25D+00,hfix=0.3271D+00,
     *           durgd=8.0D+00,unempb=0.1D+00,kgrid=278,
     *           xkbor=-2.4D+00,mgrid=30,
     *           durbd=8.0D+00,zgood=1.01D+00,zbad=0.99D+00,
     *           npbhat=2,durub=2.5D+00,ntop=9,
     *           xkglow=xkbor,xkghgh=25.0D+00,nlzpts=201,
     *           nbetas=3,nzpts=2,nrspts=132)

      common/cpr/pr(4,4),prbeta(nbetas,nbetas),prob(nbetas)
      common/cgrid/xkgpts(kgrid),xmgpts(mgrid)
      common/clzpts/xlzpts(nlzpts),xlzdat(nlzpts)
      common/cbetas/betas(nbetas)
      common/ctop/toppct(ntop)
      common/ccoef/coefk(kgrid-1,mgrid-1,4,nbetas,4),
     *             optkgd(kgrid,mgrid,4,nbetas)
      common/cvs/v1(nkpts,nmupts,nbetas,2),v2(nkpts,nmupts,nbetas,2),
     *           v3(nkpts,nmupts,nbetas,2),v4(nkpts,nmupts,nbetas,2),
     *           v1dp(nkpts,nmupts,nbetas,2),
     *           v2dp(nkpts,nmupts,nbetas,2),
     *           v3dp(nkpts,nmupts,nbetas,2),
     *           v4dp(nkpts,nmupts,nbetas,2)
      common/cws/w1(nkpts,nbetas,2),w2(nkpts,nbetas,2),
     *           w3(nkpts,nbetas,2),w4(nkpts,nbetas,2),
     *           w1dp(nkpts,nbetas,2),w2dp(nkpts,nbetas,2),
     *           w3dp(nkpts,nbetas,2),w4dp(nkpts,nbetas,2)

      dimension y(nkpts,nmupts,4,nbetas),y2(nkpts,nmupts,4,nbetas),
     *          optk(nkpts,nmupts,4,nbetas),
     *          optk2(nkpts,nmupts,4,nbetas),
     *          xkpts(nkpts),bhat1(npbhat,1),bhat2(npbhat,1),
     *          ndone(nkpts,nmupts,4,nbetas),
     *          grad(nkpts,nmupts,4,nbetas),
     *          vdiffs(nkpts,nmupts,4,nbetas),
     *          ddiffs(nkpts,nmupts,4,nbetas),
     *          r0(nmupts,2),w0(nmupts,2),r0cur(2),w0cur(2),
     *          xmupts(nmupts),topmn(ntop),topsd(ntop),
     *          topmx(ntop),topmi(ntop)

      seed1 = 3243090.0D+00
      seed2 = 254367.0D+00

      one = 1.0D+00

      toppct(1) = 0.005D+00
      toppct(2) = 0.01D+00
      toppct(3) = 0.05D+00
      toppct(4) = 0.1D+00
      toppct(5) = 0.2D+00
      toppct(6) = 0.3D+00
      toppct(7) = 0.4D+00
      toppct(8) = 0.6D+00
      toppct(9) = 0.8D+00

      pgg00 = (durug-one)/durug
      pbb00 = (durub-one)/durub
      pbg00 = 1.25D+00*pbb00
      pgb00 = 0.75D+00*pgg00
      pgg01 = (unempg - unempg*pgg00)/(one-unempg)
      pbb01 = (unempb - unempb*pbb00)/(one-unempb)
      pbg01 = (unempb - unempg*pbg00)/(one-unempg)
      pgb01 = (unempg - unempb*pgb00)/(one-unempb)
      pgg = (durgd-one)/durgd
      pgb = one - (durbd-one)/durbd

      pgg10 = one - (durug-one)/durug
      pbb10 = one - (durub-one)/durub
      pbg10 = one - 1.25D+00*pbb00
      pgb10 = one - 0.75D+00*pgg00
      pgg11 = one - (unempg - unempg*pgg00)/(one-unempg)
      pbb11 = one - (unempb - unempb*pbb00)/(one-unempb)
      pbg11 = one - (unempb - unempg*pbg00)/(one-unempg)
      pgb11 = one - (unempg - unempb*pgb00)/(one-unempb)
      pbg = one - (durgd-one)/durgd
      pbb = (durbd-one)/durbd

      pr(1,1) = pgg*pgg11
      pr(2,1) = pbg*pbg11
      pr(3,1) = pgg*pgg01
      pr(4,1) = pbg*pbg01
      pr(1,2) = pgb*pgb11
      pr(2,2) = pbb*pbb11
      pr(3,2) = pgb*pgb01
      pr(4,2) = pbb*pbb01
      pr(1,3) = pgg*pgg10
      pr(2,3) = pbg*pbg10
      pr(3,3) = pgg*pgg00
      pr(4,3) = pbg*pbg00
      pr(1,4) = pgb*pgb10
      pr(2,4) = pbb*pbb10
      pr(3,4) = pgb*pgb00
      pr(4,4) = pbb*pbb00

      open(1,file='betas.dat',form='formatted')
      read(1,265) brho1,brho2,brho3
265   format(3f10.5)
      read(1,265) betas(1),betas(2),betas(3)
      close(1)

      prbeta(1,1) = brho1
      prbeta(2,1) = one - brho1
      prbeta(3,1) = 0.0D+00
      prbeta(1,2) = (one-brho2)/2.0D+00
      prbeta(2,2) = brho2
      prbeta(3,2) = (one-brho2)/2.0D+00
      prbeta(1,3) = 0.0D+00
      prbeta(2,3) = one - brho3
      prbeta(3,3) = brho3
      call pinvar(prbeta,prob)
      write(6,156) brho1,brho2,brho3
      write(6,156) betas(1),betas(2),betas(3)
      write(6,156) prob(1),prob(2),prob(3)
156   format(3f15.9)

      write(6,325) ((pr(i,j),j=1,4),i=1,4)
325   format(4f10.6)

      xlow = 0.0D+00
      xhigh = one
      xinc = (xhigh-xlow)/dble(real(nlzpts-1))
      xlzpts(1) = xlow
      do 623 i = 2,nlzpts
         xlzpts(i) = xlzpts(i-1) + xinc
623   continue

      xlzdat(1) = 0.0D+00
      xlzdat(141) = 0.122D+00
      xlzdat(161) = 0.212D+00
      xlzdat(181) = 0.359D+00
      xlzdat(191) = 0.489D+00
      xlzdat(199) = 0.704D+00
      xlzdat(201) = one

      open(2,file='heter40w.out',form='formatted')

C      xmupts(1) = 10.85D+00
C      xmupts(2) = 11.35D+00
C      xmupts(3) = 11.85D+00
C      xmupts(4) = 12.35D+00

      xmupts(1) = 11.1D+00
      xmupts(2) = 11.52D+00
      xmupts(3) = 11.94D+00
      xmupts(4) = 12.36D+00
      xmglow = xmupts(1)
      xmghgh = xmupts(nmupts)

      xklow = xkbor
      xkhigh = 5500.0D+00
      xkcut1 = -2.3D+00
      xkcut2 = -2.3D+00
      ncalop = 0

      do 58 i = 1,nmupts
         agghrs = hfix*(one-unempg)
         r0(i,1) = zgood*alpha*(xmupts(i)**(alpha-one))*
     *                         (agghrs**(one-alpha))
         w0(i,1) = zgood*(one-alpha)*(xmupts(i)**alpha)*
     *                               (agghrs**(-alpha))
         agghrs = hfix*(one-unempb)
         r0(i,2) = zbad*alpha*(xmupts(i)**(alpha-one))*
     *                        (agghrs**(one-alpha))
         w0(i,2) = zbad*(one-alpha)*(xmupts(i)**alpha)*
     *                              (agghrs**(-alpha))
         write(6,598) xmupts(i),r0(i,1),w0(i,1),r0(i,2),w0(i,2)
598      format(5f15.8)
 58   continue

      nloop = 0

      open(1,file='agglaw.in',form='formatted')
      read(1,582) a0new,a1new
      read(1,582) c0new,c1new
      close(1)

  1   nloop = nloop + 1
      if (nloop .gt. mxloop) goto 2

      if (nloop .eq. 1) then
         mxiter = 1200
      elseif (nloop .lt. mxloop) then
         mxiter = 1200
      else
         mxiter = 1200
      endif

      a0 = a0new
      a1 = a1new
      c0 = c0new
      c1 = c1new

      nres = 0
      nkgd = 0
      if (nkgd .eq. 1) goto 3

      if (nloop .eq. 1) then

         open(7,file='heter40w.in',form='formatted')
         do 457 k = 1,nbetas
            do 56 i = 1,nkpts
               do 57 j = 1,nmupts
                  read(7,445) xkpts(i),xmupts(j),y2(i,j,1,k),
     *                        y2(i,j,2,k),y2(i,j,3,k),y2(i,j,4,k)
445               format(2f8.3,4f20.8)
C                  y2(i,j,1,k) = 0.01D+00
C                  y2(i,j,2,k) = 0.01D+00
C                  y2(i,j,3,k) = 0.01D+00
C                  y2(i,j,4,k) = 0.01D+00
 57            continue
 56         continue
457      continue
         close(7)

         xkbd = xkpts(nkpts)

         if (nres .eq. 0) then
            do 20 i = 1,nkpts
               do 21 j = 1,nmupts
                  do 22 k = 1,4
                     do 222 kk = 1,nbetas
                        optk2(i,j,k,kk) = xkpts(i)
222                  continue
 22               continue
 21            continue
 20         continue
         else

            open(1,file='heter40w.res',form='formatted')
            do 982 k = 1,nbetas

               do 972 i = 1,nrspts
                  do 996 j = 1,nmupts
                     read(1,523) dummy,dummy,
     *                           optk2(i,j,1,k),optk2(i,j,2,k),
     *                           optk2(i,j,3,k),optk2(i,j,4,k),
     *                           dummy,dummy,dummy,dummy
996               continue
972            continue

               do 678 i = nrspts+1,nkpts
                  do 995 j = 1,nmupts
                     optk2(i,j,1,k) = (optk2(nrspts,j,1,k)/
     *                                 xkpts(nrspts))*xkpts(i)
                     optk2(i,j,2,k) = (optk2(nrspts,j,2,k)/
     *                                 xkpts(nrspts))*xkpts(i)
                     optk2(i,j,3,k) = (optk2(nrspts,j,3,k)/
     *                                 xkpts(nrspts))*xkpts(i)
                     optk2(i,j,4,k) = (optk2(nrspts,j,4,k)/
     *                                 xkpts(nrspts))*xkpts(i)
995               continue
678            continue

982         continue
            close(1)
         endif

      endif

      do 30 k = 1,mxiter

         call interp(xkpts,xmupts,y2,a0,a1,c0,c1)

         do 11 i = 1,nkpts
            do 12 j = 1,nmupts
               do 13 kk = 1,4
                  do 213 jj = 1,nbetas
                     y(i,j,kk,jj) = y2(i,j,kk,jj)
                     optk(i,j,kk,jj) = optk2(i,j,kk,jj)
213               continue
 13            continue
 12         continue
 11      continue

         nprint = 0
         do 10 i = 1,nkpts
            do 73 j = 1,nmupts
              do 74 jj = 1,nbetas

                  xk = xkpts(i)
                  nmucur = j

                  neps = 1
                  ngood = 1
                  if ( (xk .lt. xkcut1) .or. (k .le. ncalop) ) then
                     call calopt(xk,neps,ngood,xknew1,value1,niter,
     *                           ndone1,1.0D-10,nmucur,
     *                           g1,xkpts,r0,w0,xkhigh,xklow,jj,
     *                           nprint)
C                     write(6,785) xkpts(i),xmupts(j),neps,ngood,
C     *                            xknew1,value1,ndone1,jj
C785                  format(2f10.4,2i5,2f20.10,i5,i3)
C                     pause
                  else
                     call altopt(xk,neps,ngood,xknew1,value1,niter,
     *                           ndone1,1.0D-10,nmucur,
     *                           optk2(i,j,1,jj),grad1,xkpts,r0,w0,jj)
C                     write(6,785) xkpts(i),xmupts(j),neps,ngood,
C     *                            xknew1,value1,ndone1,jj
C                     if (ndone1 .eq. 0)
C                        call calopt(xk,neps,ngood,xknew1,value1,niter,
C     *                              ndone1,1.0D-10,nmucur,
C     *                              g1,xkpts,r0,w0,xkhigh,xklow,jj,
C     *                              nprint)
C                     write(6,785) xkpts(i),xmupts(j),neps,ngood,
C     *                            xknew1,value1,ndone1,jj
C                     pause
                  endif
                  ndone(i,j,1,jj) = ndone1
                  grad(i,j,1,jj) = g1
                  optk2(i,j,1,jj) = xknew1
                  y2(i,j,1,jj) = value1

                  neps = 1
                  ngood = 2
                  if ( (xk .lt. xkcut1) .or. (k .le. ncalop) ) then
                     call calopt(xk,neps,ngood,xknew2,value2,niter,
     *                           ndone2,1.0D-10,nmucur,
     *                           g2,xkpts,r0,w0,xkhigh,xklow,jj,
     *                           nprint)
C                     write(6,785) xkpts(i),xmupts(j),neps,ngood,
C     *                            xknew2,value2,ndone2,jj
C                     pause
                  else
                     call altopt(xk,neps,ngood,xknew2,value2,niter,
     *                           ndone2,1.0D-10,nmucur,
     *                           optk2(i,j,2,jj),grad2,xkpts,r0,w0,jj)
C                     write(6,785) xkpts(i),xmupts(j),neps,ngood,
C     *                            xknew2,value2,ndone2,jj
C                     if (ndone2 .eq. 0)
C                        call calopt(xk,neps,ngood,xknew2,value2,niter,
C     *                              ndone2,1.0D-10,nmucur,
C     *                              g2,xkpts,r0,w0,xkhigh,xklow,jj,
C     *                              nprint)
C                     write(6,785) xkpts(i),xmupts(j),neps,ngood,
C     *                            xknew2,value2,ndone2,jj
C                     pause
                  endif
                  ndone(i,j,2,jj) = ndone2
                  grad(i,j,2,jj) = g2
                  optk2(i,j,2,jj) = xknew2
                  y2(i,j,2,jj) = value2

                  neps = 0
                  ngood = 1
                  if ( (xk .lt. xkcut2) .or. (k .le. ncalop) ) then
                     call calopt(xk,neps,ngood,xknew3,value3,niter,
     *                           ndone3,1.0D-10,nmucur,
     *                           g3,xkpts,r0,w0,xkhigh,xklow,jj,
     *                           nprint)
C                     write(6,785) xkpts(i),xmupts(j),neps,ngood,
C     *                            xknew3,value3,ndone3,jj
C                     pause
                  else
                     call altopt(xk,neps,ngood,xknew3,value3,niter,
     *                           ndone3,1.0D-10,nmucur,
     *                           optk2(i,j,3,jj),grad3,xkpts,r0,w0,jj)
C                     write(6,785) xkpts(i),xmupts(j),neps,ngood,
C     *                            xknew3,value3,ndone3,jj
C                     if (ndone3 .eq. 0)
C                        call calopt(xk,neps,ngood,xknew3,value3,niter,
C     *                              ndone3,1.0D-10,nmucur,
C     *                              g3,xkpts,r0,w0,xkhigh,xklow,jj,
C     *                              nprint)
C                     write(6,785) xkpts(i),xmupts(j),neps,ngood,
C     *                            xknew3,value3,ndone3,jj
C                     pause
                  endif
                  ndone(i,j,3,jj) = ndone3
                  grad(i,j,3,jj) = g3
                  optk2(i,j,3,jj) = xknew3
                  y2(i,j,3,jj) = value3

                  neps = 0
                  ngood = 2
                  if ( (xk .lt. xkcut2) .or. (k .le. ncalop) ) then
                     call calopt(xk,neps,ngood,xknew4,value4,niter,
     *                           ndone4,1.0D-10,nmucur,
     *                           g4,xkpts,r0,w0,xkhigh,xklow,jj,
     *                           nprint)
C                     write(6,785) xkpts(i),xmupts(j),neps,ngood,
C     *                            xknew4,value4,ndone4,jj
C                     pause
                  else
                     call altopt(xk,neps,ngood,xknew4,value4,niter,
     *                           ndone4,1.0D-10,nmucur,
     *                           optk2(i,j,4,jj),grad4,xkpts,r0,w0,jj)
C                     write(6,785) xkpts(i),xmupts(j),neps,ngood,
C     *                            xknew4,value4,ndone4,jj
C                     if (ndone4 .eq. 0)
C                        call calopt(xk,neps,ngood,xknew4,value4,niter,
C     *                              ndone4,1.0D-10,nmucur,
C     *                              g4,xkpts,r0,w0,xkhigh,xklow,jj,
C     *                              nprint)
C                     write(6,785) xkpts(i),xmupts(j),neps,ngood,
C     *                            xknew4,value4,ndone4,jj
C                     pause
                  endif
                  ndone(i,j,4,jj) = ndone4
                  grad(i,j,4,jj) = g4
                  optk2(i,j,4,jj) = xknew4
                  y2(i,j,4,jj) = value4
C                  pause

C                  write(6,398) i,j,jj,xkpts(i),xmupts(j),y2(i,j,1,jj),
C     *                         y2(i,j,2,jj),y2(i,j,3,jj),y2(i,j,4,jj)
C                  write(6,398) i,j,jj,xkpts(i),xmupts(j),
C     *                         optk2(i,j,1,jj),optk2(i,j,2,jj),
C     *                         optk2(i,j,3,jj),optk2(i,j,4,jj)
C398               format(i4,2i3,f10.3,f8.3,4f12.5)

 74            continue
 73         continue
 10      continue

         vmxdif = 0.0D+00
         dmxdif = 0.0D+00
         do 14 i = 1,nkpts
            do 15 j = 1,nmupts
               do 16 kk = 1,4
                  do 17 jj = 1,nbetas
                     vdiff = dabs(y2(i,j,kk,jj)-y(i,j,kk,jj))/
     *                       dabs(y(i,j,kk,jj))
                     if (vdiff .gt. vmxdif) vmxdif = vdiff
                     vdiffs(i,j,kk,jj) = vdiff
                     ddiff = dabs(optk2(i,j,kk,jj)-optk(i,j,kk,jj))/
     *                       dabs(optk(i,j,kk,jj))
                     if (ddiff .gt. dmxdif) dmxdif = ddiff
                     ddiffs(i,j,kk,jj) = ddiff
 17               continue
 16            continue
 15         continue
 14      continue

         if (multpl(k,40) .eq. 1) then
            write(6,341) nloop,k,mxiter,vmxdif,dmxdif
341         format(3i10,2f18.10)
            write(2,341) nloop,k,mxiter,vmxdif,dmxdif
C            pause
         endif

 30   continue

      if (nloop .gt. 0) then

         open(7,file='heter40w.vdp',form='formatted')
         do 656 j = 1,nbetas
            do 657 i = 1,nkpts
               write(7,560) xkpts(i),v1dp(i,3,j,1),v1dp(i,3,j,2),
     *                      v2dp(i,3,j,1),v2dp(i,3,j,2),
     *                      v3dp(i,3,j,1),v3dp(i,3,j,2),
     *                      v4dp(i,3,j,1),v4dp(i,3,j,2)
560            format(f8.3,8f20.10)
657         continue
656      continue
         close(7)

         open(7,file='heter40w.in2',form='formatted')
         do 759 k = 1,nbetas
            do 756 i = 1,nkpts
               do 757 j = 1,nmupts
                  write(7,445) xkpts(i),xmupts(j),y2(i,j,1,k),
     *                         y2(i,j,2,k),y2(i,j,3,k),y2(i,j,4,k)
757            continue
756         continue
759      continue
         close(7)

         open(1,file='heter40w.res',form='formatted')

         do 372 k = 1,nbetas
            do 72 i = 1,nkpts
               do 96 j = 1,nmupts
                  write(1,523) xkpts(i),xmupts(j),
     *                         optk2(i,j,1,k),optk2(i,j,2,k),
     *                         optk2(i,j,3,k),optk2(i,j,4,k),
     *                         y2(i,j,1,k),y2(i,j,2,k),y2(i,j,3,k),
     *                         y2(i,j,4,k)
523               format(2f8.3,8f12.6)
 96            continue
 72         continue
372      continue

         close(1)

         open(10,file='heter40w.dif',form='formatted')
         do 703 k = 1,nbetas
            do 715 i = 1,nkpts
               do 716 j = 1,nmupts
                  write(10,456) xkpts(i),xmupts(j),vdiffs(i,j,1,k),
     *                          vdiffs(i,j,2,k),vdiffs(i,j,3,k),
     *                          vdiffs(i,j,4,k),ddiffs(i,j,1,k),
     *                          ddiffs(i,j,2,k),ddiffs(i,j,3,k),
     *                          ddiffs(i,j,4,k)
456               format(2f8.3,8f10.7)
716            continue
715         continue
703      continue
         close(10)

         open(10,file='heter40w.tm4',form='formatted')
         open(11,file='heter40w.tm5',form='formatted')
         do 298 k = 1,nbetas
            do 98 i = 1,nkpts
               do 99 j = 1,nmupts
                  write(10,565) xkpts(i),xmupts(j),
     *                          ndone(i,j,1,k),ndone(i,j,2,k),
     *                          ndone(i,j,3,k),ndone(i,j,4,k)
565               format(2f8.3,4i6)
                  write(11,555) xkpts(i),xmupts(j),
     *                          grad(i,j,1,k),grad(i,j,2,k),
     *                          grad(i,j,3,k),grad(i,j,4,k)
555               format(2f8.3,4f13.8)
 99            continue
 98         continue
298      continue
         close(10)
         close(11)

      endif

C      goto 2
C      pause

      open(1,file='points.kgd',form='formatted')
      do 49 i = 1,kgrid
         read(1,149) xkgpts(i)
149      format(f8.3)
 49   continue
      close(1)

      xinc = (xmghgh-xmglow)/dble(real(mgrid-1))
      xmgpts(1) = xmglow
      do 44 i = 2,mgrid
         xmgpts(i) = xmgpts(i-1) + xinc
 44   continue

      do 45 j = 1,mgrid

         if (multpl(j,5) .eq. 1) then
            write(6,513) j,xmgpts(j)
513         format(' Working on mgrid point number ',i6,f15.8)
         endif

         xmu = xmgpts(j)

         agghrs = hfix*(one-unempg)
         r0cur(1) = zgood*alpha*(xmu**(alpha-one))*
     *                          (agghrs**(one-alpha))
         w0cur(1) = zgood*(one-alpha)*(xmu**alpha)*
     *                                (agghrs**(-alpha))
         agghrs = hfix*(one-unempb)
         r0cur(2) = zbad*alpha*(xmu**(alpha-one))*
     *                         (agghrs**(one-alpha))
         w0cur(2) = zbad*(one-alpha)*(xmu**alpha)*
     *                               (agghrs**(-alpha))

         call inter2(xkpts,xmupts,y2,a0,a1,c0,c1,xmu)

         do 46 i = 1,kgrid

           xk = xkgpts(i)

           do 246 k = 1,nbetas

               neps = 1
               ngood = 1
               if (xk .lt. xkcut1) then
                  call calop2(xk,neps,ngood,xknew,value,niter,
     *                        mdone,1.0D-10,
     *                        g,xkpts,r0cur,w0cur,xkhigh,xklow,k)
               else
                  if (i .eq. 1) then
                     startk = xk
                  else
                     startk = optkgd(i-1,j,1,k)
                  endif
                  call altop2(xk,neps,ngood,xknew,value,niter,
     *                        mdone,1.0D-10,
     *                        startk,grad,xkpts,r0cur,w0cur,k)
                  if (mdone .eq. 0)
     *               call calop2(xk,neps,ngood,xknew,value,niter,
     *                           mdone,1.0D-10,
     *                           g,xkpts,r0cur,w0cur,xkhigh,xklow,k)
               endif
               optkgd(i,j,1,k) = xknew

               neps = 1
               ngood = 2
               if (xk .lt. xkcut1) then
                  call calop2(xk,neps,ngood,xknew,value,niter,
     *                        mdone,1.0D-10,
     *                        g,xkpts,r0cur,w0cur,xkhigh,xklow,k)
               else
                  if (i .eq. 1) then
                     startk = xk
                  else
                     startk = optkgd(i-1,j,2,k)
                  endif
                  call altop2(xk,neps,ngood,xknew,value,niter,
     *                        mdone,1.0D-10,
     *                        startk,grad,xkpts,r0cur,w0cur,k)
                  if (mdone .eq. 0)
     *               call calop2(xk,neps,ngood,xknew,value,niter,
     *                           mdone,1.0D-10,
     *                           g,xkpts,r0cur,w0cur,xkhigh,xklow,k)
               endif
               optkgd(i,j,2,k) = xknew

               neps = 0
               ngood = 1
               if (xk .lt. xkcut2) then
                  call calop2(xk,neps,ngood,xknew,value,niter,
     *                        mdone,1.0D-10,
     *                        g,xkpts,r0cur,w0cur,xkhigh,xklow,k)
               else
                  if (i .eq. 1) then
                     startk = xk
                  else
                     startk = optkgd(i-1,j,3,k)
                  endif
                  call altop2(xk,neps,ngood,xknew,value,niter,
     *                        mdone,1.0D-10,
     *                        startk,grad,xkpts,r0cur,w0cur,k)
                  if (mdone .eq. 0)
     *               call calop2(xk,neps,ngood,xknew,value,niter,
     *                           mdone,1.0D-10,
     *                           g,xkpts,r0cur,w0cur,xkhigh,xklow,k)
               endif
               optkgd(i,j,3,k) = xknew

               neps = 0
               ngood = 2
               if (xk .lt. xkcut2) then
                  call calop2(xk,neps,ngood,xknew,value,niter,
     *                        mdone,1.0D-10,
     *                        g,xkpts,r0cur,w0cur,xkhigh,xklow,k)
               else
                  if (i .eq. 1) then
                     startk = xk
                  else
                     startk = optkgd(i-1,j,4,k)
                  endif
                  call altop2(xk,neps,ngood,xknew,value,niter,
     *                        mdone,1.0D-10,
     *                        startk,grad,xkpts,r0cur,w0cur,k)
                  if (mdone .eq. 0)
     *               call calop2(xk,neps,ngood,xknew,value,niter,
     *                           mdone,1.0D-10,
     *                           g,xkpts,r0cur,w0cur,xkhigh,xklow,k)
               endif
               optkgd(i,j,4,k) = xknew

246         continue
 46      continue
 45   continue

      open(9,file='heter40w.kgd',form='formatted')
      do 876 i = 1,kgrid
         do 877 j = 1,mgrid
            write(9,960) xkgpts(i),xmgpts(j),
     *                   optkgd(i,j,1,1),optkgd(i,j,2,1),
     *                   optkgd(i,j,3,1),optkgd(i,j,4,1),
     *                   optkgd(i,j,1,2),optkgd(i,j,2,2),
     *                   optkgd(i,j,3,2),optkgd(i,j,4,2),
     *                   optkgd(i,j,1,3),optkgd(i,j,2,3),
     *                   optkgd(i,j,3,3),optkgd(i,j,4,3)
960         format(2f8.3,12f13.7)
C960         format(2f8.3,12f11.6)
877      continue
876   continue
      close(9)

  3   if (nkgd .eq. 1) then
         write(6,299)
299      format(' Reading optkgd....')
         open(9,file='heter40w.kgd',form='formatted')
         do 576 i = 1,kgrid
            do 577 j = 1,mgrid
               read(9,960) xkgpts(i),xmgpts(j),
     *                     optkgd(i,j,1,1),optkgd(i,j,2,1),
     *                     optkgd(i,j,3,1),optkgd(i,j,4,1),
     *                     optkgd(i,j,1,2),optkgd(i,j,2,2),
     *                     optkgd(i,j,3,2),optkgd(i,j,4,2),
     *                     optkgd(i,j,1,3),optkgd(i,j,2,3),
     *                     optkgd(i,j,3,3),optkgd(i,j,4,3)
577         continue
576      continue
         close(9)
         write(6,297)
297      format(' Done reading optkgd.')
      endif

      call dosim(seed1,seed2,xkbd,
     *           xkmean,xkstd,xkmax,xkmin,
     *           ymean,ystd,ymax,ymin,
     *           cmean,cstd,cmax,cmin,
     *           xmean,xstd,xmax,xmin,
     *           hmean,hstd,hmax,hmin,
     *           xksdmn,xksdsd,xksdmx,xksdmi,
C     *           xkskmn,xksksd,xkskmx,xkskmi,
C     *           xkkrmn,xkkrsd,xkkrmx,xkkrmi,
     *           x1mxmn,x1mxsd,x1mxmx,x1mxmi,
     *           x2mxmn,x2mxsd,x2mxmx,x2mxmi,
     *           x3mxmn,x3mxsd,x3mxmx,x3mxmi,
     *           x1mimn,x1misd,x1mimx,x1mimi,
     *           x2mimn,x2misd,x2mimx,x2mimi,
     *           x3mimn,x3misd,x3mimx,x3mimi,
     *           xgmn,xgsd,xgmx,xgmi,
     *           topmn,topsd,topmx,topmi,
     *           pct0mn,pct0sd,pct0mx,pct0mi,
     *           pctcmn,pctcsd,pctcmx,pctcmi,
     *           bhat1,bhat2,rsqa,rsqb,
     *           sderra,sderrb)

      a0new = bhat1(1,1)
      a1new = bhat1(2,1)
      xmua = dexp(a0new/(one-a1new))
      xhrsa = hfix
      c0new = bhat2(1,1)
      c1new = bhat2(2,1)
      xmub = dexp(c0new/(one-c1new))
      xhrsb = hfix

      write(2,265) brho1,brho2,brho3
      write(2,265) betas(1),betas(2),betas(3)
      write(2,265) prob(1),prob(2),prob(3)
      write(2,581) xkmean,ymean,cmean,xmean,hmean
581   format(5f15.8)
      write(2,581) xkstd,ystd,cstd,xstd,hstd
      write(2,581) xkstd/xkmean,ystd/ymean,cstd/cmean,xstd/xmean,
     *             hstd/hmean
      write(2,581) xkmax,ymax,cmax,xmax,hmax
      write(2,581) xkmin,ymin,cmin,xmin,hmin
      write(2,582) a0,a1
583   format(4f15.8)
      write(2,582) a0new,a1new
      write(2,582) rsqa,sderra
      write(2,582) xmua,xhrsa
582   format(2f15.8)
      write(2,582) c0,c1
      write(2,582) c0new,c1new
      write(2,582) rsqb,sderrb
      write(2,582) xmub,xhrsb
      a0new = sfac*a0new + (one-sfac)*a0
      a1new = sfac*a1new + (one-sfac)*a1
      c0new = sfac*c0new + (one-sfac)*c0
      c1new = sfac*c1new + (one-sfac)*c1
      write(2,582) a0new,a1new
      write(2,582) c0new,c1new
      write(2,581) xksdmn,xksdsd,xksdmx,xksdmi
C      write(2,581) xkskmn,xksksd,xkskmx,xkskmi
C      write(2,581) xkkrmn,xkkrsd,xkkrmx,xkkrmi
      write(2,581) x1mxmn,x1mxsd,x1mxmx,x1mxmi
      write(2,581) x2mxmn,x2mxsd,x2mxmx,x2mxmi
      write(2,581) x3mxmn,x3mxsd,x3mxmx,x3mxmi
      write(2,581) x1mimn,x1misd,x1mimx,x1mimi
      write(2,581) x2mimn,x2misd,x2mimx,x2mimi
      write(2,581) x3mimn,x3misd,x3mimx,x3mimi
      write(2,581) xgmn,xgsd,xgmx,xgmi
      do 356 i = 1,ntop
         write(2,581) topmn(i),topsd(i),topmx(i),topmi(i)
356   continue
      write(2,581) pct0mn,pct0sd,pct0mx,pct0mi
      write(2,581) pctcmn,pctcsd,pctcmx,pctcmi

      write(6,582) a0,a1
      write(6,582) a0new,a1new
      write(6,582) c0,c1
      write(6,582) c0new,c1new
      write(6,584) xkmean,rsqa,rsqb,x3mxmx
584   format(4f15.8)

      open(1,file='agglaw.out',form='formatted')
      write(1,582) a0new,a1new
      write(1,582) c0new,c1new
      close(1)

      goto 1

  2   close(2)

      end

      subroutine calopt(xk,neps,ngood,ynew,value,niter,ndone,toler,
     *                  nmucur,grad,xkpts,r0,w0,ykhigh,yklow,nbeta,
     *                  nprint)
      implicit real*8 (a-h,o-z)
      parameter (mxiter=100,hfix=0.3271D+00,
     *           nkpts=132,nmupts=4,delta=0.025D+00,nbetas=3,
     *           whome=0.07D+00)
      common/cpr/pr(4,4),prbeta(nbetas,nbetas),prob(nbetas)
      common/cbetas/betas(nbetas)
      common/cvs/v1(nkpts,nmupts,nbetas,2),v2(nkpts,nmupts,nbetas,2),
     *           v3(nkpts,nmupts,nbetas,2),v4(nkpts,nmupts,nbetas,2),
     *           v1dp(nkpts,nmupts,nbetas,2),
     *           v2dp(nkpts,nmupts,nbetas,2),
     *           v3dp(nkpts,nmupts,nbetas,2),
     *           v4dp(nkpts,nmupts,nbetas,2)
      dimension xkpts(nkpts),r0(nmupts,2),w0(nmupts,2),
     *          y1s(nbetas),y2s(nbetas),y3s(nbetas),y4s(nbetas),
     *          y1ps(nbetas),y2ps(nbetas),y3ps(nbetas),y4ps(nbetas)

      one = 1.0D+00

      beta = betas(nbeta)

      ndone = 0
      niter = 0

      if (neps .eq. 1) then
         y = r0(nmucur,ngood)*xk + w0(nmucur,ngood)*hfix
      else
         y = r0(nmucur,ngood)*xk + whome
      endif

      xklow = yklow
      xkhigh = ykhigh

      do 10 i = 1,mxiter

         xkcur = (xkhigh+xklow)/2.0D+00

         x = xkcur - (one-delta)*xk
         if (x .ge. y) then
            xkhigh = xkcur
            goto 10
         endif

         c = y - x

         do 87 j = 1,nbetas
            call splint(xkpts,v1,v1dp,xkcur,y1,y1p,y1dp,nmucur,ngood,
     *                  j,khi,klo,0,0,1,0)
            y1ps(j) = y1p
            call splint(xkpts,v2,v2dp,xkcur,y2,y2p,y2dp,nmucur,ngood,
     *                  j,khi,klo,1,0,1,0)
            y2ps(j) = y2p
            call splint(xkpts,v3,v3dp,xkcur,y3,y3p,y3dp,nmucur,ngood,
     *                  j,khi,klo,1,0,1,0)
            y3ps(j) = y3p
            call splint(xkpts,v4,v4dp,xkcur,y4,y4p,y4dp,nmucur,ngood,
     *                  j,khi,klo,1,0,1,0)
            y4ps(j) = y4p
 87      continue

         y1p = 0.0D+00
         y2p = 0.0D+00
         y3p = 0.0D+00
         y4p = 0.0D+00
         do 81 j = 1,nbetas
            y1p = y1p + prbeta(j,nbeta)*y1ps(j)
            y2p = y2p + prbeta(j,nbeta)*y2ps(j)
            y3p = y3p + prbeta(j,nbeta)*y3ps(j)
            y4p = y4p + prbeta(j,nbeta)*y4ps(j)
 81      continue

         upc = -one/c
         if ( (neps .eq. 1) .and. (ngood .eq. 1) ) then
            grad = upc + beta*(y1p*pr(1,1) + y2p*pr(2,1) +
     *                         y3p*pr(3,1) + y4p*pr(4,1))
         elseif ( (neps .eq. 1) .and. (ngood .eq. 2) ) then
            grad = upc + beta*(y1p*pr(1,2) + y2p*pr(2,2) +
     *                         y3p*pr(3,2) + y4p*pr(4,2))
         elseif ( (neps .eq. 0) .and. (ngood .eq. 1) ) then
            grad = upc + beta*(y1p*pr(1,3) + y2p*pr(2,3) +
     *                         y3p*pr(3,3) + y4p*pr(4,3))
         else
            grad = upc + beta*(y1p*pr(1,4) + y2p*pr(2,4) +
     *                         y3p*pr(3,4) + y4p*pr(4,4))
         endif

         if (grad .ge. 0.0D+00) then
            xklow = xkcur
         else
            xkhigh = xkcur
         endif

         if (nprint .eq. 1) then
            write(6,102) i,xkcur,xklow,xkhigh,grad
102         format(i5,4f15.8,'    grad')
            pause
         endif

         if ((xkhigh-xklow) .le. toler) then
            niter = i
            if (dabs(grad) .le. 1.0D-02) ndone = 1
            goto 999
         endif

 10   continue

999   xkcur = (xklow + xkhigh)/2.0D+00
      x = xkcur - (one-delta)*xk
      c = y - x
      if (c .le. 0.0D+00) then
         write(6,777) c
777      format(' Failure in calopt: negative c',f15.8)
         pause
      endif
      utilc = dlog(c)
      do 78 j = 1,nbetas
         call splint(xkpts,v1,v1dp,xkcur,y1,y1p,y1dp,nmucur,ngood,
     *               j,khi,klo,0,1,0,0)
         y1s(j) = y1
         call splint(xkpts,v2,v2dp,xkcur,y2,y2p,y2dp,nmucur,ngood,
     *               j,khi,klo,1,1,0,0)
         y2s(j) = y2
         call splint(xkpts,v3,v3dp,xkcur,y3,y3p,y3dp,nmucur,ngood,
     *               j,khi,klo,1,1,0,0)
         y3s(j) = y3
         call splint(xkpts,v4,v4dp,xkcur,y4,y4p,y4dp,nmucur,ngood,
     *               j,khi,klo,1,1,0,0)
         y4s(j) = y4
 78   continue
      y1 = 0.0D+00
      y2 = 0.0D+00
      y3 = 0.0D+00
      y4 = 0.0D+00
      do 79 j = 1,nbetas
         y1 = y1 + prbeta(j,nbeta)*y1s(j)
         y2 = y2 + prbeta(j,nbeta)*y2s(j)
         y3 = y3 + prbeta(j,nbeta)*y3s(j)
         y4 = y4 + prbeta(j,nbeta)*y4s(j)
 79   continue
      if ( (neps .eq. 1) .and. (ngood .eq. 1) ) then
         fval = utilc + beta*(y1*pr(1,1) + y2*pr(2,1) +
     *                        y3*pr(3,1) + y4*pr(4,1))
      elseif ( (neps .eq. 1) .and. (ngood .eq. 2) ) then
         fval = utilc + beta*(y1*pr(1,2) + y2*pr(2,2) +
     *                        y3*pr(3,2) + y4*pr(4,2))
      elseif ( (neps .eq. 0) .and. (ngood .eq. 1) ) then
         fval = utilc + beta*(y1*pr(1,3) + y2*pr(2,3) +
     *                        y3*pr(3,3) + y4*pr(4,3))
      else
         fval = utilc + beta*(y1*pr(1,4) + y2*pr(2,4) +
     *                        y3*pr(3,4) + y4*pr(4,4))
      endif

      ynew = xkcur
      value = fval

      return
      end

      subroutine altopt(xk,neps,ngood,ynew,value,niter,ndone,toler,
     *                  nmucur,startk,grad,xkpts,r0,w0,nbeta)
      implicit real*8 (a-h,o-z)
      parameter (mxiter=100,hfix=0.3271D+00,
     *           nkpts=132,nmupts=4,delta=0.025D+00,nbetas=3,
     *           whome=0.07D+00)
      common/cpr/pr(4,4),prbeta(nbetas,nbetas),prob(nbetas)
      common/cbetas/betas(nbetas)
      common/cvs/v1(nkpts,nmupts,nbetas,2),v2(nkpts,nmupts,nbetas,2),
     *           v3(nkpts,nmupts,nbetas,2),v4(nkpts,nmupts,nbetas,2),
     *           v1dp(nkpts,nmupts,nbetas,2),
     *           v2dp(nkpts,nmupts,nbetas,2),
     *           v3dp(nkpts,nmupts,nbetas,2),
     *           v4dp(nkpts,nmupts,nbetas,2)
      dimension xkpts(nkpts),r0(nmupts,2),w0(nmupts,2),
     *          y1s(nbetas),y2s(nbetas),y3s(nbetas),y4s(nbetas),
     *          y1ps(nbetas),y2ps(nbetas),y3ps(nbetas),y4ps(nbetas),
     *          y1dps(nbetas),y2dps(nbetas),y3dps(nbetas),y4dps(nbetas)

      one = 1.0D+00

      beta = betas(nbeta)

      ndone = 0
      niter = 0

      if (neps .eq. 1) then
         y = r0(nmucur,ngood)*xk + w0(nmucur,ngood)*hfix
      else
         y = r0(nmucur,ngood)*xk + whome
      endif

      xk2p = startk

      do 10 i = 1,mxiter

         xk2 = xk2p
         x = xk2 - (one-delta)*xk
         c = y - x

         do 87 j = 1,nbetas
            call splint(xkpts,v1,v1dp,xk2,y1,y1p,y1dp,nmucur,ngood,
     *                  j,khi,klo,0,0,1,1)
            y1ps(j) = y1p
            y1dps(j) = y1dp
            call splint(xkpts,v2,v2dp,xk2,y2,y2p,y2dp,nmucur,ngood,
     *                  j,khi,klo,1,0,1,1)
            y2ps(j) = y2p
            y2dps(j) = y2dp
            call splint(xkpts,v3,v3dp,xk2,y3,y3p,y3dp,nmucur,ngood,
     *                  j,khi,klo,1,0,1,1)
            y3ps(j) = y3p
            y3dps(j) = y3dp
            call splint(xkpts,v4,v4dp,xk2,y4,y4p,y4dp,nmucur,ngood,
     *                  j,khi,klo,1,0,1,1)
            y4ps(j) = y4p
            y4dps(j) = y4dp
 87      continue

         y1p = 0.0D+00
         y2p = 0.0D+00
         y3p = 0.0D+00
         y4p = 0.0D+00
         y1dp = 0.0D+00
         y2dp = 0.0D+00
         y3dp = 0.0D+00
         y4dp = 0.0D+00
         do 81 j = 1,nbetas
            y1p = y1p + prbeta(j,nbeta)*y1ps(j)
            y2p = y2p + prbeta(j,nbeta)*y2ps(j)
            y3p = y3p + prbeta(j,nbeta)*y3ps(j)
            y4p = y4p + prbeta(j,nbeta)*y4ps(j)
            y1dp = y1dp + prbeta(j,nbeta)*y1dps(j)
            y2dp = y2dp + prbeta(j,nbeta)*y2dps(j)
            y3dp = y3dp + prbeta(j,nbeta)*y3dps(j)
            y4dp = y4dp + prbeta(j,nbeta)*y4dps(j)
 81      continue

         upc = -one/c
         if ( (neps .eq. 1) .and. (ngood .eq. 1) ) then
            grad = upc + beta*(y1p*pr(1,1) + y2p*pr(2,1) +
     *                         y3p*pr(3,1) + y4p*pr(4,1))
         elseif ( (neps .eq. 1) .and. (ngood .eq. 2) ) then
            grad = upc + beta*(y1p*pr(1,2) + y2p*pr(2,2) +
     *                         y3p*pr(3,2) + y4p*pr(4,2))
         elseif ( (neps .eq. 0) .and. (ngood .eq. 1) ) then
            grad = upc + beta*(y1p*pr(1,3) + y2p*pr(2,3) +
     *                         y3p*pr(3,3) + y4p*pr(4,3))
         else
            grad = upc + beta*(y1p*pr(1,4) + y2p*pr(2,4) +
     *                         y3p*pr(3,4) + y4p*pr(4,4))
         endif

         udpc = -one/(c*c)
         if ( (neps .eq. 1) .and. (ngood .eq. 1) ) then
            hess = udpc + beta*(y1dp*pr(1,1) + y2dp*pr(2,1) +
     *                          y3dp*pr(3,1) + y4dp*pr(4,1))
         elseif ( (neps .eq. 1) .and. (ngood .eq. 2) ) then
            hess = udpc + beta*(y1dp*pr(1,2) + y2dp*pr(2,2) +
     *                          y3dp*pr(3,2) + y4dp*pr(4,2))
         elseif ( (neps .eq. 0) .and. (ngood .eq. 1) ) then
            hess = udpc + beta*(y1dp*pr(1,3) + y2dp*pr(2,3) +
     *                          y3dp*pr(3,3) + y4dp*pr(4,3))
         else
            hess = udpc + beta*(y1dp*pr(1,4) + y2dp*pr(2,4) +
     *                          y3dp*pr(3,4) + y4dp*pr(4,4))
         endif

         xk2p = xk2 - grad/hess

         if (dabs(xk2p-xk2) .le. toler) then
            ndone = 1
            niter = i
            goto 999
         endif

 10   continue

999   if (c .gt. 0.0D+00) then

         utilc = dlog(c)
         do 78 j = 1,nbetas
            call splint(xkpts,v1,v1dp,xk2,y1,y1p,y1dp,nmucur,ngood,
     *                  j,khi,klo,0,1,0,0)
            y1s(j) = y1
            call splint(xkpts,v2,v2dp,xk2,y2,y2p,y2dp,nmucur,ngood,
     *                  j,khi,klo,1,1,0,0)
            y2s(j) = y2
            call splint(xkpts,v3,v3dp,xk2,y3,y3p,y3dp,nmucur,ngood,
     *                  j,khi,klo,1,1,0,0)
            y3s(j) = y3
            call splint(xkpts,v4,v4dp,xk2,y4,y4p,y4dp,nmucur,ngood,
     *                  j,khi,klo,1,1,0,0)
            y4s(j) = y4
 78      continue
         y1 = 0.0D+00
         y2 = 0.0D+00
         y3 = 0.0D+00
         y4 = 0.0D+00
         do 79 j = 1,nbetas
            y1 = y1 + prbeta(j,nbeta)*y1s(j)
            y2 = y2 + prbeta(j,nbeta)*y2s(j)
            y3 = y3 + prbeta(j,nbeta)*y3s(j)
            y4 = y4 + prbeta(j,nbeta)*y4s(j)
 79      continue
         if ( (neps .eq. 1) .and. (ngood .eq. 1) ) then
            fval = utilc + beta*(y1*pr(1,1) + y2*pr(2,1) +
     *                           y3*pr(3,1) + y4*pr(4,1))
         elseif ( (neps .eq. 1) .and. (ngood .eq. 2) ) then
            fval = utilc + beta*(y1*pr(1,2) + y2*pr(2,2) +
     *                           y3*pr(3,2) + y4*pr(4,2))
         elseif ( (neps .eq. 0) .and. (ngood .eq. 1) ) then
            fval = utilc + beta*(y1*pr(1,3) + y2*pr(2,3) +
     *                           y3*pr(3,3) + y4*pr(4,3))
         else
            fval = utilc + beta*(y1*pr(1,4) + y2*pr(2,4) +
     *                           y3*pr(3,4) + y4*pr(4,4))
         endif

         ynew = xk2
         value = fval

      else
         ndone = 0
      endif

      return
      end

      subroutine dosim(seed1,seed2,xkbd,
     *                 xkmean,xkstd,xkmax,xkmin,
     *                 ymean,ystd,ymax,ymin,
     *                 cmean,cstd,cmax,cmin,
     *                 xmean,xstd,xmax,xmin,
     *                 hmean,hstd,hmax,hmin,
     *                 xksdmn,xksdsd,xksdmx,xksdmi,
C     *                 xkskmn,xksksd,xkskmx,xkskmi,
C     *                 xkkrmn,xkkrsd,xkkrmx,xkkrmi,
     *                 x1mxmn,x1mxsd,x1mxmx,x1mxmi,
     *                 x2mxmn,x2mxsd,x2mxmx,x2mxmi,
     *                 x3mxmn,x3mxsd,x3mxmx,x3mxmi,
     *                 x1mimn,x1misd,x1mimx,x1mimi,
     *                 x2mimn,x2misd,x2mimx,x2mimi,
     *                 x3mimn,x3misd,x3mimx,x3mimi,
     *                 xgmn,xgsd,xgmx,xgmi,
     *                 topmn,topsd,topmx,topmi,
     *                 pct0mn,pct0sd,pct0mx,pct0mi,
     *                 pctcmn,pctcsd,pctcmx,pctcmi,
     *                 bhat1,bhat2,rsqa,rsqb,
     *                 sderra,sderrb)
      implicit real*8 (a-h,o-z)
      parameter (durug=1.5D+00,durub=2.5D+00,unempg=0.04D+00,
     *           unempb=0.1D+00,durgd=8.0D+00,durbd=8.0D+00,
     *           ncapt=11000,ndrop=1000,nzpts=2,nkpts=132,
     *           kgrid=278,nbetas=3,nagnt1=3000,nagnt2=24000,
     *           nagnt3=3000,nagent=nagnt1+nagnt2+nagnt3,
     *           mgrid=30,alpha=0.36D+00,delta=0.025D+00,
     *           zgood=1.01D+00,zbad=0.99D+00,hfix=0.3271D+00,
     *           xkbor=-2.4D+00,toler=1.0D-05,nlzpts=201,ntop=9,
     *           mxswch=15)
      common/cgrid/xkgpts(kgrid),xmgpts(mgrid)
      common/ctop/toppct(ntop)
      common/clzpts/xlzpts(nlzpts),xlzdat(nlzpts)
      common/ccoef/coefk(kgrid-1,mgrid-1,4,nbetas,4),
     *             optkgd(kgrid,mgrid,4,nbetas)
      common/cpr/pr(4,4),prbeta(nbetas,nbetas),prob(nbetas)
      dimension uniz(ncapt),nz(ncapt+1),xkbar(ncapt+1),
     *          preps(nzpts,nzpts),pctemp(nbetas),
     *          xkistd(ncapt+1),xx(4,4),yy(4,1),xxinv(4,4),temp(4,1),
C     *          xkiskw(ncapt+1),xkikrt(ncapt+1),
     *          hbar(ncapt),xbar(ncapt),cbar(ncapt),
     *          ybar(ncapt),tmpdat(ncapt-ndrop),xdata(ncapt,2),
     *          ydata(ncapt),bhat1(2,1),bhat2(2,1),
     *          ylzpts(nlzpts),ginis(ncapt+1),pctcon(ncapt+1),
     *          top(ntop),topdat(ncapt+1,ntop),topmn(ntop),
     *          topsd(ntop),topmx(ntop),topmi(ntop),pct0(ncapt+1),
     *          xk11(nagnt1),xk21(nagnt2),xk31(nagnt3),
     *          xk12(nagnt1),xk22(nagnt2),xk32(nagnt3),
     *          neps11(nagnt1),neps21(nagnt2),neps31(nagnt3),
     *          neps12(nagnt1),neps22(nagnt2),neps32(nagnt3),
     *          uniep1(nagnt1),uniep2(nagnt2),uniep3(nagnt3),
     *          nums1(mxswch),nums2(2*mxswch),nums3(mxswch),
     *          ntmp1(mxswch),xktmp1(mxswch),ntmp3(mxswch),
     *          xktmp3(mxswch),
     *          xkmax1(ncapt+1),xkmax2(ncapt+1),xkmax3(ncapt+1),
     *          xkmin1(ncapt+1),xkmin2(ncapt+1),xkmin3(ncapt+1),
     *          hist(kgrid-1),fractn(kgrid-1),xpt(kgrid),ypt(kgrid),
     *          a(kgrid),b(kgrid),c(kgrid),r(kgrid),ydp(kgrid)

      one = 1.0D+00

      xnagnt = dble(real(nagent))

      xmginc = xmgpts(2) - xmgpts(1)
      xmglow = xmgpts(1)
      xmghgh = xmgpts(mgrid)

      wght1 = prob(1)*xnagnt/dble(real(nagnt1))
      wght2 = prob(2)*xnagnt/dble(real(nagnt2))
      wght3 = prob(3)*xnagnt/dble(real(nagnt3))

      swtch1 = prbeta(2,1)*dble(real(nagnt1))
      swtch2 = prbeta(3,2)*dble(real(nagnt2))
      swtch3 = prbeta(2,3)*dble(real(nagnt3))

      xx(1,1) = one
      xx(2,1) = one
      xx(3,1) = one
      xx(4,1) = one

      write(6,298)
298   format(' Calculating coefficients....')

      do 330 i = 1,kgrid-1
         do 331 j = 1,mgrid-1
            do 332 k = 1,4
               do 333 kk = 1,nbetas
                  if ( (dabs(optkgd(i,j,k,kk)-xkbor)
     *                  .le. toler) .or.
     *                 (dabs(optkgd(i,j+1,k,kk)-xkbor)
     *                  .le. toler) .or.
     *                 (dabs(optkgd(i+1,j,k,kk)-xkbor)
     *                  .le. toler) .or.
     *                 (dabs(optkgd(i+1,j+1,k,kk)-xkbor)
     *                  .le. toler) )
     *            then
                     coefk(i,j,k,kk,1) = xkbor
                     coefk(i,j,k,kk,2) = 0.0D+00
                     coefk(i,j,k,kk,3) = 0.0D+00
                     coefk(i,j,k,kk,4) = 0.0D+00
                  else
                     xx(1,2) = xkgpts(i)
                     xx(1,3) = xmgpts(j)
                     xx(1,4) = xx(1,2)*xx(1,3)
                     yy(1,1) = optkgd(i,j,k,kk)
                     xx(2,2) = xkgpts(i)
                     xx(2,3) = xmgpts(j+1)
                     xx(2,4) = xx(2,2)*xx(2,3)
                     yy(2,1) = optkgd(i,j+1,k,kk)
                     xx(3,2) = xkgpts(i+1)
                     xx(3,3) = xmgpts(j)
                     xx(3,4) = xx(3,2)*xx(3,3)
                     yy(3,1) = optkgd(i+1,j,k,kk)
                     xx(4,2) = xkgpts(i+1)
                     xx(4,3) = xmgpts(j+1)
                     xx(4,4) = xx(4,2)*xx(4,3)
                     yy(4,1) = optkgd(i+1,j+1,k,kk)
                     call inv4(xx,xxinv,nsing)
                     call mult2(xxinv,yy,temp,4,4,1)
                     coefk(i,j,k,kk,1) = temp(1,1)
                     coefk(i,j,k,kk,2) = temp(2,1)
                     coefk(i,j,k,kk,3) = temp(3,1)
                     coefk(i,j,k,kk,4) = temp(4,1)
                  endif
333            continue
332         continue
331      continue
330   continue

      write(6,299)
299   format(' Done calculating coefficients.')

      pgg00 = (durug-one)/durug
      pbb00 = (durub-one)/durub
      pbg00 = 1.25D+00*pbb00
      pgb00 = 0.75D+00*pgg00
      pgg01 = (unempg - unempg*pgg00)/(one-unempg)
      pbb01 = (unempb - unempb*pbb00)/(one-unempb)
      pbg01 = (unempb - unempg*pbg00)/(one-unempg)
      pgb01 = (unempg - unempb*pgb00)/(one-unempb)

      pgg = (durgd-one)/durgd
      pgb = one - (durbd-one)/durbd

      dseed = seed1
      sseed = seed2

      call ggubs(dseed,ncapt,uniz)
      nz(1) = 1
      do 93 i = 1,ncapt-1
         if (nz(i) .eq. 1) then
            if (uniz(i) .le. pgg) then
               nz(i+1) = 1
            else
               nz(i+1) = 0
            endif
         else
            if (uniz(i) .le. pgb) then
               nz(i+1) = 1
            else
               nz(i+1) = 0
            endif
         endif
 93   continue
      nz(ncapt+1) = 1

      nswch1 = 0
      nswch2 = 0
      nswch3 = 0

      open(1,file='cross.in',form='formatted')

      do 79 i = 1,nagnt1
         read(1,855) idummy,idummy,neps11(i),xk11(i)
855      format(i6,2i3,f12.6)
 79   continue

      do 73 i = nagnt1+1,nagnt1+nagnt2
         read(1,855) idummy,idummy,neps21(i-nagnt1),xk21(i-nagnt1)
 73   continue

      do 74 i = nagnt1+nagnt2+1,nagent
         read(1,855) idummy,idummy,neps31(i-nagnt1-nagnt2),
     *               xk31(i-nagnt1-nagnt2)
 74   continue

      close(1)

      call domean(xk11,xk21,xk31,nagnt1,nagnt2,nagnt3,wght1,
     *            wght2,wght3,xksum)
      call dopct(xk11,xk21,xk31,neps11,neps21,neps31,nagnt1,nagnt2,
     *           nagnt3,wght1,wght2,wght3,pctemp,p0,pcon)

      do 10 i = 1,ncapt

         do 50 j = 1,kgrid-1
            hist(j) = 0.0D+00
            fractn(j) = 0.0D+00
 50      continue

         xkbar(i) = xksum
         xkistd(i) = dstd3(xk11,xk21,xk31,nagnt1,nagnt2,nagnt3,
     *                     xkbar(i),prob(1),prob(2),prob(3))
         pct0(i) = p0
         pctcon(i) = pcon
         xkmax1(i) = dmax(xk11,nagnt1)
         xkmax2(i) = dmax(xk21,nagnt2)
         xkmax3(i) = dmax(xk31,nagnt3)
         xkmin1(i) = dmin(xk11,nagnt1)
         xkmin2(i) = dmin(xk21,nagnt2)
         xkmin3(i) = dmin(xk31,nagnt3)

         call locid2(xkbar(i),xmglow,xmghgh,xmginc,mgrid,mcor,
     *               mlow,mhigh)

         do 70 j = 1,nagnt1
            call locidx(xk11(j),xkgpts,kgrid,kcor,klow,khigh,niter)
            if (klow .eq. 1) then
               kcor = 1
            elseif (khigh .eq. 1) then
               kcor = kgrid-1
            endif
            ns = nstate(nz(i),neps11(j))
            xk12(j) = coefk(kcor,mcor,ns,1,1) +
     *                coefk(kcor,mcor,ns,1,2)*xk11(j) +
     *                coefk(kcor,mcor,ns,1,3)*xkbar(i) +
     *                coefk(kcor,mcor,ns,1,4)*xk11(j)*xkbar(i)
C            write(6,233) j,neps11(j),xk11(j),ns,kcor,mcor,xk12(j)
C233         format(2i6,f10.6,3i6,f10.6)
            if (xk12(j) .gt. xkbd) xk12(j) = xkbd
            hist(kcor) = hist(kcor) + wght1/xnagnt
            fractn(kcor) = fractn(kcor) + wght1*xk11(j)/xnagnt
 70      continue

         do 170 j = 1,nagnt2
            call locidx(xk21(j),xkgpts,kgrid,kcor,klow,khigh,niter)
            if (klow .eq. 1) then
               kcor = 1
            elseif (khigh .eq. 1) then
               kcor = kgrid-1
            endif
            ns = nstate(nz(i),neps21(j))
            xk22(j) = coefk(kcor,mcor,ns,2,1) +
     *                coefk(kcor,mcor,ns,2,2)*xk21(j) +
     *                coefk(kcor,mcor,ns,2,3)*xkbar(i) +
     *                coefk(kcor,mcor,ns,2,4)*xk21(j)*xkbar(i)
C            write(6,233) j,neps21(j),xk21(j),ns,kcor,mcor,xk22(j)
            if (xk22(j) .gt. xkbd) xk22(j) = xkbd
            hist(kcor) = hist(kcor) + wght2/xnagnt
            fractn(kcor) = fractn(kcor) + wght2*xk21(j)/xnagnt
170      continue

         do 270 j = 1,nagnt3
            call locidx(xk31(j),xkgpts,kgrid,kcor,klow,khigh,niter)
            if (klow .eq. 1) then
               kcor = 1
            elseif (khigh .eq. 1) then
               kcor = kgrid-1
            endif
            ns = nstate(nz(i),neps31(j))
            xk32(j) = coefk(kcor,mcor,ns,3,1) +
     *                coefk(kcor,mcor,ns,3,2)*xk31(j) +
     *                coefk(kcor,mcor,ns,3,3)*xkbar(i) +
     *                coefk(kcor,mcor,ns,3,4)*xk31(j)*xkbar(i)
C            write(6,233) j,neps31(j),xk31(j),ns,kcor,mcor,xk32(j)
            if (xk32(j) .gt. xkbd) xk32(j) = xkbd
            hist(kcor) = hist(kcor) + wght3/xnagnt
            fractn(kcor) = fractn(kcor) + wght3*xk31(j)/xnagnt
270      continue

         sumx = 0.0D+00
         sumy = 0.0D+00
         xpt(1) = 0.0D+00
         ypt(1) = 0.0D+00
         nhpts = 1
         do 52 j = 2,kgrid
            sumx = sumx + hist(j-1)
            sumy = sumy + fractn(j-1)/xkbar(i)
            if (hist(j-1) .gt. 1.0D-10) then
               nhpts = nhpts + 1
               xpt(nhpts) = sumx
               ypt(nhpts) = sumy
            endif
 52      continue

         call calcsd(xpt,ypt,ydp,a,b,c,r,nhpts,-one,-one)
         call integr(xpt,ypt,ydp,nhpts,xint)
         gini = one - 2.0D+00*xint

         do 97 j = 1,ntop
            call doint(xpt,ypt,ydp,one-toppct(j),yval,yp,ydp,
     *                 khi,klo,0,nhpts,1,0,0)
            top(j) = one - yval
 97      continue

         ginis(i) = gini
         do 77 j = 1,ntop
            topdat(i,j) = top(j)
 77      continue

         if ( (i .eq. 1) .or. (multpl(i,100) .eq. 1) ) then
            write(6,255) i,nz(i),pctemp(1),pctemp(2),pctemp(3),
     *                   xkbar(i),xkistd(i),pct0(i),pctcon(i),
     *                   nswch1,nswch2,nswch3,nhpts
255         format(i6,i2,3f7.4,f9.5,f10.5,2f7.4,4i4)
C            write(6,256) xkmin1(i),xkmax1(i),xkmin2(i),xkmax2(i),
C     *                   xkmin3(i),xkmax3(i),dmean(xk11,nagnt1),
C     *                   dmean(xk21,nagnt2),dmean(xk31,nagnt3)
C256         format(9f8.3)
            write(6,256) top(1),top(2),top(3),top(4),top(5),
     *                   top(6),top(7),top(8),top(9),gini,
     *                   xmax3(xkmax1(i),xkmax2(i),xkmax3(i))
256         format(8f6.4,f7.4,f6.4,f15.5)
         endif

         if ( (nz(i) .eq. 1) .and. (nz(i+1) .eq. 1) ) then
            preps(2,2) = pgg00
            preps(2,1) = pgg01
         elseif ( (nz(i) .eq. 1) .and. (nz(i+1) .eq. 0) ) then
            preps(2,2) = pbg00
            preps(2,1) = pbg01
         elseif ( (nz(i) .eq. 0) .and. (nz(i+1) .eq. 1) ) then
            preps(2,2) = pgb00
            preps(2,1) = pgb01
         elseif ( (nz(i) .eq. 0) .and. (nz(i+1) .eq. 0) ) then
            preps(2,2) = pbb00
            preps(2,1) = pbb01
         endif

         call ggubs(sseed,nagnt1,uniep1)
         do 345 j = 1,nagnt1
            if (neps11(j) .eq. 1) then
               if (uniep1(j) .le. preps(2,1)) then
                  neps12(j) = 0
               else
                  neps12(j) = 1
               endif
            else
               if (uniep1(j) .le. preps(2,2)) then
                  neps12(j) = 0
               else
                  neps12(j) = 1
               endif
            endif
345      continue

         call ggubs(sseed,nagnt2,uniep2)
         do 346 j = 1,nagnt2
            if (neps21(j) .eq. 1) then
               if (uniep2(j) .le. preps(2,1)) then
                  neps22(j) = 0
               else
                  neps22(j) = 1
               endif
            else
               if (uniep2(j) .le. preps(2,2)) then
                  neps22(j) = 0
               else
                  neps22(j) = 1
               endif
            endif
346      continue

         call ggubs(sseed,nagnt3,uniep3)
         do 347 j = 1,nagnt3
            if (neps31(j) .eq. 1) then
               if (uniep3(j) .le. preps(2,1)) then
                  neps32(j) = 0
               else
                  neps32(j) = 1
               endif
            else
               if (uniep3(j) .le. preps(2,2)) then
                  neps32(j) = 0
               else
                  neps32(j) = 1
               endif
            endif
347      continue

         call domean(xk12,xk22,xk32,nagnt1,nagnt2,nagnt3,wght1,
     *               wght2,wght3,xksum)

         mswch1 = idnint(swtch1)
         mswch3 = idnint(swtch3)
         mswc2a = mswch1
         mswc2b = mswch3

         call nrand2(sseed,1,nagnt1,nums1,mswch1)
         call nrand2(sseed,1,nagnt2,nums2,mswch1+mswch3)
         call nrand2(sseed,1,nagnt3,nums3,mswch3)

         do 564 j = 1,mswch1
            ntmp1(j) = neps12(nums1(j))
            xktmp1(j) = xk12(nums1(j))
            neps12(nums1(j)) = neps22(nums2(j))
            xk12(nums1(j)) = xk22(nums2(j))
564      continue

         do 565 j = 1,mswch3
            ntmp3(j) = neps32(nums3(j))
            xktmp3(j) = xk32(nums3(j))
            neps32(nums3(j)) = neps22(nums2(j+mswch1))
            xk32(nums3(j)) = xk22(nums2(j+mswch1))
 565     continue

         if (mswc2a .gt. 0) then
            do 566 j = 1,mswc2a
               neps22(nums2(j)) = ntmp1(j)
               xk22(nums2(j)) = xktmp1(j)
566         continue
         endif

         if (mswc2b .gt. 0) then
            do 567 j = 1,mswc2b
               neps22(nums2(j+mswch1)) = ntmp3(j)
               xk22(nums2(j+mswch1)) = xktmp3(j)
567         continue
         endif

         call lawlrg(sseed,neps12,nz(i+1),nagnt1,nswch1)
         call lawlrg(sseed,neps22,nz(i+1),nagnt2,nswch2)
         call lawlrg(sseed,neps32,nz(i+1),nagnt3,nswch3)

         do 587 j = 1,nagnt1
            neps11(j) = neps12(j)
            xk11(j) = xk12(j)
 587     continue

         do 588 j = 1,nagnt2
            neps21(j) = neps22(j)
            xk21(j) = xk22(j)
 588     continue

         do 589 j = 1,nagnt3
            neps31(j) = neps32(j)
            xk31(j) = xk32(j)
 589     continue

         call dopct(xk11,xk21,xk31,neps11,neps21,neps31,nagnt1,nagnt2,
     *              nagnt3,wght1,wght2,wght3,pctemp,p0,pcon)

 10   continue

      open(1,file='cross.out',form='formatted')

      ntype = 1
      do 779 i = 1,nagnt1
         write(1,855) i,ntype,neps11(i),xk11(i)
779   continue

      ntype = 2
      do 773 i = 1,nagnt2
         write(1,855) i+nagnt1,ntype,neps21(i),xk21(i)
773   continue

      ntype = 3
      do 774 i = 1,nagnt3
         write(1,855) i+nagnt1+nagnt2,ntype,neps31(i),xk31(i)
774   continue

      close(1)

      do 350 j = 1,kgrid-1
         hist(j) = 0.0D+00
         fractn(j) = 0.0D+00
350   continue

      xkbar(ncapt+1) = xksum
      xkistd(ncapt+1) = dstd3(xk11,xk21,xk31,nagnt1,nagnt2,nagnt3,
     *                        xkbar(ncapt+1),prob(1),prob(2),
     *                        prob(3))
      pct0(ncapt+1) = p0
      pctcon(ncapt+1) = pcon

      do 870 j = 1,nagnt1
         call locidx(xk11(j),xkgpts,kgrid,kcor,klow,khigh,niter)
         if (klow .eq. 1) then
            kcor = 1
         elseif (khigh .eq. 1) then
            kcor = kgrid-1
         endif
         hist(kcor) = hist(kcor) + wght1/xnagnt
         fractn(kcor) = fractn(kcor) + wght1*xk11(j)/xnagnt
870   continue

      do 871 j = 1,nagnt2
         call locidx(xk21(j),xkgpts,kgrid,kcor,klow,khigh,niter)
         if (klow .eq. 1) then
            kcor = 1
         elseif (khigh .eq. 1) then
            kcor = kgrid-1
         endif
         hist(kcor) = hist(kcor) + wght2/xnagnt
         fractn(kcor) = fractn(kcor) + wght2*xk21(j)/xnagnt
871   continue

      do 872 j = 1,nagnt3
         call locidx(xk31(j),xkgpts,kgrid,kcor,klow,khigh,niter)
         if (klow .eq. 1) then
            kcor = 1
         elseif (khigh .eq. 1) then
            kcor = kgrid-1
         endif
         hist(kcor) = hist(kcor) + wght3/xnagnt
         fractn(kcor) = fractn(kcor) + wght3*xk31(j)/xnagnt
872   continue

      open(10,file='heter40w.hpt',form='formatted')
      sumx = 0.0D+00
      sumy = 0.0D+00
      xpt(1) = 0.0D+00
      ypt(1) = 0.0D+00
      write(10,736) xpt(1),ypt(1)
736   format(2f15.8)
      nhpts = 1
      do 852 j = 2,kgrid
         sumx = sumx + hist(j-1)
         sumy = sumy + fractn(j-1)/xkbar(ncapt+1)
         if (hist(j-1) .gt. 1.0D-10) then
            nhpts = nhpts + 1
            xpt(nhpts) = sumx
            ypt(nhpts) = sumy
            write(10,736) xpt(nhpts),ypt(nhpts)
         endif
852   continue
      close(10)

      call calcsd(xpt,ypt,ydp,a,b,c,r,nhpts,-one,-one)
      call integr(xpt,ypt,ydp,nhpts,xint)
      gini = one - 2.0D+00*xint

      do 897 i = 1,ntop
         call doint(xpt,ypt,ydp,one-toppct(i),yval,yp,ydp,
     *              khi,klo,0,nhpts,1,0,0)
         top(i) = one - yval
897   continue

      ginis(ncapt+1) = gini
      do 877 j = 1,ntop
         topdat(ncapt+1,j) = top(j)
877   continue

      do 597 i = 1,nlzpts
         call doint(xpt,ypt,ydp,xlzpts(i),yval,yp,ydp,
     *              khi,klo,0,nhpts,1,0,0)
         ylzpts(i) = yval
597   continue

      open(1,file='heter40w.lrz',form='formatted')
      do 56 i = 1,nlzpts
         if ( (i .eq. 141) .or. (i .eq. 161) .or.
     *        (i .eq. 181) .or. (i .eq. 191) .or.
     *        (i .eq. 199) .or.
     *        (i .eq. 1) .or. (i .eq. 201) ) then
            write(1,253) xlzpts(i),ylzpts(i),xlzdat(i)
253         format(3f12.8)
         else
            write(1,259) xlzpts(i),ylzpts(i)
259         format(2f12.8,'   .')
         endif
 56   continue
      close(1)

      do 96 i = 1,ncapt
         if (nz(i) .eq. 1) then
            hbar(i) = hfix*(one-unempg)
            ybar(i) = zgood*(xkbar(i)**alpha)*
     *                (hbar(i)**(one-alpha))
         else
            hbar(i) = hfix*(one-unempb)
            ybar(i) = zbad*(xkbar(i)**alpha)*
     *                (hbar(i)**(one-alpha))
         endif
         xbar(i) = xkbar(i+1) - (one-delta)*xkbar(i)
         cbar(i) = ybar(i) - xbar(i)
 96   continue

      open(1,file='heter40w.sim',form='formatted')
      do 68 i = ndrop+1,ncapt
         write(1,367) i-ndrop,nz(i),xkbar(i),xkbar(i+1),xkistd(i),
     *                pct0(i),xkmax1(i),xkmax2(i),xkmax3(i),
     *                xkmin1(i),xkmin2(i),xkmin3(i),topdat(i,1),
     *                topdat(i,2),topdat(i,3),topdat(i,4),
     *                topdat(i,5),topdat(i,6),topdat(i,7),
     *                topdat(i,8),topdat(i,9),ginis(i)
367      format(i6,i3,20f12.6)
 68   continue
      close(1)

      do 636 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = pct0(i)
636   continue
      pct0mn = dmean(tmpdat,ncapt-ndrop)
      pct0sd = dstd(tmpdat,ncapt-ndrop)
      pct0mx = dmax(tmpdat,ncapt-ndrop)
      pct0mi = dmin(tmpdat,ncapt-ndrop)

      do 637 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = pctcon(i)
637   continue
      pctcmn = dmean(tmpdat,ncapt-ndrop)
      pctcsd = dstd(tmpdat,ncapt-ndrop)
      pctcmx = dmax(tmpdat,ncapt-ndrop)
      pctcmi = dmin(tmpdat,ncapt-ndrop)

      do 236 j = 1,ntop
         do 237 i = ndrop+1,ncapt
            tmpdat(i-ndrop) = topdat(i,j)
237      continue
         topmn(j) = dmean(tmpdat,ncapt-ndrop)
         topsd(j) = dstd(tmpdat,ncapt-ndrop)
         topmx(j) = dmax(tmpdat,ncapt-ndrop)
         topmi(j) = dmin(tmpdat,ncapt-ndrop)
236   continue

      do 830 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = xkistd(i)
830   continue
      xksdmn = dmean(tmpdat,ncapt-ndrop)
      xksdsd = dstd(tmpdat,ncapt-ndrop)
      xksdmx = dmax(tmpdat,ncapt-ndrop)
      xksdmi = dmin(tmpdat,ncapt-ndrop)

C      do 630 i = ndrop+1,ncapt
C         tmpdat(i-ndrop) = xkiskw(i)
C630   continue
C      xkskmn = dmean(tmpdat,ncapt-ndrop)
C      xksksd = dstd(tmpdat,ncapt-ndrop)
C      xkskmx = dmax(tmpdat,ncapt-ndrop)
C      xkskmi = dmin(tmpdat,ncapt-ndrop)

C      do 730 i = ndrop+1,ncapt
C         tmpdat(i-ndrop) = xkikrt(i)
C730   continue
C      xkkrmn = dmean(tmpdat,ncapt-ndrop)
C      xkkrsd = dstd(tmpdat,ncapt-ndrop)
C      xkkrmx = dmax(tmpdat,ncapt-ndrop)
C      xkkrmi = dmin(tmpdat,ncapt-ndrop)

      do 138 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = ginis(i)
138   continue
      xgmn = dmean(tmpdat,ncapt-ndrop)
      xgsd = dstd(tmpdat,ncapt-ndrop)
      xgmx = dmax(tmpdat,ncapt-ndrop)
      xgmi = dmin(tmpdat,ncapt-ndrop)

      do 30 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = cbar(i)
 30   continue
      cmean = dmean(tmpdat,ncapt-ndrop)
      cstd = dstd(tmpdat,ncapt-ndrop)
      cmax = dmax(tmpdat,ncapt-ndrop)
      cmin = dmin(tmpdat,ncapt-ndrop)

      do 21 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = xbar(i)
 21   continue
      xmean = dmean(tmpdat,ncapt-ndrop)
      xstd = dstd(tmpdat,ncapt-ndrop)
      xmax = dmax(tmpdat,ncapt-ndrop)
      xmin = dmin(tmpdat,ncapt-ndrop)

      do 22 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = ybar(i)
 22   continue
      ymean = dmean(tmpdat,ncapt-ndrop)
      ystd = dstd(tmpdat,ncapt-ndrop)
      ymax = dmax(tmpdat,ncapt-ndrop)
      ymin = dmin(tmpdat,ncapt-ndrop)

      do 33 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = xkmax1(i)
 33   continue
      x1mxmn = dmean(tmpdat,ncapt-ndrop)
      x1mxsd = dstd(tmpdat,ncapt-ndrop)
      x1mxmx = dmax(tmpdat,ncapt-ndrop)
      x1mxmi = dmin(tmpdat,ncapt-ndrop)

      do 133 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = xkmax2(i)
133   continue
      x2mxmn = dmean(tmpdat,ncapt-ndrop)
      x2mxsd = dstd(tmpdat,ncapt-ndrop)
      x2mxmx = dmax(tmpdat,ncapt-ndrop)
      x2mxmi = dmin(tmpdat,ncapt-ndrop)

      do 203 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = xkmax3(i)
203   continue
      x3mxmn = dmean(tmpdat,ncapt-ndrop)
      x3mxsd = dstd(tmpdat,ncapt-ndrop)
      x3mxmx = dmax(tmpdat,ncapt-ndrop)
      x3mxmi = dmin(tmpdat,ncapt-ndrop)

      do 393 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = xkmin1(i)
393   continue
      x1mimn = dmean(tmpdat,ncapt-ndrop)
      x1misd = dstd(tmpdat,ncapt-ndrop)
      x1mimx = dmax(tmpdat,ncapt-ndrop)
      x1mimi = dmin(tmpdat,ncapt-ndrop)

      do 433 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = xkmin2(i)
433   continue
      x2mimn = dmean(tmpdat,ncapt-ndrop)
      x2misd = dstd(tmpdat,ncapt-ndrop)
      x2mimx = dmax(tmpdat,ncapt-ndrop)
      x2mimi = dmin(tmpdat,ncapt-ndrop)

      do 533 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = xkmin3(i)
533   continue
      x3mimn = dmean(tmpdat,ncapt-ndrop)
      x3misd = dstd(tmpdat,ncapt-ndrop)
      x3mimx = dmax(tmpdat,ncapt-ndrop)
      x3mimi = dmin(tmpdat,ncapt-ndrop)

      do 633 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = xkbar(i)
633   continue
      xkmean = dmean(tmpdat,ncapt-ndrop)
      xkstd = dstd(tmpdat,ncapt-ndrop)
      xkmax = dmax(tmpdat,ncapt-ndrop)
      xkmin = dmin(tmpdat,ncapt-ndrop)

      do 46 i = ndrop+1,ncapt
         tmpdat(i-ndrop) = hbar(i)
 46   continue
      hmean = dmean(tmpdat,ncapt-ndrop)
      hstd = dstd(tmpdat,ncapt-ndrop)
      hmax = dmax(tmpdat,ncapt-ndrop)
      hmin = dmin(tmpdat,ncapt-ndrop)

      npts = 0
      do 34 i = ndrop+1,ncapt
         if (nz(i) .eq. 1) then
            npts = npts + 1
            xdata(npts,1) = one
            xdata(npts,2) = dlog(xkbar(i))
            ydata(npts) = dlog(xkbar(i+1))
         endif
 34   continue
      call ols2(ydata,xdata,ncapt,bhat1,rsqa,sderra,npts)

      npts = 0
      do 37 i = ndrop+1,ncapt
         if (nz(i) .eq. 0) then
            npts = npts + 1
            xdata(npts,1) = one
            xdata(npts,2) = dlog(xkbar(i))
            ydata(npts) = dlog(xkbar(i+1))
         endif
 37   continue
      call ols2(ydata,xdata,ncapt,bhat2,rsqb,sderrb,npts)

      return
      end

      subroutine splint(xa,ya,y2a,x,y,yp,ydp,nycor,nzcor,nbcor,
     *                  khi,klo,kvals,ndoy,ndoyp,ndoydp)
      implicit real*8 (a-h,o-z)
      parameter (nkpts=132,nmupts=4,n=nkpts,m=nmupts,nbetas=3)
      dimension xa(n),ya(n,m,nbetas,2),y2a(n,m,nbetas,2)

      one = 1.0D+00

      if (kvals .eq. 0) then
         klo = 1
         khi = n
   1     if ((khi-klo) .gt. 1) then
            k = (khi+klo)/2
            if (xa(k) .gt. x) then
               khi = k
            else
               klo = k
            endif
            goto 1
         endif
      endif

      h = xa(khi) - xa(klo)
      a = (xa(khi) - x)/h
      b = (x - xa(klo))/h
      asq = a*a
      bsq = b*b

      if (ndoy .eq. 1)
     *   y = a*ya(klo,nycor,nbcor,nzcor) + b*ya(khi,nycor,nbcor,nzcor) +
     *       ((asq*a-a)*y2a(klo,nycor,nbcor,nzcor)
     *        +(bsq*b-b)*y2a(khi,nycor,nbcor,nzcor))*(h*h)/6.0D+00

      if (ndoyp .eq. 1)
     *   yp = (ya(khi,nycor,nbcor,nzcor)-ya(klo,nycor,nbcor,nzcor))/h -
     *        (3.0D+00*asq-one)/6.0D+00*h*y2a(klo,nycor,nbcor,nzcor) +
     *        (3.0D+00*bsq-one)/6.0D+00*h*y2a(khi,nycor,nbcor,nzcor)

      if (ndoydp .eq. 1)
     *   ydp = a*y2a(klo,nycor,nbcor,nzcor) +
     *         b*y2a(khi,nycor,nbcor,nzcor)

      return
      end

      subroutine tridag(a,b,c,r,u,n)
      implicit real*8 (a-h,o-z)
      parameter (nmax=350,toler=1.0D-12)
      dimension gam(nmax),a(n),b(n),c(n),r(n),u(n)

      bet = b(1)
      u(1) = r(1)/bet
      do 11 j = 2,n
         gam(j) = c(j-1)/bet
         bet = b(j) - a(j)*gam(j)
         if (dabs(bet) .le. toler) then
            write(6,100)
100         format(' Failure in subroutine tridag')
            pause
         endif
         u(j) = (r(j)-a(j)*u(j-1))/bet
 11   continue
      do 12 j = n-1,1,-1
         u(j) = u(j) - gam(j+1)*u(j+1)
 12   continue

      return
      end

      subroutine polint(xa,ya,n,x,y,dy)
      implicit real*8 (a-h,o-z)
      parameter (nmax=10,toler=1.0D-12)
      dimension xa(n),ya(n),c(nmax),d(nmax)

      ns = 1
      dif = dabs(x-xa(1))
      do 11 i = 1,n
         dift = dabs(x-xa(i))
         if (dift .lt. dif) then
            ns = i
            dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
 11   continue

      y = ya(ns)
      ns = ns - 1
      do 13 m = 1,n-1
         do 12 i = 1,n-m
            ho = xa(i) - x
            hp = xa(i+m) - x
            w = c(i+1) - d(i)
            den = ho - hp
            if (dabs(den) .lt. toler) then
               write(6,100)
100            format(' Error in subroutine polint')
               pause
            endif
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
 12      continue
         if ((2*ns) .lt. (n-m)) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns - 1
         endif
         y = y + dy
 13   continue

      return
      end

      subroutine interp(x1,x2,y,a0,a1,c0,c1)
      implicit real*8 (a-h,o-z)
      parameter (nkpts=132,nmupts=4,n=nkpts,m=nmupts,nbetas=3)
      common/cvs/v1(nkpts,nmupts,nbetas,2),v2(nkpts,nmupts,nbetas,2),
     *           v3(nkpts,nmupts,nbetas,2),v4(nkpts,nmupts,nbetas,2),
     *           v1dp(nkpts,nmupts,nbetas,2),
     *           v2dp(nkpts,nmupts,nbetas,2),
     *           v3dp(nkpts,nmupts,nbetas,2),
     *           v4dp(nkpts,nmupts,nbetas,2)
      dimension x1(n),x2(m),y(n,m,4,nbetas),a(n),b(n),c(n),r(n),u(n),
     *          xx(n),yy(n),xa(m),ya(m),xmup(m,2)

      do 80 i = 1,m
         xmup(i,1) = dexp(a0 + a1*dlog(x2(i)))
         xmup(i,2) = dexp(c0 + c1*dlog(x2(i)))
 80   continue

      do 70 i = 1,m
         xa(i) = x2(i)
 70   continue

      do 72 i = 1,n
         xx(i) = x1(i)
 72   continue

      a(1) = 0.0D+00
      a(n) = -1.0D+00
      b(1) = 1.0D+00
      b(n) = 1.05D+00
      c(1) = -1.05D+00
      c(n) = 0.0D+00
      r(1) = 0.0D+00
      r(n) = 0.0D+00

      do 81 k = 1,4
         do 60 kk = 1,nbetas
            do 82 ii = 1,m
               do 83 jj = 1,2

                  do 30 i = 1,n
                     do 71 j = 1,m
                        ya(j) = y(i,j,k,kk)
 71                  continue
                     call polint(xa,ya,m,xmup(ii,jj),yval,dy)
                     yy(i) = yval
                     if (k .eq. 1) then
                        v1(i,ii,kk,jj) = yval
                     elseif (k .eq. 2) then
                        v2(i,ii,kk,jj) = yval
                     elseif (k .eq. 3) then
                        v3(i,ii,kk,jj) = yval
                     elseif (k .eq. 4) then
                        v4(i,ii,kk,jj) = yval
                     endif
 30               continue

                  do 21 i = 2,n-1
                     a(i) = (xx(i)-xx(i-1))/6.0D+00
                     b(i) = (xx(i+1)-xx(i-1))/3.0D+00
                     c(i) = (xx(i+1)-xx(i))/6.0+00
                     r(i) = (yy(i+1)-yy(i))/(xx(i+1)-xx(i)) -
     *                      (yy(i)-yy(i-1))/(xx(i)-xx(i-1))
 21               continue

                  call tridag(a,b,c,r,u,n)
                  do 22 i = 1,n
                     if (k .eq. 1) then
                        v1dp(i,ii,kk,jj) = u(i)
                     elseif (k .eq. 2) then
                        v2dp(i,ii,kk,jj) = u(i)
                     elseif (k .eq. 3) then
                        v3dp(i,ii,kk,jj) = u(i)
                     elseif (k .eq. 4) then
                        v4dp(i,ii,kk,jj) = u(i)
                     endif
 22               continue

 83            continue
 82         continue
 60      continue
 81   continue

      return
      end

      subroutine calop2(xk,neps,ngood,ynew,value,niter,ndone,toler,
     *                  grad,xkpts,r0cur,w0cur,ykhigh,yklow,nbeta)
      implicit real*8 (a-h,o-z)
      parameter (mxiter=100,hfix=0.3271D+00,
     *           nkpts=132,delta=0.025D+00,nbetas=3,whome=0.07D+00)
      common/cpr/pr(4,4),prbeta(nbetas,nbetas),prob(nbetas)
      common/cbetas/betas(nbetas)
      common/cws/w1(nkpts,nbetas,2),w2(nkpts,nbetas,2),
     *           w3(nkpts,nbetas,2),w4(nkpts,nbetas,2),
     *           w1dp(nkpts,nbetas,2),w2dp(nkpts,nbetas,2),
     *           w3dp(nkpts,nbetas,2),w4dp(nkpts,nbetas,2)
      dimension xkpts(nkpts),r0cur(2),w0cur(2),
     *          y1s(nbetas),y2s(nbetas),y3s(nbetas),y4s(nbetas),
     *          y1ps(nbetas),y2ps(nbetas),y3ps(nbetas),y4ps(nbetas)

      one = 1.0D+00

      beta = betas(nbeta)

      ndone = 0
      niter = 0

      if (neps .eq. 1) then
         y = r0cur(ngood)*xk + w0cur(ngood)*hfix
      else
         y = r0cur(ngood)*xk + whome
      endif

      xklow = yklow
      xkhigh = ykhigh

      do 10 i = 1,mxiter

         xkcur = (xkhigh+xklow)/2.0D+00

         x = xkcur - (one-delta)*xk
         if (x .ge. y) then
            xkhigh = xkcur
            goto 10
         endif

         c = y - x

         do 87 j = 1,nbetas
            call splin2(xkpts,w1,w1dp,xkcur,y1,y1p,y1dp,ngood,j,
     *                  khi,klo,0,0,1,0)
            y1ps(j) = y1p
            call splin2(xkpts,w2,w2dp,xkcur,y2,y2p,y2dp,ngood,j,
     *                  khi,klo,1,0,1,0)
            y2ps(j) = y2p
            call splin2(xkpts,w3,w3dp,xkcur,y3,y3p,y3dp,ngood,j,
     *                  khi,klo,1,0,1,0)
            y3ps(j) = y3p
            call splin2(xkpts,w4,w4dp,xkcur,y4,y4p,y4dp,ngood,j,
     *                  khi,klo,1,0,1,0)
            y4ps(j) = y4p
 87      continue

         y1p = 0.0D+00
         y2p = 0.0D+00
         y3p = 0.0D+00
         y4p = 0.0D+00
         do 81 j = 1,nbetas
            y1p = y1p + prbeta(j,nbeta)*y1ps(j)
            y2p = y2p + prbeta(j,nbeta)*y2ps(j)
            y3p = y3p + prbeta(j,nbeta)*y3ps(j)
            y4p = y4p + prbeta(j,nbeta)*y4ps(j)
 81      continue

         upc = -one/c
         if ( (neps .eq. 1) .and. (ngood .eq. 1) ) then
            grad = upc + beta*(y1p*pr(1,1) + y2p*pr(2,1) +
     *                         y3p*pr(3,1) + y4p*pr(4,1))
         elseif ( (neps .eq. 1) .and. (ngood .eq. 2) ) then
            grad = upc + beta*(y1p*pr(1,2) + y2p*pr(2,2) +
     *                         y3p*pr(3,2) + y4p*pr(4,2))
         elseif ( (neps .eq. 0) .and. (ngood .eq. 1) ) then
            grad = upc + beta*(y1p*pr(1,3) + y2p*pr(2,3) +
     *                         y3p*pr(3,3) + y4p*pr(4,3))
         else
            grad = upc + beta*(y1p*pr(1,4) + y2p*pr(2,4) +
     *                         y3p*pr(3,4) + y4p*pr(4,4))
         endif

         if (grad .ge. 0.0D+00) then
            xklow = xkcur
         else
            xkhigh = xkcur
         endif

         if ((xkhigh-xklow) .le. toler) then
            niter = i
            if (dabs(grad) .le. 1.0D-02) ndone = 1
            goto 999
         endif

 10   continue

999   xkcur = (xklow + xkhigh)/2.0D+00
      x = xkcur - (one-delta)*xk
      c = y - x
      if (c .le. 0.0D+00) then
         write(6,777) c
777      format(' Failure in calop2: negative c',f15.8)
         pause
      endif
      utilc = dlog(c)
      do 78 j = 1,nbetas
         call splin2(xkpts,w1,w1dp,xkcur,y1,y1p,y1dp,ngood,j,
     *               khi,klo,0,1,0,0)
         y1s(j) = y1
         call splin2(xkpts,w2,w2dp,xkcur,y2,y2p,y2dp,ngood,j,
     *               khi,klo,1,1,0,0)
         y2s(j) = y2
         call splin2(xkpts,w3,w3dp,xkcur,y3,y3p,y3dp,ngood,j,
     *               khi,klo,1,1,0,0)
         y3s(j) = y3
         call splin2(xkpts,w4,w4dp,xkcur,y4,y4p,y4dp,ngood,j,
     *               khi,klo,1,1,0,0)
         y4s(j) = y4
 78   continue
      y1 = 0.0D+00
      y2 = 0.0D+00
      y3 = 0.0D+00
      y4 = 0.0D+00
      do 79 j = 1,nbetas
         y1 = y1 + prbeta(j,nbeta)*y1s(j)
         y2 = y2 + prbeta(j,nbeta)*y2s(j)
         y3 = y3 + prbeta(j,nbeta)*y3s(j)
         y4 = y4 + prbeta(j,nbeta)*y4s(j)
 79   continue
      if ( (neps .eq. 1) .and. (ngood .eq. 1) ) then
         fval = utilc + beta*(y1*pr(1,1) + y2*pr(2,1) +
     *                        y3*pr(3,1) + y4*pr(4,1))
      elseif ( (neps .eq. 1) .and. (ngood .eq. 2) ) then
         fval = utilc + beta*(y1*pr(1,2) + y2*pr(2,2) +
     *                        y3*pr(3,2) + y4*pr(4,2))
      elseif ( (neps .eq. 0) .and. (ngood .eq. 1) ) then
         fval = utilc + beta*(y1*pr(1,3) + y2*pr(2,3) +
     *                        y3*pr(3,3) + y4*pr(4,3))
      else
         fval = utilc + beta*(y1*pr(1,4) + y2*pr(2,4) +
     *                        y3*pr(3,4) + y4*pr(4,4))
      endif

      ynew = xkcur
      value = fval

      return
      end

      subroutine altop2(xk,neps,ngood,ynew,value,niter,ndone,toler,
     *                  startk,grad,xkpts,r0cur,w0cur,nbeta)
      implicit real*8 (a-h,o-z)
      parameter (mxiter=100,hfix=0.3271D+00,
     *           nkpts=132,delta=0.025D+00,nbetas=3,whome=0.07D+00)
      common/cpr/pr(4,4),prbeta(nbetas,nbetas),prob(nbetas)
      common/cbetas/betas(nbetas)
      common/cws/w1(nkpts,nbetas,2),w2(nkpts,nbetas,2),
     *           w3(nkpts,nbetas,2),w4(nkpts,nbetas,2),
     *           w1dp(nkpts,nbetas,2),w2dp(nkpts,nbetas,2),
     *           w3dp(nkpts,nbetas,2),w4dp(nkpts,nbetas,2)
      dimension xkpts(nkpts),r0cur(2),w0cur(2),
     *          y1s(nbetas),y2s(nbetas),y3s(nbetas),y4s(nbetas),
     *          y1ps(nbetas),y2ps(nbetas),y3ps(nbetas),y4ps(nbetas),
     *          y1dps(nbetas),y2dps(nbetas),y3dps(nbetas),y4dps(nbetas)

      one = 1.0D+00

      beta = betas(nbeta)

      ndone = 0
      niter = 0

      if (neps .eq. 1) then
         y = r0cur(ngood)*xk + w0cur(ngood)*hfix
      else
         y = r0cur(ngood)*xk + whome
      endif

      xk2p = startk

      do 10 i = 1,mxiter

         xk2 = xk2p
         x = xk2 - (one-delta)*xk
         c = y - x

         do 87 j = 1,nbetas
            call splin2(xkpts,w1,w1dp,xk2,y1,y1p,y1dp,ngood,j,
     *                  khi,klo,0,0,1,1)
            y1ps(j) = y1p
            y1dps(j) = y1dp
            call splin2(xkpts,w2,w2dp,xk2,y2,y2p,y2dp,ngood,j,
     *                  khi,klo,1,0,1,1)
            y2ps(j) = y2p
            y2dps(j) = y2dp
            call splin2(xkpts,w3,w3dp,xk2,y3,y3p,y3dp,ngood,j,
     *                  khi,klo,1,0,1,1)
            y3ps(j) = y3p
            y3dps(j) = y3dp
            call splin2(xkpts,w4,w4dp,xk2,y4,y4p,y4dp,ngood,j,
     *                  khi,klo,1,0,1,1)
            y4ps(j) = y4p
            y4dps(j) = y4dp
 87      continue

         y1p = 0.0D+00
         y2p = 0.0D+00
         y3p = 0.0D+00
         y4p = 0.0D+00
         y1dp = 0.0D+00
         y2dp = 0.0D+00
         y3dp = 0.0D+00
         y4dp = 0.0D+00
         do 81 j = 1,nbetas
            y1p = y1p + prbeta(j,nbeta)*y1ps(j)
            y2p = y2p + prbeta(j,nbeta)*y2ps(j)
            y3p = y3p + prbeta(j,nbeta)*y3ps(j)
            y4p = y4p + prbeta(j,nbeta)*y4ps(j)
            y1dp = y1dp + prbeta(j,nbeta)*y1dps(j)
            y2dp = y2dp + prbeta(j,nbeta)*y2dps(j)
            y3dp = y3dp + prbeta(j,nbeta)*y3dps(j)
            y4dp = y4dp + prbeta(j,nbeta)*y4dps(j)
 81      continue

         upc = -one/c
         if ( (neps .eq. 1) .and. (ngood .eq. 1) ) then
            grad = upc + beta*(y1p*pr(1,1) + y2p*pr(2,1) +
     *                         y3p*pr(3,1) + y4p*pr(4,1))
         elseif ( (neps .eq. 1) .and. (ngood .eq. 2) ) then
            grad = upc + beta*(y1p*pr(1,2) + y2p*pr(2,2) +
     *                         y3p*pr(3,2) + y4p*pr(4,2))
         elseif ( (neps .eq. 0) .and. (ngood .eq. 1) ) then
            grad = upc + beta*(y1p*pr(1,3) + y2p*pr(2,3) +
     *                         y3p*pr(3,3) + y4p*pr(4,3))
         else
            grad = upc + beta*(y1p*pr(1,4) + y2p*pr(2,4) +
     *                         y3p*pr(3,4) + y4p*pr(4,4))
         endif

         udpc = -one/(c*c)
         if ( (neps .eq. 1) .and. (ngood .eq. 1) ) then
            hess = udpc + beta*(y1dp*pr(1,1) + y2dp*pr(2,1) +
     *                          y3dp*pr(3,1) + y4dp*pr(4,1))
         elseif ( (neps .eq. 1) .and. (ngood .eq. 2) ) then
            hess = udpc + beta*(y1dp*pr(1,2) + y2dp*pr(2,2) +
     *                          y3dp*pr(3,2) + y4dp*pr(4,2))
         elseif ( (neps .eq. 0) .and. (ngood .eq. 1) ) then
            hess = udpc + beta*(y1dp*pr(1,3) + y2dp*pr(2,3) +
     *                          y3dp*pr(3,3) + y4dp*pr(4,3))
         else
            hess = udpc + beta*(y1dp*pr(1,4) + y2dp*pr(2,4) +
     *                          y3dp*pr(3,4) + y4dp*pr(4,4))
         endif

         xk2p = xk2 - grad/hess

         if (dabs(xk2p-xk2) .le. toler) then
            ndone = 1
            niter = i
            goto 999
         endif

 10   continue

999   if (c .gt. 0.0D+00) then

         utilc = dlog(c)
         do 78 j = 1,nbetas
            call splin2(xkpts,w1,w1dp,xk2,y1,y1p,y1dp,ngood,j,
     *                  khi,klo,0,1,0,0)
            y1s(j) = y1
            call splin2(xkpts,w2,w2dp,xk2,y2,y2p,y2dp,ngood,j,
     *                  khi,klo,1,1,0,0)
            y2s(j) = y2
            call splin2(xkpts,w3,w3dp,xk2,y3,y3p,y3dp,ngood,j,
     *                  khi,klo,1,1,0,0)
            y3s(j) = y3
            call splin2(xkpts,w4,w4dp,xk2,y4,y4p,y4dp,ngood,j,
     *                  khi,klo,1,1,0,0)
            y4s(j) = y4
 78      continue
         y1 = 0.0D+00
         y2 = 0.0D+00
         y3 = 0.0D+00
         y4 = 0.0D+00
         do 79 j = 1,nbetas
            y1 = y1 + prbeta(j,nbeta)*y1s(j)
            y2 = y2 + prbeta(j,nbeta)*y2s(j)
            y3 = y3 + prbeta(j,nbeta)*y3s(j)
            y4 = y4 + prbeta(j,nbeta)*y4s(j)
 79      continue
         if ( (neps .eq. 1) .and. (ngood .eq. 1) ) then
            fval = utilc + beta*(y1*pr(1,1) + y2*pr(2,1) +
     *                           y3*pr(3,1) + y4*pr(4,1))
         elseif ( (neps .eq. 1) .and. (ngood .eq. 2) ) then
            fval = utilc + beta*(y1*pr(1,2) + y2*pr(2,2) +
     *                           y3*pr(3,2) + y4*pr(4,2))
         elseif ( (neps .eq. 0) .and. (ngood .eq. 1) ) then
            fval = utilc + beta*(y1*pr(1,3) + y2*pr(2,3) +
     *                           y3*pr(3,3) + y4*pr(4,3))
         else
            fval = utilc + beta*(y1*pr(1,4) + y2*pr(2,4) +
     *                           y3*pr(3,4) + y4*pr(4,4))
         endif

         ynew = xk2
         value = fval

      else
         ndone = 0
      endif

      return
      end

      subroutine inter2(x1,x2,y,a0,a1,c0,c1,xmu)
      implicit real*8 (a-h,o-z)
      parameter (nkpts=132,nmupts=4,n=nkpts,m=nmupts,nbetas=3)
      common/cws/w1(nkpts,nbetas,2),w2(nkpts,nbetas,2),
     *           w3(nkpts,nbetas,2),w4(nkpts,nbetas,2),
     *           w1dp(nkpts,nbetas,2),w2dp(nkpts,nbetas,2),
     *           w3dp(nkpts,nbetas,2),w4dp(nkpts,nbetas,2)
      dimension x1(n),x2(m),y(n,m,4,nbetas),a(n),b(n),c(n),r(n),u(n),
     *          xx(n),yy(n),xa(m),ya(m),xmup(2)

      xmup(1) = dexp(a0 + a1*dlog(xmu))
      xmup(2) = dexp(c0 + c1*dlog(xmu))

      do 70 i = 1,m
         xa(i) = x2(i)
 70   continue

      do 72 i = 1,n
         xx(i) = x1(i)
 72   continue

      a(1) = 0.0D+00
      a(n) = -1.0D+00
      b(1) = 1.0D+00
      b(n) = 1.05D+00
      c(1) = -1.05D+00
      c(n) = 0.0D+00
      r(1) = 0.0D+00
      r(n) = 0.0D+00

      do 81 k = 1,4
         do 60 kk = 1,nbetas
            do 82 ii = 1,2

               do 30 i = 1,n
                  do 71 j = 1,m
                     ya(j) = y(i,j,k,kk)
 71               continue
                  call polint(xa,ya,m,xmup(ii),yval,dy)
                  yy(i) = yval
                  if (k .eq. 1) then
                     w1(i,kk,ii) = yval
                  elseif (k .eq. 2) then
                     w2(i,kk,ii) = yval
                  elseif (k .eq. 3) then
                     w3(i,kk,ii) = yval
                  elseif (k .eq. 4) then
                     w4(i,kk,ii) = yval
                  endif
 30            continue

               do 21 i = 2,n-1
                  a(i) = (xx(i)-xx(i-1))/6.0D+00
                  b(i) = (xx(i+1)-xx(i-1))/3.0D+00
                  c(i) = (xx(i+1)-xx(i))/6.0+00
                  r(i) = (yy(i+1)-yy(i))/(xx(i+1)-xx(i)) -
     *                   (yy(i)-yy(i-1))/(xx(i)-xx(i-1))
 21            continue

               call tridag(a,b,c,r,u,n)
               do 22 i = 1,n
                  if (k .eq. 1) then
                     w1dp(i,kk,ii) = u(i)
                  elseif (k .eq. 2) then
                     w2dp(i,kk,ii) = u(i)
                  elseif (k .eq. 3) then
                     w3dp(i,kk,ii) = u(i)
                  elseif (k .eq. 4) then
                     w4dp(i,kk,ii) = u(i)
                  endif
 22            continue

 82         continue
 60      continue
 81   continue

      return
      end

      subroutine splin2(xa,ya,y2a,x,y,yp,ydp,nycor,nbcor,
     *                  khi,klo,kvals,ndoy,ndoyp,ndoydp)
      implicit real*8 (a-h,o-z)
      parameter (nkpts=132,n=nkpts,nbetas=3)
      dimension xa(n),ya(n,nbetas,2),y2a(n,nbetas,2)

      one = 1.0D+00

      if (kvals .eq. 0) then
         klo = 1
         khi = n
   1     if ((khi-klo) .gt. 1) then
            k = (khi+klo)/2
            if (xa(k) .gt. x) then
               khi = k
            else
               klo = k
            endif
            goto 1
         endif
      endif

      h = xa(khi) - xa(klo)
      a = (xa(khi) - x)/h
      b = (x - xa(klo))/h
      asq = a*a
      bsq = b*b

      if (ndoy .eq. 1)
     *    y = a*ya(klo,nbcor,nycor) + b*ya(khi,nbcor,nycor) +
     *       ((asq*a-a)*y2a(klo,nbcor,nycor)
     *        +(bsq*b-b)*y2a(khi,nbcor,nycor))*(h*h)/6.0D+00

      if (ndoyp .eq. 1)
     *   yp = (ya(khi,nbcor,nycor)-ya(klo,nbcor,nycor))/h -
     *        (3.0D+00*asq-one)/6.0D+00*h*y2a(klo,nbcor,nycor) +
     *        (3.0D+00*bsq-one)/6.0D+00*h*y2a(khi,nbcor,nycor)

      if (ndoydp .eq. 1)
     *   ydp = a*y2a(klo,nbcor,nycor) + b*y2a(khi,nbcor,nycor)

      return
      end

      subroutine locid2(x,xlow,xhigh,xinc,npts,ngrid,nlow,nhigh)
      implicit real*8 (a-h,o-z)

      if (x .lt. xlow) then
         ngrid = 1
         nlow = 1
         nhigh = 0
      elseif (x .gt. xhigh) then
         ngrid = npts - 1
         nlow = 0
         nhigh = 1
      else
         ngrid = idint((x-xlow)/xinc) + 1
         nlow = 0
         nhigh = 0
      endif

      return
      end

      subroutine integr(xa,ya,y2a,n,xint)
      implicit real*8 (a-h,o-z)
      dimension xa(n),ya(n),y2a(n)

      one = 1.0D+00
      one24 = one/24.0D+00

      sum = 0.0D+00
      do 10 i = 1,n-1
         diff = xa(i+1)-xa(i)
         sum = sum + diff*( (ya(i)+ya(i+1))/2.0D+00 -
     *               one24*diff*diff*(y2a(i)+y2a(i+1)) )
 10   continue

      xint = sum

      return
      end

      subroutine calcsd(x,y,vdp,a,b,c,r,n,an,c1)
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n),vdp(n),a(n),b(n),c(n),r(n)

      a(1) = 0.0D+00
      a(n) = an
      b(1) = 1.0D+00
      b(n) = 1.0D+00
      c(1) = c1
      c(n) = 0.0D+00
      r(1) = 0.0D+00
      r(n) = 0.0D+00

      do 21 i = 2,n-1
         a(i) = (x(i)-x(i-1))/6.0D+00
         b(i) = (x(i+1)-x(i-1))/3.0D+00
         c(i) = (x(i+1)-x(i))/6.0+00
         r(i) = (y(i+1)-y(i))/(x(i+1)-x(i)) -
     *          (y(i)-y(i-1))/(x(i)-x(i-1))
 21   continue

      call tridag(a,b,c,r,vdp,n)

      return
      end

      subroutine doint(xa,ya,y2a,x,y,yp,ydp,khi,klo,kvals,n,ndoy,
     *                 ndoyp,ndoydp)
      implicit real*8 (a-h,o-z)
      dimension xa(n),ya(n),y2a(n)

      one = 1.0D+00

      if (kvals .eq. 0) then
         klo = 1
         khi = n
   1     if ((khi-klo) .gt. 1) then
            k = (khi+klo)/2
            if (xa(k) .gt. x) then
               khi = k
            else
               klo = k
            endif
            goto 1
         endif
      endif

      h = xa(khi) - xa(klo)
      a = (xa(khi) - x)/h
      b = (x - xa(klo))/h
      asq = a*a
      bsq = b*b

      if (ndoy .eq. 1)
     *   y = a*ya(klo) + b*ya(khi) +
     *       ((asq*a-a)*y2a(klo)
     *        +(bsq*b-b)*y2a(khi))*(h*h)/6.0D+00

      if (ndoyp .eq. 1)
     *   yp = (ya(khi)-ya(klo))/h -
     *        (3.0D+00*asq-one)/6.0D+00*h*y2a(klo) +
     *        (3.0D+00*bsq-one)/6.0D+00*h*y2a(khi)

      if (ndoydp .eq. 1)
     *   ydp = a*y2a(klo) + b*y2a(khi)

      return
      end

      subroutine pinvar(p,prob)
      implicit real*8 (a-h,o-z)
      parameter (mxiter=100,toler=1.0D-10,nbetas=3,n=nbetas)
      dimension p(n,n),ptemp1(n,n),ptemp2(n,n),prob(n)

      call copmat(p,ptemp1,n,n)

      do 10 i = 1,mxiter
         call mult2(ptemp1,ptemp1,ptemp2,n,n,n)
         diffmx = dabs(ptemp2(1,1)-ptemp1(1,1))/ptemp1(1,1)
         do 11 j = 2,n
            diff = dabs(ptemp2(j,j)-ptemp1(j,j))/ptemp1(j,j)
            if (diff .gt. diffmx) diffmx = diff
 11      continue
         if (diffmx .lt. toler) then
            goto 999
         else
            call copmat(ptemp2,ptemp1,n,n)
         endif
 10   continue

999   do 12 i = 1,n
         prob(i) = ptemp2(i,i)
 12   continue

      return
      end

      integer function nstate(nz,neps)
      implicit real*8 (a-h,o-z)

      if ( (nz .eq. 1) .and. (neps .eq. 1) ) then
         nstate = 1
      elseif ( (nz .eq. 0) .and. (neps .eq. 1) ) then
         nstate = 2
      elseif ( (nz .eq. 1) .and. (neps .eq. 0) ) then
         nstate = 3
      elseif ( (nz .eq. 0) .and. (neps .eq. 0) ) then
         nstate = 4
      endif

      return
      end

      subroutine nrand2(sseed,nlow,nhigh,nums,n)
      implicit real*8 (a-h,o-z)
      parameter (mxiter=10000,mxswch=15)
      dimension nums(2*mxswch)

      nums(1) = nrand(sseed,nlow,nhigh)
      do 10 i = 2,n
         do 11 j = 1,mxiter
            kk = nrand(sseed,nlow,nhigh)
            match = 0
            do 12 k = 1,i-1
               if (kk .eq. nums(k)) then
                  match = 1
                  goto 998
               endif
 12         continue
998         if (match .eq. 0) goto 999
 11      continue
999      nums(i) = kk
 10   continue

      return
      end

      subroutine lawlrg(sseed,neps,nz,nagent,nswtch)
      implicit real*8 (a-h,o-z)
      parameter (unempg=0.04D+00,unempb=0.1D+00,mxswch=100000)
      dimension neps(nagent)

      one = 1.0D+00

      if (nz .eq. 1) then
         nemptg = idnint(dble(real(nagent))*(one-unempg))
      else
         nemptg = idnint(dble(real(nagent))*(one-unempb))
      endif

      numemp = 0
      do 10 i = 1,nagent
         if (neps(i) .eq. 1) numemp = numemp + 1
 10   continue

      nswtch = numemp - nemptg

      if (nswtch .gt. 0) then
         nchang = 0
         do 26 j = 1,mxswch
            k = nrand(sseed,1,nagent)
            if (neps(k) .eq. 1) then
               neps(k) = 0
               nchang = nchang + 1
            endif
            if (nchang .eq. nswtch) goto 999
 26      continue
      elseif (nswtch .lt. 0) then
         nchang = 0
         do 27 j = 1,mxswch
            k = nrand(sseed,1,nagent)
            if (neps(k) .eq. 0) then
               neps(k) = 1
               nchang = nchang + 1
            endif
            if (nchang .eq. (-nswtch)) goto 999
 27      continue
      endif

999   return
      end

      double precision function dstd3(xk1,xk2,xk3,nagnt1,nagnt2,
     *                                nagnt3,xkmean,p1,p2,p3)
      implicit real*8 (a-h,o-z)
      dimension xk1(nagnt1),xk2(nagnt2),xk3(nagnt3)

      sum1 = 0.0D+00
      do 10 i = 1,nagnt1
         sum1 = sum1 + xk1(i)*xk1(i)
 10   continue
      eksq1 = sum1/dble(real(nagnt1))

      sum2 = 0.0D+00
      do 11 i = 1,nagnt2
         sum2 = sum2 + xk2(i)*xk2(i)
 11   continue
      eksq2 = sum2/dble(real(nagnt2))

      sum3 = 0.0D+00
      do 12 i = 1,nagnt3
         sum3 = sum3 + xk3(i)*xk3(i)
 12   continue
      eksq3 = sum3/dble(real(nagnt3))

      dstd3 = dsqrt(p1*eksq1 + p2*eksq2 + p3*eksq3 - xkmean*xkmean)

      return
      end

      subroutine domean(xk1,xk2,xk3,nagnt1,nagnt2,nagnt3,wght1,
     *                  wght2,wght3,xkmean)
      implicit real*8 (a-h,o-z)
      dimension xk1(nagnt1),xk2(nagnt2),xk3(nagnt3)

      xnagnt = dble(real(nagnt1+nagnt2+nagnt3))

      sum = 0.0D+00
      do 82 j = 1,nagnt1
         sum = sum + xk1(j)
 82   continue
      sumtot = wght1*sum

      sum = 0.0D+00
      do 182 j = 1,nagnt2
         sum = sum + xk2(j)
182   continue
      sumtot = sumtot + wght2*sum

      sum = 0.0D+00
      do 282 j = 1,nagnt3
         sum = sum + xk3(j)
282   continue
      sumtot = sumtot + wght3*sum

      xkmean = sumtot/xnagnt

      return
      end

      subroutine dopct(xk1,xk2,xk3,neps1,neps2,neps3,nagnt1,nagnt2,
     *                 nagnt3,wght1,wght2,wght3,pctemp,p0,pcon)
      implicit real*8 (a-h,o-z)
      parameter (nbetas=3,xkbor=-2.4D+00,toler=1.0D-06)
      dimension xk1(nagnt1),xk2(nagnt2),xk3(nagnt3),pctemp(nbetas),
     *          neps1(nagnt1),neps2(nagnt2),neps3(nagnt3)

      xnagnt = dble(real(nagnt1+nagnt2+nagnt3))

      sum1 = 0.0D+00
      sum2 = 0.0D+00

      msum = 0
      do 82 j = 1,nagnt1
         if (xk1(j) .le. 0.0D+00) sum1 = sum1 + wght1
         if (dabs(xk1(j)-xkbor) .le. toler) sum2 = sum2 + wght1
         msum = msum + neps1(j)
 82   continue
      pctemp(1) = dble(real(msum))/dble(real(nagnt1))

      msum = 0
      do 182 j = 1,nagnt2
         if (xk2(j) .le. 0.0D+00) sum1 = sum1 + wght2
         if (dabs(xk2(j)-xkbor) .le. toler) sum2 = sum2 + wght2
         msum = msum + neps2(j)
182   continue
      pctemp(2) = dble(real(msum))/dble(real(nagnt2))

      msum = 0
      do 282 j = 1,nagnt3
         if (xk3(j) .le. 0.0D+00) sum1 = sum1 + wght3
         if (dabs(xk3(j)-xkbor) .le. toler) sum2 = sum2 + wght3
         msum = msum + neps3(j)
282   continue
      pctemp(3) = dble(real(msum))/dble(real(nagnt3))

      p0 = sum1/xnagnt
      pcon = sum2/xnagnt

      return
      end

      double precision function xmax3(x1,x2,x3)
      implicit real*8 (a-h,o-z)

      if ( (x1 .gt. x2) .and. (x1 .gt. x3) ) then
         xmax = x1
      elseif ( (x2 .gt. x1) .and. (x2 .gt. x3) ) then
         xmax = x2
      else
         xmax = x3
      endif

      xmax3 = xmax

      return
      end

      subroutine ludcmp(a,n,np,indx,d,nsing)
      implicit real*8 (a-h,o-z)
      parameter (nmax=1000,tiny=1.0D-15)
      dimension a(np,np),indx(n),vv(nmax)

      d = 1.0D+00
      nsing = 0

      do 12 i = 1,n
         aamax = 0.0D+00
         do 11 j = 1,n
            if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
 11      continue
C         if (aamax .eq. 0.0D+00) pause 'Singular matrix'
         if (dabs(aamax) .le. 1.0D-13) then
            nsing = 1
            goto 9999
         endif
         vv(i) = 1.0D+00/aamax
 12   continue

      do 19 j = 1,n
         do 14 i = 1,j-1
            sum = a(i,j)
            do 13 k = 1,i-1
               sum = sum - a(i,k)*a(k,j)
 13         continue
            a(i,j) = sum
 14      continue
         aamax = 0.0D+00
         do 16 i = j,n
            sum = a(i,j)
            do 15 k = 1,j-1
               sum = sum - a(i,k)*a(k,j)
 15         continue
            a(i,j) = sum
            dum = vv(i)*abs(sum)
            if (dum .ge. aamax) then
               imax = i
               aamax = dum
            endif
 16      continue
         if (j .ne. imax) then
            do 17 k = 1,n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
 17         continue
            d = -d
            vv(imax) = vv(j)
         endif
         indx(j) = imax
         if (a(j,j) .eq. 0.0D+00) a(j,j) = tiny
         if (j .ne. n) then
            dum = 1.0D+00/a(j,j)
            do 18 i = j+1,n
               a(i,j) = a(i,j)*dum
 18         continue
         endif
 19   continue

9999  return
      end

      subroutine lubksb(a,n,np,indx,b)
      implicit real*8 (a-h,o-z)
      dimension a(np,np),indx(n),b(n)

      ii = 0

      do 12 i = 1,n
         ll = indx(i)
         sum = b(ll)
         b(ll) = b(i)
         if (ii .ne. 0) then
            do 11 j = ii,i-1
               sum = sum - a(i,j)*b(j)
 11         continue
         else if (sum .ne. 0.0D+00) then
            ii = i
         endif
         b(i) = sum
 12   continue

      do 14 i = n,1,-1
         sum = b(i)
         if (i .lt. n) then
            do 13 j = i+1,n
               sum = sum - a(i,j)*b(j)
 13         continue
         endif
         b(i) = sum/a(i,i)
 14   continue

      return
      end

      subroutine inv4(a,b,nsing)
      implicit real*8 (a-h,o-z)
      parameter (n=4)
      dimension a(n,n),b(n,n),atemp(n,n),eye(n,n),indx(n)

      call inimat(eye,n,n)
      do 10 i = 1,n
         eye(i,i) = 1.0D+00
 10   continue

      call copmat(a,atemp,n,n)

      call ludcmp(atemp,n,n,indx,d,nsing)
      if (nsing .eq. 1) goto 9999
      do 11 j = 1,n
         call lubksb(atemp,n,n,indx,eye(1,j))
 11   continue

      call copmat(eye,b,n,n)

9999  return
      end

      subroutine inimat(a,n,k)
      implicit real*8 (a-h,o-z)
      dimension a(n,k)

      do 10 i = 1,n
         do 20 j = 1,k
            a(i,j) = 0.0D+00
 20      continue
 10   continue

      return
      end

      subroutine copmat(a,b,n,k)
      implicit real*8 (a-h,o-z)
      dimension a(n,k),b(n,k)

      do 10 i = 1,n
         do 20 j = 1,k
            b(i,j) = a(i,j)
 20      continue
 10   continue

      return
      end

      subroutine mult2(a,b,c,n,k,l)
      implicit real*8 (a-h,o-z)
      dimension a(n,k),b(k,l),c(n,l)

      if (l .eq. 1) then
         do 20 i = 1,n
            c(i,1) = 0.0D+00
            do 22 m = 1,k
               c(i,1) = c(i,1) + a(i,m)*b(m,1)
 22         continue
 20      continue
      else
         do 10 i = 1,n
            do 11 j = 1,l
               c(i,j) = 0.0D+00
               do 12 m = 1,k
                  c(i,j) = c(i,j) + a(i,m)*b(m,j)
 12            continue
 11         continue
 10      continue
      endif

      return
      end

      double precision function dmean(x,n)
      implicit real*8 (a-h,o-z)
      dimension x(n)

      sum = 0.0D+00
      do 10 i = 1,n
         sum = sum + x(i)
 10   continue

      dmean = sum/dble(real(n))

      return
      end

      double precision function dstd(x,n)
      implicit real*8 (a-h,o-z)
      dimension x(n)

      xmean = dmean(x,n)

      sum = 0.0D+00
      do 10 i = 1,n
         sum = sum + (x(i)-xmean)*(x(i)-xmean)
 10   continue

      dstd = dsqrt(sum/dble(real(n)))

      return
      end

      double precision function dmax(x,n)
      implicit real*8 (a-h,o-z)
      dimension x(n)

      xmax = x(1)
      do 10 i = 2,n
         if (x(i) .gt. xmax) xmax = x(i)
 10   continue

      dmax = xmax

      return
      end

      double precision function dmin(x,n)
      implicit real*8 (a-h,o-z)
      dimension x(n)

      xmin = x(1)
      do 10 i = 2,n
         if (x(i) .lt. xmin) xmin = x(i)
 10   continue

      dmin = xmin

      return
      end

      integer function multpl(i,m)
      implicit real*8 (a-h,o-z)
      parameter (toler=1.0D-12)

      remain = dble(real(i))/dble(real(m)) - i/m

      if (dabs(remain) .lt. toler) then
         multpl = 1
      else
         multpl = 0
      endif

      return
      end

C dseed: inital seed (dseed is automatically updated)
C nr: dimension of the double precision array r
C r: real*8 r(nr)
C GGUBS generates nr independent uniform random variables
C GGNML generates nr independent normal variables


      SUBROUTINE GGUBS(DSEED,NR,R)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   GGUBS GENERATES NR SINGLE PRECISION RANDOM VARIATES UNIFORM C
C ON (0,1) BY A LINEAR CONGRUENTIAL SCHEME.  THIS ROUTINE IS    C
C DEPENDENT ON MACHINE WORD SIZE.                               C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   LAST MODIFIED                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C   ON ENTRY                                                    C
C       DSEED   DOUBLE PRECISION                                C
C               SEED FOR GENERATOR                              C
C       NR      INTEGER                                         C
C               NUMBER OF VARIATES TO GENERATE                  C
C   ON RETURN                                                   C
C       R       REAL (NR)                                       C
C               SINGLE PRECISION ARRAY CONTAINING THE VARIATES  C
C   GGUBS CALLS                                                 C
C               DMOD                                            C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
      INTEGER            NR
      REAL*8             R(NR)
      DOUBLE PRECISION   DSEED
C
C                              LOCAL
C
      INTEGER            I
      DOUBLE PRECISION   D2P31M,D2P31,DMULTX
C
C                              MACHINE CONSTANTS
C                              D2P31M=(2**31) - 1
C                              D2P31 =(2**31)(OR AN ADJUSTED VALUE)
C
      DATA               D2P31M/2147483647.D0/
      DATA               D2P31/2147483711.D0/
      DATA               DMULTX/16807.0D+00/
C
      DO 5 I=1,NR
         DSEED=DMOD(DMULTX*DSEED,D2P31M)
         R(I) =DSEED / D2P31
  5   CONTINUE
C
C                               END OF GGUBS
C
      RETURN
      END
C===============================================================C
C===============================================================C
      SUBROUTINE GGNML(DSEED,NR,R)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   GGNML GENERATES NR SINGLE PRECISION N(0,1) RANDOM VARIATES  C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   LAST MODIFIED                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C   ON ENTRY                                                    C
C       DSEED   DOUBLE PRECISION                                C
C               SEED FOR GENERATOR                              C
C       NR      INTEGER                                         C
C               NUMBER OF VARIATES TO GENERATE                  C
C   ON RETURN                                                   C
C       R       REAL (NR)                                       C
C               SINGLE PRECISION ARRAY CONTAINING THE VARIATES  C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
      INTEGER            NR
      REAL*8             R(NR)
      DOUBLE PRECISION   DSEED
C                              LOCAL
      INTEGER             IER
C
C                              GET NR RANDOM NUMBERS
C                              UNIFORM (0,1)
C
      CALL GGUBS(DSEED,NR,R)
C
C                              TRANSFORM EACH UNIFORM DEVIATE
C
      DO 5 I=1,NR
         CALL MDNRIS(R(I),R(I),IER)
    5 CONTINUE
C
C                               END OF GGNML
C
      RETURN
      END


      SUBROUTINE MDNRIS (P,Y,IER)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   LAST MODIFIED                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL*8             P,Y
      INTEGER            IER
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL*8             EPS,G0,G1,G2,G3,H0,H1,H2,A,W,WI,SN,SD
      REAL*8             SIGMA,SQRT2,X,XINF
      DATA               XINF/1.7014E+38/
      DATA               SQRT2/1.414214/
      DATA               EPS/1.1921E-07/
      DATA               G0/.1851159E-3/,G1/-.2028152E-2/
      DATA               G2/-.1498384/,G3/.1078639E-1/
      DATA               H0/.9952975E-1/,H1/.5211733/
      DATA               H2/-.6888301E-1/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (P .GT. 0.0 .AND. P .LT. 1.0) GO TO 5
      IER = 129

      if (p .lt. 0.0D+00) then
         sigma = -1.0D+00
      else
         sigma = 1.0D+00
      endif

C      SIGMA = SIGN(1.0,P)
C      write(6,666) p,sigma
C666   format(' Sign #1: ',2f15.8)
C      pause

      Y = SIGMA * XINF
      GO TO 20
    5 IF(P.LE.EPS) GO TO 10
      X = 1.0 -(P + P)
      CALL MERFI (X,Y,IER)
      Y = -SQRT2 * Y
      GO TO 20
C                                  P TOO SMALL, COMPUTE Y DIRECTLY
   10 A = P+P
      W = SQRT(-dLOG(A+(A-A*A)))
C                                  USE A RATIONAL FUNCTION IN 1./W
      WI = 1./W
      SN = ((G3*WI+G2)*WI+G1)*WI
      SD = ((WI+H2)*WI+H1)*WI+H0
      Y = W + W*(G0+SN/SD)
      Y = -Y*SQRT2
C                               END OF MDNRIS
  20  RETURN
      END
C===============================================================C
      SUBROUTINE MERFI (P,Y,IER)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   LAST MODIFIED                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL*8             P,Y
      INTEGER            IER
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL*8             A,B,X,Z,W,WI,SN,SD,F,Z2,RINFM,A1,A2,A3,B0,B1,
     *                   B2,B3,C0,C1,C2,C3,D0,D1,D2,E0,E1,E2,E3,F0,F1,
     *                   F2,G0,G1,G2,G3,H0,H1,H2,SIGMA
      DATA               A1/-.5751703/,A2/-1.896513/,A3/-.5496261E-1/
      DATA               B0/-.1137730/,B1/-3.293474/,B2/-2.374996/
      DATA               B3/-1.187515/
      DATA               C0/-.1146666/,C1/-.1314774/,C2/-.2368201/
      DATA               C3/.5073975E-1/
      DATA               D0/-44.27977/,D1/21.98546/,D2/-7.586103/
      DATA               E0/-.5668422E-1/,E1/.3937021/,E2/-.3166501/
      DATA               E3/.6208963E-1/
      DATA               F0/-6.266786/,F1/4.666263/,F2/-2.962883/
      DATA               G0/.1851159E-3/,G1/-.2028152E-2/
      DATA               G2/-.1498384/,G3/.1078639E-1/
      DATA               H0/.9952975E-1/,H1/.5211733/
      DATA               H2/-.6888301E-1/
      DATA               RINFM/1.7014E+38/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      X = P

      if (x .lt. 0.0D+00) then
         sigma = -1.0D+00
      else
         sigma = 1.0D+00
      endif

C      SIGMA = SIGN(1.0,X)
C      write(6,666) x,sigma
C666   format(' Sign #2: ',2f15.8)
C      pause

C                                  TEST FOR INVALID ARGUMENT
      IF (.NOT.(X.GT.-1. .AND. X.LT.1.)) GO TO 30
      Z = ABS(X)
      IF (Z.LE. .85) GO TO 20
      A = 1.-Z
      B = Z
C                                  REDUCED ARGUMENT IS IN (.85,1.),
C                                     OBTAIN THE TRANSFORMED VARIABLE
    5 W = SQRT(-dLOG(A+A*B))
      IF (W.LT.2.5) GO TO 15
      IF (W.LT.4.) GO TO 10
C                                  W GREATER THAN 4., APPROX. F BY A
C                                     RATIONAL FUNCTION IN 1./W
      WI = 1./W
      SN = ((G3*WI+G2)*WI+G1)*WI
      SD = ((WI+H2)*WI+H1)*WI+H0
      F = W + W*(G0+SN/SD)
      GO TO 25
C                                  W BETWEEN 2.5 AND 4., APPROX. F
C                                     BY A RATIONAL FUNCTION IN W
   10 SN = ((E3*W+E2)*W+E1)*W
      SD = ((W+F2)*W+F1)*W+F0
      F = W + W*(E0+SN/SD)
      GO TO 25
C                                  W BETWEEN 1.13222 AND 2.5, APPROX.
C                                     F BY A RATIONAL FUNCTION IN W
   15 SN = ((C3*W+C2)*W+C1)*W
      SD = ((W+D2)*W+D1)*W+D0
      F = W + W*(C0+SN/SD)
      GO TO 25
C                                  Z BETWEEN 0. AND .85, APPROX. F
C                                     BY A RATIONAL FUNCTION IN Z
   20 Z2 = Z*Z
      F = Z+Z*(B0+A1*Z2/(B1+Z2+A2/(B2+Z2+A3/(B3+Z2))))
C                                  FORM THE SOLUTION BY MULT. F BY
C                                     THE PROPER SIGN
   25 Y = SIGMA*F
      IER = 0
      GO TO 40
C                                  ERROR EXIT. SET SOLUTION TO PLUS
C                                     (OR MINUS) INFINITY
   30 IER = 129
      Y = SIGMA * RINFM
C                               END OF MERFI
   40 RETURN
      END
C===============================================================C

      SUBROUTINE DDNOR(Y,GAUSS)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 P(6), Q(5), A(9), B(8), C(5), D(4)
      DATA P(1)/-6.58749161529837803157D-04/,
     1     P(2)/-1.60837851487422766278D-02/,
     2     P(3)/-1.25781726111229246204D-01/,
     3     P(4)/-3.60344899949804439429D-01/,
     4     P(5)/-3.05326634961232344035D-01/,
     5     P(6)/-1.63153871373020978498D-02/
      DATA Q(1)/2.33520497626869185443D-03/,
     1     Q(2)/6.05183413124413191178D-02/,
     2     Q(3)/5.27905102951428412248D-01/,
     3     Q(4)/1.87295284992346047209D00/,
     4     Q(5)/2.56852019228982242072D00/
      DATA A(1)/1.23033935479799725272D03/,
     1     A(2)/2.05107837782607146532D03/,
     2     A(3)/1.71204761263407058314D03/,
     3     A(4)/8.81952221241769090411D02/,
     4     A(5)/2.98635138197400131132D02/,
     5     A(6)/6.61191906371416294775D01/,
     6     A(7)/8.88314979438837594118D00/,
     7     A(8)/5.64188496988670089180D-01/,
     8     A(9)/2.15311535474403846343D-08/
      DATA B(1)/1.23033935480374942043D03/,
     1     B(2)/3.43936767414372163696D03/,
     2     B(3)/4.36261909014324715820D03/,
     3     B(4)/3.29079923573345962678D03/,
     4     B(5)/1.62138957456669018874D03/,
     5     B(6)/5.37181101862009857509D02/,
     6     B(7)/1.17693950891312499305D02/,
     7     B(8)/1.57449261107098347253D01/
      DATA C(1)/3.209377589138469472562D03/,
     1     C(2)/3.774852376853020208137D02/,
     2     C(3)/1.138641541510501556495D02/,
     3     C(4)/3.161123743870565596947D00/,
     4     C(5)/1.857777061846031526730D-01/
      DATA D(1)/2.844236833439170622273D03/,
     1     D(2)/1.282616526077372275645D03/,
     2     D(3)/2.440246379344441733056D02/,
     3     D(4)/2.360129095234412093499D01/
      DATA ORPI/.5641895835477562869483D0/,
     1   ROOT2/.70710678118654752440083D0/
C  THIS SUBROUTINE USES CODY'S METHOD TO EVALUATE THE CUMULATIVE
C NORMAL DISTRIBUTION. IT IS PROBABLY ACCURATE TO 19 OR 20
C SIGNIFICANT DIGITS. IT WAS WRITTEN BY JAMES MACKINNON LATE IN
C 1977, BASED ON THE CODY ARTICLE REFERRED TO IN THE DOCUMENTATION
C FOR IMSL SUBROUTINE MDNOR.
      ISW = 1
      IF (Y.LT.-16.D0) Y = -16.D0
      IF (Y.GT.16.D0) Y = 16.D0
      X = -Y*ROOT2
      IF(X.GT.0.D0) GO TO 1
      IF(X.LT.0.D0) GO TO 2
      GAUSS = .5D0
      RETURN
    2 CONTINUE
      X = - X
      ISW = -1
    1 CONTINUE
      IF(X.LT..477D0) GO TO 10
      IF(X.LE.4.D0) GO TO 20
C  EVALUATE ERFC FOR X.GT.4.0
      X2 = X*X
      XM2 = 1.D0/X2
      XM4 = XM2*XM2
      XM6 = XM4*XM2
      XM8 = XM4*XM4
      XM10 = XM6*XM4
      TOP = P(1) + P(2)*XM2 + P(3)*XM4 + P(4)*XM6 + P(5)*XM8 + P(6)*XM10
      BOT = Q(1) + Q(2)*XM2 + Q(3)*XM4 + Q(4)*XM6 + Q(5)*XM8 + XM10
      CRAP = ORPI + TOP/(BOT*X2)
      ERFC = DEXP(-X2)*CRAP/X
C
      IF(ISW.EQ.-1) ERFC = 2.D0 - ERFC
      GAUSS = ERFC*.5D0
      RETURN
   20 CONTINUE
C  EVALUATE ERFC FOR .477.LT.X.LE.4.0
      X2 = X*X
      X3 = X2*X
      X4 = X2*X2
      X5 = X3*X2
      X6 = X3*X3
      X7 = X3*X4
      X8 = X4*X4
      TOP = A(1) + A(2)*X + A(3)*X2 + A(4)*X3 + A(5)*X4 + A(6)*X5 +
     1 A(7)*X6 + A(8)*X7 + A(9)*X8
      BOT = B(1) + B(2)*X + B(3)*X2 + B(4)*X3 + B(5)*X4 + B(6)*X5 +
     1 B(7)*X6 + B(8)*X7 + X8
      ERFC = DEXP(-X2)*TOP/BOT
C
      IF(ISW.EQ.-1) ERFC = 2.D0 - ERFC
      GAUSS = ERFC*.5D0
      RETURN
   10 CONTINUE
C  EVALUATE ERF FOR X.LT..477
      X2 = X*X
      X4 = X2*X2
      X6 = X4*X2
      X8 = X4*X4
      TOP = C(1) + C(2)*X2 + C(3)*X4 + C(4)*X6 + C(5)*X8
      BOT = D(1) + D(2)*X2 + D(3)*X4 + D(4)*X6 + X8
      ERF = X*TOP/BOT
C
      ERF = ERF*ISW
      ERFC = 1.D0 - ERF
      GAUSS = ERFC*.5D0
      RETURN
      END

      integer function nrand(sseed,nlow,nhigh)
      implicit real*8 (a-h,o-z)

      call ggubs(sseed,1,uni)

      nrand = nlow + idint(uni*dble(real(nhigh-nlow+1)))

      return
      end

      subroutine ols2(ydata,xdata,ndata,bhat,rsq,stderr,npts)
      implicit real*8 (a-h,o-z)
      parameter (np=2)
      dimension xdata(ndata,np),ydata(ndata),xpx(np,np),
     *          xpxinv(np,np),xpy(np,1),bhat(np,1),temp1(np,1),
     *          temp2(np,1),xt(np,1),xtp(1,np),xtxtp(np,np),
     *          temp(np,np)

      call inimat(xpx,np,np)
      do 51 i = 1,npts
         do 52 j = 1,np
            xt(j,1) = xdata(i,j)
 52      continue
         call tr(xt,xtp,np,1)
         call mult2(xt,xtp,xtxtp,np,1,np)
         call add2(xpx,xtxtp,temp,np,np)
         call copmat(temp,xpx,np,np)
 51   continue
      call inv2(xpx,xpxinv,nsing)

      call inimat(xpy,np,1)
      do 10 i = 1,npts
         do 20 j = 1,np
            temp1(j,1) = xdata(i,j)*ydata(i)
 20      continue
         call add2(xpy,temp1,temp2,np,1)
         call copmat(temp2,xpy,np,1)
 10   continue

      call mult2(xpxinv,xpy,bhat,np,np,1)

      sse = 0.0D+00
      sumy = 0.0D+00
      do 12 i = 1,npts
         diff = ydata(i)
         do 14 j = 1,np
            diff = diff - bhat(j,1)*xdata(i,j)
 14      continue
         sse = sse + diff*diff
         sumy = sumy + ydata(i)
 12   continue
      ymean = sumy/dble(real(npts))
      stderr = dsqrt(sse/dble(real(npts)))

      tss = 0.0D+00
      do 13 i = 1,npts
         tss = tss + (ydata(i)-ymean)*(ydata(i)-ymean)
 13   continue

      rsq = 1.0D+00 - sse/tss

      return
      end

      subroutine locidx(x,xpts,npts,idx,nlow,nhigh,niter)
      implicit real*8 (a-h,o-z)
      parameter (mxiter=50)
      dimension xpts(npts)

      nhigh = 0
      nlow = 0

      if (x .ge. xpts(npts)) then
         nhigh = 1
         idxtmp = 0
         niter = 0
      elseif (x .lt. xpts(1)) then
         nlow = 1
         idxtmp = 0
         niter = 0
      else

         ilow = 1
         ihigh = npts
         do 10 i = 1,mxiter
            icur = (ilow+ihigh)/2
            if (x .ge. xpts(icur+1)) then
               ilow = icur
            elseif (x .lt. xpts(icur)) then
               ihigh = icur
            else
               idxtmp = icur
               niter = i
               goto 999
            endif
 10      continue

      endif

999   idx = idxtmp

      return
      end

      subroutine add2(a,b,c,n,k)
      implicit real*8 (a-h,o-z)
      dimension a(n,k),b(n,k),c(n,k)

      do 10 i = 1,n
         do 20 j = 1,k
            c(i,j) = a(i,j) + b(i,j)
 20      continue
 10   continue

      return
      end

      subroutine inv2(a,b,nsing)
      implicit real*8 (a-h,o-z)
      parameter (n=2)
      dimension a(n,n),b(n,n),atemp(n,n),eye(n,n),indx(n)

      call inimat(eye,n,n)
      do 10 i = 1,n
         eye(i,i) = 1.0D+00
 10   continue

      call copmat(a,atemp,n,n)

      call ludcmp(atemp,n,n,indx,d,nsing)
      if (nsing .eq. 1) goto 9999
      do 11 j = 1,n
         call lubksb(atemp,n,n,indx,eye(1,j))
 11   continue

      call copmat(eye,b,n,n)

9999  return
      end

      subroutine tr(a,b,n,k)
      implicit real*8 (a-h,o-z)
      dimension a(n,k),b(k,n)

      do 10 i = 1,n
         do 20 j = 1,k
            b(j,i) = a(i,j)
 20      continue
 10   continue

      return
      end
