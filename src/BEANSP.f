c subroutine to include censored observations.
      subroutine mbeansp(jj,nx,nn,zl,zr,x,y,a,b,sigma,day,enc,tha,
     +   covar,n0,ntotal,q,del,DIC,Dbar,pd,trace1,trace2,trace3,trace4,
     +   venc,vday,vsigma,vcovar,vdel,vq,asr,vasr,casr,vcasr,sr1,vsr,
     +   mutha,vartha,trace5,mi,WAIC)
      implicit   real*8 (a-h,o-z)
      implicit   integer*4 (i-n)

      integer ig,n0,ntotal,nn,nx,ntha
      integer jj,iu,idown,iup,j1,nn1
      integer iseed, ir(1) ! for randomization
      integer nnmax
      parameter (nnmax=20000)
      parameter (ntha=2)
      integer mi1(nnmax), tm1(nnmax), mi(nn)
      integer u(nn), y(nn), t(nn),zl(nn),zr(nn)
      integer ul(nn),ur(nn),tl(nn),tr(nn)

      real*8 outs(1), ms,mf,pys,py1,py2
      real*8 scale1,temp,temp1,temp2,logliks,loglikm,slogliks
      real*8 del(jj),q(nn,jj+1),covar(nx),p(jj)
      real*8 sdel(jj),sdel2(jj),vdel(jj)
      real*8 sq(nn,jj+1),sq2(nn,jj+1),vq(nn,jj+1)
      real*8 stha(ntha),stha2(ntha),vtha(ntha)
      real*8 day(jj),sday(jj),sday2(jj),vday(jj)
      real*8 enc(jj-1),senc(jj-1),senc2(jj-1),venc(jj-1)
      real*8 AA(jj,jj), EE(jj-1,jj-1)
      real*8 sr1(nn,jj),ssr(nn,jj),ssr2(nn,jj),vsr(nn,jj)
      real*8 csr(nn,jj), scsr(nn,jj),scsr2(nn,jj),vcsr(nn,jj)
      real*8 asr(jj),sasr(jj),sasr2(jj),vasr(jj)
      real*8 casr(jj),scasr(jj),scasr2(jj),vcasr(jj)
      real*8 scovar(nx),scovar2(nx),vcovar(nx)
      real*8 sigma(2),a(2),b(2),x(nn,nx)
      real*8 ssigma(2),vsigma(2),ssigma2(2)
      real*8 x2(nnmax), x3(nnmax), x4(nnmax)
      real*8 te1(nnmax), te2(nnmax)
      real*8 td1(nnmax), td2(nnmax)
      real*8 tha(ntha), mutha(ntha), vartha(ntha)
      real*8 tha1(ntha),mutha1(ntha),vartha1(ntha)
      real*8 trace1(ntotal,jj-1),trace2(ntotal,jj)
      real*8 trace3(ntotal,nx),trace4(ntotal,2),trace5(ntotal,ntha)
      real*8 WAIC
      real*8 lpd,pwaic
      real*8 spost(nn),logpost(nn),logpost2(nn),vlogpost(nn)
      real*8 post(nn)


      integer trace6(ntotal,nn)



      COMMON /torandom/IXR,IYR,IZR
      common /toencm/ emu,sigmae,ne
      common /todaym/ dmu,sigmad,nd
      common /toxm/ varx,x1
      common /tosub1m/ te1,te2
      common /tosub2m/ td1,td2
      common /tosub3m/ x2,x3,x4
      common /tosub4m/ mi1,tm1
      common /tosub5m/ tha1,mutha1,vartha1
      common /tosub6m/ tm0
      common /inputm/ nn1


      external RANDOM, eval_enc_m, smpl_enc_m, eval_day_m, smpl_day_m,
     +   eval_tha, smpl_tha, eval_tha2, smpl_tha2
      varx=10.0d0
      nn1=nn
      tha1=tha
      mutha1=mutha
      vartha1=vartha
      iseed=1976629
      IXR=3951
      IYR=1234
      IZR=6345

   	  AA=0.0D0                !!! IAR(2) precision matrix for day
      DO I=1,jj
         DO J=1,jj
         IF ((I==J).AND.((I==1).OR.(I==jj))) AA(I,J)=1.0D0
         IF ((I==J).AND.((I==2).OR.(I==jj-1))) AA(I,J)=5.0D0
         IF ((I==J).AND.((I>2).AND.(I<jj-1))) AA(I,J)=6.0D0
         IF (ABS(I-J)==1) AA(I,J)=-4.0D0
         IF (ABS(I-J)==2) AA(I,J)=1.0D0
         ENDDO
      ENDDO
      AA(1,2)=-2.0D0
      AA(2,1)=-2.0D0
      AA(jj-1,jj)=-2.0D0
      AA(jj,jj-1)=-2.0D0

      EE=0.0D0                !!! IAR(2) precision matrix for enc
      DO I=1,jj-1
         DO J=1,jj-1
         IF ((I==J).AND.((I==1).OR.(I==jj-1))) EE(I,J)=1.0D0
         IF ((I==J).AND.((I==2).OR.(I==jj-2))) EE(I,J)=5.0D0
         IF ((I==J).AND.((I>2).AND.(I<jj-2))) EE(I,J)=6.0D0
         IF (ABS(I-J)==1) EE(I,J)=-4.0D0
         IF (ABS(I-J)==2) EE(I,J)=1.0D0
         ENDDO
      ENDDO
      EE(1,2)=-2.0D0
      EE(2,1)=-2.0D0
      EE(jj-2,jj-1)=-2.0D0
      EE(jj-1,jj-2)=-2.0D0

      mi1(1:nn)=mi
      do j=1,nx
         xmean=sum(x(1:nn,j))/(nn*1.0d0)
         x(1:nn,j)=x(1:nn,j)-xmean
      enddo
c**************************************************************
c     the following accumulate in the Gibbs sampling
c     they have to be initialized for every new estimation
c     specify the prior hyperparameters
c***************************************************************

   	  sdel=0.0d0; sdel2=0.0d0; sq=0.0d0; sq2=0.0d0
      ssr=0.0d0; ssr2=0.0d0; ssigma=0.0d0; ssigma2=0.0d0
      sasr=0.0d0; sasr2=0.0d0; scasr=0.0d0; scasr2=0.0d0
      sday=0.0d0; sday2=0.0d0; senc=0.0d0; senc2=0.0d0
      scovar=0.0d0; scovar2=0.0d0; stha=0.0d0; stha2=0.0d0;
      slogliks=0.0d0;
c***********************  assign initial value
      q=1.0d0/(jj+1)
      del(1:jj)=1.0d0/jj
      ur(1:nn1)=int(10)
      do k=1,nn
c*************** updating Y_k
         if(mi(k).eq.1) then
            y(k)=0
            ur(k)=jj-zl(k)+1
            ul(k)=jj-zr(k)+1     !ul,ur are generated with S fate
            py1=q(k,jj+1)*sum(del(ul(k):ur(k)))*
     +      dexp(tha(1)+tha(2))/(1+dexp(tha(1)+tha(2)))
            py2=0.0d0
  		    do iu=1,ur(k)        !ul=1,ur are with F fate
               idown=iu+zl(k)-1    !idown,iup are t_k interval given u_k=iu
               iup=iu+zr(k)-1
               if (iup .gt. jj) then
                  iup=jj
                  if (idown .gt. jj) idown=jj
               endif
               py2=py2+del(iu)*sum(q(k,idown:iup))
            enddo
            py2=dexp(tha(1))/(1+dexp(tha(1)))*py2
            pys=py1/(py1+py2)
            rn=ran1(iseed)
            if(rn .lt. pys) y(k)=1
   	     endif
c*************** updating U_k
         ur(k)=jj-zl(k)+1
         ul(k)=jj-zr(k)+1
         if (y(k).eq.0)  ul(k)=1
         if (ul(k) < 1)   ul(k)=1
         if (ur(k) < 1)   ur(k)=1
         if (ur(k) .eq. ul(k)) then
            u(k)=ul(k)
         else
            ir = myund(ur(k)-ul(k), iseed)
            u(k)=ul(k)+ir(1)
         endif
c*************** updating T_k
         p=0.0d0
         if (y(k) .eq. 0) then
            idown=u(k)+zl(k)-1
            iup=u(k)+zr(k)-1
            if (iup .gt. jj) then
               iup=jj
               if (idown .gt. jj) idown=jj
            endif
            temp=0.0d0
            do j1=idown,iup
               temp=temp+q(k,j1)
               p(j1)=temp
            enddo
            do j1=idown,iup-1
               p(j1)=p(j1)/temp
            enddo
            p(iup)=1.0d0
   			rn=ran1(iseed)
            if (rn .lt. p(idown)) then
               t(k)=idown
            else
               do j1=idown,iup-1
                  if((rn .ge. p(j1)).and.(rn .lt. p(j1+1))) then
                    t(k)=j1+1
                  endif
               enddo
            endif
         else
		    t(k)=jj+1
         endif
      enddo
c*************** initial values for computing waic
      lpd=0
      pwaic=0
      do k=1,nn
          post(k)=0
          spost(k)=0
          logpost(k)=0
          logpost2(k)=0
          vlogpost(k)=0
      enddo




      do 16 ig=1,ntotal
!*****************  updating encounter effect, enc()
       do j=1,jj-1
          ne=0
          do k=1,nn
             if(u(k).eq.j) ne=ne+1
   			 te1(k)=0.0d0
	         do i=1,jj
                if (i.eq.j) then
				   do j1=i,jj+1
                      te1(k)=te1(k)+q(k,j1)
				   enddo
   				endif
             enddo
   			 te2(k)=0.0d0
  			 do j1=jj,jj+1
			    te2(k)=te2(k)+q(k,j1)
	         enddo
             do i=1,jj-1
			    temp=0.0d0
			    do j1=i,jj+1
  				temp=temp+q(k,j1)
	            enddo
                te2(k)=te2(k)+dexp(enc(i))*temp
	         enddo
   			 te2(k)=te2(k)-dexp(enc(j))*te1(k)
          enddo

          sigmae=sigma(1)/EE(j,j)
          emu=0.0d0
          do je=1,jj-1
             emu=emu+EE(j,je)*enc(je)
          enddo
          emu=-(emu-EE(j,j)*enc(j))/EE(j,j)
          call smpl_enc_m(oute)
          enc(j)=oute
          trace1(ig,j)=enc(j)
       enddo
       deno=1.0d0
       do j=1,jj-1
          deno=deno+dexp(enc(j))
       enddo
       do j=1,jj-1
          del(j)=dexp(enc(j))/deno
       enddo
       del(jj)=1.0-sum(del(1:(jj-1)))
c*************** updating q(j), failure rate
       do j=1,jj
          nd=0
          do k=1,nn
             if(t(k).eq.j) nd=nd+1
   			 td1(k)=0.0d0
	         do i=1,jj
                if (i.eq.j) then
				   do i1=1,i
                      td1(k)=td1(k)+del(i1)
				   enddo
   				endif
             enddo
             td1(k)=td1(k)*dexp(sum(x(k,1:nx)*covar))
             td2(k)=0.0d0
  			 do i=1,jj
			    td2(k)=td2(k)+del(i)
	         enddo
             do j1=1,jj
			    temp=0.0d0
			    do i=1,j1
  				temp=temp+del(i)
	            enddo
   				td2(k)=td2(k)+dexp(sum(x(k,1:nx)*covar)+day(j1))*temp
	         enddo
  			 td2(k)=td2(k)-dexp(day(j))*td1(k)
          enddo
          sigmad=sigma(2)/AA(j,j)
          dmu=0.0d0
          do jd=1,jj
             dmu=dmu+AA(j,jd)*day(jd)
          enddo
          dmu=-(dmu-AA(j,j)*day(j))/AA(j,j)
          call smpl_day_m(outd)
          day(j)=outd
          trace2(ig,j)=day(j)
       enddo
       do j=1,nx
          x1=sum((1.0-y)*x(1:nn,j))
          do k=1,nn
              x3(k)=x(k,j)
          enddo
          do k=1,nn
		     delage=0.0d0
		     do j1=1,jj
			    temp=0.0d0
  				do i1=1,j1
	  			temp=temp+del(i1)
		   		enddo
                delage=delage+dexp(day(j1))*temp
             enddo
             x2(k)=delage*dexp(sum(x(k,1:nx)*covar)-x(k,j)*covar(j))
             x4(k)=0.0d0
             do i1=1,jj
                x4(k)=x4(k)+del(i1)
             enddo
          enddo
          call smpl_x_m(outx)
          covar(j)=outx
          trace3(ig,j)=covar(j)
       enddo
       do k=1,nn
          deno=1.0d0
          do j=1,jj
             deno=deno+dexp(day(j)+sum(x(k,1:nx)*covar))
          enddo
          do j=1,jj
             q(k,j)=dexp(day(j)+sum(x(k,1:nx)*covar))/deno
          enddo
          q(k,jj+1)=1.0-sum(q(k,1:jj))
       enddo
c*************** updating sigma
ccc    updating sigma(1)    !! encounter variation
       scale=0.0d0
       do j=3,jj-1
          scale=scale+(enc(j)-2*enc(j-1)+enc(j-2))**2
       enddo
       scale=scale/2.0d0+b(1)
       shape=((jj-3)*1.0d0)/2.0d0+a(1)
       scale1 = 1.0
       outs = gamdev(shape,scale1,iseed)
       sigma(1)=scale/outs(1)

ccc    updating sigma()     !! age variation
       scale=0.0d0
       do j=3,jj
          scale=scale+(day(j)-2*day(j-1)+day(j-2))**2
       enddo
       scale=scale/2.0d0+b(2)
       shape=((jj-2)*1.0d0)/2.0d0+a(2)
       scale1 = 1.0
       outs = gamdev(shape,scale1,iseed)
       sigma(2)=scale/outs(1)
       trace4(ig,1:2)=sigma(1:2)
c*************** updating tha
ccc    updating tha(1)
       tm0=sum(mi*y)
       tm1=y
       call smpl_tha(outtha)
       tha(1)=outtha
ccc	   updating tha(2)
       call smpl_tha2(outtha)
  	 tha(2)=outtha
       tha1=tha
       trace5(ig,1:2)=tha(1:2)
!!!************* updating U_k,T_k,Y_k
         do k=1,nn
!*************** updating Y_k
            if(mi(k) .eq. 1) then
              y(k)=0
              ur(k)=jj-zl(k)+1
              ul(k)=jj-zr(k)+1
           py1=dexp(tha(1)+tha(2))/(1+dexp(tha(1)+tha(2)))
     +         *q(k,jj+1)*sum(del(ul(k):ur(k)))
           py2=0.0d0
  		      do iu=1,ur(k)
               idown=iu+zl(k)-1
                 iup=iu+zr(k)-1
               if (iup .gt. jj) then
                  iup=jj
                  if (idown .gt. jj) idown=jj
               endif
                 py2=py2+del(iu)*sum(q(k,idown:iup))
              enddo
              py2=dexp(tha(1))/(1+dexp(tha(1)))*py2
              pys=py1/(py1+py2)
              rn=ran1(iseed)
              if(rn.lt.pys) y(k)=1
            endif
            trace6(ig,k)=y(k)
!*************** updating U_k (mi=1,ul and ur are updated; mi=0, ul and ur are created at the initial part.)
              ur(k)=jj-zl(k)+1       ! after updating fate, updating ur(k)
              ul(k)=jj-zr(k)+1       ! after updating fate, updating ul(k)
              if (y(k).eq.0)  ul(k)=1
              if (ul(k) < 1)   ul(k)=1
              if (ur(k) < 1)   ur(k)=1
                iup=ur(k)
              idown=ul(k)
    			  p=0.0d0
               if (y(k) .eq. 1) then
                  temp=0.0d0
                  do j1=idown,iup
                     temp=temp+del(j1)
                     p(j1)=temp
                  enddo
                  do j1=idown,iup-1
                     p(j1)=p(j1)/temp
                  enddo
               else
                  temp1=0.0d0
                  do i1=idown,iup
                     iup1=i1+zr(k)-1
                     idown1=i1+zl(k)-1
                     if (iup1 .gt. jj) then
                        iup1=jj
                        if (idown1 .gt. jj) then
                           idown1=jj
                        endif
                     endif
                     temp2=0.0d0
                     do j1=idown1,iup1
                        temp2=temp2+q(k,j1)
                     enddo
                     temp1=temp1+del(i1)*temp2
                     p(i1)=temp1
                  enddo
                  do i1=idown,iup-1
                     p(i1)=p(i1)/temp1
                  enddo
               endif
               p(iup)=1.0d0
               rn = ran1(iseed)
               if (rn .lt. p(idown)) then
                  u(k)=idown
               else
                  do  j1=idown,iup-1
                      if ((rn .ge. p(j1)) .and. (rn .lt. p(j1+1))) then
                         u(k)=j1+1
                      endif
                  enddo
               endif
!*************** updating T_k
            p=0.0d0
            if (y(k) .eq. 0) then
               idown=u(k)+zl(k)-1
               iup=u(k)+zr(k)-1
               if (iup .gt. jj) then
                  iup=jj
                  if (idown .gt. jj) idown=jj
               endif
               temp=0.0d0
               do j1=idown,iup
                  temp=temp+q(k,j1)
                  p(j1)=temp
               enddo
               do j1=idown,iup-1
                  p(j1)=p(j1)/temp
               enddo
               p(iup)=1.0d0
               rn=ran1(iseed)
               if (rn .lt. p(idown)) then
                  t(k)=idown
               else
                  do j1=idown,iup-1
                   if((rn .ge. p(j1)).and.(rn .lt. p(j1+1))) then
                       t(k)=j1+1
                   endif
                  enddo
               endif
            else
		       t(k)=jj+1
            endif
         enddo
! --------- calculate survival rate
       do k=1,nn
  		  do j=1,jj
             sr1(k,j)=sum(q(k,(j+1):(jj+1)))/sum(q(k,j:(jj+1)))
  	         if ( j .eq. 1) then
		        csr(k,j) = sr1(k,j)
		     else
		        csr(k,j) = sr1(k,j)*csr(k,j-1)
		     endif
	      enddo
       enddo

   	   do j=1,jj
	      asr(j)=sum(sr1(1:nn,j))/nn
	      if (j .eq. 1) then
		     casr(j) = asr(j)
	      else
	         casr(j) = asr(j)*casr(j-1)
	      endif
       enddo
!---------- calculate log-likelihood for update data
       logliks=0.0d0     ! logliks is for loglik for single cycle
       do k = 1,nn
           post(k)=0
           if (y(k)==1) then
              logliks=logliks+dlog(q(k,jj+1)*sum(del(ul(k):ur(k))))
              post(k)=q(k,jj+1)*sum(del(ul(k):ur(k)))
           else
              part=0.0d0
	        if(ul(k)==ur(k)) then
 			  ii=ul(k)
              part=part+del(ii)*sum(q(k,(ii+zl(k)-1):(ii+zr(k)-1)))
              else
                  do ii=1,jj-zr(k)+1
                   part=part+del(ii)*sum(q(k,(ii+zl(k)-1):(ii+zr(k)-1)))
                  enddo
                  do ii=jj-zr(k)+2, ur(k)
                      part=part+del(ii)*sum(q(k,(ii+zl(k)-1):jj))
                  enddo
              endif
              logliks=logliks+dlog(part)
              post(k)=part

          endif
       enddo
! --------- summary
         if (ig .gt. n0) then
             sdel=sdel+del
             sdel2=sdel2+del**2
             sq=sq+q
             sq2=sq2+q**2
 	       ssr=ssr+sr1
	       ssr2=ssr2+sr1**2
             scsr=scsr+csr
             scsr2=scsr2+csr**2
	       sasr=sasr+asr
	       sasr2=sasr2+asr**2
	       scasr=scasr+casr
             scasr2=scasr2+casr**2
             slogliks=slogliks+logliks
	       stha=stha+tha
	       stha2=stha2+tha**2
	       scovar=scovar+covar
	       scovar2=scovar2+covar**2
             senc=senc+enc
             senc2=senc2+enc**2
             sday=sday+day
             sday2=sday2+day**2
             ssigma=ssigma+sigma
             ssigma2=ssigma2+sigma**2
             do k=1,nn
                 spost(k)=spost(k)+post(k)
                 logpost(k)=logpost(k)+dlog(post(k))
                 logpost2(k)=logpost2(k)+dlog(post(k))**2
              enddo
         endif
 16    continue

cc******* averages of the Gibbs samplers
      ddd=dble(ntotal-n0)
      del=sdel/ddd
      vdel=sqrt(sdel2/ddd - del**2)
      q=sq/ddd
      vq=sqrt(sq2/ddd - q**2)
      sr1=ssr/ddd
      vsr=sqrt(ssr2/ddd - sr1**2)
      csr=scsr/ddd
      vcsr=sqrt(scsr2/ddd - csr**2)
      asr=sasr/ddd
      vasr=sqrt(sasr2/ddd - asr**2)
      casr=scasr/ddd
      vcasr=sqrt(scasr2/ddd - casr**2)
      tha=stha/ddd
      vtha=sqrt(stha2/ddd - tha**2)
      covar=scovar/ddd
      vcovar=sqrt(scovar2/ddd - covar**2)
      enc=senc/ddd
      venc=sqrt(senc2/ddd-enc**2)
      day=sday/ddd
      vday=sqrt(sday2/ddd-day**2)
      sigma=ssigma/ddd
      vsigma=sqrt(ssigma2/ddd-sigma**2)
      do k=1,nn
          spost(k)=spost(k)/ddd
          lpd=lpd+dlog(spost(k))
          logpost(k)=logpost(k)/ddd
          vlogpost(k)=logpost2(k)/(ddd-1)-logpost(k)**2*ddd/(ddd-1)
          pwaic=pwaic+vlogpost(k)
      enddo
!---------- calculate deviance of the mean for DIC
      loglikm=0.0d0     ! loglikm is for loglik-mean
       do k = 1,nn
          if (y(k)==1) then
          loglikm=loglikm+dlog(q(k,jj+1)*sum(del(ul(k):ur(k))))
          else
            part=0.0d0
	        if(ul(k)==ur(k)) then
			  ii=ul(k)
              part=part+del(ii)*sum(q(k,(ii+zl(k)-1):(ii+zr(k)-1)))
   			else
              do ii=1,jj-zr(k)+1
              part=part+del(ii)*sum(q(k,(ii+zl(k)-1):(ii+zr(k)-1)))
              enddo
              do ii=jj-zr(k)+2, ur(k)
              part=part+del(ii)*sum(q(k,(ii+zl(k)-1):jj))
              enddo
              endif
            loglikm=loglikm+dlog(part)
          endif
       enddo
!---------------  DIC
         Dbar=-2*slogliks/ddd
         pd=Dbar+2*loglikm
         DIC=Dbar+pd
         WAIC=-2*(lpd-pwaic)






      end








      subroutine eval_enc_m(xs,hxx,hpxx) ! hxx=Ln(f(xs)); hpxx=(Ln(f(xs)))' ;
      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      integer nnmax
      parameter (nnmax=20000)
      real*8 te1(nnmax),te2(nnmax)
      real*8 emu, sigmae
      common /toencm/ emu,sigmae, ne
      common /tosub1m/ te1,te2
      common /inputm/ nn1
      hxx = ne*xs-(xs-emu)**2/2.0/sigmae
      do k=1,nn1
         hxx=hxx-dlog(te1(k)*dexp(xs)+te2(k))
      enddo
      hpxx = ne-(xs-emu)/sigmae
      do k=1,nn1
         hpxx=hpxx-te1(k)*dexp(xs)/(te1(k)*dexp(xs)+te2(k))
      enddo
      return
      end
ccc*************************************************************************************
      subroutine smpl_enc_m(oute)
      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      logical lb,ub
      parameter(ns=100,m=3,emax=600.0d0)
      parameter(lb=.TRUE.,ub=.TRUE.,xlb=-4.0d0,xub=2.0d0)
      real*8 rwv(6*ns+15),x(m),hx(m),hpx(m)
      integer*4 iwv(ns+7)
      external eval_enc_m
      x(1)= -3.0d0                      ! lower starting point
      x(3)= 1.0d0                       ! upper starting point
      x(2)=(x(1)+x(3))/2.0d0
      do j=1,m
         xs=x(j)
         call eval_enc_m(xs,hxout,hpxout)
         hx(j)=hxout
         hpx(j)=hpxout
      enddo
      call initial(ns,m,emax,x,hx,hpx,lb,xlb,ub,xub,ifault,iwv,rwv)
      call sample(iwv,rwv,eval_enc_m,beta,ifault)
      oute=beta
      return
      end
ccc***********************************************************************************
      subroutine eval_day_m(xs,hxx,hpxx) ! hxx=Ln(f(xs)); hpxx=(Ln(f(xs)))' ;
cccc  hxx=Ln(f(xs)); hpxx=(Ln(f(xs)))' ;  jnow: xi(now)
      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      integer nnmax
      parameter (nnmax=20000)
      real*8 td1(nnmax),td2(nnmax)
      real*8  dmu, sigmad
      common /todaym/ dmu,sigmad,nd
      common /tosub2m/ td1,td2
      common /inputm/ nn1

      hxx = nd*xs-(xs-dmu)**2/2.0/sigmad
      do k=1,nn1
         hxx=hxx-dlog(td1(k)*dexp(xs)+td2(k))
      enddo
      hpxx = nd-(xs-dmu)/sigmad
      do k=1,nn1
         hpxx=hpxx-td1(k)*dexp(xs)/(td1(k)*dexp(xs)+td2(k))
      enddo
      return
      end
ccc*************************************************************************************
      subroutine smpl_day_m(outd)
      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      logical lb,ub
      parameter(ns=100,m=3,emax=600.0d0)
      parameter(lb=.TRUE.,ub=.TRUE.,xlb=-4.0d0,xub=0.0d0)
      real*8 rwv(6*ns+15),x(m),hx(m),hpx(m)
      integer*4 iwv(ns+7)
      external eval_day_m
      x(1)= -3.0d0                      ! lower starting point
      x(3)= -1.0d0                       ! upper starting point
      x(2)=(x(1)+x(3))/2.0d0
      do j=1,m
         xs=x(j)
         call eval_day_m(xs,hxout,hpxout)
         hx(j)=hxout
         hpx(j)=hpxout
      enddo
      call initial(ns,m,emax,x,hx,hpx,lb,xlb,ub,xub,ifault,iwv,rwv)
      call sample(iwv,rwv,eval_day_m,beta,ifault)
      outd=beta
      return
      end
ccc***********************************************************************************
      subroutine eval_x_m(xs,hxx,hpxx) ! hxx=Ln(f(xs)); hpxx=(Ln(f(xs)))' ;
      implicit  integer*4 (i-n), real*8 (a-h,o-z)
      integer   nnmax
      parameter (nnmax=20000)
      real*8 x2(nnmax),x3(nnmax),x4(nnmax)
      real*8 varx,x1
      common /toxm/ varx,x1
      common /tosub3m/ x2,x3,x4
      common /inputm/ nn1
      hxx = x1*xs-xs**2/2.0/varx
      do k=1,nn1
         hxx=hxx-dlog(x2(k)*dexp(x3(k)*xs)+x4(k))
      enddo
      hpxx = x1-xs/varx
      do k=1,nn1
         hpxx=hpxx-x3(k)*x2(k)*dexp(x3(k)*xs)/(x2(k)*dexp(x3(k)*xs)
     +   +x4(k))
      enddo
      return
      end
ccc*************************************************************************************
      subroutine smpl_x_m(outx)
      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      logical lb,ub
      parameter(ns=100,m=3,emax=600.0d0)
      parameter(lb=.TRUE.,ub=.TRUE.,xlb=-2.0d0,xub=2.0d0)
      real*8 rwv(6*ns+15),x(m),hx(m),hpx(m)
      integer*4 iwv(ns+7)
      external eval_x_m
      x(1)= -1.0d0                      ! lower starting point
      x(3)= 1.0d0                       ! upper starting point
      x(2)=(x(1)+x(3))/2.0d0
      do j=1,m
         xs=x(j)
         call eval_x_m(xs,hxout,hpxout)
         hx(j)=hxout
         hpx(j)=hpxout
      enddo
      call initial(ns,m,emax,x,hx,hpx,lb,xlb,ub,xub,ifault,iwv,rwv)
      call sample(iwv,rwv,eval_x_m,beta,ifault)
      outx=beta
      return
      end

ccc***********************************************************************************
      subroutine eval_tha(xs,hxx,hpxx) ! hxx=Ln(f(xs)); hpxx=(Ln(f(xs)))' ;
      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      integer   nnmax
      parameter (nnmax=20000)
      integer ntha
      parameter (ntha=2)
      real*8 tha1(ntha), mutha1(ntha), vartha1(ntha)
      integer*4 mi1(nnmax), tm1(nnmax)
      common /tosub4m/ mi1,tm1
      common /tosub5m/ tha1,mutha1,vartha1
      common /tosub6m/ tm0
      common /inputm/ nn1
      hxx = tha1(2)*tm0-(xs-mutha1(1))**2/2.0/vartha1(1)
      do k=1,nn1
         hxx=hxx+mi1(k)*xs-dlog(1+dexp(xs+tha1(2)*tm1(k)))
      enddo
      hpxx =-(xs-mutha1(1))/vartha1(1)
      do k=1,nn1
      hpxx=hpxx+mi1(k)-dexp(xs+tha1(2)*tm1(k))/(1+dexp(xs+tha1(2)
     +   *tm1(k)))
      enddo
      return
      end
ccc*************************************************************************************
      subroutine smpl_tha(outtha)
      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      logical lb,ub
      parameter(ns=100,m=3,emax=600.0d0)
      parameter(lb=.TRUE.,ub=.TRUE.,xlb=-3.0d0,xub=-1.0d0)
      real*8 rwv(6*ns+15),x(m),hx(m),hpx(m)
      integer*4 iwv(ns+7)
      external eval_tha
      x(1)= -2.5d0                      ! lower starting point
      x(3)= -1.5d0                       ! upper starting point
      x(2)=(x(1)+x(3))/2.0d0
      do j=1,m
         xs=x(j)
         call eval_tha(xs,hxout,hpxout)
         hx(j)=hxout
         hpx(j)=hpxout
      enddo
      call initial(ns,m,emax,x,hx,hpx,lb,xlb,ub,xub,ifault,iwv,rwv)
      call sample(iwv,rwv,eval_tha,beta,ifault)
      outtha=beta
      return
      end
ccc***********************************************************************************
      subroutine eval_tha2(xs,hxx,hpxx) ! hxx=Ln(f(xs)); hpxx=(Ln(f(xs)))' ;
      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      integer   nnmax
      parameter (nnmax=20000)
      integer ntha
      parameter (ntha=2)
      real*8 tha1(ntha), mutha1(ntha), vartha1(ntha)
   	  integer*4 mi1(nnmax), tm1(nnmax)
      common /tosub4m/ mi1,tm1
      common /tosub5m/ tha1,mutha1,vartha1
      common /tosub6m/ tm0
      common /inputm/ nn1
      hxx = xs*tm0+tha1(1)*sum(mi1)-(xs-mutha1(2))**2/2.0/vartha1(2)
      do k=1,nn1
         hxx=hxx-dlog(1+dexp(tha1(1)+xs*tm1(k)))
      enddo
      hpxx = tm0 - (xs-mutha1(2))/vartha1(2)
      do k=1,nn1
         hpxx=hpxx-tm1(k)*dexp(tha1(1)+xs*tm1(k))/(1+dexp(tha1(1)+xs
     +   *tm1(k)))
      enddo
      return
      end
!!!*************************************************************************************
      subroutine smpl_tha2(outtha)
      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      logical lb,ub
      parameter(ns=100,m=3,emax=600.0d0)
      parameter(lb=.TRUE.,ub=.TRUE.,xlb=0.0d0,xub=2.0d0)
      real*8 rwv(6*ns+15),x(m),hx(m),hpx(m)
      integer*4 iwv(ns+7)
      external eval_tha2
      x(1)= 0.5d0                      ! lower starting point
      x(3)= 1.5d0                       ! upper starting point
      x(2)=(x(1)+x(3))/2.0d0
      do j=1,m
         xs=x(j)
         call eval_tha2(xs,hxout,hpxout)
         hx(j)=hxout
         hpx(j)=hpxout
      enddo
      call initial(ns,m,emax,x,hx,hpx,lb,xlb,ub,xub,ifault,iwv,rwv)
      call sample(iwv,rwv,eval_tha2,beta,ifault)
      outtha=beta
      return
      end







c***************************************************************
c subroutine to include noncensored observations.

      subroutine BEANSP(jj,nx,nn,ul,ur,zl,zr,x,y,a,b,sigma,day,enc,
     +   covar,n0,ntotal,q,del,DIC,Dbar,pd,trace1,trace2,trace3,trace4,
     +   venc,vday,vsigma,vcovar,vdel,vq,asr,vasr,casr,vcasr,sr,vsr)
      implicit   real*8 (a-h,o-z)
      implicit   integer*4 (i-n)
      integer n0, ntotal, nx, jj, nn1
      integer ir(1)
      integer y(nn),ul(nn),ur(nn),zl(nn),zr(nn)
      integer u(nn),t(nn),tl(nn),tr(nn)
      integer nn
      integer nnmax
      parameter (nnmax=20000)

      real*8  logliks, loglikm
      real*8  varx,x1,DIC,Dbar,pd,part
      real*8  x(nn,nx)
      real*8  outs(1),outu(1),rr(1)
      real*8  scale1
      real*8  sigma(2),a(2), b(2)
      real*8  day(jj),enc(jj-1),covar(nx)
      real*8  del(jj),q(nn,jj+1)
      real*8  p(jj)
      real*8  sdel(jj), sdel2(jj), vdel(jj)
      real*8  sq(nn,jj+1), sq2(nn,jj+1),vq(nn,jj+1)

      real*8  sday(jj),sday2(jj),vday(jj)
      real*8  senc(jj-1),senc2(jj-1),venc(jj-1)
      real*8  scovar(nx),scovar2(nx),vcovar(nx)
      real*8  sr(nn,jj), ssr(nn,jj),ssr2(nn,jj),vsr(nn,jj)
      real*8  csr(nn,jj), scsr(nn,jj),scsr2(nn,jj),vcsr(nn,jj)
      real*8  asr(jj), sasr(jj),sasr2(jj),vasr(jj)
      real*8  casr(jj),scasr(jj),scasr2(jj),vcasr(jj)

      real*8  ssigma(2), vsigma(2), ssigma2(2)
      real*8  AA(jj,jj), EE(jj-1,jj-1)
      real*8  x2(nnmax), x3(nnmax)
      real*8  te1(nnmax), te2(nnmax)
      real*8  td1(nnmax), td2(nnmax)
      real*8  probs(jj,jj),probm(jj,jj)
      real*8  trace1(ntotal,jj-1),trace2(ntotal,jj)
      real*8  trace3(ntotal,nx),trace4(ntotal,2)

      COMMON /torandom/IXR,IYR,IZR
      common /toenc/ emu,sigmae,ne
      common /today/ dmu,sigmad,nd
      common /tox/ varx,x1
      common /tosub1/ te1,te2
      common /tosub2/ td1,td2
      common /tosub3/ x2,x3
      common /input/ nn1

      external   RANDOM, eval_enc, smpl_enc, eval_day, smpl_day
      varx=10.0d0
      nn1=nn
      iseed=1976629
      IXR=3951
      IYR=1234
      IZR=6345

      AA=0.0D0                !!! IAR(2) precision matrix for day
      DO I=1,jj
         DO J=1,jj
         IF ((I==J).AND.((I==1).OR.(I==jj))) AA(I,J)=1.0D0
         IF ((I==J).AND.((I==2).OR.(I==jj-1))) AA(I,J)=5.0D0
         IF ((I==J).AND.((I>2).AND.(I<jj-1))) AA(I,J)=6.0D0
         IF (ABS(I-J)==1) AA(I,J)=-4.0D0
         IF (ABS(I-J)==2) AA(I,J)=1.0D0
         ENDDO
      ENDDO
      AA(1,2)=-2.0D0
      AA(2,1)=-2.0D0
      AA(jj-1,jj)=-2.0D0
      AA(jj,jj-1)=-2.0D0


      EE=0.0D0                !!! IAR(2) precision matrix for enc
      DO I=1,jj-1
         DO J=1,jj-1
         IF ((I==J).AND.((I==1).OR.(I==jj-1))) EE(I,J)=1.0D0
         IF ((I==J).AND.((I==2).OR.(I==jj-2))) EE(I,J)=5.0D0
         IF ((I==J).AND.((I>2).AND.(I<jj-2))) EE(I,J)=6.0D0
         IF (ABS(I-J)==1) EE(I,J)=-4.0D0
         IF (ABS(I-J)==2) EE(I,J)=1.0D0
         ENDDO
      ENDDO
      EE(1,2)=-2.0D0
      EE(2,1)=-2.0D0
      EE(jj-2,jj-1)=-2.0D0
      EE(jj-1,jj-2)=-2.0D0
      do j=1,nx
         xmean=sum(x(1:nn,j))/(nn*1.0d0)
         x(1:nn,j)=x(1:nn,j)-xmean
      enddo
c**************************************************************
c     the following accumulate in the Gibbs sampling
c     they have to be initialized for every new estimation
c     specify the prior hyperparameters
c***************************************************************
      sdel=0.0d0; sdel2=0.0d0; sq=0.0d0; sq2=0.0d0
      ssr=0.0d0; ssr2=0.0d0; ssigma=0.0d0; ssigma2=0.0d0
      sasr=0.0d0; sasr2=0.0d0; scasr=0.0d0; scasr2=0.0d0
      sday=0.0d0; sday2=0.0d0; senc=0.0d0; senc2=0.0d0
      scovar=0.0d0; scovar2=0.0d0
      slogliks=0.0d0;  sqb=0.0d0

c***********************  assign initial value

      q=0.04;    q(1:nn,(1+jj))=1.0-sum(q(1,1:jj))
      do k=1,nn
c              -------  update the upper bound of u_k
         if (ur(k) .gt. (jj-zl(k)+1)) then
            ur(k)=jj-zl(k)+1
            if (ul(k) .gt. (jj-zl(k)+1)) then
               ul(k)=jj-zl(k)+1
            endif
         endif
         if (y(k) .eq. 1) then
c              -------  update the lower bound of u_k
            if (ul(k) .lt. (jj-zr(k)+1)) then
                ul(k)=jj-zr(k)+1
            endif
         endif
c              -------  initial value for latent variable u_k
         if (ur(k) .eq. ul(k)) then
            u(k)=ul(k)
         else
            ir = myund(ur(k)-ul(k), iseed)
            u(k)=ul(k)+ir(1)
         endif
c              -------  initial value for latent variable t_k
         if (y(k) .eq. 1) then
            t(k)=jj+1
         else
            tl(k)=u(k)+zl(k)-1
            if ( tl(k) .gt. jj) then
               tl(k)=jj
            endif
            tr(k)=u(k)+zr(k)-1
            if ( tr(k) .gt. jj ) then
               tr(k)=jj
            endif
            if (tr(k) .eq. tl(k)) then
               t(k)=tl(k)
            else
                ir = myund(tr(k)-tl(k),iseed)
               t(k)=tl(k)+ir(1)
            endif
         endif	
!         print *,u(k),t(k)		 
      enddo
c**********************    big Gibbs sampling loop
      do 15 i=1,ntotal
c*****************  updating encounter effect, enc()
       do j=1,jj-1
          ne=0
          do k=1,nn
             if(u(k)==j) ne=ne+1
             te1(k)=sum(q(k,j:(jj+1)))
             te2(k)=sum(q(k,jj:(jj+1)))
             do je=1,jj-1
                te2(k)=te2(k)+sum(q(k,je:(jj+1)))*dexp(enc(je))
             enddo
             te2(k)=te2(k)-te1(k)*dexp(enc(j))
          enddo
          sigmae=sigma(1)/EE(j,j)
          emu=0.0d0
          do je=1,jj-1
             emu=emu+EE(j,je)*enc(je)
          enddo
          emu=-(emu-EE(j,j)*enc(j))/EE(j,j)
 ! 		  print *,j,emu,sigma(1)
          call smpl_enc(oute)
          enc(j)=oute
          trace1(i,j)=enc(j)
       enddo
       deno=1.0d0
       do j=1,jj-1
          deno=deno+dexp(enc(j))
       enddo
       do j=1,jj-1
          del(j)=dexp(enc(j))/deno
       enddo
       del(jj)=1.0-sum(del(1:(jj-1)))
c*************** updating q(j), failure rate
       do j=1,jj
          nd=0
          do k=1,nn
             if(t(k)==j) nd=nd+1
c             td1(k)=sum(del(1:j))*dexp(sum(x(k,1:nx)*covar))
        td1(k)=sum(del(1:j))*dexp(sum(x(k,1:nx)*covar))
             td2(k)=1.0+dexp(day(jj)+sum(x(k,1:nx)*covar))
             do jd=1,jj-1
         td2(k)=td2(k)+sum(del(1:jd))*dexp(day(jd)+sum(x(k,1:nx)*covar))
             enddo
             td2(k)=td2(k)-td1(k)*dexp(day(j))
          enddo

          sigmad=sigma(2)/AA(j,j)
          dmu=0.0d0
          do jd=1,jj
             dmu=dmu+AA(j,jd)*day(jd)
          enddo
          dmu=-(dmu-AA(j,j)*day(j))/AA(j,j)
          call smpl_day(outd)
          day(j)=outd
          trace2(i,j)=day(j)
       enddo
       delage=0.0d0
       do j=1,jj
          delage=delage+sum(del(1:j))*dexp(day(j))
       enddo

       do j=1,nx
          x1=sum((1.0-y)*x(1:nn,j))
          do k=1,nn
             x3(k)=x(k,j)
             x2(k)=delage*dexp(sum(x(k,1:nx)*covar)-x(k,j)*covar(j))
          enddo
          call smpl_x(outx)
          covar(j)=outx
          trace3(i,j)=covar(j)
       enddo

       do k=1,nn
          deno=1.0d0
          do j=1,jj
             deno=deno+dexp(day(j)+sum(x(k,1:nx)*covar))
          enddo
          do j=1,jj
             q(k,j)=dexp(day(j)+sum(x(k,1:nx)*covar))/deno
          enddo
          q(k,jj+1)=1.0-sum(q(k,1:jj))
       enddo
c*************** updating sigma   
ccc    updating sigma(1)    !! encounter variation

       scale=0.0d0
       do j=3,jj-1
          scale=scale+(enc(j)-2*enc(j-1)+enc(j-2))**2
       enddo
       scale=scale/2.0d0+b(1)
       shape=((jj-3)*1.0d0)/2.0d0+a(1)

       scale1 = 1.0
       outs = gamdev(shape,scale1,iseed)       
       sigma(1)=scale/outs(1)

ccc    updating sigma()     !! age variation

       scale=0.0d0
       do j=3,jj
          scale=scale+(day(j)-2*day(j-1)+day(j-2))**2
       enddo
       scale=scale/2.0d0+b(2)
       shape=((jj-2)*1.0d0)/2.0d0+a(2)

       scale1 = 1.0
       outs = gamdev(shape,scale1,iseed)
       sigma(2)=scale/outs(1)
       trace4(i,1:2)=sigma(1:2)

c*************** updating u_k

            do k=1,nn
               iup=ur(k)
               idown=ul(k)
               if (y(k) .eq. 1) then
                  temp=0.0d0
                  do j1=idown,iup
                     temp=temp+del(j1)
                     p(j1)=temp
                  enddo
                  do j1=idown,iup-1
                     p(j1)=p(j1)/temp
                  enddo
              else
                  temp1=0.0d0
                  do i1=idown,iup
                     iup1=i1+zr(k)-1
                     idown1=i1+zl(k)-1
                     if (iup1 .gt. jj) then
                        iup1=jj
                        if (idown1 .gt. jj) then
                           idown1=jj
                        endif
                     endif
                     temp2=0.0d0
                     do j1=idown1,iup1
                        temp2=temp2+q(k,j1)
                     enddo
                     temp1=temp1+del(i1)*temp2
                     p(i1)=temp1
                  enddo
                  do i1=idown,iup-1
                     p(i1)=p(i1)/temp1
                  enddo
               endif
               p(iup)=1.0d0
               rn = ran1(iseed)
               if (rn .lt. p(idown)) then
                  u(k)=idown

               else
                  do  j1=idown,iup-1
                      if ((rn .ge. p(j1)) .and. (rn .lt. p(j1+1))) then
                         u(k)=j1+1
                      endif
                  enddo
               endif
c*************** updating t_k
              if (y(k) .eq. 0) then
                  idown=u(k)+zl(k)-1
                  iup=u(k)+zr(k)-1
                  if (iup .gt. jj) then
                     iup=jj
                     if (idown .gt. jj) then
                        idown=jj
                     endif
                  endif
                  temp=0.0d0
                  do j1=idown,iup
                     temp=temp+q(k,j1)
                     p(j1)=temp
                  enddo
                  do j1=idown,iup-1
                     p(j1)=p(j1)/temp
                  enddo
                  p(iup)=1.0d0
                  rn = ran1(iseed)
                  if (rn .lt. p(idown)) then
                     t(k)=idown

                  else
                     do  j1=idown,iup-1
                        if((rn .ge. p(j1)).and.(rn .lt. p(j1+1))) then
                           t(k)=j1+1
                        endif
                     enddo
                  endif
               endif
            enddo

c --------- calculate survival rate
          do k=1,nn
             do j=1,jj
                sr(k,j)=sum(q(k,(j+1):(jj+1)))/sum(q(k,j:(jj+1)))
  	            if ( j .eq. 1) then
		               csr(k,j) = sr(k,j)
		           else
		               csr(k,j) = sr(k,j)*csr(k,j-1)
		           endif
             enddo
          enddo

	     do j=1,jj
	        asr(j)=sum(sr(1:nn,j))/nn
	        if (j .eq. 1) then
		         casr(j) = asr(j)
	        else
		         casr(j) = asr(j)*casr(j-1)
	        endif
	     enddo

c---------- calculate log-likelihood for DIC
              logliks=0.0d0     ! logliks is for loglik-singlecycle
              do k = 1,nn
                 if (y(k)==1) then
                    logliks=logliks+dlog(q(k,jj+1)
     &                      *sum(del(ul(k):ur(k))))
                 else
                   do ie=1,jj     ! ie i_encounter
                      do io=1,jj    ! io i_outcome
                         probs(ie,io)=del(ie)*q(k,io)
                      enddo
                   enddo

                   part=0.0d0
                   do ii=1, jj-zr(k)+1
                  part=part+sum(probs(ii,(ii+zl(k)-1):(ii+zr(k)-1)))
                   enddo
                   do ii=jj-zr(k)+2, ur(k)
                      part=part+sum(probs(ii,(ii+zl(k)-1):jj))
                   enddo

                   logliks=logliks+dlog(part)
                 endif
              enddo  
c --------- summary
         
         if (i .gt. n0) then

             sdel=sdel+del
             sdel2=sdel2+del**2

             sq=sq+q
             sq2=sq2+q**2

             ssr=ssr+sr
             ssr2=ssr2+sr**2

             scsr=scsr+csr
             scsr2=scsr2+csr**2

	           sasr=sasr+asr
	           sasr2=sasr2+asr**2

	           scasr=scasr+casr
	           scasr2=scasr2+casr**2

             ssigma=ssigma+sigma
             ssigma2=ssigma2+sigma**2

             senc=senc+enc
             senc2=senc2+enc**2
             sday=sday+day
             sday2=sday2+day**2
             scovar=scovar+covar
             scovar2=scovar2+covar**2

             sqb=sqb+qb
             slogliks=slogliks+logliks
          endif
c       --------------- end of the big Gibbs sampling loop
 15    continue

c******* averages of the Gibbs samplers
      ddd=dble(ntotal-n0)

      del=sdel/ddd
      vdel=sqrt(sdel2/ddd - del**2)

      q=sq/ddd
      vq=sqrt(sq2/ddd - q**2)

      sr=ssr/ddd
      vsr=sqrt(ssr2/ddd - sr**2)

      csr=scsr/ddd
      vcsr=sqrt(scsr2/ddd - csr**2)

      asr=sasr/ddd
      vasr=sqrt(sasr2/ddd - asr**2)

      casr=scasr/ddd
      vcasr=sqrt(scasr2/ddd - casr**2)

      sigma=ssigma/ddd
      vsigma=sqrt(ssigma2/ddd-sigma**2)

      day=sday/ddd
      vday=sqrt(sday2/ddd-day**2)

      covar=scovar/ddd
      vcovar=sqrt(scovar2/ddd-covar**2)

      enc=senc/ddd
      venc=sqrt(senc2/ddd-enc**2)
c---------- calculate deviance of the mean for DIC
              loglikm=0.0d0     ! logliks is for loglik-singlecycle
              do k = 1,nn

                 if (y(k)==1) then
                    loglikm=loglikm+dlog(q(k,jj+1)
     &                      *sum(del(ul(k):ur(k))))
                 else

                   do ie=1,jj     ! ie i_encounter
                      do io=1,jj    ! io i_outcome
                         probm(ie,io)=del(ie)*q(k,io)
                      enddo
                   enddo

                   part=0.0d0
                   do ii=1, jj-zr(k)+1
                  part=part+sum(probm(ii,(ii+zl(k)-1):(ii+zr(k)-1)))
                   enddo
                   do ii=jj-zr(k)+2, ur(k)
                      part=part+sum(probm(ii,(ii+zl(k)-1):jj))
                   enddo

                   loglikm=loglikm+dlog(part)

                 endif
              enddo

c---------------  DIC

         Dbar=-2*slogliks/ddd
         pd=Dbar+2*loglikm

         DIC=Dbar+pd
c       write(*,*) 'Completed'

      end
ccc***********************************************************************************
      subroutine eval_enc(xs,hxx,hpxx) 
cccc  hxx=Ln(f(xs));hpxx=(Ln(f(xs)))';jnow:xi(now)
      implicit  integer*4 (i-n), real*8 (a-h,o-z)
      integer nnmax
      parameter (nnmax=20000)
      real*8 te1(nnmax),te2(nnmax)
      real*8 emu, sigmae
      common /toenc/ emu,sigmae, ne
      common /tosub1/ te1,te2
      common /input/ nn1
      hxx = ne*xs-(xs-emu)**2/2.0/sigmae
      do k=1,nn1
         hxx=hxx-dlog(te1(k)*dexp(xs)+te2(k))
      enddo

      hpxx = ne-(xs-emu)/sigmae
      do k=1,nn1
         hpxx=hpxx-te1(k)*dexp(xs)/(te1(k)*dexp(xs)+te2(k))
      enddo
      return
      end
ccc*************************************************************************************
      subroutine smpl_enc(oute)
      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      logical lb,ub
      parameter(ns=100,m=3,emax=600.0d0)
      parameter(lb=.TRUE.,ub=.TRUE.,xlb=-4.0d0,xub=2.0d0)
      real*8 rwv(6*ns+15),x(m),hx(m),hpx(m)
      integer*4 iwv(ns+7)
      external eval_enc

      x(1)= -3.0d0                      ! lower starting point
      x(3)= 1.0d0                       ! upper starting point
      x(2)=(x(1)+x(3))/2.0d0
      do j=1,m
         xs=x(j)
         call eval_enc(xs,hxout,hpxout)
         hx(j)=hxout
         hpx(j)=hpxout
      enddo
      call initial(ns,m,emax,x,hx,hpx,lb,xlb,ub,xub,ifault,iwv,rwv)

      call sample(iwv,rwv,eval_enc,beta,ifault)

      oute=beta

      return
      end
ccc***********************************************************************************

      subroutine eval_day(xs,hxx,hpxx) 
cccc  hxx=Ln(f(xs)); hpxx=(Ln(f(xs)))' ;  jnow: xi(now)
      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      integer nnmax
      parameter (nnmax=20000)
      real*8 td1(nnmax),td2(nnmax)
      real*8  dmu, sigmad
      common /today/ dmu,sigmad,nd
      common /tosub2/ td1,td2
      common /input/ nn1

      hxx = nd*xs-(xs-dmu)**2/2.0/sigmad

      do k=1,nn1
         hxx=hxx-dlog(td1(k)*dexp(xs)+td2(k))
      enddo

      hpxx = nd-(xs-dmu)/sigmad

      do k=1,nn1
         hpxx=hpxx-td1(k)*dexp(xs)/(td1(k)*dexp(xs)+td2(k))
      enddo

      return
      end


ccc*************************************************************************************

      subroutine smpl_day(outd)

      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      logical lb,ub
      parameter(ns=100,m=3,emax=600.0d0)
      parameter(lb=.TRUE.,ub=.TRUE.,xlb=-4.0d0,xub=0.0d0)
      real*8 rwv(6*ns+15),x(m),hx(m),hpx(m)
      integer*4 iwv(ns+7)
      external eval_day


      x(1)= -3.0d0                      ! lower starting point
      x(3)= -1.0d0                       ! upper starting point
      x(2)=(x(1)+x(3))/2.0d0
      do j=1,m
         xs=x(j)
         call eval_day(xs,hxout,hpxout)
         hx(j)=hxout
         hpx(j)=hpxout
      enddo
      call initial(ns,m,emax,x,hx,hpx,lb,xlb,ub,xub,ifault,iwv,rwv)

      call sample(iwv,rwv,eval_day,beta,ifault)

      outd=beta

      return
      end
ccc***********************************************************************************

      subroutine eval_x(xs,hxx,hpxx) 
cccc  hxx=Ln(f(xs)); hpxx=(Ln(f(xs)))' ;  jnow: xi(now)
      implicit  integer*4 (i-n), real*8 (a-h,o-z)
      integer   nnmax
      parameter (nnmax=20000)
      real*8 x2(nnmax),x3(nnmax)
      real*8 varx,x1
      common /tox/ varx,x1
      common /tosub3/ x2,x3
      common /input/ nn1

      hxx = x1*xs-xs**2/2.0/varx

      do k=1,nn1
         hxx=hxx-dlog(x2(k)*dexp(x3(k)*xs)+1.0)
      enddo

      hpxx = x1-xs/varx

      do k=1,nn1
         hpxx=hpxx-x3(k)*x2(k)*dexp(x3(k)*xs)/(x2(k)*dexp(x3(k)*xs)+1.0)
      enddo

      return
      end


ccc*************************************************************************************

      subroutine smpl_x(outx)

      implicit      integer*4 (i-n), real*8 (a-h,o-z)
      logical lb,ub
      parameter(ns=100,m=3,emax=600.0d0)
      parameter(lb=.TRUE.,ub=.TRUE.,xlb=-2.0d0,xub=2.0d0)
      real*8 rwv(6*ns+15),x(m),hx(m),hpx(m)
      integer*4 iwv(ns+7)
      external eval_x


      x(1)= -1.0d0                      ! lower starting point
      x(3)= 1.0d0                       ! upper starting point
      x(2)=(x(1)+x(3))/2.0d0
      do j=1,m
         xs=x(j)
         call eval_x(xs,hxout,hpxout)
         hx(j)=hxout
         hpx(j)=hpxout
      enddo
      call initial(ns,m,emax,x,hx,hpx,lb,xlb,ub,xub,ifault,iwv,rwv)

      call sample(iwv,rwv,eval_x,beta,ifault)

      outx=beta

      return
      end

c--------------End of written program-----------------------------
      DOUBLE PRECISION FUNCTION EXPON(X,EMAX)
C
C Performs an exponential without underflow
C
      IMPLICIT REAL*8 (A-H,O-Z)
      ZEXP(X)=DEXP(X)
      IF (X .LT. -EMAX) THEN
        EXPON=0.0D0
      ELSE
        EXPON=ZEXP(X)
      ENDIF
      RETURN
      END
C
C***********************************************************************

C
      DOUBLE PRECISION FUNCTION RANDOM(L)
C
C Algorithm AS 183 APPL STATIST. (1982) Vol 31 n 2 Modified according to
C AS R58.  Returns a pseudo-random number rectangularly distributed between
C 0 and 1. IXR, IYR, AND IZR should be set to integer values between 1 and
C 30000 before first entry. Integer arithmetic up to 30323 required.
C
       IMPLICIT NONE
       INTEGER*4 L
       COMMON /TORANDOM/IXR,IYR,IZR
       INTEGER*4        IXR,IYR,IZR

       IXR = 171 * MOD(IXR,177) -2 *(IXR/177)
       IYR = 172 * MOD(IYR,176) -35*(IYR/176)
       IZR = 170 * MOD(IZR,178) -63*(IZR/178)
       IF (IXR.LT.0) IXR=IXR+30269
       IF (IYR.LT.0) IYR=IYR+30307
       IF (IZR.LT.0) IZR=IZR+30323
       RANDOM= AMOD(FLOAT(IXR)/30269.0 +FLOAT(IYR) /30307.0 +
     +              FLOAT(IZR)/30323.0,1.0)
       IF (RANDOM.GT.0.0)  RETURN
       RANDOM=DMOD(DBLE(FLOAT(IXR))/30269.0D0 +
     +  DBLE(FLOAT(IYR))/30307.0D0+DBLE(FLOAT(IZR))/30323.0D0,1.0D0)
       IF (RANDOM.GE.1) RANDOM=0.999999
       RETURN
       END
C
C***********************************************************************

C
      SUBROUTINE INITIAL(ns,m,emax,x,hx,hpx,lb,xlb,ub,xub,ifault,
     +                                                          iwv,rwv)

C
C This subroutine takes as input the number of starting values m
C and the starting values x(i),hx(i),hpx(i)  i=1,m
C As output we have pointer ipt along with ilow and ihigh and the lower
C and upper hulls defined  by z,hz,scum,cu,hulb,huub stored in working
C vectors iwv and rwv
C Ifault detects wrong starting points or non-concavity
C
      INTEGER ns,nn,m,ilow,ihigh,ifault,i,iwv(*)
      INTEGER iipt,iz,ihuz,iscum,ix,ihx,ihpx
      LOGICAL ub,lb,horiz
      DOUBLE PRECISION xlb,xub,emax,x(*),hx(*),hpx(*),rwv(*),expon
      DOUBLE PRECISION hulb,huub,eps,cu,alcu,huzmax,zlog,zmax
      DOUBLE PRECISION a,b
      zlog(a)=dlog(a)
      zmax(a,b)=dmax1(a,b)
C
C DESCRIPTION OF PARAMETERS and place of storage
C
C lb   iwv(5) : boolean indicating if there is a lower bound to the domain
C ub   iwv(6) : boolean indicating if there is an upper bound
C xlb  rwv(8) : value of the lower bound
C xub  rwv(9) : value of the upper bound
C emax rwv(3) : large value for which it is possible to compute an
C               exponential, eps=exp(-emax) is taken as a small value used to
C               test for numerical unstability
C m    iwv(4) : number of starting points
C ns   iwv(3) : maximum number of points defining the hulls
C x    rwv(ix+1)  : vector containing the abscissae of the starting points
C hx   rwv(ihx+1) : vector containing the ordinates
C hpx  rwv(ihpx+1): vector containing the derivatives
C ifault      : diagnostic
C iwv,rwv     : integer and real working vectors
C
      eps=expon(-emax,emax)
      ifault=0
      ilow=1
      ihigh=1
      nn=ns+1
C at least one starting point
      IF (m.LT.1) ifault=1
      huzmax=hx(1)
      IF (.NOT.ub) xub=0.0
      IF (.NOT.lb) xlb=0.0
      hulb=(xlb-x(1))*hpx(1)+hx(1)
      huub=(xub-x(1))*hpx(1)+hx(1)
C if bounded on both sides
      IF ((ub).AND.(lb)) THEN
        huzmax=zmax(huub,hulb)
        horiz=(abs(hpx(1)).LT.eps)
        IF (horiz) THEN
          cu=expon((huub+hulb)*0.5-huzmax,emax)*(xub-xlb)
        ELSE
          cu=expon(huub-huzmax,emax)*(1-expon(hulb-huub,emax))/hpx(1)
        ENDIF
      ELSE IF ((ub).AND.(.NOT.lb))THEN
C if bounded on the right and unbounded on the left
        huzmax=huub
        cu=1.0/hpx(1)
      ELSE IF ((.NOT.ub).AND.(lb))THEN
C if bounded on the left and unbounded on the right
        huzmax=hulb
        cu=-1.0/hpx(1)
C if unbounded at least 2 starting points
      ELSE
        cu=0.0
        IF (m.LT.2) ifault=1
      ENDIF
      IF (cu.GT.0.0) alcu=zlog(cu)
C set pointers
      iipt=6
      iz=9
      ihuz=nn+iz
      iscum=nn+ihuz
      ix=nn+iscum
      ihx=nn+ix
      ihpx=nn+ihx
C store values in working vectors
      iwv(1)=ilow
      iwv(2)=ihigh
      iwv(3)=ns
      iwv(4)=1
      IF (lb) THEN
        iwv(5)=1

      ELSE
        iwv(5)=0
      ENDIF
      IF (ub) THEN
        iwv(6)=1
      ELSE
        iwv(6)=0
      ENDIF
      IF (ns.LT.m) ifault=2
      iwv(iipt+1)=0
      rwv(1)=hulb
      rwv(2)=huub
      rwv(3)=emax
      rwv(4)=eps
      rwv(5)=cu
      rwv(6)=alcu
      rwv(7)=huzmax
      rwv(8)=xlb
      rwv(9)=xub
      rwv(iscum+1)=1.0
      DO 9 i=1,m
        rwv(ix+i)=x(i)
        rwv(ihx+i)=hx(i)
        rwv(ihpx+i)=hpx(i)
  9   CONTINUE
C create lower and upper hulls
      i=1
 10   IF (i.LT.m)THEN
         CALL update(iwv(4),iwv(1),iwv(2),iwv(iipt+1),rwv(iscum+1),
     +           rwv(5),rwv(ix+1),rwv(ihx+1),rwv(ihpx+1),rwv(iz+1),
     +           rwv(ihuz+1),rwv(7),rwv(3),lb,rwv(8),rwv(1),ub,rwv(9),
     +           rwv(2),ifault,rwv(4),rwv(6))
         i=iwv(4)
         IF (ifault.NE.0) RETURN
      GOTO 10
      ENDIF
C test for wrong starting points
      IF ((.NOT.lb).AND.(hpx(iwv(1)).LT.eps)) ifault=3
      IF ((.NOT.ub).AND.(hpx(iwv(2)).GT.-eps)) ifault=4
      RETURN
      END

C
C***********************************************************************

C
      SUBROUTINE SAMPLE(iwv,rwv,eval,beta,ifault)

      INTERFACE
         SUBROUTINE eval(A,B,C)
         REAL*8 :: A,B,C
         END SUBROUTINE eval
      END INTERFACE
C
      DOUBLE PRECISION beta,rwv(*)
      INTEGER iipt,iz,ns,nn,ihuz,iscum,ix,ihx,ihpx,ifault,iwv(*)
      LOGICAL ub,lb
C     EXTERNAL eval
C
C     set pointers
c      print*,'in sample, stage 0'
      iipt=6
      iz=9
      ns=iwv(3)
      nn=ns+1
      ihuz=nn+iz
      iscum=nn+ihuz
      ix=nn+iscum
      ihx=nn+ix
      ihpx=nn+ihx
      lb=.FALSE.
      ub=.FALSE.
      IF (iwv(5).EQ.1) lb=.TRUE.
      IF (iwv(6).EQ.1) ub=.TRUE.
C     call sampling subroutine
      call spl1(ns,iwv(4),iwv(1),iwv(2),iwv(iipt+1),rwv(iscum+1),
     +  rwv(5),rwv(ix+1),rwv(ihx+1),rwv(ihpx+1),rwv(iz+1),rwv(ihuz+1),
     +  rwv(7),lb,rwv(8),rwv(1),ub,rwv(9),rwv(2),eval,beta,ifault,
     +  rwv(3),rwv(4),rwv(6))
c      print*,'in sample, stage 1,beta=',beta
      RETURN
      END

C
C***********************************************************************

C
      SUBROUTINE SPL1(ns,n,ilow,ihigh,ipt,scum,cu,x,hx,hpx,z,huz,
     +  huzmax,lb,xlb,hulb,ub,xub,huub,eval,beta,ifault,emax,eps,alcu)

      INTERFACE
         SUBROUTINE eval(A,B,C)
         REAL*8 :: A,B,C
         END SUBROUTINE eval
      END INTERFACE

C
C this subroutine performs the adaptive rejection sampling, it calls
C subroutine splhull to sample from the upper hull ,if the sampling
C involves a function evaluation it calls the updating subroutine
C ifault is a diagnostic of any problem: non concavity, 0 random number
C or numerical imprecision
C
      INTEGER ns,n,ilow,ihigh,ifault,i,j,ipt(*),n1,l
      LOGICAL ub,lb,sampld
      DOUBLE PRECISION z(*),huz(*),x(*),hx(*),hpx(*),scum(*)
      DOUBLE PRECISION random,xlb,xub,emax,u11,u22,beta,alu1,fx
      DOUBLE PRECISION hulb,huub,eps,cu,alcu,huzmax,alhl,alhu
      DOUBLE PRECISION zlog,a
C      EXTERNAL eval
      zlog(a)=dlog(a)
C
      sampld=.FALSE.
 10   IF (.NOT.sampld) THEN
        u22=random(L)
C test for zero random number
        IF (u22.EQ.0.0) THEN
          ifault=6
          RETURN
        ENDIF
      CALL splhull(u22,ipt,ilow,lb,xlb,hulb,huzmax,alcu,x,hx,hpx,
     +                              z,huz,scum,eps,emax,beta,i,j)
C sample u11 to compute rejection
        u11=random(l)
        IF (u11.EQ.0.0) ifault=6
        alu1=zlog(u11)
C compute alhu: upper hull at point u11
        alhu=hpx(i)*(beta-x(i))+hx(i)-huzmax
        IF ((beta.GT.x(ilow)).AND.(beta.LT.x(ihigh))) THEN
C compute alhl: value of the lower hull at point u11
          IF (beta.GT.x(i)) THEN
            j=i
            i=ipt(i)
          ENDIF
          alhl=hx(i)+(beta-x(i))*(hx(i)-hx(j))/(x(i)-x(j))-huzmax

C squeezing test
          IF ((alhl-alhu).GT.alu1) THEN
             sampld=.TRUE.
          ENDIF
        ENDIF
C if not sampled evaluate the function ,do the rejection test and update

        IF (.NOT.sampld) THEN
          n1=n+1
          x(n1)=beta
          CALL eval(x(n1),hx(n1),hpx(n1))
          fx=hx(n1)-huzmax
          IF (alu1.lt.(fx-alhu)) sampld=.TRUE.
C update while the number of points defining the hulls is lower than ns
          IF (n.LT.ns)  THEN
            CALL update(n,ilow,ihigh,ipt,scum,cu,x,hx,hpx,z,huz,huzmax,
     +                 emax,lb,xlb,hulb,ub,xub,huub,ifault,eps,alcu)
          ENDIF
          IF (ifault.NE.0) RETURN
        ENDIF
        GOTO 10
      ENDIF
      RETURN
      END

C
C***********************************************************************

C
      SUBROUTINE SPLHULL(u22,ipt,ilow,lb,xlb,hulb,huzmax,alcu,x,hx,hpx,
     +                              z,huz,scum,eps,emax,beta,i,j)
C
C this subroutine samples beta from the normalised upper hull
C
      DOUBLE PRECISION z(*),huz(*),x(*),hx(*),hpx(*),scum(*)
      DOUBLE PRECISION expon,xlb,emax,eh,beta,u22,a,logdu,logtg,sign
      DOUBLE PRECISION hulb,eps,alcu,huzmax,zlog
      INTEGER ilow,i,j,ipt(*)
      LOGICAL lb,horiz
      zlog(a)=dlog(a)
C
        i=ilow
C
C find from which exponential piece you sample
 20   IF (u22.gt.scum(i))THEN
        j=i
        i=ipt(i)
        GOTO 20
      ENDIF
      IF (i.eq.ilow) THEN
C sample below z(ilow),depending on the existence of a lower bound
        IF (lb) THEN
          eh=hulb-huzmax-alcu
          horiz=(abs(hpx(ilow)).LT.eps)
          IF (horiz) THEN
            beta=xlb+u22*expon(-eh,emax)
          ELSE
           sign=abs(hpx(i))/hpx(i)
           logtg=zlog(abs(hpx(i)))
           logdu=zlog(u22)
           eh=logdu+logtg-eh
           IF (eh.LT.emax) THEN
            beta=xlb+zlog(1.0+sign*expon(eh,emax))/hpx(i)
           ELSE
            beta=xlb+eh/hpx(i)
           ENDIF
          ENDIF
        ELSE

C hpx(i) must be positive , x(ilow) is left of the mode
          beta=(zlog(hpx(i)*u22)+alcu-hx(i)+x(i)*hpx(i)+huzmax)/hpx(i)
        ENDIF
      ELSE
C sample above(j)
        eh=huz(j)-huzmax-alcu
        horiz=(abs(hpx(i)).LT.eps)
        IF (horiz) THEN
         beta=z(j)+(u22-scum(j))*expon(-eh,emax)
        ELSE
         sign=abs(hpx(i))/hpx(i)
         logtg=zlog(abs(hpx(i)))
         logdu=zlog(u22-scum(j))
         eh=logdu+logtg-eh
         IF (eh.LT.emax) THEN
          beta=z(j)+(zlog(1.0+sign*expon(eh,emax)))/hpx(i)
         ELSE
          beta=z(j)+eh/hpx(i)
         ENDIF
        ENDIF
      ENDIF
      RETURN
      END

C
C***********************************************************************

C
      SUBROUTINE INTERSECTION(x1,y1,yp1,x2,y2,yp2,z1,hz1,eps,ifault)
C
C computes the intersection (z1,hz1) between 2 tangents defined by
C   x1,y1,yp1 and x2,y2,yp2
C
      DOUBLE PRECISION x1,y1,yp1,x2,y2,yp2,z1,hz1,eps,y12,y21,dh
      INTEGER ifault
C
C first test for non-concavity
      y12=y1+yp1*(x2-x1)
      y21=y2+yp2*(x1-x2)
      IF ((y21.lt.y1).OR.(y12.LT.y2)) THEN
         ifault=5
         RETURN
      ENDIF
      dh=yp2-yp1
C IF the lines are nearly parallel,the intersection is taken at the midpoint
      IF (abs(dh).LE.eps)THEN
        z1=0.5*(x1+x2)
        hz1=0.5*(y1+y2)
C Else compute from the left or the right for greater numerical precision
      ELSE IF (abs(yp1).LT.abs(yp2)) THEN
        z1=x2+(y1-y2+yp1*(x2-x1))/dh
        hz1=yp1*(z1-x1)+y1
      ELSE
        z1=x1+(y1-y2+yp2*(x2-x1))/dh
        hz1=yp2*(z1-x2)+y2
      ENDIF
C Test for misbehaviour due to numerical imprecision
      IF ((z1.LT.x1).OR.(z1.GT.x2)) ifault=7
      RETURN
      END

C
C***********************************************************************

C
      SUBROUTINE UPDATE(n,ilow,ihigh,ipt,scum,cu,x,hx,hpx,z,huz,
     +          huzmax,emax,lb,xlb,hulb,ub,xub,huub,ifault,eps,alcu)

      INTERFACE
         SUBROUTINE eval(A,B,C)
         REAL*8 :: A,B,C
         END SUBROUTINE eval
      END INTERFACE
C
C This subroutine increments n and updates all the parameters which
C define the lower and the upper hull
C
      INTEGER n,ilow,ihigh,ifault,i,j,ipt(*)
      LOGICAL ub,lb,horiz
      DOUBLE PRECISION xlb,xub,x(*),hx(*),hpx(*),z(*),huz(*),scum(*)
      DOUBLE PRECISION hulb,huub,eps,cu,alcu,huzmax,emax,expon,dh,u
      DOUBLE PRECISION a,b,zlog,zmax
C      EXTERNAL eval
      zlog(a)=dlog(a)
      zmax(a,b)=dmax1(a,b)
C
C DESCRIPTION OF PARAMETERS and place of storage
C ilow iwv(1)    : index of the smallest x(i)
C ihigh iwv(2)   : index of the largest x(i)
C n    iwv(4)    : number of points defining the hulls
C ipt  iwv(iipt) : pointer array:  ipt(i) is the index of the x(.)
C                  immediately larger than x(i)
C hulb rwv(1)    : value of the upper hull at xlb
C huub rwv(2)    : value of the upper hull at xub
C cu   rwv(5)    : integral of the exponentiated upper hull divided
C                  by exp(huzmax)
C alcu rwv(6)    : logarithm of cu
C huzmax rwv(7)  : maximum of huz(i); i=1,n
C z    rwv(iz+1) : z(i) is the abscissa of the intersection between
C                  the tangents at x(i) and x(ipt(i))
C huz  rwv(ihuz+1): huz(i) is the ordinate of the intersection defined above
C scum rwv(iscum): scum(i) is the cumulative probability of the normalised
C                  exponential of the upper hull calculated at z(i)
C eps  rwv(4)    : =exp(-emax) a very small number
C
       n=n+1
C update z,huz and ipt
      IF (x(n).LT.x(ilow)) THEN

C insert x(n) below x(ilow)
C test for non-concavity
        IF (hpx(ilow).GT.hpx(n)) ifault=5
        ipt(n)=ilow
        CALL intersection(x(n),hx(n),hpx(n),x(ilow),hx(ilow),hpx(ilow),
     +         z(n),huz(n),eps,ifault)
        IF (ifault.NE.0) RETURN
        IF (lb) hulb=hpx(n)*(xlb-x(n))+hx(n)
        ilow=n
      ELSE
        i=ilow
        j=i
C find where to insert x(n)
  10    IF ((x(n).GE.x(i)).AND.(ipt(i).ne.0))THEN
          j=i
          i=ipt(i)
          GOTO 10
         ENDIF
        IF (x(n).GE.x(i)) THEN
C insert above x(ihigh)
C test for non-concavity
          IF (hpx(i).LT.hpx(n)) ifault=5
          ihigh=n
          ipt(i)=n
          ipt(n)=0
          CALL intersection(x(i),hx(i),hpx(i),x(n),hx(n),hpx(n),
     +         z(i),huz(i),eps,ifault)
          IF (ifault.NE.0) RETURN
          huub=hpx(n)*(xub-x(n))+hx(n)
          z(n)=0.0
          huz(n)=0.0
        ELSE
C insert x(n) between x(j) and x(i)
C test for non-concavity
          IF ((hpx(j).LT.hpx(n)).OR.(hpx(i).GT.hpx(n))) ifault=5
          ipt(j)=n
          ipt(n)=i
C insert z(j) between x(j) and x(n)
          CALL intersection(x(j),hx(j),hpx(j),x(n),hx(n),hpx(n),
     +         z(j),huz(j),eps,ifault)
          IF (ifault.NE.0) RETURN

C insert z(n) between x(n) and x(i)
          CALL intersection(x(n),hx(n),hpx(n),x(i),hx(i),hpx(i),
     +         z(n),huz(n),eps,ifault)
          IF (ifault.NE.0) RETURN
        ENDIF
      ENDIF
C update huzmax
      j=ilow
      i=ipt(j)
      huzmax=huz(j)
 20   IF ((huz(j).LT.huz(i)).AND.(ipt(i).NE.0))THEN
        j=i
        i=ipt(i)
        huzmax=zmax(huzmax,huz(j))
        GOTO 20
      ENDIF
      IF (lb) huzmax=zmax(huzmax,hulb)
      IF (ub) huzmax=zmax(huzmax,huub)
C update cu
C scum receives area below exponentiated upper hull left of z(i)
      i=ilow
      horiz=(abs(hpx(ilow)).LT.eps)
      IF ((.NOT.lb).AND.(.NOT.horiz)) THEN
        cu=expon(huz(i)-huzmax,emax)/hpx(i)
      ELSE IF (lb.AND.horiz) THEN
        cu=(z(ilow)-xlb)*expon(hulb-huzmax,emax)
      ELSE IF (lb.AND.(.NOT.horiz)) THEN
        dh=hulb-huz(i)
        IF (dh.GT.emax) THEN
          cu=-expon(hulb-huzmax,emax)/hpx(i)
        ELSE
          cu=expon(huz(i)-huzmax,emax)*(1-expon(dh,emax))/hpx(i)
        ENDIF
      ELSE
        cu=0
      ENDIF
      scum(i)=cu
      j=i
      i=ipt(i)
 30   IF (ipt(i).NE.0)THEN
        dh=huz(j)-huz(i)
        horiz=(abs(hpx(i)).LT.eps)
        IF (horiz) THEN
          cu=cu+(z(i)-z(j))*expon((huz(i)+huz(j))*0.5-huzmax,emax)
        ELSE

          IF (dh.LT.emax) THEN
            cu=cu+expon(huz(i)-huzmax,emax)*(1-expon(dh,emax))/hpx(i)
          ELSE
            cu=cu-expon(huz(j)-huzmax,emax)/hpx(i)
          ENDIF
        ENDIF
        j=i
        i=ipt(i)
        scum(j)=cu
        GOTO 30
      ENDIF
      horiz=(abs(hpx(i)).LT.eps)
C if the derivative is very small the tangent is nearly horizontal
      IF (.NOT.(ub.OR.horiz)) THEN
        cu=cu-expon(huz(j)-huzmax,emax)/hpx(i)
      ELSE IF (ub.AND.horiz) THEN
        cu=cu+(xub-x(i))*expon((huub+hx(i))*0.5-huzmax,emax)
      ELSE IF (ub.AND.(.NOT.horiz)) THEN
        dh=huz(j)-huub
        IF (dh.GT.emax) THEN
          cu=cu-expon(huz(j)-huzmax,emax)/hpx(i)
        ELSE
          cu=cu+expon(huub-huzmax,emax)*(1-expon(dh,emax))/hpx(i)
        ENDIF
      ENDIF
      scum(i)=cu
      IF (cu.GT.0) alcu=zlog(cu)
C normalize scum to obtain a cumulative probability while excluding
C unnecessary points
      i=ilow
      u=(cu-scum(i))/cu
      IF ((u.EQ.1.0).AND.(hpx(ipt(i)).GT.0)) THEN
        ilow=ipt(i)
        scum(i)=0.0
      ELSE
        scum(i)=1.0-u
      ENDIF
      j=i
      i=ipt(i)
 40   IF (ipt(i).NE.0) THEN
        j=i
        i=ipt(i)
        u=(cu-scum(j))/cu
        IF ((u.EQ.1.0).AND.(hpx(i).GT.0)) THEN
          ilow=i

        ELSE
          scum(j)=1.0-u
        ENDIF
        GOTO 40
      ENDIF
      scum(i)=1.0
      IF (ub) huub=hpx(ihigh)*(xub-x(ihigh))+hx(ihigh)
      IF (lb) hulb=hpx(ilow)*(xlb-x(ilow))+hx(ilow)
      RETURN
      END



c-----------------------------------------------------
c  This function gascum returns the cumulative distribution of 
c  Normal(0,1) at point x by the rational fraction approximation.
c  Reference: Statistical Computing, pp.92.
c  By William J. Kennedy, Jr., and James E. Gentle. 1980.

      real*8 function gascum(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter(mm=1000)
      real*8 w(0:mm)
      if(x.le.-1.d1)gascum=0.0d0
      if(x.ge.1.d1)gascum=1.0d0

      if(-1.d1.lt.x.and.x.le.-1.5d0)then
         m=iabs(int(4.0d0-5.2d1/x+3.45d2/x/x))
         w(4*m+3)=x*x
         do i=m,1,-1
            w(4*i-1)=x*x+dble(4*i-1)-dble(2*i*(2*i+1))/w(4*i+3)
         end do
         gascum=gasden(x)*(1.0d0/w(3)-1.0d0)/x
      endif
      if(3.0d0.lt.x.and.x.lt.1.d1)then
         m=iabs(int(-10+8.5d1/x+2.55d2/x/x))
         w(4*m+3)=x*x
         do i=m,1,-1
            w(4*i-1)=x*x+dble(4*i-1)-dble(2*i*(2*i+1))/w(4*i+3)
         end do
         gascum=1.0d0+gasden(x)*(1.0d0/w(3)-1.0d0)/x
      endif
      if(-1.5d0.lt.x.and.x.le.3.0d0)then
         m=int(1.d1+1.d1*abs(x))
         w(2*m+1)=x
         do i=m-1,1,-1
            w(2*i-1)=w(2*i+1)*x*x/dble(2*i+1)+x
         end do
         gascum=gasden(x)*w(1)+5.d-1
      endif
      return
      end

c  The following function returns Normal N(0,1) density function.
      real*8 function gasden(x)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter(pi=3.14159265358979323846d0)

      if(dabs(x).gt.1.d1)then
         gasden=0.0d0
      else
         gasden=dexp(-x*x/2.0d0)/dsqrt(2.0d0*pi)
      endif
      return
      end

c   this subroutine returns the inverse of the cumulative distribution 
c   of Normal(0,1) at probability p by the rational fraction approximation.
c   Reference: Statistical Computing, pp.95.
c   By William J. Kennedy, Jr., and James E. Gentle. 1980.

      subroutine gauinv(pp,xp,error)
      implicit real*8 (a-h,o-z),integer*4 (i-n)
      dlim=1.d-20
      p0=-3.22232431088d-1
      p1=-1.0d0
      p2=-3.42242088547d-1
      p3=-2.04231210245d-2
      p4=-4.53642210148d-5
      q0=9.93484626060d-2
      q1=5.88581570495d-1
      q2=5.31103462366d-1
      q3=1.03537752850d-1
      q4=3.8560700634d-3

      error=1.0d0
      xp=0.0d0

      p=pp
      if(p.gt.5.d-1)p=1.0d0-pp
      if(p.lt.dlim)return

      error=0.0d0
      if(p.eq.5.d-1)return

      y=dsqrt(-dlog(p*p))
      xp=y+((((y*p4+p3)*y+p2)*y+p1)*y+p0)/((((y*q4+q3)*y+q2)*
     +   y+q1)*y+q0)
      if(pp.lt.5.d-1)xp=-xp

      return
      end

c  This function uses ran1 to returns a normally N(0,1) r.v.
c  (using ran1(idum) as the source of uniform deviates).

      real*8 function gasdev(idum)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      save iset,gset
      data iset/0/
      if (iset.eq.0) then
1        v1=2.0d0*ran1(idum)-1.0d0
         v2=2.0d0*ran1(idum)-1.0d0
         rsq=v1**2.0d0+v2**2.0d0
         if(rsq.ge.1.0d0.or.rsq.eq.0.0d0)goto 1
         fac=dsqrt(-2.0d0*dlog(rsq)/rsq)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
      else
         gasdev=gset
         iset=0
      endif

      return
      end

c  The followng function uses ran1,gam1dev,gamfdev
c  Returns a deviate distributed as a Gamma distribution of parameters a and
c  b, using random ran1(idum),gam1dev(ia,idum),gamfdev(a,idum).

      real*8 function gamdev(a,b,idum)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

c      if(a.le.0.0d0.or.b.le.0.0d0)pause 'bad argument in gammadev'
      n=int(a)
      aa=a-dble(n)
      data x1/0.0d0/,x2/0.0d0/
      if(aa.gt.0.0d0)x1=gamfdev(aa,idum)
      if(n.gt.0)x2=gamidev(n,idum)
      gamdev=b*(x1+x2)
      return
      end

      real*8 function gamidev(ia,idum)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
c  The following function uses ran1 to returns a deviate distributed as a 
c  Gamma distribution of integer order, i.e., a
c  waiting time to the iath event in a Poisson process of unit mean, 
c  using ran1(idum)  as the source of uniform deviates.
c     if(ia.lt.1)pause 'bad argument in gammadev'

      if(ia.lt.6)then
        x=1.0d0
	      do j=1,ia
	         x=x*ran1(idum)
	      enddo
	      x=-dlog(x)
      else
1     v1=2.0d0*ran1(idum)-1.0d0
      v2=2.0d0*ran1(idum)-1.0d0
	    if(v1**2.0d0+v2**2.0d0.gt.1.0d0)goto 1
	    y=v2/v1
	    am=dble(ia-1)
	    s=dsqrt(2.0d0*am+1.0d0)
	    x=s*y+am
	    if(x.le.0)goto 1
	    e=(1.0d0+y**2.0d0)*dexp(am*dlog(x/am)-s*y)
	    if(ran1(idum).gt.e)goto 1
      endif
      gamidev=x
      return
      end

c  The following function uses ran1 to
c  returns a deviate distributed as a Gamma distribution of a real number
c  between 0 and 1, using the rejection method and Weibull distribution.

      real*8 function gamfdev(a,idum)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

c     if(a.le.0.0d0.or.a.ge.1.0d0)pause 'bad argument in gammadev'

        h=dexp(a**(a/(1.0d0-a))-a**(1.0d0/(1.0d0-a)))
1       v1=ran1(idum)
      	v2=h*ran1(idum)
	      x=(-dlog(v1))**(1.0d0/a)
	      y=dexp(x**a-x)
	      if(v2.ge.y)goto 1

	      gamfdev=x
       	return
	     end
      
c  The following function generates U(0,1) random numbers by the iterative 
c  process.  The program is copied from "Numerical Recipes, 2nd edition" on
c  page 271.  The method is using the recursive formula
c         x(n+1)=16807 x(n) mod (2147483647).
c
      DOUBLE PRECISION FUNCTION ran1(idum)      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      PARAMETER (IA=16807, IM=2147483647, AM=1.0d0/DBLE(REAL(IM)),
     +     IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2d-10,
     +     RNMX=1.0d0-EPS)
      INTEGER iv(NTAB)
      SAVE iv,iy
      DATA iv /NTAB*0/, iy/0/
      if (idum.le.0.or.iy.eq.0) then
	 idum=max(-idum,1)
	 do j=NTAB+8,1,-1
	     k=idum/IQ
	     idum=IA*(idum-k*IQ)-IR*k
	     if (idum.lt.0) idum=idum+IM
	     if (j.le.NTAB) iv(j)=idum
	 end do
	 iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*DBLE(REAL(iy)),RNMX)
      return
      end

c  The function gammln returns the value ln[Gamma(xx)] for xx>0
      real*8 function gammln(xx)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*8 cof(6)
c  Internal arithmetic will be done in double precision, a nicety that you can
c  omit if five-figure accuracy is good enough.
      save cof,stp
      data cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *  24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *  -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+5.d-1)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
         y=y+1.0d0
         ser=ser+cof(j)/y
      end do
      gammln=tmp+dlog(stp*ser/x)
      return
      end

C  The following subroutine sorts an array arr(1:n) into ascending numerical 
c  order using the Quicksort algorithm.  n is input; 
c  arr is replaced on output by its sorted rearrangement.
C  Parameters: M is the size of subarray sorted by straught insertion and 
C  NSTACK is the required auxiliary storage.

      subroutine sort(n,arr)
      integer n,M,NSTACK
      REAL arr(n)
      PARAMETER(M=7,NSTACK=50000)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
	do j=l+1,ir
	  a=arr(j)
	  do i=j-1,1,-1
	    if(arr(i).le.a)goto 2
	    arr(i+1)=arr(i)
        end do
	  i=0
2       arr(i+1)=a
      end do
	if(jstack.eq.0)return
	   ir=istack(jstack)
	   l=istack(jstack-1)
	   jstack=jstack-2
      else
	   k=(l+ir)/2
	   temp=arr(k)
	   arr(k)=arr(l+1)
	   arr(l+1)=temp
	   if(arr(l+1).gt.arr(ir))then
	     temp=arr(l+1)
	     arr(l+1)=arr(ir)
	     arr(ir)=temp
	   endif
	   if(arr(l).gt.arr(ir))then
	     temp=arr(l)
	     arr(l)=arr(ir)
	     arr(ir)=temp
	   endif
	   if(arr(l+1).gt.arr(l))then
	     temp=arr(l+1)
	     arr(l+1)=arr(l)
	     arr(l)=temp
	   endif
	   i=l+1
	   j=ir
	   a=arr(l)
3        continue
	     i=i+1
	   if(arr(i).lt.a)goto 3
4          continue
	     j=j-1
	   if(arr(j).gt.a)goto 4
	   if(j.lt.i)goto 5
         temp=arr(i)
	   arr(i)=arr(j)
	   arr(j)=temp
	   goto 3
5        arr(l)=arr(j)
	   arr(j)=a
	   jstack=jstack+2

C  Push pointers to larger subarray on stack, process smaller subarray
C  immediately.
c	   if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'

	   if(ir-i+1.ge.j-1)then
	     istack(jstack)=ir
	     istack(jstack-1)=i
	     ir=j-1
	   else
	     istack(jstack)=j-1
	     istack(jstack-1)=l
	     l=i
	   endif
      endif
      goto 1
      end

c  The following function generates pseudorandom numbers from a discrete 
c  uniform distribution (1,k) using ran1(idum). idum is the seed
      integer*4 function myund(k,idum)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      
      myund = AINT(k*ran1(idum)+1.0)
      
      return
      end	  
