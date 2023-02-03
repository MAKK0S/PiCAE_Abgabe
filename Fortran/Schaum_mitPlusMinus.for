c      Titel: Isotrope, viskoelastische UMAT im Zeitbereich
c      zu beachten: gleiche Relaxationszeiten für die
c                   unterschiedlichen Materialkennwerte müssen
c                   verwendet werden
c                   Es werden 6+6 Zustandsvariablen benötigt
c
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,knc)
c------Deklaration ABAQUS
      implicit none
      integer :: kspt,layer,npt,noel,nprops,nstatv,ntens,
     1 nshr,ndi,knc
      double precision :: sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,
     1 pnewdt,celent,dfgrd0(3,3),dfgrd1(3,3),time(2),stress(ntens),
     2 statev(nstatv),ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     3 stran(ntens),dstran(ntens),predef(1),dpred(1),props(nprops),
     4 coords(3),drot(3,3)
      integer :: jstep(4)
      character*80 cmname
c------Lokale Deklarationen
      integer :: i, j, l, AnzahlMaxwellElem
      real :: stressM(ntens),stressE(ntens),
     1 zero,one,two,three,four,half,ten,expTerm
      real :: E,nu,G0,
     1 k0,G,k,Gm,km,f0(ntens,ntens),
     2 f(ntens,ntens)
      real :: gi((nprops-2)/3),taui((nprops-2)/3),ki((nprops-2)/3)
      
c
c-----Nuetzliche Zahlen
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, four=4.d0,
     1  half=0.5d0,ten=10.d0)
c
c-----Initialisieren mit null
      AnzahlMaxwellElem=(nprops-2)/3
      gi=zero; ki=zero; taui=zero
c-----Materialparameter aus props einlesen
      E=props(1); nu=props(2)
      gi(1:AnzahlMaxwellElem)=props(3:nprops-2:3)
      ki(1:AnzahlMaxwellElem)=props(4:nprops-1:3)
      taui(1:AnzahlMaxwellElem)=props(5:nprops:3)
c-----Zuweisen der Materialparamter
      G0=E/(two*(one+nu)); k0=E/(three*(one-two*nu))
      G=G0
      k=k0
      do l=1,AnzahlMaxwellElem
        G=G-G0*gi(l)
        k=k-k0*ki(l)
      end do
c-----elastische Steifigkeitsmatrix
     f0=zero
      do i=1,2
       f0(i,i)= four*G*(G+three*k)/(four*G+three*k)
      end do
      f0(2,1)= -two*G*(two*G-three*k)/(four*G+three*k)
      f0(1,2)=f0(2,1)
      f0(3,3)=G
      km=zero;Gm=zero;f=zero
      Gm=gi*G0
      km=ki*k0
c
c-----Berechnung des elastischen Spannungsanteils
      stressE=zero
      do i=1,3
        stressE(i)=statev(i)
          do j=1,3
            stressE(i)=stressE(i)+f0(i,j)*dstran(j)
          end do
      end do
  
c 
c----Aktualisieren des elastischen Anteils der Statusvariablen, ddsdde und Spannung
      statev(1:6)=stressE(1:6)
      ddsdde=zero
      ddsdde(1,1)=(k+four/three*G)
      ddsdde(1,2)=(k-two/three*G)
      
      stress=zero
      stress(1:6)=stressE(1:6)
      WRITE(*,*)'stressE',stressE
      WRITE(*,*)'stress',stress
c
      stressM=zero
      do l=1,AnzahlMaxwellElem
c----Berechnung des Relaxationsanteils
        stressM=zero
        do j=1,ntens
          stressM(j)=exp(-dtime/taui(l))*statev(6+j+6*(l-1))
        end do
c-------Steifigkeitsmatrix für n-tes Maxwell-Element
        km=zero;Gm=zero;f=zero
        Gm=gi(l)*G0
        km=ki(l)*k0
        do i=1,3
          do j=1,3
            f(i,j)=km-two/three*Gm
          end do
          f(i,i)=km+four/three*Gm
        end do
    
c------Berechnung des übrigen Anteils
        expTerm=taui(l)*(one-exp(-dtime/taui(l)))
        do i=1,3
          do j=1,3
            stressM(i)=stressM(i)+f(i,j)*dstran(j)/dtime*expTerm
          end do
        end do
       
c-----Aktualisieren der viskosen Zustandsvariable
        do i=1,ntens
          statev(6+i+6*(l-1))=stressM(i)
        end do
c-----Aktualisieren der Spannung
        do i=1,3
          stress(i)=stress(i)+stressM(i)
        end do
c
c-----Aktualisieren des rests von ddsdde
        Gm=gi(l)*G0
        km=ki(l)*k0
        ddsdde(1,1)=ddsdde(1,1)+(km+four/three*Gm)*taui(l)/dtime*(one-exp(-(dtime/taui(l))))
        ddsdde(1,2)=ddsdde(1,2)+(km-two/three*Gm)*taui(l)/dtime*(one-exp(-(dtime/taui(l))))
       
      end do
      ddsdde(2,2)=ddsdde(1,1)
      ddsdde(3,3)=ddsdde(1,1)
      ddsdde(1,3)=ddsdde(1,2)
      ddsdde(2,1)=ddsdde(1,2)
      ddsdde(2,3)=ddsdde(1,2)
      ddsdde(3,1)=ddsdde(1,2)
      ddsdde(3,2)=ddsdde(1,2)
    
c
c-----Berechnen der sphärischen Spannung (bzw. tr(sigma)) und speichern als zusätzliche innere Variable
      IF ((stress(1)+stress(2)+stress(3)) >= zero) THEN
        statev(6*AnzahlMaxwellElem+7) = one
      ELSE
        statev(6*AnzahlMaxwellElem+7) = -1.d0
      END IF
      end subroutine umat
