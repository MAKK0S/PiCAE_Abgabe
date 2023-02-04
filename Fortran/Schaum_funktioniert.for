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
      integer :: i, j, l, nmaxwell
      real :: stressm(ntens),stressE(ntens),
     1 zero,one,two,three,four,half,ten,expTerm,trace_sigma
      real :: E,nu,gi(nprops),taui(nprops),ki(nprops),G0,
     1 k0,G,k,Gm,km,f0(ntens,ntens),
     2 f(ntens,ntens)
c
c-----Nuetzliche Zahlen
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, four=4.d0,
     1  half=0.5d0,ten=10.d0)
c-----Anzahl MW-Elemente bestimmen
      nmaxwell = (nprops-2)/3
c-----Initialisieren mit null
      gi=zero; ki=zero; taui=zero
c-----Materialparameter aus props einlesen 
      E=props(1); nu=props(2)
      gi(1:nmaxwell)=props(3:nprops-2:3)
      ki(1:nmaxwell)=props(4:nprops-1:3)
      taui(1:nmaxwell)=props(5:nprops:3)

c-----Zuweisen der Materialparamter
      G0=E/(two*(one+nu)); k0=E/(three*(one-two*nu))
      G=G0
      k=k0
      do l=1,nmaxwell
        G=G-G0*gi(l)
        k=k-k0*ki(l)
      end do
c-----Fallunterscheidung Schale-Solid
      IF ( ntens > three ) THEN
c-----elastische Steifigkeitsmatrix
        f0=zero
        do i=1,3
            do j=1,3
            f0(i,j)=k-two/three*G
            end do
            f0(i,i)=k+four/three*G
        end do
        do i=4,6
            f0(i,i)=G
        end do
c-----Elastischer Anteil Steifigkeit zuweisen
        ddsdde=zero
        ddsdde(1,1)=(k+four/three*G)
        ddsdde(1,2)=(k-two/three*G)
        ddsdde(4,4)=G
c-----Berechnung des elastischen Spannungsanteils
        stressE=zero
        do i=1,3
            stressE(i)=statev(i)
            do j=1,3
                stressE(i)=stressE(i)+f0(i,j)*dstran(j)
            end do
        end do
        do j=4,ntens
            stressE(j)=statev(j)+f0(j,j)*dstran(j)
        end do
c-----Update der Spannung
        stress = zero
        stress(1:6)=stressE(1:6)
c-----Aktualisieren der gespeicherten Zustandsvariablen
        statev(1:6)=stressE(1:6)

c-----Fuer alle Maxwell-Elemente iterieren      
        do l=1,nmaxwell
            km=zero;Gm=zero;f=zero
            Gm=gi(l)*G0
            km=ki(l)*k0
c-----Komponenten der Steifigkeitsmatrix für jedes Prony-Element
            do i=1,3
                do j=1,3
                f(i,j)=km-two/three*Gm
                end do
                f(i,i)=km+four/three*Gm
            end do
            do i=4,6
                f(i,i)=Gm
            end do
c-----Relaxieren der Spannung
            stressm=zero
            do j=1,ntens
                stressm(j)=exp(-dtime/taui(l))*statev(6+j+6*(l-1))
            end do
c-----Berechnung der Spannung am Ende des Inkrementes
            expTerm=taui(l)*(one-exp(-dtime/taui(l)))
c
            do i=1,3
                do j=1,3
                stressm(i)=stressm(i)+f(i,j)*dstran(j)/dtime*expTerm
                end do
            end do
            do j=4,ntens
                stressm(j)=stressm(j)+f(j,j)*dstran(j)/dtime*expTerm
            end do
c-----Berechnung der Gesamtspannung am Ende des Inkrementes
            do i=1,6
                stress(i)=stress(i)+stressm(i)
            end do
c
c-----Jacobi-Matrix
            ddsdde(1,1)=ddsdde(1,1)+(km+four/three*Gm)*taui(l)/dtime*(one-exp(-(dtime/taui(l))))
            ddsdde(1,2)=ddsdde(1,2)+(km-two/three*Gm)*taui(l)/dtime*(one-exp(-(dtime/taui(l))))
            ddsdde(4,4)=ddsdde(4,4)+Gm*taui(l)/dtime*(one-exp(-(dtime/taui(l))))

            ddsdde(2,2)=ddsdde(1,1)
            ddsdde(3,3)=ddsdde(1,1)
            ddsdde(1,3)=ddsdde(1,2)
            ddsdde(2,1)=ddsdde(1,2)
            ddsdde(2,3)=ddsdde(1,2)
            ddsdde(3,1)=ddsdde(1,2)
            ddsdde(3,2)=ddsdde(1,2)
            ddsdde(5,5)=ddsdde(4,4)
            ddsdde(6,6)=ddsdde(4,4)
c
c-----Aktualisieren der gespeicherten Zustandsvariablen
            do i=1,ntens
                statev(6+i+6*(l-1))=stressm(i)
            end do
c
        end do
    
      ELSE
c-----elastische Steifigkeitsmatrix
        f0=zero
        f0(1,1)=four*G*(G+three*K) / (four*G+3*K)
        f0(1,2)=-two*G*(two*G-three*K) / (four*G+3*K)
        f0(2,1)=f0(1,2)
        f0(2,2)=f0(1,1)
        f0(3,3)=G
c-----Elastischer Anteil Steifigkeit zuweisen
        ddsdde=zero
        ddsdde(1,1)=f0(1,1)
        ddsdde(1,2)=f0(1,2)
        ddsdde(3,3)=f0(3,3)
c-----Berechnung des elastischen Spannungsanteils
        stressE=zero
        do i=1,2
            stressE(i)=statev(i)
            do j=1,2
                stressE(i)=stressE(i)+f0(i,j)*dstran(j)
            end do
        end do
        do j=3,ntens
            stressE(j)=statev(j)+f0(j,j)*dstran(j)
        end do
c-----Update der Spannung
        stress=zero
        stress(1:3)=stressE(1:3)
c-----Aktualisieren der gespeicherten Zustandsvariablen
        statev(1:3)=stressE(1:3)

c-----Fuer alle Maxwell-Elemente iterieren      
        do l=1,nmaxwell
            km=zero;Gm=zero;f=zero
            Gm=gi(l)*G0
            km=ki(l)*k0
c-----Komponenten der Steifigkeitsmatrix für jedes Prony-Element
            f(1,1)=four*Gm*(Gm+three*km) / (four*Gm+3*km)
            f(1,2)=-two*Gm*(two*Gm-three*km) / (four*Gm+3*km)
            f(2,1)=f(1,2)
            f(2,2)=f(1,1)
            f(3,3)=Gm
c-----Relaxieren der Spannung
            stressm=zero
            do j=1,ntens
                stressm(j)=exp(-dtime/taui(l))*statev(ntens+j+ntens*(l-1))
            end do
c-----Berechnung der Spannung am Ende des Inkrementes
            expTerm=taui(l)*(one-exp(-dtime/taui(l)))
c
            do i=1,2
                do j=1,2
                stressm(i)=stressm(i)+f(i,j)*dstran(j)/dtime*expTerm
                end do
            end do
            do j=3,ntens
                stressm(j)=stressm(j)+f(j,j)*dstran(j)/dtime*expTerm
            end do
c-----Berechnung der Gesamtspannung am Ende des Inkrementes
            do i=1,6
                stress(i)=stress(i)+stressm(i)
            end do
c
c-----Jacobi-Matrix
            ddsdde(1,1)=ddsdde(1,1)+f(1,1)*expTerm/dtime
            ddsdde(1,2)=ddsdde(1,2)+f(1,2)*expTerm/dtime
            ddsdde(3,3)=ddsdde(3,3)+f(3,3)*expTerm/dtime            
c-----Aktualisieren der gespeicherten Zustandsvariablen
            do i=1,ntens
                statev(ntens+i+ntens*(l-1))=stressm(i)
            end do
c
        end do        
        ddsdde(2,1)=ddsdde(1,2)
        ddsdde(2,2)=ddsdde(1,1)
      END IF


      end subroutine umat
