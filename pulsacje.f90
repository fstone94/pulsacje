!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!                                -------??CO ROBI TEN PROGRAM??-------

!Zadaniem programu jest wygenerowanie animacji *.gif przedstawiajacej drgania powiezchni gwiazdy w czasie 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

PROGRAM pulsacje
IMPLICIT NONE

!-------------------------------------------------------------------------
!SPECYFIKACJA ZMIENNYCH GLOBALNYCH
!-------------------------------------------------------------------------

integer, parameter     :: dp=selected_real_kind(15)       !specyfikacja dokladnosci typy REAL
integer,parameter      :: pp=selected_real_kind(15,80)
real,parameter         :: pi=4*atan(1.0)                  !liczba pi  
integer                :: ilosc_arg                       !ilosc argumentow przekazanych z lini komend 
character(len=30)      :: parametr                        !zmienna przechowujaca parametry z lini komend
character(len=30)      :: filein                          !nazwa pliku zawierajacego parametry modow 
character(len=30)      :: fileout                         !nazwa pliku wyjsciowego zawierajacego animacje 
real(kind=pp)          :: t_total                         !calkowita dlugosc czasowa animacji wyrazona w dniach 
real(kind=pp)          :: delta_t                         !krok czasowy animacji wyrazony w dniach 
integer                :: n_teta                          !ilosc dyskretnych punktow rozlozonych w dlugosci asterograficznej gwiazdy 
integer                :: n_fi                            !ilosc dyskretnych punktow rozlozonych w szerokosci asterograficznej gwiazdy         
real(kind=pp)          :: t_i 

!-------------------------------------------------------------------------
!POBRANIE PARAMETROW Z LINI KOMEND ORAZ WYSWIETLANIE HELPA
!-------------------------------------------------------------------------

!patrzymy ile argumentow zostal‚o przekazanych z lini komend
ilosc_arg = command_argument_count()

!jesli ilosc argumentow wynosi 0 lub jest inna niz rzadana, wowczas wyswietlamy helpa
IF (ilosc_arg == 0 .OR. ilosc_arg /= 6) THEN
   write(*,*)
   write(*,*) '###############################################################################'
   write(*,*) '                            *** PROGRAM PULSACJE ***                    '
   write(*,*) '###############################################################################'
   write(*,*) 
   write(*,*) 'AUTOR    : MARTA KAMINSKA'
   write(*,*) 
   write(*,*) '-------------------------------------------------------------------------------'
   write(*,*) 
   write(*,*) 'WYWOLANIE: >./pulsacje <filein> <fileout> <t_total> <delta_t> <n_teta> <n_fi>'
   write(*,*)
   write(*,*) '-------------------------------------------------------------------------------'
   write(*,*) 
   write(*,*) '<filein>   --->  nazwa pliku zawierajacego parametry modow '
   write(*,*) '<fileout>  --->  nazwa pliku wyjsciowego zawierajacego animacje'
   write(*,*) '<t_total>  --->  calkowita dlugosc czasowa animacji wyrazona w dniach '
   write(*,*) '<delta_t>  --->  krok czasowy animacji wyrazony w dniach '
   write(*,*) '<n_teta>   --->  ilosc punktow w dlugosci asterograficznej gwiazdy '
   write(*,*) '<n_fi>     --->  ilosc punktow w szerokosci asterograficznej gwiazdy '
   write(*,*) '-------------------------------------------------------------------------------' 
   STOP
END IF

!jesli ilosc argumentow sie zgadza, wowczas pobieramy ich wartosci do odpowiednich zmiennych
call get_command_argument(1,parametr)
read(parametr,*) filein
call get_command_argument(2,parametr)
read(parametr,*) fileout
call get_command_argument(3,parametr)
read(parametr,*) t_total 
call get_command_argument(4,parametr) 
read(parametr,*) delta_t
call get_command_argument(5,parametr)
read(parametr,*) n_teta
call get_command_argument(6,parametr) 
read(parametr,*) n_fi 


!przywoluje subrutyne obliczajaca wartosc wektorow wodzacych dla kazdego z dyskretnych punktow
!oraz zapisujaca wartosci do plikow, rysujaca obrazki, tworzaca gifa  

call wektory(filein,t_total,delta_t,n_teta,n_fi)

!-------------------------------------------------------------------------
!KOMENTARZE KONCOWE
!-------------------------------------------------------------------------

!jesli wszysto przebiegnie pomyslnie, wyswietlony zostanie komunikat 
write(*,*)'---------------------------------------------------------------'
write(*,*)'           !! WSZYSTKO POSZLO JAK TRZEBA !!                    '
write(*,*)'---------------------------------------------------------------'
write(*,*)
write(*,*)'   WYNIK DZIALANIA PROGRAMU ZOSTAL ZAPISANY NA DYSKU           '
write(*,*)
write(*,*)'---------------------------------------------------------------' 


!#########################################################################
CONTAINS
!#########################################################################

!-------------------------------------------------------------------------
!SUBRUTYNA OBLICZAJACA STOWARZYSZONE WIELOMIANY LEGENDRA 	
!-------------------------------------------------------------------------

SUBROUTINE pm_polynomial_value ( mm, n, m, x, y )

IMPLICIT NONE 

  integer ( kind = 4 ) mm
  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(mm,0:n)
  real ( kind = 8 ) y(mm)
  real ( kind = 8 ) fact
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(mm)

  cx(1:mm,0:n) = 0.0D+00
!
!  J = M is the first nonzero function.
!
  if ( m <= n ) then
    cx(1:mm,m) = 1.0D+00

    fact = 1.0D+00
    do j = 1, m
      cx(1:mm,m) = - cx(1:mm,m) * fact * sqrt ( 1.0D+00 - x(1:mm)**2 )
      fact = fact + 2.0D+00
    end do

  end if
!
!  J = M + 1 is the second nonzero function.
!
  if ( m + 1 <= n ) then
    cx(1:mm,m+1) = x(1:mm) * real ( 2 * m + 1, kind = 8 ) * cx(1:mm,m)
  end if
!
!  Now we use a three term recurrence.
!
  do j = m + 2, n
    cx(1:mm,j) = ( real ( 2 * j     - 1, kind = 8 ) * x(1:mm) * cx(1:mm,j-1) &
                 + real (   - j - m + 1, kind = 8 ) *           cx(1:mm,j-2) ) &
                 / real (     j - m,     kind = 8 )
  end do

  y = cx(1:mm,n)

  return
  
END SUBROUTINE pm_polynomial_value

!----------------------------------------------------------------------------------
!SUBRUTYNA OBLICZAJACA WARTOSC WEKTORA WODZACEGO 
!----------------------------------------------------------------------------------

SUBROUTINE wektory(filein,t_total,delta_t,n_teta,n_fi) 
IMPLICIT NONE 

integer, parameter                        :: dp=selected_real_kind(15)              !specyfikacja dokladnosci typy REAL
integer,parameter                         :: pp=selected_real_kind(15,80)           !specyfikacja dokladnosci dla silnii 
integer,parameter                         :: mm=1 
real,parameter         :: pi=4*atan(1.0)  
character(len=30)                         :: filein,nazwapng                        !plik wejsciowy i nazwa obrazka w png                           
character(len=100)                        :: title                                  !zmienna przechowujaca tytul kazdego obrazka 
real(kind=pp)                             :: t_total                                !calkowity czas dla animacji 
real(kind=pp)                             :: argument                               !zmienna posrednia przechowujace dane do obliczenia skladowych wzoru na r_i
real(kind=pp)                             :: delta_t,delta_fi,delta_teta            !krok dla czasu, dla kata phi i dla kata theta 
real(kind=pp)                             :: r_i,phi_i,theta_i,theta_gnuplot,teta_i,fi_i !wektor wodzacy, kat phi, kat theta, theta przeksztalcona dla gnuplota, chwila czasu t_i 
real,dimension(:,:),allocatable           :: mody                                   !tablica przechowujaca parametry modow wczytywane z pliku wejsciowego
real(kind=8),dimension(:),allocatable     :: cosinus_teta,cosinus_m_fi,cosinus_t_i  !wektory przechowujace skladowe cosinusy wzoru na r_i
real(kind=8),dimension(:),allocatable     :: wielomiany,amplituda                   !wektory przechowujace wielomiany Legendra i wartosci amplitudy z pliku wejsciowego
integer                                   :: n_teta,l_obr,i,j,k,q,n_fi              !liczba punktow theta i phi; liczba obrazkow a animacji; numeratory petli                 
integer                                   :: liczba_wierszy,liczenie_wierszy        !liczba wierszy w pliku wejsciowym policzona przez funkcje liczenie_wierszy    
 



!liczba wierszy obliczona przez zewnetrzna funkcje 
liczba_wierszy=liczenie_wierszy(filein)

!------------------------------------------------------------------------
!WCZYTANIE WARTOSCI PARAMETROW Z PLIKU WEJSCIOWEGO
!------------------------------------------------------------------------

!alokuje tablice o 5 kolumnach i liczbie wierszy obliczonej powyzej 
ALLOCATE(mody(5,liczba_wierszy)) 
OPEN(unit=2,file=filein) 

!wczytuje dane z pliku do tablicy 
READ(2,*,END=2)mody
2 Close(2)

!------------------------------------------------------------------------
!WYZNACZENIE WARTOSCI POTRZEBNYCH DO OBLICZENIA WEKTOROW WODZACYCH
!------------------------------------------------------------------------

!liczba generowanych obrazkow odpowiada ilorazowi calkowitego czasu i delty t 
!l_obr jest liczba calkowita
l_obr=t_total/delta_t 

!kat theta przebiega od 0 do pi
!kat phi przebiega od 0 do 2pi 
!delta dla obu katow odpowiada wiec ilorazowi dlugosci zakresu i wpisanych z linii komend punktow dyskretnych 

delta_teta= pi/n_teta
delta_fi=(2*pi)/n_fi

!-------------------------------------------------------------------------
!PETLA ODPOWIEDZIALNA ZA OBLICZENIE WIELOMIANOW LEGENDRA I COSINUSOW  
!-------------------------------------------------------------------------

!pierwszy wektor bedzie przechowywal wartosci cosinusow 
!do drugiego wektora zapisane zostana obliczone wartosci wielomianow dla kazdego cosinusa 
ALLOCATE(cosinus_teta(n_teta)) 
ALLOCATE(wielomiany(n_teta)) 

DO i=1,n_teta                                    !iteracja po liczbie punktow dyskretnych 

	teta_i=(i-1)*delta_teta
	cosinus_teta(i)=cos(teta_i) 

		DO j=1,liczba_wierszy                    !iteracja po liczbie modow w pliku wejsciowym 
		
			call pm_polynomial_value (mm,int(mody(4,j)),int(mody(5,j)),cosinus_teta(i),wielomiany(i))
			
		END DO 
END DO 
		

!obliczenie cosinusa m_phi i wrzucenie wartosci do wektora
ALLOCATE(cosinus_m_fi(n_fi))

DO i=1,n_fi
 
	fi_i=(i-1)*delta_fi
		
		DO j=1,liczba_wierszy 	
				
			cosinus_m_fi(i)=cos(fi_i*mody(5,j))
		
		END DO 
		
END DO 


!obliczenie cosinusa (2*pi*...) i wrzucenie do wektora 
ALLOCATE(cosinus_t_i(l_obr))

DO i=1,l_obr

	t_i=(i-1)*delta_t 
	argument=2*pi*mody(1,1)*t_i

	cosinus_t_i(i)=cos(argument)
	
END DO 
 
!tworze wektor przechowujacy wartosci amplitudy, co przyda sie przy ustalaniu zakresu osi  
ALLOCATE(amplituda(liczba_wierszy))
 
DO i=1,liczba_wierszy 

	amplituda(i)=mody(2,1)


END DO 
 
 
!--------------------------------------------------------------------------
!CLUE PROGRAMU CZYLI WYLICZENIE WARTOSCI r_i 
!--------------------------------------------------------------------------

DO i=1,l_obr                                            !petla po chwilach czasu, ktore odpowiadaja ilosci obrazow

	t_i=(i-1)*delta_t 
	!otwieram pik do ktorego wpisze komendy do gnuplopta
	OPEN(unit=3,file="skrypt.plt") 

!kazdy obrazek bedzie mial inna nazwe tak, zeby sie nie nadpisywaly 
!na kazdym obrazku umieszczona bedzie inna chwila czasu t_i
		nazwapng='set output  "'//'sfera'//trim(numer(i))//'.png'//'"'
		WRITE(title,222) t_i
	    
write(3,*) 'set terminal pngcairo size 1600,1000 enhanced'
write(3,*) trim(nazwapng) 
write(3,*) 'set pm3d'
write(3,*) 'set hidden3d'
write(3,*) 'set palette'
write(3,*) 'set xrange [',-1-sum(amplituda),':',1+sum(amplituda),']'
write(3,*) 'set yrange [',-1-sum(amplituda),':',1+sum(amplituda),']'
write(3,*) 'set zrange [',-1-sum(amplituda),':',1+sum(amplituda),']'
write(3,*) 'set cbrange [',1-0.7*sum(amplituda),':',1+0.7*sum(amplituda),']'
write(3,*) 'set mapping spherical'
write(3,*) 'set border linewidth 2'
write(3,*) 'unset key'
write(3,*) 'set mxtics'
write(3,*) 'set mytics'
write(3,*) 'set mztics'
write(3,*) 'set mcbtics'
write(3,*) 'set xtics font "Helvetica,17" offset 0,-1,0'
write(3,*) 'set ytics font "Helvetica,17" offset 2,0,0'
write(3,*) 'set ztics font "Helvetica,17" offset 0,0,0'
write(3,*) 'set cbtics font "Helvetica,17"'
write(3,*) 'set format x "%.1f"'
write(3,*) 'set format y "%.1f"'
write(3,*) 'set format z "%.1f"'
write(3,*) 'set format cb "%.2f"'
write(3,*) 'set xlabel "(X/R_{star})" font "Helvetica,20"'
write(3,*) 'set ylabel "(Y/R_{star})" font "Helvetica,20"'
write(3,*) 'set zlabel "(Z/R_{star})" font "Helvetica,20" offset -6'
write(3,*) 'set cblabel "Radial Displacement (1+{/Symbol D}r/R_{star})" font "Helvetica,20" offset 5'
write(3,*)  trim(title)
write(3,*) 'set view 65,45,1.1'
write(3,*) 'set view equal xyz'
write(3,*) 'set ticslevel 0.1'
write(3,*) 'splot "sfera.dat" u 1:2:3:3 w p lt 7 ps 0.9 palette'

222 format('set title "Elapsed Time = ',F9.3,' d" font "Helvetica,30" offset 0,0,-2')

	Close(3) 
	
	!tworze plik do ktorego zapisane zostana wyliczone parametry 
	OPEN(unit=4,file='sfera.dat')

	DO j=1,n_teta                   !petla po kacie theta
	 
		DO k=1,n_fi                 !petla po kacie phi 
		
			r_i=1 
		
			DO q=1,liczba_wierszy   !petla po wszystkich modach 
			
			
			!obliczam theta i zmieniam na czytelne dla gnuplota 
			theta_i=(j-1)*delta_teta
			theta_gnuplot=((pi/2)-theta_i)
			
		    !obliczam phi 
			phi_i=(k-1)*delta_fi

			!obliczam promien wodzacy 
			r_i=1+mody(2,q)*sqrt(((2*mody(4,q)+1)/(4*pi))*&
			&((factorial(int(mody(4,q))-int(mody(5,q))))/(factorial(int(mody(4,q))+int(mody(5,q))))))&
			&*wielomiany(j)*cosinus_m_fi(k)*cosinus_t_i(i)
			

			END DO
			
			!zapisuje obliczone wartosci do pliku 
			!plik zostaje nadpisany przez kolejna petle, ale jako ze rysowane obrazki maja inne nazwy to jest ok 
			WRITE(4,*)phi_i,theta_gnuplot,r_i
			
		END DO 
		
	END DO 
	
	!wywoluje gnuplota i wrzucam do niego utworzony skrypt do rysowania obrazka 
	call execute_command_line('gnuplot < skrypt.plt')
	
	!zamykam plik do ktorego wpisywaly sie dane 
	CLOSE(4)
	
END DO

    !lacze wszystkie obrazki w jednego gifa 
    call execute_command_line('convert -delay 10 -loop 0 *.png '//trim(fileout)//'.gif')
    !usuwam zarowno plik ze skryptem jak i wszystkie skladowe z formatem .png
    call execute_command_line('rm *.png *.plt')
    
END SUBROUTINE wektory

!-------------------------------------------------------------------------
!INTEGER NA CHARACTER
!-------------------------------------------------------------------------

!funkcja zamieniajaca integer na character 
character(len=30)FUNCTION numer(s)

integer,intent(in) ::s

WRITE(numer,*)s
numer=adjustl(numer)

END FUNCTION numer
	
!------------------------------------------------------------------------
!SILNIA
!------------------------------------------------------------------------	
	
!funkcja obliczajaca silnie 
real(kind=pp) FUNCTION factorial(n) !beda wychodzic bardzo duze liczby, potrzebne jest wiec zeby 
IMPLICIT NONE 

real(kind=pp)                     :: j
integer                           :: n

factorial=1.0 
DO j=2,n,1.0  
 
   factorial=factorial*j 
   
END DO 

END FUNCTION factorial 	
		
				
END PROGRAM pulsacje 

!------------------------------------------------------------------------
!LICZENIE WIERSZY
!------------------------------------------------------------------------

!funkcja, ktora zwraca liczbe wierszy w pliku wejsciowym 
integer FUNCTION liczenie_wierszy(filein)
IMPLICIT NONE   
character(len=30)                :: filein 

OPEN(unit=1,file=trim(filein)) 

liczenie_wierszy=0 
!czyta plik zwiekszajac licznik n o 1 az do konca pliku 
DO 

		READ(1,*,END=1)
		liczenie_wierszy=liczenie_wierszy+1 
		
END DO 
1 CLOSE(1) 	!kiedy skoncza sie dane plik zostanie zamkniety 

END FUNCTION liczenie_wierszy	

