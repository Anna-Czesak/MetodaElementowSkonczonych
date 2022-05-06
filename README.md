# Symulacja nieustalonych procesów cieplnych
Anna Czesak 
20.01.2022r. 

Metoda elementów skończonych

Problem:
Rozwiązanie zagadnienia nieustalonej wymiany ciepła w układzie 
dwuwymiarowym.

Cel:
Obliczenie temperatur wewnątrz materiału w określonej chwili.

Grid – klasa zawierająca, najważniejsze dane dla siatki oraz jej tworzenie. W obiekcie tej 
klasy utworzyłam większość metod do obliczania macierzy globalnych potrzebnych do 
układu równań. Tu również znajduje się metoda wykonująca symulację.

Element – klasa zawierająca tablicę z Id elementów, a także macierze lokalne H, C, Hbc oraz 
lokalny wektor P. Na potrzeby obliczeń zawarłam tu również jakobian, pochodne funkcji 
kształtu względem x oraz y.

Element_uniwersalny – klasa jest odpowiedzialna za przechowywanie danych dotyczących 
schematu całkowania: punktów całkowania, ich ilości, wag. W tym miejscu obliczam również 
funkcje kształtu, ich pochodne po ksi oraz eta. Element uniwersalny zawiera 4 ściany.

Node – klasa zawiera współrzędne węzłów, ich temperaturę oraz informację o obecności 
warunku brzegowego.

Sciana – klasa przedstawiająca ścianę elementu, jej punkty całkowania, wagi, funkcje 
kształtu, a także macierz hbc.

GaussElimination – klasa pomocnicza, w której znajdują się metody pozwalające obliczyć 
układ równań


MES
Metoda elementów skończonych polega na przedstawieniu danego obszaru materiału jako 
siatkę podzieloną na elementy, a temperaturę jako wartości funkcji kształtu w kolejnych 
węzłach tych elementów.



