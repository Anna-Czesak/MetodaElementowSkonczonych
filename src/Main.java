import java.util.Scanner;


public class Main {


    static public void main(String [] args){


   Scanner scanner = new Scanner(System.in);

        double H, B;
        int nH,nB;

      /*
      H=0.025;
      B=0.025;
      nH=2;
      nB=2;*/

      H=0.1;
      B=0.1;
      nH=31;
      nB=31;

       // H=0.1;
       // B=0.1;
       // nH=4;
       // nB=4;

        Grid siatka = new Grid(H,  B,  nH,  nB);

        siatka.addElement();
        siatka.addNode();

//2 lub 3 punktowy schemat całkowania
    int k=3;

//pochodne funkcji kształtu po dKsi i dEta znajduja sie w elemencie uniwersalnym
    Element_uniwersalny elementUniwersalny=new Element_uniwersalny(k);
    elementUniwersalny.printdNdKsi();
    elementUniwersalny.printdNdEta();

//obliczanie jakobianu
        for(int i=0; i<siatka.nE;i++){
            siatka.elements[i].obliczJakobian(elementUniwersalny,siatka,i);
            siatka.elements[i].obliczWyznacznik(elementUniwersalny);
            siatka.elements[i].obliczJakobianOdwrotny(elementUniwersalny);
            siatka.elements[i].obliczdNdXinNnY(elementUniwersalny);
        }
        for(int i=0; i<siatka.nE;i++){
            System.out.println("\nElement "+(i+1));
            siatka.elements[i].printJakobian(elementUniwersalny);
            siatka.elements[i].printPochodneDNdxDNdy(elementUniwersalny);
        }
//obliczanie macierzy H

        int kt=25; //współczynnik przewodnosci cieplnej
        for(int i=0; i<siatka.nE;i++){
            System.out.println("\n\nElement nr: "+(i+1));
            siatka.elements[i].obliczMacierzH(kt,elementUniwersalny);
            siatka.elements[i].printMacierzH();
        }

//obliczanie macierzy Hbc oraz wektora P

        int alpha=300; //współczynnik konwekcji [W/m2K] laby -25 ostatecznie 300
        int t=1200; //temperatura otoczenia do wektora P

        elementUniwersalny.setSciany();

        System.out.println("\nFunkcje kształtu elementu uniwersalnego:");
        for (int i=0; i<4; i++){
            System.out.println("\nsciana "+(i+1));
            elementUniwersalny.sciany[i].printSciana();
        }

        System.out.println("\nMacierz Hbc elementu uniwersalnego:");
        for (int i=0; i<4; i++){
            System.out.println("Sciana "+(i+1));
            elementUniwersalny.sciany[i].obliczHbcSciany(alpha);
            elementUniwersalny.sciany[i].printHbcSciany();
            System.out.println();
        }

        siatka.obliczMacierzHbc(elementUniwersalny, alpha, t);
        for (int i=0; i<siatka.nE; i++){
            System.out.println("\nElement nr: "+(i+1));
            System.out.println("Macierz Hcb:");
            siatka.elements[i].printMacierzHBC();
            System.out.println("Wektor P ");
            siatka.elements[i].printP();
        }

        siatka.obliczHglobal();
        siatka.printHglobal();

        siatka.obliczHiHbcglobal();
        siatka.printHiHbcglobal();


        siatka.obliczPglobal();
        siatka.printPglobal();

//Obliczanie macierzy C

        int ro=7800; //gęstość
        int c=700;  //ciepło właściwe
        elementUniwersalny.printFunkcjeKsztaltu();
        for (int i=0; i<siatka.nE;i++){
            System.out.println("\nElement "+(i+1));
            siatka.elements[i].obliczMacierzC(elementUniwersalny, ro, c);
            siatka.elements[i].printMacierzC();
        }

        siatka.obliczCglobal();
        siatka.printCglobal();

        //symulacja

        /*
        matrix H = [H]+[C]/krok
        matrix P= [C]/krok * {t0} + {P}
        [matrix H]*{t}=[matrix P]
        */

        double t0=100.0;        //poczatkowa temperatura w węzłach
        int czasSymulacji =500; //100;
        int krok = 50;//1;      //po jakim kroku czasowym otrzymamy wynik?

        for (int i =0; i<siatka.nN; i++){
            siatka.nodes[i].setTempertature(t0); //wypełnienie siatki temperaturami początkowymi
        }

        /*
        siatka.matrixH(krok);
        siatka.printMatrixH();
        siatka.matrixP(krok);
        siatka.printMatrixP();
*/
        siatka.symulacja(czasSymulacji, krok);

        /*
        // GAUSS
        //dane do schematu dwupunktowego
        double wagi2[]={1.0,1.0};
        double pkt2[]={-1.0/Math.sqrt(3.0),1 /Math.sqrt(3.0) };
        Schemat dwupunktowy=new Schemat(wagi2,pkt2,2);
        dwupunktowy.schematPrint();

        //dane do schematu trzypunktowego
        double wagi3[]={(5.0/9.0),(8.0/9.0), (5.0/9.0)};
        double pkt3[]={-Math.sqrt(3.0/5.0), 0, Math.sqrt(3.0/5.0)};
        Schemat trzypunktowy=new Schemat(wagi3,pkt3,3);
        trzypunktowy.schematPrint();

        //Gauss

        Gauss dwupunktowe =new Gauss(dwupunktowy);
        Gauss trzypunktowe =new Gauss(trzypunktowy);

        System.out.format("Gauss 1d, schemat 2 punktowy: %.3f\t",dwupunktowe.calculateGauss1d());

        System.out.format("\nGauss 1d, schemat 3 punktowy: %.3f\t",trzypunktowe.calculateGauss1d());

        System.out.format("\nGauss 2d, schemat 2 punktowy: %.3f\t",dwupunktowe.calculateGauss2d());

        System.out.format("\nGauss 2d, schemat 3 punktowy: %.3f\t",trzypunktowe.calculateGauss2d());
        */


    }
}

