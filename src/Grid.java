import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

public class Grid {
   double H;                        //wysokość siatki
   double B;                        //szerokość siatki
   int nH;                          //liczba węzłów na wysokość
   int nB;                          //liczba węzłów na szerokość
   public int nN;                   //liczba węzłów
   public  int nE;                  //liczba elementów
   public Element [] elements;      //tablica z elementami
   public Node [] nodes;            //tablica z węzłami
   Element_uniwersalny uniwersalnyElement;

   //macierze globalne
   double Hglobal[][];
   double HiHbcGlobal[][];
   double Pglobal[];
   double Cglobal[][];

   //macierze do układu równań
   double matrixH[][];  //matrix H = C/krok + H
   double matrixP[];    //matrix P = C/krok * tot + P

   public Grid(double H, double B, int nH, int nB){
       this.H=H;
       this.B=B;
       this.nH=nH;
       this.nB=nB;
       nN=nH*nB;
       nE=(nH-1)*(nB-1);
       elements = new Element[nE];
       nodes = new Node[nN];

       Hglobal = new double [nN][nN];
       HiHbcGlobal = new double [nN][nN];
       Pglobal=new double [nN];
       Cglobal= new double[nN][nN];
       matrixH=new double[nN][nN];
       matrixP=new double[nN];
   }
   //budowa siatki
   public void addElement(){
       int j=1; //numeracja zaczyna się od 1
           for (int i=0; i<=nE-1;i++) {
               if(j%nH==0) j++;
               elements[i]=new Element(j, nH);
               System.out.println("ID Elementu nr "+(i+1)+": ");
               elements[i].print();
               j++;
           }
   }
   public void addNode(){
       double x=0; //ile odleglosci na osi x
       double y=0; //ile odleglosci na osi y
       double dX=B/(nB-1);
       double dY=H/(nH-1);
       for (int i=0; i<=(nN-1);i++) { //iterator tablicy z węzłami
           if(y%nH==0&&y!=0) {
               x++;
               y=0;
           }
           nodes[i]=new Node(x*dX, y*dY);

        if(nodes[i].x==0||nodes[i].x==B||nodes[i].y==0||nodes[i].y==H){ //sprwdzenie obecności warunku brzegowego w węzłach
            nodes[i].setBC(1);
        }
           int a=i+1;
           System.out.println("Współrzędne wezła nr "+a+": ");
           nodes[i].print();
           y++;
       }
   }
   //jakobian do 1d Hbc, P
   public double obliczDetJ(Node a, Node b){
       return (sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2)))/2.0;
   }
   //zapisywanie HBC w elementach
   public void obliczMacierzHbc(Element_uniwersalny uniwersalnyElement, int alpha, int tot){
       double DetJ=0;
       for(int i=0; i<nE;i++) {
           if (nodes[(elements[i].ID[0]) - 1].BC == 1 && nodes[(elements[i].ID[1]) - 1].BC == 1) { //sprawdzamy czy na danej ścianie elementu jest warunek brzegowy
               DetJ = obliczDetJ(nodes[(elements[i].ID[1]) - 1], nodes[(elements[i].ID[0]) - 1]); //jesli tak obliczamy DetJ czyli długość ściany podzielona przez 2.

               for (int j = 0; j < 4; j++) {
                   for (int k = 0; k < 4; k++) {
                       elements[i].macierzhbc[j][k] += uniwersalnyElement.sciany[0].macierzHbcSciany[j][k] * DetJ;
                   }
               }
               for(int k=0;k<uniwersalnyElement.k;k++){ //dla wszystkich punktow całkowania
                   for(int x=0; x<4;x++) {
                       elements[i].P[x]+=alpha*uniwersalnyElement.sciany[0].wagi[k]*uniwersalnyElement.sciany[0].funkcjeKsztaltu[k][x]* tot *DetJ;
                   }
               }
           }
           if (nodes[(elements[i].ID[1]) - 1].BC == 1 && nodes[(elements[i].ID[2]) - 1].BC == 1) {
               DetJ = obliczDetJ(nodes[(elements[i].ID[2]) - 1], nodes[(elements[i].ID[1]) - 1]);

               for (int j = 0; j < 4; j++) {
                   for (int k = 0; k < 4; k++) {
                       elements[i].macierzhbc[j][k] += uniwersalnyElement.sciany[1].macierzHbcSciany[j][k] * DetJ;
                   }
               }
               for(int k=0;k<uniwersalnyElement.k;k++){ //dla dwoch punktow całkowania
                   for(int x=0; x<4;x++) {
                       elements[i].P[x]+=alpha*uniwersalnyElement.sciany[1].wagi[k]*uniwersalnyElement.sciany[1].funkcjeKsztaltu[k][x]* tot *DetJ;
                   }
               }
           }
           if (nodes[(elements[i].ID[2]) - 1].BC == 1 && nodes[(elements[i].ID[3]) - 1].BC == 1) {
               DetJ = obliczDetJ(nodes[(elements[i].ID[3]) - 1], nodes[(elements[i].ID[2]) - 1]);

               for (int j = 0; j < 4; j++) {
                   for (int k = 0; k < 4; k++) {
                       elements[i].macierzhbc[j][k] += uniwersalnyElement.sciany[2].macierzHbcSciany[j][k] * DetJ;
                   }
               }
               for(int k=0;k<uniwersalnyElement.k;k++){ //dla dwoch punktow całkowania
                   for(int x=0; x<4;x++) {
                       elements[i].P[x]+=alpha*uniwersalnyElement.sciany[2].wagi[k]*uniwersalnyElement.sciany[2].funkcjeKsztaltu[k][x]* tot *DetJ;
                   }
               }
           }
           if (nodes[(elements[i].ID[3]) - 1].BC == 1 && nodes[(elements[i].ID[0]) - 1].BC == 1) {
               DetJ = obliczDetJ(nodes[(elements[i].ID[0]) - 1], nodes[(elements[i].ID[3]) - 1]);

               for (int j = 0; j < 4; j++) {
                   for (int k = 0; k < 4; k++) {
                       elements[i].macierzhbc[j][k] += uniwersalnyElement.sciany[3].macierzHbcSciany[j][k] * DetJ;
                   }
               }
               for(int k=0;k<uniwersalnyElement.k;k++){ //dla dwoch punktow całkowania
                   for(int x=0; x<4;x++) {
                       elements[i].P[x]+=alpha*uniwersalnyElement.sciany[3].wagi[k]*uniwersalnyElement.sciany[3].funkcjeKsztaltu[k][x]* tot *DetJ;
                   }
               }
           }
       }
   }
   //agregacje
   public void obliczHglobal(){
       for (int i = 0; i < nE; i++) {
            for (int x = 0; x < 4; x++) {
                for (int y = 0; y < 4; y++) {
                    Hglobal[(elements[i].ID[x])-1][(elements[i].ID[y])-1] += elements[i].macierzH[x][y];
                }
            }
        }
   }
   public void obliczHiHbcglobal(){ //H + Hbc
        for (int i = 0; i < nE; i++) {
            for (int x = 0; x < 4; x++) {
                for (int y = 0; y < 4; y++) {
                    HiHbcGlobal[(elements[i].ID[x])-1][(elements[i].ID[y])-1] += elements[i].macierzH[x][y] +elements[i].macierzhbc[x][y];
                }
            }
        }
    }
   public void obliczPglobal(){
       for(int i = 0; i < nE; i++){
           for(int x = 0; x < 4; x++){
               Pglobal[(elements[i].ID[x])-1] += elements[i].P[x];
           }
       }
   }
   public void obliczCglobal(){
        for (int i = 0; i < nE; i++) {
            for (int x = 0; x < 4; x++) {
                for (int y = 0; y < 4; y++) {
                    Cglobal[(elements[i].ID[x])-1][(elements[i].ID[y])-1] += elements[i].macierzC[x][y] ;
                }
            }
        }
   }
   //wypisywanie
   public void printCglobal(){

        System.out.println("\nMacierz globalna C\n");
        for (int i = 0; i < nN; i++) {
            for (int j = 0; j < nN; j++) {
                System.out.format("%.3f\t", Cglobal[i][j]);
            }
            System.out.println();
        }
    }
   public void printHglobal(){

       System.out.println("\nMacierz globalna H \n");
           for (int i = 0; i < nN; i++) {
               for (int j = 0; j < nN; j++) {
                   System.out.format("%.3f\t", Hglobal[i][j]);
               }
               System.out.println();
           }
   }
   public void printHiHbcglobal(){

        System.out.println("\nMacierz globalna H+Hbc\n");
        for (int i = 0; i < nN; i++) {
            for (int j = 0; j < nN; j++) {
                System.out.format("%.3f\t", HiHbcGlobal[i][j]);
            }
            System.out.println();
        }
    }
   public void printPglobal(){
        System.out.println("\nWektor globalny P\n");
        for (int i = 0; i < nN; i++) {
                System.out.format("%.3f\t", Pglobal[i]);
            System.out.println();
        }
    }
   //obliczanie elementów układu równań
   public void matrixH(int krok){ //matrix H=[H]+[C]/krok
        for (int i=0; i<nN;i++){
            for(int j=0; j<nN;j++){
                matrixH[i][j]=HiHbcGlobal[i][j]+Cglobal[i][j]/ krok;
            }
        }
    }
   public void matrixP(int krok){ //matrix P={P}+{[C]/krok} * {T0}
       double [] CprzezKrokRazyT0 =new double[nN];

       for (int i=0; i<nN; i++){
           for (int j=0; j<nN;j++) {
               CprzezKrokRazyT0[i] += (Cglobal[i][j] / krok)*nodes[j].tempertature;
           }
       }
        for (int i=0; i<nN;i++){
            matrixP[i]=Pglobal[i]+ CprzezKrokRazyT0[i];
        }
    }
   public void printMatrixH(){ //matrix H=h+C/dtał
        System.out.println("\nMatrix H\n");
        for (int i=0; i<nN;i++){
            for(int j=0; j<nN;j++){
                System.out.format("%.3f\t", matrixH[i][j]);
            }
            System.out.println();
        }
    }
   public void printMatrixP(){
        System.out.println("\nMatrix P\n");
        for (int i=0; i<nN;i++){
            System.out.format("%.1f\t", matrixP[i]);
        }
        System.out.println("\n");
    }
   //symulacja
   public void symulacja(int czas, int krok){
       int iloscIteracji= czas / krok;
        for(int i = 0; i < iloscIteracji; i++){
            System.out.format("Iteracja %d, czas %d",i,((i* krok)+ krok) );

            matrixP(krok);
            matrixH(krok);

            double []temp =new double[nN];
            GaussElimination temperatures=new GaussElimination();
            double[][] A = new double [nN][nN];
            A=matrixH;
            double[] B = new double [nN];
            B=matrixP;
            temp=temperatures.solve(A,B);
            /*for (int j=0; j<nN;j++){
                System.out.format("%.2f\t", temp[j]);
            }*/
            for(int j = 0; j < nN; j++){
                nodes[j].tempertature = temp[j];
            }
            double min,max;
            min = temp[0];
            max = temp[0];

            for(int j = 0; j < nN; j++){
                if(max < temp[j])  max = temp[j];
                if(min > temp[j])  min = temp[j];
            }
            System.out.format("\nmin: %.2f \tmax: %.2f\n", min, max);
        }
    }
}
