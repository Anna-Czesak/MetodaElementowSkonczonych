import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

public class Element_uniwersalny {

    double []Ksi;
    double []Eta;
    int iloscPktCalkowania;
    double []wagi;
    int k;
    double funkcjeKsztaltu[][];

    public double[][]dNdKsi; //pochodne funkcji kształtu po ksi i eta
    public double[][]dNdEta;

    public double []jakobianPrzeksztalcenia=new double [2];
    public Sciana [] sciany=new Sciana[4];

    public Element_uniwersalny(int k){ //dwu lub trzypunktowy schemat
        if(k==2){ //wspolrzedne punktow calkowania 2D 2punktowy
            Ksi= new double[]{(-1.0 / sqrt(3.0)), (1.0 / sqrt(3.0)), (1.0 / sqrt(3.0)), (-1.0 / sqrt(3.0))};
            Eta= new double[]{(-1.0/sqrt(3.0)), (-1.0/sqrt(3.0)), (1.0/sqrt(3.0)), (1.0/sqrt(3.0))};
            wagi=new double[]{1.0,1.0};
            iloscPktCalkowania = 4;
            this.k=k;
        }
        if(k==3){ //wspolrzedne punktow calkowania 2D 3punktowy
            Ksi= new double[]{(-sqrt(3.0 / 5.0)), 0.0, sqrt(3.0 / 5.0), (-sqrt(3.0 / 5.0)), 0.0, sqrt(3.0 / 5.0), (-sqrt(3.0 / 5.0)), 0.0, sqrt(3.0 / 5.0)};
            Eta= new double[]{(-sqrt(3.0 / 5.0)), -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0.0 , 0.0, 0.0, sqrt(3.0 / 5.0), sqrt(3.0 / 5.0), sqrt(3.0 / 5.0)};
            wagi=new double[]{5.0/9.0, 8.0/9.0, 5.0/9.0 };
            iloscPktCalkowania =9;
            this.k=k;

        }

        funkcjeKsztaltu=new double[iloscPktCalkowania][4];

        dNdKsi=new double[iloscPktCalkowania][4];//-1/4(1-Eta)
        dNdEta=new double[iloscPktCalkowania][4];//-1/4(1+Ksi)


        for (int i = 0; i < iloscPktCalkowania; i++) {
            for (int j = 0; j < 4; j++) {
                if (j == 0) dNdEta[i][j] = (-1.0 / 4.0) * (1.0 - (Ksi[i]));
                if (j == 1) dNdEta[i][j] = (-1.0 / 4.0) * (1.0 + (Ksi[i]));
                if (j == 2) dNdEta[i][j] = (1.0 / 4.0) * (1.0 + (Ksi[i]));
                if (j == 3) dNdEta[i][j] = (1.0 / 4.0) * (1.0 - (Ksi[i]));
            }
        }
        for (int i = 0; i < iloscPktCalkowania; i++) {
            for (int j=0; j<4;j++){
                if (j==0) dNdKsi[i][j]=(-1.0 / 4.0) * (1.0 - (Eta[i]));
                if (j==1) dNdKsi[i][j]=(1.0 / 4.0) * (1.0 - (Eta[i]));
                if (j==2) dNdKsi[i][j]= (1.0 / 4.0) * (1.0 + (Eta[i]));
                if (j==3) dNdKsi[i][j]=(-1.0 / 4.0) * (1.0 + (Eta[i]));
            }
        }
        //funkcje kształtu do macierzy C
        for (int i = 0; i < iloscPktCalkowania; i++) {
            funkcjeKsztaltu[i][0] = 0.25 * (1.0 - Ksi[i]) * (1.0 - Eta[i]);
            funkcjeKsztaltu[i][1] = 0.25 * (1.0 + Ksi[i]) * (1.0 - Eta[i]);
            funkcjeKsztaltu[i][2] = 0.25 * (1.0 + Ksi[i]) * (1.0 + Eta[i]);
            funkcjeKsztaltu[i][3] = 0.25 * (1.0 - Ksi[i]) * (1.0 + Eta[i]);
        }
    }

    public void printdNdEta(){
        System.out.println("pochodna funkcji kształtu względem ETA: \n\t\tN1\t\t\tN2\t\t\tN3\t\t\tN4");
        for (int i=0; i <iloscPktCalkowania;i++){
            System.out.print("pc "+(i+1)+"\t");
            for (int j=0; j<4; j++){
                System.out.format("%.5f\t ",dNdEta[i][j]);
            }
            System.out.println() ;
        }
    }
    public void printdNdKsi(){
    System.out.println("\npochodna funkcji kształtu względem KSI: \n\t\tN1\t\t\tN2\t\t\tN3\t\t\tN4");
        for (int i=0; i <iloscPktCalkowania;i++){
            System.out.print("pc "+(i+1)+"\t");
            for (int j=0; j<4; j++){
            System.out.format("%.5f\t ",dNdKsi[i][j]) ;
            }
        System.out.println() ;
        }
    }
    public void setSciany(){

        Sciana sciana1;
        Sciana sciana2;
        Sciana sciana3;
        Sciana sciana4;

        if (k==2) {
        //schemat całkowania dla ściany 1
        double wagi1[] = {1.0, 1.0};
        double ksi1[] = {-1.0 / Math.sqrt(3.0), 1 / Math.sqrt(3.0)};
        double eta1[] = {-1.0, -1.0};
        sciana1 = new Sciana(ksi1, eta1, wagi1, this.k);
        //schemat całkowania dla ściany 2
        double wagi2[] = {1.0, 1.0};
        double ksi2[] = {1.0, 1.0};
        double eta2[] = {-1.0 / Math.sqrt(3.0), 1 / Math.sqrt(3.0)};
        sciana2 = new Sciana(ksi2, eta2, wagi2,this.k);
        //schemat całkowania dla ściany 3
        double wagi3[] = {1.0, 1.0};
        double ksi3[] = {1.0 / Math.sqrt(3.0), -1 / Math.sqrt(3.0)};
        double eta3[] = {1.0, 1.0};
        sciana3 = new Sciana(ksi3, eta3, wagi3,this.k);
        //schemat całkowania dla ściany 4
        double wagi4[] = {1.0, 1.0};
        double ksi4[] = {-1.0, -1.0};
        double eta4[] = {1.0 / Math.sqrt(3.0), -1 / Math.sqrt(3.0)};
        sciana4 = new Sciana(ksi4, eta4, wagi4,this.k);

        sciany[0]=sciana1;
        sciany[1]=sciana2;
        sciany[2]=sciana3;
        sciany[3]=sciana4;
        }
        if (k==3) {
            //schemat całkowania dla ściany 1
            double wagi1[] = {5.0/9.0, 8.0/9.0,5.0/9.0};
            double ksi1[] = {(-sqrt(3.0 / 5.0)), 0.0, sqrt(3.0 / 5.0)};
            double eta1[] = {-1.0, -1.0,-1.0};
            sciana1 = new Sciana(ksi1, eta1, wagi1,k);
            //schemat całkowania dla ściany 2
            double wagi2[] = {5.0/9.0, 8.0/9.0,5.0/9.0};
            double ksi2[] = {1.0, 1.0,1.0};
            double eta2[] = {(-sqrt(3.0 / 5.0)), 0.0, sqrt(3.0 / 5.0)};
            sciana2 = new Sciana(ksi2, eta2, wagi2,k);
            //schemat całkowania dla ściany 3
            double wagi3[] = {5.0/9.0, 8.0/9.0,5.0/9.0};
            double ksi3[] = {(sqrt(3.0 / 5.0)), 0.0, -sqrt(3.0 / 5.0)};
            double eta3[] = {1.0, 1.0, 1.0};
            sciana3 = new Sciana(ksi3, eta3, wagi3,k);
            //schemat całkowania dla ściany 4
            double wagi4[] = {5.0/9.0, 8.0/9.0,5.0/9.0};
            double ksi4[] = {-1.0, -1.0, -1.0};
            double eta4[] = {(sqrt(3.0 / 5.0)), 0.0, -sqrt(3.0 / 5.0)};
            sciana4 = new Sciana(ksi4, eta4, wagi4,k);

            sciany[0]=sciana1;
            sciany[1]=sciana2;
            sciany[2]=sciana3;
            sciany[3]=sciana4;
        }
    }
    public void printFunkcjeKsztaltu(){
        System.out.println("\nFunkcje kształtu elemenru uniwersalnego\n");
        for (int i = 0; i < iloscPktCalkowania; i++) {
            System.out.print("pc: " + (i + 1) + "\t ");
            System.out.format("ksi: %.4f\teta: %.4f \n", Ksi[i], Eta[i]);
            for (int j=0; j<4;j++){
                System.out.format("%.4f\t", funkcjeKsztaltu[i][j]);
            }
            System.out.println();
        }
    }
}
