public class Sciana {
    int k; //schemat
    public double[] ksi; //pierwsza wspolrzedna pkt calkowania
    public double[] eta; //druga wspolrzedna pkt całkowania
    public double[] wagi;
    public double[][] funkcjeKsztaltu;// = new double[k][4]; //po 4 funkcje kształtu dla kazdego punktu całkowania
    public double[][] macierzHbcSciany =new double [4][4];

    public Sciana(double[] ksi, double[] eta, double[] wagi, int k) {
        this.k=k;
        this.ksi = ksi;
        this.eta = eta;
        this.wagi = wagi;
        funkcjeKsztaltu = new double[this.k][4];
        for (int i = 0; i < this.k; i++) {
            funkcjeKsztaltu[i][0] = 0.25 * (1.0 - this.ksi[i]) * (1.0 - this.eta[i]);
            funkcjeKsztaltu[i][1] = 0.25 * (1.0 + this.ksi[i]) * (1.0 - this.eta[i]);
            funkcjeKsztaltu[i][2] = 0.25 * (1.0 + this.ksi[i]) * (1.0 + this.eta[i]);
            funkcjeKsztaltu[i][3] = 0.25 * (1.0 - this.ksi[i]) * (1.0 + this.eta[i]);
        }
    }

    public void printSciana() {
        for (int i = 0; i < this.k; i++) {
            System.out.format("Pc %d (%.3f, %.3f) \t", (i + 1), ksi[i], eta[i]);
            for (int j = 0; j < 4; j++) {
                System.out.format(" %.3f \t", funkcjeKsztaltu[i][j]);

            }
            System.out.println();
        }
    }

    public void obliczHbcSciany(int alpha){ //mnozenie * DetJ dopiero pezy zapisywaniu do elementu
        double [][]hbcLocal=new double[4][4];
        double [][]hbc=new double[4][4];
       for (int i=0; i<k; i++){ //obliczamy hbc w kazdym punkcie całkowania
           for (int j = 0; j < 4; j++)
               for (int k = 0; k < 4; k++) {
                   hbcLocal[j][k] = wagi[i]*funkcjeKsztaltu[i][j] * funkcjeKsztaltu[i][k];
                   hbc[j][k] += hbcLocal[j][k]; //sumowanie hbc w kazdym punkcie całkowania
               }
       }
       for (int i=0; i<4;i++){
           for(int j=0; j<4;j++){
               hbc[i][j]*=alpha; //mnożenie ostatecznej macierzy hbc ściany przez współczynnik konwekcji cieplnej
           }
       }
        macierzHbcSciany =hbc;
    }

    public void printHbcSciany(){
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                System.out.format("%.3f\t", macierzHbcSciany[i][j]); //DetJ=L/2
            }
            System.out.println();
        }
    }
}
