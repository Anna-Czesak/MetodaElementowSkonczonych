public class Gauss {

    public static double f1(double x){
        return 5*x*x+3*x+6;
    }
    public static double f2(double x, double y){
        return 5*x*x*y*y+3*x*y+6;
    }

    Schemat schemat;

    public Gauss(Schemat schemat ){ //k-2 dla dwupunktowego,
        this.schemat = new Schemat();
        this.schemat=schemat;
    }

    public double calculateGauss1d(){
        double calka =0;
        for (int i = 0; i < schemat.k; i++){
            calka += schemat.wagi[i] * f1(schemat.punkty[i]);
        }

        return calka;
    }

    public double calculateGauss2d(){
            double calka =0;
            for (int i = 0; i < schemat.k; i++){
                for (int j = 0; j < schemat.k; j++){
                    calka += schemat.wagi[i]* schemat.wagi[j] * f2(schemat.punkty[i], schemat.punkty[j]);
                }
            }

            return calka;
    }
}

