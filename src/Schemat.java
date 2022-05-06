public class Schemat {
    public double wagi[];
    public double punkty[];
    public int k;
    public Schemat(){}
    public Schemat(double []wagi, double[]pkt, int k){
        this.wagi=new double[k];
        this.wagi=wagi;

        this.punkty=new double[k];
        this.punkty=pkt;

        this.k=k;
    }
    public void schematPrint(){
        System.out.println(k+"-punktowy schemat calkowania\nwagi:");
        for (int i=0; i<k;i++){
            System.out.format("%.3f\t",this.wagi[i]);
        }
        System.out.println("\npunkty calkowania:");
        for (int i=0; i<k;i++){
            System.out.format("%.3f\t",this.punkty[i]);
        }
        System.out.println("\n");
    }
}
