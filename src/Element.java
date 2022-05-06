import java.util.Arrays;

public class Element {

    int []ID=new int[4];

    double jakobian[][];
    double wyznacznik[];
    double odwrotnoscWyznacznika[];
    double jakobianOdwrotny[][];

    double dNdX[][];//[ilosc pkt calkowania][4];
    double dNdY[][];//[ilosc pkt calkowania][4];

    double macierzH[][];//new double[4][4];
    double macierzC[][]=new double[4][4];
    public double macierzhbc[][]=new double [4][4];
    double P[]=new double[4];

    public Element(int i, int nH){ //id to numer elementu w tablicy węzłów
        ID[0]=i;
        ID[1]=ID[0]+nH;
        ID[2]=ID[1]+1;
        ID[3]=ID[0]+1;
    }
    public void obliczJakobian(Element_uniwersalny e, Grid siatka, int ktoryelement){

        double dxDksi[]=new double[e.iloscPktCalkowania];
        double dyDksi[]=new double[e.iloscPktCalkowania];
        double dxDeta[]=new double[e.iloscPktCalkowania];
        double dyDeta[]=new double[e.iloscPktCalkowania];

        jakobian=new double[e.iloscPktCalkowania][4]; //dx/dksi, dy/dksi, dx/deta, dy/deta - w wierszu dla kazdego punktu całkowania

        Node pkt[]=new Node[4];
        pkt[0]=siatka.nodes[((siatka.elements[ktoryelement]).ID[0])-1];
        pkt[1]=siatka.nodes[((siatka.elements[ktoryelement]).ID[1])-1];
        pkt[2]=siatka.nodes[((siatka.elements[ktoryelement]).ID[2])-1];
        pkt[3]=siatka.nodes[((siatka.elements[ktoryelement]).ID[3])-1];

        for(int i=0; i<4;i++){
            pkt[i].print();
        }

        for (int j = 0; j < e.iloscPktCalkowania; j++) { //w każdym punkcie całkowania liczymy jakobian
            for (int i = 0; i < 4; i++) {
                dxDksi[j] += e.dNdKsi[j][i] * pkt[i].x;
                dyDksi[j] += e.dNdKsi[j][i] * pkt[i].y;
                dxDeta[j] += e.dNdEta[j][i] * pkt[i].x;
                dyDeta[j] += e.dNdEta[j][i] * pkt[i].y;
            } //zapisanie poszczególnych elementów do jakobianów
            jakobian[j][0]=dxDksi[j];
            jakobian[j][1]=dyDksi[j];
            jakobian[j][2]=dxDeta[j];
            jakobian[j][3]=dyDeta[j];
        }
    }
    public void obliczWyznacznik(Element_uniwersalny e){
        wyznacznik=new double [e.iloscPktCalkowania];
        odwrotnoscWyznacznika=new double[e.iloscPktCalkowania];
        for(int i=0; i<e.iloscPktCalkowania;i++){ //dla kazdego punktu calkowania
            wyznacznik[i]=jakobian[i][0]*jakobian[i][3]-jakobian[i][1]*jakobian[i][2];
            odwrotnoscWyznacznika[i]=1/wyznacznik[i];
        }
    }
    public void obliczJakobianOdwrotny(Element_uniwersalny e){
        jakobianOdwrotny=new double[e.iloscPktCalkowania][4];
        for (int i=0; i<e.iloscPktCalkowania;i++){
            jakobianOdwrotny[i][0] = jakobian[i][3] * odwrotnoscWyznacznika[i];
            jakobianOdwrotny[i][1] = jakobian[i][1] * odwrotnoscWyznacznika[i] * (-1);
            jakobianOdwrotny[i][2] = jakobian[i][2] * odwrotnoscWyznacznika[i] * (-1);
            jakobianOdwrotny[i][3] = jakobian[i][0] * odwrotnoscWyznacznika[i];
        }
    }
    public void obliczdNdXinNnY(Element_uniwersalny e){
        dNdX=new double[e.iloscPktCalkowania][4];
        for (int j=0; j<e.iloscPktCalkowania;j++) { //dla 4/9 punktów całkowania
            for (int i = 0; i < 4; i++) {
                dNdX[j][i] = jakobianOdwrotny[j][0] * e.dNdKsi[j][i] + jakobianOdwrotny[j][1] * e.dNdEta[j][i];
            }
        }
        dNdY=new double[e.iloscPktCalkowania][4];
        for (int j=0; j<e.iloscPktCalkowania;j++) { //dla 4/9 punktów całkowania
            for (int i = 0; i < 4; i++) {
                dNdY[j][i] = jakobianOdwrotny[j][2] * e.dNdKsi[j][i] + jakobianOdwrotny[j][3] * e.dNdEta[j][i];
            }
        }
    }

    public void printJakobian(Element_uniwersalny e){
        System.out.println("\njakobian przekształcenia:");
        for (int i=0; i<e.iloscPktCalkowania;i++){ //w kazdym punkcie całkowania
            System.out.println("W punkcie calkowania nr: "+(i+1));
            System.out.println("Jakobian [J]: ");
            System.out.format("%.4f\t %.4f\n", jakobian[i][0],jakobian[i][1]);
            System.out.format("%.4f\t %.4f\n", jakobian[i][2],jakobian[i][3]);
            System.out.format("Det[J]: %f\n", wyznacznik[i]);
            System.out.format("1/Det[J]: %.2f\n", odwrotnoscWyznacznika[i]);
            System.out.println("Jakobian odwrotny [J]: ");
            System.out.format("%.4f\t %.4f\n",jakobianOdwrotny[i][0],jakobianOdwrotny[i][1]);
            System.out.format("%.4f\t %.4f\n\n", jakobianOdwrotny[i][2],jakobianOdwrotny[i][3]);
        }
    }

    // [H]=kt*(dN/dx*dN/dx^T+dN/dy*dN/dy^T)*detJ
    public void obliczMacierzH(int kt, Element_uniwersalny e){
        macierzH=new double[4][4];
        double macierzHlocal[][]=new double[4][4];
        double iloczynx[][]=new double[4][4];
        double iloczyny[][]=new double[4][4];
        double suma[][]=new double[4][4];
        for(int i=0; i<e.iloscPktCalkowania;i++){
            for(int j=0;j<4;j++){
                for(int k=0;k<4;k++){
                    iloczynx[j][k]=dNdX[i][j]*dNdX[i][k];
                    iloczyny[j][k]=dNdY[i][j]*dNdY[i][k];
                }
            }
            for(int j=0;j<4;j++){
                for(int k=0;k<4;k++){
                    suma[j][k]=iloczynx[j][k]+iloczyny[j][k];
                }
            }
            for(int j=0;j<4;j++){
                for(int k=0;k<4;k++){
                     int y = i / e.k;
                     int x = i % e.k;
                    macierzHlocal[j][k] =  (e.wagi[x] * e.wagi[y])*kt * suma[j][k] * wyznacznik[i];
                }
            }
            for(int j=0;j<4;j++){
                for(int k=0;k<4;k++){
                    macierzH[j][k]+=macierzHlocal[j][k];
                }
            }
        }
    }

    // [C]=ro*cp*N*N^T*detJ,
    public void obliczMacierzC(Element_uniwersalny ele, int ro, int c){
        double cLocal[][]=new double[4][4];
        System.out.println("Obliczanie macierzy C:");
        for (int k=0; k< ele.iloscPktCalkowania; k++){
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    int x=k % ele.k;
                    int y=k / ele.k;

                    cLocal[i][j] += ro*c*ele.wagi[x]*ele.wagi[y] * ele.funkcjeKsztaltu[k][i] * ele.funkcjeKsztaltu[k][j]*wyznacznik[k];
                    macierzC[i][j] = cLocal[i][j];
                }
            }
        }
    }

    //wypisywanie
    public void printMacierzC(){
        System.out.println("\nMacierz [C]");
        for (int i=0; i< 4;i++){
            for(int j=0; j<4;j++){
                System.out.format("%.5f\t", macierzC[i][j]);
            }
            System.out.println();
        }
    }
    public void printPochodneDNdxDNdy(Element_uniwersalny e){
        System.out.println("Pochodne dN/dx: ");
        for (int i=0; i<e.iloscPktCalkowania;i++){
            System.out.print("\npc "+(i+1)+"\t");
            for (int j=0;j<4;j++){
                System.out.format("%.6f\t", dNdX[i][j]);
            }
        }
        System.out.println("\n\nPochodne dN/dy: ");
        for (int i=0; i<e.iloscPktCalkowania;i++){
            System.out.print("\npc "+(i+1)+"\t ");
            for (int j=0;j<4;j++){
                System.out.format("%.6f\t", dNdY[i][j]);
            }
        }
    }
    public void print(){
        System.out.println(Arrays.toString(ID));
    }
    public void printMacierzH(){
        System.out.println("\nMacierz [H]");
        for (int i=0; i< 4;i++){
            for(int j=0; j<4;j++){
                System.out.format("%.3f\t", macierzH[i][j]);
            }
            System.out.println();
        }
    }
    public void printMacierzHBC(){
       for (int i = 0; i < 4; ++i) {
           for (int j = 0; j < 4; ++j) {
               System.out.format("%.3f\t", macierzhbc[i][j]);
           }
           System.out.println();
       }
   }
    public void printP(){
        for (int i=0; i<4; i++){
            System.out.format("%.3f\t", P[i]);
        }
        System.out.println();
}

}
