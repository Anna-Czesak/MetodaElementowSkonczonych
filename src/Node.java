public class Node {
    double x, y;
    int BC;
    double tempertature;

    public Node(double x, double y) {
        this.x = x;
        this.y = y;
        this.BC=0; //informacja o obecności warunku brzegowego
        this.tempertature=0;
    }
    public void setBC(int BC){
        this.BC=BC;
    }
    public void print(){
        System.out.format("x: %.3f y: %.3f BC: %d", x,y,BC);
        System.out.println();
    }
    public void setTempertature(double t){
        tempertature=t;
    }

}
