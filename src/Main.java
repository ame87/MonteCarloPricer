import java.util.Arrays;
import java.util.Scanner;


public class Main {

    public static void main(String[] args) {
        System.out.println("This is a Monte Carlo Pricer please enter the varaibles needed");
        Scanner myObj = new Scanner(System.in);

        System.out.println("Enter initial stock price:");
        double So = myObj.nextDouble();

        System.out.println("Enter Strike price:");
        double K = myObj.nextDouble();

        System.out.println("Enter time to maturity (must be int):");
        int T = myObj.nextInt();

        System.out.println("Enter interest rate:");
        double r = myObj.nextDouble();

        System.out.println("Enter volatility:");
        double sigma = myObj.nextDouble();

        System.out.println("Enter number of iterations:");
        int N = myObj.nextInt();

        System.out.println("Enter timestep (int):");
        int nT = myObj.nextInt();

        System.out.println("Enter type as a char, can be c (call) or p(put):");
        char type = myObj.next().charAt(0);

        Scanner scanner = new Scanner(System.in);
        System.out.println("Enter Barrier type, can be (UO,UI,DO,DI) or you anything else if Vanilla option:");
        String barrierType = scanner.nextLine();

        System.out.println("Enter Barrier:");
        double barrier = scanner.nextDouble();




        MC_pricer thing=new MC_pricer(So,K,T,r,sigma,N,nT,type,barrierType,barrier);
        double X=thing.computeMCPrice();
        System.out.printf("The price is  %f\n", X);

        double[] greeks= thing.computeBSGreeks();
        System.out.printf("The greeks are Delta, Vega and Gamma respectively:\n");
        System.out.println(Arrays.toString(greeks));

    }


}
