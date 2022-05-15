
import java.util.*;
import static java.lang.Math.sqrt;

import java.lang.Math;
import java.util.Random;
import java.util.List;
import java.util.ArrayList;


public class MC_pricer {

    // variables
    double So;// underlying
    double K;// Strike_price
    int T;// Maturity
    double r;// rate
    double sigma;//volatility

    int N;//number of trajectories

    int nT;// Timestep

    char type;

    String barrierType;

    double barrier;

    public MC_pricer(double So, double K, int T, double r, double sigma,int N,int nT,char type,String barrierType,double barrier) {
        this.So = So;
        this.K = K;
        this.T = T;
        this.r = r;
        this.sigma = sigma;
        this.N = N;
        this.nT = nT;
        this.type=type;
        this.barrierType=barrierType;
        this.barrier=barrier;

    }

    public static List<double[]> generateGBM(double sigma, double r, double T, double So, int N, int nT) {
        Random rnd = new Random();
        double Deltat = T/nT;
        double sqrtDeltat = Math.sqrt(Deltat);
        double drift=r;
        double vol=sigma;
        List<double[]> res = new ArrayList<double[]>();
        for (int i=0; i < N; i++) {
            // create vector of simulations with first value So;
            double[] vectorS = new double[nT+1];
            vectorS[0] = So;
            for (int j=1; j < nT+1; j++) {
                vectorS[j] = vectorS[j-1]*(1+drift*Deltat+vol*rnd.nextGaussian()*sqrtDeltat);
            }
            res.add(vectorS);
        }
        return res;
    }


    public static void printList(List<double[]> l) {
        int n = l.size();
        int m = l.get(0).length;

        // loop print
        for (int i=0; i < n; i++) {
            double[] vectorS = l.get(i);
            System.out.printf("Array %d: ", i+1);
            for (int j=0; j < m; j++) {
                System.out.printf("%f ", vectorS[j]);
            }
            System.out.printf("\n");
        }
    }

    private static List<OptionalDouble> getMean(List<double[]> stock){
        int n = stock.size();
        int m = stock.get(0).length;

        double[][] array = stock.toArray(new double[0][0]);

        double[][] transpose = new double[m][n];
        for(int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                transpose[j][i] = array[i][j];
            }
        }
        List<OptionalDouble> temp = new ArrayList<OptionalDouble>();
        for(int i = 0; i < m; i++) {
            temp.add(Arrays.stream(transpose[i]).average());
        }

        return  temp;

    }
    private static Double barrierAnalysis(List<OptionalDouble> stock, String barrierType, double barrier) {
        boolean crossedUp = false;
        boolean crossedDown=false;


        for (OptionalDouble optionalDouble : stock) {

            if (optionalDouble.orElseThrow(IllegalStateException::new) >= barrier) {
                crossedUp = true;
                break;
            }

            }
        for (OptionalDouble optionalDouble : stock) {

            if (optionalDouble.orElseThrow(IllegalStateException::new) < barrier) {
                crossedDown = true;
                break;
            }

        }

        switch(barrierType) {
            case "UI":
                // UP-AND-IN
                if(crossedUp){
                    return stock.get(stock.size() - 1).orElseThrow(IllegalStateException::new);
            }
                break;
            case "UO":
                // UP-AND-OUT
                if(crossedUp){
                    return 0.0;
                }
                break;
            case "DI":
                // DOWN-AND-IN
                if(crossedDown){
                    return stock.get(stock.size() - 1).orElseThrow(IllegalStateException::new);
                }
                break;
            case "DO":
                // DOWN-AND-OUT
                if(crossedDown){
                    return 0.0;
                }
                break;
            default:
                return stock.get(stock.size() - 1).orElseThrow(IllegalStateException::new);

        }
        return stock.get(stock.size() - 1).orElseThrow(IllegalStateException::new);
    }


    public double computeMCPrice() {
        List<double[]> stock=generateGBM(sigma,r,T,So,N,nT);
        List<OptionalDouble> S=getMean(stock);
        double final_S=barrierAnalysis(S,barrierType, barrier);
        double price=0;
        if (type == 'c') price = Math.max(final_S - K* Math.exp(-r*T), 0);// price of a call option
        else price =  Math.max(K*Math.exp(-r*T)-final_S, 0);// price of a put option
        return price;
    }

    public double[] computeBSGreeks() {

        double dOne = (Math.log(So / K) + (r  + Math.pow(sigma, 2) / 2) * T) / (sigma * Math.sqrt(T));
        double dTwo = dOne - sigma * Math.sqrt(T);
        // compute Delta, Vega, Psi, Theta, Rho, Gamma and Volga
        double[] greeks = new double[3];
        if (type == 'c') {
            // Delta
            greeks[0] =  CUMNORMSDIST(dOne);
            // Vega
            greeks[1] = So *  normalDist(dOne) * Math.sqrt(T);
            // Gamma
            greeks[2] =  normalDist(dOne) / (So * sigma * Math.sqrt(T));

        } else {
            // Delta
            greeks[0] = -1 * CUMNORMSDIST(-dOne);
            // Vega
            greeks[1] = So  * normalDist(dOne) * Math.sqrt(T);
            // Gamma
            greeks[2] =  normalDist(dOne) / (So * sigma * Math.sqrt(T));

        }

        return greeks;
    }

    public static double normalDist(double x) {
        // probability density function of normal distribution
        return 1/Math.sqrt(2 * Math.PI) * Math.exp(-Math.pow(x,2)/2);
    }

    public static long factorial(int n) {
        long fact = 1;
        for (int i = 2; i <= n; i++) {
            fact = fact * i;
        }
        return fact;
    }

    public static double cumlativeNormalDist(double x, int n) {

        // CUMULATIVE NORMAL DISTRIBUTION
        double res=0.5+1/Math.sqrt(2 * Math.PI);
        double taylor_series_approx= x;

        for (int i = n - 1; i > 0; --i )
            taylor_series_approx = taylor_series_approx + (Math.pow(-1,n)*Math.pow(x,2*n+1))/(Math.pow(2,n)* factorial(n)*(2*n+1));

        return res + taylor_series_approx;
    }

    private static double erf(double x)
    {
        //A&S formula 7.1.26
        double a1 = 0.254829592;
        double a2 = -0.284496736;
        double a3 = 1.421413741;
        double a4 = -1.453152027;
        double a5 = 1.061405429;
        double p = 0.3275911;
        x = Math.abs(x);
        double t = 1 / (1 + p * x);
        //Direct calculation using formula 7.1.26 is absolutely correct
        //But calculation of nth order polynomial takes O(n^2) operations
        //return 1 - (a1 * t + a2 * t * t + a3 * t * t * t + a4 * t * t * t * t + a5 * t * t * t * t * t) * Math.Exp(-1 * x * x);

        //Horner's method, takes O(n) operations for nth order polynomial
        return 1 - ((((((a5 * t + a4) * t) + a3) * t + a2) * t) + a1) * t * Math.exp(-1 * x * x);
    }
    public static double CUMNORMSDIST(double z)
    {
        double sign = 1;
        if (z < 0) sign = -1;

        double result=0.5 * (1.0 + sign * erf(Math.abs(z)/Math.sqrt(2)));
        return result;
    }


}
