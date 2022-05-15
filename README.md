# MonteCarloPricer

this an assignment to find the monte carlo option prices, whether theyre vanilla or not.
 I created an executable jar folder under the artifacts folder.
 The code takes the following as input:
 
 
    double So;// underlying
    
    double K;// Strike_price
    
    int T;// Maturity
    
    double r;// rate
    
    double sigma;//volatility

    int N;//number of trajectories

    int nT;// Timestep

    char type; // c for call and p for. put 

    String barrierType;// UI,UO,DI,DO taken as string, if the user wishes for the option to be vanilla then he has to input any other string 

    double barrier;// the barrier needed if the option is knock in or knock out 
