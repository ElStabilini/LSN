int main(){
  J = 1;
  int npoint = 10;
  double T [npoint  + 1];

  // Caso h = 0
  h = 0;
  //Misuro le osservabili per 10 valori di temperatura nell'intervallo [0.5, 2.0]
  double dt = 1.5/npoint;
  for(int i = 0; i <= npoint; i++){
    T[i] = 0.5 + i*dt;
  }

  // i = 0 --> Gibbs
  // i = 1 --> Metropolis
  for(int i = 0; i < 2; i++){
    metro = i;
    //Per ogni valore di temperatura, misuro le grandezze richieste, dopo aver equilibrato la simulazione
    for(auto t:T){
      Input(h, t);
      // Equilibrazione di 10000 passi
      nblk = 10000;
      for(int i = 0; i < nblk; i++){
        Move(metro);
      }

      //Eseguo l'algoritmo per 20 blocchi di 1000 step
      nblk = 20;
      nstep = 40000;
      for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
      {
        Reset(iblk);
        for(int istep=1; istep <= nstep; ++istep)
        {
          Move(metro);
          Measure();
          Accumulate();
        }
        Averages(iblk);
      }
      Averages(t, h);
      
    }
  }
  
  cout << endl;
  
}


