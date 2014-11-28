inline void countEandJ (MainType& Enew, MainType& Em, MainType* Jnew, MainType* Jm, MainType& rotdt, DispStruct& CfArr, int xc, int md=Md) {
  MainType Etmp = Enew;
  MainType sum=0.;
  for(int i=0; i<md; i++) sum+= CfArr.kEJ[xc][i][0]*Jnew[i] + CfArr.kEJ[xc][i][1]*Jm[i];  // в статье слагаемые также!!!
  Enew = CfArr.kE[xc][0]*Em + CfArr.kE[xc][1]*Enew + CfArr.kE[xc][2]*rotdt - sum;  //тут rotor*dt = (Dnew-Dold) !!!

  MainType Jtmp[md];
  for(int i=0; i<md; i++) {
    Jtmp[i] = Jnew[i];
    Jnew[i] = CfArr.kJ[xc][i][0]*Jnew[i] + CfArr.kJ[xc][i][1]*Jm[i] + CfArr.kJ[xc][i][2]*Enew + CfArr.kJ[xc][i][3]*Em + CfArr.kJ[xc][i][4]*Etmp;
  }

  Em = Etmp;
  for(int i=0; i<md; i++) Jm[i] = Jtmp[i];

  if(isnan(Enew)) {
  printf("\nEnew=%g\n",Enew);
  printf("Em=%g\n",Em);
  for(int i=0; i<md; i++) printf("Jnew=%g  i=%d\n",Jnew[i],i);
  for(int i=0; i<md; i++) printf("Jm=%g  i=%d\n",Jm[i],i);
  printf("rotdt=%g\n", rotdt);

  for(int i=0; i<md; i++) printf("i=%d  kEJ=%g  %g\n",i,CfArr.kEJ[xc][i][0], CfArr.kEJ[xc][i][1] );
  printf("kE=%g  %g  %g\n",CfArr.kE[xc][0], CfArr.kE[xc][1], CfArr.kE[xc][2] );
  for(int i=0; i<md; i++) printf("i=%d  kJ=%g  %g  %g  %g  %g\n",i,CfArr.kJ[xc][i][0], CfArr.kJ[xc][i][1], CfArr.kJ[xc][i][2], CfArr.kJ[xc][i][3], CfArr.kJ[xc][i][4]);

  exit(0);
  }
}
