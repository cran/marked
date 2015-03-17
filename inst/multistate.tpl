// Fixed-effect Multi-State Cormack-Jolly-Seber model with unobservable states based on work of Jessica Ford in MEE Dec 2012
// Jeff Laake; 12 Jan 2013
DATA_SECTION 
    init_int n;                            // number of capture histories
    init_int m;                            // number of capture occasions
	init_int nS;                           // number of states excluding death state
	init_int uS;                           // number of unobserved states
    init_imatrix ch(1,n,1,m);              // capture history matrix; uses numeric values for states
    init_ivector frst(1,n);                // occasion first seen for each history
	init_vector freq(1,n);                 // frequency of each capture history
    init_matrix tint(1,n,1,m-1);           // time interval between occasions for each history-interval
    init_int kphi;                         // number of columns in the design matrix for Phi - survival
    int nrows;                             // number of entries in design matrix m-1 values for each of nS states
    int all_nrows;                         // number of rows in design matrix for all histories
    !! nrows=nS*(m-1);                
    !! all_nrows=n*nrows;                
    init_matrix phidm(1,all_nrows,1,kphi); // design matrix for Phi
    init_int kp;                           // number of columns in the design matrix for p - capture probability
    init_matrix pdm(1,all_nrows,1,kp);     // design matrix for p
	int nT;                                // number of transitions excluding death
	int all_nT;                            // number of transitions for all histories
	!! nT=nS*nS*(m-1);                       
	!! all_nT=n*nT;                       
	init_int kpsi;                         // number of columns in Psi design matrix
	init_matrix psidm(1,all_nT,1,kpsi);    // design matrix for transition probs
		
PARAMETER_SECTION
    init_vector phibeta(1,kphi);       // parameter vector for Phi
    init_vector pbeta(1,kp);           // parameter vector for p
    init_vector psibeta(1,kpsi);       // parameter vector for p
	objective_function_value g; 

PROCEDURE_SECTION
    int i;                             // index over observations
    for(i=1;i<=n;i++)                  // loop over capture histories - one per capture history
        ll_i(i,phibeta,pbeta,psibeta);
 
FUNCTION void ll_i(const int i, const dvar_vector& phibeta, const dvar_vector& pbeta, const dvar_vector& psibeta )
//  Code implements an algorithm described on pg 47 in Zucchini and MacDonald (Z&M) for computing likelihood value for a 
//  Hidden Markov Model. The algorithm is recursive over the encounter occasions. The algorithm involves a simpler recursion 
//  shown on pages 38 and 45 but to avoid numerical underflows it is scaled by the sum of the state probabilities which makes
//  it appear a little more complex.
    dvar_vector phi(1,nrows);          // temp vector for Phis for an occasion
    dvar_vector p(1,nrows);            // temp vector for ps for an occasion
    dvar_vector psiexp(1,nT);          // temp vector for psis 
    dvariable psisum;                  // sum of psiexp for each state
    dvariable u;                       // sum of state probabilities
    dvar_matrix psi(1,nS,1,nS);        // matrix for psis 
	int ni=m-frst(i)+1;                // length of capture history for individual i
	dvar_vector pS(1,nS+1);            // update vector for prob of being in state j=1,nS + death       
	dvar_vector S(1,nS+1);             // prob of being in state j=1,nS + death for each occasion
	dvar_vector pD(1,nS+1);            // observation probability; equivalent to diagonal matrix
	dvar_vector v(1,nS+1);             // temporary update vector
	dvar_vector pseen(1,nS+1);         // capture probability for each state 
	dvar_matrix gamma(1,nS+1,1,nS+1);  // transition matrix including death
	dvariable Lglki;                   // log-likelihood accumulator
    int j,k,k2;                        // lop variables
	int index,index2,bindex;           // indices
    pS.initialize();                   // initialize values to 0
    Lglki.initialize();	
	bindex=(i-1)*nrows;                                 // initialize index into phi,p for ith history
    for(j=1;j<=nrows;j++)
	   phi(j)=1/(1+exp(-phidm(bindex+j)*phibeta));      // compute phi for the interval
    for (j=1;j<=m-1;j++)                                // adjust phi for the time interval length if needed
   	for (k=1;k<=nS;k++)
       phi((j-1)*nS+k)=pow(phi((j-1)*nS+k),tint(i,j));    
    for(j=1;j<=nrows;j++)
	   p(j)=1/(1+exp(-pdm(bindex+j)*pbeta));  // compute p 

	bindex=(i-1)*nT;                                    // initialize index into psi for ith history
    for(j=1;j<=nT;j++)
       psiexp(j)=exp(psidm(bindex+j)*psibeta);          // compute exp of psidm*psibeta; these are components of psi

	S.initialize();                                     // initialize S to all 0
	S(ch(i,frst(i)))=1;                                 // set state prob to 1 for observed state at first observation

    for(j=2;j<=ni;j++)                                  // loop over remaining possible occasions from 2 to ni
    {
 	    int j2=j+frst(i)-1;
 	    index=(j2-2)*nS+1;                              // index into phi and p vectors for each occasion

	    gamma.initialize();                             // initialize all transitions to zero
 	    pseen.initialize();                             // initialize p to all 0
 	    index2=(j2-2)*nS*nS;                            // index in psi vector
	    for(k=1;k<=nS;k++)                              // loop over states computing psi matrix 
	    {
	       psisum=0;
	       for(k2=(k-1)*nS+1;k2<=k*nS;k2++)              
		     psisum=psisum+psiexp(index2+k2);
  	       for(k2=1;k2<=nS;k2++)
	          psi(k,k2)=psiexp((k-1)*nS+k2+index2)/psisum;
	    }
	    for(k=1;k<=nS;k++)                              // loop over states creating p and gamma values
		{
		   if(k <= (nS-uS))                             // get and assign p
		      pseen(k)=p(index);                        // observable state
           else
		      pseen(k)=0;                               // unobservable state
	       for(k2=1;k2<=nS;k2++)
		     gamma(k,k2)=psi(k,k2)*phi(index);          // adjust psi for survival
		   gamma(k,nS+1)=1-phi(index);                  // add death state value for each state
		   index++;
		}
		gamma(nS+1,nS+1)=1;                             // death is an absorbing state
        pS=S*gamma;        	                            // previous scaled state probability*transition matrix	
        pD.initialize();                                // probability of observation; equivalent to diagonal matrix in Z&M
		if(ch(i,j2)>0)                          
           pD(ch(i,j2))=pseen(ch(i,j2));
		else
		   pD=(1-pseen); 
        v=elem_prod(pS,pD);                             // v is temp state vector alpha in Z&M
    	u=sum(v);                                       // sum across states
        S=v/u;                                          // update S;S is normalized alpha(j) (phi) in Z&M
	    Lglki+=log(u);    	                            // accumulate log-likelihood value
	}
	g-=freq(i)*Lglki;
